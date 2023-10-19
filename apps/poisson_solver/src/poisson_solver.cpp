/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <iostream>
#include <vector>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>
#include <map>

#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/hho"
#include "mumps.hpp"
#include "diskpp/output/silo.hpp"

#include "diskpp/methods/implementation_hho/curl.hpp"

#include "diskpp/common/timecounter.hpp"

#include "lua_interface.h"
#include "poisson_solver.h"

#include "sgr.hpp"



template<disk::mesh_2D Mesh>
void
adjust_stabfree_recdeg(const Mesh& msh, const typename Mesh::cell_type& cl,
    disk::hho_degree_info& hdi)
{
    size_t cd = hdi.cell_degree();
    size_t fd = hdi.face_degree();
    bool is_mixed_high = (hdi.cell_degree() > hdi.face_degree());
    size_t n = faces(msh, cl).size();   
    size_t rpd = cd+2;

    /* HHO space dofs */
    size_t from = ((cd+2)*(cd+1))/2 + n*(fd+1);
    /* Reconstruction dofs */
    size_t to = ((rpd+2)*(rpd+1))/2;

    if (from <= to) {
        hdi.reconstruction_degree(rpd);
    }
    else {
        /* Every harmonic degree provides 2 additional dofs, therefore
         * we need an increment that it is sufficient to accomodate
         * (from-to) dofs => ((from - to) + (2-1))/2 */
        size_t incr = (from - to + 1)/2;
        hdi.reconstruction_degree(rpd+incr);
    }
}

template<typename Mesh, typename ProblemData>
auto
compute_laplacian_operator(hho_poisson_solver_state<Mesh>& state,
    typename Mesh::cell_type& cl, const ProblemData& pd)
{
    using mesh_type = Mesh;
    using T = typename mesh_type::coordinate_type;

    auto& msh = state.msh;
    auto& hdi = state.hdi;

    bool is_mixed_high = (hdi.cell_degree() > hdi.face_degree());

    disk::dynamic_matrix<T> A;
    disk::dynamic_matrix<T> GR;

    auto subdomain_id = msh.domain_info(cl).tag(); 
    constexpr size_t DIM = Mesh::dimension;
    disk::static_matrix<T, DIM, DIM> diff_tens
        = pd.diffusion_coefficient(subdomain_id, barycenter(msh, cl)).tensor();

    if (state.use_stabfree)
    {
        if (is_mixed_high) {
            auto oper = make_shl_face_proj_harmonic(msh, cl, hdi, diff_tens);
            A = oper.second;
            GR = oper.first;
        }
        else {
            auto oper = make_sfl(msh, cl, hdi, diff_tens);
            A = oper.second;
            GR = oper.first;
        }
        state.A_norm += A.norm();
    } else {
        auto oper = make_scalar_hho_laplacian(msh, cl, hdi, diff_tens);
        A = oper.second;
        GR = oper.first;
        state.A_norm += A.norm();
        if (is_mixed_high) {
            if (state.use_dt_stab)
                A = A + make_scalar_hdg_stabilization2(msh, cl, hdi, diff_tens);
            else
                A = A + make_scalar_hdg_stabilization(msh, cl, hdi);
        }
        else {
            if (state.use_dt_stab)
                A = A + make_scalar_hho_stabilization2(msh, cl, GR, hdi, diff_tens);
            else
                A = A + make_scalar_hho_stabilization(msh, cl, GR, hdi);
        }
        state.AS_norm += A.norm();
    }
    
    return std::pair(GR, A);
}

template<typename Mesh, typename ProblemData>
auto
compute_element_contribution(hho_poisson_solver_state<Mesh>& state,
    typename Mesh::cell_type& cl, const ProblemData& pd)
{
    using mesh_type = Mesh;
    using T = typename mesh_type::coordinate_type;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using point_type = typename Mesh::point_type;

    auto& msh = state.msh;
    auto& hdi = state.hdi;

    if (state.use_stabfree)
        adjust_stabfree_recdeg(msh, cl, hdi);

    auto bar = barycenter(msh, cl);
    auto subdomain_id = msh.domain_info(cl).tag();
    constexpr size_t DIM = Mesh::dimension;
    disk::static_matrix<T, DIM, DIM> diff_tens
        = pd.diffusion_coefficient(subdomain_id, barycenter(msh, cl)).tensor();

    auto [GR, lhs] = compute_laplacian_operator(state, cl, pd);

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();
    auto cbs = disk::scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    auto fbs = disk::scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    Eigen::Matrix<T, Eigen::Dynamic, 1> strong_data =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(num_faces*fbs);

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> robin_lhs =
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(cbs + num_faces*fbs, cbs + num_faces*fbs);

    Eigen::Matrix<T, Eigen::Dynamic, 1> weak_data =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(cbs + num_faces*fbs);

    Eigen::Matrix<T, Eigen::Dynamic, 1> wddofs =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(cbs + num_faces*fbs);

    size_t fnum = 0;
    for (auto& fc : fcs)
    {
        auto n = normal(msh, cl, fc);
        auto gfnum = offset(msh, fc);
        auto b = state.boundary_info[gfnum];
        auto fb = make_scalar_monomial_basis(msh, fc, hdi.face_degree());

        auto dirichlet_fun = [&](const point_type& pt) {
            return pd.dirichlet_data(b.number, pt);
        };

        auto neumann_fun = [&](const point_type& pt) {
            return pd.neumann_data(b.number, pt, disk::point<T,DIM>(n));
        };

        if (b.type == boundary_type::dirichlet)
        {
            strong_data.segment(fnum*fb.size(), fb.size()) = 
                disk::project_function(msh, fc, fb, dirichlet_fun);
        }

        if (b.type == boundary_type::neumann) {
            weak_data.segment(cbs + fnum*fb.size(), fb.size()) =
                disk::make_rhs(msh, fc, fb, neumann_fun); 
        }

        if (b.type == boundary_type::robin) {
            auto fofs = cbs + fnum*fb.size();
            T bnd_val = b.value.value_or(n.dot(diff_tens*n));
            robin_lhs.block(fofs, fofs, fb.size(), fb.size()) +=
                bnd_val * disk::make_mass_matrix(msh, fc, fb);
        }

        if (b.type == boundary_type::jump /* and cell is jump source */) {
            wddofs.segment(cbs + fnum*fb.size(), fb.size()) +=
                disk::project_function(msh, fc, fb, dirichlet_fun);
            weak_data.segment(cbs + fnum*fb.size(), fb.size()) +=
                disk::make_rhs(msh, fc, fb, neumann_fun);
        }

        fnum++;
    }

    lhs = lhs + robin_lhs;

    auto rhs_fun = [&](const point_type& pt) {
        return pd.right_hand_side(subdomain_id, pt);
    };

    auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
    Eigen::Matrix<T, Eigen::Dynamic, 1> rhs =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(lhs.rows());
    rhs.head(cb.size()) = make_rhs(msh, cl, cb, rhs_fun);

    rhs += -lhs*wddofs + weak_data;

    return std::tuple(lhs, rhs, strong_data);
}


template<typename Mesh, typename ProblemData>
int
assemble_stabfree(hho_poisson_solver_state<Mesh>& state, const ProblemData& pd)
{
    using mesh_type = Mesh;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using point_type = typename Mesh::point_type;

    auto& msh = state.msh;
    auto& hdi = state.hdi;
    auto& assm = state.assm;

    state.recdegs.resize(msh.cells_size());

    size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto cbs = disk::scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
        auto [lhs, rhs, dd] = compute_element_contribution(state, cl, pd);
        auto [lhsC, rhsC] = disk::static_condensation(lhs, rhs, cbs);
        assm.assemble(msh, cl, lhsC, rhsC, dd);
        state.recdegs[cell_i] = hdi.reconstruction_degree();
        cell_i++;
    }
    assm.finalize();

    return 0;
}

template<typename Mesh, typename ProblemData>
int
assemble_standard(hho_poisson_solver_state<Mesh>& state, const ProblemData& pd)
{
    using mesh_type = Mesh;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using point_type = typename Mesh::point_type;

    auto& msh = state.msh;
    auto& hdi = state.hdi;
    auto& assm = state.assm;

    state.recdegs.resize(msh.cells_size(), hdi.reconstruction_degree());

    for (auto& cl : msh)
    {
        auto cbs = disk::scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
        auto [lhs, rhs, dd] = compute_element_contribution(state, cl, pd);
        auto [lhsC, rhsC] = disk::static_condensation(lhs, rhs, cbs);
        assm.assemble(msh, cl, lhsC, rhsC, dd);
    }
    assm.finalize();

    return 0;
}

template<typename Mesh, typename ProblemData>
void
assemble(hho_poisson_solver_state<Mesh>& state, const ProblemData& pd)
{
    auto& msh = state.msh;
    auto& hdi = state.hdi;
    auto& assm = state.assm;

    state.A_norm = 0.0;
    state.AS_norm = 0.0;

    std::vector<bool> is_dirichlet;
    is_dirichlet.resize( msh.faces_size() );

    for (size_t i = 0; i < msh.faces_size(); i++) {
        const auto& b = state.boundary_info[i];
        is_dirichlet[i] = false;

        if (b.type == boundary_type::dirichlet)
            is_dirichlet[i] = true;

        if (b.type == boundary_type::undefined and not b.internal)
            is_dirichlet[i] = true;
    }

    assm.initialize(state.msh, state.hdi, is_dirichlet);

    if (state.use_stabfree) {
        std::cout << "Assembling stabilization-free problem" << std::endl;
        assemble_stabfree(state, pd);
    }
    else {
        std::cout << "Assembling standard HHO problem" << std::endl;
        assemble_standard(state, pd);
    }

    double S_norm = state.AS_norm - state.A_norm;
    std::cout << "A norm  = " << state.A_norm << std::endl;
    std::cout << "S norm  = " << S_norm << " (" << 100.0*S_norm/state.A_norm << "%)" << std::endl;
    std::cout << "AS norm = " << state.AS_norm << std::endl;
}

template<typename Mesh, typename ProblemData>
void
solve(hho_poisson_solver_state<Mesh>& state, const ProblemData& pd)
{
    using mesh_type = Mesh;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using scalar_type = typename Mesh::coordinate_type;

    auto& msh = state.msh;
    auto& assm = state.assm;
    auto& sol = state.sol;
    auto& hdi = state.hdi;

    sol = disk::dynamic_vector<scalar_type>::Zero(assm.syssz);
    std::cout << "Running MUMPS..." << std::flush;
    sol = mumps_lu(assm.LHS, assm.RHS);
    std::cout << "done" << std::endl;

    std::cout << "Expanding solution..." << std::flush;
    auto cd = hdi.cell_degree();
    auto fd = hdi.face_degree();
    auto cbs = disk::scalar_basis_size(cd, Mesh::dimension);
    auto fbs = disk::scalar_basis_size(fd, Mesh::dimension-1);
    auto fullsz = cbs*msh.cells_size() + fbs*msh.faces_size();

    state.sol_full = disk::dynamic_vector<scalar_type>::Zero(fullsz);

    size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto [lhs, rhs, dd] = compute_element_contribution(state, cl, pd);
        auto cb = disk::make_scalar_monomial_basis(msh, cl, hdi.cell_degree());

        auto edofs = assm.get_element_dofs(msh, cl, sol);
        edofs += dd;

        Matrix<scalar_type, Dynamic, 1> esol = disk::static_decondensation(lhs, rhs, edofs);
    
        state.sol_full.segment(cbs*cell_i, cbs) = esol.head(cbs);

        auto fcs = faces(msh, cl);
        for (size_t iF = 0; iF < fcs.size(); iF++) {
            auto& fc = fcs[iF];
            auto lofs = cbs + fbs*iF;
            auto gofs = cbs*msh.cells_size() + fbs*offset(msh, fc);
            state.sol_full.segment(gofs, fbs) = esol.segment(lofs, fbs);
        }
        cell_i++;
    }
    std::cout << "done" << std::endl;
}

template<typename Mesh, typename ProblemData, typename SolutionData>
auto
check(hho_poisson_solver_state<Mesh>& state, const ProblemData& pd,
    const SolutionData& sd)
{
    using mesh_type = Mesh;
    using point_type = typename mesh_type::point_type;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using scalar_type = typename Mesh::coordinate_type;
    using dv = disk::dynamic_vector<scalar_type>;
    using dm = disk::dynamic_matrix<scalar_type>;

    auto& msh = state.msh;
    auto& assm = state.assm;
    auto& sol = state.sol;
    auto& hdi = state.hdi;

    auto cd = hdi.cell_degree();
    auto fd = hdi.face_degree();
    auto cbs = disk::scalar_basis_size(cd, Mesh::dimension);
    auto fbs = disk::scalar_basis_size(fd, Mesh::dimension-1);
    auto fullsz = cbs*msh.cells_size() + fbs*msh.faces_size();

    scalar_type L2normsq = 0.0;
    scalar_type L2errsq = 0.0;
    scalar_type H1normsq = 0.0;
    scalar_type H1errsq = 0.0;
    scalar_type Aerrsq = 0.0;
    scalar_type Anormsq = 0.0;

    

    size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto domain_num = msh.domain_info(cl).tag();
        auto cb = disk::make_scalar_monomial_basis(msh, cl, hdi.cell_degree());

        auto sol_fun = [&](const point_type& pt) {
            return sd.sol(domain_num, pt);
        };

        auto sol_grad = [&](const point_type& pt) {
            return sd.grad(domain_num, pt);
        };

        auto fcs = faces(msh, cl);
        size_t num_dofs = cbs + fbs*fcs.size();
        
        dv num_loc_sol = dv::Zero(num_dofs);

        num_loc_sol.head(cbs) = state.sol_full.segment(cbs*cell_i, cbs);
        for (size_t iF = 0; iF < fcs.size(); iF++) {
            auto& fc = fcs[iF];
            auto lofs = cbs + fbs*iF;
            auto gofs = cbs*msh.cells_size() + fbs*offset(msh, fc);
            num_loc_sol.segment(lofs, fbs) = state.sol_full.segment(gofs, fbs);
        }

        dv ana_loc_sol = disk::project_function(msh, cl, hdi, sol_fun, 2);
        dv err_loc_sol = ana_loc_sol - num_loc_sol;

        dv ana_uT = ana_loc_sol.head(cbs);
        dv err_uT = err_loc_sol.head(cbs);

        dm MM = disk::make_mass_matrix(msh, cl, cb);

        L2errsq  += err_uT.dot(MM*err_uT);
        L2normsq += ana_uT.dot(MM*ana_uT);

        auto [GR, A] = compute_laplacian_operator(state, cl, pd);
        Aerrsq  += err_loc_sol.dot(A*err_loc_sol);
        Anormsq += ana_loc_sol.dot(A*ana_loc_sol);

        dm uR = GR*num_loc_sol;
        
        auto rb = disk::make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());

        auto subdomain_id = msh.domain_info(cl).tag();
        constexpr size_t DIM = Mesh::dimension;
        disk::static_matrix<scalar_type, DIM, DIM> diff_tens
            = pd.diffusion_coefficient(subdomain_id, barycenter(msh, cl)).tensor();

        auto qps = disk::integrate(msh, cl, 2*hdi.reconstruction_degree());
        for (auto& qp : qps) {
            disk::static_vector<scalar_type,DIM> ana_grad = sol_grad(qp.point());

            auto dphi = rb.eval_gradients(qp.point());
            disk::static_vector<scalar_type,DIM> num_grad =
                dphi.block(1,0,rb.size()-1,DIM).transpose()*uR;

            disk::static_vector<scalar_type,DIM> diff_grad = num_grad - ana_grad;

            H1normsq += qp.weight() * ana_grad.dot(diff_tens*ana_grad);
            H1errsq += qp.weight() * diff_grad.dot(diff_tens*diff_grad);
        }

        cell_i++;
    }

    std::cout << "Norm: " << std::sqrt(L2normsq) << std::endl;

    return std::tuple(std::sqrt(L2errsq/L2normsq),
        std::sqrt(H1errsq/H1normsq),
        std::sqrt(Aerrsq/Anormsq)
    );
}

template<typename Mesh, typename ProblemData>
void
export_to_visit(hho_poisson_solver_state<Mesh>& state,
    const ProblemData& pd, const char *filename)
{
    auto& msh = state.msh;
    auto& sol = state.sol_full;
    auto& hdi = state.hdi;

    auto cbs = disk::scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    std::vector<double> u;
    std::vector<double> rhs;
    size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto bar = barycenter(msh, cl);
        u.push_back( sol(cbs*cell_i) );
        rhs.push_back( pd.right_hand_side(1, bar) );
        cell_i++;
    }

    disk::silo_database db;
    db.create(filename);
    db.add_mesh(msh, "mesh");
    db.add_variable("mesh", "u", u, disk::zonal_variable_t);
    db.add_variable("mesh", "rhs", rhs, disk::zonal_variable_t);
    db.add_variable("mesh", "rd", state.recdegs, disk::zonal_variable_t);
}

template<typename Mesh>
void
collect_boundary_info(sol::state& lua, hho_poisson_solver_state<Mesh>& state)
{
    auto& msh = state.msh;
    auto& boundary_info = state.boundary_info;
    boundary_info.reserve( msh.faces_size() );

    for (auto& fc : faces(msh))
    {
        auto bi = msh.boundary_info(fc);

        using T = typename hho_poisson_solver_state<Mesh>::scalar_type;
        boundary_condition_descriptor<T> bcd;
        bcd.type = boundary_type::not_a_boundary;
        bcd.internal = false;
        bcd.number = 0;

        if (bi.is_boundary()) {
            bcd = lua_get_boundary_condition<T>(lua, bi.tag());
            bcd.internal = bi.is_internal();
        }

        boundary_info.push_back(bcd);
    }
}

#include "diskpp/mesh/point_lua_binding.hpp"

template<typename Mesh>
int
setup_and_run(sol::state& lua, hho_poisson_solver_state<Mesh>& state)
{
    constexpr size_t DIM = Mesh::dimension;
    lua[NODE_NAME_SIM]["dimension"] = DIM;

    collect_boundary_info(lua, state);
    lua_problem_data lpd(lua);
    lua_solution_data lsd(lua);

    register_point_usertype(lua, typename Mesh::point_type{});

    /* Ok, here the assumption is that state and lua are alive during the
     * execution of setup_and_run, so all the following lambdas are safe.
     * For extra caution, we unregister all the functions after calling
     * the user code, so they don't get called after a call to run() on
     * the Lua side. */

    /* Bind functions */
    lua["assemble"] = [&]() {
        return assemble(state, lpd);
    };

    lua["solve"] = [&]() {
        return solve(state, lpd);
    };

    lua["export_to_visit"] = [&](const char *fn) {
        return export_to_visit(state, lpd, fn);
    };

    lua["check"] = [&]() {
        return check(state, lpd, lsd);
    };

    lua["mesh_h"] = [&]() {
        return disk::average_diameter(state.msh);
    };

    /* Call user code */
    int err = lua_call_user_code(lua);
    
    /* Unbind everything to avoid unsafe calls */
    lua["assemble"] = nullptr;
    lua["solve"] = nullptr;
    lua["export_to_visit"] = nullptr;
    lua["check"] = nullptr;
    lua["mesh_h"] = nullptr;
    
    return err;
}

template<typename Mesh>
int
init_solver_state(sol::state& lua, hho_poisson_solver_state<Mesh>& state)
{
    hho_variant variant = lua_get_hho_variant(lua);
    int order = lua_get_hho_order(lua);

    if (order < 0 and variant != hho_variant::mixed_order_low) {
        std::cout << "Invalid HHO order " << order << ", reverting ";
        std::cout << "to order 0" << std::endl;
        order = 0;
    }

    if (order < 1 and variant == hho_variant::mixed_order_low) {
        std::cout << "Invalid HHO order " << order << ", reverting ";
        std::cout << "to order 1" << std::endl;
        order = 1;
    }

    /* Reconstruction degree will be changed if stabfree is used. */
    state.hdi.reconstruction_degree(order+1);
    state.hdi.face_degree(order);
    switch(variant) {
        case hho_variant::mixed_order_low:
            state.hdi.cell_degree(order-1);
            break;
        case hho_variant::equal_order:
            state.hdi.cell_degree(order);
            break;
        case hho_variant::mixed_order_high:
            state.hdi.cell_degree(order+1);
            break;
    }

    state.use_stabfree = lua_use_stabfree_hho(lua);
    state.use_dt_stab = lua_use_diffusion_tensor_in_stab(lua);

    return 0;
}

template<typename T>
int run(sol::state& lua)
{
    /* The user has called run() on the Lua side, so we got called.
     * Here we just prepare the solver state according to the mesh
     * requested by the user, nothing more.
     * The solver state is invisible from the Lua side (maybe this
     * will change), but the important thing is that it will be alive
     * (together with lua) during setup_and_run(). The binding to the
     * other functions happens in setup_and_run(), when we know exactly
     * the actual mesh type. */
    auto mps = lua_get_mesh_parameters(lua);

    if (mps.source == mesh_source::internal)
    {
        if (mps.type == internal_mesh_type::triangles)
        {
            using mesh_type = disk::simplicial_mesh<T,2>;
            hho_poisson_solver_state<mesh_type> state;
            init_solver_state(lua, state);
            auto mesher = disk::make_simple_mesher(state.msh);
            for (size_t i = 0; i < mps.level; i++)
                mesher.refine();
            /*
            std::random_device rd; 
            std::mt19937 gen(rd());
            auto h = disk::average_diameter(state.msh);
            std::uniform_real_distribution<> dis(-0.1*h, 0.1*h);

            auto tr = [&](const disk::point<T,2>& pt) {
                return disk::point<T,2>(pt.x()+dis(gen), pt.y()+dis(gen));
            };
            state.msh.transform(tr);
            */
            return setup_and_run(lua, state);
        }

        if (mps.type == internal_mesh_type::quadrangles)
        {
            using mesh_type = disk::cartesian_mesh<T,2>;
            hho_poisson_solver_state<mesh_type> state;
            init_solver_state(lua, state);
            auto mesher = disk::make_simple_mesher(state.msh);
            for (size_t i = 0; i < mps.level; i++)
                mesher.refine();
            return setup_and_run(lua, state);
        }

        if (mps.type == internal_mesh_type::hexagons)
        {
            using mesh_type = disk::generic_mesh<T,2>;
            hho_poisson_solver_state<mesh_type> state;
            init_solver_state(lua, state);
            auto mesher = disk::make_fvca5_hex_mesher(state.msh);
            mesher.make_level(mps.level);
            return setup_and_run(lua, state);
        }
        /*
        if (mps.type == internal_mesh_type::tetrahedra)
        {
            using mesh_type = disk::simplicial_mesh<T,3>;
            hho_poisson_solver_state<mesh_type> state;
            init_solver_state(lua, state);
            auto mesher = disk::make_simple_mesher(state.msh);
            for (size_t i = 0; i < mps.level; i++)
                mesher.refine();
            return setup_and_run(lua, state);
        }
        */
    }

    if (mps.source == mesh_source::file)
    {
        std::cout << "file" << std::endl;        
        auto mesh_filename = lua_mesh_filename(lua);
        
#ifdef HAVE_GMSH
        /* GMSH 2D simplicials */
        if (std::regex_match(mesh_filename, std::regex(".*\\.geo2s$") ))
        {
            using mesh_type = disk::simplicial_mesh<T,2>;
            hho_poisson_solver_state<mesh_type> state;
            init_solver_state(lua, state);
            std::cout << "Guessed mesh format: GMSH 2D simplicials" << std::endl;
            disk::gmsh_geometry_loader<mesh_type> loader;
        
            loader.read_mesh(mesh_filename);
            loader.populate_mesh(state.msh);
            return setup_and_run(lua, state);
        }
#else
        std::cout << "GMSH support not compiled. Exiting." << std::endl;
        return 1;
#endif
    }

    return 1;
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <param filename>" << std::endl;
        return 1;
    }

    using T = double;

    sol::state lua;
    lua_init_environment(lua);

    using dp2d = diffusion_parameter<T,2>;
    using dp3d = diffusion_parameter<T,3>;

    lua.new_usertype<dp2d>("tensor_2D",
		sol::constructors<dp2d(), dp2d(T)>(),
        "entry",
        sol::overload(
            sol::resolve<void(int, int, T)>(&dp2d::entry),
            sol::resolve<T(int, int) const>(&dp2d::entry)
        )
    );

    lua.new_usertype<dp3d>("tensor_3D",
		sol::constructors<dp3d(), dp3d(T)>(),
        "entry",
        sol::overload(
            sol::resolve<void(int, int, T)>(&dp3d::entry),
            sol::resolve<T(int, int) const>(&dp3d::entry)
        )
    );

    lua.new_usertype<timecounter>("timecounter",
		sol::constructors<timecounter()>(),
        "tic", &timecounter::tic,
        "toc", &timecounter::toc,
        "elapsed", &timecounter::elapsed
    );

    /* The code is going to be a bit spaghetti in the initial part.
     * Here we bind the run function above to the Lua environment. */
    lua["run"] = [&](){ return run<T>(lua); };

    /* and here we load the Lua script. If the user does not call run(),
     * nothing happens. However, when run() is called on the Lua side,
     * the run() on the C++ side is called and the ball starts rolling. */
    lua_load_script(lua, argv[1]);

    return 0;
}