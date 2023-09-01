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

#include "sgr.hpp"

















template<typename Mesh>
class hho_assembler
{
    using coordinate_type = typename Mesh::coordinate_type;
    using scalar_type = coordinate_type;
    using mesh_type = Mesh;
    using face_type = typename Mesh::face_type;
    using cell_type = typename Mesh::cell_type;

    //Mesh                                msh;
    disk::hho_degree_info               hdi;

    std::vector<Triplet<scalar_type>>   triplets;
    std::vector<bool>                   is_dirichlet;

    std::vector<size_t>                 compress_table;
    std::vector<size_t>                 expand_table;

    const size_t INVALID_OFFSET = (size_t) ~0;

    size_t face_basis_size(void) const {
        return disk::scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    }

    /* Get the offset of a face in the linear system */
    size_t get_system_offset(const Mesh& msh,
        const typename Mesh::face_type& fc)
    {
        auto face_num = offset(msh, fc);
        assert(face_num < msh.faces_size());
        assert(compress_table.size() == msh.faces_size());
        auto cnum = compress_table[face_num];
        assert(cnum != INVALID_OFFSET);
        auto fbs = face_basis_size();
        return cnum*fbs;
    }

    /* Determine if a face should be assembled */
    bool is_in_system(const Mesh& msh,
        const typename Mesh::face_type& fc)
    {
        auto bi = msh.boundary_info(fc);
        auto face_num = offset(msh, fc);
        assert(face_num < msh.faces_size());
        assert(is_dirichlet.size() == msh.faces_size());
        return not (bi.is_boundary() and is_dirichlet[face_num]);
    }

    void make_tables(mesh_type& msh)
    {
        compress_table.resize( msh.faces_size(), INVALID_OFFSET);
        expand_table.resize( sysfcs );

        size_t face_i = 0;
        size_t compressed_ofs = 0;
        for (auto& fc : faces(msh))
        {
            assert(compressed_ofs <= face_i);
            if ( is_in_system(msh, fc) )
            {
                assert(face_i < compress_table.size());
                compress_table[face_i] = compressed_ofs;
                assert(compressed_ofs < expand_table.size());
                expand_table[compressed_ofs] = face_i;
                compressed_ofs++;
            }

            face_i++;
        }

        assert(face_i == msh.faces_size());
    }

public:
    
    SparseMatrix<scalar_type>           LHS;
    Matrix<scalar_type, Dynamic, 1>     RHS;

    size_t                              syssz;
    size_t                              sysfcs;


    hho_assembler()
    {}

    hho_assembler(Mesh& msh, const disk::hho_degree_info& hdi,
                  const std::vector<bool>& is_dirichlet)
    {
        initialize(msh, hdi, is_dirichlet);
    }

    void clear()
    {
        triplets.clear();
        is_dirichlet.clear();
        compress_table.clear();
        expand_table.clear();
        syssz = 0;
        sysfcs = 0;
    }

    void initialize(Mesh& msh, const disk::hho_degree_info& p_hdi,
                    const std::vector<bool>& p_is_dirichlet)
    {
        clear();

        is_dirichlet = p_is_dirichlet;
        hdi = p_hdi;

        auto fbs = face_basis_size();

        auto in_system = [&](const face_type& fc) -> bool {
            auto ofs = offset(msh, fc);
            assert(ofs < is_dirichlet.size());
            return not (msh.is_boundary(fc) and is_dirichlet[ofs]);
        };

        sysfcs = std::count_if(msh.faces_begin(), msh.faces_end(), in_system);
        syssz = fbs*sysfcs;

        make_tables(msh);

        LHS = SparseMatrix<scalar_type>(syssz, syssz);
        RHS = Matrix<scalar_type, Dynamic, 1>::Zero(syssz);

        std::cout << "Assembler initialized: " << sysfcs << " faces in system, ";
        std::cout << syssz << " DoFs" << std::endl;
    }

    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<scalar_type, Dynamic, Dynamic>& lhsc,
             const Matrix<scalar_type, Dynamic, 1>& rhs,
             const Matrix<scalar_type, Dynamic, 1>& dirichlet_data)
    {
        auto fbs = face_basis_size();

        auto fcs = faces(msh, cl);
        for (size_t fi = 0; fi < fcs.size(); fi++)
        {
            if ( not is_in_system(msh, fcs[fi]) )
                continue;

            auto cofsi = get_system_offset(msh, fcs[fi]);
            for (size_t fj = 0; fj < fcs.size(); fj++)
            {
                auto lofsi = fi*fbs;
                auto lofsj = fj*fbs;

                if ( not is_in_system(msh, fcs[fj]) ) 
                {
                    RHS.segment(cofsi, fbs) += -lhsc.block(lofsi, lofsj, fbs, fbs)*dirichlet_data.segment(lofsj, fbs);
                    continue;
                }

                auto cofsj = get_system_offset(msh, fcs[fj]);
                for (size_t i = 0; i < fbs; i++)
                    for(size_t j = 0; j < fbs; j++)
                        triplets.push_back( Triplet<scalar_type>(cofsi+i, cofsj+j, lhsc(lofsi+i, lofsj+j)) );
            }

            RHS.segment(cofsi, fbs) += rhs.segment(fi*fbs, fbs);
        }
    }

    void finalize(void)
    {
        LHS.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();
        std::cout << "Matrix has " << LHS.nonZeros() << " nonzeros." << std::endl; 
    }

    disk::dynamic_vector<scalar_type>
    get_expanded_solution(const Mesh& msh, disk::dynamic_vector<scalar_type>& sol)
    {
        auto fbs = face_basis_size();

        disk::dynamic_vector<scalar_type> ret = 
            disk::dynamic_vector<scalar_type>::Zero( fbs*msh.faces_size() );

        for (size_t i = 0; i < sysfcs; i++)
        {
            auto in_offset = i*fbs;
            auto out_offset = expand_table.at(i)*fbs;
            ret.segment(out_offset, fbs) = sol.segment(in_offset, fbs);
        }

        return ret;
    }

    disk::dynamic_vector<scalar_type>
    get_element_dofs(const Mesh& msh, const typename Mesh::cell& cl, disk::dynamic_vector<scalar_type>& sol)
    {
        auto fbs = face_basis_size();
        auto fcs = faces(msh, cl);
        disk::dynamic_vector<scalar_type> ret = 
            disk::dynamic_vector<scalar_type>::Zero( fbs*fcs.size() );

        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto& fc = fcs[i];
            if ( not is_in_system(msh, fc) )
                continue;

            auto ofs = get_system_offset(msh, fc);
            ret.segment(i*fbs, fbs) = sol.segment(ofs, fbs);
        }

        return ret;
    }
};




















template<typename Mesh>
struct hho_poisson_solver_state
{
    using mesh_type = Mesh;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using scalar_type = typename Mesh::coordinate_type;
    using assembler_type = hho_assembler<Mesh>;
    using vector_type = disk::dynamic_vector<scalar_type>;

    mesh_type               msh;
    assembler_type          assm;
    vector_type             sol;
    vector_type             sol_full;

    std::vector<boundary>   boundary_info;

    disk::hho_degree_info   hdi;
    hho_variant             variant;
    bool                    use_stabfree;
    bool                    use_dt_stab;

    std::vector<double>     recdegs;
};

template<disk::mesh_2D Mesh>
void
adjust_stabfree_recdeg(const Mesh& msh, const typename Mesh::cell_type& cl,
    disk::hho_degree_info& hdi)
{
    size_t cd = hdi.cell_degree();
    size_t fd = hdi.face_degree();
    size_t n = faces(msh, cl).size();

    /* HHO space dofs */
    size_t from = ((cd+2)*(cd+1))/2 + n*(fd+1);
    /* Reconstruction dofs, polynomial part (degree is cd+2) */
    size_t to = ((cd+4)*(cd+3))/2;

    if (from <= to) {
        hdi.reconstruction_degree(cd+2);
    }
    else {
        /* Every harmonic degree provides 2 additional dofs, therefore
         * we need an increment that it is sufficient to accomodate
         * (from-to) dofs => ((from - to) + (2-1))/2 */
        size_t incr = (from - to + 1)/2;
        hdi.reconstruction_degree(cd+2+incr);
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

    constexpr size_t DIM = Mesh::dimension;
    disk::static_matrix<T, DIM, DIM> diff_tens
        //= disk::static_matrix<T, DIM, DIM>::Identity();
        = pd.diffusion_coefficient(1, barycenter(msh, cl)).tensor();

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
    } else {
        auto oper = make_scalar_hho_laplacian(msh, cl, hdi, diff_tens);
        A = oper.second;
        GR = oper.first;
    
        if (is_mixed_high) {
            A = A + make_scalar_hdg_stabilization(msh, cl, hdi);
        }
        else {
            if (state.use_dt_stab)
                A = A + make_scalar_hho_stabilization2(msh, cl, GR, hdi, diff_tens);
            else
                A = A + make_scalar_hho_stabilization(msh, cl, GR, hdi);
        }
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
    auto domain_num = msh.domain_info(cl).tag(); 

    auto [GR, A] = compute_laplacian_operator(state, cl, pd);

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();
    auto fbs = disk::scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    Eigen::Matrix<T, Eigen::Dynamic, 1> dirichlet_data =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(num_faces*fbs);

    size_t fnum = 0;
    for (auto& fc : fcs)
    {
        if (not msh.is_boundary(fc)) {
            fnum++;
            continue;
        }

        auto bnd_num = msh.boundary_info(fc).tag();
        auto dirichlet_fun = [&](const point_type& pt) {
            return pd.dirichlet_data(bnd_num, pt);
        };

        auto fb = make_scalar_monomial_basis(msh, fc, hdi.face_degree());
        Eigen::Matrix<T, Eigen::Dynamic, 1> dd =
            disk::project_function(msh, fc, fb, dirichlet_fun);
        
        dirichlet_data.segment(fnum*fb.size(), fb.size()) = dd; 

        fnum++;
    }

    auto rhs_fun = [&](const point_type& pt) {
        return pd.right_hand_side(domain_num, pt);
    };

    auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
    Eigen::Matrix<T, Eigen::Dynamic, 1> rhs = make_rhs(msh, cl, cb, rhs_fun);
    
    Eigen::Matrix<T, Eigen::Dynamic, 1> true_rhs =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(A.rows());
    true_rhs.head(cb.size()) = rhs;

    return std::tuple(A, true_rhs, dirichlet_data);
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

    std::vector<bool> is_dirichlet;
    is_dirichlet.resize( msh.faces_size() );

    for (size_t i = 0; i < msh.faces_size(); i++) {
        const auto& b = state.boundary_info[i];
        is_dirichlet[i] = false;

        if (b.type == boundary_type::dirichlet)
            is_dirichlet[i] = true;

        if (b.type == boundary_type::dirichlet_zero)
            is_dirichlet[i] = true;

        if (b.type == boundary_type::undefined)
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
        for (size_t iF = 0; iF < fcs.size(); iF++)
        {
            auto& fc = fcs[iF];
            auto lofs = cbs + fbs*iF;
            auto gofs = cbs*msh.cells_size() + fbs*offset(msh, fc);
            state.sol_full.segment(gofs, fbs) = esol.segment(lofs, fbs);
            /* XXX: adjust dirichlet bcs */
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

    auto& msh = state.msh;
    auto& assm = state.assm;
    auto& sol = state.sol;
    auto& hdi = state.hdi;

    auto cd = hdi.cell_degree();
    auto fd = hdi.face_degree();
    auto cbs = disk::scalar_basis_size(cd, Mesh::dimension);
    auto fbs = disk::scalar_basis_size(fd, Mesh::dimension-1);
    auto fullsz = cbs*msh.cells_size() + fbs*msh.faces_size();

    scalar_type L2errsq = 0.0;
    scalar_type Aerrsq = 0.0;

    size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto domain_num = msh.domain_info(cl).tag();
        auto cb = disk::make_scalar_monomial_basis(msh, cl, hdi.cell_degree());

        auto sol_fun = [&](const point_type& pt) {
            return sd.sol(domain_num, pt);
        };

        auto fcs = faces(msh, cl);
        size_t num_dofs = cbs + fbs*fcs.size();
        disk::dynamic_vector<scalar_type> num_loc_sol =
            disk::dynamic_vector<scalar_type>::Zero(num_dofs);

        num_loc_sol.head(cbs) = state.sol_full.segment(cbs*cell_i, cbs);
        num_loc_sol.tail(fbs*fcs.size()) =
            state.assm.get_element_dofs(msh, cl, state.sol);

        disk::dynamic_vector<scalar_type> ana_loc_sol =
            disk::project_function(msh, cl, hdi, sol_fun, 2);
            
        disk::dynamic_vector<scalar_type> diff =
            ana_loc_sol - num_loc_sol;

        disk::dynamic_vector<scalar_type> diff_uT = diff.head(cbs);

        disk::dynamic_matrix<scalar_type> MM = 
            disk::make_mass_matrix(msh, cl, cb);

        L2errsq += diff_uT.dot(MM*diff_uT);

        auto [GR, A] = compute_laplacian_operator(state, cl, pd);
        Aerrsq += diff.dot(A*diff);
        cell_i++;
    }

    return std::pair(sqrt(L2errsq), sqrt(Aerrsq));
}

template<typename Mesh, typename ProblemData>
void
export_to_visit(hho_poisson_solver_state<Mesh>& state,
    const ProblemData& pd)
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
        u.push_back( sol(cbs*cell_i) );
        rhs.push_back( pd.right_hand_side(1, barycenter(msh, cl)) );
        cell_i++;
    }

    disk::silo_database db;
    db.create("poisson.silo");
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
        boundary b;
        b.type = boundary_type::none;
        b.internal = false;
        b.number = 0;

        auto bi = msh.boundary_info(fc);

        if (bi.is_boundary()) {
            b.type = lua_get_boundary_type(lua, bi.tag());
            b.internal = false;
            b.number = bi.tag();
        }

        if (bi.is_internal()) {
            b.type = boundary_type::none;
            b.internal = true;
            b.number = bi.tag();
        }

        boundary_info.push_back(b);
    }
}

template<typename Mesh>
int
setup_and_run(sol::state& lua, hho_poisson_solver_state<Mesh>& state)
{
    constexpr size_t DIM = Mesh::dimension;
    lua[NODE_NAME_SIM]["dimension"] = DIM;

    collect_boundary_info(lua, state);
    lua_problem_data lpd(lua);
    lua_solution_data lsd(lua);


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

    lua["export_to_visit"] = [&]() {
        return export_to_visit(state, lpd);
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
        
        const char *mesh_filename = "../../../refactor_old_diskpp_code/meshes/2D_triangles/fvca5/mesh1_4.typ1";
        using mesh_type = disk::generic_mesh<T, 2>;
        hho_poisson_solver_state<mesh_type> state;
        disk::load_mesh_fvca5_2d<T>(mesh_filename, state.msh);
        
        //const char *mesh_filename = "../../../refactor_old_diskpp_code/meshes/2D_triangles/netgen/tri03.mesh2d";
        //using mesh_type = disk::simplicial_mesh<T, 2>;
        //hho_poisson_solver_state<mesh_type> state;
        //disk::load_mesh_netgen<T>(mesh_filename, state.msh);
        
        return setup_and_run(lua, state);
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