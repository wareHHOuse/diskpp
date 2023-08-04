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

#include "lua_interface.h"

#include "sgr.hpp"




template<typename T, size_t DIM>
class diffusion_parameter
{
    using tens_type = Eigen::Matrix<T, DIM, DIM>;
    tens_type diff_tens;

public:
    diffusion_parameter()
    {
        diff_tens = tens_type::Identity();
    }

    diffusion_parameter(T value)
    {
        diff_tens = tens_type::Identity() * value;
    }
    
    void entry(int row, int col, T data) {
        bool row_ok = row >= 0 and row <= DIM;
        bool col_ok = col >= 0 and col <= DIM;
        if (row_ok and col_ok)
            diff_tens(row, col) = data;
    }

    T entry(int row, int col) const {
        bool row_ok = row >= 0 and row <= DIM;
        bool col_ok = col >= 0 and col <= DIM;
        if (not row_ok or not col_ok)
            throw std::invalid_argument("tensor access out of bounds");
        return diff_tens(row, col);
    }

    tens_type tensor(void) const {
        return diff_tens;
    }
};

template<typename T, size_t DIM>
std::ostream&
operator<<(std::ostream& os, const diffusion_parameter<T, DIM>& dp)
{
    os << dp.tensor();
    return os;
}

template<typename T>
using diffusion_parameter_2D = diffusion_parameter<T,2>;

template<typename T>
using diffusion_parameter_3D = diffusion_parameter<T,3>;

template<typename Mesh>
struct hho_poisson_solver_state
{
    using mesh_type = Mesh;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using scalar_type = typename Mesh::coordinate_type;
    using assembler_type = disk::diffusion_condensed_assembler<Mesh>;
    using vector_type = disk::dynamic_vector<scalar_type>;

    mesh_type               msh;
    assembler_type          assm;
    vector_type             sol;
    vector_type             sol_full;

    std::vector<boundary>   boundary_info;

    disk::hho_degree_info   hdi;
    hho_variant             variant;
    bool                    use_stabfree;
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
compute_element_contribution_stabfree(hho_poisson_solver_state<Mesh>& state,
    typename Mesh::cell_type& cl, const ProblemData& pd)
{
    using mesh_type = Mesh;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using point_type = typename Mesh::point_type;

    auto& msh = state.msh;
    auto& hdi = state.hdi;

    adjust_stabfree_recdeg(msh, cl, hdi);
}

template<typename Mesh, typename ProblemData>
auto
compute_element_contribution_standard(hho_poisson_solver_state<Mesh>& state,
    typename Mesh::cell_type& cl, const ProblemData& pd)
{
    using mesh_type = Mesh;
    using T = typename mesh_type::coordinate_type;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using point_type = typename Mesh::point_type;

    auto& msh = state.msh;
    auto& hdi = state.hdi;

    auto bar = barycenter(msh, cl);
    auto domain_num = msh.domain_info(cl).tag();
    //auto diff_tens = 

    bool is_mixed_high = (hdi.cell_degree() > hdi.face_degree());

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> GR;
    
    Eigen::Matrix<T, Mesh::dimension, Mesh::dimension> diff_tens =
        Eigen::Matrix<T, Mesh::dimension, Mesh::dimension>::Identity();

    auto oper = make_scalar_hho_laplacian(msh, cl, hdi, diff_tens);
    A = oper.second;
    GR = oper.first;

    if (is_mixed_high)
        A = A + make_scalar_hdg_stabilization(msh, cl, hdi);
    else
        A = A + make_scalar_hho_stabilization(msh, cl, GR, hdi);
    
    auto rhs_data = [&](const point_type& pt) {
        return pd.right_hand_side(domain_num, pt);
    };

    auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
    Eigen::Matrix<T, Eigen::Dynamic, 1> rhs = make_rhs(msh, cl, cb, rhs_data);
    //auto sc = make_scalar_static_condensation(msh, cl, hdi, A, rhs);
    //assm.assemble(msh, cl, sc.first, sc.second, sol_fun);

    return std::pair(GR, A);
    //return std::tuple(lhs, rhs, dirichlet_data);
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

    for (auto& cl : msh)
    {}

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

    for (auto& cl : msh)
    {
        auto [GR, A] = compute_element_contribution_standard(state, cl, pd);
    }

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

    std::cout << "Assembler initialized: " << assm.LHS.rows() << " DOFs" << std::endl;
    
    if (state.use_stabfree)
        assemble_stabfree(state, pd);
    else
        assemble_standard(state, pd);
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
    //lua[NODE_NAME_SIM]["dimension"] = Mesh::dimension;

    collect_boundary_info(lua, state);
    lua_problem_data lpd(lua);

    /* Ok, here the assumption is that state and lua are alive during the
     * execution of setup_and_run, so all the following lambdas are safe.
     * For extra caution, we unregister all the functions after calling
     * the user code, so they don't get called after a call to run() on
     * the Lua side. */

    /* Bind functions */
    lua["assemble"] = [&]() {
        return assemble(state, lpd);
    };

    /* Call user code */
    int err = lua_call_user_code(lua);
    
    /* Unbind everything to avoid unsafe calls */
    lua["assemble"] = nullptr;

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
    }

    if (order < 1 and variant == hho_variant::mixed_order_low) {
        std::cout << "Invalid HHO order " << order << ", reverting ";
        std::cout << "to order 1" << std::endl;
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

    state.use_stabfree = false;

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
            return setup_and_run(lua, state);
        }

        if (mps.type == internal_mesh_type::hexahedra)
        {
            using mesh_type = disk::generic_mesh<T,2>;
            hho_poisson_solver_state<mesh_type> state;
            init_solver_state(lua, state);
            auto mesher = disk::make_fvca5_hex_mesher(state.msh);
            //mesher.make_level(5);
            return setup_and_run(lua, state);
        }

        if (mps.type == internal_mesh_type::tetrahedra)
        {
            using mesh_type = disk::simplicial_mesh<T,3>;
            hho_poisson_solver_state<mesh_type> state;
            init_solver_state(lua, state);
            auto mesher = disk::make_simple_mesher(state.msh);
            return setup_and_run(lua, state);
        }
    }

    if (mps.source == mesh_source::file)
    {
        std::cout << "file" << std::endl;
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

    lua.new_usertype<dp2d>("diffusion_parameter_2D",
		sol::constructors<dp2d(), dp2d(T)>(),
        "entry",
        sol::overload(
            sol::resolve<void(int, int, T)>(&dp2d::entry),
            sol::resolve<T(int, int) const>(&dp2d::entry)
        )
    );

    lua.new_usertype<dp3d>("diffusion_parameter_3D",
		sol::constructors<dp3d(), dp3d(T)>(),
        "entry",
        sol::overload(
            sol::resolve<void(int, int, T)>(&dp3d::entry),
            sol::resolve<T(int, int) const>(&dp3d::entry)
        )
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