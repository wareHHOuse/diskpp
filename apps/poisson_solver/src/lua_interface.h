/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#pragma once

#include "diskpp/mesh/point.hpp"

#define SOL_ALL_SAFETIES_ON 1
#include "sol/sol.hpp"

#define NODE_NAME_SIM           "sim"


#define NODE_NAME_SILO          "silo"


#define NODE_NAME_HHO           "hho"
/**/
#define HHO_FIELD_VARIANT       "variant"
#define HHO_FIELD_ORDER         "order"
#define HHO_MOL_VARIANT_NAME    "mixed_order_low"
#define HHO_MOH_VARIANT_NAME    "mixed_order_high"
#define HHO_EO_VARIANT_NAME     "equal_order"
/**/
#define HHO_FIELD_STABFREE      "use_stabfree"

#define NODE_NAME_MESH          "mesh"
/**/
#define MESH_FIELD_SOURCE       "source"
#define MESH_SOURCE_FILE        "file"
#define MESH_SOURCE_INTERNAL    "internal"
/**/
#define MESH_FIELD_TYPE         "type"
#define MESH_TYPE_TRIANGLES     "triangles"
#define MESH_TYPE_HEXAGONS      "hexagons"
#define MESH_TYPE_TETRAHEDRA    "tetrahedra"
/**/
#define MESH_FIELD_LEVEL        "refinement_level"
#define MESH_FIELD_FILENAME     "filename"


#define NODE_NAME_DOMAIN        "domain"
#define NODE_NAME_BOUNDARY      "boundary"

enum class mesh_source {
    internal,
    file,
    invalid,
};

enum class internal_mesh_type {
    triangles,
    hexagons,
    tetrahedra,
    invalid,
};

enum class hho_variant {
    mixed_order_low,
    equal_order,
    mixed_order_high,
};

enum class boundary_type {
    undefined,
    none,
    internal_interface,
    dirichlet,
    dirichlet_zero,
    neumann,
    neumann_zero,
};

struct boundary
{
    boundary_type       type;
    bool                internal;
    size_t              number;
};

struct mesh_parameters {
    mesh_source         source;
    internal_mesh_type  type;
    std::string         filename;
    size_t              level;
};


void lua_init_environment(sol::state&);
bool lua_load_script(sol::state&, const std::string&);
mesh_parameters lua_get_mesh_parameters(sol::state&);
int lua_get_hho_order(sol::state&);
hho_variant lua_get_hho_variant(sol::state&);
bool lua_use_stabfree_hho(sol::state&);
boundary_type lua_get_boundary_type(sol::state&, size_t);
int lua_call_user_code(sol::state&);

template<typename T>
T
lua_eval_rhs(sol::state& lua, size_t domain_num, const disk::point<T,2>& pt)
{
    return lua["right_hand_side"](domain_num, pt.x(), pt.y());
}

template<typename T>
T
lua_eval_rhs(sol::state& lua, size_t domain_num, const disk::point<T,3>& pt)
{
    return lua["right_hand_side"](domain_num, pt.x(), pt.y(), pt.z());
}



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




struct lua_problem_data
{
    sol::state& lua;

    lua_problem_data(sol::state& plua)
        : lua(plua)
    {}

    template<typename T>
    T dirichlet_data(size_t boundary_num, const disk::point<T,2>& pt) const
    {
        return lua["dirichlet_data"](boundary_num, pt.x(), pt.y());
    }

    template<typename T>
    T dirichlet_data(size_t boundary_num, const disk::point<T,3>& pt) const
    {
        return lua["dirichlet_data"](boundary_num, pt.x(), pt.y(), pt.z());
    }

    template<typename T>
    T right_hand_side(size_t domain_num, const disk::point<T,2>& pt) const
    {
        return lua["right_hand_side"](domain_num, pt.x(), pt.y());
    }

    template<typename T>
    T right_hand_side(size_t domain_num, const disk::point<T,3>& pt) const
    {
        return lua["right_hand_side"](domain_num, pt.x(), pt.y(), pt.z());
    }

    template<typename T>
    diffusion_parameter<T,2>
    diffusion_coefficient(size_t domain_num, const disk::point<T,2>& pt) const
    {
        return lua["diffusion_coefficient"](domain_num, pt.x(), pt.y());
    }
};

struct lua_solution_data
{
    sol::state& lua;

    lua_solution_data(sol::state& plua)
        : lua(plua)
    {}

    template<typename T>
    T sol(size_t domain_num, const disk::point<T,2>& pt) const
    {
        return lua["solution"](domain_num, pt.x(), pt.y());
    }

    template<typename T>
    T sol(size_t domain_num, const disk::point<T,3>& pt) const
    {
        return lua["solution"](domain_num, pt.x(), pt.y(), pt.z());
    }

    template<typename T>
    Eigen::Matrix<T,2,1>
    grad(size_t domain_num, const disk::point<T,2>& pt) const
    {
        Eigen::Matrix<T,2,1> ret;
        return lua["grad_solution"](domain_num, pt.x(), pt.y());
    }

    template<typename T>
    Eigen::Matrix<T,3,1>
    grad(size_t domain_num, const disk::point<T,3>& pt) const
    {
        Eigen::Matrix<T,3,1> ret;
        return lua["grad_solution"](domain_num, pt.x(), pt.y(), pt.z());
    }
};