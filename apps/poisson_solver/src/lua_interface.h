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
#define MESH_TYPE_HEXAHEDRA     "hexahedra"
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
    hexahedra,
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
};