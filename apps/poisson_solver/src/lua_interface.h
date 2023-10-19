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
#include "poisson_solver.h"

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
#define MESH_TYPE_QUADRANGLES   "quadrangles"
#define MESH_TYPE_HEXAGONS      "hexagons"
#define MESH_TYPE_TETRAHEDRA    "tetrahedra"
/**/
#define MESH_FIELD_LEVEL        "refinement_level"
#define MESH_FIELD_FILENAME     "filename"


#define NODE_NAME_DOMAIN        "domain"
#define NODE_NAME_BOUNDARY      "boundary"




void lua_init_environment(sol::state&);
bool lua_load_script(sol::state&, const std::string&);
mesh_parameters lua_get_mesh_parameters(sol::state&);
int lua_get_hho_order(sol::state&);
hho_variant lua_get_hho_variant(sol::state&);
bool lua_use_stabfree_hho(sol::state&);
bool lua_use_diffusion_tensor_in_stab(sol::state&);
std::string lua_mesh_filename(sol::state&);
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







template<typename T>
boundary_condition_descriptor<T>
lua_get_boundary_condition(sol::state& lua, size_t bnd)
{
    boundary_condition_descriptor<T> ret;
    ret.type = boundary_type::undefined;
    ret.number = bnd;

    sol::optional<std::string> bnd_type_opt = lua[NODE_NAME_BOUNDARY][bnd]["type"];
    if (bnd_type_opt)
    {
        ret.type = boundary_type::invalid;
        std::string bnd_type_name = bnd_type_opt.value();
        if (bnd_type_name == "dirichlet")   ret.type = boundary_type::dirichlet;
        if (bnd_type_name == "neumann")     ret.type = boundary_type::neumann;
        if (bnd_type_name == "robin")       ret.type = boundary_type::robin;
        if (bnd_type_name == "jump")        ret.type = boundary_type::jump;

        if (ret.type == boundary_type::invalid) {
            std::cout << "Invalid type '" << bnd_type_name;
            std::cout << "' set on boundary " << bnd << std::endl;
        }
    }

    if (ret.type == boundary_type::robin) {
        sol::optional<T> bnd_value_opt = lua[NODE_NAME_BOUNDARY][bnd]["value"];
        if (bnd_value_opt)
            ret.value = bnd_value_opt.value();
        else
            std::cout << "Warning: value not set on Robin boundary " << bnd << std::endl;
    }

    return ret;
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
    {
        auto dd_fun = lua["dirichlet_data"];
        if (not dd_fun.valid()) {
            std::cout << "dirichlet_data() undefined, defaulting to zero" << std::endl;
        }

        auto nd_fun = lua["neumann_data"];
        if (not nd_fun.valid()) {
            std::cout << "neumann_data() undefined, defaulting to zero" << std::endl;
        }

        auto rhs_fun = lua["right_hand_side"];
        if (not rhs_fun.valid()) {
            std::cout << "right_hand_side() undefined, defaulting to zero" << std::endl;
        }
    }

    template<typename T, size_t N>
    T dirichlet_data(size_t boundary_num, const disk::point<T,N>& pt) const
    {
        auto dd_fun = lua["dirichlet_data"];
        if (not dd_fun.valid())
            return 0.0;

        return dd_fun(boundary_num, pt);
    }

    template<typename T, size_t N>
    T neumann_data(size_t boundary_num, const disk::point<T,N>& pt,
        const disk::point<T,N>& normal) const
    {
        auto nd_fun = lua["neumann_data"];
        if (not nd_fun.valid())
            return 0.0;

        return nd_fun(boundary_num, pt, normal);
    }

    template<typename T, size_t N>
    T right_hand_side(size_t domain_num, const disk::point<T,N>& pt) const
    {
        auto rhs_fun = lua["right_hand_side"];
        if (not rhs_fun.valid())
            return 0.0;

        return rhs_fun(domain_num, pt);
    }

    template<typename T>
    diffusion_parameter<T,2>
    diffusion_coefficient(size_t domain_num, const disk::point<T,2>& pt) const
    {
        auto dc = lua["diffusion_coefficient"](domain_num, pt.x(), pt.y());

        sol::optional<double> dc_dbl = dc;
        if (dc_dbl) {
            diffusion_parameter<T,2> ret;
            ret.entry(0, 0, dc_dbl.value());
            ret.entry(1, 1, dc_dbl.value());
            return ret;
        }

        sol::optional<diffusion_parameter<T,2>> dc_tens = dc;
        if (dc_tens) {
            return dc_tens.value();
        }

        throw std::invalid_argument("Can't convert diffusion coefficient to a valid type.");
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
        T gx, gy;
        sol::tie(gx, gy) = lua["solution_gradient"](domain_num, pt.x(), pt.y());
        ret(0) = gx;
        ret(1) = gy;
        return ret;
    }

    template<typename T>
    Eigen::Matrix<T,3,1>
    grad(size_t domain_num, const disk::point<T,3>& pt) const
    {
        Eigen::Matrix<T,3,1> ret;
        T gx, gy, gz;
        sol::tie(gx, gy, gz) = lua["solution_gradient"](domain_num, pt.x(), pt.y(), pt.z());
        ret(0) = gx;
        ret(1) = gy;
        ret(2) = gz;
        return ret;
    }
};