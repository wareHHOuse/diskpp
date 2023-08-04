/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include "lua_interface.h"



void
lua_init_environment(sol::state& lua)
{
    lua.open_libraries(sol::lib::base, sol::lib::math, sol::lib::io);
    lua[NODE_NAME_SIM] = lua.create_table();
    lua[NODE_NAME_MESH] = lua.create_table();
    lua[NODE_NAME_MESH][MESH_FIELD_SOURCE] = MESH_SOURCE_INTERNAL;
    lua[NODE_NAME_MESH][MESH_FIELD_TYPE] = MESH_TYPE_TRIANGLES;
    lua[NODE_NAME_DOMAIN] = lua.create_table();
    lua[NODE_NAME_BOUNDARY] = lua.create_table();
    lua[NODE_NAME_SILO] = lua.create_table();
    lua[NODE_NAME_HHO] = lua.create_table();
}

bool
lua_load_script(sol::state& lua, const std::string& fn)
{
    auto sf = lua.script_file(fn);
    return sf.valid();
}


mesh_parameters
lua_get_mesh_parameters(sol::state& lua)
{
    mesh_parameters ret;

    auto fn = lua[NODE_NAME_MESH][MESH_FIELD_FILENAME];
    if (fn.valid())
        ret.filename = fn;
    
    if (ret.filename == "")
        ret.source = mesh_source::internal;
    else
        ret.source = mesh_source::file;

    auto ms = lua[NODE_NAME_MESH][MESH_FIELD_SOURCE];
    if (ms.valid()) {
        std::string mss = ms;
        if (mss == MESH_SOURCE_INTERNAL) ret.source = mesh_source::internal;
        else if (mss == MESH_SOURCE_FILE) ret.source = mesh_source::file; 
        else {
            std::cout << "Invalid value '" << mss << "' for " << NODE_NAME_MESH;
            std::cout << "." << MESH_FIELD_SOURCE << std::endl;
            ret.source = mesh_source::invalid; 
        } 
    }

    auto mt = lua[NODE_NAME_MESH][MESH_FIELD_TYPE];
    if (mt.valid()) {
        std::string mts = mt;
        if (mts == MESH_TYPE_TRIANGLES) ret.type = internal_mesh_type::triangles;
        else if (mts == MESH_TYPE_HEXAHEDRA) ret.type = internal_mesh_type::hexahedra; 
        else if (mts == MESH_TYPE_TETRAHEDRA) ret.type = internal_mesh_type::tetrahedra;
        else {
            std::cout << "Invalid value '" << mts << "' for " << NODE_NAME_MESH;
            std::cout << "." << MESH_FIELD_TYPE << std::endl;
            ret.type = internal_mesh_type::invalid;
        } 
    }

    return ret;
}

int
lua_get_hho_order(sol::state& lua)
{
    auto order = lua[NODE_NAME_HHO][HHO_FIELD_ORDER];
    if (not order.valid())
        return 1;

    return order;
}

hho_variant
lua_get_hho_variant(sol::state& lua)
{
    enum hho_variant ret = hho_variant::equal_order;
    auto varname = lua[NODE_NAME_HHO][HHO_FIELD_VARIANT];
    if (varname.valid()) {
        std::string vn = varname;
        if (vn == HHO_EO_VARIANT_NAME) ret = hho_variant::equal_order;
        else if (vn == HHO_MOL_VARIANT_NAME) ret = hho_variant::mixed_order_low;
        else if (vn == HHO_MOH_VARIANT_NAME) ret = hho_variant::mixed_order_high;
        else {
            std::cout << "Invalid value '" << vn << "' for " << NODE_NAME_HHO;
            std::cout << "." << HHO_FIELD_VARIANT << std::endl;
        }
    }
    return ret;
}

boundary_type
lua_get_boundary_type(sol::state& lua, size_t bnd)
{
    auto bndnode = lua[NODE_NAME_BOUNDARY][bnd];
    if (not bndnode.valid())
        return boundary_type::undefined;
    
    std::string k = bndnode;
    if (k == "dirichlet") return boundary_type::dirichlet;
    else if (k == "dirichlet_zero") return boundary_type::dirichlet_zero;
    else if (k == "neumann_zero") return boundary_type::neumann_zero;
    else {
        std::cout << "Boundary " << bnd << ": invalid kind '" << k;
        std::cout << "'" << std::endl;
        return boundary_type::undefined;
    }
}


int
lua_call_user_code(sol::state& lua)
{
    auto solproc = lua["solution_process"];
    if (not solproc.valid()) {
        std::cout << "The 'solution_process()' function is undefined ";
        std::cout << "in the user code. I have nothing to do." << std::endl;
        return 1;
    }

    return solproc();
}
