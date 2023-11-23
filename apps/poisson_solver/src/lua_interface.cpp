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
    lua.open_libraries(sol::lib::base, sol::lib::math, sol::lib::io, sol::lib::table, sol::lib::string);
    lua[NODE_NAME_SIM] = lua.create_table();
    lua[NODE_NAME_MESH] = lua.create_table();
    lua[NODE_NAME_MESH][MESH_FIELD_SOURCE] = MESH_SOURCE_INTERNAL;
    lua[NODE_NAME_MESH][MESH_FIELD_TYPE] = MESH_TYPE_TRIANGLES;
    lua[NODE_NAME_DOMAIN] = lua.create_table();
    lua[NODE_NAME_BOUNDARY] = lua.create_table();
    lua[NODE_NAME_SILO] = lua.create_table();
    lua[NODE_NAME_HHO] = lua.create_table();
    lua[NODE_NAME_HHO][HHO_FIELD_STABFREE] = false;
}

bool
lua_load_script(sol::state& lua, const std::string& fn)
{
    auto sf = lua.safe_script_file(fn);
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

    sol::optional<std::string> element_type_opt = lua[NODE_NAME_MESH][MESH_FIELD_TYPE];
    if (element_type_opt) {
        std::string et = element_type_opt.value();
        if (et == MESH_TYPE_TRIANGLES) ret.type = internal_mesh_type::triangles;
        else if (et == MESH_TYPE_QUADRANGLES) ret.type = internal_mesh_type::quadrangles; 
        else if (et == MESH_TYPE_HEXAGONS) ret.type = internal_mesh_type::hexagons; 
        else if (et == MESH_TYPE_TETRAHEDRA) ret.type = internal_mesh_type::tetrahedra;
        else {
            std::cout << "Invalid value '" << et << "' for " << NODE_NAME_MESH;
            std::cout << "." << MESH_FIELD_TYPE << std::endl;
            ret.type = internal_mesh_type::invalid;
        } 
    }
    else {
        std::cout << "Internal mesh: element type not specified" << std::endl;
    }

    ret.level = 0;
    auto ml = lua[NODE_NAME_MESH][MESH_FIELD_LEVEL];
    if (ml.valid()) {
        int mli = ml;
        if (mli < 0)
            mli = 0;
        ret.level = mli;
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

bool
lua_use_stabfree_hho(sol::state& lua)
{
    sol::optional<bool> use_stabfree_opt = lua[NODE_NAME_HHO][HHO_FIELD_STABFREE];
    if (use_stabfree_opt)
        return use_stabfree_opt.value();

    return false;
}

bool
lua_use_diffusion_tensor_in_stab(sol::state& lua)
{
    sol::optional<bool> dts = lua[NODE_NAME_HHO]["dt_in_stab"];
    if (dts)
        return dts.value();

    return false;
}

std::string
lua_mesh_filename(sol::state& lua)
{
    std::optional<std::string> mesh_filename_opt = lua[NODE_NAME_MESH]["filename"];
    if (mesh_filename_opt)
        return mesh_filename_opt.value();

    return "";
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
