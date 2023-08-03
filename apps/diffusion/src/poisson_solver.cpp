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

#include "sol/sol.hpp"

#include "sgr.hpp"

#define NODE_NAME_SIM           "sim"
#define SIM_MESH_FILENAME       "mesh_filename"
#define SIM_MESH_SOURCE         "mesh_source"
#define SIM_MESH_TYPE           "mesh_type"

#define NODE_NAME_SILO          "silo"

#define NODE_NAME_HHO           "hho"

#define NODE_NAME_DOMAIN        "domain"
#define NODE_NAME_BOUNDARY      "boundary"

#define MESH_SOURCE_FILE        "file"
#define MESH_SOURCE_INTERNAL    "internal"
#define MESH_TYPE_TRIANGLES     "triangles"
#define MESH_TYPE_HEXAHEDRA     "hexahedra"
#define MESH_TYPE_TETRAHEDRA    "tetrahedra"

enum class boundary_type {
    UNDEFINED,
    DIRICHLET,
    NEUMANN,
};

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

template<typename T>
class config_loader
{
    sol::state lua;

public:
    config_loader()
    {
        lua.open_libraries(sol::lib::base, sol::lib::math, sol::lib::io);
        lua[NODE_NAME_SIM] = lua.create_table();
        lua[NODE_NAME_SIM][SIM_MESH_SOURCE] = MESH_SOURCE_INTERNAL;
        lua[NODE_NAME_SIM][SIM_MESH_TYPE] = MESH_TYPE_TRIANGLES;
        lua[NODE_NAME_DOMAIN] = lua.create_table();
        lua[NODE_NAME_BOUNDARY] = lua.create_table();
        lua[NODE_NAME_SILO] = lua.create_table();
        lua[NODE_NAME_HHO] = lua.create_table();
    }

    sol::state& state() {
        return lua;
    }

    bool load(const std::string& fn)
    {
        bool success = true;
        lua.script_file(fn);

        return success;
    }

    struct mesh_params {
        mesh_source         source;
        internal_mesh_type  type;
        std::string         filename;
        size_t              level;
    };

    mesh_params mesh_parameters() const {
        mesh_params ret;

        auto fn = lua[NODE_NAME_SIM][SIM_MESH_FILENAME];
        if (fn.valid())
            ret.filename = fn;
        
        if (ret.filename == "")
            ret.source = mesh_source::internal;
        else
            ret.source = mesh_source::file;

        auto ms = lua[NODE_NAME_SIM][SIM_MESH_SOURCE];
        if (ms.valid()) {
            std::string mss = ms;
            if (mss == MESH_SOURCE_INTERNAL) ret.source = mesh_source::internal;
            else if (mss == MESH_SOURCE_FILE) ret.source = mesh_source::file; 
            else {
                std::cout << "Invalid value '" << mss << "' for " << NODE_NAME_SIM;
                std::cout << "." << SIM_MESH_SOURCE << std::endl;
                ret.source = mesh_source::invalid; 
            } 
        }

        auto mt = lua[NODE_NAME_SIM][SIM_MESH_TYPE];
        if (mt.valid()) {
            std::string mts = mt;
            if (mts == MESH_TYPE_TRIANGLES) ret.type = internal_mesh_type::triangles;
            else if (mts == MESH_TYPE_HEXAHEDRA) ret.type = internal_mesh_type::hexahedra; 
            else if (mts == MESH_TYPE_TETRAHEDRA) ret.type = internal_mesh_type::tetrahedra;
            else {
                std::cout << "Invalid value '" << mts << "' for " << NODE_NAME_SIM;
                std::cout << "." << SIM_MESH_TYPE << std::endl;
                ret.type = internal_mesh_type::invalid;
            } 
        }

        return ret;
    }

    std::string mesh_filename() const
    {
        auto mfn = lua[NODE_NAME_SIM]["mesh_filename"];
        if (not mfn.valid())
        {
            std::cout << "[CONFIG]: Mesh file name not specified ";
            std::cout << "(sim.mesh_filename)" << std::endl;
            throw std::invalid_argument("sim.mesh_filename");
        }

        return mfn;
    }
};

template<typename Mesh>
struct hho_poisson_solver_state
{
    using mesh_type = Mesh;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using scalar_type = typename Mesh::coordinate_type;

    mesh_type                                       msh;
    disk::hho_degree_info                           hdi;
    //maxwell_hho_assembler<mesh_type, scalar_type>   assm;
    disk::dynamic_vector<scalar_type>               sol;
    disk::dynamic_vector<scalar_type>               sol_full;
    disk::dynamic_vector<scalar_type>               reco;
};

template<typename Mesh>
auto
compute_element_contribution(hho_poisson_solver_state<Mesh>& state,
                             typename Mesh::cell_type& cl)
{
    using mesh_type = Mesh;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using scalar_type = typename hho_poisson_solver_state<Mesh>::scalar_type;
    using point_type = typename Mesh::point_type;

    auto& msh = state.msh;
    auto& hdi = state.hdi;
}

template<typename T>
int run(config_loader<T>& cfg)
{
    auto mps = cfg.mesh_parameters();

    if (mps.source == mesh_source::internal)
    {
        if (mps.type == internal_mesh_type::triangles)
        {
            std::cout << "tri" << std::endl;
        }

        if (mps.type == internal_mesh_type::hexahedra)
        {
            std::cout << "hex" << std::endl;
        }

        if (mps.type == internal_mesh_type::tetrahedra)
        {
            std::cout << "tet" << std::endl;
        }
    }

    if (mps.source == mesh_source::file)
    {
        std::cout << "file" << std::endl;
    }

    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <param filename>" << std::endl;
        return 1;
    }

    using T = double;

    config_loader<T> cfg;
    sol::state &lua = cfg.state();
    lua["run"] = [&](){ return run(cfg); };

    cfg.load(argv[1]);

    return 0;
}