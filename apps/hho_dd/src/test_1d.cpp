/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2026
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <cstddef>
#include <iostream>
#include <regex>
#include <set>
#include <string>
#include <filesystem>

#include "diskpp/common/util.h"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/loaders/loader_gmsh.hpp"
#include "diskpp/loaders/loader_utils.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/output/silo.hpp"
#include "diskpp/methods/hho_slapl.hpp"
#include "diskpp/methods/hho"
#include "diskpp/methods/dg"
#include "rasdd.hpp"
#include "common.hpp"
#include "solvers.hpp"
#include "diskpp_git_revision.h"


int main(void)
{
    using T = double;

    //disk::generic_mesh<T,1> msh;
    //disk::uniform_mesh_loader<T,1> loader(0,1,100);
    //loader.populate_mesh(msh);

    disk::generic_mesh<T, 2> msh;
    auto mesher = disk::make_fvca5_hex_mesher(msh);
    mesher.make_level(2);

    disk::silo_database db;
    db.create("test_1d.silo");
    db.add_mesh(msh, "srcmesh");
    db.add_mesh(msh, "dstmesh");
    hho_diffusion_solver(msh, 7, db);
    dg_diffusion_solver(msh, 8, 10.0, db);


    return 0;
}