/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2024, 2025
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <iostream>
#include "diskpp/mesh/mesh.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/loaders/loader_gmsh.hpp"
#include "diskpp/loaders/loader_utils.hpp"
#include "diskpp/output/silo.hpp"
#include "common.hpp"

int main(void)
{
    using T = double;
    disk::simplicial_mesh<T,2> finemsh;
    auto mesher = make_simple_mesher(finemsh);
    for (size_t l = 0; l < 3; l++)
        mesher.refine();
        
    partition_unit_square_mesh(finemsh, 3);

    disk::generic_mesh<T,2> coarsemsh;

    agglomerate_by_subdomain(finemsh, coarsemsh);

    disk::silo_database db;
    db.create("test_agglo_by_subdomain.silo");
    db.add_mesh(finemsh, "finemsh");
    db.add_mesh(coarsemsh, "coarsemsh");

    return 0;
}