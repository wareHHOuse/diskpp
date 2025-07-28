#include <iostream>
#include "diskpp/loaders/loader.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/output/silo.hpp"
#include "gmsh.h"




int main(void)
{
    using T = double;
    using coarse_mesh_type = disk::generic_mesh<T,2>;
    using fine_mesh_type = disk::simplicial_mesh<T, 2>;

    coarse_mesh_type coarse_msh;
    auto mesher = disk::make_fvca5_hex_mesher(coarse_msh);
    mesher.make_level(1);

    fine_mesh_type fine_msh;
    submesh_via_gmsh(coarse_msh, fine_msh, 0.1);

    disk::silo_database db;
    db.create("out.silo");
    db.add_mesh(coarse_msh, "coarse_mesh");
    db.add_mesh(fine_msh, "fine_mesh");

    std::vector<double> dn;
    for (const auto& cl : fine_msh) {
        auto di = fine_msh.domain_info(cl);
        dn.push_back( di.tag() );
    }

    db.add_variable("fine_mesh", "parent", dn, disk::zonal_variable_t);

    return 0;
}