/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2021
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group
 */

#include <iostream>

#include "mesh/mesh.hpp"
#include "core/loaders/loader.hpp"
#include "core/output/silo.hpp"

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cout << "Please specify GMSH .geo file" << std::endl;
        return 1;
    }

    using T = double;
    disk::generic_mesh<T,3>  msh;
    disk::gmsh_geometry_loader< disk::generic_mesh<T,3> > loader;

    //loader.verbose(true);
    loader.read_mesh(argv[1]);
    loader.populate_mesh(msh);

    std::vector<double> dn;
    for (auto& cl : cells(msh))
    {
        auto di = msh.domain_info(cl);
        dn.push_back( di.tag() );
    }

    std::map<size_t, size_t> bnd_elemcount;
    for (auto& fc : faces(msh))
    {
        auto bi = msh.boundary_info(fc);
        bnd_elemcount[bi.tag()]++;
    }

    for (auto& [bnd, count] : bnd_elemcount)
        std::cout << bnd << ": " << count << " elements" << std::endl;

    disk::silo_zonal_variable<T> dns("subdomain_number", dn);

    disk::silo_database silo_db;
    silo_db.create("gmsh.silo");
    silo_db.add_mesh(msh, "mesh");
    silo_db.add_variable("mesh", dns);

    return 0;
}