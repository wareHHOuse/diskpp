/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *  
 * Matteo Cicuttin (C) 2020,2021
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group  
 */

#include <iostream>
#include <iomanip>
#include <regex>

#include <sstream>
#include <string>

#include <unistd.h>

#include "bases/bases.hpp"
#include "quadratures/quadratures.hpp"
#include "methods/hho"
#include "methods/implementation_hho/curl.hpp"
#include "core/loaders/loader.hpp"
#include "output/silo.hpp"
#include "solvers/solver.hpp"
#include "solvers/mumps.hpp"

int main(int argc, const char **argv)
{
    using T = double;

    if (argc != 2)
        return 1;
    
    const char *mesh_filename = argv[1];

    /* GMSH */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo2s$") ))
    {
        std::cout << "Guessed mesh format: GMSH 2D simplicials" << std::endl;
        disk::simplicial_mesh<T,2> msh;
        disk::gmsh_geometry_loader< disk::simplicial_mesh<T,2> > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        disk::silo_database silo;
        silo.create("test.silo");
        silo.add_mesh(msh, "mesh");
    }

    return 0;
}