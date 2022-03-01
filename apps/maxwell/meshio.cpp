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

#include "compinfo.h"

#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include "paramloader.hpp"


int main(int argc, const char *argv[])
{
    using T = double;
    disk::simplicial_mesh<T,3> msh;
    disk::gmsh_geometry_loader< disk::simplicial_mesh<T,3> > loader;
        
    loader.read_mesh(argv[1]);
    loader.populate_mesh(msh);

    disk::silo_database silo_db;
    silo_db.create("mesh.silo");
    silo_db.add_mesh(msh, "mesh");
    silo_db.close();

    return 0;
}

