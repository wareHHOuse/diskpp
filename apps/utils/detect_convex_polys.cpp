#include <iostream>
#include <iomanip>
#include <regex>

#include <sstream>
#include <string>

#include <unistd.h>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/output/silo.hpp"

template<typename Mesh>
int detect_convex(Mesh& msh, const char *silo_filename)
{
    std::vector<double> cvx;

    for (auto& cl : msh) {
        cvx.push_back( double(is_convex(msh, cl)) );
    }

    disk::silo_database silo_db;
    silo_db.create(silo_filename);
    silo_db.add_mesh(msh, "mesh");
    silo_db.add_variable("mesh", "isconvex", cvx, disk::zonal_variable_t);
    silo_db.close();

    return 0;
}

int main(int argc, const char *argv[])
{
    if (argc != 3)
    {
        std::cout << argv[0] << " <mesh_file> <silo_file>" << std::endl;
        return 1;
    }

    const char *mesh_filename = argv[1];
    const char *silo_filename = argv[2];

    using T = double;

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        disk::generic_mesh<T,2> msh;
        disk::load_mesh_fvca5_2d<T>(mesh_filename, msh);
        detect_convex(msh, silo_filename);
        return 0;
    }

    std::cout << "Only FVCA5 file format" << std::endl;
    return 1;
}