#include <iostream>
#include <iomanip>
#include <regex>

#include <sstream>
#include <string>

#include <unistd.h>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/output/silo.hpp"

template<typename Mesh>
int export_mesh_to_silo(Mesh& msh, const char *silo_filename)
{
    disk::silo_database silo_db;
    silo_db.create(silo_filename);
    silo_db.add_mesh(msh, "mesh");
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
        auto msh = disk::load_fvca5_2d_mesh<T>(mesh_filename);
        export_mesh_to_silo(msh, silo_filename);
        return 0;
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<T>(mesh_filename);
        export_mesh_to_silo(msh, silo_filename);
        return 0;
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        auto msh = disk::load_cartesian_2d_mesh<T>(mesh_filename);
        export_mesh_to_silo(msh, silo_filename);
        return 0;
    }


    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        auto msh = disk::load_netgen_3d_mesh<T>(mesh_filename);
        export_mesh_to_silo(msh, silo_filename);
        return 0;
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        auto msh = disk::load_cartesian_3d_mesh<T>(mesh_filename);
        export_mesh_to_silo(msh, silo_filename);
        return 0;
    }

    /* FVCA6 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
        disk::generic_mesh<T,3> msh;
        disk::load_mesh_fvca6_3d<T>(mesh_filename, msh);
        export_mesh_to_silo(msh, silo_filename);
        return 0;
    }
    
    return 0;
}

