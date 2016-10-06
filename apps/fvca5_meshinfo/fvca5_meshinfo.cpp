#include <iostream>
#include <fstream>
#include <regex>

#include "MeshWidget.h"

#include "geometry/geometry.hpp"
#include "loaders/loader.hpp"

using RealType = double;

template<typename Mesh>
int
run_qt(int argc, char **argv, const Mesh& msh)
{
    QApplication app(argc, argv);

    MeshWidget<Mesh> w;
    w.setMesh(msh);
    w.show();
    return app.exec();
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cout << "Please specify the path of the mesh" << std::endl;
        return -1;
    }

    char *filename = argv[1];

    if (std::regex_match(filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;

        typedef disk::generic_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::fvca5_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        return run_qt(argc, argv, msh);
    }

    if (std::regex_match(filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;

        typedef disk::simplicial_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::netgen_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        return run_qt(argc, argv, msh);
    }
}
