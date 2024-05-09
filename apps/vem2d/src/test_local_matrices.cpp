

#include <iostream>
#include <iomanip>
#include <regex>

#include <unistd.h>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/vem"
#include "diskpp/mesh/meshgen.hpp"

int main(int argc, char **argv)
{
    using T = double;
    size_t degree = 2;

    std::string mesh_filename;
    if (argc > 1)
        mesh_filename = argv[1];

    if ( mesh_filename.length() > 0 )
    {
        /* Single element CSV 2D */
        if (std::regex_match(mesh_filename, std::regex(".*\\.csv$") ))
        {
            std::cout << "Guessed mesh format: CSV 2D" << std::endl;
            disk::generic_mesh<T,2> msh;
            disk::load_single_element_csv(msh, mesh_filename);
    
            for(const auto& cl : msh)
                std::cout<< "B matrix only in cell T: \n"
                    << disk::vem_2d::matrix_BT(msh, cl, degree) << "\n";

            for(const auto& cl : msh)
                std::cout<< "B matrix only in cell F: \n"
                    << disk::vem_2d::matrix_BF(msh, cl, degree) << "\n";

            return 0;
        }
    }
    else {
        double radius = 1.0;
        size_t num_faces = 4;

        disk::generic_mesh<double,2> msh;
        disk::make_single_element_mesh(msh, radius, num_faces);

        for(const auto& cl : msh)
            std::cout<< "B matrix only in cell T: \n"
                    << disk::vem_2d::matrix_BT(msh, cl, degree) << "\n";

        for(const auto& cl : msh)
            std::cout<< "B matrix only in cell F: \n"
                    << disk::vem_2d::matrix_BF(msh, cl, 2) << "\n";
    }

    return 0;
}