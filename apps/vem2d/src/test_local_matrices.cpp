

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
    disk::generic_mesh<T,2> msh;

    std::string mesh_filename;
    if (argc > 1 && mesh_filename.length() > 0 )
    {
        mesh_filename = argv[1];

        /* Single element CSV 2D */
        if (std::regex_match(mesh_filename, std::regex(".*\\.csv$") ))
        {
            std::cout << "Guessed mesh format: CSV 2D" << std::endl;
            disk::load_single_element_csv(msh, mesh_filename);
        }    
    }
    else 
    {    
        double radius = 1.0;
        size_t num_faces = 4;

        disk::make_single_element_mesh(msh, radius, num_faces);
    }

    std::cout << "LEX ordering" << std::endl;
    for (auto& cl : msh) {
        auto fcs = faces(msh, cl);
        for (auto& fc : fcs)
            std::cout << fc << std::endl;
    }

    std::cout << "CCW ordering" << std::endl;
    for (auto& cl : msh) {
        auto fcs_ccw = faces_ccw(msh, cl);
        for (auto& [fc, flip] : fcs_ccw)
            std::cout << fc << " " << flip << std::endl;
    }

    for(const auto& cl : msh)
        std::cout<< "G matrix: \n"
            << disk::vem_2d::matrix_BT(msh, cl, degree) << "\n";

    for(const auto& cl : msh)
        std::cout<< "B matrix: \n"
            << disk::vem_2d::matrix_B(msh, cl, 2) << "\n";
    
    for(const auto& cl : msh)
        std::cout<< "D matrix: \n"
            << disk::vem_2d::matrix_D(msh, cl, degree) << "\n";
    
    return 0;
}