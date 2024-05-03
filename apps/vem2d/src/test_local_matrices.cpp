

#include <iostream>
#include <iomanip>
#include <regex>

#include <unistd.h>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/vem"
#include "diskpp/mesh/meshgen.hpp"

int main()
{
    double radius = 1.0;
    size_t num_faces = 4;
    size_t degree = 5;

    disk::generic_mesh<double,2> msh;
    disk::make_single_element_mesh(msh, radius, num_faces);

    for(const auto& cl : msh)
        std::cout<< "B matrix only in cell T: \n"
                 << disk::vem_2d::matrix_BT(msh, cl, degree) << "\n";

    return 0;
}