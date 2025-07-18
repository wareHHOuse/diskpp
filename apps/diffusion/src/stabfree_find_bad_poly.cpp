#include <iostream>
#include <complex>

#include "stabfree_explorer.hpp"

int main(void)
{
    using T = double;
    using mesh_type = disk::generic_mesh<T,2>;
    using point_type = mesh_type::point_type;

    int numCPU = sysconf(_SC_NPROCESSORS_ONLN);
    std::cout << numCPU << std::endl;

    for (size_t i = 6; i < 11; i++) {
        std::cout << " **** Faces = " << i << " ****" << std::endl;
        mesh_type msh_gen;
        double scale = 1.0;
        disk::make_single_element_mesh(msh_gen, scale, i);

        disk::silo_database silo;
        std::string fn = "badpolys_" + std::to_string(i) + ".silo";
        silo.create(fn);
        silo.add_mesh(msh_gen, "mesh");

        minimize_step(msh_gen);

        explore(msh_gen, 0);

        break;
    }

    

    return 0;
}