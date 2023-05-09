/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <iostream>
#include <fstream>

#include "diskpp/bases/bases_new.hpp"
#include "diskpp/bases/bases_operations.hpp"
#include "diskpp/quadratures/quadratures.hpp"

template<typename Mesh>
void
eval_conditioning(const Mesh& msh)
{
    for (auto& cl : msh)
    {
        
    }
}

int main(int argc, char **argv)
{

    if (argc != 2)
    {
        std::cout << argv[0] << " <mesh file>" << std::endl;
        return 1;
    }

    auto mesh_filename = argv[1];

    return disk::dispatch_all_meshes(mesh_filename,
              [](auto ...args) { eval_conditioning(args...); }
            );
}
