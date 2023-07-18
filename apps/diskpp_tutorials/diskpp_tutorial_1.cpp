/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

/*
 *       /\        Matteo Cicuttin (C) 2016-2020
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */

/* DiSk++ tutorial: create a mesh on the unit square/unit cube and iterate
 * on its elements. Print some properties of the elements.
 */

#include <iostream>
#include <unistd.h>

// For the mesh data structure
#include "diskpp/mesh/mesh.hpp"
// For the unit square/cube mesh generators
#include "diskpp/mesh/meshgen.hpp"

enum class element_type {
    triangle,
    hexagon,
    tetrahedron,
};

/* This function is generic and will work for any kind of mesh */
template<typename Mesh>
void
print_element_properties(const Mesh& msh)
{
    for (auto& cl : msh)
    {
        std::cout << cl << std::endl;
        std::cout << "  Diameter: " << diameter(msh, cl) << std::endl;
        std::cout << "  Measure: " << measure(msh, cl) << std::endl;
        std::cout << "  Barycenter: " << barycenter(msh, cl) << std::endl;
        auto fcs = faces(msh, cl);
        for (auto& fc : fcs)
        {
            std::cout << "  " <<fc << std::endl;
            std::cout << "    Diameter: " << diameter(msh, fc) << std::endl;
            std::cout << "    Measure: " << measure(msh, fc) << std::endl;
            std::cout << "    Barycenter: " << barycenter(msh, fc) << std::endl;
        }
    }
}

/* If your function is appropriate only for 1D, 2D or 3D meshes, you
 * can constrain the mesh type via the appropriate concept */
template<disk::mesh_2D Mesh>
void
print_space_dimension(const Mesh& msh)
{
    std::cout << "This is a 2D mesh" << std::endl;
}

template<disk::mesh_3D Mesh>
void
print_space_dimension(const Mesh& msh)
{
    std::cout << "This is a 3D mesh" << std::endl;
}

int main(int argc, char **argv)
{
    int num_refs = 0;
    element_type elem_type = element_type::triangle;

    /* Parse the command line arguments */
    int ch;
    while ( (ch = getopt(argc, argv, "r:thT")) != -1 )
    {
        switch(ch)
        {
            case 'r':
                num_refs = std::stoi(optarg);
                if (num_refs < 0)
                    num_refs = 0;
                break;

            case 't':
                elem_type = element_type::triangle;
                break;

            case 'h':
                elem_type = element_type::hexagon;
                break;

            case 'T':
                elem_type = element_type::tetrahedron;
                break;

            case '?':
            default:
                std::cout << "Invalid option" << std::endl;
                return 1;
        }
    }

    /* Coordinate type is double */
    using T = double;
    
    switch (elem_type)
    {
        case element_type::triangle: {
            /* The mesh type is a 2D simplicial mesh */
            using mesh_type = disk::simplicial_mesh<T,2>;
            /* Declare the mesh object... */
            mesh_type msh;
            /* ...and instantiate a mesh generator over it */
            auto mesher = disk::make_simple_mesher(msh);
            /* Refine the mesh */
            for (auto nr = 0; nr < num_refs; nr++)
                mesher.refine();
            /* Call the code. */
            print_space_dimension(msh);
            print_element_properties(msh);
        } break;

        case element_type::hexagon: {
            /* The mesh type is a 2D polygonal mesh */
            using mesh_type = disk::generic_mesh<T,2>;
            /* Declare the mesh object... */
            mesh_type msh;
            /* ...and instantiate a mesh generator over it */
            auto mesher = disk::make_fvca5_hex_mesher(msh);
            /* Make the specified level of the mesh */
            mesher.make_level(num_refs);
            /* Call the code. */
            print_space_dimension(msh);
            print_element_properties(msh);
        } break;

        case element_type::tetrahedron: {
            /* The mesh type is a 3D simplicial mesh */
            using mesh_type = disk::simplicial_mesh<T,3>;
            /* Declare the mesh object... */
            mesh_type msh;
            /* ...and instantiate a mesh generator over it */
            auto mesher = disk::make_simple_mesher(msh);
            /* Refine the mesh */
            for (auto nr = 0; nr < num_refs; nr++)
                mesher.refine();
            /* Call the code. */
            print_space_dimension(msh);
            print_element_properties(msh);
        } break;
    }

    return 0;
}


