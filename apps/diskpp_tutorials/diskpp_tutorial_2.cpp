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

/* DiSk++ tutorial: create a mesh on the unit square/unit cube, iterate on
 * its elements and compute the values of some function. Output it to VisIt.
 */



#include <iostream>
#include <unistd.h>

// For the mesh data structure
#include "diskpp/mesh/mesh.hpp"
// For the unit square/cube mesh generators
#include "diskpp/mesh/meshgen.hpp"
// For the output to SILO/VisIt
#include "diskpp/output/silo.hpp"

enum class element_type {
    triangle,
    hexagon,
    tetrahedron,
};

template<typename T>
auto fun(const disk::point<T,2>& pt)
{
    auto sx = std::sin(M_PI*pt.x());
    auto sy = std::sin(M_PI*pt.y());
    return sx*sy;
}

template<typename T>
auto fun(const disk::point<T,3>& pt)
{
    auto sx = std::sin(M_PI*pt.x());
    auto sy = std::sin(M_PI*pt.y());
    auto sz = std::sin(M_PI*pt.z());
    return sx*sy*sz;
}

template<typename Mesh>
void
plot_function(const Mesh& msh, const std::string& output_filename)
{
    using T = typename Mesh::coordinate_type;
    std::vector<T> fun_vals;

    for (auto& cl : msh)
    {
        auto pt = barycenter(msh, cl);
        auto val = fun(pt);
        fun_vals.push_back( fun(pt) ); 
    }

    disk::silo_database db; /* SILO database object*/
    db.create(output_filename); /* Create a new database */
    db.add_mesh(msh, "mesh"); /* Put the mesh in it */

    /* The data is zonal, meaning that for each mesh cell there is one
     * single value. SILO/VisIt do not support high-order visualization. */
    db.add_variable("mesh", "fun", fun_vals, disk::zonal_variable_t);
    /* Nodal variables also exist, see diskpp/output/silo.hpp for the details. */
}

int main(int argc, char **argv)
{
    int num_refs = 0;
    element_type elem_type = element_type::triangle;
    std::string output_filename = "diskpp_tutorial_2.silo";

    /* Parse the command line arguments */
    int ch;
    while ( (ch = getopt(argc, argv, "r:thTf:")) != -1 )
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

            case 'f':
                output_filename = optarg;
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
            plot_function(msh, output_filename);
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
            plot_function(msh, output_filename);
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
            plot_function(msh, output_filename);
        } break;
    }

    return 0;
}

