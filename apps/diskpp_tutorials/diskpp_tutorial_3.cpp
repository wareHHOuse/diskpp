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

/* This DiSk++ tutorial shows how to iterate on the elements of a mesh
 * previously loaded.
 */



#include <iostream>

// For the mesh data structure
#include "diskpp/mesh/mesh.hpp"

// For the loaders and related helper functions
#include "diskpp/loaders/loader.hpp"

#include "diskpp/output/silo.hpp"

using disk::cells;
using disk::faces;

template<typename Mesh>
void run(Mesh& msh)
{
    /* Iterate on the mesh cells */
    for (auto& cl : cells(msh))
    {
        auto h = measure(msh, cl);
        auto bar = barycenter(msh, cl);
        std::cout << cl << ". h = " << h << ", barycenter: " << bar << std::endl;
    }

    /* iterate on faces */
    for (auto& fc : faces(msh))
    {
        auto h = measure(msh, fc);
        auto bar = barycenter(msh, fc);
        std::cout << fc << ". h = " << h << ", barycenter: " << bar <<  std::endl;
    }

    disk::silo_database silo;
    silo.create("mesh.silo");
    silo.add_mesh(msh, "mesh");
    silo.close();
}

int main(int argc, const char **argv)
{
    if (argc < 2)
    {
        std::cout << "Please specify a mesh" << std::endl;
        return 1;
    }

    const char *mesh_filename = argv[1];

    using T = double;

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        disk::generic_mesh<T, 2> msh;
        auto success = disk::load_mesh_fvca5_2d<T>(mesh_filename, msh);
        if (!success)
            return 1;
        run(msh);
        return 0;
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        disk::simplicial_mesh<T, 2> msh;
        auto success = disk::load_mesh_netgen<T>(mesh_filename, msh);
        if (!success)
            return 1;
        run(msh);
        return 0;
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        disk::cartesian_mesh<T, 2> msh;
        auto success = disk::load_mesh_diskpp_cartesian<T>(mesh_filename, msh);
        if (!success)
            return 1;
        run(msh);
        return 0;
    }

    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        disk::simplicial_mesh<T, 3> msh;
        auto success = disk::load_mesh_netgen<T>(mesh_filename, msh);
        if (!success)
            return 1;
        run(msh);
        return 0;
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        disk::cartesian_mesh<T, 3> msh;
        auto success = disk::load_mesh_diskpp_cartesian<T>(mesh_filename, msh);
        if (!success)
            return 1;
        run(msh);
        return 0;
    }

    /* FVCA6 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
        disk::generic_mesh<T,3> msh;
        auto success = disk::load_mesh_fvca6_3d<T>(mesh_filename, msh);
        if (!success)
            return 1;
        run(msh);
        return 0;
    }


    return 0;

}
