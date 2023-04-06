/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016 - 2020
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet (C) 2020        nicolas.pignet@enpc.fr
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

#include <iomanip>
#include <iostream>
#include <sstream>
#include <regex>
#include <unistd.h>

#include "diskpp/loaders/loader.hpp"

#include "diskpp/output/gmshConvertMesh.hpp"
#include "diskpp/output/gmshDisk.hpp"
#include "diskpp/output/postMesh.hpp"

template<typename MeshType>
void split_2D_polygons(MeshType& msh)
{
    std::cout << "Split 2D polygons" << std::endl;

    // compute mesh for post-processing
    const disk::PostMesh<MeshType>  post_mesh = disk::PostMesh<MeshType>(msh);
    const gmsh::Gmesh               gmsh      = disk::convertMesh(post_mesh);

    gmsh.writeGmesh("splittedMesh.msh", 2);
}

  int
  main(int argc, char** argv)
{
    using RealType = double;

    char* mesh_filename = nullptr;

    int ch;

    while ((ch = getopt(argc, argv, "r:")) != -1)
    {
        switch (ch)
        {
            case 'r':
                break;
            default: std::cout << "wrong arguments" << std::endl; exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    if (argc == 0)
    {
        std::cout << "Error" << std::endl;
        return 0;
    }

    mesh_filename = argv[0];

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$")))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        auto msh = disk::load_fvca5_2d_mesh<RealType>(mesh_filename);
        split_2D_polygons(msh);
        return 0;
    }

    /* Medit 2d*/
    if (std::regex_match(mesh_filename, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        auto msh = disk::load_medit_2d_mesh<RealType>(mesh_filename);
        split_2D_polygons(msh);
        return 0;
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$")))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);
        split_2D_polygons(msh);
        return 0;
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$")))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        auto msh = disk::load_cartesian_2d_mesh<RealType>(mesh_filename);
        split_2D_polygons(msh);
        return 0;
    }

}