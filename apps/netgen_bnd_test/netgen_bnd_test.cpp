/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
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

#include <iostream>
#include <regex>

#include "loaders/loader.hpp"


template<typename MeshType>
void
process_mesh(const MeshType& msh)
{
    for (auto& cl : msh)
    {
        std::cout << "This is cell " << cl << std::endl;
        auto fcs = faces(msh, cl);
        for (auto& fc : fcs)
        {
            std::cout << "  This is face " << fc;
            if ( msh.is_boundary(fc) )
                std::cout << " [boundary, id = " << msh.boundary_id(fc) << "]";
            std::cout << std::endl;
        }
    }
}


int main(int argc, char **argv)
{
    using RealType = double;

    char    *filename       = nullptr;
    int     elems_1d        = 8;
    int ch;


    if (argc != 2)
    {
        std::cout << "Please specify filename" << std::endl;
        return 1;
    }

    filename = argv[1];

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

        process_mesh(msh);
    }

    return 0;
}
