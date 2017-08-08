/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#include <iostream>
#include <fstream>
#include <regex>
#include "loaders/loader.hpp"

int main(int argc, char **argv)
{

    char *filename = argv[1];

    if (std::regex_match(filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;

        typedef disk::generic_mesh<double, 3>   mesh_type;

        mesh_type msh;
        disk::fvca6_mesh_loader<double, 3> loader;


        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }

        loader.populate_mesh(msh);

        std::ofstream ofs(argv[2]);

        ofs << msh.points_size() << std::endl;

        for (auto itor = msh.points_begin(); itor != msh.points_end(); itor++)
            ofs << itor->x() << " " << itor->y() << " " << itor->z() << std::endl;

        ofs << msh.cells_size() << std::endl;
        for (auto& cl : msh)
        {
            auto ptids = cl.point_ids();
            if (ptids.size() != 4)
                throw std::invalid_argument("mesh is not tetrahedral");

            ofs << 1 << " ";
            ofs << ptids[0]+1 << " " << ptids[1]+1 << " ";
            ofs << ptids[2]+1 << " " << ptids[3]+1 << std::endl;

        }

        ofs << msh.boundary_faces_size() << std::endl;
        for (auto itor = msh.boundary_faces_begin();
             itor != msh.boundary_faces_end(); itor++)
        {
            auto fc = *itor;
            auto ptids = fc.point_ids();
            if (ptids.size() != 3)
                throw std::invalid_argument("mesh is not tetrahedral");

            ofs << 1 << " ";
            ofs << ptids[0]+1 << " " << ptids[1]+1 << " ";
            ofs << ptids[2]+1 << " " << std::endl;

        }

        ofs.close();
    }
}
