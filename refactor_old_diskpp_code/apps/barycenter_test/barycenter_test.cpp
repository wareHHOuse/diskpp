/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
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


template<template<typename, size_t, typename> class Mesh,
typename T, typename Storage>
void
dump_the_fucking_barycenters_to_matlab(const Mesh<T, 2, Storage>& msh,
                                       const std::string& filename)
{
    std::ofstream ofs(filename);

    ofs << "clf;\nhold on;" << std::endl;

    for (auto cl : msh)
    {
        auto bar = barycenter(msh, cl);

        ofs << "plot(" << bar.x() << ", " << bar.y() << ", 'r+');" << std::endl;

        auto fcs = faces(msh, cl);
        for (auto fc : fcs)
        {
            auto ptids = fc.point_ids();
            auto pts = points(msh, fc);
            assert(ptids.size() == pts.size());

            ofs << "line([" << pts[0].x() << " " << pts[1].x() << "], [";
            ofs << pts[0].y() << " " << pts[1].y() << "], 'Color', 'b');";
            ofs << std::endl;
        }
    }

    ofs.close();
}

int main(int argc, char **argv)
{
    using RealType = double;


    if (argc < 2)
    {
        std::cout << "Filename please!" << std::endl;
        return 1;
    }

    char * filename = argv[1];

    if (std::regex_match(filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;

        typedef disk::generic_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::fvca5_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        dump_the_fucking_barycenters_to_matlab(msh, "barycenters.m");


    }

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

        dump_the_fucking_barycenters_to_matlab(msh, "barycenters.m");
    }


    return 0;
}
