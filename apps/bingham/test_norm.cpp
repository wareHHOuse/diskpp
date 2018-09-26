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
 * Karol Cascavita (C) 2018                     karol.cascavita@enpc.fr
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
#include <iomanip>
#include <regex>

#include <unistd.h>

#include "revolution/bases"
#include "revolution/quadratures"
#include "revolution/methods/hho"

#include "core/loaders/loader.hpp"


int main(void)
{
    using RealType = double;

    typedef Matrix<RealType, Dynamic, Dynamic>         matrix_type;
    typedef Matrix<RealType, Dynamic, 1>               vector_type;
    typedef disk::cartesian_mesh<RealType, 2>  mesh_type;

    auto input_mesh ="../../../diskpp/meshes/2D_quads/diskpp/testmesh-256-256.quad";

    //#if 0
    mesh_type  msh;
    disk::cartesian_mesh_loader<RealType, 2> loader;
    if (!loader.read_mesh(input_mesh))
    {
        std::cout << "Problem loading mesh." << std::endl;
        return 0;
    }
    loader.populate_mesh(msh);
    //#endif

    auto face_degree = 0;
    auto dim = 2;


    std::ofstream ofs("bars_x_driven_256_256.data");
    if (!ofs.is_open())
        std::cout << "Error opening file"<<std::endl;


    for(auto& cl : msh)
    {
        auto bar = barycenter(msh, cl);
        ofs << bar.x() << " "<< bar.y()<<std::endl;
    }

    ofs.close();

    return 0;
}
