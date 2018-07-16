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
#include <fstream>

#include "geometry/geometry.hpp"
#include "loaders/loader.hpp"
#include "mesh_debug.hpp"

template<template<typename, size_t, typename> class Mesh, typename T, size_t DIM, typename Storage>
bool test_mesh(Mesh<T, DIM, Storage>& msh)
{
    for (auto& cell : msh)
    {
        std::cout << "* Cell of type " << cell << std::endl;
        std::cout << "  Measure: " << measure(msh, cell) << std::endl;
        std::cout << "  Barycenter: " << barycenter(msh, cell) << std::endl;

        auto fcs = faces(msh, cell);

        for (auto& f : fcs)
        {
            std::cout << "  * Face of type " << f << std::endl;
            std::cout << "    Measure: " << measure(msh, f) << std::endl;
            std::cout << "    Barycenter: " << barycenter(msh, f) << std::endl;
            std::cout << "    Is boundary? " << msh.is_boundary(f) << std::endl;
        }
    }

    std::cout << "Num. of boundary faces: " << msh.boundary_faces_size() << std::endl;
    std::cout << "Num. of internal faces: " << msh.internal_faces_size() << std::endl;
    //std::cout << "Mesh h: " << average_diameter(msh) << std::endl;

    size_t bcount = 0;
    for(auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
        bcount++;
    std::cout << "bcount: " << bcount << std::endl;

    size_t icount = 0;
    for(auto itor = msh.internal_faces_begin(); itor != msh.internal_faces_end(); itor++)
        icount++;
    std::cout << "icount: " << icount << std::endl;

    return true;
}

template<typename PointT>
struct point_transformer
{
    PointT operator()(const PointT& pt)
    {
        return PointT{pt.x()*cos(.1) - pt.y()*sin(.1), pt.x()*sin(.1) + pt.y()*cos(.1)};
    }
};

int main(void)
{
    using T = double;
    const size_t DIM_s = 2;

    typedef typename disk::simplicial_mesh<T, DIM_s>::point_type point_type;
    typedef point_transformer<point_type> transformer_type;
    transformer_type my_transformer;

    disk::cartesian_mesh_loader<T, DIM_s> loader_s;
    disk::cartesian_mesh<T, DIM_s> msh_s;

    //disk::medit_mesh_loader<T, 2> loader_s;
    //disk::generic_mesh<T, 2> msh_s;

    loader_s.read_mesh("../../../diskpp/meshes/2D_quads/diskpp/cube-50-50.quad");
    loader_s.populate_mesh(msh_s);
    //mesh_to_postscript(msh_s, "cube_diff_1.ps");

    test_mesh(msh_s);
    dump_to_matlab(msh_s, "test_mesh.m");

    return 0;
}
