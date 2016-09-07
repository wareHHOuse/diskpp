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

#include "geometry/geometry_generic.hpp"
#include "geometry/geometry_simplicial.hpp"
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
    //std::cout << "Mesh h: " << mesh_h(msh) << std::endl;

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
    const size_t DIM_g = 2;
    const size_t DIM_s = 3;

    typedef typename disk::generic_mesh<T, DIM_g>::point_type point_type;
    typedef point_transformer<point_type> transformer_type;
    transformer_type my_transformer;

    disk::generic_mesh<T, DIM_g> msh_g;
    disk::fvca5_mesh_loader<T, DIM_g> loader_g;
    loader_g.read_mesh("../../../meshes/mesh2_1.typ1");
    //loader_g.read_mesh("./data/pi6_tiltedhexagonal_1.typ1");
    loader_g.populate_mesh(msh_g);
    msh_g.transform(my_transformer);
    mesh_to_postscript(msh_g, "mesh2_1.ps");

//    hho::simplicial_mesh<T, DIM_s> msh_s;
//    hho::netgen_mesh_loader<T, DIM_s> loader_s;
//    loader_s.read_mesh("../../hho/tests/data/simplecube.mesh");
//    loader_s.populate_mesh(msh_s);

//    hho::generic_mesh<T, 1> msh_u;
//    hho::uniform_mesh_loader<T, 1> loader_u;
//    loader_u.populate_mesh(msh_u);

    test_mesh(msh_g);
//    test_mesh(msh_s);
//    test_mesh(msh_u);

    return 0;
}
