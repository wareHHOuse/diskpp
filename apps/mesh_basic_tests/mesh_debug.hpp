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

#pragma once

#include <sstream>

#include "mesh/mesh.hpp"
#include "simpleps.h"
#include "geometry/geometry.hpp"

template<typename T>
T transform(T point, T min, T max, T scale, T offset)
{
    //return offset + (point - min)*scale/(max-min);
    return (point+min)*scale + offset;
}

template<typename T, typename Storage>
void mesh_to_postscript(const disk::mesh<T, 2, Storage>& msh, const std::string& filename)
{
    auto storage = msh.backend_storage();
    point<T, 2> min, max;

    for (auto& pt : storage->points)
    {
        if ( pt.x() < min.x() ) min.x() = pt.x();
        if ( pt.y() < min.y() ) min.y() = pt.y();
        if ( max.x() < pt.x() ) max.x() = pt.x();
        if ( max.y() < pt.y() ) max.y() = pt.y();
    }

    std::cout << "bounding box: " << min << " " << max << std::endl;


    sps::simpleps ps;

    for (auto& e : storage->edges)
    {
        auto e_center = barycenter(msh, e);
        auto n1 = *e.subelement_id_begin();
        auto n2 = *std::next(e.subelement_id_begin());

        auto pt1 = *(storage->points.begin() + n1);
        auto pt2 = *(storage->points.begin() + n2);

        T pt1_newx = transform(pt1.x(), min.x(), max.x(), 500.0, 10.0);
        T pt1_newy = transform(pt1.y(), min.y(), max.y(), 500.0, 10.0);
        T pt2_newx = transform(pt2.x(), min.x(), max.x(), 500.0, 10.0);
        T pt2_newy = transform(pt2.y(), min.y(), max.y(), 500.0, 10.0);

        T e_center_newx = transform(e_center.x(), min.x(), max.x(), 500.0, 10.0);
        T e_center_newy = transform(e_center.y(), min.y(), max.y(), 500.0, 10.0);

        sps::path *p = new sps::path( {sps::ps_point( {pt1_newx, pt1_newy} ),
                                       sps::ps_point( {pt2_newx, pt2_newy} ) } );

        std::stringstream ss;
        ss << e;
        sps::text *t = new sps::text(ss.str(), sps::ps_point( {e_center_newx, e_center_newy} ) );

        ps.add(p);
        ps.add(t);
    }

    ps.write(filename);
}
