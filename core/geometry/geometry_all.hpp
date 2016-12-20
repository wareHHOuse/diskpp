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

 #ifndef _GEOMETRY_HPP_WAS_INCLUDED_
     #error "You must NOT include this file directly. Include geometry.hpp."
 #endif

 #ifndef _GEOMETRY_ALL_HPP_
 #define _GEOMETRY_ALL_HPP_

/*
 * Here we have all the queries that don't depend on the particular mesh
 * type. For example 'barycenter()' is computed in the same way for all the
 * elements, so it goes here. It can happen that there are queries that are
 * applicable to all kinds of elements, but for only some classes of elements
 * more efficient ways to do the computation exist. Put the general one here
 * and put the specialization in the right geometry_<whatever>.hpp file.
 */

#include <algorithm>
#include <vector>

#include "common/util.h"

namespace disk {

/* Compute an estimate of the mesh discretization step 'h' */
template<typename Mesh>
typename Mesh::scalar_type
mesh_h(const Mesh& msh)
{
    typename Mesh::scalar_type h{};
    for (auto itor = msh.cells_begin(); itor != msh.cells_end(); itor++)
    {
        auto cell = *itor;
        auto cell_measure = measure(msh, cell);

        auto fcs = faces(msh, cell);
        typename Mesh::scalar_type face_sum{};
        for (auto& f : fcs)
        {
            auto m = measure(msh, f);
            face_sum += m;
        }
        h = std::max(h, cell_measure/face_sum);
    }

    return h;
}

template<typename Mesh, typename Element>
std::vector<typename Mesh::point_type>
points(const Mesh& msh, const Element& elem)
{
    auto ptids = elem.point_ids();

    auto points_begin = msh.points_begin();
    auto ptid_to_point = [&](const point_identifier<Mesh::dimension>& pi) -> auto {
        return *std::next(points_begin, pi);
    };

    std::vector<typename Mesh::point_type> pts(ptids.size());
    std::transform(ptids.begin(), ptids.end(), pts.begin(), ptid_to_point);

    return pts;
}

/* Compute the barycenter of a cell */
template<typename Mesh, typename Element>
point<typename Mesh::coordinate_type, Mesh::dimension>
barycenter(const Mesh& msh, const Element& elm)
{
    auto pts = points(msh, elm);
    auto bar = std::accumulate(std::next(pts.begin()), pts.end(), pts.front());
    return bar / typename Mesh::coordinate_type( pts.size() );
}

template<typename Mesh, typename Element>
typename Mesh::scalar_type
diameter(const Mesh& msh, const Element& elem)
{
    auto pts = points(msh, elem);

    typename Mesh::scalar_type diam = 0.;

    for (size_t i = 0; i < pts.size(); i++)
        for (size_t j = i+1; j < pts.size(); j++)
            diam = std::max((pts[i] - pts[j]).to_vector().norm(), diam);

    return diam;
}

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
static_vector<T, 2>
normal(const Mesh<T,2,Storage>& msh,
       const typename Mesh<T,2,Storage>::cell& cl,
       const typename Mesh<T,2,Storage>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 2);

    auto v = pts[1] - pts[0];
    auto n = (point<T,2>({-v.y(), v.x()})).to_vector();

    auto cell_bar = barycenter(msh, cl);
    auto face_bar = barycenter(msh, fc);
    auto outward_vector = (face_bar - cell_bar).to_vector();

    if ( n.dot(outward_vector) < T(0) )
        return -n/n.norm();

    return n/n.norm();
}


} // namespace disk

#endif /* _GEOMETRY_ALL_HPP_ */
