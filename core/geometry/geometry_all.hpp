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
[[deprecated("Use average_diameter(). This is not the correct way to compute h.")]]
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

template<typename Mesh>
typename Mesh::scalar_type
average_diameter(const Mesh& msh)
{
    typename Mesh::scalar_type h{};
    for (auto& cl : msh)
    {
        h += diameter(msh, cl);
    }

    return h/msh.cells_size();
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
    const auto pts = points(msh, elem);

    typename Mesh::scalar_type diam = 0.;

    for (size_t i = 0; i < pts.size(); i++)
        for (size_t j = i+1; j < pts.size(); j++)
            diam = std::max((pts[i] - pts[j]).to_vector().norm(), diam);

    return diam;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
std::array<T, 3>
diameter_boundingbox(const Mesh<T, 3, Storage>& msh, const typename Mesh<T, 3, Storage>::cell& cl)
{
    const auto pts = points(msh, cl);

    T xmin = 0;
    T xmax = 0;
    T ymin = 0;
    T ymax = 0;
    T zmin = 0;
    T zmax = 0;

    for (auto& pt : pts)
    {
        if (pt.x() < xmin)
        {
            xmin = pt.x();
        }
        else if (pt.x() > xmax)
        {
            xmax = pt.x();
        }

        if (pt.y() < ymin)
        {
            ymin = pt.y();
        }
        else if (pt.y() > ymax)
        {
            ymax = pt.y();
        }

        if (pt.z() < zmin)
        {
            zmin = pt.z();
        }
        else if (pt.z() > zmax)
        {
            zmax = pt.z();
        }
    }

    return {std::abs(xmax - xmin), std::abs(ymax - ymin), std::abs(zmax - zmin)};
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
std::array<T,2>
diameter_boundingbox(const Mesh<T, 2, Storage>&                      msh,
                     const typename Mesh<T, 2, Storage>::cell&       cl)
{
    const auto pts = points(msh, cl);

    T xmin = 0;
    T xmax = 0;
    T ymin = 0;
    T ymax = 0;

    for (auto& pt : pts)
    {
        if(pt.x() < xmin)
        {
            xmin = pt.x();
        }
        else if (pt.x() > xmax)
        {
            xmax = pt.x();
        }

        if (pt.y() < ymin)
        {
            ymin = pt.y();
        }
        else if (pt.y() > ymax)
        {
            ymax = pt.y();
        }
    }

    return {std::abs(xmax - xmin), std::abs(ymax - ymin)};
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
bool
is_inside(const Mesh<T, 2, Storage>&                      msh,
          const typename Mesh<T, 2, Storage>::cell&       cl,
          const typename Mesh<T, 2, Storage>::point_type& pt)
{
    /* Nodes MUST be ordered COUNTERCLOCKWISE and the polygon must be CONVEX */
    auto pts = points(msh, cl);

    for (size_t i = 0; i < pts.size(); i++)
    {
        auto p0 = pts[i];
        auto p1 = pts[i % pts.size()];

        auto x  = pt.x();
        auto y  = pt.y();
        auto x0 = p0.x();
        auto y0 = p0.y();
        auto x1 = p1.x();
        auto y1 = p1.y();

        if ((y - y0) * (x1 - x0) - (x - x0) * (y1 - y0) < 0.0)
            return false;
    }

    return true;
}

template<typename Mesh>
bool
has_faces_on_boundary(const Mesh& msh, const typename Mesh::cell& cl)
{
    auto fcs = faces(msh, cl);
    bool has_bnd = false;
    for (auto& fc : fcs)
        if ( msh.is_boundary(fc) )
            return true;

    return has_bnd;
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

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
static_vector<T, 3>
normal(const Mesh<T, 3, Storage>& msh,
       const typename Mesh<T, 3, Storage>::cell& cl,
       const typename Mesh<T, 3, Storage>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() >= 3);

    auto v0 = (pts[1] - pts[0]).to_vector();
    auto v1 = (pts[2] - pts[1]).to_vector();
    auto n = v0.cross(v1);

    auto cell_bar = barycenter(msh, cl);
    auto face_bar = barycenter(msh, fc);
    auto outward_vector = (face_bar - cell_bar).to_vector();

    if ( n.dot(outward_vector) < T(0) )
        return -n/n.norm();

    return n/n.norm();
}

} // namespace disk

#endif /* _GEOMETRY_ALL_HPP_ */
