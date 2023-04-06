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
 * Nicolas Pignet  (C) 2018, 2019               nicolas.pignet@enpc.fr
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

#ifndef _QUADRATURES_HPP_WAS_INCLUDED_
    #error "You must NOT include this file. Include quadratures.hpp"
#endif

#ifndef _QUAD_GENERIC_HPP_
#define _QUAD_GENERIC_HPP_

#include "raw_simplices.hpp"

namespace disk {

namespace priv {

template<typename Iterator>
auto
barycenter(Iterator begin, Iterator end)
{
    typedef typename std::iterator_traits<Iterator>::value_type point_type;
    typedef typename point_type::value_type                     T;

    point_type ret;
    T          den = 0.0;

    auto numpts = std::distance(begin, end);
    auto p0     = *begin;

    for (size_t i = 2; i < numpts; i++)
    {
        auto pprev = *std::next(begin, i - 1) - p0;
        auto pcur  = *std::next(begin, i) - p0;
        auto d     = det(pprev, pcur) / 2.0;
        ret        = ret + (pprev + pcur) * d;
        den += d;
    }

    return p0 + ret / (den * 3);
}

//#define OPTIMAL_TRIANGLE_NUMBER

/* The 'optimal triangle number' version gives almost the same results
 * of the other version. In bigger meshes there are some advantages in
 * assembly time. The problem with this version is that it could generate
 * triangles with a very bad aspect ratio and at the moment I don't know
 * if and how this can affect computations, so I leave it turned off.
 */
template<typename T>
std::vector<disk::quadrature_point<T, 2>>
integrate_convex_polygon(const size_t degree, const std::vector<point<T, 2>>& pts)
{
    using quadpoint_type = disk::quadrature_point<T, 2>;

    const auto qps = disk::triangle_quadrature(degree);

    std::vector<quadpoint_type> ret;

#ifdef OPTIMAL_TRIANGLE_NUMBER
    /* Break the cell in triangles, compute the transformation matrix and
     * map quadrature data in the physical space. Edges of the triangle as
     * column vectors in the transformation matrix. */
    ret.resize(qps.size() * pts.size() - 2);
    for (size_t i = 1; i < pts.size() - 1; i++)
    {
        const auto pt1  = pts[i];
        const auto pt2  = pts[i + 1];
        // Compute the integration basis
        const auto [p0, v0, v1] = integration_basis(pts[0], pt1, pt2);

        /* Compute the area of the sub-triangle */
        const auto tm = area_triangle_kahan(pts[0], pt1, pt2);

        auto tr = [ p0 = p0, v0 = v0, v1 = v1, tm ](const std::pair<point<T, 2>, T>& qd) -> auto
        {
            const auto point  = p0 + v0 * qd.first.x() + v1 * qd.first.y();
            const auto weight = qd.second * tm;
            return make_qp(point, weight);
        };

        auto retbegin = ret.begin();
        std::advance(retbegin, qps.size() * (i - 1));

        std::transform(qps.begin(), qps.end(), retbegin, tr);
    }
#else
    const auto c_center = priv::barycenter(pts.begin(), pts.end());

    /* Break the cell in triangles, compute the transformation matrix and
     * map quadrature data in the physical space. Edges of the triangle as
     * column vectors in the transformation matrix. */
    ret.resize(qps.size() * pts.size());
    for (size_t i = 0; i < pts.size(); i++)
    {
        const auto pt1  = pts[i];
        const auto pt2  = pts[(i + 1) % pts.size()];
        // Compute the integration basis
        const auto [p0, v0, v1] = integration_basis(c_center, pt1, pt2);

        /* Compute the area of the sub-triangle */
        const auto tm = area_triangle_kahan(c_center, pt1, pt2);

        auto tr = [ p0 = p0, v0 = v0, v1 = v1, tm ](const std::pair<point<T, 2>, T>& qd) -> auto
        {
            const auto point  = p0 + v0 * qd.first.x() + v1 * qd.first.y();
            const auto weight = qd.second * tm;
            return make_qp(point, weight);
        };

        auto retbegin = ret.begin();
        std::advance(retbegin, qps.size() * i);

        std::transform(qps.begin(), qps.end(), retbegin, tr);
    }
#endif
    return ret;
}

/* Integrate non-convex polygon. More costly than convex polygon
 * due to spliting
 */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
std::vector<disk::quadrature_point<T, 2>>
integrate_nonconvex_polygon(const Mesh<T, 2, Storage>&                msh,
                            const typename Mesh<T, 2, Storage>::cell& cl,
                            const size_t                              degree)
{
    using quadpoint_type = disk::quadrature_point<T, 2>;

    const auto qps = disk::triangle_quadrature(degree);

    /* Break the cell in triangles, compute the transformation matrix and
     * map quadrature data in the physical space. Edges of the triangle as
     * column vectors in the transformation matrix. */

    const auto rss = split_in_raw_triangles(msh, cl);

    std::vector<quadpoint_type> ret;
    ret.reserve(qps.size() * rss.size());

    for (auto& rs : rss)
    {
        const auto pts = rs.points();
        assert(pts.size() == 3);
        const auto meas = area_triangle_kahan(pts[0], pts[1], pts[2]);

        // Compute the integration basis
        const auto [p0, v0, v1] = integration_basis(pts[0], pts[1], pts[2]);

        for (auto& qd : qps)
        {
            const auto point  = p0 + v0 * qd.first.x() + v1 * qd.first.y();
            const auto weight = qd.second * meas;
            ret.push_back(disk::make_qp(point, weight));
        }
    }

    return ret;
}

}

template<typename T>
std::vector<disk::quadrature_point<T, 2>>
integrate(const disk::generic_mesh<T, 2>& msh, const typename disk::generic_mesh<T, 2>::cell& cl, const size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, cl);
    }

    const auto pts = points(msh, cl);
    switch (pts.size())
    {
        case 0:
        case 1:
        case 2:
            throw std::invalid_argument("A 2D cell cannot have less than three points. "
                                        "This looks like a nice bug.");

        case 3: return priv::integrate_triangle<T>(degree, pts);

        // case 4: return priv::integrate_quadrangle_tens(degree, pts);

        default: return priv::integrate_nonconvex_polygon(msh, cl, degree);
        // default: return priv::integrate_convex_polygon(degree, pts);
    }
}

template<typename T>
std::vector<disk::quadrature_point<T, 2>>
integrate(const disk::generic_mesh<T, 2>& msh, const typename disk::generic_mesh<T, 2>::face& fc, const size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, fc);
    }

    return priv::integrate_2D_face(msh, fc, degree);
}


namespace priv {

template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate_polyhedron(const disk::generic_mesh<T, 3>&                msh,
                     const typename disk::generic_mesh<T, 3>::cell& cl,
                     const size_t                                   degree)
{
    using quadpoint_type = disk::quadrature_point<T, 3>;

    const auto m_quadrature_data = disk::tetrahedron_quadrature(degree);

    const auto rss = split_in_raw_tetrahedra(msh, cl);

    std::vector<quadpoint_type> ret;
    ret.reserve(m_quadrature_data.size() * rss.size());
    for (auto& rs : rss)
    {
        const auto pts = rs.points();
        assert(pts.size() == 4);

        const auto meas = measure(rs);
        // Compute the integration basis
        const auto [p0, v0, v1, v2] = integration_basis(pts[0], pts[1], pts[2], pts[3]);

        for (auto& qd : m_quadrature_data)
        {
            auto point  = p0 + v0 * qd.first.x() + v1 * qd.first.y() + v2 * qd.first.z();
            auto weight = qd.second * meas;
            ret.push_back(make_qp(point, weight));
        }
    }

    return ret;
}

template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate_polyhedron_face(const disk::generic_mesh<T, 3>&                msh,
                          const typename disk::generic_mesh<T, 3>::face& fc,
                          const size_t                                         degree)
{
    using quadpoint_type = disk::quadrature_point<T, 3>;

    const auto m_quadrature_data = disk::triangle_quadrature(degree);

    const auto rss = split_in_raw_triangles(msh, fc);

    std::vector<quadpoint_type> ret;
    ret.reserve(m_quadrature_data.size() * rss.size());
    for (auto& rs : rss)
    {
        const auto pts = rs.points();
        assert(pts.size() == 3);
        const auto meas = measure(rs);

        // Compute the integration basis
        const auto [p0, v0, v1] = integration_basis(pts[0], pts[1], pts[2]);

        for (auto& qd : m_quadrature_data)
        {
            const auto point  = p0 + v0 * qd.first.x() + v1 * qd.first.y();
            const auto weight = qd.second * meas;
            ret.push_back(disk::make_qp(point, weight));
        }
    }

    return ret;
}

} // end priv

template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate(const disk::generic_mesh<T, 3>& msh, const typename disk::generic_mesh<T, 3>::cell& cl, const size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, cl);
    }

    const auto pts = points(msh, cl);
    switch (pts.size())
    {
        case 0:
        case 1:
        case 2:
        case 3:
            throw std::invalid_argument("A 3D cell cannot have less than four points. "
                                        "This looks like a nice bug.");

        default: return priv::integrate_polyhedron(msh, cl, degree);
    }
}

template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate(const disk::generic_mesh<T, 3>& msh, const typename disk::generic_mesh<T, 3>::face& fc, const size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, fc);
    }

    const auto pts = points(msh, fc);
    switch (pts.size())
    {
        case 0:
        case 1:
        case 2:
            throw std::invalid_argument("A 3D face cannot have less than three points. "
                                        "This looks like a nice bug.");

        default: return priv::integrate_polyhedron_face(msh, fc, degree);
    }
}

template<typename T>
std::vector<disk::quadrature_point<T, 1>>
integrate(const disk::generic_mesh<T,1>& msh, const typename generic_mesh<T, 1>::cell& cl, size_t degree)
{
    const auto qps = disk::edge_quadrature<T>(degree);
    const auto pts = points(msh, cl);

    assert(pts.size() == 2);

    const auto scale = (pts[1] - pts[0]);
    const auto meas  = scale.to_vector().norm();

    std::vector<disk::quadrature_point<T, 1>> ret;
    ret.reserve(qps.size());

    for (auto itor = qps.begin(); itor != qps.end(); itor++)
    {
        const auto qp = *itor;
        const auto t  = qp.first.x();
        const auto p  = 0.5 * (1 - t) * pts[0] + 0.5 * (1 + t) * pts[1];
        const auto w  = qp.second * meas * 0.5;

        ret.push_back(disk::make_qp(p, w));
    }

    return ret;
}

template<typename T>
std::vector<disk::quadrature_point<T, 1>>
integrate(const disk::generic_mesh<T, 1>& msh, const typename generic_mesh<T, 1>::face& fc, size_t degree)
{
    using mesh_type = generic_mesh<T, 1>;
    std::vector<disk::quadrature_point<typename mesh_type::coordinate_type, 1>> ret;

    const auto bar  = barycenter(msh, fc);
    const auto meas = 1.0;

    ret.push_back(disk::make_qp(bar, meas));

    return ret;
}

} // namespace disk

#endif /* _QUAD_GENERIC_HPP_ */
