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

#pragma once

#include "quad_raw_dunavant.hpp"
#include "quad_raw_gauss.hpp"


namespace disk {
namespace quadrature {

namespace priv {

/* Integrate a convex element: determine a rough center and build triangles
 * between it and all the edges. Then use a triangle quadrature. */
template<typename T>
auto
integrate_convex(const disk::generic_mesh<T,2>& msh,
    const typename disk::generic_mesh<T,2>::cell_type& cl, size_t degree)
{
    auto pts = points(msh, cl);
    assert(pts.size() > 2);
    auto center = std::accumulate(pts.begin(), pts.end(), point<T,2>(0,0));
    center = center/T(pts.size());

    std::vector<quadrature_point<T,2>> ret;

    for (size_t i = 0; i < pts.size(); i++)
    {
        auto p0 = pts[i];
        auto p1 = pts[(i+1)%pts.size()];
        auto qps = dunavant(degree, center, p0, p1);
        ret.insert( ret.end(), qps.begin(), qps.end() );
    }

    return ret;
}

/* Integrate a non-convex element: triangulate by calling a mesh generator.
 * Then use a triangle quadrature. */
template<typename T>
auto
integrate_nonconvex(const disk::generic_mesh<T,2>& msh,
    const typename disk::generic_mesh<T,2>::cell_type& cl, size_t degree)
{
    auto tris = triangulate_nonconvex_polygon(msh, cl);

    std::vector<quadrature_point<T,2>> ret;
    for (auto& tri : tris)
    {
        auto qps = dunavant(degree, tri.p0, tri.p1, tri.p2);
        ret.insert( ret.end(), qps.begin(), qps.end() );
    }

    return ret;
}

} // namespace priv

} // namespace quadrature


template<typename T>
std::vector<disk::quadrature_point<T, 2>>
integrate(const disk::generic_mesh<T, 2>& msh, const typename disk::generic_mesh<T, 2>::cell& cl, size_t degree)
{
    const auto pts = points(msh, cl);

    assert( (pts.size() > 2) && "Insufficient points for a 2D cell" );

    if (pts.size() == 3) {
        return disk::quadrature::dunavant(degree, pts[0], pts[1], pts[2]);
    }

    bool convex = is_convex(msh, cl);
    
    if (pts.size() == 4 and convex) {
        return disk::quadrature::tensorized_gauss_legendre(degree, pts[0], pts[1], pts[2], pts[3]);
    }

    if (convex) {
        return quadrature::priv::integrate_convex(msh, cl, degree);
    }

    return quadrature::priv::integrate_nonconvex(msh, cl, degree);
}

template<typename T>
auto
integrate(const disk::generic_mesh<T, 2>& msh,
    const typename disk::generic_mesh<T, 2>::face& fc, size_t degree)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 2);
    return disk::quadrature::gauss_legendre(degree, pts[0], pts[1]);
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

