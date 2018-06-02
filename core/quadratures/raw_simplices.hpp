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

#pragma once

#include "geometry/geometry.hpp"
#include "mesh/point.hpp"


template<typename PtT, size_t N>
class raw_simplex
{
    typedef PtT                                 point_type;
    typedef typename point_type::value_type     coordinate_type;

    std::array<point_type, N+1>  m_points;

public:
    raw_simplex()
    {}

    raw_simplex(std::initializer_list<point_type> l)
    {
        if (l.size() != N+1)
            throw std::invalid_argument("Wrong initializer list size");
        std::copy(l.begin(), l.end(), m_points.begin());
    }

    auto points() const { return m_points; }

};

template<typename T, size_t N>
auto
points(const raw_simplex<T, N>& rs)
{
    return rs.points();
}

template<typename T>
auto
measure(const raw_simplex<T,1>& rs)
{
    auto pts = rs.points();
    assert(pts.size() == 2);

    return (pts[1] - pts[0]).to_vector().norm();
}

template<typename T>
auto
measure(const raw_simplex<T,2>& rs)
{
    auto pts = rs.points();
    assert(pts.size() == 3);

    auto v0 = (pts[1] - pts[0]).to_vector();
    auto v1 = (pts[2] - pts[0]).to_vector();

    return v0.cross(v1).norm()/2.0;
}

template<typename T>
auto
measure(const raw_simplex<point<T,3>,3>& rs)
{
    auto pts = rs.points();
    assert(pts.size() == 4);

    auto v0 = (pts[1] - pts[0]).to_vector();
    auto v1 = (pts[2] - pts[0]).to_vector();
    auto v2 = (pts[3] - pts[0]).to_vector();

    return std::abs( v0.dot(v1.cross(v2))/T(6) );
}

template<typename T, size_t N>
T
barycenter(const raw_simplex<T,N>& rs)
{
    auto pts = rs.points();
    auto bar = std::accumulate(std::next(pts.begin()), pts.end(), pts.front());
    return bar / pts.size();
}

template<typename T>
point<T,3>
map_from_reference(const point<T,3>& p, const raw_simplex<point<T,3>, 3>& rs)
{
    auto pts = rs.points();
    assert(pts.size() == 4);

    auto pp = (pts[1] - pts[0]) * p.x() +
              (pts[2] - pts[0]) * p.y() +
              (pts[3] - pts[0]) * p.z() +
              pts[0];

    return pp;
}

template<typename T>
point<T,3>
map_from_reference(const point<T,2>& p, const raw_simplex<point<T,3>, 2>& rs)
{
    auto pts = rs.points();
    assert(pts.size() == 3);

    auto pp = (pts[1] - pts[0]) * p.x() +
              (pts[2] - pts[0]) * p.y() +
              pts[0];

    return pp;
}

template<typename Mesh>
std::vector<raw_simplex<typename Mesh::point_type, 2>>
split_in_raw_triangles(const Mesh& msh, const typename Mesh::face& fc)
{
    typedef raw_simplex<typename Mesh::point_type, 2>   raw_simplex_type;

    auto pts    = points(msh, fc);

    std::vector<raw_simplex_type>   raw_simplices;

    if (pts.size() == 3)
        raw_simplices.push_back( raw_simplex_type({pts[0], pts[1], pts[2]}) );
    else if (pts.size() == 4)
    {
        raw_simplices.push_back( raw_simplex_type({pts[0], pts[1], pts[2]}) );
        raw_simplices.push_back( raw_simplex_type({pts[2], pts[3], pts[0]}) );
    }
    else
    {
        auto bar = barycenter(msh, fc);
        for (size_t i = 0; i < pts.size(); i++)
        {
            auto rs = raw_simplex_type( {bar, pts[i], pts[(i+1)%pts.size()]} );
            raw_simplices.push_back( rs );
        }
    }

    return raw_simplices;
}

template<typename Mesh>
std::vector<raw_simplex<typename Mesh::point_type, 3>>
split_in_raw_tetrahedra(const Mesh& msh, const typename Mesh::cell& cl)
{
    typedef raw_simplex<typename Mesh::point_type, 3>   raw_simplex_type;

    auto pts    = points(msh, cl);
    assert(pts.size() >= 4);

    std::vector<raw_simplex_type> raw_simplices;
    if(pts.size() == 4)
    {
        raw_simplices.push_back(raw_simplex_type({pts[0], pts[1], pts[2], pts[3]}));
    }
    else
    {
        auto bar = barycenter(msh, cl);
        auto fcs = faces(msh, cl);

        for (auto& fc : fcs)
        {
            auto tris = split_in_raw_triangles(msh, fc);

            for (auto& tri : tris)
            {
                auto tpts = points(tri);
                auto rs   = raw_simplex_type({bar, tpts[0], tpts[1], tpts[2]});
                raw_simplices.push_back(rs);
            }
        }
    }

    return raw_simplices;
}
