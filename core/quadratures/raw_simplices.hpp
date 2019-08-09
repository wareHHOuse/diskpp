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

#include "geometry/geometry.hpp"
#include "mesh/point.hpp"

namespace disk{


/**
 * @brief This class contains the basic informations to define a simplex, i.e. a triangle
 * in 2D and a tetrahedra in 3D.
 *
 * @tparam PtT type of point
 * @tparam N dimention of the simplex
 */
template<typename PtT, size_t N>
class raw_simplex
{
    typedef PtT                                 point_type;
    typedef typename point_type::value_type     coordinate_type;

    std::array<point_type, N+1>  m_points;

public:
    /**
     * @brief Construct a new raw simplex object
     *
     * Default contructor
     */
    raw_simplex()
    {}

    /**
     * @brief Construct a new raw simplex object
     *
     * Create a simplex from the list of its vertex
     *
     * @param l list of points
     */
    raw_simplex(std::initializer_list<point_type> l)
    {
        if (l.size() != N+1)
            throw std::invalid_argument("Wrong initializer list size");
        std::copy(l.begin(), l.end(), m_points.begin());
    }

    /**
     * @brief Return the list of the vertex
     *
     * @return std::array<point_type, N+1> vertex
     */
    std::array<point_type, N+1> points() const { return m_points; }

};

/**
 * @brief Compute the measure, i.e. the area of the specified triangle
 *
 * @tparam PtT type of point
 * @param rs specified triangle
 * @return T area of the triangle
 */
template<typename PtT>
auto
measure(const raw_simplex<PtT,2>& rs)
{
    const auto pts = rs.points();
    assert(pts.size() == 3);

    const auto v0 = (pts[1] - pts[0]).to_vector();
    const auto v1 = (pts[2] - pts[0]).to_vector();

    return v0.cross(v1).norm()/2.0;
}

/**
 * @brief Compute the measure, i.e. the volume of the specified tetrahedra
 *
 * @tparam PtT type of point
 * @param rs specified tetrahedra
 * @return T area of the tetrahedra
 */
template<typename PtT>
auto
measure(const raw_simplex<PtT,3>& rs)
{
    using T        = typename PtT::value_type;
    const auto pts = rs.points();
    assert(pts.size() == 4);

    const auto v0 = (pts[1] - pts[0]).to_vector();
    const auto v1 = (pts[2] - pts[0]).to_vector();
    const auto v2 = (pts[3] - pts[0]).to_vector();

    return std::abs( v0.dot(v1.cross(v2))/T(6) );
}

/**
 * @brief Compute the barycenter of the specfied simplex
 *
 * @tparam T scalar type
 * @tparam N dimension of the simplex
 * @param rs specified simplex
 * @return point<T,N> barycenter
 */
template<typename T, size_t N>
point<T,N>
barycenter(const raw_simplex<point<T,N>,N>& rs)
{
    const auto pts = rs.points();
    const auto bar = std::accumulate(std::next(pts.begin()), pts.end(), pts.front());
    return bar / T(pts.size());
}

/**
 * @brief Map a point from the reference tetrahedra to the physical tetrahedra
 *
 * @tparam T scalar type
 * @param rs specified tetrahedra
 * @param p point in the reference tetrahedra
 * @return point<T, 3> point in the physical tetrahedra
 */
template<typename T>
point<T,3>
map_from_reference(const raw_simplex<point<T,3>, 3>& rs, const point<T,3>& p)
{
    const auto pts = rs.points();
    assert(pts.size() == 4);

    auto pp = (pts[1] - pts[0]) * p.x() +
              (pts[2] - pts[0]) * p.y() +
              (pts[3] - pts[0]) * p.z() +
              pts[0];

    return pp;
}

/**
 * @brief Map a point from the reference triangle to the physical triangle
 *
 * @tparam T scalar type
 * @param rs specified triangle
 * @param pm point in the reference triangle
 * @return point<T, 3> point in the physical triangle
 */
template<typename T>
point<T,3>
map_from_reference(const point<T,2>& p, const raw_simplex<point<T,3>, 2>& rs)
{
    const auto pts = rs.points();
    assert(pts.size() == 3);

    auto pp = (pts[1] - pts[0]) * p.x() +
              (pts[2] - pts[0]) * p.y() +
              pts[0];

    return pp;
}

/**
 * @brief Split a 3D-face in triangles. There is several possibilities
 *
 * If the 3D-faces has:
 * - 3 vertices, this is a triangle and only one triangle is created
 * - 4 verticies, this is a quadrangle which is splitted in two triangles
 * - otherwise, a triangle is created from the two vertices of an edge of the triangle
 * and the barycenter of the traingle
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param fc specified 3D-face
 * @return std::vector<raw_simplex<typename Mesh::point_type, 2>> list of triangles
 */
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
        const auto bar = barycenter(msh, fc);
        for (size_t i = 0; i < pts.size(); i++)
        {
            auto rs = raw_simplex_type( {bar, pts[i], pts[(i+1)%pts.size()]} );
            raw_simplices.push_back( rs );
        }
    }

    return raw_simplices;
}

/**
 * @brief Split a 3D-cell in tetrahedra. There is several possibilities
 *
 * If the 3D-cell has:
 * - 4 vertices, this is a tetrahedra and only one tetrahedra is created
 * - otherwise, a tetrahedra is created from the two vertices of an edge of the face
 * and the barycenter of the cell and the face
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl specified cell
 * @return std::vector<raw_simplex<typename Mesh::point_type, 3>> list of tetrahedra
 */
template<typename Mesh>
std::vector<raw_simplex<typename Mesh::point_type, 3>>
split_in_raw_tetrahedra(const Mesh& msh, const typename Mesh::cell& cl)
{
    typedef raw_simplex<typename Mesh::point_type, 3>   raw_simplex_type;

    const auto pts    = points(msh, cl);
    assert(pts.size() >= 4);

    std::vector<raw_simplex_type> raw_simplices;
    if(pts.size() == 4)
    {
        raw_simplices.push_back(raw_simplex_type({pts[0], pts[1], pts[2], pts[3]}));
    }
    else
    {
        const auto bar = barycenter(msh, cl);
        const auto fcs = faces(msh, cl);

        for (auto& fc : fcs)
        {
            const auto tris = split_in_raw_triangles(msh, fc);

            for (auto& tri : tris)
            {
                const auto tpts = tri.points();
                const auto rs   = raw_simplex_type({bar, tpts[0], tpts[1], tpts[2]});
                raw_simplices.push_back(rs);
            }
        }
    }

    return raw_simplices;
}

} // end disk
