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

#include "diskpp/geometry/geometry.hpp"
#include "diskpp/common/simplicial_formula.hpp"
#include "diskpp/contrib/earcut.hpp"

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

    return area_triangle_kahan(pts[0], pts[1], pts[2]);
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

    return volume_tetrahedron_kahan(pts[0], pts[1], pts[2], pts[3]);
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

    // Compute the integration basis
    const auto [p0, v0, v1, v2] = integration_basis(pts[0], pts[1], pts[2], pts[3]);

    return p0 + v0 * p.x() + v1 * p.y() + v2 * p.z();
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

    // Compute the integration basis
    const auto [p0, v0, v1] = integration_basis(pts[0], pts[1], pts[2]);

    return v0 * p.x() +
              v1 * p.y() +
              p0;
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
    static_assert(Mesh::dimension == 3);
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
 * @brief Split a 2D-cell in triangles. There is several possibilities
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl specified 2D-cell
 * @return std::vector<raw_simplex<typename Mesh::point_type, 2>> list of triangles
 */
template<typename Mesh>
std::vector<raw_simplex<typename Mesh::point_type, 2>>
split_in_raw_triangles(const Mesh& msh, const typename Mesh::cell& cl)
{
    static_assert(Mesh::dimension == 2);
    typedef raw_simplex<typename Mesh::point_type, 2> raw_simplex_type;

    auto pts = points(msh, cl);

    std::vector<typename Mesh::point_type> vpts;
    vpts.insert(vpts.begin(), std::begin(pts), std::end(pts));

    const auto triangles = triangulation_earcut_algorithm(vpts);
    std::vector<raw_simplex_type> raw_simplices;
    raw_simplices.reserve(triangles.size());

    for(auto& tri : triangles)
    {
        raw_simplices.push_back(raw_simplex_type({pts[tri[0]], pts[tri[1]], pts[tri[2]]}));
    }

    return raw_simplices;
}

/**
 * @brief Split a 2D-cell in triangles.
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl specified 2D-cell
 * @return std::vector<std::array<point_identifier<2>, 3>> list of triangles (index of points)
 */
template<typename Mesh>
std::vector<std::array<point_identifier<2>, 3>>
split_in_raw_triangles_index(const Mesh& msh, const typename Mesh::cell& cl)
{
    static_assert(Mesh::dimension == 2);
    using Triangle = std::array<point_identifier<2>, 3>;

    const auto pts = points(msh, cl);
    const auto pts_id = cl.point_ids();

    std::vector<typename Mesh::point_type> vpts;
    vpts.insert(vpts.begin(), std::begin(pts), std::end(pts));

    const auto triangles = triangulation_earcut_algorithm(vpts);

    std::vector<Triangle> raw_simplices;
    raw_simplices.reserve(triangles.size());

    for (auto& tri : triangles)
    {
        raw_simplices.push_back(Triangle{
            point_identifier<2>(pts_id[tri[0]]),
            point_identifier<2>(pts_id[tri[1]]),
            point_identifier<2>(pts_id[tri[2]])});
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
    static_assert(Mesh::dimension == 3);
    typedef raw_simplex<typename Mesh::point_type, 3> raw_simplex_type;

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

/**
 * @brief Split a 2D-cell in triangles using earcut algorithm
 *
 * Return the list of triangles (local index of nodes)
 *
 * There is (n-2) triangles for n points
 *
 * @tparam scalar type
 * @param std::vector<point<T, 2>>& list of points of the cell
 * @return std::vector<std::array<std::size_t, 3>> list of triangles
 */
template<typename T>
std::vector<std::array<std::size_t, 3>>
triangulation_earcut_algorithm(const std::vector<point<T, 2>>& pts)
{
    // Create array
    using Point = std::array<T, 2>;
    using Triangle = std::array<std::size_t, 3>;
    std::vector<std::vector<Point>> polygon;

    std::vector<Point> list_pts;
    list_pts.reserve(pts.size());

    for( auto& pt : pts)
    {
        list_pts.push_back(Point{{pt.x(), pt.y()}});
    }

    // Fill polygon structure with actual data. Any winding order works.
    // The first polyline defines the main polygon.
    polygon.push_back(list_pts);
    // Following polylines define holes.
    // polygon.push_back({{75, 25}, {75, 75}, {25, 75}, {25, 25}});

    // Run tessellation
    // Returns array of indices that refer to the vertices of the input polygon.
    // e.g: the index 6 would refer to {25, 75} in this example.
    // Three subsequent indices form a triangle. Output triangles are clockwise.
    const auto indices = mapbox::earcut<std::size_t>(polygon);
    const auto nb_ind  = indices.size();

    std::vector<Triangle> triangles;
    triangles.reserve(nb_ind / 3);

    for (std::size_t i = 0; i < nb_ind; i+=3)
    {
        triangles.push_back(Triangle{{indices[i], indices[i + 1], indices[i + 2]}});
    }

    // std::cout << "pts: ";
    // for (auto& pt : pts)
    //     std ::cout << "(" << pt.x() << "," << pt.y() << ")";
    // std::cout << std::endl;
    // std::cout << "ind: ";

    // for (auto& i : indices)
    //     std ::cout << i << ", ";
    // std::cout << std::endl;

    return triangles;
}

    } // end disk
