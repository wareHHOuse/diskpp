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

#include "quad_raw_tetra.hpp"
#include "quad_raw_triangle.hpp"

namespace disk
{
namespace quadrature
{

namespace priv
{

/* Integrate a convex element: determine a rough center and build triangles
 * between it and all the edges. Then use a triangle quadrature. */
template<typename T>
auto
integrate_convex(const disk::generic_mesh<T, 2>&                     msh,
                 const typename disk::generic_mesh<T, 2>::cell_type& cl,
                 size_t                                              degree)
{
    auto pts = points(msh, cl);
    assert(pts.size() > 2);
    std::vector< quadrature_point< T, 2 > > ret;

    if ( pts.size() == 3 ) {
        auto p0 = pts[0];
        auto p1 = pts[1];
        auto p2 = pts[2];
        auto qps = triangle_gauss( degree, p0, p1, p2 );
        ret.insert( ret.end(), qps.begin(), qps.end() );
        return ret;
    }

    auto center = std::accumulate(pts.begin(), pts.end(), point<T, 2>(0, 0));
    center      = center / T(pts.size());

    for (size_t i = 0; i < pts.size(); i++)
    {
        auto p0  = pts[i];
        auto p1  = pts[(i + 1) % pts.size()];
        auto qps = triangle_gauss(degree, center, p0, p1);
        ret.insert(ret.end(), qps.begin(), qps.end());
    }

    return ret;
}

/* Integrate a non-convex element: triangulate by calling a mesh generator.
 * Then use a triangle quadrature. */
template<typename T>
auto
integrate_nonconvex(const disk::generic_mesh<T, 2>&                     msh,
                    const typename disk::generic_mesh<T, 2>::cell_type& cl,
                    size_t                                              degree)
{
    auto tris = triangulate_nonconvex_polygon(msh, cl);

    std::vector<quadrature_point<T, 2>> ret;
    for (auto& tri : tris)
    {
        auto qps = triangle_gauss(degree, tri.p0, tri.p1, tri.p2);
        ret.insert(ret.end(), qps.begin(), qps.end());
    }

    return ret;
}

template < typename T >
bool is_ortho_quad( const disk::generic_mesh< T, 2 > &msh,
                    const typename disk::generic_mesh< T, 2 >::cell_type &cl ) {
    const auto pts = points( msh, cl );

    if ( pts.size() != 4 ) {
        return false;
    };

    const T thrs = 1e-8;

    const auto v02 = ( pts[2] - pts[0] ).to_vector();

    const auto p2_test = pts[1] + ( pts[3] - pts[0] );

    const auto dist = ( pts[2] - p2_test ).to_vector().norm();

    return dist < thrs * v02.norm();
}

template < typename T >
bool is_hexa( const disk::generic_mesh< T, 3 > &msh,
              const typename disk::generic_mesh< T, 3 >::cell_type &cl ) {

    if ( cl.point_ids().size() != 8 ) {
        return false;
    }

    const auto fcs = faces( msh, cl );
    if ( fcs.size() != 6 ) {
        return false;
    }

    for ( auto &fc : fcs ) {
        if ( fc.point_ids().size() != 4 ) {
            return false;
        }
    }

    return true;
}

template < typename T >
bool is_ortho_hexa( const disk::generic_mesh< T, 3 > &msh,
                    const typename disk::generic_mesh< T, 3 >::cell_type &cl ) {

    const auto fcs = faces( msh, cl );
    if ( fcs.size() != 6 ) {
        return false;
    }

    const auto bT = barycenter( msh, cl );

    for ( auto &fc : fcs ) {
        const auto n = normal( msh, cl, fc );
        const auto bF = barycenter( msh, fc );

        const auto vFT = ( bF - bT ).to_vector().normalized();

        if ( n.dot( vFT ) < ( 1. - 1e-8 ) ) {
            return false;
        }
    }

    return true;
}

template < typename T >
bool is_extruted_hexa( const std::array< point< T, 3 >, 8 > &pts ) {

    const T thrs = 1e-8;

    const auto n = pts[4] - pts[0];
    const T n_norm = n.to_vector().norm();

    for ( int i = 0; i < 3; i++ ) {
        const auto p_test = pts[i] + n;
        const auto dist = ( pts[i + 4] - p_test ).to_vector().norm();

        if ( dist > thrs * n_norm ) {
            return false;
        }
    }

    return true;
}

template < typename T >
std::array< point< T, 3 >, 8 >
renumber_hexa( const disk::generic_mesh< T, 3 > &msh,
               const typename disk::generic_mesh< T, 3 >::cell_type &cl ) {

    using point_id = point_identifier< 3 >;

    const auto pts = points( msh, cl );
    const auto pts_ids = cl.point_ids();

    std::map< point_id, point< T, 3 > > map_cl;

    int i = 0;
    for ( auto &pt : pts ) {
        map_cl[pts_ids[i++]] = pt;
    }

    const auto fcs = faces( msh, cl );
    assert( fcs.size() == 6 );

    std::map< point_id, std::vector< point_id > > map_fc0;
    auto pts_fc0 = fcs[0].point_ids();
    auto pts2_fc0 = points( msh, fcs[0] );

    const auto n0 = normal( msh, cl, fcs[0] );

    auto v0 = ( pts2_fc0[1] - pts2_fc0[0] ).to_vector();
    auto v1 = ( pts2_fc0[2] - pts2_fc0[1] ).to_vector();
    auto n = v0.cross( v1 );

    if ( n.dot( n0 ) <= 0 ) {
        point_id pt = pts_fc0[1];
        pts_fc0[1] = pts_fc0[3];
        pts_fc0[3] = pt;
    }

    for ( auto &pt : pts_fc0 ) {
        map_fc0[pt] = std::vector< point_id >();
    }

    for ( int i = 1; i < fcs.size(); i++ ) {
        const auto fc = fcs[i];
        const auto pts_fc = fc.point_ids();

        for ( auto &pt : pts_fc ) {
            if ( map_fc0.contains( pt ) ) {
                for ( auto &pt2 : pts_fc ) {
                    if ( !map_fc0.contains( pt2 ) ) {
                        map_fc0[pt].push_back( pt2 );
                    }
                }
            }
        }
    }

    std::array< point< T, 3 >, 8 > new_num;
    int j = 0;
    for ( auto &pt : pts_fc0 ) {
        const auto vect = map_fc0[pt];
        point_id pt_corr = -1;
        for ( auto &val : vect ) {
            const auto nb_elem = std::count( vect.begin(), vect.end(), val );
            if ( nb_elem == 2 ) {
                pt_corr = val;
                break;
            }
        }
        if ( pt_corr < 0 ) {
            throw std::runtime_error( "Error" );
        };

        new_num[j] = map_cl[pt];
        new_num[j + 4] = map_cl[pt_corr];
        j++;
    }

    // std::cout << "(" << new_num[0] << ", " << new_num[1] << ", " << new_num[2] << ", " <<
    // new_num[3]
    //           << ", " << new_num[4] << ", " << new_num[5] << ", " << new_num[6] << ", "
    //           << new_num[7] << ")" << std::endl;

    return new_num;
}

} // namespace priv

} // namespace quadrature

template<typename T>
std::vector<disk::quadrature_point<T, 2>>
integrate(const disk::generic_mesh<T, 2>& msh, const typename disk::generic_mesh<T, 2>::cell& cl, size_t degree)
{
    const auto pts = points(msh, cl);

    assert((pts.size() > 2) && "Insufficient points for a 2D cell");

    if (pts.size() == 3)
    {
        return disk::quadrature::triangle_gauss(degree, pts[0], pts[1], pts[2]);
    }

    /* The coordinate transformation could become non-linear, so use tensorized Gauss
     * points only on quadrilaterals which are almost square. */
    if ( pts.size() == 4 and quadrature::priv::is_ortho_quad( msh, cl ) ) {
        return disk::quadrature::tensorized_gauss_legendre(degree, pts[0], pts[1], pts[2], pts[3]);
    }

    if ( is_convex( msh, cl ) ) {
        return quadrature::priv::integrate_convex(msh, cl, degree);
    }

    return quadrature::priv::integrate_nonconvex(msh, cl, degree);
}

template<typename T>
auto
integrate(const disk::generic_mesh<T, 2>& msh, const typename disk::generic_mesh<T, 2>::face& fc, size_t degree)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 2);
    return disk::quadrature::gauss_legendre(degree, pts[0], pts[1]);
}

namespace priv
{

template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate_polyhedron(const disk::generic_mesh<T, 3>&                msh,
                     const typename disk::generic_mesh<T, 3>::cell& cl,
                     const size_t                                   degree)
{
    using quadpoint_type = disk::quadrature_point<T, 3>;

    const auto rss = split_in_raw_tetrahedra(msh, cl);

    std::vector< quadpoint_type > ret;
    for (auto& rs : rss)
    {
        const auto pts = rs.points();
        assert(pts.size() == 4);
        const auto quad_tet = disk::quadrature::grundmann_moeller(degree, pts[0], pts[1], pts[2], pts[3]);
        ret.insert(ret.end(), quad_tet.begin(), quad_tet.end());
    }

    return ret;
}

template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate_polyhedron_face(const disk::generic_mesh<T, 3>&                msh,
                          const typename disk::generic_mesh<T, 3>::face& fc,
                          const size_t                                   degree)
{
    using quadpoint_type = disk::quadrature_point<T, 3>;

    const auto rss = split_in_raw_triangles(msh, fc);

    std::vector<quadpoint_type> ret;
    ret.reserve(quadrature::priv::triangle_gauss_rules[degree].num_points * rss.size());
    for (auto& rs : rss)
    {
        const auto pts = rs.points();
        assert(pts.size() == 3);

        const auto quad_tri = disk::quadrature::triangle_gauss(degree, pts[0], pts[1], pts[2]);
        ret.insert(ret.end(), quad_tri.begin(), quad_tri.end());
    }

    return ret;
}

template < typename T >
std::vector< disk::quadrature_point< T, 3 > >
integrate_hexahedron_extruded( const std::array< point< T, 3 >, 8 > &pts, const size_t degree ) {
    using quadpoint_type = disk::quadrature_point< T, 3 >;

    //  compute quadrature on the basis
    std::vector< quadpoint_type > quad_basis;
    quad_basis.reserve( quadrature::priv::triangle_gauss_rules[degree].num_points * 2 );

    const auto quad_tri0 = disk::quadrature::triangle_gauss( degree, pts[0], pts[1], pts[2] );
    quad_basis.insert( quad_basis.end(), quad_tri0.begin(), quad_tri0.end() );

    const auto quad_tri1 = disk::quadrature::triangle_gauss( degree, pts[2], pts[3], pts[0] );
    quad_basis.insert( quad_basis.end(), quad_tri1.begin(), quad_tri1.end() );

    // compute quad on axis
    const auto axis = pts[4] - pts[0];
    const T dist = distance( pts[0], pts[4] );

    const auto quad_axis = disk::quadrature::gauss_legendre( degree, 0.0, 1.0 );

    std::vector< quadpoint_type > quad_extr;
    quad_extr.reserve( quad_axis.size() * quad_basis.size() );

    for ( auto &qp_a : quad_axis ) {
        const auto weight = dist * qp_a.weight();
        const auto a = qp_a.point().x() * axis;
        for ( auto &qp_b : quad_basis ) {
            quad_extr.push_back( { qp_b.point() + a, qp_b.weight() * weight } );
        }
    }

    return quad_extr;
}

template < typename T >
std::vector< disk::quadrature_point< T, 3 > >
integrate_hexahedron( const disk::generic_mesh< T, 3 > &msh,
                      const typename disk::generic_mesh< T, 3 >::cell &cl, const size_t degree ) {
    using quadpoint_type = disk::quadrature_point< T, 3 >;
    using raw_simplex_type = raw_simplex< typename disk::generic_mesh< T, 3 >::point_type, 3 >;

    const auto pts = quadrature::priv::renumber_hexa( msh, cl );

    if ( quadrature::priv::is_ortho_hexa( msh, cl ) ) {
        return disk::quadrature::tensorized_gauss_legendre( degree, pts );
    } else if ( quadrature::priv::is_extruted_hexa( pts ) ) {
        return priv::integrate_hexahedron_extruded( pts, degree );
    }

    std::vector< raw_simplex_type > rss;
    rss.push_back( raw_simplex_type( { pts[0], pts[1], pts[3], pts[4] } ) );
    rss.push_back( raw_simplex_type( { pts[1], pts[2], pts[3], pts[6] } ) );
    rss.push_back( raw_simplex_type( { pts[1], pts[3], pts[4], pts[6] } ) );
    rss.push_back( raw_simplex_type( { pts[1], pts[4], pts[5], pts[6] } ) );
    rss.push_back( raw_simplex_type( { pts[3], pts[4], pts[6], pts[7] } ) );

    std::vector< quadpoint_type > ret;
    for ( auto &rs : rss ) {
        const auto pts_rs = rs.points();
        assert( pts_rs.size() == 4 );
        const auto quad_tet = disk::quadrature::grundmann_moeller( degree, pts_rs[0], pts_rs[1],
                                                                   pts_rs[2], pts_rs[3] );
        ret.insert( ret.end(), quad_tet.begin(), quad_tet.end() );
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

    if ( pts.size() < 4 ) {
        throw std::invalid_argument( "A 3D cell cannot have less than four points. "
                                     "This looks like a nice bug." );
    }

    if ( quadrature::priv::is_hexa( msh, cl ) ) {
        return priv::integrate_hexahedron( msh, cl, degree );
    }

    return priv::integrate_polyhedron( msh, cl, degree );
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
integrate(const disk::generic_mesh<T, 1>& msh, const typename generic_mesh<T, 1>::cell& cl, size_t degree)
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
