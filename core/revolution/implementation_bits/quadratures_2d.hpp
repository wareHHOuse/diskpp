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

#include "geometry/geometry.hpp"
#include "quadratures/quad_bones.hpp"

namespace revolution {
namespace proton {

template<typename T>
std::vector<std::pair<point<T,1>, T>>
golub_welsch(size_t degree)
{
    std::vector<std::pair<point<T,1>, T>> ret;

    using namespace Eigen;

    if ( degree % 2 == 0)
        degree++;

    size_t reqd_nodes = (degree+1)/2;

    if (reqd_nodes == 1)
    {
        auto qp = std::make_pair(point<T,1>({0.}), 2.0);
        ret.push_back(qp);
        return ret;
    }

    Matrix<T, Dynamic, Dynamic> M = Matrix<T, Dynamic, Dynamic>::Zero(reqd_nodes, reqd_nodes);
    for (size_t i = 1; i < reqd_nodes; i++)
    {
        T p = 4.0 - 1.0/(i*i);
        M(i, i-1) = sqrt(1.0/p);
    }

    SelfAdjointEigenSolver<Matrix<T, Dynamic, Dynamic>> solver;
    solver.compute(M);

    Matrix<T, Dynamic, Dynamic> weights = solver.eigenvectors().row(0);
                                weights = weights.array().square();
    Matrix<T, Dynamic, Dynamic> nodes = solver.eigenvalues();

    assert( weights.size() == nodes.size() );

    for (size_t i = 0; i < nodes.size(); i++)
    {
        auto qp = std::make_pair(point<T,1>({nodes(i)}), 2*weights(i));
        ret.push_back(qp);
    }

    return ret;
}


template<typename T>
std::vector<std::pair<point<T,1>, T>>
gauss_legendre(size_t degree)
{
    auto comp_degree = degree;

    if ( degree % 2 == 0 )
        comp_degree = degree + 1;

    size_t reqd_nodes = (comp_degree+1)/2;

    std::vector< std::pair<point<T,1>, T> > ret;
    ret.reserve(reqd_nodes);

    point<T,1>  qp;
    T           qw;
    T           a1, a2;

    switch(reqd_nodes)
    {
        case 1:
            qp = point<T,1>({0.0});
            qw = 2.0;
            ret.push_back( std::make_pair(qp, qw) );
            return ret;

        case 2:
            qp = point<T,1>({ 1.0/std::sqrt(3.0) });
            qw = 1.0;
            ret.push_back( std::make_pair(-qp, qw ) );
            ret.push_back( std::make_pair( qp, qw ) );
            return ret;

        case 3:
            qp = point<T,1>({ std::sqrt(3.0/5.0) });
            qw = 5.0/9.0;
            ret.push_back( std::make_pair(-qp, qw ) );
            ret.push_back( std::make_pair( qp, qw ) );

            qp = point<T,1>({0.0});
            qw = 8.0/9.0;
            ret.push_back( std::make_pair( qp, qw ) );

            return ret;

        case 4:
            a1 = 15;
            a2 = 2.0*std::sqrt(30.0);

            qp = point<T,1>({ std::sqrt( (a1 - a2)/35.0 ) });
            qw = (18.0 + std::sqrt(30.0))/36.0;
            ret.push_back( std::make_pair(-qp, qw ) );
            ret.push_back( std::make_pair( qp, qw ) );

            qp = point<T,1>({ std::sqrt( (a1 + a2)/35.0 ) });
            qw = (18.0 - std::sqrt(30.0))/36.0;
            ret.push_back( std::make_pair(-qp, qw ) );
            ret.push_back( std::make_pair( qp, qw ) );

            return ret;

        case 5:
            qp = point<T,1>({ 0.0 });
            qw = 128.0/225.0;
            ret.push_back( std::make_pair( qp, qw ) );

            a1 = 5.0;
            a2 = 2.0*std::sqrt(10.0/7.0);
            qp = point<T,1>({ std::sqrt(a1 - a2)/3.0 });
            qw = (322 + 13.0*std::sqrt(70.0))/900.0;
            ret.push_back( std::make_pair(-qp, qw ) );
            ret.push_back( std::make_pair( qp, qw ) );

            qp = point<T,1>({ std::sqrt(a1 + a2)/3.0 });
            qw = (322 - 13.0*std::sqrt(70.0))/900.0;
            ret.push_back( std::make_pair(-qp, qw ) );
            ret.push_back( std::make_pair( qp, qw ) );
            return ret;

        default:
            return golub_welsch<T>(degree);

    }

    return ret;
}

template<typename T>
std::vector<std::pair<point<T,1>, T>>
edge_quadrature(size_t doe)
{
    return gauss_legendre<T>(doe);
}

template<typename T>
std::vector<std::pair<point<T, 2>, T>>
quadrangle_quadrature(size_t degree)
{
    auto qps = proton::edge_quadrature<T>(degree);

    std::vector<std::pair<point<T, 2>, T>> ret; //;
    ret.reserve(qps.size() * qps.size());

    for (auto jtor = qps.begin(); jtor != qps.end(); jtor++)
    {
        auto qp_y = *jtor;
        auto eta  = qp_y.first.x();

        for (auto itor = qps.begin(); itor != qps.end(); itor++)
        {
            auto qp_x = *itor;
            auto xi  = qp_x.first.x();

            auto w = qp_x.second * qp_y.second;

            ret.push_back(std::make_pair(point<T, 2>({xi, eta}), w));
        }
    }

    return ret;
}
} //proton

namespace priv {

template<typename T, typename PtA>
std::vector< disk::quadrature_point<T,2> >
integrate_triangle(size_t degree, const PtA& pts)
{
	assert(pts.size() == 3);

    static_assert( std::is_same<typename PtA::value_type, point<T,2>>::value,
                   "This function is for 2D points" );

	using quadpoint_type = disk::quadrature_point<T,2>;

    auto qps = disk::triangle_quadrature(degree);

        std::vector<quadpoint_type> ret;

        ret.resize(qps.size());

        auto col1 = pts[1] - pts[0];
        auto col2 = pts[2] - pts[0];

        /* Compute the area of the sub-triangle */
        auto tm = (col1.x() * col2.y() - col2.x() * col1.y()) / 2.;

        auto tr = [&](const std::pair<point<T, 2>, T>& qd) -> auto
        {
            auto point  = col1 * qd.first.x() + col2 * qd.first.y() + pts[0];
            auto weight = qd.second * std::abs(tm);
            return disk::make_qp(point, weight);
    };

    std::transform(qps.begin(), qps.end(), ret.begin(), tr);

    return ret;
}



template<typename T>
std::vector< disk::quadrature_point<T,2> >
integrate_quadrangle_tens(size_t degree, const std::vector< point<T,2> >& pts)
{
    auto qps = proton::quadrangle_quadrature<T>(degree);

    std::vector<disk::quadrature_point<T, 2>> ret;
    ret.reserve(qps.size());

    auto P = [&](T xi, T eta) -> T {
        return 0.25 * pts[0].x() * (1-xi)*(1-eta) +
               0.25 * pts[1].x() * (1+xi)*(1-eta) +
               0.25 * pts[2].x() * (1+xi)*(1+eta) +
               0.25 * pts[3].x() * (1-xi)*(1+eta);
    };

    auto Q = [&](T xi, T eta) -> T {
        return 0.25 * pts[0].y() * (1-xi)*(1-eta) +
               0.25 * pts[1].y() * (1+xi)*(1-eta) +
               0.25 * pts[2].y() * (1+xi)*(1+eta) +
               0.25 * pts[3].y() * (1-xi)*(1+eta);
    };

    auto J = [&](T xi, T eta) -> T {
        auto j11 = 0.25*((pts[1].x() - pts[0].x())*(1-eta) + (pts[2].x() - pts[3].x())*(1+eta));
        auto j12 = 0.25*((pts[1].y() - pts[0].y())*(1-eta) + (pts[2].y() - pts[3].y())*(1+eta));
        auto j21 = 0.25*((pts[3].x() - pts[0].x())*(1-xi) + (pts[2].x() - pts[1].x())*(1+xi));
        auto j22 = 0.25*((pts[3].y() - pts[0].y())*(1-xi) + (pts[2].y() - pts[1].y())*(1+xi));

        return std::abs(j11*j22 - j12*j21);
    };

    for(auto& qp : qps)
    {
            T xi  = qp.first.x();
            T eta = qp.first.y();

            T px = P(xi, eta);
            T py = Q(xi, eta);

            T weight = qp.second * J(xi, eta);
            ret.push_back(disk::make_qp(point<T, 2>({px, py}), weight));
    }

    return ret;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename PtA>
std::vector< disk::quadrature_point<T,2> >
integrate_quadrangle(const Mesh<T,2,Storage>& msh,
				     const typename Mesh<T,2,Storage>::cell& cl,
				     size_t degree,
                     const PtA& pts)
{
	assert(pts.size() == 4);

	using mesh_type = Mesh<T,2,Storage>;
	using quadpoint_type = disk::quadrature_point<T,2>;

	auto qps = disk::triangle_quadrature(degree);

    std::vector<quadpoint_type> ret;

    ret.resize( qps.size() * 2 );

	for (size_t i = 1; i < 3; i++)
    {
        auto pt1 = pts[i];
        auto pt2 = pts[i+1];
        auto col1 = pt1 - pts[0];
        auto col2 = pt2 - pts[0];

        /* Compute the area of the sub-triangle */
        auto tm = (col1.x()*col2.y() - col2.x()*col1.y())/2.;

        auto tr = [&](const std::pair<point<T,2>, T>& qd) -> auto {
            auto point = col1*qd.first.x() + col2*qd.first.y() + pts[0];
            auto weight = qd.second * std::abs(tm);
            return disk::make_qp(point, weight);
        };

        auto retbegin = ret.begin();
        std::advance(retbegin, qps.size()*(i-1));

        std::transform(qps.begin(), qps.end(), retbegin, tr);
    }

    return ret;
}

namespace priv {

template<typename Iterator>
auto
barycenter(Iterator begin, Iterator end)
{
    typedef typename std::iterator_traits<Iterator>::value_type point_type;
    typedef typename point_type::value_type T;

    point_type  ret;
    T           den = 0.0;

    auto numpts = std::distance(begin, end);
    auto p0 = *begin;

    for (size_t i = 2; i < numpts; i++)
    {
        auto pprev  = *std::next(begin, i-1) - p0;
        auto pcur   = *std::next(begin, i) - p0;
        auto d      = det(pprev, pcur) / 2.0;
        ret         = ret + (pprev + pcur) * d;
        den         += d;
    }

    return p0 + ret/(den*3);
}
} // priv

//#define OPTIMAL_TRIANGLE_NUMBER

    /* The 'optimal triangle number' version gives almost the same results
     * of the other version. In bigger meshes there are some advantages in
     * assembly time. The problem with this version is that it could generate
     * triangles with a very bad aspect ratio and at the moment I don't know
     * if and how this can affect computations, so I leave it turned off.
     */
template<typename T>
std::vector< disk::quadrature_point<T,2> >
integrate_polygon(size_t degree, const std::vector< point<T,2> >& pts)
{
	using quadpoint_type = disk::quadrature_point<T,2>;

	auto qps = disk::triangle_quadrature(degree);

    std::vector<quadpoint_type> ret;

#ifdef OPTIMAL_TRIANGLE_NUMBER
    /* Break the cell in triangles, compute the transformation matrix and
     * map quadrature data in the physical space. Edges of the triangle as
     * column vectors in the transformation matrix. */
    ret.resize( qps.size() * pts.size()-2 );
    for (size_t i = 1; i < pts.size()-1; i++)
    {
        auto pt1 = pts[i];
        auto pt2 = pts[i+1];
        auto col1 = pt1 - pts[0];
        auto col2 = pt2 - pts[0];

        /* Compute the area of the sub-triangle */
        auto tm = (col1.x()*col2.y() - col2.x()*col1.y())/2.;

        auto tr = [&](const std::pair<point<T,2>, T>& qd) -> auto {
            auto point = col1*qd.first.x() + col2*qd.first.y() + pts[0];
            auto weight = qd.second * std::abs(tm);
            return disk::make_qp(point, weight);
        };

        auto retbegin = ret.begin();
        std::advance(retbegin, qps.size()*(i-1));

        std::transform(qps.begin(), qps.end(), retbegin, tr);
    }
#else
    auto c_center   = priv::barycenter(pts.begin(), pts.end());

    /* Break the cell in triangles, compute the transformation matrix and
     * map quadrature data in the physical space. Edges of the triangle as
     * column vectors in the transformation matrix. */
    ret.resize( qps.size() * pts.size() );
    for (size_t i = 0; i < pts.size(); i++)
    {
        auto pt1 = pts[i];
        auto pt2 = pts[(i+1)%pts.size()];
        auto col1 = pt1 - c_center;
        auto col2 = pt2 - c_center;

        /* Compute the area of the sub-triangle */
        auto tm = (col1.x()*col2.y() - col2.x()*col1.y())/2.;

        auto tr = [&](const std::pair<point<T,2>, T>& qd) -> auto {
            auto point = col1*qd.first.x() + col2*qd.first.y() + c_center;
            auto weight = qd.second * std::abs(tm);
            return disk::make_qp(point, weight);
        };

        auto retbegin = ret.begin();
        std::advance(retbegin, qps.size()*i);

        std::transform(qps.begin(), qps.end(), retbegin, tr);
    }
#endif
    return ret;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
std::vector< disk::quadrature_point<T,2> >
integrate_2D_face(const Mesh<T,2,Storage>& msh,
		          const typename Mesh<T,2,Storage>::face& fc,
		          size_t degree)
{
	using mesh_type 	= disk::simplicial_mesh<T,2>;
	using point_type 	= typename mesh_type::point_type;

    auto qps = proton::edge_quadrature<T>(degree);
    auto pts = points(msh, fc);

    auto scale = (pts[1] - pts[0]);
    auto meas = scale.to_vector().norm();

    std::vector< disk::quadrature_point<T,2> >   ret;
    ret.reserve(qps.size());

    for (auto itor = qps.begin(); itor != qps.end(); itor++)
    {
        auto qp = *itor;
        auto t = qp.first.x();
        auto p = 0.5*(1-t)*pts[0] + 0.5*(1+t)*pts[1];
        auto w = qp.second * meas * 0.5;

        ret.push_back( disk::make_qp(p, w) );
    }

    return ret;
}

template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate_hexahedron_tens(int degree, const std::array<point<T, 3>, 8>& pts)
{
    assert(degree >= 0);
    auto qps = proton::edge_quadrature<T>(degree);

    std::vector<disk::quadrature_point<T, 3>> ret;

    auto P = [&](T xi, T eta, T zeta) -> T {
        return 0.125 * pts[0].x() * (1 - xi) * (1 - eta) * (1 - zeta) +
               0.125 * pts[1].x() * (1 + xi) * (1 - eta) * (1 - zeta) +
               0.125 * pts[2].x() * (1 + xi) * (1 + eta) * (1 - zeta) +
               0.125 * pts[3].x() * (1 - xi) * (1 + eta) * (1 - zeta) +
               0.125 * pts[4].x() * (1 - xi) * (1 - eta) * (1 + zeta) +
               0.125 * pts[5].x() * (1 + xi) * (1 - eta) * (1 + zeta) +
               0.125 * pts[6].x() * (1 + xi) * (1 + eta) * (1 + zeta) +
               0.125 * pts[7].x() * (1 - xi) * (1 + eta) * (1 + zeta);
    };

    auto Q = [&](T xi, T eta, T zeta) -> T {
        return 0.125 * pts[0].y() * (1 - xi) * (1 - eta) * (1 - zeta) +
               0.125 * pts[1].y() * (1 + xi) * (1 - eta) * (1 - zeta) +
               0.125 * pts[2].y() * (1 + xi) * (1 + eta) * (1 - zeta) +
               0.125 * pts[3].y() * (1 - xi) * (1 + eta) * (1 - zeta) +
               0.125 * pts[4].y() * (1 - xi) * (1 - eta) * (1 + zeta) +
               0.125 * pts[5].y() * (1 + xi) * (1 - eta) * (1 + zeta) +
               0.125 * pts[6].y() * (1 + xi) * (1 + eta) * (1 + zeta) +
               0.125 * pts[7].y() * (1 - xi) * (1 + eta) * (1 + zeta);
    };

    auto R = [&](T xi, T eta, T zeta) -> T {
        return 0.125 * pts[0].z() * (1 - xi) * (1 - eta) * (1 - zeta) +
               0.125 * pts[1].z() * (1 + xi) * (1 - eta) * (1 - zeta) +
               0.125 * pts[2].z() * (1 + xi) * (1 + eta) * (1 - zeta) +
               0.125 * pts[3].z() * (1 - xi) * (1 + eta) * (1 - zeta) +
               0.125 * pts[4].z() * (1 - xi) * (1 - eta) * (1 + zeta) +
               0.125 * pts[5].z() * (1 + xi) * (1 - eta) * (1 + zeta) +
               0.125 * pts[6].z() * (1 + xi) * (1 + eta) * (1 + zeta) +
               0.125 * pts[7].z() * (1 - xi) * (1 + eta) * (1 + zeta);
    };

    auto J = [&](T xi, T eta, T zeta) -> T {
        static_matrix<T, 3, 3> Jac = static_matrix<T, 3, 3>::Zero();
        Jac(0, 0) = -0.125 * pts[0].x() * (1 - eta) * (1 - zeta) + 0.125 * pts[1].x() * (1 - eta) * (1 - zeta) +
                    0.125 * pts[2].x() * (1 + eta) * (1 - zeta) - 0.125 * pts[3].x() * (1 + eta) * (1 - zeta) -
                    0.125 * pts[4].x() * (1 - eta) * (1 + zeta) + 0.125 * pts[5].x() * (1 - eta) * (1 + zeta) +
                    0.125 * pts[6].x() * (1 + eta) * (1 + zeta) - 0.125 * pts[7].x() * (1 + eta) * (1 + zeta);

        Jac(0, 1) = -0.125 * pts[0].y() * (1 - eta) * (1 - zeta) + 0.125 * pts[1].y() * (1 - eta) * (1 - zeta) +
                    0.125 * pts[2].y() * (1 + eta) * (1 - zeta) - 0.125 * pts[3].y() * (1 + eta) * (1 - zeta) -
                    0.125 * pts[4].y() * (1 - eta) * (1 + zeta) + 0.125 * pts[5].y() * (1 - eta) * (1 + zeta) +
                    0.125 * pts[6].y() * (1 + eta) * (1 + zeta) - 0.125 * pts[7].y() * (1 + eta) * (1 + zeta);

        Jac(0, 2) = -0.125 * pts[0].z() * (1 - eta) * (1 - zeta) + 0.125 * pts[1].z() * (1 - eta) * (1 - zeta) +
                    0.125 * pts[2].z() * (1 + eta) * (1 - zeta) - 0.125 * pts[3].z() * (1 + eta) * (1 - zeta) -
                    0.125 * pts[4].z() * (1 - eta) * (1 + zeta) + 0.125 * pts[5].z() * (1 - eta) * (1 + zeta) +
                    0.125 * pts[6].z() * (1 + eta) * (1 + zeta) - 0.125 * pts[7].z() * (1 + eta) * (1 + zeta);

        Jac(1, 0) = -0.125 * pts[0].x() * (1 - xi) * (1 - zeta) + 0.125 * pts[1].x() * (1 + xi) * (1 - zeta) +
                    0.125 * pts[2].x() * (1 + xi) * (1 - zeta) - 0.125 * pts[3].x() * (1 - xi) * (1 - zeta) -
                    0.125 * pts[4].x() * (1 - xi) * (1 + zeta) + 0.125 * pts[5].x() * (1 + xi) * (1 + zeta) +
                    0.125 * pts[6].x() * (1 + xi) * (1 + zeta) - 0.125 * pts[7].x() * (1 - xi) * (1 + zeta);

        Jac(1, 1) = -0.125 * pts[0].y() * (1 - xi) * (1 - zeta) + 0.125 * pts[1].y() * (1 + xi) * (1 - zeta) +
                    0.125 * pts[2].y() * (1 + xi) * (1 - zeta) - 0.125 * pts[3].y() * (1 - xi) * (1 - zeta) -
                    0.125 * pts[4].y() * (1 - xi) * (1 + zeta) + 0.125 * pts[5].y() * (1 + xi) * (1 + zeta) +
                    0.125 * pts[6].y() * (1 + xi) * (1 + zeta) - 0.125 * pts[7].y() * (1 - xi) * (1 + zeta);

        Jac(1, 2) = -0.125 * pts[0].z() * (1 - xi) * (1 - zeta) + 0.125 * pts[1].z() * (1 + xi) * (1 - zeta) +
                    0.125 * pts[2].z() * (1 + xi) * (1 - zeta) - 0.125 * pts[3].z() * (1 - xi) * (1 - zeta) -
                    0.125 * pts[4].z() * (1 - xi) * (1 + zeta) + 0.125 * pts[5].z() * (1 + xi) * (1 + zeta) +
                    0.125 * pts[6].z() * (1 + xi) * (1 + zeta) - 0.125 * pts[7].z() * (1 - xi) * (1 + zeta);

        Jac(2, 0) = -0.125 * pts[0].x() * (1 - xi) * (1 - eta) - 0.125 * pts[1].x() * (1 + xi) * (1 - eta) -
                    0.125 * pts[2].x() * (1 + xi) * (1 + eta) - 0.125 * pts[3].x() * (1 - xi) * (1 + eta) -
                    0.125 * pts[4].x() * (1 - xi) * (1 - eta) + 0.125 * pts[5].x() * (1 + xi) * (1 - eta) +
                    0.125 * pts[6].x() * (1 + xi) * (1 + eta) + 0.125 * pts[7].x() * (1 - xi) * (1 + eta);

        Jac(2, 1) = -0.125 * pts[0].y() * (1 - xi) * (1 - eta) - 0.125 * pts[1].y() * (1 + xi) * (1 - eta) -
                    0.125 * pts[2].y() * (1 + xi) * (1 + eta) - 0.125 * pts[3].y() * (1 - xi) * (1 + eta) -
                    0.125 * pts[4].y() * (1 - xi) * (1 - eta) + 0.125 * pts[5].y() * (1 + xi) * (1 - eta) +
                    0.125 * pts[6].y() * (1 + xi) * (1 + eta) + 0.125 * pts[7].y() * (1 - xi) * (1 + eta);

        Jac(2, 2) = -0.125 * pts[0].z() * (1 - xi) * (1 - eta) - 0.125 * pts[1].z() * (1 + xi) * (1 - eta) -
                    0.125 * pts[2].z() * (1 + xi) * (1 + eta) - 0.125 * pts[3].z() * (1 - xi) * (1 + eta) -
                    0.125 * pts[4].z() * (1 - xi) * (1 - eta) + 0.125 * pts[5].z() * (1 + xi) * (1 - eta) +
                    0.125 * pts[6].z() * (1 + xi) * (1 + eta) + 0.125 * pts[7].z() * (1 - xi) * (1 + eta);

        return std::abs(Jac.determinant());
    };

    for (auto jtor = qps.begin(); jtor != qps.end(); jtor++)
    {
        for (auto itor = qps.begin(); itor != qps.end(); itor++)
        {
            for (auto ktor = qps.begin(); ktor != qps.end(); ktor++)
            {
                auto qp_x = *itor;
                auto qp_y = *jtor;
                auto qp_z = *ktor;

                auto xi   = qp_x.first.x();
                auto eta  = qp_y.first.x();
                auto zeta = qp_z.first.x();

                auto px = P(xi, eta, zeta);
                auto py = Q(xi, eta, zeta);
                auto pz = R(xi, eta, zeta);

                auto w = qp_x.second * qp_y.second * qp_z.second * J(xi, eta, zeta);

                ret.push_back(disk::make_qp(point<T, 3>({px, py, pz}), w));
            }
        }
    }

    return ret;
}

template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate_polyhedron(const disk::generic_mesh<T, 3>&                msh,
                     const typename disk::generic_mesh<T, 3>::cell& cl,
                     size_t                                         degree)
{
    using quadpoint_type = disk::quadrature_point<T, 3>;

        const auto m_quadrature_data = disk::tetrahedron_quadrature(degree);

        auto rss = split_in_raw_tetrahedra(msh, cl);

        std::vector<quadpoint_type> ret;
        ret.reserve(m_quadrature_data.size() * rss.size());
        for (auto& rs : rss)
        {
            auto pts = rs.points();
            assert(pts.size() == 4);

            auto meas = measure(rs);
            auto col1 = pts[1] - pts[0];
            auto col2 = pts[2] - pts[0];
            auto col3 = pts[3] - pts[0];

            for (auto& qd : m_quadrature_data)
            {
                auto point  = col1 * qd.first.x() + col2 * qd.first.y() + col3 * qd.first.z() + pts[0];
                auto weight = qd.second * meas;
                ret.push_back(disk::make_qp(point, weight));
            }
    }

    return ret;
}

template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate_polyhedron_face(const disk::generic_mesh<T, 3>&                msh,
                     const typename disk::generic_mesh<T, 3>::face& fc,
                     size_t                                         degree)
{
    using quadpoint_type = disk::quadrature_point<T, 3>;

    const auto m_quadrature_data = disk::triangle_quadrature(degree);

    auto rss = split_in_raw_triangles(msh, fc);

    std::vector<quadpoint_type> ret;
    ret.reserve(m_quadrature_data.size() * rss.size());
    for (auto& rs : rss)
    {
        auto pts = rs.points();
        assert(pts.size() == 3);
        auto meas = measure(rs);

        auto col1 = pts[1] - pts[0];
        auto col2 = pts[2] - pts[0];

        for (auto& qd : m_quadrature_data)
        {
            auto point  = col1 * qd.first.x() + col2 * qd.first.y() + pts[0];
            auto weight = qd.second * meas;
            ret.push_back(disk::make_qp(point, weight));
        }
    }

    return ret;
}

template<typename MeshType, typename Element>
std::vector<disk::quadrature_point<typename MeshType::coordinate_type, MeshType::dimension>>
integrate_degree0(const MeshType& msh, const Element& elem)
{
    std::vector<disk::quadrature_point<typename MeshType::coordinate_type, MeshType::dimension>> ret;

    const auto bar  = barycenter(msh, elem);
    const auto meas = measure(msh, elem);

    ret.push_back(disk::make_qp(bar, meas));

    return ret;
}

} //priv


template<typename T>
std::vector< disk::quadrature_point<T,2> >
integrate(const disk::simplicial_mesh<T,2>& msh,
		  const typename disk::simplicial_mesh<T,2>::cell& cl,
		  size_t degree)
{
    if(degree == 0)
    {
        return priv::integrate_degree0(msh, cl);
    }

    auto pts = points(msh, cl);
    return priv::integrate_triangle<T>(degree, pts);
}


template<typename T>
std::vector< disk::quadrature_point<T,2> >
integrate(const disk::simplicial_mesh<T,2>& msh,
		  const typename disk::simplicial_mesh<T,2>::face& fc,
		  size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, fc);
    }

    return priv::integrate_2D_face(msh, fc, degree);
}

template<typename T>
std::vector< disk::quadrature_point<T,2> >
integrate(const disk::cartesian_mesh<T,2>& msh,
		  const typename disk::cartesian_mesh<T,2>::cell& cl,
		  size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, cl);
    }

    auto pts = points(msh, cl);
    // transform arry to vector
    std::vector< point<T,2> > ptsv{pts[0], pts[1], pts[3], pts[2]};
    return priv::integrate_quadrangle_tens(degree, ptsv);
}



template<typename T>
std::vector< disk::quadrature_point<T,2> >
integrate(const disk::cartesian_mesh<T,2>& msh,
		  const typename disk::cartesian_mesh<T,2>::face& fc,
		  size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, fc);
    }

    return priv::integrate_2D_face(msh, fc, degree);
}

template<typename T>
std::vector< disk::quadrature_point<T,2> >
integrate(const disk::generic_mesh<T,2>& msh,
		  const typename disk::generic_mesh<T,2>::cell& cl,
		  size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, cl);
    }

    auto pts = points(msh, cl);
    switch (pts.size())
    {
        case 0:
        case 1:
        case 2:
            throw std::invalid_argument("A 2D cell cannot have less than three points. "
                                        "This looks like a nice bug.");

        case 3: return priv::integrate_triangle<T>(degree, pts);

        case 4: return priv::integrate_quadrangle_tens(degree, pts);

        default: return priv::integrate_polygon(degree, pts);
	}
}

template<typename T>
std::vector< disk::quadrature_point<T,2> >
integrate(const disk::generic_mesh<T,2>& msh,
		  const typename disk::generic_mesh<T,2>::face& fc,
		  size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, fc);
    }

    return priv::integrate_2D_face(msh, fc, degree);
}


namespace priv {

template<typename T>
point<T,3>
map_to_physical(const disk::simplicial_mesh<T,3>& msh,
                const typename disk::simplicial_mesh<T,3>::face& face,
                const point<T,2>& pm)
{
    auto pts = points(msh, face);
    return pts[0] + (pts[1]-pts[0])*pm.x() + (pts[2]-pts[0])*pm.y();
}

template<typename T>
point<T,3>
map_to_physical(const disk::simplicial_mesh<T,3>& msh,
                const typename disk::simplicial_mesh<T,3>::cell& cell,
                const point<T,3>& p)
{
    auto pts = points(msh, cell);

    auto pp = (pts[1] - pts[0]) * p.x() +
              (pts[2] - pts[0]) * p.y() +
              (pts[3] - pts[0]) * p.z() +
              pts[0];

    return pp;
}



} // namespace priv

template<typename T>
std::vector< disk::quadrature_point<T,3> >
integrate(const disk::simplicial_mesh<T,3>& msh,
          const typename disk::simplicial_mesh<T,3>::cell& cl,
          size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, cl);
    }

    auto m_quadrature_data = disk::tetrahedron_quadrature(degree);

    auto pts  = points(msh, cl);
    assert(pts.size() == 4);
    auto meas = measure(msh, cl);

    auto col1 = pts[1] - pts[0];
    auto col2 = pts[2] - pts[0];
    auto col3 = pts[3] - pts[0];

    auto tr = [&](const std::pair<point<T,3>, T>& qd) -> auto {
        auto point  = col1 * qd.first.x() + col2 * qd.first.y() + col3 * qd.first.z() + pts[0];
        auto weight = qd.second * meas;
        return disk::make_qp(point, weight);
    };

    std::vector<disk::quadrature_point<T,3>> ret(m_quadrature_data.size());
    std::transform(m_quadrature_data.begin(), m_quadrature_data.end(),
                   ret.begin(), tr);

    return ret;
}

template<typename T>
std::vector< disk::quadrature_point<T,3> >
integrate(const disk::simplicial_mesh<T,3>& msh,
          const typename disk::simplicial_mesh<T,3>::face& fc,
          size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, fc);
    }

    auto m_quadrature_data = disk::triangle_quadrature(degree);

    auto pts = points(msh, fc);
    assert(pts.size() == 3);
    auto meas = measure(msh, fc);

    auto col1 = pts[1] - pts[0];
    auto col2 = pts[2] - pts[0];

    auto tr = [&](const std::pair<point<T,2>, T>& qd) -> auto {
        auto point  = col1 * qd.first.x() + col2 * qd.first.y() + pts[0];
        auto weight = qd.second * meas;
        return disk::make_qp(point, weight);
    };

    std::vector<disk::quadrature_point<T,3>> ret(m_quadrature_data.size());
    std::transform(m_quadrature_data.begin(), m_quadrature_data.end(),
                   ret.begin(), tr);

    return ret;
}

template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate(const disk::generic_mesh<T, 3>& msh, const typename disk::generic_mesh<T, 3>::cell& cl, size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, cl);
    }

    auto pts = points(msh, cl);
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
integrate(const disk::generic_mesh<T, 3>& msh, const typename disk::generic_mesh<T, 3>::face& fc, size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, fc);
    }

    auto pts = points(msh, fc);
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
std::vector< disk::quadrature_point<T,3> >
integrate(const disk::cartesian_mesh<T,3>& msh,
		  const typename disk::cartesian_mesh<T,3>::cell& cl,
		  size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, cl);
    }

    auto pts = points(msh, cl);
    // transform arry to vector
    std::array< point<T,3>, 8 > ptsv{pts[0], pts[1], pts[3], pts[2], pts[4], pts[5], pts[7], pts[6]};
    return priv::integrate_hexahedron_tens(degree, ptsv);
}

template<typename T>
std::vector< disk::quadrature_point<T,3> >
integrate(const disk::cartesian_mesh<T,3>& msh,
          const typename disk::cartesian_mesh<T,3>::face& fc,
          size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, fc);
    }

    auto pts = points(msh, fc);
    auto meas              = measure(msh, fc);
    auto m_quadrature_data = proton::quadrangle_quadrature<T>(degree);

    std::vector<disk::quadrature_point<T, 3>> ret(m_quadrature_data.size());

    for(auto& qp : m_quadrature_data)
    {
        auto point  = (pts[1] - pts[0]) * qp.first.x() + (pts[3] - pts[0]) * qp.first.y() + pts[0];
        auto weight = qp.second * meas;
        ret.push_back(disk::make_qp(point, weight));
    }

    return ret;
}

} // revolution


