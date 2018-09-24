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
 * Nicolas Pignet  (C) 2018                     nicolas.pignet@enpc.fr
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

#include <vector>
#include <utility>

#include "mesh/point.hpp"
#include "jburkardt/triangle_dunavant_rule.hpp"

#define USE_ARBQ

#ifdef USE_ARBQ
    #include "jburkardt/tetrahedron_arbq_rule.hpp"
#else
    #include "jburkardt/simplex_gm_rule.hpp"  /* Has negative weights */
#endif

namespace disk {

/* Golub-Welsch algorithm to get quadrature points on a line. Pay attention
 * that this functions solves an eigenvalue problem each time you call it,
 * so it is quite expensive. You should use the Gauss-Legendre quadrature
 * defined below, which calls Golub-Welsch only when really needed. */
template<typename T>
std::vector<std::pair<point<T, 1>, T>>
golub_welsch(const size_t degree)
{
    using namespace Eigen;
    typedef Matrix<T, Dynamic, Dynamic>    matrix_type;
    std::vector<std::pair<point<T, 1>, T>> ret;

    size_t comp_degree = degree;
    if (comp_degree % 2 == 0)
        comp_degree++;

    size_t reqd_nodes = (comp_degree + 1) / 2;

    if (reqd_nodes == 1)
    {
        auto qp = std::make_pair(point<T, 1>({0.}), 2.0);
        ret.push_back(qp);
        return ret;
    }

    matrix_type M = matrix_type::Zero(reqd_nodes, reqd_nodes);
    for (size_t i = 1; i < reqd_nodes; i++)
    {
        T p         = 4.0 - 1.0 / (i * i);
        M(i, i - 1) = sqrt(1.0 / p);
    }

    SelfAdjointEigenSolver<matrix_type> solver;
    solver.compute(M);

    matrix_type weights  = solver.eigenvectors().row(0);
    weights              = weights.array().square();
    matrix_type nodes    = solver.eigenvalues();

    assert(weights.size() == nodes.size());

    for (size_t i = 0; i < nodes.size(); i++)
    {
        auto qp = std::make_pair(point<T, 1>({nodes(i)}), 2 * weights(i));
        ret.push_back(qp);
    }

    return ret;
}

/* Gauss-Legendre 1D quadrature. Fall-back to Golub-Welsch if the required
 * degree is too high. */
template<typename T>
std::vector<std::pair<point<T, 1>, T>>
gauss_legendre(size_t degree)
{
    auto comp_degree = degree;

    if (degree % 2 == 0)
        comp_degree = degree + 1;

    size_t reqd_nodes = (comp_degree + 1) / 2;

    std::vector<std::pair<point<T, 1>, T>> ret;
    ret.reserve(reqd_nodes);

    point<T, 1> qp;
    T           qw;
    T           a1, a2;

    switch (reqd_nodes)
    {
        case 1:
            qp = point<T, 1>({0.0});
            qw = 2.0;
            ret.push_back(std::make_pair(qp, qw));
            return ret;

        case 2:
            qp = point<T, 1>({1.0 / std::sqrt(3.0)});
            qw = 1.0;
            ret.push_back(std::make_pair(-qp, qw));
            ret.push_back(std::make_pair(qp, qw));
            return ret;

        case 3:
            qp = point<T, 1>({std::sqrt(3.0 / 5.0)});
            qw = 5.0 / 9.0;
            ret.push_back(std::make_pair(-qp, qw));
            ret.push_back(std::make_pair(qp, qw));

            qp = point<T, 1>({0.0});
            qw = 8.0 / 9.0;
            ret.push_back(std::make_pair(qp, qw));

            return ret;

        case 4:
            a1 = 15;
            a2 = 2.0 * std::sqrt(30.0);

            qp = point<T, 1>({std::sqrt((a1 - a2) / 35.0)});
            qw = (18.0 + std::sqrt(30.0)) / 36.0;
            ret.push_back(std::make_pair(-qp, qw));
            ret.push_back(std::make_pair(qp, qw));

            qp = point<T, 1>({std::sqrt((a1 + a2) / 35.0)});
            qw = (18.0 - std::sqrt(30.0)) / 36.0;
            ret.push_back(std::make_pair(-qp, qw));
            ret.push_back(std::make_pair(qp, qw));

            return ret;

        case 5:
            qp = point<T, 1>({0.0});
            qw = 128.0 / 225.0;
            ret.push_back(std::make_pair(qp, qw));

            a1 = 5.0;
            a2 = 2.0 * std::sqrt(10.0 / 7.0);
            qp = point<T, 1>({std::sqrt(a1 - a2) / 3.0});
            qw = (322 + 13.0 * std::sqrt(70.0)) / 900.0;
            ret.push_back(std::make_pair(-qp, qw));
            ret.push_back(std::make_pair(qp, qw));

            qp = point<T, 1>({std::sqrt(a1 + a2) / 3.0});
            qw = (322 - 13.0 * std::sqrt(70.0)) / 900.0;
            ret.push_back(std::make_pair(-qp, qw));
            ret.push_back(std::make_pair(qp, qw));
            return ret;

        default:
            return golub_welsch<T>(degree);
    }

    return ret;
}

/* The client code should call this when a quadrature on the reference edge
 * is required. In our case the reference edge is [-1,1]. */
template<typename T>
std::vector<std::pair<point<T, 1>, T>>
edge_quadrature(const size_t doe)
{
    return gauss_legendre<T>(doe);
}

/* Tetrahedron quadrature. This is just a driver to call either ABRQ or GM
 * quadrature rules in the jburkardt's code. */
std::vector<std::pair<point<double,3>, double>>
tetrahedron_quadrature(size_t degree)
{
    int dimension = 3;

#ifdef USE_ARBQ
    int rule = degree;
    if (rule == 0)
        rule = 1;
#else
    int rule = degree/2;
#endif
    int point_num;

#ifdef USE_ARBQ
    point_num = tetrahedron_arbq_size(rule);
#else
    point_num = gm_rule_size(rule, dimension);
#endif
    std::vector<double> ws(point_num);
    std::vector<double> pts(3*point_num);

#ifdef USE_ARBQ
    tetrahedron_arbq(rule, point_num, &pts[0], &ws[0]);
#else
    gm_rule_set(rule, dimension, point_num, &ws[0], &pts[0]);
#endif

    std::vector<std::pair<point<double,3>, double>> ret;
    ret.reserve(point_num);

#ifdef USE_ARBQ
    static_matrix<double,3,3> invA;
    invA(0,0) = 1./2.;
    invA(0,1) = -1./(2.*std::sqrt(3.));
    invA(0,2) = -1./(2.*std::sqrt(6.));
    invA(1,0) = 0.;
    invA(1,1) = 1./std::sqrt(3.);
    invA(1,2) = -1./(2.*std::sqrt(6.));
    invA(2,0) = 0.;
    invA(2,1) = 0.;
    invA(2,2) = 3./(2.*std::sqrt(6.));

    static_vector<double,3> b;
    b(0) = -1.;
    b(1) = -1./std::sqrt(3.);
    b(2) = -1./std::sqrt(6.);
#endif

    for (size_t i = 0; i < point_num; i++)
    {
#ifdef USE_ARBQ
        double w = ws[i] / ( std::sqrt(8.) / 3. );
#else
        double w = ws[i];
#endif

        auto r = pts[3*i];
        auto s = pts[3*i+1];
        auto t = pts[3*i+2];

        auto p = point<double, 3>{r,s,t};

#ifdef USE_ARBQ
        static_vector<double,3> q = invA*(p.to_vector() - b);
        auto qq = point<double, 3>{q(0), q(1), q(2)};
        ret.push_back( std::make_pair(qq, w) );
#else
        ret.push_back( std::make_pair(p, w) );
#endif
    }

    return ret;
}


/* Triangle quadrature. This is just a driver to call Dunavant
 * quadrature rule in the jburkardt's code. */
std::vector<std::pair<point<double,2>, double>>
triangle_quadrature(size_t degree)
{
    std::vector<std::pair<point<double,2>, double>> ret;

    if ( degree == 0 )
        degree++;

    int max_degree = dunavant_rule_num();

    if ( degree > max_degree )
        throw std::invalid_argument("Degree too high");

    int num_points = dunavant_order_num(degree);
    std::vector<double> xytab(2*num_points), wtab(num_points);

    dunavant_rule(degree, num_points, &xytab[0], &wtab[0]);

    ret.reserve( num_points );
    for(int iQN = 0; iQN < num_points; iQN++)
    {
        point<double,2> xQN;
        xQN.at(0) = xytab[0+iQN*2];
        xQN.at(1) = xytab[1+iQN*2];

        ret.push_back( std::make_pair(xQN, wtab[iQN]) );
    }

    return ret;
}

/* Quadrature for cartesian quadrangles, it is just tensorized Gauss points. */
template<typename T>
std::vector<std::pair<point<T, 2>, T>>
quadrangle_quadrature(const size_t degree)
{
    const auto qps = disk::edge_quadrature<T>(degree);

    std::vector<std::pair<point<T, 2>, T>> ret;
    ret.reserve(qps.size() * qps.size());

    for (auto jtor = qps.begin(); jtor != qps.end(); jtor++)
    {
        const auto qp_y = *jtor;
        const auto eta  = qp_y.first.x();

        for (auto itor = qps.begin(); itor != qps.end(); itor++)
        {
            const auto qp_x = *itor;
            const auto xi   = qp_x.first.x();

            const auto w = qp_x.second * qp_y.second;

            ret.push_back(std::make_pair(point<T, 2>({xi, eta}), w));
        }
    }

    return ret;
}

/* This class represents a quadrature point, which is composed by the
 * coordinates (method .point()) and the weight (method .weight()) */
template<typename T, size_t DIM>
class quadrature_point : private std::pair<point<T,DIM>, T>
{
    typedef std::pair<point<T,DIM>, T> base_type;

public:
    quadrature_point()
    {}

    quadrature_point(const point<T, DIM>& p, const T& w)
        : base_type(p, w)
    {}

    point<T, DIM>
    point() const
    { return this->first; }

    T
    weight() const
    { return this->second; }
};

template<typename T, size_t DIM>
quadrature_point<T, DIM>
make_qp(const point<T, DIM>& point, const T& weight)
{
    return quadrature_point<T, DIM>(point, weight);
}

template<typename T, size_t DIM>
std::ostream&
operator<<(std::ostream& os, const quadrature_point<T,DIM>& qp)
{
    os << qp.point() << " " << qp.weight();
    return os;
}


/* Private stuff to support quadratures. User should not rely on this stuff. */
namespace priv {

/* See Ern & Guermond - Theory and practice of FEM, pag 360. */
/*
template<typename T>
std::vector<quadrature_point<T,2>>
triangle_quadrature_low_order(const point<T,2>& p0,
                              const point<T,2>& p1, 
                              const point<T,2>& p2, size_t deg)
{
    std::vector<quadrature_point<T,2>>   ret;
    auto v0 = p1 - p0;
    auto v1 = p2 - p0;
    auto area = std::abs( (v0.x() * v1.y() - v0.y() * v1.x())/2.0 );
    point<T,2>      qp;
    T               qw;
    T               a1 = (6. - std::sqrt(15.)) / 21;
    T               a2 = (6. + std::sqrt(15.)) / 21;
    T               w1 = (155. - std::sqrt(15.)) / 1200;
    T               w2 = (155. + std::sqrt(15.)) / 1200;
    switch(deg)
    {
        case 0:
        case 1:
            qw = area;
            qp = (p0 + p1 + p2)/3;      ret.push_back( make_qp(qp, qw) );
            return ret;
        case 2:
            qw = area/3;
            qp = p0/6 + p1/6 + 2*p2/3;  ret.push_back( make_qp(qp, qw) );
            qp = p0/6 + 2*p1/3 + p2/6;  ret.push_back( make_qp(qp, qw) );
            qp = 2*p0/3 + p1/6 + p2/6;  ret.push_back( make_qp(qp, qw) );
            return ret;
        case 3:
            qw = 9*area/20;
            qp = (p0 + p1 + p2)/3;      ret.push_back( make_qp(qp, qw) );
            qw = 2*area/15;
            qp = (p0 + p1)/2;           ret.push_back( make_qp(qp, qw) );
            qp = (p0 + p2)/2;           ret.push_back( make_qp(qp, qw) );
            qp = (p1 + p2)/2;           ret.push_back( make_qp(qp, qw) );
            qw = area/20;
            qp = p0;                    ret.push_back( make_qp(qp, qw) );
            qp = p1;                    ret.push_back( make_qp(qp, qw) );
            qp = p2;                    ret.push_back( make_qp(qp, qw) );
            return ret;
        case 4:
        case 5:
            qw = 9*area/40;
            qp = (p0 + p1 + p2)/3;      ret.push_back( make_qp(qp, qw) );
            qw = w1 * area;
            qp = a1*p0 + a1*p1 + (1-2*a1)*p2;   ret.push_back( make_qp(qp, qw) );
            qp = a1*p0 + (1-2*a1)*p1 + a1*p2;   ret.push_back( make_qp(qp, qw) );
            qp = (1-2*a1)*p0 + a1*p1 + a1*p2;   ret.push_back( make_qp(qp, qw) );
            qw = w2 * area;
            qp = a2*p0 + a2*p1 + (1-2*a2)*p2;   ret.push_back( make_qp(qp, qw) );
            qp = a2*p0 + (1-2*a2)*p1 + a2*p2;   ret.push_back( make_qp(qp, qw) );
            qp = (1-2*a2)*p0 + a2*p1 + a2*p2;   ret.push_back( make_qp(qp, qw) );
            return ret;
            
        default:
            throw std::invalid_argument("Triangle quadrature: requested order too high");
    }
    return ret;
}
*/
/* Get quadrature points for a triangle specified as a list of points. The
 * list of points is contained in a STL random-access containter (PtA) */
template<typename T, typename PtA>
std::vector<disk::quadrature_point<T, 2>>
integrate_triangle(size_t degree, const PtA& pts)
{
    assert(pts.size() == 3);
    static_assert(std::is_same<typename PtA::value_type, point<T, 2>>::value, "This function is for 2D points");

    using quadpoint_type = disk::quadrature_point<T, 2>;

    const auto qps = disk::triangle_quadrature(degree);

    std::vector<quadpoint_type> ret;

    ret.resize(qps.size());

    const auto col1 = pts[1] - pts[0];
    const auto col2 = pts[2] - pts[0];

    /* Compute the area of the sub-triangle */
    const auto tm = (col1.x() * col2.y() - col2.x() * col1.y()) / 2.;

    auto tr = [&](const std::pair<point<T, 2>, T>& qd) -> auto
    {
        const auto point  = col1 * qd.first.x() + col2 * qd.first.y() + pts[0];
        const auto weight = qd.second * std::abs(tm);
        return disk::make_qp(point, weight);
    };

    std::transform(qps.begin(), qps.end(), ret.begin(), tr);

    return ret;
}

/* Integrate using tensorized Gauss points on any quadrangle (not only
 * cartesian quadrangles) */
template<typename T>
std::vector<disk::quadrature_point<T, 2>>
integrate_quadrangle_tens(size_t degree, const std::vector<point<T, 2>>& pts)
{
    const auto qps = disk::quadrangle_quadrature<T>(degree);

    std::vector<disk::quadrature_point<T, 2>> ret;
    ret.reserve(qps.size());

    auto P = [&](T xi, T eta) -> T {
        return 0.25 * pts[0].x() * (1 - xi) * (1 - eta) +
               0.25 * pts[1].x() * (1 + xi) * (1 - eta) +
               0.25 * pts[2].x() * (1 + xi) * (1 + eta) +
               0.25 * pts[3].x() * (1 - xi) * (1 + eta);
    };

    auto Q = [&](T xi, T eta) -> T {
        return 0.25 * pts[0].y() * (1 - xi) * (1 - eta) +
               0.25 * pts[1].y() * (1 + xi) * (1 - eta) +
               0.25 * pts[2].y() * (1 + xi) * (1 + eta) +
               0.25 * pts[3].y() * (1 - xi) * (1 + eta);
    };

    auto J = [&](T xi, T eta) -> T {
        auto j11 = 0.25 * ( (pts[1].x() - pts[0].x()) * (1 - eta) +
                            (pts[2].x() - pts[3].x()) * (1 + eta) );
        
        auto j12 = 0.25 * ( (pts[1].y() - pts[0].y()) * (1 - eta) +
                            (pts[2].y() - pts[3].y()) * (1 + eta) );
        
        auto j21 = 0.25 * ( (pts[3].x() - pts[0].x()) * (1 - xi) +
                            (pts[2].x() - pts[1].x()) * (1 + xi) );
        
        auto j22 = 0.25 * ( (pts[3].y() - pts[0].y()) * (1 - xi) +
                            (pts[2].y() - pts[1].y()) * (1 + xi) );

        return std::abs(j11 * j22 - j12 * j21);
    };

    for (auto& qp : qps)
    {
        const T xi  = qp.first.x();
        const T eta = qp.first.y();

        const T px = P(xi, eta);
        const T py = Q(xi, eta);

        const T weight = qp.second * J(xi, eta);
        ret.push_back(disk::make_qp(point<T, 2>({px, py}), weight));
    }

    return ret;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename PtA>
std::vector<disk::quadrature_point<T, 2>>
integrate_quadrangle(const Mesh<T, 2, Storage>&                msh,
                     const typename Mesh<T, 2, Storage>::cell& cl,
                     size_t                                    degree,
                     const PtA&                                pts)
{
    assert(pts.size() == 4);

    using mesh_type      = Mesh<T, 2, Storage>;
    using quadpoint_type = disk::quadrature_point<T, 2>;

    const auto qps = disk::triangle_quadrature(degree);

    std::vector<quadpoint_type> ret;

    ret.resize(qps.size() * 2);

    for (size_t i = 1; i < 3; i++)
    {
        const auto pt1  = pts[i];
        const auto pt2  = pts[i + 1];
        const auto col1 = pt1 - pts[0];
        const auto col2 = pt2 - pts[0];

        /* Compute the area of the sub-triangle */
        const auto tm = (col1.x() * col2.y() - col2.x() * col1.y()) / 2.;

        auto tr = [&](const std::pair<point<T, 2>, T>& qd) -> auto
        {
            const auto point  = col1 * qd.first.x() + col2 * qd.first.y() + pts[0];
            const auto weight = qd.second * std::abs(tm);
            return disk::make_qp(point, weight);
        };

        auto retbegin = ret.begin();
        std::advance(retbegin, qps.size() * (i - 1));

        std::transform(qps.begin(), qps.end(), retbegin, tr);
    }

    return ret;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
std::vector<disk::quadrature_point<T, 2>>
integrate_2D_face(const Mesh<T, 2, Storage>& msh, const typename Mesh<T, 2, Storage>::face& fc, size_t degree)
{
    using mesh_type  = disk::simplicial_mesh<T, 2>;
    using point_type = typename mesh_type::point_type;

    const auto qps = disk::edge_quadrature<T>(degree);
    const auto pts = points(msh, fc);

    const auto scale = (pts[1] - pts[0]);
    const auto meas  = scale.to_vector().norm();

    std::vector<disk::quadrature_point<T, 2>> ret;
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

} // end priv


} // namespace disk



///
