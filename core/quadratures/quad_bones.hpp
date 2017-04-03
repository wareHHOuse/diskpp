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

std::vector<std::pair<point<double,3>, double>>
tetrahedron_quadrature(size_t degree)
{
    int dimension = 3;

#ifdef USE_ARBQ
    int rule = degree+1;
    if (rule == 0)
        rule = 1;
#else
    int rule = degree/2 + 1;
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


std::vector<std::pair<point<double,2>, double>>
triangle_quadrature(size_t order)
{
    std::vector<std::pair<point<double,2>, double>> ret;

    int rule_num = dunavant_rule_num();
    int rule, order_num, degree;

    for(rule = 1; rule <= rule_num; rule++)
    {
        degree = dunavant_degree(rule);
        if(degree > order)
            break;
    }
    assert(rule != rule_num or degree >= order);

    order_num = dunavant_order_num(rule);
    std::vector<double> xytab(2*order_num), wtab(order_num);

    dunavant_rule(rule, order_num, &xytab[0], &wtab[0]);

    ret.reserve( order_num );
    for(int iQN = 0; iQN < order_num; iQN++)
    {
        point<double,2> xQN;
        xQN.at(0) = xytab[0+iQN*2];
        xQN.at(1) = xytab[1+iQN*2];

        ret.push_back( std::make_pair(xQN, wtab[iQN]) );
    }

    return ret;
}

template<typename T>
std::vector<std::pair<point<T,1>, T>>
edge_quadrature(size_t doe)
{
    std::vector<std::pair<point<T,1>, T>> ret;

    using namespace Eigen;

    if (doe%2 == 0)
        doe++;

    size_t num_nodes = (doe+1)/2;

    if (num_nodes == 1)
    {
        auto qp = std::make_pair(point<T,1>({0.5}), 1.0);
        ret.push_back(qp);
        return ret;
    }

    dynamic_matrix<T> M = dynamic_matrix<T>::Zero(num_nodes, num_nodes);
    for (size_t i = 1; i < num_nodes; i++)
    {
        T p = 4.0 - 1.0/(i*i);
        M(i, i-1) = sqrt(1.0/p);
    }

    SelfAdjointEigenSolver<dynamic_matrix<T>> solver;
    solver.compute(M);

    dynamic_vector<T> weights = solver.eigenvectors().row(0);
                      weights = weights.array().square();
    dynamic_vector<T> nodes = solver.eigenvalues();

    assert( weights.size() == nodes.size() );

    for (size_t i = 0; i < nodes.size(); i++)
    {
        auto qp = std::make_pair(point<T,1>({(nodes(i) + 1.)/2.}), weights(i));
        ret.push_back(qp);
    }

    return ret;
}


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


} // namespace disk



///
