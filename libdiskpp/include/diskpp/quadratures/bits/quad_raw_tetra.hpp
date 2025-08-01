/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#pragma once

#include <vector>

#include "diskpp/common/simplicial_formula.hpp"
#include "diskpp/mesh/point.hpp"
#include "diskpp/quadratures/quadrature_point.hpp"

#include "simplex_gm_rule.hpp"
#include "tetrahedron_arbq_rule.hpp"

namespace disk
{

namespace quadrature
{

namespace priv
{

static const double arbq_to_ref_A[3][3] = {{0.500000000000000, -0.288675134594813, -0.204124145231932},
                                           {0.000000000000000, 0.577350269189626, -0.204124145231932},
                                           {0.000000000000000, 0.000000000000000, 0.612372435695795}};

static const double arbq_to_ref_b[3] = {-1.000000000000000, -0.577350269189626, -0.408248290463863};

static const double inv_arbq_vol = 1.060660171779821;

} // namespace priv

template<typename T>
std::vector<quadrature_point<T, 3>>
arbq(size_t degree, const point<T, 3>& p0, const point<T, 3>& p1, const point<T, 3>& p2, const point<T, 3>& p3)
{
    int rule = degree;
    if (rule == 0)
        rule = 1;

    int point_num = tetrahedron_arbq_size(rule);

    std::vector<double> ws(point_num);
    std::vector<double> pts(3 * point_num);

    auto v0   = p1 - p0;
    auto v1   = p2 - p0;
    auto v2   = p3 - p0;
    auto tvol = volume_tetrahedron_kahan(p0, p1, p2, p3);

    tetrahedron_arbq(rule, point_num, pts.data(), ws.data());

    std::vector<quadrature_point<T, 3>> ret;
    ret.reserve(point_num);

    for (size_t pn = 0; pn < point_num; pn++)
    {
        double qpr[3];
        for (size_t i = 0; i < 3; i++)
        {
            qpr[i] = 0.0;
            for (size_t j = 0; j < 3; j++)
                qpr[i] += priv::arbq_to_ref_A[i][j] * (pts[3 * pn + j] - priv::arbq_to_ref_b[j]);
        }
        point<T, 3> qp = p0 + qpr[0] * v0 + qpr[1] * v1 + qpr[2] * v2;
        T           qw = ws[pn] * priv::inv_arbq_vol * tvol;
        ret.push_back(make_qp(qp, qw));
    }

    return ret;
}

/* Call the Grundmann-Moeller quadrature into the jburkardt's code.
 * The order of this quadrature is in theory arbitrary, but in practice
 * it works ok until order ~40. */
template<typename T>
std::vector<quadrature_point<T, 3>>
grundmann_moeller(size_t             degree,
                  const point<T, 3>& p0,
                  const point<T, 3>& p1,
                  const point<T, 3>& p2,
                  const point<T, 3>& p3)
{
    int rule      = (degree % 2) == 0 ? degree / 2 : (degree - 1) / 2;
    int dimension = 3;
    int point_num = gm_rule_size(rule, dimension);

    std::vector<double> ws(point_num);
    std::vector<double> pts(3 * point_num);

    gm_rule_set(rule, dimension, point_num, ws.data(), pts.data());

    std::vector<quadrature_point<T, 3>> ret;
    ret.reserve(point_num);

    auto v0   = p1 - p0;
    auto v1   = p2 - p0;
    auto v2   = p3 - p0;
    auto tvol = volume_tetrahedron_kahan(p0, p1, p2, p3);

    for (size_t i = 0; i < point_num; i++)
    {
        auto        r  = pts[3 * i];
        auto        s  = pts[3 * i + 1];
        auto        t  = pts[3 * i + 2];
        point<T, 3> qp = p0 + r * v0 + s * v1 + t * v2;
        auto        qw = ws[i] * tvol;
        ret.push_back(make_qp(qp, qw));
    }

    return ret;
}

template<typename T>
std::vector<quadrature_point<T, 3>>
default_tetrahedron_quadrature(size_t degree, point<T, 3> p0, point<T, 3> p1, point<T, 3> p2, point<T, 3> p3)
{
    return arbq(degree, p0, p1, p2, p3);
}

} // namespace quadrature

} // namespace disk
