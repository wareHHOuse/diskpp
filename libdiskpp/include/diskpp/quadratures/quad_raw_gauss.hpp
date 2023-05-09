#pragma once

#include <vector>

#include "diskpp/quadratures/quadrature_point.hpp"

namespace disk::quadrature {

namespace priv {

struct gauss_point {
    double point;
    double weight;
};

/* 1-point Gauss rule */
static gauss_point gauss_rule_1[] = {
    {                        0,                        2 }
};

/* 2-point Gauss rule */
static gauss_point gauss_rule_2[] = {
    {       -0.577350269189626,                        1 },
    {        0.577350269189626,                        1 }
};

/* 3-point Gauss rule */
static gauss_point gauss_rule_3[] = {
    {       -0.774596669241483,        0.555555555555556 },
    {                        0,        0.888888888888889 },
    {        0.774596669241483,        0.555555555555556 }
};

/* 4-point Gauss rule */
static gauss_point gauss_rule_4[] = {
    {       -0.861136311594053,        0.347854845137454 },
    {       -0.339981043584856,        0.652145154862546 },
    {        0.339981043584856,        0.652145154862546 },
    {        0.861136311594052,        0.347854845137453 }
};

/* 5-point Gauss rule */
static gauss_point gauss_rule_5[] = {
    {       -0.906179845938664,        0.236926885056189 },
    {       -0.538469310105683,        0.478628670499366 },
    {                        0,        0.568888888888889 },
    {        0.538469310105683,        0.478628670499367 },
    {        0.906179845938664,        0.236926885056189 }
};

/* 6-point Gauss rule */
static gauss_point gauss_rule_6[] = {
    {       -0.932469514203152,         0.17132449237917 },
    {       -0.661209386466264,        0.360761573048138 },
    {       -0.238619186083197,        0.467913934572691 },
    {        0.238619186083197,        0.467913934572691 },
    {        0.661209386466264,        0.360761573048138 },
    {        0.932469514203152,        0.171324492379171 }
};

/* 7-point Gauss rule */
static gauss_point gauss_rule_7[] = {
    {       -0.949107912342759,         0.12948496616887 },
    {       -0.741531185599394,        0.279705391489276 },
    {       -0.405845151377397,        0.381830050505119 },
    {                        0,        0.417959183673469 },
    {        0.405845151377397,        0.381830050505119 },
    {        0.741531185599394,        0.279705391489276 },
    {        0.949107912342758,         0.12948496616887 }
};

/* 8-point Gauss rule */
static gauss_point gauss_rule_8[] = {
    {       -0.960289856497536,        0.101228536290376 },
    {       -0.796666477413627,        0.222381034453375 },
    {       -0.525532409916329,        0.313706645877887 },
    {        -0.18343464249565,        0.362683783378362 },
    {         0.18343464249565,        0.362683783378362 },
    {        0.525532409916329,        0.313706645877888 },
    {        0.796666477413627,        0.222381034453374 },
    {        0.960289856497536,        0.101228536290377 }
};

/* 9-point Gauss rule */
static gauss_point gauss_rule_9[] = {
    {       -0.968160239507626,       0.0812743883615746 },
    {       -0.836031107326636,        0.180648160694858 },
    {        -0.61337143270059,        0.260610696402936 },
    {       -0.324253423403809,        0.312347077040003 },
    {                        0,         0.33023935500126 },
    {        0.324253423403809,        0.312347077040003 },
    {        0.613371432700591,        0.260610696402936 },
    {        0.836031107326637,        0.180648160694857 },
    {        0.968160239507627,       0.0812743883615752 }
};

/* 10-point Gauss rule */
static gauss_point gauss_rule_10[] = {
    {       -0.973906528517171,       0.0666713443086882 },
    {       -0.865063366688984,         0.14945134915058 },
    {       -0.679409568299024,        0.219086362515982 },
    {       -0.433395394129247,        0.269266719309996 },
    {       -0.148874338981631,        0.295524224714753 },
    {        0.148874338981631,        0.295524224714753 },
    {        0.433395394129247,        0.269266719309997 },
    {        0.679409568299025,        0.219086362515982 },
    {        0.865063366688985,        0.149451349150582 },
    {        0.973906528517172,       0.0666713443086878 }
};

struct gauss_rule {
    size_t          num_entries;
    gauss_point     *points;
};


static struct gauss_rule gauss_rules[] = {
    {  sizeof(gauss_rule_1)/(sizeof(gauss_point)), gauss_rule_1  },
    {  sizeof(gauss_rule_2)/(sizeof(gauss_point)), gauss_rule_2  },
    {  sizeof(gauss_rule_3)/(sizeof(gauss_point)), gauss_rule_3  },
    {  sizeof(gauss_rule_4)/(sizeof(gauss_point)), gauss_rule_4  },
    {  sizeof(gauss_rule_5)/(sizeof(gauss_point)), gauss_rule_5  },
    {  sizeof(gauss_rule_6)/(sizeof(gauss_point)), gauss_rule_6  },
    {  sizeof(gauss_rule_7)/(sizeof(gauss_point)), gauss_rule_7  },
    {  sizeof(gauss_rule_8)/(sizeof(gauss_point)), gauss_rule_8  },
    {  sizeof(gauss_rule_9)/(sizeof(gauss_point)), gauss_rule_9  },
    { sizeof(gauss_rule_10)/(sizeof(gauss_point)), gauss_rule_10 },
    { 0, NULL }
};

} // namespace priv

template<typename T>
std::vector<quadrature_point<T,1>>
gauss_legendre(size_t degree, const T& a, const T& b)
{
    const auto num_rules = sizeof(priv::gauss_rules)/sizeof(priv::gauss_rule) - 1;
    const auto rule_num = degree/2;
    if (rule_num >= num_rules)
        throw std::invalid_argument("gauss_legendre: order too high");

    auto npts = priv::gauss_rules[rule_num].num_entries;
    std::vector<quadrature_point<T,1>> qps;
    qps.reserve(npts);
    for (size_t i = 0; i < npts; i++) {
        const auto &qp = priv::gauss_rules[rule_num].points[i];
        auto tr_qp = 0.5*(a+b) + 0.5*(b-a)*qp.point;
        auto tr_qw = (b-a)*qp.weight;
        qps.push_back({tr_qp, tr_qw});
    }

    return qps;
}

template<typename T, size_t DIM>
std::vector<quadrature_point<T,DIM>>
gauss_legendre(size_t degree, const point<T,DIM>& a, const point<T,DIM>& b)
{
    const auto num_rules = sizeof(priv::gauss_rules)/sizeof(priv::gauss_rule) - 1;
    const auto rule_num = degree/2;
    if (rule_num >= num_rules)
        throw std::invalid_argument("gauss_legendre: order too high");

    const auto& qrule = priv::gauss_rules[rule_num];
    auto npts = qrule.num_entries;
    std::vector<quadrature_point<T,DIM>> qps;
    qps.reserve(npts);
    for (size_t i = 0; i < npts; i++) {
        const auto &qp = qrule.points[i];
        auto tr_qp = 0.5 * (a + b) + 0.5 * (b - a) * qp.point;
        auto tr_qw = 0.5 * distance(a,b) * qp.weight;
        qps.push_back({tr_qp, tr_qw});
    }

    return qps;
}

template<typename T>
std::vector<quadrature_point<T,2>>
tensorized_gauss_legendre(size_t degree, const point<T,2>& p0,
    const point<T,2>& p1, const point<T,2>& p2, const point<T,2>& p3)
{
    const auto num_rules = sizeof(priv::gauss_rules)/sizeof(priv::gauss_rule) - 1;
    const auto rule_num = degree/2;
    if (rule_num >= num_rules)
        throw std::invalid_argument("tensorized_gauss_legendre: order too high");

    const auto& qrule = priv::gauss_rules[rule_num];

    auto npts = qrule.num_entries;
    std::vector<quadrature_point<T,2>> qps;
    qps.reserve(npts*npts); /* Collect the points in the [-1, 1]^2 square */
    for (size_t i = 0; i < npts; i++) {
        const auto &qpy = qrule.points[i];
        for (size_t j = 0; j < npts; j++) {
            const auto &qpx = qrule.points[j];
            point<T,2> qp(qpx.point, qpy.point);
            auto qw = qpx.weight * qpy.weight;
            qps.push_back({qp, qw});
        }
    }

    /* points in pts must be in _counterclockwise_ order */
    auto X = [&](T xi, T eta) -> T {
        return 0.25 * p0.x() * (1 - xi) * (1 - eta) +
               0.25 * p1.x() * (1 + xi) * (1 - eta) +
               0.25 * p2.x() * (1 + xi) * (1 + eta) +
               0.25 * p3.x() * (1 - xi) * (1 + eta);
    };

    auto Y = [&](T xi, T eta) -> T {
        return 0.25 * p0.y() * (1 - xi) * (1 - eta) +
               0.25 * p1.y() * (1 + xi) * (1 - eta) +
               0.25 * p2.y() * (1 + xi) * (1 + eta) +
               0.25 * p3.y() * (1 - xi) * (1 + eta);
    };

    auto J = [&](T xi, T eta) -> T {
        auto j11 = 0.25 * ( (p1.x() - p0.x()) * (1 - eta) +
                            (p2.x() - p3.x()) * (1 + eta) );

        auto j12 = 0.25 * ( (p1.y() - p0.y()) * (1 - eta) +
                            (p2.y() - p3.y()) * (1 + eta) );

        auto j21 = 0.25 * ( (p3.x() - p0.x()) * (1 - xi) +
                            (p2.x() - p1.x()) * (1 + xi) );

        auto j22 = 0.25 * ( (p3.y() - p0.y()) * (1 - xi) +
                            (p2.y() - p1.y()) * (1 + xi) );

        return std::abs(j11 * j22 - j12 * j21);
    };

    /* do ref->phys transform in place */
    for (auto& qp : qps)
    {
        const T xi  = qp.point().x();
        const T eta = qp.point().y();

        point<T,2> phys_point(X(xi, eta), Y(xi, eta));
        const T phys_weight = qp.weight() * J(xi, eta);
        qp = disk::make_qp(phys_point, phys_weight);
    }

    return qps;
}

template<typename T>
std::vector<quadrature_point<T,2>>
tensorized_gauss_legendre(size_t degree, const std::array<point<T,2>, 4>& pts)
{
    return tensorized_gauss_legendre(degree, pts[0], pts[1], pts[2], pts[3]);
}

} // namespace disk::quadrature