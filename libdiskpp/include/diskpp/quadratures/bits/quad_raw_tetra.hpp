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

#error "Not ready yet"

namespace disk {

namespace quadrature {

namespace priv {

template<typename T>
inline T 
tetra_volume(const Eigen::Matrix<T,3,1>& v0, const Eigen::Matrix<T,3,1>& v1,
    const Eigen::Matrix<T,3,1>& v2 )
{
    return std::abs(v0.dot(v1.cross(v2))) / 6.;
}

} // namespace priv

template<typename T>
std::vector<quadrature_point<T,3>>
arbq(size_t degree, point<T,3> p0, point<T,3> p1, point<T,3> p2, point<T,3> p3)
{
    int rule = degree;
    if (rule == 0)
        rule = 1;
    
    int point_num = tetrahedron_arbq_size(rule);

    std::vector<double> ws(point_num);
    std::vector<double> pts(3*point_num);

    tetrahedron_arbq(rule, point_num, pts.data(), ws.data());

}


template<typename T>
std::vector<quadrature_point<T,3>>
grundmann_moeller(size_t degree, point<T,3> p0, point<T,3> p1, point<T,3> p2, point<T,3> p3)
{
    int rule = degree/2;
    int dimension = 3;
    int point_num = gm_rule_size(rule, dimension);

    std::vector<double> ws(point_num);
    std::vector<double> pts(3*point_num);

    gm_rule_set(rule, dimension, point_num, ws.data(), pts.data());

    std::vector<quadrature_point<T,3>> ret;
    ret.reserve(point_num);

    auto v0 = p1 - p0;
    auto v1 = p2 - p0;
    auto v2 = p3 - p0;

    for (size_t i = 0; i < point_num; i++)
    {
        auto r = pts[3*i];
        auto s = pts[3*i+1];
        auto t = pts[3*i+2];
        point<T,3> qp = p0 + r*v0 + s*v1 + t*v2;
        
        auto tvol = priv::tetra_volume(v0.to_vector(), v1.to_vector(), v2.to_vector());
        auto qw = ws[i] * tvol;
        ret.push_back( make_qp(qp, qw) );
    }

    return ret;
}

template<typename T>
std::vector<quadrature_point<T,3>>
default_tetrahedron_quadrature(size_t degree, point<T,3> p0, point<T,3> p1,
    point<T,3> p2, point<T,3> p3)
{
    return arbq(degree, p0, p1, p2, p3);
}

/* Tetrahedron quadrature. This is just a driver to call either ABRQ or GM
 * quadrature rules in the jburkardt's code. */
std::vector<std::pair<point<double,3>, double>>
tetrahedron_quadrature(size_t degree)
{


    
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

} // namespace quadrature

} // namespace disk

