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

#include <limits>

#include <xmmintrin.h>

#include "bases/bases.hpp"
#include "quadratures/quadratures.hpp"
#include "methods/hho"

#include "common.hpp"

enum class quadtype
{
    GOLUB_WELSCH,
    GAUSS_LEGENDRE,
    DRIVER
};

/* Test 1D quadratures on [0,1] interval by integrating monomials */
template<typename T>
bool
test_1d_quadrature(quadtype qt)
{
    size_t max_test_degree = 20; /* max degree to test */
    T ULP_max = 50; /* max error in ULP */
    
    size_t failed_tests = 0;
    size_t total_tests = 0;
    
    auto transform = [](const std::pair<point<T, 1>, T>& qp) -> auto {
        point<T,1> p({0.5 + qp.first.x()/2.0});
        T w = 0.5*qp.second;
        return std::make_pair(p,w);
    };
    
    /* This is \int_0^1 x^n dx */
    auto analytic_integral = [](size_t degree) -> T {
        return 1.0/(degree+1);
    };
    
    auto monomial = [](const point<T,1>& p, size_t degree) -> T {
        return iexp_pow(p.x(), degree);
    };
    
    auto get_quadrature = [](quadtype qt, size_t degree) -> auto {
        /* These are the functions we are testing: see quad_bones.hpp */
        if ( qt == quadtype::GOLUB_WELSCH )
            return disk::golub_welsch<T>(degree);
        
        if ( qt == quadtype::GAUSS_LEGENDRE )
            return disk::gauss_legendre<T>(degree);
            
        if ( qt == quadtype::DRIVER )
            return disk::edge_quadrature<T>(degree);
            
        throw std::invalid_argument("Unknown quadrature to test");
    };
    
    for (size_t qd = 0; qd <= max_test_degree; qd++)
    {        
        auto qps = get_quadrature(qt, qd);
        
        for (size_t md = 0; md <= qd; md++)
        {
            T int_ana = analytic_integral(md);
            T int_num = 0.0;
            for (auto& qp : qps)
            { 
                auto real_qp = transform(qp);
                int_num += real_qp.second * monomial(real_qp.first, md);
            }
            
            if ( not almost_equal(int_ana, int_num, ULP_max) )
            {
                std::cout << "FAIL: md = " << md << ", qd = " << qd;
                std::cout << " analytical value = " << std::setprecision(16) << int_ana;
                std::cout << " numerical value = " << std::setprecision (16) << int_num;
                std::cout << std::endl;
                
                failed_tests++;
            }
            
            total_tests++;
        }
    }
    
    std::string qt_string;
    switch (qt)
    {
        case quadtype::GOLUB_WELSCH:
            qt_string = "Golub-Welsch quadrature test:   ";
            break;
            
        case quadtype::GAUSS_LEGENDRE:
            qt_string = "Gauss-Legendre quadrature test: ";
            break;
            
        case quadtype::DRIVER:
            qt_string = "Driver function test:           ";
            break;
    }
    
    std::cout << qt_string;
    
    if (failed_tests != 0)
    {
        std::cout << "[\x1b[31mFAIL\x1b[0m] (";
        std::cout << failed_tests << "/" << total_tests << " tests failed)";
        std::cout << std::endl;
        return false;
    }
    
    std::cout << "[ \x1b[32mOK\x1b[0m ] (" << total_tests << " tests passed)";
    std::cout << std::endl;
    return true;
}

/* Test 2D quadrature for triangles */
template<typename T>
bool
test_2d_triangle_quadrature(void)
{
    /* max degree to test: pay attention that setting degree k
     * means testing up to degree 2*k, because we consider the
     * monomial x^m y^n where m and n go from 0 to k. */
    size_t max_test_degree = 8;
    
    /* max error in ULP */
    T ULP_max = 100;
    
    size_t failed_tests = 0;
    size_t total_tests = 0;
    
    /* triangle on which we integrate */
    std::array<point<T,2>,3> tp;
    tp[0] = point<T,2>({0,0});
    tp[1] = point<T,2>({1,0});
    tp[2] = point<T,2>({1,1});
    
    /* This is \int_0^1 \int_0^x x^m y^n dy dx */
    auto analytic_integral = [](size_t degree_x, size_t degree_y) -> T {
        return 1.0/((degree_y+1)*(degree_x+degree_y+2));
    };
    
    auto monomial = [](const point<T,2>& pt, size_t m, size_t n) -> T {
        return iexp_pow(pt.x(), m)*iexp_pow(pt.y(), n);
    };
    
    for (size_t m = 0; m <= max_test_degree; m++)
    {
        for (size_t n = 0; n <= max_test_degree; n++)
        {
            /* This is the function we are testing: see quad_bones.hpp */
            auto qps = disk::priv::integrate_triangle<T>(m+n, tp);
            T int_num = 0.0;
            for (auto& qp : qps)
                int_num += qp.weight()*monomial(qp.point(),m,n);
            
            T int_ana = analytic_integral(m,n);
            
            if ( not almost_equal(int_ana, int_num, ULP_max) )
            {
                std::cout << "FAIL: m = " << m << ", n = " << n;
                std::cout << " analytical value = " << std::setprecision(16) << int_ana;
                std::cout << " numerical value = " << std::setprecision (16) << int_num;
                std::cout << std::endl;
                failed_tests++;
            }
            
            total_tests++;
        }
    }
    
    std::cout << "Triangle quadrature test:       ";
    if (failed_tests != 0)
    {
        std::cout << "[\x1b[31mFAIL\x1b[0m] (";
        std::cout << failed_tests << "/" << total_tests << " tests failed)";
        std::cout << std::endl;
        return false;
    }
    
    std::cout << "[ \x1b[32mOK\x1b[0m ] (" << total_tests << " tests passed)";
    std::cout << std::endl;
    return true;
}

/* Test 2D quadrature for triangles */
template<typename T>
bool
test_2d_quad_quadrature(void)
{
    /* max degree to test: pay attention that setting degree k
     * means testing up to degree 2*k, because we consider the
     * monomial x^m y^n where m and n go from 0 to k. */
    size_t max_test_degree = 10;
    
    /* max error in ULP */
    T ULP_max = 50;
    
    size_t failed_tests = 0;
    size_t total_tests = 0;
    
    /* triangle on which we integrate */
    std::vector<point<T,2>> tp;
    tp.resize(4);
    tp[0] = point<T,2>({0,0});
    tp[1] = point<T,2>({1,0});
    tp[2] = point<T,2>({1,0.5});
    tp[3] = point<T,2>({1,1});
    
    /* This is \int_0^1 \int_0^x x^m y^n dy dx */
    auto analytic_integral = [](size_t degree_x, size_t degree_y) -> T {
        return 1.0/((degree_y+1)*(degree_x+degree_y+2));
    };
    
    auto monomial = [](const point<T,2>& pt, size_t m, size_t n) -> T {
        return iexp_pow(pt.x(), m)*iexp_pow(pt.y(), n);
    };
    
    for (size_t m = 0; m <= max_test_degree; m++)
    {
        for (size_t n = 0; n <= max_test_degree; n++)
        {
            /* This is the function we are testing: see quad_bones.hpp */
            auto qps = disk::priv::integrate_quadrangle_tens<T>(m+n+2, tp);
            T int_num = 0.0;
            for (auto& qp : qps)
                int_num += qp.weight()*monomial(qp.point(),m,n);
            
            T int_ana = analytic_integral(m,n);
            
            if ( not almost_equal(int_ana, int_num, ULP_max) )
            {
                std::cout << "FAIL: m = " << m << ", n = " << n;
                std::cout << " analytical value = " << std::setprecision(16) << int_ana;
                std::cout << " numerical value = " << std::setprecision (16) << int_num;
                std::cout << std::endl;
                failed_tests++;
            }
            total_tests++;
        }
    }

    std::cout << "Quad tensorized test:           ";
    if (failed_tests != 0)
    {
        std::cout << "[\x1b[31mFAIL\x1b[0m] (";
        std::cout << failed_tests << "/" << total_tests << " tests failed)";
        std::cout << std::endl;
        return false;
    }
    
    std::cout << "[ \x1b[32mOK\x1b[0m ] (" << total_tests << " tests passed)";
    std::cout << std::endl;
    return true;
}

/* Test 2D quadrature for triangles */
template<typename T>
bool
test_tetrahedron_quadrature(void)
{
    /* max degree to test: pay attention that setting degree k
     * means testing up to degree 3*k, because we consider the
     * monomial x^m y^n z^k where m, n and k go from 0 to k. */
    size_t max_test_degree = 5;
    
    /* max error in ULP */
    T ULP_max = 12; /* ARBQ is fuckin' precise...use 12 here for it. */
    
    size_t failed_tests = 0;
    size_t total_tests = 0;
    
    /* tetrahedron on which we integrate */
    std::array<point<T,3>,4> tp;
    tp[0] = point<T,3>({0,0,0});
    tp[1] = point<T,3>({0,0,1});
    tp[2] = point<T,3>({1,0,1});
    tp[3] = point<T,3>({1,1,1});
    
    auto v0 = (tp[1] - tp[0]).to_vector();
    auto v1 = (tp[2] - tp[0]).to_vector();
    auto v2 = (tp[3] - tp[0]).to_vector();
    
    auto tv = std::abs( v0.dot(v1.cross(v2))/T(6) );
    
    /* This is \int_0^1 \int_0^x \int_x^1 x^m y^n z^k dz dy dx */
    auto analytic_integral = [](size_t m, size_t n, size_t k) -> T {
        T a = 1./((k+1)*(n+1));
        T b = 1./(m+n+2);
        T c = 1./(m+n+k+3);
        return a*(b-c);
    };
    
    auto monomial = [](const point<T,3>& pt, size_t m, size_t n, size_t k) -> T {
        return iexp_pow(pt.x(), m)*iexp_pow(pt.y(), n)*iexp_pow(pt.z(), k);
    };
    
    auto transform = [&](const std::pair<point<T,3>, T>& rqp) -> auto {
        auto p = rqp.first;
        auto w = rqp.second;
        auto rp = (tp[1] - tp[0]) * p.x() +
                  (tp[2] - tp[0]) * p.y() +
                  (tp[3] - tp[0]) * p.z() + tp[0];
        
        auto rw = w*tv;
        
        return std::make_pair(rp, rw);
    };
    
    for (size_t m = 0; m <= max_test_degree; m++)
    {
        for (size_t n = 0; n <= max_test_degree; n++)
        {
            for (size_t k = 0; k <= max_test_degree; k++)
            {
                /* This is the function we are testing: see quad_bones.hpp */
                auto rqps = disk::tetrahedron_quadrature(m+n+k);
                
                T int_num = 0.0;
                for (auto& rqp : rqps)
                {
                    auto qp = transform(rqp);
                    int_num += qp.second*monomial(qp.first,m,n,k);
                }
                T int_ana = analytic_integral(m,n,k);
            
                if ( not almost_equal(int_ana, int_num, ULP_max) )
                {
                    std::cout << "FAIL: m = " << m << ", n = " << n << ", k = " << k;
                    std::cout << " analytical value = " << std::setprecision(16) << int_ana;
                    std::cout << " numerical value = " << std::setprecision (16) << int_num;
                    std::cout << std::endl;
                    failed_tests++;
                }
                total_tests++;
            }
        }
    }
    
    std::cout << "Tetrahedron quadrature test:    ";
    if (failed_tests != 0)
    {
        std::cout << "[\x1b[31mFAIL\x1b[0m] (";
        std::cout << failed_tests << "/" << total_tests << " tests failed)";
        std::cout << std::endl;
        return false;
    }
    
    std::cout << "[ \x1b[32mOK\x1b[0m ] (" << total_tests << " tests passed)";
    std::cout << std::endl;
    return true;
}

int main(void)
{
    using T = double;
    
    int ret = EXIT_SUCCESS;

    if ( !test_1d_quadrature<T>(quadtype::GOLUB_WELSCH) )
        ret = EXIT_FAILURE;
    
    if ( !test_1d_quadrature<T>(quadtype::GAUSS_LEGENDRE) )
        ret = EXIT_FAILURE;
        
    if ( !test_1d_quadrature<T>(quadtype::DRIVER) )
        ret = EXIT_FAILURE;

    if (!test_2d_triangle_quadrature<T>())
        ret = EXIT_FAILURE;
    
    if (!test_2d_quad_quadrature<T>())
        ret = EXIT_FAILURE;
    
    if (!test_tetrahedron_quadrature<T>())
        ret = EXIT_FAILURE;
    
    return ret;
}
