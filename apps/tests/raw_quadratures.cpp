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

template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
    almost_equal(T x, T y, T ulp)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x-y) <= std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
    // unless the result is subnormal
           || std::abs(x-y) < std::numeric_limits<T>::min();
}

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
    
    auto transform = [](const std::pair<point<T, 1>, T>& qp) -> auto {
        point<T,1> p({0.5 + qp.first.x()/2.0});
        T w = 0.5*qp.second;
        return std::make_pair(p,w);
    };
    
    auto analytic_integral = [](size_t degree) -> T {
        return 1.0/(degree+1);
    };
    
    auto monomial = [](const point<T,1>& p, size_t degree) -> T {
        return iexp_pow(p.x(), degree);
    };
    
    auto get_quadrature = [](quadtype qt, size_t degree) -> auto {
        if ( qt == quadtype::GOLUB_WELSCH )
            return disk::golub_welsch<T>(degree);
        
        if ( qt == quadtype::GAUSS_LEGENDRE )
            return disk::gauss_legendre<T>(degree);
            
        if ( qt == quadtype::DRIVER )
            return disk::edge_quadrature<T>(degree);
            
        throw std::invalid_argument("Unknown quadrature to test");
    };

    bool success = true;
    
    for (size_t qd = 0; qd < max_test_degree; qd++)
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
            
            //auto relerr = std::abs(int_ana - int_num)/int_ana;
            //std::cout << relerr << std::endl;
            
            if ( not almost_equal(int_ana, int_num, ULP_max) )
            {
                std::cout << "FAIL: md = " << md << ", qd = " << qd;
                std::cout << " analytical value = " << std::setprecision(16) << int_ana;
                std::cout << " numerical value = " << std::setprecision (16) << int_num;
                std::cout << std::endl;
                success = false;
            }
        }
    } 
    
    return success;
}

void
test_1d_driver(const char *name, quadtype qt, int& status)
{
    std::string str_success = "[ \x1b[32mOK\x1b[0m ]";
    std::string str_failure = "[\x1b[31mFAIL\x1b[0m]";
    
    std::cout << "Testing 1D " << name << ": ";
    if ( not test_1d_quadrature<double>(qt) )
    {
        std::cout << str_failure << std::endl;
        status = EXIT_FAILURE; 
    }
    else
    {
        std::cout << str_success << std::endl;
    }
}

int main(void)
{
    int ret = EXIT_SUCCESS;

    test_1d_driver("Golub-Welsch", quadtype::GOLUB_WELSCH, ret);
    test_1d_driver("Gauss-Legendre", quadtype::GAUSS_LEGENDRE, ret);
    test_1d_driver("general driver", quadtype::DRIVER, ret);
    
    return ret;
}
