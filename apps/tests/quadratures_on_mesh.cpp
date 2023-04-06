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

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/quadratures/quadratures.hpp"
#include "common.hpp"

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
bool
test_integration_on_mesh(const Mesh<T, 2, Storage>& msh)
{
    /* max degree to test: pay attention that setting degree k
     * means testing up to degree 2*k, because we consider the
     * monomial x^m y^n where m and n go from 0 to k. */
    size_t max_test_degree = 8;

    /* max error in ULP */
    T ULP_max = 100;

    /* ASSUMES THAT THE MESH REPRESENTS THE DOMAIN [0,1]^d */
    auto analytic_integral = [](size_t m, size_t n) -> auto {
        return 1./((m+1)*(n+1));
    };

    auto monomial = [](const disk::point<T,2>& pt, size_t m, size_t n) -> T {
        return iexp_pow(pt.x(), m)*iexp_pow(pt.y(), n);
    };

    T area = 0.0;
    for (auto& cl : msh)
        area += measure(msh, cl);

    if ( !almost_equal(area, 1.0, 2.0/* Max 2ULP for area */) )
    {
        std::cout << bold << magenta << "  Area not computed accurately: ";
        std::cout << std::setprecision(16) << area << reset << std::endl;
        //return false;
    }

    for (size_t m = 0; m < max_test_degree; m++)
    {
        for (size_t n = 0; n < max_test_degree; n++)
        {
            T int_num = 0.0;
            for (auto& cl : msh)
            {
                auto qps = integrate(msh, cl, m+n);
                for (auto& qp : qps)
                    int_num += qp.weight()*monomial(qp.point(),m,n);
            }

            T int_ana = analytic_integral(m,n);

            if ( not almost_equal(int_ana, int_num, ULP_max) )
            {
                std::cout << "  Test FAIL for monomial x^" << m << " y^";
                std::cout << n << std::endl;
                std::cout << "   Analytical value = ";
                std::cout << std::setprecision(16) << int_ana << std::endl;
                std::cout << "   Numerical value  = ";
                std::cout << std::setprecision(16) << int_num << std::endl;
            }
        }
    }

    return true;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
bool
test_integration_on_mesh(const Mesh<T, 3, Storage>& msh)
{
    /* max degree to test: pay attention that setting degree k
     * means testing up to degree 3*k, because we consider the
     * monomial x^m y^n z^k where m, n and k go from 0 to k. */
    size_t max_test_degree = 5;

    /* max error in ULP */
    T ULP_max = 50;

    /* ASSUMES THAT THE MESH REPRESENTS THE DOMAIN [0,1]^d */
    auto analytic_integral = [](size_t m, size_t n, size_t k) -> auto {
        return 1./((m+1)*(n+1)*(k+1));
    };

    auto monomial = [](const disk::point<T,3>& pt, size_t m, size_t n, size_t k) -> T {
        return iexp_pow(pt.x(), m)*iexp_pow(pt.y(), n)*iexp_pow(pt.z(), k);
    };

    T volume = 0.0;
    for (auto& cl : msh)
        volume += measure(msh, cl);

    if ( !almost_equal(volume, 1.0, 2.0/* Max 2ULP for volume */) )
    {
        std::cout << bold << magenta << "  Volume not computed accurately: ";
        std::cout << std::setprecision(16) << volume << reset << std::endl;
        //return false;
    }

    for (size_t m = 0; m < max_test_degree; m++)
    {
        for (size_t n = 0; n < max_test_degree; n++)
        {
            for (size_t k = 0; k < max_test_degree; k++)
            {
                T int_num = 0.0;
                for (auto& cl : msh)
                {
                    auto qps = integrate(msh, cl, m+n+k);
                    for (auto& qp : qps)
                        int_num += qp.weight()*monomial(qp.point(),m,n,k);
                }

                T int_ana = analytic_integral(m,n,k);

                if ( not almost_equal(int_ana, int_num, ULP_max) )
                {
                std::cout << "  Test FAIL for monomial x^" << m << " y^";
                std::cout << n << " z^" << k << std::endl;
                std::cout << "   Analytical value = ";
                std::cout << std::setprecision(16) << int_ana << std::endl;
                std::cout << "   Numerical value  = ";
                std::cout << std::setprecision(16) << int_num << std::endl;
                }
            }
        }
    }

    return true;
}

template<typename T>
void test(const T& meshes)
{
    size_t num_meshes = meshes.size();
    size_t cur_mesh = 1;
    for (auto& msh : meshes)
    {
        std::cout << cyan << " Mesh " << cur_mesh << "/";
        std::cout << num_meshes << nocolor << std::endl;
        test_integration_on_mesh(msh);
        cur_mesh++;
    }
}

int main(void)
{
    std::cout << bold << yellow << "Testing triangular meshes (generic)" << reset << std::endl;
    test( get_triangle_generic_meshes<double>() );
    std::cout << bold << yellow << "Testing triangular meshes (netgen)" << reset << std::endl;
    test( get_triangle_netgen_meshes<double>() );
    std::cout << bold << yellow << "Testing quad meshes" << reset << std::endl;
    test( get_quad_generic_meshes<double>() );
    std::cout << bold << yellow << "Testing polygonal meshes" << reset << std::endl;
    test( get_polygonal_generic_meshes<double>() );
    std::cout << bold << yellow << "Testing tetrahedral meshes" << reset << std::endl;
    test( get_tetrahedra_netgen_meshes<double>() );
}


