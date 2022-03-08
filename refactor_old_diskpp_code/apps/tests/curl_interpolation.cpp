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
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
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

#include <iostream>
#include <iomanip>
#include <regex>

#include "contrib/colormanip.h"
#include <unistd.h>

#include <xmmintrin.h>
//#define EIGEN_USE_MKL_ALL
#include "bases/bases.hpp"
#include "quadratures/quadratures.hpp"
#include "methods/hho"

#include "core/loaders/loader.hpp"

#include "common.hpp"


template<typename Mesh>
struct test_functor_equal_order
{
    /* Expect k+1 convergence */
    typename Mesh::coordinate_type
    operator()(const Mesh& msh, size_t degree) const
    {
        typedef Mesh mesh_type;
        typedef typename mesh_type::cell        cell_type;
        typedef typename mesh_type::face        face_type;
        typedef typename mesh_type::coordinate_type scalar_type;
        typedef typename mesh_type::point_type  point_type;

        using T = scalar_type;
        static const size_t DIM = Mesh::dimension;


        auto f = [](const disk::point<T, DIM>& pt) -> Matrix<T,DIM,1> {
            Matrix<T,DIM,1> ret;
            ret(0) = std::sin(M_PI*pt.x());
            ret(1) = std::sin(M_PI*pt.y());
            if constexpr (DIM==3)
                ret(2) = std::sin(M_PI*pt.z());
            return ret;
        };

        typename disk::hho_degree_info hdi(degree);

        scalar_type error = 0.0;
        for (auto& cl : msh)
        {
            Matrix<scalar_type, Dynamic, 1> proj = disk::project_tangent(msh, cl, hdi, f, 1);

            auto cb = disk::make_vector_monomial_basis(msh, cl, hdi.cell_degree());

            Matrix<scalar_type, Dynamic, 1> cp = proj.segment(0, cb.size()); 
            auto qps = integrate(msh, cl, 2*hdi.cell_degree());
            for (auto& qp : qps)
            {
                if constexpr (DIM==3){
                    Matrix<T, Dynamic, 3> phi = cb.eval_functions(qp.point());
                    disk::dynamic_vector<T> ediff = disk::eval(cp, phi) - f(qp.point());
                    error += qp.weight() * ediff.dot(ediff);
                }
            }
        }

        return std::sqrt(error);
    }

    size_t
    expected_rate(size_t k)
    {
        return k+1;
    }
};


template<typename Mesh>
struct test_functor_equal_order_face
{
    /* Expect k+1 convergence */
    typename Mesh::coordinate_type
    operator()(const Mesh& msh, size_t degree) const
    {
        typedef Mesh mesh_type;
        typedef typename mesh_type::cell        cell_type;
        typedef typename mesh_type::face        face_type;
        typedef typename mesh_type::coordinate_type scalar_type;
        typedef typename mesh_type::point_type  point_type;

        using T = scalar_type;
        static const size_t DIM = Mesh::dimension;


        auto f = [](const disk::point<T, DIM>& pt) -> Matrix<T,DIM,1> {
            Matrix<T,DIM,1> ret;
            ret(0) = std::sin(M_PI*pt.x());
            ret(1) = std::sin(M_PI*pt.y());
            if constexpr (DIM==3)
                ret(2) = std::sin(M_PI*pt.z());
            return ret;
        };

        typename disk::hho_degree_info hdi(degree);

        scalar_type error = 0.0;
        for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
        {
            if constexpr (DIM==3)
            {
                auto fc = *itor;
                auto fb = disk::make_vector_monomial_tangential_basis(msh, fc, hdi.face_degree());

                Matrix<scalar_type, Dynamic, 1> fp = disk::project_tangent(msh, fc, hdi.face_degree(), f, 1);
                auto qps = integrate(msh, fc, 2*hdi.face_degree());
                for (auto& qp : qps)
                {
                    Matrix<T, Dynamic, 3> phi = fb.eval_functions(qp.point());
                    auto n = normal(msh, fc);
                    Matrix<T,DIM,1> fv = n.cross( f(qp.point()).cross(n) );
                    disk::dynamic_vector<T> ediff = disk::eval(fp, phi) - fv;
                    error += qp.weight() * ediff.dot(ediff);
                }
            }
        }

        return std::sqrt(error);
    }

    size_t
    expected_rate(size_t k)
    {
        return k+1;
    }
};


int main(void)
{
    // face order: k, cell order: k
    std::cout << cyan << "Cell order: k" << std::endl;
    tester<test_functor_equal_order> tstr1;
    tstr1.run();
    // face order: k, cell order: k+1
    std::cout << cyan << "Face order: k" << std::endl;
    tester<test_functor_equal_order_face> tstr2;
    tstr2.run();
    return 0;
}
