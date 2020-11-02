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
    typename Mesh::coordinate_type
    operator()(const Mesh& msh, size_t degree) const
    {
        return 0;
    }

    size_t
    expected_rate(size_t k)
    {
        return k + 50;
    } 
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct test_functor_equal_order<Mesh<T,3,Storage>>
{
    /* Expect k+1 convergence */
    using mesh_type = Mesh<T,3,Storage>;
    typename mesh_type::coordinate_type
    operator()(const mesh_type& msh, size_t degree) const
    {
        typedef typename mesh_type::cell        cell_type;
        typedef typename mesh_type::face        face_type;
        typedef typename mesh_type::coordinate_type scalar_type;
        typedef typename mesh_type::point_type  point_type;


        auto f = [](const point_type& pt) -> double
        {
            return std::sin(M_PI*pt.x());
        };

        auto sf = [](const point_type& pt) -> Matrix<scalar_type, 3, 1>
        {
            Matrix<scalar_type, 3, 1> ret;
            ret(0) = M_PI * std::cos(M_PI*pt.x());
            ret(1) = 0.0;
            ret(2) = 0.0;
            return ret;
        };

        typename disk::hho_degree_info hdi({ .rd = degree,
                                             .cd = degree,
                                             .fd = degree });

        scalar_type error = 0.0;
        for (auto& cl : msh)
        {
            Matrix<scalar_type, Dynamic, 1> proj = disk::project_function(msh, cl, hdi, f, 2);
            auto gr = disk::lapl_mm(msh, cl, hdi);

            size_t rec_size = disk::vector_basis_size(hdi.reconstruction_degree(), 3, 3);

            Matrix<scalar_type, Dynamic, 1> reconstr = Matrix<scalar_type, Dynamic, 1>::Zero(rec_size);
            reconstr = gr * proj;

            auto cb = disk::make_vector_monomial_basis(msh, cl, hdi.reconstruction_degree());
            //Matrix<scalar_type, Dynamic, Dynamic> mass = disk::make_mass_matrix(msh, cl, cb);
            //Matrix<scalar_type, Dynamic, 1> rhs = disk::make_rhs(msh, cl, cb, sf);
            //Matrix<scalar_type, Dynamic, 1> exp_reconstr = mass.llt().solve(rhs);


            Matrix<scalar_type, Dynamic, 1> exp_reconstr = project_function(msh, cl, cb, sf);
            Matrix<scalar_type, Dynamic, 1> diff = reconstr - exp_reconstr;

            Matrix<scalar_type, Dynamic, Dynamic> mass = disk::make_mass_matrix(msh, cl, cb);

            error += diff.dot(mass*diff);
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
    std::cout << blue << "Face order: k and Cell order: k" << std::endl;
    tester<test_functor_equal_order> tstr1;
    tstr1.run();

}
