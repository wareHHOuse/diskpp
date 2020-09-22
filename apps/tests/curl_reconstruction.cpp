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

#include <iomanip>
#include <iostream>
#include <regex>

#include <unistd.h>

#include "bases/bases.hpp"
#include "methods/hho"
#include "methods/implementation_hho/curl.hpp"
#include "quadratures/quadratures.hpp"

#include "core/loaders/loader.hpp"

#include "common.hpp"

template<typename Mesh, bool mixed>
struct test_functor_curl_reconstruction
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

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, bool mixed>
struct test_functor_curl_reconstruction<Mesh<T,3,Storage>, mixed>
{
    using mesh_type = Mesh<T,3,Storage>;
    /* Expect k+1 convergence (hho stabilization) */
    typename mesh_type::coordinate_type
    operator()(const mesh_type& msh, size_t degree) const
    {
        typedef typename mesh_type::cell            cell_type;
        typedef typename mesh_type::face            face_type;
        typedef typename mesh_type::coordinate_type scalar_type;
        typedef typename mesh_type::point_type      point_type;

        typedef Matrix<scalar_type, mesh_type::dimension, 1> ret_type;

        auto f = make_vector_testing_data(msh);

        size_t fd = degree;
        size_t cd = mixed ? degree+1 : degree;
        size_t rd = degree+1;

        disk::hho_degree_info hdi( { .rd = rd, .cd = cd, .fd = fd } );

        scalar_type error = 0.0;
        for (auto& cl : msh)
        {
            auto CR = disk::make_vector_hho_curl_impl(msh, cl, hdi);
            auto proj = disk::project_tangent(msh, cl, hdi, f, 1);
            Matrix<T, Dynamic, 1> rf = CR.first * proj;
            auto rb = disk::make_vector_monomial_basis(msh, cl, rd);

            Matrix<T, Dynamic, Dynamic> mass         = disk::make_mass_matrix(msh, cl, rb, 1);
            Matrix<T, Dynamic, 1>       rhs          = disk::make_rhs(msh, cl, rb, f, 1);
            Matrix<T, Dynamic, 1>       exp_reconstr = mass.ldlt().solve(rhs);
            Matrix<T, Dynamic, 1>       diff = exp_reconstr;
            diff.head(3) = exp_reconstr.head(3);
            diff.tail(rb.size()-3) -= rf;
        
            Matrix<T, Dynamic, Dynamic> CC = make_curl_curl_matrix(msh, cl, rb);
            error += diff.dot(CC * diff);

        }

        return std::sqrt(error);
    }

    size_t
    expected_rate(size_t k)
    {
        return k + 1;
    }
};

template<typename Mesh>
using test_functor_curl_reconstruction_eo = test_functor_curl_reconstruction<Mesh, false>;

template<typename Mesh>
using test_functor_curl_reconstruction_mo = test_functor_curl_reconstruction<Mesh, true>;

int
main(void)
{
    std::cout << red << "Test HHO curl reconstruction operator" << std::endl;
    // face order: k, cell order: k
    std::cout << cyan << "Face order: k and Cell order: k" << std::endl;
    tester<test_functor_curl_reconstruction_eo> tstr1;
    tstr1.run();

    std::cout << cyan << "Face order: k and Cell order: k+1" << std::endl;
    tester<test_functor_curl_reconstruction_mo> tstr2;
    tstr2.run();

    return 0;
}
