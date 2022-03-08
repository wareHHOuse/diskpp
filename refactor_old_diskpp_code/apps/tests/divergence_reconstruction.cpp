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
#include "quadratures/quadratures.hpp"

#include "core/loaders/loader.hpp"

#include "common.hpp"

template<typename Mesh>
struct test_functor_equal_order
{
    /* Expect k+1 convergence */
    typename Mesh::coordinate_type
    operator()(const Mesh& msh, size_t degree) const
    {
        typedef Mesh                                mesh_type;
        typedef typename mesh_type::cell            cell_type;
        typedef typename mesh_type::face            face_type;
        typedef typename mesh_type::coordinate_type scalar_type;
        typedef typename mesh_type::point_type      point_type;

        typedef Matrix<scalar_type, Mesh::dimension, 1> ret_type;

        auto f    = make_vector_testing_data(msh);
        auto divf = make_vector_testing_data_div(msh);

        typename disk::hho_degree_info hdi(degree);

        scalar_type error = 0.0;
        for (auto& cl : msh)
        {
            Matrix<scalar_type, Dynamic, 1> proj = disk::project_function(msh, cl, hdi, f, 2);
            const auto                      dr   = make_hho_divergence_reconstruction(msh, cl, hdi);

            Matrix<scalar_type, Dynamic, 1> div = dr.first * proj;

            const auto                            cb   = disk::make_scalar_monomial_basis(msh, cl, hdi.face_degree());
            Matrix<scalar_type, Dynamic, Dynamic> mass = disk::make_mass_matrix(msh, cl, cb);
            Matrix<scalar_type, Dynamic, 1>       rhs  = disk::make_rhs(msh, cl, cb, divf);
            Matrix<scalar_type, Dynamic, 1>       exp_div = mass.llt().solve(rhs);

            Matrix<scalar_type, Dynamic, 1> diff = div - exp_div;

            error += diff.dot(mass * diff);
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
struct test_functor_mixed_order1
{
    /* Expect k+1 convergence */
    typename Mesh::coordinate_type
    operator()(const Mesh& msh, size_t degree) const
    {
        typedef Mesh                                mesh_type;
        typedef typename mesh_type::cell            cell_type;
        typedef typename mesh_type::face            face_type;
        typedef typename mesh_type::coordinate_type scalar_type;
        typedef typename mesh_type::point_type      point_type;

        typedef Matrix<scalar_type, Mesh::dimension, 1> ret_type;

        auto f    = make_vector_testing_data(msh);
        auto divf = make_vector_testing_data_div(msh);

        typename disk::hho_degree_info hdi(degree + 1, degree);

        scalar_type error = 0.0;
        for (auto& cl : msh)
        {
            Matrix<scalar_type, Dynamic, 1> proj = disk::project_function(msh, cl, hdi, f, 2);
            const auto                      dr   = make_hho_divergence_reconstruction(msh, cl, hdi);

            Matrix<scalar_type, Dynamic, 1> div = dr.first * proj;

            const auto                            cb   = disk::make_scalar_monomial_basis(msh, cl, hdi.face_degree());
            Matrix<scalar_type, Dynamic, Dynamic> mass = disk::make_mass_matrix(msh, cl, cb);
            Matrix<scalar_type, Dynamic, 1>       rhs  = disk::make_rhs(msh, cl, cb, divf);
            Matrix<scalar_type, Dynamic, 1>       exp_div = mass.llt().solve(rhs);

            Matrix<scalar_type, Dynamic, 1> diff = div - exp_div;

            error += diff.dot(mass * diff);
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
struct test_functor_mixed_order2
{
    /* Expect k+1 convergence */
    typename Mesh::coordinate_type
    operator()(const Mesh& msh, size_t degree) const
    {
        typedef Mesh                                mesh_type;
        typedef typename mesh_type::cell            cell_type;
        typedef typename mesh_type::face            face_type;
        typedef typename mesh_type::coordinate_type scalar_type;
        typedef typename mesh_type::point_type      point_type;

        typedef Matrix<scalar_type, Mesh::dimension, 1> ret_type;

        auto f    = make_vector_testing_data(msh);
        auto divf = make_vector_testing_data_div(msh);

        typename disk::hho_degree_info hdi(degree - 1, degree);

        scalar_type error = 0.0;
        for (auto& cl : msh)
        {
            Matrix<scalar_type, Dynamic, 1> proj = disk::project_function(msh, cl, hdi, f, 2);
            const auto                      dr   = make_hho_divergence_reconstruction(msh, cl, hdi);

            Matrix<scalar_type, Dynamic, 1> div = dr.first * proj;

            const auto                            cb   = disk::make_scalar_monomial_basis(msh, cl, hdi.face_degree());
            Matrix<scalar_type, Dynamic, Dynamic> mass = disk::make_mass_matrix(msh, cl, cb);
            Matrix<scalar_type, Dynamic, 1>       rhs  = disk::make_rhs(msh, cl, cb, divf);
            Matrix<scalar_type, Dynamic, 1>       exp_div = mass.llt().solve(rhs);

            Matrix<scalar_type, Dynamic, 1> diff = div - exp_div;

            error += diff.dot(mass * diff);
        }

        return std::sqrt(error);
    }

    size_t
    expected_rate(size_t k)
    {
        return k + 1;
    }
};

int
main(void)
{
    // face order: k, cell order: k
    std::cout << blue << "Face order: k and Cell order: k" << std::endl;
    tester<test_functor_equal_order> tstr1;
    tstr1.run();
    // face order: k, cell order: k+1
    std::cout << blue << "Face order: k and Cell order: k+1" << std::endl;
    tester<test_functor_mixed_order1> tstr2;
    tstr2.run();
    // face order: k, cell order: k-1
    std::cout << blue << "Face order: k and Cell order: k-1" << std::endl;
    tester<test_functor_mixed_order2> tstr3;
    tstr3.run(1, 3);
    return 0;
}
