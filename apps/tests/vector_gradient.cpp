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

#include "diskpp/bases/bases.hpp"
#include "diskpp/methods/hho"
#include "diskpp/quadratures/quadratures.hpp"
#include "diskpp/loaders/loader.hpp"

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

        auto f  = make_scalar_testing_data(msh);
        auto gf = make_scalar_testing_data_grad(msh);

        typename disk::hho_degree_info hdi(degree, degree, degree+1);

        scalar_type error = 0.0;
        for (auto& cl : msh)
        {
            Matrix<scalar_type, Dynamic, 1> proj = disk::project_function(msh, cl, hdi, f, 2);
            const auto                      gr   = disk::make_vector_hho_gradrec(msh, cl, hdi);

            Matrix<scalar_type, Dynamic, 1> gradc = gr.first * proj;
            Matrix<scalar_type, Dynamic, 1> grad = disk::project_function(msh, cl, hdi.grad_degree(), gf, 2);
            Matrix<scalar_type, Dynamic, 1> diff = grad - gradc;

            const auto                            gb   = disk::make_vector_monomial_basis(msh, cl, hdi.grad_degree());
            Matrix<scalar_type, Dynamic, Dynamic> mass = disk::make_mass_matrix(msh, cl, gb);

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

        auto f = make_vector_testing_data(msh);

        typename disk::hho_degree_info hdi(degree + 1, degree);

        scalar_type error = 0.0;
        for (auto& cl : msh)
        {
            Matrix<scalar_type, Dynamic, 1> proj = disk::project_function(msh, cl, hdi, f, 2);
            auto                            gr   = disk::make_vector_hho_laplacian(msh, cl, hdi);

            size_t rec_size = disk::vector_basis_size(hdi.reconstruction_degree(), Mesh::dimension, Mesh::dimension);

            Matrix<scalar_type, Dynamic, 1> reconstr  = Matrix<scalar_type, Dynamic, 1>::Zero(rec_size);
            reconstr.tail(rec_size - Mesh::dimension) = gr.first * proj;

            reconstr.head(Mesh::dimension) = proj.head(Mesh::dimension);

            auto cb = disk::make_vector_monomial_basis(msh, cl, hdi.reconstruction_degree());
            Matrix<scalar_type, Dynamic, Dynamic> mass         = disk::make_mass_matrix(msh, cl, cb);
            Matrix<scalar_type, Dynamic, 1>       rhs          = disk::make_rhs(msh, cl, cb, f);
            Matrix<scalar_type, Dynamic, 1>       exp_reconstr = mass.llt().solve(rhs);

            Matrix<scalar_type, Dynamic, 1> diff = reconstr - exp_reconstr;

            Matrix<scalar_type, Dynamic, Dynamic> stiffness = disk::make_stiffness_matrix(msh, cl, cb);

            error += diff.dot(stiffness * diff);
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

        auto f = make_vector_testing_data(msh);

        typename disk::hho_degree_info hdi(degree - 1, degree);

        scalar_type error = 0.0;
        for (auto& cl : msh)
        {
            Matrix<scalar_type, Dynamic, 1> proj = disk::project_function(msh, cl, hdi, f, 2);
            auto                            gr   = disk::make_vector_hho_laplacian(msh, cl, hdi);

            size_t rec_size = disk::vector_basis_size(hdi.reconstruction_degree(), Mesh::dimension, Mesh::dimension);

            Matrix<scalar_type, Dynamic, 1> reconstr  = Matrix<scalar_type, Dynamic, 1>::Zero(rec_size);
            reconstr.tail(rec_size - Mesh::dimension) = gr.first * proj;

            reconstr.head(Mesh::dimension) = proj.head(Mesh::dimension);

            auto cb = disk::make_vector_monomial_basis(msh, cl, hdi.reconstruction_degree());
            Matrix<scalar_type, Dynamic, Dynamic> mass         = disk::make_mass_matrix(msh, cl, cb);
            Matrix<scalar_type, Dynamic, 1>       rhs          = disk::make_rhs(msh, cl, cb, f);
            Matrix<scalar_type, Dynamic, 1>       exp_reconstr = mass.llt().solve(rhs);

            Matrix<scalar_type, Dynamic, 1> diff = reconstr - exp_reconstr;

            Matrix<scalar_type, Dynamic, Dynamic> stiffness = disk::make_stiffness_matrix(msh, cl, cb);

            error += diff.dot(stiffness * diff);
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
