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


#include <iostream>
#include <iomanip>
#include <regex>

#include <unistd.h>

#include "contrib/colormanip.h"
#include <xmmintrin.h>

#include "bases/bases.hpp"
#include "quadratures/quadratures.hpp"
#include "methods/hho"

#include "core/loaders/loader.hpp"

#include "common.hpp"

template<typename Mesh>
struct test_functor_hho
{
    /* Expect k+1 convergence (hho stabilization) */
    typename Mesh::coordinate_type
    operator()(const Mesh& msh, size_t degree) const
    {
        typedef Mesh mesh_type;
        typedef typename mesh_type::cell        cell_type;
        typedef typename mesh_type::face        face_type;
        typedef typename mesh_type::coordinate_type scalar_type;
        typedef typename mesh_type::point_type  point_type;


        auto f = make_scalar_testing_data(msh);

        typename disk::hho_degree_info hdi(degree);

        scalar_type error = 0.0;
        for (auto& cl : msh)
        {
            auto gr = disk::make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab = disk::make_hho_scalar_stabilization(msh, cl, gr.first, hdi);

            Matrix<scalar_type, Dynamic, 1> proj = disk::project_function(msh, cl, hdi, f, 2);

            error += proj.dot(stab*proj);
        }

        return std::sqrt( error );
    }

    size_t
    expected_rate(size_t k)
    {
        return k+1;
    }
};

template<typename Mesh>
struct test_functor_hdg
{
    /* Expect k+1 convergence (hdg stabilization) */
    typename Mesh::coordinate_type
    operator()(const Mesh& msh, size_t degree) const
    {
        typedef Mesh                                mesh_type;
        typedef typename mesh_type::cell            cell_type;
        typedef typename mesh_type::face            face_type;
        typedef typename mesh_type::coordinate_type scalar_type;
        typedef typename mesh_type::point_type      point_type;

        auto f = make_scalar_testing_data(msh);

        typename disk::hho_degree_info hdi(degree+1, degree);

        scalar_type error = 0.0;
        for (auto& cl : msh)
        {
            auto stab = disk::make_hdg_scalar_stabilization(msh, cl, hdi);

            Matrix<scalar_type, Dynamic, 1> proj = disk::project_function(msh, cl, hdi, f, 2);

            error += proj.dot(stab * proj);
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
struct test_functor_dg
{
    /* Expect k convergence (dg stabilization) */
    typename Mesh::coordinate_type
    operator()(const Mesh& msh, size_t degree) const
    {
        typedef Mesh                                mesh_type;
        typedef typename mesh_type::cell            cell_type;
        typedef typename mesh_type::face            face_type;
        typedef typename mesh_type::coordinate_type scalar_type;
        typedef typename mesh_type::point_type      point_type;

        auto f = make_scalar_testing_data(msh);

        typename disk::hho_degree_info hdi(degree);

        scalar_type error = 0.0;
        for (auto& cl : msh)
        {
            auto stab = disk::make_dg_scalar_stabilization(msh, cl, hdi);

            Matrix<scalar_type, Dynamic, 1> proj = disk::project_function(msh, cl, hdi, f, 2);

            error += proj.dot(stab * proj);
        }

        return std::sqrt(error);
    }

    size_t
    expected_rate(size_t k)
    {
        return k;
    }
};

int main(void)
{
    std::cout << red << "Test HHO-stabilization operator" << std::endl;
    tester<test_functor_hho> tstr_hho;
    tstr_hho.run();

    std::cout << red << "Test HDG-stabilization operator" << std::endl;
    tester<test_functor_hdg> tstr_hdg;
    tstr_hdg.run();

    std::cout << red << "Test dG-stabilization operator" << std::endl;
    tester<test_functor_dg> tstr_dg;
    tstr_dg.run();
    return 0;
}
