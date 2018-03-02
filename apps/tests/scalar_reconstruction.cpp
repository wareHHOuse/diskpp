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

#include <xmmintrin.h>

#include "revolution/bases"
#include "revolution/quadratures"
#include "revolution/methods/hho"

#include "core/loaders/loader.hpp"

#include "common.hpp"


template<typename Mesh>
struct test_functor
{
    typename Mesh::scalar_type
    operator()(const Mesh& msh, size_t degree) const
    {
        typedef Mesh mesh_type;
        typedef typename mesh_type::cell        cell_type;
        typedef typename mesh_type::face        face_type;
        typedef typename mesh_type::scalar_type scalar_type;
        typedef typename mesh_type::point_type  point_type;


        auto f = make_scalar_testing_data(msh);

        typename revolution::hho_degree_info hdi(degree);

        scalar_type error = 0.0;
        for (auto& cl : msh)
        {
            Matrix<scalar_type, Dynamic, 1> proj = revolution::project_function(msh, cl, hdi, f);
            auto gr = revolution::make_hho_scalar_laplacian(msh, cl, hdi);

            size_t rec_size = revolution::scalar_basis_size(hdi.reconstruction_degree(), Mesh::dimension);

            Matrix<scalar_type, Dynamic, 1> reconstr = Matrix<scalar_type, Dynamic, 1>::Zero(rec_size);
            reconstr.tail(rec_size-1) = gr.first * proj;
            reconstr(0) = proj(0);

            auto cb = revolution::make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
            Matrix<scalar_type, Dynamic, Dynamic> mass = revolution::make_mass_matrix(msh, cl, cb);
            Matrix<scalar_type, Dynamic, 1> rhs = revolution::make_rhs(msh, cl, cb, f);
            Matrix<scalar_type, Dynamic, 1> exp_reconstr = mass.llt().solve(rhs);

            Matrix<scalar_type, Dynamic, 1> diff = reconstr - exp_reconstr;

            Matrix<scalar_type, Dynamic, Dynamic> stiffness = revolution::make_stiffness_matrix(msh, cl, cb);

            error += diff.dot(stiffness*diff);
        }

        return std::sqrt(error);
    }
};

template<typename Mesh>
test_functor<Mesh>
get_test_functor(const std::vector<Mesh>& meshes)
{
    return test_functor<Mesh>();
}

void test_triangles(void)
{
    std::cout << "*** TESTING TRIANGLES ON GENERIC MESH ***" << std::endl;
    using T = double;

    auto meshes = get_triangle_generic_meshes<T>();
    auto tf = get_test_functor(meshes);
    do_testing(meshes, tf);
}

void test_quads(void)
{
    std::cout << "*** TESTING QUADS ON GENERIC MESH ***" << std::endl;
    using T = double;

    auto meshes = get_quad_generic_meshes<T>();
    auto tf = get_test_functor(meshes);
    do_testing(meshes, tf);
}

int main(void)
{
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    
    test_triangles();
    test_quads();

    return 0;
}
