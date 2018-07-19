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

#include <xmmintrin.h>

#include "bases/bases.hpp"
#include "quadratures/quadratures.hpp"
#include "methods/hho"

#include "common.hpp"

template<typename Mesh>
struct test_functor
{
    /* Expect k+1 convergence on the cells and k+0.5 on the faces. */
    typename Mesh::scalar_type
    operator()(const Mesh& msh, size_t degree) const
    {
        typedef Mesh mesh_type;
        typedef typename mesh_type::cell        cell_type;
        typedef typename mesh_type::face        face_type;
        typedef typename mesh_type::scalar_type scalar_type;
        typedef typename mesh_type::point_type  point_type;


        auto f = make_scalar_testing_data(msh);

        scalar_type error = 0.0;

        for (auto itor = msh.cells_begin(); itor != msh.cells_end(); itor++)
        {
            auto cl = *itor;
            auto basis = disk::make_scalar_monomial_basis(msh, cl, degree);

            Matrix<scalar_type, Dynamic, 1> proj = disk::project_function(msh, cl, degree, f);

            auto qps = disk::integrate(msh, cl, 2*degree+2);
            for (auto& qp : qps)
            {
                auto tv = f(qp.point());

                auto phi = basis.eval_functions(qp.point());
                auto pv = proj.dot(phi);

                error += qp.weight() * (tv - pv) * (tv - pv);
            }
        }

        return std::sqrt( error );
    }

    size_t
    expected_rate(size_t k)
    {
        return k+1;
    }
};

int main(void)
{
    tester<test_functor> tstr;
    tstr.run();
    return 0;
}
