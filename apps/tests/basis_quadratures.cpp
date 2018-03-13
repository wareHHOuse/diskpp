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

#include "revolution/bases"
#include "revolution/quadratures"
//#include "revolution/methods/hho"

#include "core/loaders/loader.hpp"

int main(void)
{
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    using T = double;


    typedef disk::simplicial_mesh<T, 3>  mesh_type;


    mesh_type msh;
    disk::netgen_mesh_loader<T, 3> loader;

    if (!loader.read_mesh("../../../diskpp/meshes/3D_tetras/netgen/basiccube.mesh"))
    {
        std::cout << "Problem loading mesh." << std::endl;
        return 1;
    }
    loader.populate_mesh(msh);

    for (auto& cl : msh)
    {
    	auto cb = revolution::make_scalar_monomial_basis(msh, cl, 2);
    	auto qps = revolution::integrate(msh, cl, 4);

        Matrix<T, Dynamic, 1> phi = Matrix<T, Dynamic, 1>::Zero(cb.size());

        for (auto& qp : qps)
            phi += qp.weight() * cb.eval_functions(qp.point());

        std::cout << phi.transpose() << std::endl;

    	auto fcs = faces(msh, cl);
    	for (auto& fc : fcs)
    	{
    		auto fb = revolution::make_scalar_monomial_basis(msh, fc, 2);
    		auto qps = revolution::integrate(msh, fc, 4);
    	}
    }

    return 0;
}
