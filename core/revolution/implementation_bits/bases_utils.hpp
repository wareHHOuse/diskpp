/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
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

#include "../quadratures"

namespace revolution {

template<typename Mesh, typename Element, typename Basis>
Matrix<typename Basis::scalar_type, Dynamic, Dynamic>
make_mass_matrix(const Mesh& msh, const Element& elem, const Basis& basis, size_t di = 0)
{
	auto degree		= basis.degree();
    auto basis_size = basis.size();

    using T = typename Basis::scalar_type;

    Matrix<T, Dynamic, Dynamic> ret = Matrix<T, Dynamic, Dynamic>::Zero(basis_size, basis_size);

    auto qps = integrate(msh, elem, 2*(degree+di));

    for (auto& qp : qps)
    {
        auto phi = basis.eval_functions(qp.point());
        ret += qp.weight() * phi * phi.transpose();
    }

    return ret;
}

template<typename Mesh, typename Element, typename Basis>
Matrix<typename Basis::scalar_type, Dynamic, Dynamic>
make_stiffness_matrix(const Mesh& msh, const Element& elem, const Basis& basis, size_t di = 0)
{
	auto degree		= basis.degree();
    auto basis_size = basis.size();

    using T = typename Basis::scalar_type;

    Matrix<T, Dynamic, Dynamic> ret = Matrix<T, Dynamic, Dynamic>::Zero(basis_size, basis_size);

    auto qps = integrate(msh, elem, 2*(degree+di));

    for (auto& qp : qps)
    {
        auto dphi = basis.eval_gradients(qp.point());
        ret += qp.weight() * dphi * dphi.transpose();
    }

    return ret;
}

template<typename Mesh, typename Element, typename Basis, typename Function>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
make_rhs(const Mesh& msh, const Element& elem,
         const Basis& basis, const Function& f, size_t di = 0)
{
    auto degree		= basis.degree();
    auto basis_size = basis.size();

    using T = typename Basis::scalar_type;

    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(basis_size);

    auto qps = integrate(msh, elem, 2*(degree+di));

    for (auto& qp : qps)
    {
        auto phi = basis.eval_functions(qp.point());
        ret += qp.weight() * phi * f(qp.point());
    }

    return ret;
}


} //revolution



