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
 * Karol Cascavita (C) 2018                     karol.cascavita@enpc.fr
 * Nicolas Pignet  (C) 2018, 2019               nicolas.pignet@enpc.fr
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

#pragma once

#include "common/eigen.hpp"
#include "quadratures/quadratures.hpp"

namespace disk {

namespace priv
{
// inner_product : generalisation of dot_product, scalar_product in l^p, matrix-vector product
// multiplication by a scalar

// dot-product in R x R
template<typename T>
T
inner_product(const T& a, const T& b)
{
    return a * b;
}

// dot-product in R^n x R^n
template<typename T, int N>
T
inner_product(const Matrix<T, N, 1>& a, const Matrix<T, N, 1>& b)
{
    return a.dot(b);
}

// dot-product in R^(nxm) x R^(nxm)
template<typename T, int N, int M>
T
inner_product(const Matrix<T, N, M>& b, const Matrix<T, N, M>& a)
{
    return a.cwiseProduct(b).sum();
}

// dot-product of a list of vector a with a vector b
template<typename T, int N>
Matrix<T, Dynamic, 1>
inner_product(const Matrix<T, Dynamic, N>& a, const Matrix<T, N, 1>& b)
{
    return a * b;
}

// matrix-vector product
template<typename T, int N>
Matrix<T, N, 1>
inner_product(const Matrix<T, N, N>& a, const Matrix<T, N, 1>& b)
{

    return a * b;
}

// matrix-vector product
template<typename T, int N>
Matrix<T, Dynamic, N>
inner_product(const Matrix<T, N, N> a, const Matrix<T, Dynamic, N>& b)
{
    Matrix<T, Dynamic, N> ret(b.rows(), N);

    for (size_t i = 0; i < b.rows(); i++)
    {
        ret.row(i) = (a * (b.row(i).transpose())).transpose();
    }

    return ret;
}

// matrix-vector product of a list of matrix a with a vector b
template<typename T, int N>
Matrix<T, Dynamic, N>
inner_product(const eigen_compatible_stdvector<Matrix<T, N, N>>& a, const Matrix<T, N, 1>& b)
{
    Matrix<T, Dynamic, N> ret(a.size(), N);

    for (size_t i = 0; i < a.size(); i++)
        ret.row(i) = inner_product(a[i], b).transpose();

    return ret;
}

// multiplication by a scalar
template<typename T, int N>
Matrix<T, Dynamic, N>
inner_product(const T& a, const Matrix<T, Dynamic, N>& b)
{
    return a * b;
}

// multiplication by a scalar
template<typename T, int N, int M>
Matrix<T, N, M>
inner_product(const T& a, const Matrix<T, N, M>& b)
{
    return a * b;
}

// multiplication by a scalar
template<typename T, int N, int M>
eigen_compatible_stdvector<Matrix<T, N, M>>
inner_product(const T& a, const eigen_compatible_stdvector<Matrix<T, N, M>>& b)
{

    eigen_compatible_stdvector<Matrix<T, N, M>> ret = b;

    for (size_t i = 0; i < b.size(); i++)
        ret[i] *= a;

    return ret;
}

// outer_product : equivalent to a scalar-product in L^p (.,.)_L^p or duality product <.,.>
template<typename T, int N>
Matrix<T, Dynamic, Dynamic>
outer_product(const Matrix<T, Dynamic, N>& a, const Matrix<T, Dynamic, N>& b)
{
    return a * b.transpose();
}

template<typename T, int N, int M>
Matrix<T, Dynamic, Dynamic>
outer_product(const eigen_compatible_stdvector<Matrix<T, N, M>>& a,
              const eigen_compatible_stdvector<Matrix<T, N, M>>& b)
{
    Matrix<T, Dynamic, Dynamic> ret(a.size(), b.size());

    for (size_t i = 0; i < a.size(); i++)
        for (size_t j = 0; j < b.size(); j++)
            ret(i, j) = inner_product(a[i], b[j]);

    return ret;
}

template<typename T>
T
outer_product(const T& a, const T& b)
{
    return a * b;
}

// outer_product (b is a constant in this case) return a vector
template<typename T>
Matrix<T, Dynamic, 1>
outer_product(const Matrix<T, Dynamic, 1>& a, const T& b)
{
    return a * b;
}

// outer_product (b is a constant in this case) return a vector
template<typename T, int N>
Matrix<T, Dynamic, 1>
outer_product(const Matrix<T, Dynamic, N>& a, const Matrix<T, N, 1>& b)
{
    return a * b;
}

// outer_product (b is a constant in this case) return a vector
template<typename T, int N>
Matrix<T, Dynamic, 1>
outer_product(const eigen_compatible_stdvector<Matrix<T, N, N>>& a, const Matrix<T, N, N>& b)
{
    Matrix<T, Dynamic, 1> ret(a.size());

    for (size_t i = 0; i < a.size(); i++)
        ret(i) = inner_product(a[i], b);

    return ret;
}

} // priv


template<typename Mesh, typename Element, typename Basis>
Matrix<typename Basis::scalar_type, Dynamic, Dynamic>
make_mass_matrix(const Mesh& msh, const Element& elem, const Basis& basis, size_t di = 0)
{
    const auto degree     = basis.degree();
    const auto basis_size = basis.size();

    using T = typename Basis::scalar_type;

    Matrix<T, Dynamic, Dynamic> ret = Matrix<T, Dynamic, Dynamic>::Zero(basis_size, basis_size);

    const auto qps = integrate(msh, elem, 2 * (degree+di));

    for (auto& qp : qps)
    {
        const auto phi    = basis.eval_functions(qp.point());
        const auto qp_phi = priv::inner_product(qp.weight(), phi);
        ret += priv::outer_product(qp_phi, phi);
    }

    return ret;
}

/* We have a problem here: this definition could be ambiguous with the previous
 * one. */
//#if 0
template<typename Mesh, typename Element, typename Basis, typename MaterialField>
[[deprecated("DiSk++ issue: this declaration is ambiguous, fix is needed")]]
Matrix<typename Basis::scalar_type, Dynamic, Dynamic>
make_mass_matrix(const Mesh& msh, const Element& elem, const Basis& basis, const MaterialField& material_tensor, size_t di = 0)
{
    const auto degree     = basis.degree();
    const auto basis_size = basis.size();

    using T = typename Basis::scalar_type;

    Matrix<T, Dynamic, Dynamic> ret = Matrix<T, Dynamic, Dynamic>::Zero(basis_size, basis_size);

    const auto qps = integrate(msh, elem, 2 * (degree + di));

    for (auto& qp : qps)
    {
        const auto phi      = basis.eval_functions(qp.point());
        const auto qp_m     = priv::inner_product(qp.weight(), material_tensor(qp.point()));
        const auto qp_phi_m = priv::inner_product(qp_m, phi);
        ret += priv::outer_product(phi, qp_phi_m);
    }

    return ret;
}
//#endif

template<typename Mesh, typename Element, typename Basis>
Matrix<typename Basis::scalar_type, Dynamic, Dynamic>
make_stiffness_matrix(const Mesh& msh, const Element& elem, const Basis& basis)
{
    const auto degree     = basis.degree();
    const auto basis_size = basis.size();

    using T = typename Basis::scalar_type;

    Matrix<T, Dynamic, Dynamic> ret = Matrix<T, Dynamic, Dynamic>::Zero(basis_size, basis_size);

    if(degree > 0)
    {
        const auto qps = integrate(msh, elem, 2 * (degree - 1));

        for (auto& qp : qps)
        {
            const auto dphi    = basis.eval_gradients(qp.point());
            const auto qp_dphi = priv::inner_product(qp.weight(), dphi);
            ret += priv::outer_product(qp_dphi, dphi);
        }
    }


    return ret;
}

template<typename Mesh, typename Element, typename Basis, typename Function>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
make_rhs(const Mesh& msh, const Element& elem, const Basis& basis, const Function& rhs_fun, size_t di = 0)
{
    const auto degree     = basis.degree();
    const auto basis_size = basis.size();

    using T = typename Basis::scalar_type;

    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(basis_size);

    const auto qps = integrate(msh, elem, 2 * (degree + di));

    for (auto& qp : qps)
    {
        const auto phi  = basis.eval_functions(qp.point());
        const auto qp_f = priv::inner_product(qp.weight(), rhs_fun(qp.point()));
        ret += priv::outer_product(phi, qp_f);
    }

    return ret;
}

template<typename Mesh>
using SRT = typename Mesh::coordinate_type;

template<typename Mesh>
using VRT = Matrix<typename Mesh::coordinate_type, Mesh::dimension, 1>;

template<typename Mesh>
using MRT = Matrix<typename Mesh::coordinate_type, Mesh::dimension, Mesh::dimension>;

template<typename Mesh>
using PT = typename Mesh::point_type;

template<typename Mesh>
using scalar_rhs_function = std::function<SRT<Mesh>(PT<Mesh>)>;

template<typename Mesh>
using vector_rhs_function = std::function<VRT<Mesh>(PT<Mesh>)>;

template<typename Mesh>
using matrix_rhs_function = std::function<MRT<Mesh>(PT<Mesh>)>;

template<typename Mesh, typename Element>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
project_function(const Mesh& msh, const Element& elem, const size_t degree, const scalar_rhs_function<Mesh>& f, size_t di = 0)
{
    using T = typename Mesh::coordinate_type;

    const auto                  basis       = make_scalar_monomial_basis(msh, elem, degree);
    const Matrix<T, Dynamic, Dynamic> mass  = make_mass_matrix(msh, elem, basis);
    const Matrix<T, Dynamic, 1>       rhs   = make_rhs(msh, elem, basis, f, di);

    return mass.ldlt().solve(rhs);
}

template<typename Mesh, typename Element, typename Basis>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
project_function(const Mesh&                      msh,
                 const Element&                   elem,
                 const Basis&                     basis,
                 const scalar_rhs_function<Mesh>& f,
                 size_t                           di = 0)
{
    using T = typename Mesh::coordinate_type;

    const Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, elem, basis);
    const Matrix<T, Dynamic, 1>       rhs  = make_rhs(msh, elem, basis, f, di);

    return mass.ldlt().solve(rhs);
}

template<typename Mesh, typename Element>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
project_function(const Mesh& msh, const Element& elem, const size_t degree, const vector_rhs_function<Mesh>& f, size_t di = 0)
{
    using T = typename Mesh::coordinate_type;

    const auto                        basis = make_vector_monomial_basis(msh, elem, degree);
    const Matrix<T, Dynamic, Dynamic> mass  = make_mass_matrix(msh, elem, basis);
    const Matrix<T, Dynamic, 1>       rhs   = make_rhs(msh, elem, basis, f, di);

    return mass.ldlt().solve(rhs);
}

template<typename Mesh, typename Element, typename Basis>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
project_function(const Mesh&    msh,
                 const Element& elem,
                 const Basis& basis,
                 const vector_rhs_function<Mesh>& f,
                 size_t                           di = 0)
{
    using T = typename Mesh::coordinate_type;

    const Matrix<T, Dynamic, Dynamic> mass  = make_mass_matrix(msh, elem, basis);
    const Matrix<T, Dynamic, 1>       rhs   = make_rhs(msh, elem, basis, f, di);

    return mass.ldlt().solve(rhs);
}

template<typename Mesh, typename Element>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
project_gradient(const Mesh& msh, const Element& elem, size_t degree, const matrix_rhs_function<Mesh>& f, size_t di = 0)
{
    using T = typename Mesh::coordinate_type;

    const auto basis = make_vector_monomial_basis(msh, elem, degree);
    const size_t N     = Mesh::dimension;

    const auto basis_size = vector_basis_size(degree, Mesh::dimension, Mesh::dimension);

    const Matrix<T, Dynamic, Dynamic> stiff = make_stiffness_matrix(msh, elem, basis);
    Matrix<T, Dynamic, 1>             rhs   = Matrix<T, Dynamic, 1>::Zero(basis_size);

    const auto qps = integrate(msh, elem, 2 * (degree + di));

    for (auto& qp : qps)
    {
        const auto dphi        = basis.eval_gradients(qp.point());
        const auto qp_rhs_func = priv::inner_product(qp.weight(), f(qp.point()));

        rhs += priv::outer_product(dphi, qp_rhs_func);
    }

    return stiff.block(N, N, basis_size - N, basis_size - N).ldlt().solve(rhs.tail(basis_size - N));
}

template<typename T>
T
eval(const dynamic_vector<T>& tab_coeff, const dynamic_vector<T>& base)
{
    assert(tab_coeff.rows() == base.rows());

    return tab_coeff.dot(base);
}

template<typename T, int DIM>
static_vector<T, DIM>
eval(const dynamic_vector<T>& tab_coeff, const Matrix<T, Dynamic, DIM>& base)
{
    static_vector<T, DIM> ret = static_vector<T, DIM>::Zero();
    assert(tab_coeff.size() == (base.rows()));

    for (int i = 0; i < tab_coeff.size(); i++)
    {
        ret += tab_coeff(i) * (base.row(i)).transpose();
    }

    return ret;
}

template<typename T, int DIM>
static_matrix<T, DIM, DIM>
eval(const dynamic_vector<T>&                                      tab_coeff,
     const eigen_compatible_stdvector<static_matrix<T, DIM, DIM>>& base,
     const size_t                                                  begin = 0)
{
    static_matrix<T, DIM, DIM> ret = static_matrix<T, DIM, DIM>::Zero();
    assert(tab_coeff.size() == (base.size() - begin));

    for (int i = 0; i < tab_coeff.size(); i++)
    {
        ret += tab_coeff(i) * base[i + begin];
    }

    return ret;
}

} // namespace disk
