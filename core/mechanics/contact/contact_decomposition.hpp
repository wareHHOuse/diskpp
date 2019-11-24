/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018, 2019
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework.
 * M. Abbas, A. Ern, N. Pignet.
 * International Journal of Numerical Methods in Engineering (2019)
 * 120(3), 303-327
 * DOI: 10.1002/nme.6137
 */

#pragma once

#include "bases/bases.hpp"
#include "common/eigen.hpp"

namespace disk
{
namespace mechanics
{

template<typename T, size_t DIM>
T
normal_part(const static_vector<T, DIM>& vec, const static_vector<T, DIM>& n)
{
    return vec.dot(n);
}

template<typename T, size_t DIM>
static_vector<T, DIM>
tangential_part(const static_vector<T, DIM>& vec, const static_vector<T, DIM>& n)
{
    return vec - vec.dot(n) * n;
}

template<typename Basis, typename T, size_t DIM>
dynamic_vector<T>
normal_part(const Basis& basis, const static_vector<T, DIM>& n, const point<T,DIM>& pt)
{
    const auto basis_pt = basis.eval_functions(pt);

    //(basis . n)
    return priv::inner_product(basis_pt, n);
}

template<typename Basis, typename T, size_t DIM>
Eigen::Matrix<T, Eigen::Dynamic, DIM>
tangential_part(const Basis& basis, const static_vector<T, DIM>& n, const point<T, DIM>& pt)
{
    const auto basis_pt = basis.eval_functions(pt);
    const auto u_n   = normal_part(basis, n, pt);

    // phi_T - (phi_T . n)n
    return basis_pt - priv::inner_product(u_n, n);
}

template<typename T, size_t DIM>
static_vector<T, DIM>
stress_n(const static_matrix<T, DIM, DIM>& stress, const static_vector<T, DIM>& n)
{
    return stress * n;
}

template<typename T, size_t DIM>
T
stress_nn(const static_matrix<T, DIM, DIM>& stress, const static_vector<T, DIM>& n)
{
    return normal_part(stress * n, n);
}

template<typename T, size_t DIM>
T
stress_nn(const static_vector<T, DIM>& stress_n, const static_vector<T, DIM>& n)
{
    return normal_part(stress_n, n);
}

template<typename T, size_t DIM>
static_vector<T, DIM>
stress_nt(const static_matrix<T, DIM, DIM>& stress, const static_vector<T, DIM>& n)
{
    const auto stress_n = stress_n(stress, n);
    return tangential_part(stress_n, n);
}

template<typename T, size_t DIM>
static_vector<T, DIM>
stress_nt(const static_vector<T, DIM>& stress_n, const static_vector<T, DIM>& n)
{
    return tangential_part(stress_n, n);
}

template<typename StressBasis, typename T>
auto
stress(const StressBasis& sb, const dynamic_matrix<T>& stress_coeff, const point_type pt) const
{
    const auto sphi = sb.eval_functions(pt);

    return eval(stress_coeff, sphi);
}


template<typename StressBasis, typename T, size_t DIM>
static_vector<T, DIM>
stress_n(const StressBasis&             sb,
           const dynamic_matrix<T>&     coeff_stress,
           const static_vector<T, DIM>& n,
           const point<T, DIM>&         pt
{
    const auto Stress = stress(sb, coeff_stress, pt);

    // sigma_n . n
    return stress_n(Stress, n);
}

template<typename StressBasis, typename T, size_t DIM>
T
stress_nn(const StressBasis&             sb,
           const dynamic_matrix<T>&     coeff_stress,
           const static_vector<T, DIM>& n,
           const point<T, DIM>&         pt
{
    const auto Stress = stress(sb, coeff_stress, pt);

    // sigma_n . n
    return stress_nn(Stress, n);
}

template<typename T, size_t DIM>
T
u_n(const static_vector<T, DIM>& u, const static_vector<T, DIM>& n)
{
    return normal_part(u, n);
}

template<typename T, size_t DIM>
static_vector<T, DIM>
u_t(const static_vector<T, DIM>& u, const static_vector<T, DIM>& n)
{
    return tangential_part(u, n);
}

template<typename Basis, typename T, size_t DIM>
T
u_n(const Basis& basis, const dynamic_vector<T>& coeff_u, const static_vector<T, DIM>& n, const point<T,DIM>& pt)
{
    const auto basis_pt = basis.eval_function(pt);
    const auto u        = basis_pt.dot(coeff_u);
    return normal_part(u, n);
}

template<typename Basis, typename T, size_t DIM>
static_vector<T, DIM>
u_t(const Basis& basis, const dynamic_vector<T>& coeff_u, const static_vector<T, DIM>& n, const point<T, DIM>& pt)
{
    const auto basis_pt = basis.eval_function(pt);
    const auto u        = basis_pt.dot(coeff_u);
    return tangential_part(u, n);
}

template<typename Basis, typename T, size_t DIM>
dynamic_vector<T>
v_n(const Basis& basis, const static_vector<T, DIM>& n, const point<T, DIM>& pt)
{
    return normal_part(basis, n, pt);
}

template<typename Basis, typename T, size_t DIM>
Eigen::Matrix<T, Eigen::Dynamic, DIM>
v_t(const Basis& basis, const static_vector<T, DIM>& n, const point<T, DIM>& pt)
{
    return tangential_part(basis, n, pt);
}

// cauchy traction : sigma_n = sigma * n
template<typename GradBasis, typename T, size_t DIM, typename Material>
Eigen::Matrix<T, Eigen::Dynamic, DIM>
sigma_v_n(const GradBasis&             gb,
          const dynamic_matrix<T>&     ET,
          const static_vector<T, DIM>& n,
          const point<T, DIM>&         pt,
          const Material&              material)
{
    Eigen::Matrix<T, Eigen::Dynamic, DIM> sigma_n = Eigen::Matrix<T, Eigen::Dynamic, DIM>::Zero(ET.cols(), DIM);

    const auto gphi   = gb.eval_functions(pt);
    const auto gphi_n = normal_part(gb, n, n);

    sigma_n = 2.0 * material.getMu() * ET.transpose() * gphi_n;

    const auto gphi_trace_n = priv::inner_product(trace(gphi), n);

    sigma_n += material.getLambda() * ET.transpose() * gphi_trace_n;

    // sigma. n
    return sigma_n;
}

template<typename GradBasis>
dynamic_vector<T>
sigma_v_nn(const GradBasis&             gb,
           const dynamic_matrix<T>&     ET,
           const static_vector<T, DIM>& n,
           const point<T, DIM>&         pt,
           const Material&              material)
{
    const auto sigma_n = sigma_v_n(gb, ET, n, pt, material);

    // sigma_n . n
    return priv::inner_product(sigma_n, n);
}

// cauchy traction : sigma_n = sigma * n
template<typename GradBasis, typename T, size_t DIM, typename Material>
Eigen::Matrix<T, Eigen::Dynamic, DIM>
sigma_v_nt(const GradBasis&             gb,
           const dynamic_matrix<T>&     ET,
           const static_vector<T, DIM>& n,
           const point<T, DIM>&         pt,
           const Material&              material)
{
    const auto sigma_n = sigma_v_n(gb, ET, n, pt, material);
    const auto sigma_nn = sigma_v_nn(gb, ET, n, pt, material);

    // sigma_n - sigma_nn * n
    return sigma_n - priv::inner_product(sigma_nn, n);
}

} // end mechanics
} // end diskpp