/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2018                     nicolas.pignet@enpc.fr
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
#include "mechanics/behaviors/maths_tensor.hpp"
#include "mechanics/behaviors/maths_utils.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

namespace disk
{

/* Law for Hencky-von Mises model in small deformations
Input : symetric stain tensor (Gs)

dev = normL2(Gs - trace(Gs)/dim * Id)

Stress: sigma = 2*\tilde{mu}(dev(Gs)) *Gs + \tilde{lambda}(dev(Gs)) * trace(Gs) * Id
\tilde{mu}(dev(Gs)) = mu*(1 +(1+ dev(Gs))^{-1/2})
\tilde{lambda}(dev(Gs)) = ((lambda + mu/2)-mu/2*(1+ dev(Gs))^{-1/2})

Tangent Moduli: C = 2*mu * I4 + lambda * prod_Kronecker(Id , Id) /it is the elastic moduli

*/

template<typename scalar_type>
class HenckyMisesLaw
{
    scalar_type m_lambda;
    scalar_type m_mu;

    scalar_type
    tildemu(const scalar_type& normL2dev) const
    {
        return m_mu * (1 + 1.0 / std::sqrt(1 + normL2dev * normL2dev));
    }

    scalar_type
    tildelambda(const scalar_type& normL2dev) const
    {
        return (m_lambda + m_mu / 2.0) - m_mu / (2.0 * std::sqrt(1 + normL2dev * normL2dev));
    }

    scalar_type
    derivativetildemu(const scalar_type& normL2dev) const
    {
        const auto term3_2 = (1 + normL2dev * normL2dev) * std::sqrt(1 + normL2dev * normL2dev);

        return -m_mu * normL2dev / term3_2;
    }

    scalar_type
    derivativetildelambda(const scalar_type& normL2dev) const
    {
        return derivativetildemu(normL2dev) / scalar_type(2);
    }

  public:
    HenckyMisesLaw() : m_lambda(1.0), m_mu(1.0) {}

    HenckyMisesLaw(const scalar_type lambda, const scalar_type mu) : m_lambda(lambda), m_mu(mu) {}

    void
    setLambda(const scalar_type lambda)
    {
        m_lambda = lambda;
    }

    scalar_type
    giveLambda() const
    {
        return m_lambda;
    }

    void
    setMu(const scalar_type mu)
    {
        m_mu = mu;
    }

    scalar_type
    giveMu() const
    {
        return m_mu;
    }

    template<int DIM>
    static_matrix<scalar_type, DIM, DIM>
    compute_stress(const static_matrix<scalar_type, DIM, DIM>& Gs) const
    {
        const static_matrix<scalar_type, DIM, DIM> Id = static_matrix<scalar_type, DIM, DIM>::Identity();

        const scalar_type normFrodev = deviator(Gs).norm();

        return tildemu(normFrodev) * Gs + tildelambda(normFrodev) * Gs.trace() * Id;
    }

    template<int DIM>
    static_tensor<scalar_type, DIM>
    compute_tangent_moduli(const static_matrix<scalar_type, DIM, DIM>& Gs) const
    {
        return 2 * m_mu * compute_IdentityTensor<scalar_type, DIM>() + m_lambda * compute_IxI<scalar_type, DIM>();

        // const auto dev        = deviator(Gs);
        // const auto normFrodev = dev.norm();

        // if (normFrodev == scalar_type(0)) {
        //    return 2 * m_mu * compute_IdentityTensor<scalar_type, DIM>() +
        //           m_lambda * compute_IxI<scalar_type, DIM>();
        // }

        // const auto devFirstDerivative = computeNormFroFirstDerivate(dev);
        // // const static_matrix<scalar_type, DIM, DIM> Id =
        // //   static_matrix<scalar_type, DIM, DIM>::Identity();

        // const auto term1 = tildemu(normFrodev) * compute_IdentityTensor<scalar_type, DIM>();
        // const auto term2 = tildelambda(normFrodev) * compute_IxI<scalar_type, DIM>();

        // const auto term3 =
        //   derivativetildemu(normFrodev) * computeKroneckerProduct(devFirstDerivative, Gs);
        // const auto term4 = Gs.trace() * derivativetildelambda(normFrodev) *
        //                    computeKroneckerProduct(devFirstDerivative, Id);

        // return term1 + term2 + term3 + term4;
    }

    template<int DIM>
    std::pair<static_matrix<scalar_type, DIM, DIM>, static_tensor<scalar_type, DIM>>
    compute_whole(const static_matrix<scalar_type, DIM, DIM>& Gs) const
    {
        const static_matrix<scalar_type, DIM, DIM> sigma = compute_stress(Gs);
        const static_tensor<scalar_type, DIM>      C     = compute_tangent_moduli<DIM>(Gs);

        return std::make_pair(sigma, C);
    }
};
}
