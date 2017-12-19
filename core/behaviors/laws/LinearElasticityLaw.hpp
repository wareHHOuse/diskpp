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

#pragma once

#include "behaviors/maths_tensor.hpp"
#include "common/eigen.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

namespace disk {

/* Law for linear elasticity in small deformations
Input : symetric stain tensor (Gs)

Stress: sigma = 2*mu *Gs + lambda * trace(Gs) * Id
Tangent Moduli: C = 2*mu * I4 + lambda * prod_Kronecker(Id , Id)

*/

template<typename scalar_type>
class LinearElasticityLaw
{
   scalar_type m_lambda;
   scalar_type m_mu;

 public:
   LinearElasticityLaw() : m_lambda(1.0), m_mu(1.0) {}

   LinearElasticityLaw(const scalar_type lambda, const scalar_type mu) : m_lambda(lambda), m_mu(mu)
   {}

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
      const static_matrix<scalar_type, DIM, DIM> Id =
        static_matrix<scalar_type, DIM, DIM>::Identity();

      return 2 * m_mu * Gs + m_lambda * Gs.trace() * Id;
   }

   template<int DIM>
   static_tensor<scalar_type, DIM>
   compute_tangent_moduli() const
   {
      return 2 * m_mu * compute_IdentityTensor<scalar_type, DIM>() +
             m_lambda * compute_IxI<scalar_type, DIM>();
   }

   template<int DIM>
   std::pair<static_matrix<scalar_type, DIM, DIM>, static_tensor<scalar_type, DIM>>
   compute_whole(const static_matrix<scalar_type, DIM, DIM>& Gs) const
   {
      const static_matrix<scalar_type, DIM, DIM> sigma = compute_stress(Gs);
      const static_tensor<scalar_type, DIM>      C     = compute_tangent_moduli<DIM>();

      return std::make_pair(sigma, C);
   }
};
}
