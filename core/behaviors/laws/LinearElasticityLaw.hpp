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

#include "common/eigen.hpp"
#include "BehaviorLaws/maths_tensor.hpp"


#define _USE_MATH_DEFINES
#include <cmath>



template<typename scalar_type>
class LinearElasticityLaw
{
   scalar_type m_lambda;


public:
   LinearElasticityLaw()
   : m_lambda(0.0)
   {}

   LinearElasticityLaw(const scalar_type lambda)
   :m_lambda(lambda)
   {}

   void
   setLambda(const scalar_type lambda)
   {
      m_lambda = lambda;}

   scalar_type
   giveLambda() const {return m_lambda;}


   template<int DIM>
   static_matrix<scalar_type, DIM, DIM>
   compute_PK1(const static_matrix<scalar_type, DIM, DIM>& F) const
   {
      const static_matrix<scalar_type, DIM, DIM> Id = static_matrix<scalar_type, DIM , DIM>::Identity();

      return m_lambda * (F- Id);
   }

   template<int DIM>
   static_tensor<scalar_type, DIM>
   compute_tangent_moduli(const static_matrix<scalar_type,DIM,DIM>& F) const
   {
      return m_lambda * compute_IdentityTensor<scalar_type, DIM>();
   }

   template<int DIM>
   std::pair<static_matrix<scalar_type, DIM, DIM>, static_tensor<scalar_type, DIM> >
   compute_whole_PK1(const static_matrix<scalar_type, DIM, DIM>& F) const
   {
      const static_matrix<scalar_type, DIM, DIM> PK1 = compute_PK1(F);
      const static_tensor<scalar_type, DIM> A = compute_tangent_moduli(F);

      return std::make_pair(PK1, A);
   }

};
