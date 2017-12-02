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

#include <iostream>

#include <sstream>
#include <string>
#include <math.h>

#include "common/eigen.hpp"
#include "BehaviorLaws/maths_tensor.hpp"
#include <Eigen/LU>

#define _USE_MATH_DEFINES
#include <cmath>

/* Material: Neo-nookean
 * Energy :  W(F) = Wiso(F) + Wvol(F)
 *   - Wiso(F) =    mu / 2 *[tr(F^T * F) - d]
 *   - Wvol(F) =  lambda/2 * U(J)**2 - mu * ln(J)
 * ** We set T1(J) = J * U(J) * U'(J) and T2(J) =  U(J) * J *( U''(J) * J + U'(J)) + ( J * U'(J))^2
 * Stress :  PK1(F) = Piso(F) + Pvol(F)
 *   - Piso(F) = mu * F
 *   - Pvol(F) = (lambda * T1(J) - mu) * F^{-T}
 * Module :  A(F) = PK1(F) = Aiso(F) + Avol(F)
 *   - Aiso(F) = mu *  I4
 *   - Avol(F) = (mu -lambda * T1(J)) * F^{-T} \time_inf F^{-1}
 *                  + lambda * T2(J) * F^{-T} \kronecker F^{-T}
 */


template<typename scalar_type>
class MoonleyLaw
{
   scalar_type m_mu1;
   scalar_type m_mu2;
   scalar_type m_lambda;


   static_matrix<scalar_type, 2, 2>
   PK3(const static_matrix<scalar_type, 2, 2>& F) const
   {
      static_matrix<scalar_type, 2, 2> P = static_matrix<scalar_type, 2, 2>::Zero();

      P(0,0) = 4.0 * F(0,0) * (F(0,0)*F(0,0) + F(0,1)*F(0,1) + F(1,0)*F(1,0) - 1.0)
               + 2.0 * F(0,1)*F(1,0)*F(1,1);

      P(0,1) = 4.0 * F(0,1) * (F(0,1)*F(0,1) + F(0,0)*F(0,0) + F(1,1)*F(1,1) - 1.0)
               + 2.0 * F(0,0)*F(1,0)*F(1,1);

      P(1,0) = 4.0 * F(1,0) * (F(1,0)*F(1,0) + F(0,0)*F(0,0) + F(1,1)*F(1,1) - 1.0)
               + 2.0 * F(0,0)*F(0,1)*F(1,1);

      P(1,1) = 4.0 * F(1,1) * (F(1,1)*F(1,1) + F(0,1)*F(0,1) + F(1,0)*F(1,0) - 1.0)
               + 2.0 * F(0,0)*F(0,1)*F(1,0);

      return  P;
   }


   static_matrix<scalar_type, 3, 3>
   PK3(const static_matrix<scalar_type, 3, 3>& F) const
   {
      static_matrix<scalar_type, 3, 3> P = static_matrix<scalar_type, 3, 3>::Zero();

      P(0,0) = 4.0 * F(0,0) * (F(0,0)*F(0,0) + F(0,1)*F(0,1) + F(1,0)*F(1,0) - 1.0)
               + 2.0 * F(0,1)*F(1,0)*F(1,1);

      P(0,1) = 4.0 * F(0,1) * (F(0,1)*F(0,1) + F(0,0)*F(0,0) + F(1,1)*F(1,1) - 1.0)
               + 2.0 * F(0,0)*F(1,0)*F(1,1);

      P(1,0) = 4.0 * F(1,0) * (F(1,0)*F(1,0) + F(0,0)*F(0,0) + F(1,1)*F(1,1) - 1.0)
               + 2.0 * F(0,0)*F(0,1)*F(1,1);

      P(1,1) = 4.0 * F(1,1) * (F(1,1)*F(1,1) + F(0,1)*F(0,1) + F(1,0)*F(1,0) - 1.0)
               + 2.0 * F(0,0)*F(0,1)*F(1,0);

      return  P;
   }

   static_tensor<scalar_type, 2>
   AK3(const static_matrix<scalar_type, 2, 2>& F) const
   {
      static_tensor<scalar_type, 2> A = static_tensor<scalar_type, 2>::Zero();

      //block A11
      A(0,0) = 3.0 * F(0,0)*F(0,0) + F(0,1)*F(0,1) + F(1,0)*F(1,0);
      A(0,1) = 2.0 * F(0,0)*F(0,1) + F(1,0)*F(1,1);
      A(1,0) = 2.0 * F(0,0)*F(1,0) + F(0,1)*F(1,1);
      A(1,1) = F(0,1)*F(1,0);

      //block A12
      A(0,2) = 2.0 * F(0,1)*F(0,0) + F(1,0)*F(1,1);
      A(0,3) = 3.0 * F(0,1)*F(0,1) + F(0,0)*F(0,0) + F(1,1)*F(1,1);
      A(1,2) = F(0,0)*F(1,1);
      A(1,3) = 2.0 * F(0,1)*F(1,1) + F(0,0)*F(1,0);

      //block A21
      A(2,0) = 2.0 * F(1,0)*F(0,0) + F(0,1)*F(1,1);
      A(2,1) = F(0,0)*F(1,1);
      A(3,0) = 3.0 * F(1,0)*F(1,0) + F(0,0)*F(0,0) + F(1,1)*F(1,1);
      A(3,1) = 2.0 * F(1,0)*F(1,1) + F(0,0)*F(0,1);

      //block A22
      A(2,2) = F(0,1)*F(1,0);
      A(2,3) = 2.0 * F(1,1)*F(0,1) + F(0,0)*F(1,0);
      A(3,2) = 2.0 * F(1,1)*F(1,0) + F(0,0)*F(0,1);
      A(3,3) = 3.0 * F(1,1)*F(1,1) + F(0,1)*F(0,1) + F(1,0)*F(1,0);
      return A;
   }


   static_tensor<scalar_type, 3>
   AK3(const static_matrix<scalar_type, 3, 3>& F) const
   {
      static_tensor<scalar_type, 3> A = static_tensor<scalar_type, 3>::Zero();

      //block A11
      A(0,0) = 12.0 * F(0,0)*F(0,0) + 4.0 * (F(0,1)*F(0,1) + F(1,0)*F(1,0) - 1.0);
      A(0,1) = 8.0 * F(0,0)*F(0,1) + 2.0 * F(1,0)*F(1,1);
      A(1,0) = 8.0 * F(0,0)*F(1,0) + 2.0 * F(0,1)*F(1,1);
      A(1,1) = 2.0 * F(0,1)*F(1,0);

      //block A12
      A(0,2) = 8.0 * F(0,1)*F(0,0) + 2.0 * F(1,0)*F(1,1);
      A(0,3) = 12.0 * F(0,1)*F(0,1) + 4.0 * (F(0,0)*F(0,0) + F(1,1)*F(1,1) - 1.0);
      A(1,2) = 2.0 * F(0,0)*F(1,1);
      A(1,3) = 8.0 * F(0,1)*F(1,1) + 2.0 * F(0,0)*F(1,0);

      //block A21
      A(2,0) = 8.0 * F(1,0)*F(0,0) + 2.0 * F(0,1)*F(1,1);
      A(2,1) = 2.0 * F(0,0)*F(1,1);
      A(3,0) = 12.0 * F(1,0)*F(1,0) + 4.0 * (F(0,0)*F(0,0) + F(1,1)*F(1,1) - 1.0);
      A(3,1) = 8.0 * F(1,0)*F(1,1) + 2.0 * F(0,0)*F(0,1);

      //block A22
      A(2,2) = 2.0 * F(0,1)*F(1,0);
      A(2,3) = 8.0 * F(1,1)*F(0,1) + 2.0 * F(0,0)*F(1,0);
      A(3,2) = 8.0 * F(1,1)*F(1,0) + 2.0 * F(0,0)*F(0,1);
      A(3,3) = 12.0 * F(1,1)*F(1,1) + 4.0 * (F(0,1)*F(0,1) + F(1,0)*F(1,0) - 1.0);

      return A;
   }


public:
   MoonleyLaw()
   : m_mu1(0.0), m_mu2(0.0), m_lambda(0.0)
   {}

   MoonleyLaw(const scalar_type mu1, const scalar_type mu2, const scalar_type lambda)
   : m_mu1(mu1), m_mu2(mu2), m_lambda(lambda)
   {}

   void
   setMu1(const scalar_type mu1)
   {
      m_mu1 = mu1;
   }

   void
   setMu2(const scalar_type mu2)
   {
      m_mu2 = mu2;
   }

   void
   setLambda(const scalar_type lambda)
   {
      m_lambda = lambda;
   }

   scalar_type
   giveMu1() const {return m_mu1;}

   scalar_type
   giveMu2() const {return m_mu2;}

   scalar_type
   giveLambda() const {return m_lambda;}


   template< int DIM>
   scalar_type
   compute_energy(const static_matrix<scalar_type, DIM, DIM>& F) const
   {
      const scalar_type J = F.determinant();
      if(J <=0.0){
         const std::string mess = "J= " + std::to_string(J) + " <= 0";
         throw std::invalid_argument(mess);
      }

      const static_matrix<scalar_type, DIM, DIM> C = compute_CauchyGreenRightTensor(F);
      const static_matrix<scalar_type, DIM, DIM> EGl = compute_GreenLagrangeTensor(C);

      const scalar_type W1 = m_lambda / J;
      const scalar_type W2 = m_mu1 /2.0 * C.trace();
      const scalar_type W3 = m_mu2 /2.0 * (EGl.transpose()*EGl).trace();
      return  W1 + W2 + W3;
   }

   template< int DIM>
   static_matrix<scalar_type, DIM, DIM>
   compute_PK1(const static_matrix<scalar_type, DIM, DIM>& F) const
   {
      const scalar_type J = F.determinant();
      if(J <=0.0){
         const std::string mess = "J= " + std::to_string(J) + " <= 0";
         throw std::invalid_argument(mess);
      }

      const static_matrix<scalar_type, DIM, DIM> invF = F.inverse();
      const static_matrix<scalar_type, DIM, DIM> invFT = invF.transpose();

      const auto P1 = -m_lambda / J *invFT;
      const auto P2 = m_mu1 * F;
      const auto P3 = m_mu2/2.0*(F*F.transpose()*F - F);

      return  P1 + P2 + P3;
   }


   template<int DIM>
   static_tensor<scalar_type, DIM>
   compute_tangent_moduli_A(const static_matrix<scalar_type, DIM, DIM>& F) const
   {
      const scalar_type J = F.determinant();
      if(J <=0.0){
         const std::string mess = "J= " + std::to_string(J) + " <= 0";
         throw std::invalid_argument(mess);
      }

      const static_matrix<scalar_type, DIM, DIM> invF = F.inverse();
      const static_matrix<scalar_type, DIM, DIM> invFt = invF.transpose();

      const static_tensor<scalar_type, DIM> I4 = compute_IdentityTensor<scalar_type,DIM>();
      const static_tensor<scalar_type, DIM> invFt_invF = computeProductInf(invFt, invF);
      const static_tensor<scalar_type, DIM> invFt_invFt = computeKroneckerProduct(invFt, invFt);

      const auto A1 = m_lambda/J * (invFt_invFt + invFt_invF);
      const auto A2 = m_mu1 * I4;
      const auto A3 = m_mu2 / 2.0 * (AK3(F) - I4);

      return A1 + A2 + A3;
   }


   template<int DIM>
   std::pair<static_matrix<scalar_type, DIM, DIM>, static_tensor<scalar_type, DIM> >
   compute_whole_PK1(const static_matrix<scalar_type, DIM, DIM>& F) const
   {
      const scalar_type J = F.determinant();
      if(J <=0.0){
         const std::string mess = "J= " + std::to_string(J) + " <= 0";
         throw std::invalid_argument(mess);
      }

      const static_matrix<scalar_type, DIM, DIM> invF = F.inverse();
      const static_matrix<scalar_type, DIM, DIM> invFt = invF.transpose();

      const static_tensor<scalar_type, DIM> I4 = compute_IdentityTensor<scalar_type,DIM>();
      const static_tensor<scalar_type, DIM> invFt_invF = computeProductInf(invFt, invF);
      const static_tensor<scalar_type, DIM> invFt_invFt = computeKroneckerProduct(invFt, invFt);

      const auto P1 = -m_lambda / J * invFt;
      const auto P2 = m_mu1 * F;
      const auto P3 = m_mu2/2.0*(F*F.transpose()*F - F);

      const auto A1 = m_lambda / J * (invFt_invFt + invFt_invF);
      const auto A2 = m_mu1 * I4;
      const auto A3 = m_mu2 / 2.0 * (AK3(F) - I4);

      return std::make_pair(P1 + P2 + P3, A1 + A2 + A3);
   }

};
