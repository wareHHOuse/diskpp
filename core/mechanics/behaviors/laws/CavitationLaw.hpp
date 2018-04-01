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

#include <iostream>

#include <math.h>
#include <sstream>
#include <string>

#include "mechanics/behaviors/maths_tensor.hpp"
#include "common/eigen.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

namespace disk {

/* Material: Neo-nookean for cavitation
 * Energy :  W(F) = Wiso(F) + Wvol(F)
 *   - Wiso(F) =    2*mu/3^{5/4}  *[tr(F^T * F)]^{3/4}
 *   - Wvol(F) =  lambda/2 * U(J)**2 - mu * ln(J)
 * ** We set T1(J) = J * U(J) * U'(J) and T2(J) =  U(J) * J *( U''(J) * J + U'(J)) + ( J * U'(J))^2
 * Stress :  PK1(F) = Piso(F) + Pvol(F)
 *   - Piso(F) =  mu * (3*[tr(F^T * F)])^{-1/4}
 *   - Pvol(F) = (lambda * T1(J) -mu) * F^{-T}
 * Module :  A(F) = PK1(F) = Aiso(F) + Avol(F)
 *   - Aiso(F) = mu * 3^{-1/4} * ( [tr(F^T * F)])^{-1/4} * I4 -1/4 * [tr(F^T * F)])^{-5/4} F
 * \kronecker F)
 *   - Avol(F) = (mu-lambda * T1(J)) * F^{-T} \time_inf F^{-1}
 *                  + lambda * T2(J) * F^{-T} \kronecker F^{-T}
 */

/* Laws:
 * 1- U(J) = ln(J)
 * 2- U(J) = J -1
 * 3- U(J) = log10(J)
 * 4- U(J) = 1 -1 /J
 * 5- U(J) = J^2 -1
 * 6- U(J) = sqrt( ( J^2 -1 - 2 *ln(J)) /2)
 * */

template<typename scalar_type>
class CavitationLaw
{
   scalar_type m_mu;
   scalar_type m_lambda;
   size_t      m_type;

   const size_t maxtype = 6;

   scalar_type
   compute_U(scalar_type J) const
   {
      if (m_type == 1)
         return log(J);
      else if (m_type == 2)
         return (J - 1.0);
      else if (m_type == 3)
         return log10(J);
      else if (m_type == 4)
         return 1.0 / (1.0 - J);
      else if (m_type == 5) {
         scalar_type J2 = J * J;
         return (J2 - 1.0);
      } else if (m_type == 6)
         return sqrt((J * J - 1.0 - 2.0 * log(J)) / 2.0);
      else
         throw std::invalid_argument("CavitationLaw: m_type have to be <= 6");
   }

   scalar_type
   compute_T1(scalar_type J) const
   {
      if (m_type == 1)
         return log(J);
      else if (m_type == 2)
         return J * (J - 1.0);
      else if (m_type == 3)
         return log(J) / (log(10) * log(10));
      else if (m_type == 4)
         return (J - 1.0) / (J * J);
      else if (m_type == 5) {
         scalar_type J2 = J * J;
         return 2 * J2 * (J2 - 1.0);
      } else if (m_type == 6)
         return (J * J - 1.0) / 2.0;
      else
         throw std::invalid_argument("CavitationLaw: m_type have to be <= 6");
   }

   scalar_type
   compute_T2(scalar_type J) const
   {
      if (m_type == 1)
         return 1.0;
      else if (m_type == 2)
         return J * (2.0 * J - 1.0);
      else if (m_type == 3)
         return 1.0 / (log(10) * log(10));
      else if (m_type == 4)
         return (2.0 - J) / (J * J);
      else if (m_type == 5) {
         scalar_type J2 = J * J;
         return J2 * (8.0 * J2 - 4.0);
      } else if (m_type == 6)
         return J * J;
      else
         throw std::invalid_argument("CavitationLaw: m_type have to be <= 6");
   }

 public:
   CavitationLaw() : m_mu(0.0), m_lambda(0.0), m_type(1)
   {
      if (m_type <= 0 || m_type > maxtype) {
         std::cout << "Unknown option for NeoNookean material" << '\n';
         std::cout << "We use U(J) = ln(J)" << '\n';
         m_type = 1;
      }
   }

   CavitationLaw(const scalar_type mu, const scalar_type lambda, const size_t type) :
     m_mu(mu), m_lambda(lambda), m_type(type)
   {
      if (m_type <= 0 || m_type > maxtype) {
         std::cout << "Unknown option for NeoNookean material" << '\n';
         std::cout << "We use U(J) = ln(J)" << '\n';
         m_type = 1;
      }
   }

   void
   setMu(const scalar_type mu)
   {
      m_mu = mu;
   }

   void
   setLambda(const scalar_type lambda)
   {
      m_lambda = lambda;
   }

   void
   setType(const size_t type)
   {
      m_type = type;
      if (m_type <= 0 || m_type > maxtype) {
         std::cout << "Unknown option for NeoNookean material" << '\n';
         std::cout << "We use U(J) = ln(J)" << '\n';
         m_type = 1;
      }
   }

   scalar_type
   giveMu() const
   {
      return m_mu;
   }

   scalar_type
   giveLambda() const
   {
      return m_lambda;
   }

   template<int DIM>
   scalar_type
   compute_energy(const static_matrix<scalar_type, DIM, DIM>& F) const
   {
      const scalar_type J = F.determinant();
      if (J <= 0.0) {
         const std::string mess = "J= " + std::to_string(J) + " <= 0";
         throw std::invalid_argument(mess);
      }

      const static_matrix<scalar_type, DIM, DIM> C = F.transpose() * F;

      const scalar_type Wiso =
        2.0 * m_mu / std::pow(3.0, 5.0 / 4.0) * std::pow(C.trace(), 3.0 / 4.0);
      const scalar_type Wvol = m_lambda / 2.0 * compute_U(J) * compute_U(J) - m_mu * log(J);

      return Wiso + Wvol;
   }

   template<int DIM>
   static_matrix<scalar_type, DIM, DIM>
   compute_PK1(const static_matrix<scalar_type, DIM, DIM>& F) const
   {
      const scalar_type J = F.determinant();
      if (J <= 0.0) {
         const std::string mess = "J= " + std::to_string(J) + " <= 0";
         throw std::invalid_argument(mess);
      }

      const static_matrix<scalar_type, DIM, DIM> invF = F.inverse();
      const scalar_type                          T1   = compute_T1(J);
      const static_matrix<scalar_type, DIM, DIM> C    = F.transpose() * F;

      const auto Piso = m_mu * std::pow(3.0 * C.trace(), -1.0 / 4.0) * F;
      const auto Pvol = (m_lambda * T1 - m_mu) * invF.transpose();

      return Piso + Pvol;
   }

   template<int DIM>
   static_tensor<scalar_type, DIM>
   compute_tangent_moduli_A(const static_matrix<scalar_type, DIM, DIM>& F) const
   {
      const scalar_type J = F.determinant();
      if (J <= 0.0) {
         const std::string mess = "J= " + std::to_string(J) + " <= 0";
         throw std::invalid_argument(mess);
      }

      const static_matrix<scalar_type, DIM, DIM> invF  = F.inverse();
      const static_matrix<scalar_type, DIM, DIM> invFt = invF.transpose();
      const static_matrix<scalar_type, DIM, DIM> C     = F.transpose() * F;

      const scalar_type trace_C = C.trace();
      const scalar_type T1      = compute_T1(J);
      const scalar_type T2      = compute_T2(J);

      const static_tensor<scalar_type, DIM> I4         = compute_IdentityTensor<scalar_type, DIM>();
      const static_tensor<scalar_type, DIM> invFt_invF = computeProductInf(invFt, invF);
      const static_tensor<scalar_type, DIM> invFt_invFt = computeKroneckerProduct(invFt, invFt);
      const static_tensor<scalar_type, DIM> F_F         = computeKroneckerProduct(F, F);

      const auto Aiso = m_mu * std::pow(3.0, -0.25) *
                        (std::pow(trace_C, -0.25) * I4 - 0.5 * std::pow(trace_C, -5.0 / 4.0) * F_F);
      const auto Avol = m_lambda * T2 * invFt_invFt + (m_mu - m_lambda * T1) * invFt_invF;

      return Aiso + Avol;
   }

   template<int DIM>
   std::pair<static_matrix<scalar_type, DIM, DIM>, static_tensor<scalar_type, DIM>>
   compute_whole_PK1(const static_matrix<scalar_type, DIM, DIM>& F) const
   {
      const scalar_type J = F.determinant();
      if (J <= 0.0) {
         const std::string mess = "J= " + std::to_string(J) + " <= 0";
         throw std::invalid_argument(mess);
      }

      const static_matrix<scalar_type, DIM, DIM> invF  = F.inverse();
      const static_matrix<scalar_type, DIM, DIM> invFt = invF.transpose();
      const static_matrix<scalar_type, DIM, DIM> C     = F.transpose() * F;

      const scalar_type trace_C = C.trace();
      const scalar_type T1      = compute_T1(J);
      const scalar_type T2      = compute_T2(J);

      const static_tensor<scalar_type, DIM> I4         = compute_IdentityTensor<scalar_type, DIM>();
      const static_tensor<scalar_type, DIM> invFt_invF = computeProductInf(invFt, invF);
      const static_tensor<scalar_type, DIM> invFt_invFt = computeKroneckerProduct(invFt, invFt);
      const static_tensor<scalar_type, DIM> F_F         = computeKroneckerProduct(F, F);

      const auto Piso = m_mu * std::pow(3.0 * trace_C, -1.0 / 4.0) * F;
      const auto Pvol = (m_lambda * T1 - m_mu) * invF.transpose();

      const auto Aiso = m_mu * std::pow(3.0, -0.25) *
                        (std::pow(trace_C, -0.25) * I4 - 0.5 * std::pow(trace_C, -5.0 / 4.0) * F_F);
      const auto Avol = m_lambda * T2 * invFt_invFt + (m_mu - m_lambda * T1) * invFt_invF;

      return std::make_pair(Piso + Pvol, Aiso + Avol);
   }
};
}