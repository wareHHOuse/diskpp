/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
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
#include "core/mechanics/behaviors/laws/materialData.hpp"
#include "core/mechanics/behaviors/maths_tensor.hpp"
#include "core/mechanics/behaviors/maths_utils.hpp"
#include "core/mechanics/behaviors/tensor_conversion.hpp"
#include "core/mechanics/deformation_tensors.hpp"
#include "mesh/point.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

namespace disk
{

// Law for Linear Isotropic and Kinematic Hardening model with von Mises Criteria  in small

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

template<typename T, int DIM>
class Cavitation_qp
{
  public:
    typedef T                                    scalar_type;
    typedef static_matrix<scalar_type, DIM, DIM> static_matrix_type;
    typedef static_matrix<scalar_type, 3, 3>     static_matrix_type3D;
    const static size_t dimension = DIM;
    typedef MaterialData<scalar_type>            data_type;

  private:
    // coordinat and weight of considered gauss point.
    point<scalar_type, DIM> m_point;
    scalar_type             m_weight;

    // internal variables at previous step
    static_matrix_type m_F_prev; // elastic strain

    // internal variables at current step
    static_matrix_type m_F_curr; // elastic strain

    static_tensor<scalar_type, DIM>
    elastic_modulus(const data_type& data) const
    {

        return 2 * data.getMu() * compute_IdentitySymTensor<scalar_type, DIM>() +
               data.getLambda() * compute_IxI<scalar_type, DIM>();
    }

    scalar_type
    compute_U(const data_type& data, scalar_type J) const
    {
        switch (data.getType())
        {
            case 1: return log(J);
            case 2: return (J - 1.0);
            case 3: return log10(J);
            case 4: return 1.0 / (1.0 - J);
            case 5: return (J * J - 1.0);
            case 6: return sqrt((J * J - 1.0 - 2.0 * log(J)) / 2.0);

            default: throw std::invalid_argument("NeoHookeanLaw: m_type have to be <= 6");
        }
    }

    scalar_type
    compute_T1(const data_type& data, scalar_type J) const
    {
        switch (data.getType())
        {
            case 1: return log(J);
            case 2: return J * (J - 1.0);
            case 3: return log(J) / (log(10) * log(10));
            case 4: return (J - 1.0) / (J * J);
            case 5: return 2 * J * J * (J * J - 1.0);
            case 6: return (J * J - 1.0) / 2.0;

            default: throw std::invalid_argument("NeoHookeanLaw: m_type have to be <= 6");
        }
    }

    scalar_type
    compute_T2(const data_type& data, scalar_type J) const
    {
        switch (data.getType())
        {
            case 1: return 1.0;
            case 2: return J * (2.0 * J - 1.0);
            case 3: return 1.0 / (log(10) * log(10));
            case 4: return (2.0 - J) / (J * J);
            case 5: return J * J * (8.0 * J * J - 4.0);
            case 6: return J * J;

            default: throw std::invalid_argument("NeoHookeanLaw: m_type have to be <= 6");
        }
    }

    static_tensor<scalar_type, DIM>
    compute_tangent_moduli_A(const data_type& data) const
    {
        const scalar_type J = m_F_curr.determinant();
        if (J <= 0.0)
        {
            const std::string mess = "J= " + std::to_string(J) + " <= 0";
            throw std::invalid_argument(mess);
        }

        const static_matrix_type invF  = m_F_curr.inverse();
        const static_matrix_type invFt = invF.transpose();
        const static_matrix_type C     = convertFtoCauchyGreenRight(m_F_curr);

        const scalar_type trace_C = C.trace();
        const scalar_type T1      = compute_T1(J);
        const scalar_type T2      = compute_T2(J);

        const static_tensor<scalar_type, DIM> I4          = compute_IdentityTensor<scalar_type, DIM>();
        const static_tensor<scalar_type, DIM> invFt_invF  = computeProductInf(invFt, invF);
        const static_tensor<scalar_type, DIM> invFt_invFt = computeKroneckerProduct(invFt, invFt);
        const static_tensor<scalar_type, DIM> F_F         = computeKroneckerProduct(m_F_curr, m_F_curr);

        const auto Aiso = data.getMu() * std::pow(3.0, -0.25) *
                          (std::pow(trace_C, -0.25) * I4 - 0.5 * std::pow(trace_C, -5.0 / 4.0) * F_F);
        const auto Avol = data.getLambda() * T2 * invFt_invFt + (data.getMu() - data.getLambda() * T1) * invFt_invF;

        return Aiso + Avol;
    }

  public:
    Cavitation_qp() : m_weight(0), m_F_prev(static_matrix_type::Zero()), m_F_curr(static_matrix_type::Zero()) {}

    Cavitation_qp(const point<scalar_type, DIM>& point, const scalar_type& weight) :
      m_point(point), m_weight(weight), m_F_prev(static_matrix_type::Zero()), m_F_curr(static_matrix_type::Zero())
    {
    }

    point<scalar_type, DIM>
    point() const
    {
        return m_point;
    }

    scalar_type
    weight() const
    {
        return m_weight;
    }

    bool
    is_plastic() const
    {
        return false;
    }

    static_matrix_type3D
    getElasticStrain() const
    {
        return convertMatrix3DwithOne(m_F_curr);
    }

    static_matrix_type3D
    getPlasticStrain() const
    {
        return static_matrix_type3D::Zero();
    }

    static_matrix_type
    getTotalStrain() const
    {
        return m_F_curr;
    }

    static_matrix_type
    getTotalStrainPrev() const
    {
        return m_F_prev;
    }

    scalar_type
    getAccumulatedPlasticStrain() const
    {
        return scalar_type(0);
    }

    void
    update()
    {
        m_F_prev = m_F_curr;
    }

    static_matrix_type
    compute_stress(const data_type& data) const
    {
        const scalar_type J = m_F_curr.determinant();
        if (J <= 0.0)
        {
            const std::string mess = "J= " + std::to_string(J) + " <= 0";
            throw std::invalid_argument(mess);
        }

        const static_matrix_type invF = m_F_curr.inverse();
        const scalar_type        T1   = compute_T1(J);
        const static_matrix_type C    = m_F_curr.transpose() * m_F_curr;

        const auto Piso = data.getMu() * std::pow(3.0 * C.trace(), -1.0 / 4.0) * m_F_curr;
        const auto Pvol = (data.getLambda() * T1 - data.getMu()) * invF.transpose();

        return Piso + Pvol;
    }

    static_matrix_type
    compute_stressPrev(const data_type& data) const
    {
        const scalar_type J = m_F_prev.determinant();
        if (J <= 0.0)
        {
            const std::string mess = "J= " + std::to_string(J) + " <= 0";
            throw std::invalid_argument(mess);
        }

        const static_matrix_type invF = m_F_prev.inverse();
        const scalar_type        T1   = compute_T1(J);
        const static_matrix_type C    = m_F_prev.transpose() * m_F_prev;

        const auto Piso = data.getMu() * std::pow(3.0 * C.trace(), -1.0 / 4.0) * m_F_prev;
        const auto Pvol = (data.getLambda() * T1 - data.getMu()) * invF.transpose();

        return Piso + Pvol;
    }

    scalar_type
    compute_energy(const data_type& data) const
    {
        const scalar_type J = m_F_curr.determinant();
        if (J <= 0.0)
        {
            const std::string mess = "J= " + std::to_string(J) + " <= 0";
            throw std::invalid_argument(mess);
        }

        const static_matrix_type C = m_F_curr.transpose() * m_F_curr;

        const scalar_type Wiso = 2.0 * data.getMu() / std::pow(3.0, 5.0 / 4.0) * std::pow(C.trace(), 3.0 / 4.0);
        const scalar_type Wvol = data.getLambda() / 2.0 * compute_U(J) * compute_U(J) - data.getMu() * log(J);

        return Wiso + Wvol;
    }

    std::pair<static_matrix_type, static_tensor<scalar_type, DIM>>
    compute_whole(const static_matrix_type& F_curr, const data_type& data, bool tangentmodulus = true)
    {
        // is always elastic
        m_F_curr = F_curr;

        const auto PK1 = this->compute_stress(data);
        const auto A   = this->compute_tangent_moduli_A(data);

        return std::make_pair(PK1, A);
    }
};
}
