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
#include "core/mechanics/behaviors/laws/law_qp_bones.hpp"
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

/* Laws:
 * 1- U(J) = ln(J)
 * 2- U(J) = J -1
 * 3- U(J) = log10(J)
 * 4- U(J) = 1 -1 /J
 * 5- U(J) = J^2 -1
 * 6- U(J) = sqrt( ( J^2 -1 - 2 *ln(J)) /2)
 * */

template<typename T, int DIM>
class Neohookean_qp : public law_qp_bones<T, DIM>
{
  public:
    typedef T                                    scalar_type;
    typedef static_matrix<scalar_type, DIM, DIM> static_matrix_type;
    typedef static_matrix<scalar_type, 3, 3>     static_matrix_type3D;
    const static size_t                          dimension = DIM;
    typedef MaterialData<scalar_type>            data_type;

  private:
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

    static_tensor<scalar_type, 3>
    compute_tangent_moduli_A(const data_type& data) const
    {
        const scalar_type J = this->m_estrain_curr.determinant();
        if (J <= 0.0)
        {
            const std::string mess = "J= " + std::to_string(J) + " <= 0";
            throw std::invalid_argument(mess);
        }

        const static_matrix_type3D invF  = this->m_estrain_curr.inverse();
        const static_matrix_type3D invFt = invF.transpose();

        const scalar_type T1 = compute_T1(data, J);
        const scalar_type T2 = compute_T2(data, J);

        const static_tensor<scalar_type, 3> I4          = IdentityTensor4<scalar_type, 3>();
        const static_tensor<scalar_type, 3> invFt_invF  = ProductInf(invFt, invF);
        const static_tensor<scalar_type, 3> invFt_invFt = Kronecker(invFt, invFt);

        return data.getMu() * (I4 + invFt_invF) + data.getLambda() * (T2 * invFt_invFt - T1 * invFt_invF);
    }

  public:
    Neohookean_qp() : law_qp_bones<T, DIM>() {}

    Neohookean_qp(const point<scalar_type, DIM>& point, const scalar_type& weight) : law_qp_bones<T, DIM>(point, weight)
    {
    }

    static_matrix_type3D
    compute_stress3D(const data_type& data) const
    {
        const scalar_type J = this->m_estrain_curr.determinant();
        if (J <= 0.0)
        {
            const std::string mess = "J= " + std::to_string(J) + " <= 0";
            throw std::invalid_argument(mess);
        }

        const static_matrix_type3D invF = this->m_estrain_curr.inverse();
        const scalar_type          T1   = compute_T1(data, J);

        return data.getMu() * this->m_estrain_curr + (data.getLambda() * T1 - data.getMu()) * invF.transpose();
    }

    static_matrix_type
    compute_stress(const data_type& data) const
    {
        return convertMatrix<scalar_type, DIM>(compute_stress3D(data));
    }

    static_matrix_type3D
    compute_stressPrev3D(const data_type& data) const
    {
        const scalar_type J = this->m_estrain_prev.determinant();
        if (J <= 0.0)
        {
            const std::string mess = "J= " + std::to_string(J) + " <= 0";
            throw std::invalid_argument(mess);
        }

        const static_matrix_type3D invF = this->m_estrain_prev.inverse();
        const scalar_type          T1   = compute_T1(data, J);

        return data.getMu() * this->m_estrain_prev + (data.getLambda() * T1 - data.getMu()) * invF.transpose();
    }

    scalar_type
    compute_energy(const data_type& data) const
    {
        const scalar_type J = this->m_estrain_curr.determinant();
        if (J <= 0.0)
        {
            const std::string mess = "J= " + std::to_string(J) + " <= 0";
            throw std::invalid_argument(mess);
        }

        const static_matrix_type3D C = convertFtoCauchyGreenRight(this->m_estrain_curr);

        const scalar_type Wiso = data.getMu() / 2.0 * (C.trace() - 3);
        const scalar_type Wvol =
          data.getLambda() / 2.0 * compute_U(data, J) * compute_U(data, J) - data.getMu() * log(J);

        return Wiso + Wvol;
    }

    std::pair<static_matrix_type3D, static_tensor<scalar_type, 3>>
    compute_whole3D(const static_matrix_type3D& F_curr, const data_type& data, bool tangentmodulus = true)
    {
        // is always elastic
        this->m_estrain_curr = F_curr;

        const auto PK1 = this->compute_stress3D(data);
        const auto A   = this->compute_tangent_moduli_A(data);

        return std::make_pair(PK1, A);
    }

    std::pair<static_matrix_type, static_tensor<scalar_type, DIM>>
    compute_whole(const static_matrix_type& F_curr, const data_type& data, bool tangentmodulus = true)
    {
        const static_matrix_type3D F3D         = convertMatrix3DwithOne(F_curr);
        const auto                 behaviors3D = compute_whole3D(F3D, data, tangentmodulus);

        const static_matrix_type              stress = convertMatrix<scalar_type, DIM>(behaviors3D.first);
        const static_tensor<scalar_type, DIM> Cep    = convertTensor<scalar_type, DIM>(behaviors3D.second);

        return std::make_pair(stress, Cep);
    }
};
}
