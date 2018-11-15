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
#include "mesh/point.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

namespace disk
{

// Law for Isotropic Hardening model with von Mises Criteria in small deformation
// where the curve R(p) is given point by point
// see https://www.code-aster.org/doc/default/en/man_r/r5/r5.03.02.pdf section 3.1.2

template<typename T, int DIM>
class IsotropicHardeningVMis_qp : public law_qp_bones<T, DIM>
{
  public:
    const static size_t                          dimension = DIM;
    typedef T                                    scalar_type;
    typedef static_matrix<scalar_type, DIM, DIM> static_matrix_type;
    typedef static_matrix<scalar_type, 3, 3>     static_matrix_type3D;
    typedef MaterialData<scalar_type>            data_type;

  private:
    // internal variables at previous step
    static_matrix_type3D m_pstrain_prev; // plastic strain
    scalar_type          m_p_prev;       // cumulate plastic strain

    // internal variables at current step
    static_matrix_type3D m_pstrain_curr;    // plastic strain
    scalar_type          m_p_curr;          // cumulate plastic strain
    bool                 m_is_plastic_curr; // the gauss point is plastic ?

  public:
    IsotropicHardeningVMis_qp() :
      law_qp_bones<T, DIM>(), m_pstrain_prev(static_matrix_type3D::Zero()), m_p_prev(scalar_type(0)),
      m_pstrain_curr(static_matrix_type3D::Zero()), m_p_curr(scalar_type(0)), m_is_plastic_curr(false)
    {
    }

    IsotropicHardeningVMis_qp(const point<scalar_type, DIM>& point, const scalar_type& weight) :
      law_qp_bones<T, DIM>(point, weight), m_pstrain_prev(static_matrix_type3D::Zero()), m_p_prev(scalar_type(0)),
      m_pstrain_curr(static_matrix_type3D::Zero()), m_p_curr(scalar_type(0)), m_is_plastic_curr(false)
    {
    }

    bool
    is_plastic() const
    {
        return m_is_plastic_curr;
    }

    static_matrix_type3D
    getPlasticStrain() const
    {
        return m_pstrain_curr;
    }

    static_matrix_type
    getTotalStrain() const
    {
        return convertMatrix<scalar_type, DIM>(this->m_estrain_curr + m_pstrain_curr);
    }

    static_matrix_type
    getTotalStrainPrev() const
    {
        return convertMatrix<scalar_type, DIM>(this->m_estrain_prev + m_pstrain_prev);
    }

    scalar_type
    getAccumulatedPlasticStrain() const
    {
        return m_p_curr;
    }

    void
    update()
    {
        law_qp_bones<T, DIM>::update();
        m_pstrain_prev = m_pstrain_curr;
        m_p_prev       = m_p_curr;
    }

    static_matrix_type3D
    compute_stress3D(const data_type& data) const
    {
        const static_matrix_type3D Id = static_matrix_type3D::Identity();

        const auto stress =
          2 * data.getMu() * this->m_estrain_curr + data.getLambda() * this->m_estrain_curr.trace() * Id;

        return stress;
    }

    static_matrix_type
    compute_stress(const data_type& data) const
    {
        return convertMatrix<scalar_type, DIM>(compute_stress3D(data));
    }

    static_matrix_type3D
    compute_stress3DPrev(const data_type& data) const
    {
        const static_matrix_type3D Id = static_matrix_type3D::Identity();

        const auto stress =
          2 * data.getMu() * this->m_estrain_prev + data.getLambda() * this->m_estrain_prev.trace() * Id;

        return stress;
    }

    static_matrix_type
    compute_stressPrev(const data_type& data) const
    {
        const static_matrix_type3D Id = static_matrix_type3D::Identity();

        const auto stress =
          2 * data.getMu() * this->m_estrain_prev + data.getLambda() * this->m_estrain_prev.trace() * Id;

        return convertMatrix<scalar_type, DIM>(stress);
    }

    std::pair<static_matrix_type3D, static_tensor<scalar_type, 3>>
    compute_whole3D(const static_matrix_type3D& strain_curr, const data_type& data, bool tangentmodulus = true)
    {
        static_tensor<scalar_type, 3> Cep      = this->elastic_modulus3D(data);
        const auto                    RpCurve  = data.getRpCurve();
        const auto                    nb_point = RpCurve.size();

        const static_matrix_type3D incr_strain   = strain_curr - convertMatrix3D(this->getTotalStrainPrev());
        const static_matrix_type3D estrain_trial = this->m_estrain_prev + incr_strain; // elastic strain trial

        if (nb_point == 0)
        {
            // We don't have points for the traction curve
            // we suppose that we are in linear elasticity
            // elastic evolution
            m_is_plastic_curr = false;
            // update
            this->m_estrain_curr = estrain_trial;
            m_pstrain_curr       = m_pstrain_prev;
            m_p_curr             = m_p_prev;
        }
        else
        {
            // we search i0 such that p_prev \in [p_i0, p_{i0+1}]
            size_t i0 = nb_point - 2;
            for (size_t i = 0; i < nb_point - 1; i++)
            {
                if (m_p_prev < RpCurve[i + 1].getP())
                {
                    i0 = i;
                    break;
                }
            }

            const scalar_type Rp0 = RpCurve[i0].getRp();
            const scalar_type p0  = RpCurve[i0].getP();
            const scalar_type H0  = (RpCurve[i0 + 1].getRp() - Rp0) / (RpCurve[i0 + 1].getP() - p0);

            // prediction
            const static_matrix_type3D se         = 2 * data.getMu() * deviator(estrain_trial);
            const scalar_type          se_eq      = this->sigmaeq(se);
            const scalar_type          Phi_trial0 = se_eq - Rp0 - H0 * (m_p_prev - p0);

            // debug informations
            // std::cout << "point: " << i0 << " on " << nb_point << std::endl;
            // std::cout << "Rp0: " << Rp0 << ", p0: " << p0 << ", H0: " << H0 << std::endl;
            // std::cout << "se_eq: " << se_eq << ", Rp " << Rp0 + H0 * (m_p_prev - p0) << ", Phi: " << Phi_trial0 <<
            // std::endl;

            if (std::abs(se.trace()) > 1E-8)
            {
                const std::string mess = "Se_trace= " + std::to_string(se.trace()) + " <= 0";
                throw std::invalid_argument(mess);
            }

            // check
            if (Phi_trial0 < scalar_type(0))
            {
                // elastic evolution
                m_is_plastic_curr = false;
            }
            else
            {
                // plastic evolution
                m_is_plastic_curr = true;
            }

            // corection
            if (m_is_plastic_curr)
            {
                const scalar_type troismu = 3 * data.getMu();
                size_t            i1      = nb_point - 2;
                for (size_t i = i0 + 1; i < nb_point - 1; i++)
                {
                    const scalar_type eq = RpCurve[i].getRp() - troismu * (m_p_prev - RpCurve[i].getRp()) - se_eq;

                    if (eq > scalar_type(0))
                    {
                        i1 = i - 1;
                        break;
                    }
                }

                const scalar_type Rp1 = RpCurve[i1].getRp();
                const scalar_type p1  = RpCurve[i1].getP();
                const scalar_type H1  = (RpCurve[i1 + 1].getRp() - Rp1) / (RpCurve[i1 + 1].getP() - p1);

                const scalar_type Phi_trial1 = se_eq - Rp1 - H1 * (m_p_prev - p1);

                const scalar_type          dem     = troismu + H1;
                const static_matrix_type3D normal  = scalar_type(3.) * se / (scalar_type(2.) * se_eq);
                const scalar_type          delta_p = Phi_trial1 / dem;

                // update
                m_p_curr             = m_p_prev + delta_p;
                this->m_estrain_curr = estrain_trial - delta_p * normal;
                m_pstrain_curr       = m_pstrain_prev + delta_p * normal;

                if (tangentmodulus)
                {
                    // compute cep coherent
                    const static_tensor<scalar_type, 3> nxn  = Kronecker(normal, normal);
                    const static_tensor<scalar_type, 3> Is   = IdentitySymTensor4<scalar_type, 3>();
                    const static_tensor<scalar_type, 3> Pdev = Is - IxI<scalar_type, 3>() / scalar_type(3);
                    const scalar_type                   mu2  = data.getMu() * data.getMu();

                    Cep += 4 * mu2 * (delta_p / se_eq - 1.0 / dem) * nxn - 6.0 * mu2 * delta_p / se_eq * Pdev;
                }
            }
            else
            {
                // update
                this->m_estrain_curr = estrain_trial;
                m_pstrain_curr       = m_pstrain_prev;
                m_p_curr             = m_p_prev;
            }
        }

        if (std::abs(m_pstrain_curr.trace()) > 1E-8)
        {
            const std::string mess = "eps_p= " + std::to_string(m_pstrain_curr.trace()) + " <= 0";
            throw std::invalid_argument(mess);
        }

        // std::cout << "ep:" << std::endl;
        // std::cout << m_pstrain_curr << std::endl;
        // std::cout << "ee:" << std::endl;
        // std::cout << this->m_estrain_curr << std::endl;
        // std::cout << "e:" << std::endl;
        // std::cout << this->m_estrain_curr + m_pstrain_curr << std::endl;
        // std::cout << "p:" << std::endl;
        // std::cout << m_p_curr << std::endl;

        // compute Cauchy stress
        const static_matrix_type3D stress = this->compute_stress3D(data);

        return std::make_pair(stress, Cep);
    }

    std::pair<static_matrix_type, static_tensor<scalar_type, DIM>>
    compute_whole(const static_matrix_type& strain_curr, const data_type& data, bool tangentmodulus = true)
    {
        const static_matrix_type3D strain3D_curr = convertMatrix3D(strain_curr);
        const auto                 behaviors3D   = this->compute_whole3D(strain3D_curr, data, tangentmodulus);

        const static_matrix_type              stress = convertMatrix<scalar_type, DIM>(behaviors3D.first);
        const static_tensor<scalar_type, DIM> Cep    = convertTensor<scalar_type, DIM>(behaviors3D.second);

        return std::make_pair(stress, Cep);
    }
};
}
