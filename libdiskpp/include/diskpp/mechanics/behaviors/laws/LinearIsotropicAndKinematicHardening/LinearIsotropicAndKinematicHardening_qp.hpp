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

#include "diskpp/common/eigen.hpp"
#include "diskpp/mechanics/behaviors/laws/law_qp_bones.hpp"
#include "diskpp/mechanics/behaviors/laws/materialData.hpp"
#include "diskpp/mechanics/behaviors/maths_tensor.hpp"
#include "diskpp/mechanics/behaviors/maths_utils.hpp"
#include "diskpp/mechanics/behaviors/tensor_conversion.hpp"
#include "diskpp/mesh/point.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

namespace disk
{

// Law for Linear Isotropic and Kinematic Hardening model with von Mises Criteria  in small

template<typename T, int DIM>
class LinearIsotropicAndKinematicHardening_qp : public law_qp_bones<T, DIM>
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
    LinearIsotropicAndKinematicHardening_qp() :
      law_qp_bones<T, DIM>(), m_pstrain_prev(static_matrix_type3D::Zero()), m_p_prev(scalar_type(0)),
      m_pstrain_curr(static_matrix_type3D::Zero()), m_p_curr(scalar_type(0)), m_is_plastic_curr(false)
    {
    }

    LinearIsotropicAndKinematicHardening_qp(const point<scalar_type, DIM>& point, const scalar_type& weight) :
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
        return convertMatrix<scalar_type, DIM>(this->m_estrain_curr + this->m_pstrain_curr);
    }

    static_matrix_type
    getTotalStrainPrev() const
    {
        return convertMatrix<scalar_type, DIM>(this->m_estrain_prev + this->m_pstrain_prev);
    }

    scalar_type
    getEquivalentPlasticStrain() const
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
        static_tensor<scalar_type, 3> Cep = this->elastic_modulus3D(data);

        // prediction
        const static_matrix_type3D incr_strain   = strain_curr - convertMatrix3D(this->getTotalStrainPrev());
        const static_matrix_type3D estrain_trial = this->m_estrain_prev + incr_strain; // elastic strain trial
        const static_matrix_type3D X_prev        = data.getK() * m_pstrain_prev;       // back-stress previous
        const static_matrix_type3D se            = 2 * data.getMu() * deviator(estrain_trial) - X_prev;
        const scalar_type          se_eq         = this->sigmaeq(se);
        const scalar_type          Phi_trial     = se_eq - data.getSigma_y0() - data.getH() * m_p_prev;

        if ((std::abs(X_prev.trace()) / X_prev.norm()) > 1E-8)
        {
            const std::string mess = "X_trace= " + std::to_string(X_prev.trace()) + " <= 0";
            throw std::invalid_argument(mess);
        }

        if ((std::abs(se.trace()) / se.norm()) > 1E-8)
        {
            const std::string mess = "Se_trace= " + std::to_string(se.trace()) + " <= 0";
            throw std::invalid_argument(mess);
        }

        // check
        if (Phi_trial < scalar_type(0))
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
            const scalar_type dem = 3 * data.getMu() + data.getH() + scalar_type(3.) * data.getK() / scalar_type(2.);
            const static_matrix_type3D normal  = scalar_type(3.) * se / (scalar_type(2.) * se_eq);
            const scalar_type          delta_p = Phi_trial / dem;

            //  std::cout << "n:" << std::endl;
            //  std::cout << normal << std::endl;

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

        if ((std::abs(m_pstrain_curr.trace()) / m_pstrain_curr.norm()) > 1E-8)
        {
            const std::string mess = "eps_p= " + std::to_string(m_pstrain_curr.trace()) + " <= 0";
            throw std::invalid_argument(mess);
        }

        // std::cout << "ep:" << std::endl;
        // std::cout << m_pstrain_curr << std::endl;
        // std::cout << "ee:" << std::endl;
        // std::cout << m_estrain_curr << std::endl;
        // std::cout << "e:" << std::endl;
        // std::cout << m_estrain_curr + m_pstrain_curr << std::endl;

        // compute Cauchy stress
        const static_matrix_type3D stress = this->compute_stress3D(data);

        // std::cout << "stress:" << std::endl;
        // std::cout << stress << std::endl;

        // std::cout << 2 * data.getMu() * m_estrain_curr(2, 2) +
        //                data.getLambda() * (m_estrain_curr.trace())
        //           << std::endl;

        // std::cout << "p:" << std::endl;
        // std::cout << m_p_curr << std::endl;

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
