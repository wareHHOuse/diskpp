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
#include "core/mechanics/behaviors/maths_tensor.hpp"
#include "core/mechanics/behaviors/maths_utils.hpp"
#include "core/mechanics/behaviors/tensor_conversion.hpp"
#include "mesh/point.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

namespace disk
{

// Law for Linear Isotropic and Kinematic Hardening model with von Mises Criteria  in small

template<typename scalar_type>
class LinearIsotropicAndKinematicHardening_Data
{
  private:
    scalar_type m_lambda;
    scalar_type m_mu;
    scalar_type m_H;
    scalar_type m_K;
    scalar_type m_sigma_y0;

  public:
    LinearIsotropicAndKinematicHardening_Data() :
      m_lambda(1.0), m_mu(1.0), m_H(0), m_K(0), m_sigma_y0(std::numeric_limits<scalar_type>::max())
    {
    }

    LinearIsotropicAndKinematicHardening_Data(const scalar_type& lambda,
                                              const scalar_type& mu,
                                              const scalar_type& H,
                                              const scalar_type& K,
                                              const scalar_type& sigma_y0) :
      m_lambda(lambda),
      m_mu(mu), m_H(H), m_K(K), m_sigma_y0(sigma_y0)
    {
    }

    LinearIsotropicAndKinematicHardening_Data(const scalar_type& lambda, const scalar_type& mu) :
      m_lambda(lambda), m_mu(mu), m_H(0), m_K(0), m_sigma_y0(std::numeric_limits<scalar_type>::max())
    {
    }

    scalar_type
    getE() const
    {
        return m_mu * (3 * m_lambda + 2 * m_mu) / (m_lambda + m_mu);
    }

    scalar_type
    getNu() const
    {
        return m_lambda / (2 * (m_lambda + m_mu));
    }

    scalar_type
    getET() const
    {
        const scalar_type E = getE();
        return E * (m_H + 1.5 * m_K) / (m_H + 1.5 * m_K + E);
    }

    scalar_type
    getLambda() const
    {
        return m_lambda;
    }

    scalar_type
    getMu() const
    {
        return m_mu;
    }

    scalar_type
    getH() const
    {
        return m_H;
    }

    scalar_type
    getK() const
    {
        return m_K;
    }

    scalar_type
    getSigma_y0() const
    {
        return m_sigma_y0;
    }

    void
    print() const
    {
        std::cout << "Material parameters: " << std::endl;
        std::cout << "* E: " << getE() << std::endl;
        std::cout << "* Nu: " << getNu() << std::endl;
        std::cout << "* ET: " << getET() << std::endl;
        std::cout << "* H: " << getH() << std::endl;
        std::cout << "* K: " << getK() << std::endl;
        std::cout << "* Sy0: " << getSigma_y0() << std::endl;
        std::cout << "* Lambda: " << getLambda() << std::endl;
        std::cout << "* Mu: " << getMu() << std::endl;
    }
};

template<typename T, int DIM>
class LinearIsotropicAndKinematicHardening_qp
{
  public:
    typedef T                                                      scalar_type;
    typedef static_matrix<scalar_type, DIM, DIM>                   static_matrix_type;
    typedef static_matrix<scalar_type, 3, 3>                       static_matrix_type3D;
    typedef LinearIsotropicAndKinematicHardening_Data<scalar_type> data_type;

    const static size_t dimension = DIM;

  private:
    // coordinat and weight of considered gauss point.
    point<scalar_type, DIM> m_point;
    scalar_type             m_weight;

    // internal variables at previous step
    static_matrix_type3D m_estrain_prev; // elastic strain
    static_matrix_type3D m_pstrain_prev; // plastic strain
    scalar_type          m_p_prev;       // cumulate plastic strain

    // internal variables at current step
    static_matrix_type3D m_estrain_curr;    // elastic strain
    static_matrix_type3D m_pstrain_curr;    // plastic strain
    scalar_type          m_p_curr;          // cumulate plastic strain
    bool                 m_is_plastic_curr; // the gauss point is plastic ?

    static_tensor<scalar_type, 3>
    elastic_modulus(const data_type& data) const
    {

        return 2 * data.getMu() * compute_IdentitySymTensor<scalar_type, 3>() +
               data.getLambda() * compute_IxI<scalar_type, 3>();
    }

    scalar_type
    sigmaeq(const static_matrix_type3D& dev) const
    {
        return sqrt(scalar_type(1.5) * dev.squaredNorm());
    }

  public:
    LinearIsotropicAndKinematicHardening_qp() :
      m_weight(0), m_estrain_prev(static_matrix_type3D::Zero()), m_pstrain_prev(static_matrix_type3D::Zero()),
      m_p_prev(scalar_type(0)), m_estrain_curr(static_matrix_type3D::Zero()),
      m_pstrain_curr(static_matrix_type3D::Zero()), m_p_curr(scalar_type(0)), m_is_plastic_curr(false)
    {
    }

    LinearIsotropicAndKinematicHardening_qp(const point<scalar_type, DIM>& point, const scalar_type& weight) :
      m_point(point), m_weight(weight), m_estrain_prev(static_matrix_type3D::Zero()),
      m_pstrain_prev(static_matrix_type3D::Zero()), m_p_prev(scalar_type(0)),
      m_estrain_curr(static_matrix_type3D::Zero()), m_pstrain_curr(static_matrix_type3D::Zero()),
      m_p_curr(scalar_type(0)), m_is_plastic_curr(false)
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
        return m_is_plastic_curr;
    }

    static_matrix_type3D
    getElasticStrain() const
    {
        return m_estrain_curr;
    }

    static_matrix_type3D
    getPlasticStrain() const
    {
        return m_pstrain_curr;
    }

    static_matrix_type
    getTotalStrain() const
    {
        return convertMatrix<scalar_type, DIM>(m_estrain_curr + m_pstrain_curr);
    }

    static_matrix_type
    getTotalStrainPrev() const
    {
        return convertMatrix<scalar_type, DIM>(m_estrain_prev + m_pstrain_prev);
    }

    scalar_type
    getAccumulatedPlasticStrain() const
    {
        return m_p_curr;
    }

    void
    update()
    {
        m_estrain_prev = m_estrain_curr;
        m_pstrain_prev = m_pstrain_curr;
        m_p_prev       = m_p_curr;
    }

    static_matrix_type3D
    compute_stress3D(const data_type& data) const
    {
        const static_matrix_type3D Id = static_matrix_type3D::Identity();

        const auto stress = 2 * data.getMu() * m_estrain_curr + data.getLambda() * m_estrain_curr.trace() * Id;

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

        const auto stress = 2 * data.getMu() * m_estrain_prev + data.getLambda() * m_estrain_prev.trace() * Id;

        return stress;
    }

    static_matrix_type
    compute_stressPrev(const data_type& data) const
    {
        const static_matrix_type3D Id = static_matrix_type3D::Identity();

        const auto stress = 2 * data.getMu() * m_estrain_prev + data.getLambda() * m_estrain_prev.trace() * Id;

        return convertMatrix<scalar_type, DIM>(stress);
    }

    std::pair<static_matrix_type3D, static_tensor<scalar_type, 3>>
    compute_whole3D(const static_matrix_type3D& strain_curr, const data_type& data, bool tangentmodulus = true)
    {
        static_tensor<scalar_type, 3> Cep = elastic_modulus(data);

        // prediction
        const static_matrix_type3D incr_strain   = strain_curr - convertMatrix3D(this->getTotalStrainPrev());
        const static_matrix_type3D estrain_trial = m_estrain_prev + incr_strain; // elastic strain trial
        const static_matrix_type3D X_prev        = data.getK() * m_pstrain_prev; // back-stress previous
        const static_matrix_type3D se            = 2 * data.getMu() * deviator(estrain_trial) - X_prev;
        const scalar_type          se_eq         = sigmaeq(se);
        const scalar_type          Phi_trial     = se_eq - data.getSigma_y0() - data.getH() * m_p_prev;

        assert(std::abs(X_prev.trace()) <= 1E-8);
        assert(std::abs(se.trace()) <= 1E-8);

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
            m_p_curr       = m_p_prev + delta_p;
            m_estrain_curr = estrain_trial - delta_p * normal;
            m_pstrain_curr = m_pstrain_prev + delta_p * normal;

            assert(std::abs(m_pstrain_curr.trace()) <= 1E-8);

            if (tangentmodulus)
            {
                // compute cep coherent
                const static_tensor<scalar_type, 3> nxn  = computeKroneckerProduct(normal, normal);
                const static_tensor<scalar_type, 3> IxI  = compute_IxI<scalar_type, 3>();
                const static_tensor<scalar_type, 3> Is   = compute_IdentityTensor<scalar_type, 3>();
                const static_tensor<scalar_type, 3> Pdev = Is - IxI / scalar_type(3);
                const scalar_type                   mu2  = data.getMu() * data.getMu();

                Cep += 4 * mu2 * (delta_p / se_eq - 1.0 / dem) * nxn - 6.0 * mu2 * delta_p / se_eq * Pdev;
            }
        }
        else
        {
            // update
            m_estrain_curr = estrain_trial;
            m_pstrain_curr = m_pstrain_prev;
            m_p_curr       = m_p_prev;
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
