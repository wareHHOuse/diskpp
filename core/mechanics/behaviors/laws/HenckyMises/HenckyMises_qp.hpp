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
class HenckyMises_Data
{
  private:
    scalar_type m_lambda;
    scalar_type m_mu;

  public:
    HenckyMises_Data() : m_lambda(1.0), m_mu(1.0) {}

    HenckyMises_Data(const scalar_type& lambda, const scalar_type& mu) : m_lambda(lambda), m_mu(mu) {}

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
    getLambda() const
    {
        return m_lambda;
    }

    scalar_type
    getMu() const
    {
        return m_mu;
    }

    void
    print() const
    {
        std::cout << "Material parameters: " << std::endl;
        std::cout << "* E: " << getE() << std::endl;
        std::cout << "* Nu: " << getNu() << std::endl;
        std::cout << "* Lambda: " << getLambda() << std::endl;
        std::cout << "* Mu: " << getMu() << std::endl;
    }
};

// Input : symetric stain tensor(Gs)

//           dev = normL2(Gs - trace(Gs) / dim * Id)

//             Stress : sigma = 2 *\tilde{mu}(dev(Gs)) * Gs + \tilde{lambda}(dev(Gs)) * trace(Gs) * Id
// \tilde{mu}(dev(Gs)) = mu * (1 + (1 + dev(Gs)) ^ {-1 / 2})
// \tilde{lambda}(dev(Gs)) = ((lambda + mu / 2) - mu / 2 * (1 + dev(Gs)) ^ {-1 / 2})

//                          Tangent Moduli : C = 2 * mu * I4 + lambda * prod_Kronecker(Id, Id) /
//                                                               it is the elastic moduli

template<typename T, int DIM>
class HenckyMises_qp
{
  public:
    typedef T                                    scalar_type;
    typedef static_matrix<scalar_type, DIM, DIM> static_matrix_type;
    typedef static_matrix<scalar_type, 3, 3>     static_matrix_type3D;
    typedef HenckyMises_Data<scalar_type>        data_type;

    const static size_t dimension = DIM;

  private:

    // coordinat and weight of considered gauss point.
    point<scalar_type, DIM> m_point;
    scalar_type             m_weight;

    // internal variables at previous step
    static_matrix_type3D m_estrain_prev; // elastic strain

    // internal variables at current step
    static_matrix_type3D m_estrain_curr; // elastic strain

    scalar_type
    tildemu(const data_type& data, const scalar_type& normL2dev) const
    {
        return data.getMu() * (1.0 + 1.0 / std::sqrt(1 + normL2dev * normL2dev));
    }

    scalar_type
    tildelambda(const data_type& data, const scalar_type& normL2dev) const
    {
        return (data.getLambda() + data.getMu() / 2.0) - data.getMu() / (2.0 * std::sqrt(1.0 + normL2dev * normL2dev));
    }

    scalar_type
    derivativetildemu(const data_type& data, const scalar_type& normL2dev) const
    {
        const auto term3_2 = (1.0 + normL2dev * normL2dev) * std::sqrt(1.0 + normL2dev * normL2dev);

        return -data.getMu() * normL2dev / term3_2;
    }

    scalar_type
    derivativetildelambda(const data_type& data, const scalar_type& normL2dev) const
    {
        return derivativetildemu(normL2dev) / scalar_type(2);
    }

    static_tensor<scalar_type, DIM>
    elastic_modulus(const data_type& data) const
    {

        return 2 * data.getMu() * compute_IdentitySymTensor<scalar_type, DIM>() +
               data.getLambda() * compute_IxI<scalar_type, DIM>();
    }


  public:
    HenckyMises_qp() :
      m_weight(0), m_estrain_prev(static_matrix_type3D::Zero()), m_estrain_curr(static_matrix_type3D::Zero())
    {
    }

    HenckyMises_qp(const point<scalar_type, DIM>& point, const scalar_type& weight) :
      m_point(point), m_weight(weight), m_estrain_prev(static_matrix_type3D::Zero()), m_estrain_curr(static_matrix_type3D::Zero())
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
        return m_estrain_curr;
    }

    static_matrix_type3D
    getPlasticStrain() const
    {
        return static_matrix_type3D::Zero();
    }

    static_matrix_type
    getTotalStrain() const
    {
        return m_estrain_curr.block(0, 0, DIM, DIM);
    }

    static_matrix_type
    getTotalStrainPrev() const
    {
        return m_estrain_prev.block(0, 0, DIM, DIM);
    }

    scalar_type
    getAccumulatedPlasticStrain() const
    {
        return 0.0;
    }

    void
    update()
    {
        m_estrain_prev = m_estrain_curr;
    }

    static_matrix_type
    compute_stress(const data_type& data) const
    {
        const static_matrix_type Id = static_matrix_type::Identity();
        const static_matrix_type Gs = m_estrain_curr.block(0, 0, DIM, DIM);

        const scalar_type normFrodev = deviator(Gs).norm();

        return tildemu(data, normFrodev) * Gs + tildelambda(data, normFrodev) * Gs.trace() * Id;
    }

    std::pair<static_matrix_type, static_tensor<scalar_type, DIM>>
    compute_whole(const static_matrix_type& strain_curr, const data_type& data, bool tangentmodulus = true)
    {
        static_tensor<scalar_type, DIM> Cep = elastic_modulus(data);

        // is always elastic
        m_estrain_curr = convertMatrix3D(strain_curr);

        // compute Cauchy stress
        const static_matrix_type stress = this->compute_stress(data);

        return std::make_pair(stress, Cep);
    }
};
}
