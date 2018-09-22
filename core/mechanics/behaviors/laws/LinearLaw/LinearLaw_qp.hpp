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
#include "mesh/point.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

namespace disk
{

// Law for LinearLaw (test of finite deformations)

// Input : symetric stain tensor(Gs)

//           dev = normL2(Gs - trace(Gs) / dim * Id)

//             Stress : sigma = 2 *\tilde{mu}(dev(Gs)) * Gs + \tilde{lambda}(dev(Gs)) * trace(Gs) * Id
// \tilde{mu}(dev(Gs)) = mu * (1 + (1 + dev(Gs)) ^ {-1 / 2})
// \tilde{lambda}(dev(Gs)) = ((lambda + mu / 2) - mu / 2 * (1 + dev(Gs)) ^ {-1 / 2})

//                          Tangent Moduli : C = 2 * mu * I4 + lambda * prod_Kronecker(Id, Id) /
//                                                               it is the elastic moduli

template<typename T, int DIM>
class LinearLaw_qp
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
    static_matrix_type3D m_estrain_prev; // elastic strain

    // internal variables at current step
    static_matrix_type3D m_estrain_curr; // elastic strain

    static_tensor<scalar_type, 3>
    elastic_modulus(const data_type& data) const
    {

        return 2 * data.getMu() * compute_IdentitySymTensor<scalar_type, 3>() +
               data.getLambda() * compute_IxI<scalar_type, 3>();
    }

  public:
    LinearLaw_qp() :
      m_weight(0), m_estrain_prev(static_matrix_type3D::Zero()), m_estrain_curr(static_matrix_type3D::Zero())
    {
    }

    LinearLaw_qp(const point<scalar_type, DIM>& point, const scalar_type& weight) :
      m_point(point), m_weight(weight), m_estrain_prev(static_matrix_type3D::Zero()),
      m_estrain_curr(static_matrix_type3D::Zero())
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
        return convertMatrix<scalar_type, DIM>(m_estrain_curr);
    }

    static_matrix_type
    getTotalStrainPrev() const
    {
        return convertMatrix<scalar_type, DIM>(m_estrain_prev);
    }

    scalar_type
    getAccumulatedPlasticStrain() const
    {
        return scalar_type(0);
    }

    void
    update()
    {
        m_estrain_prev = m_estrain_curr;
    }

    static_matrix_type3D
    compute_stress3D(const data_type& data) const
    {
        return data.getLambda() * m_estrain_curr;
    }

    static_matrix_type
    compute_stress(const data_type& data) const
    {
        return convertMatrix<scalar_type, DIM>(compute_stress3D(data));
    }

    std::pair<static_matrix_type3D, static_tensor<scalar_type, 3>>
    compute_whole3D(const static_matrix_type3D& F_curr, const data_type& data, bool tangentmodulus = true)
    {
        static_tensor<scalar_type, 3> Cep = data.getLambda() * compute_IdentityTensor<scalar_type, 3>();

        // is always elastic
        m_estrain_curr = F_curr;

        // compute Cauchy stress
        const static_matrix_type3D stress = this->compute_stress3D(data);

        return std::make_pair(stress, Cep);
    }

    std::pair<static_matrix_type, static_tensor<scalar_type, DIM>>
    compute_whole(const static_matrix_type& F_curr, const data_type& data, bool tangentmodulus = true)
    {
        const static_matrix_type3D F3D         = convertMatrix3D(F_curr);
        const auto                 behaviors3D = compute_whole3D(F3D, data, tangentmodulus);

        const static_matrix_type              stress = convertMatrix<scalar_type, DIM>(behaviors3D.first);
        const static_tensor<scalar_type, DIM> Cep    = convertTensor<scalar_type, DIM>(behaviors3D.second);

        return std::make_pair(stress, Cep);
    }
};
}
