/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework.
 * M. Abbas, A. Ern, N. Pignet.
 * International Journal of Numerical Methods in Engineering (2019)
 * 120(3), 303-327
 * DOI: 10.1002/nme.6137
 */

#pragma once

#include "diskpp/common/eigen.hpp"
#include "diskpp/mechanics/behaviors/laws/law_qp_bones.hpp"
#include "diskpp/mechanics/behaviors/laws/materialData.hpp"
#include "diskpp/mechanics/behaviors/maths_tensor.hpp"
#include "diskpp/mechanics/behaviors/maths_utils.hpp"
#include "diskpp/mechanics/behaviors/tensor_conversion.hpp"
#include "diskpp/mesh/point.hpp"

namespace disk
{

// Law for linear elasticity

template<typename T, int DIM>
class LinearElasticity_qp : public law_qp_bones<T, DIM>
{
  public:
    const static size_t                          dimension = DIM;
    typedef T                                    scalar_type;
    typedef static_matrix<scalar_type, DIM, DIM> static_matrix_type;
    typedef static_matrix<scalar_type, 3, 3>     static_matrix_type3D;
    typedef MaterialData<scalar_type>            data_type;

    LinearElasticity_qp() : law_qp_bones<T, DIM>() {}

    LinearElasticity_qp(const point<scalar_type, DIM>& point, const scalar_type& weight) :
      law_qp_bones<T, DIM>(point, weight)
    {
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
        this->m_estrain_curr                  = strain_curr;
        const static_tensor<scalar_type, 3> C = this->elastic_modulus3D(data);

        // compute Cauchy stress
        const static_matrix_type3D stress = this->compute_stress3D(data);

        return std::make_pair(stress, C);
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
