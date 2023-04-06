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

// Law for LinearLaw (test of finite deformations)

template<typename T, int DIM>
class LinearLaw_qp : public law_qp_bones<T, DIM>
{
  public:
    typedef T                                    scalar_type;
    typedef static_matrix<scalar_type, DIM, DIM> static_matrix_type;
    typedef static_matrix<scalar_type, 3, 3>     static_matrix_type3D;
    const static size_t                          dimension = DIM;
    typedef MaterialData<scalar_type>            data_type;

  public:
    LinearLaw_qp() : law_qp_bones<T, DIM>() {}

    LinearLaw_qp(const point<scalar_type, DIM>& point, const scalar_type& weight) : law_qp_bones<T, DIM>(point, weight)
    {
    }

    static_matrix_type3D
    compute_stress3D(const data_type& data) const
    {
        return data.getLambda() * this->m_estrain_curr;
    }

    static_matrix_type
    compute_stress(const data_type& data) const
    {
        return convertMatrix<scalar_type, DIM>(compute_stress3D(data));
    }

    std::pair<static_matrix_type3D, static_tensor<scalar_type, 3>>
    compute_whole3D(const static_matrix_type3D& F_curr, const data_type& data, bool tangentmodulus = true)
    {
        static_tensor<scalar_type, 3> Cep = data.getLambda() * IdentityTensor4<scalar_type, 3>();

        // is always elastic
        this->m_estrain_curr = F_curr;

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
