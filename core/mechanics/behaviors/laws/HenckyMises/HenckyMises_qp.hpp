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

// Law for Linear Isotropic and Kinematic Hardening model with von Mises Criteria  in small

// Input : symetric stain tensor(Gs)

//           dev = normL2(Gs - trace(Gs) / dim * Id)

//             Stress : sigma = 2 *\tilde{mu}(dev(Gs)) * Gs + \tilde{lambda}(dev(Gs)) * trace(Gs) * Id
// \tilde{mu}(dev(Gs)) = mu * (1 + (1 + dev(Gs)) ^ {-1 / 2})
// \tilde{lambda}(dev(Gs)) = ((lambda + mu / 2) - mu / 2 * (1 + dev(Gs)) ^ {-1 / 2})

//                          Tangent Moduli : C = 2 * mu * I4 + lambda * prod_Kronecker(Id, Id) /
//                                                               it is the elastic moduli

template<typename T, int DIM>
class HenckyMises_qp : public law_qp_bones<T, DIM>
{
  public:
    typedef T                                    scalar_type;
    typedef static_matrix<scalar_type, DIM, DIM> static_matrix_type;
    typedef static_matrix<scalar_type, 3, 3>     static_matrix_type3D;
    const static size_t                          dimension = DIM;
    typedef MaterialData<scalar_type>            data_type;

  private:
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

  public:
    HenckyMises_qp() : law_qp_bones<T, DIM>() {}

    HenckyMises_qp(const point<scalar_type, DIM>& point, const scalar_type& weight) :
      law_qp_bones<T, DIM>(point, weight)
    {
    }

    static_matrix_type
    compute_stress(const data_type& data) const
    {
        const static_matrix_type Id = static_matrix_type::Identity();
        const static_matrix_type Gs = this->m_estrain_curr.block(0, 0, DIM, DIM);

        const scalar_type normFrodev = deviator(Gs).norm();

        return tildemu(data, normFrodev) * Gs + tildelambda(data, normFrodev) * Gs.trace() * Id;
    }

    std::pair<static_matrix_type, static_tensor<scalar_type, DIM>>
    compute_whole(const static_matrix_type& strain_curr, const data_type& data, bool tangentmodulus = true)
    {
        static_tensor<scalar_type, DIM> Cep = this->elastic_modulus(data);

        // is always elastic
        this->m_estrain_curr = convertMatrix3D(strain_curr);

        // compute Cauchy stress
        const static_matrix_type stress = this->compute_stress(data);

        return std::make_pair(stress, Cep);
    }
};
}
