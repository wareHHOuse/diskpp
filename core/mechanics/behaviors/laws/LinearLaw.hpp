/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
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
#include "mechanics/behaviors/maths_tensor.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

namespace disk
{

/* Material: LinearLaw
 * Stress :  Stress(Gs) = lambda * Gs
 * Module :  A(Gs) = lambda * I4
 */

template<typename scalar_type>
class LinearLaw
{
    scalar_type m_lambda;

  public:
    LinearLaw() : m_lambda(1.0) {}

    LinearLaw(const scalar_type lambda) : m_lambda(lambda) {}

    void
    setLambda(const scalar_type lambda)
    {
        m_lambda = lambda;
    }

    scalar_type
    giveLambda() const
    {
        return m_lambda;
    }

    template<int DIM>
    static_matrix<scalar_type, DIM, DIM>
    compute_stress(const static_matrix<scalar_type, DIM, DIM>& Gs) const
    {
        return m_lambda * Gs;
    }

    template<int DIM>
    static_tensor<scalar_type, DIM>
    compute_tangent_moduli(const static_matrix<scalar_type, DIM, DIM>& Gs) const
    {
        return m_lambda * compute_IdentityTensor<scalar_type, DIM>();
    }

    template<int DIM>
    std::pair<static_matrix<scalar_type, DIM, DIM>, static_tensor<scalar_type, DIM>>
    compute_whole(const static_matrix<scalar_type, DIM, DIM>& Gs) const
    {
        const static_matrix<scalar_type, DIM, DIM> stress = compute_stress(Gs);
        const static_tensor<scalar_type, DIM>      C      = compute_tangent_moduli(Gs);

        return std::make_pair(stress, C);
    }
};
}
