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
#include "core/mechanics/behaviors/laws/behaviorlaws.hpp"
// #include "mechanics/behaviors/maths_tensor.hpp"
// #include "mechanics/behaviors/maths_utils.hpp"
#include "mesh/point.hpp"



namespace disk
{

namespace mechanics
{

// Routine for Logarithmic Stain

/* For details see the paper:
 *   Anisotropic additive plasticity in the logaritmic strain: modular
 *   kinematic formulation and implementation based on incremental minimization
 *   principles for standard materials
 *   C. Miehe, N. Apel, M. Lambrecht
 *   Comput. Methods Appl. Mech. Engrg. (2002)
 */

template<typename LawTypeQp>
class LogarithmicStrain_qp
{
  public:
    typedef LawTypeQp                             law_hpp_qp_type;
    typedef typename law_hpp_qp_type::data_type   data_type;
    typedef typename law_hpp_qp_type::scalar_type scalar_type;

    const static size_t DIM = LawTypeQp::dimension;

    typedef static_matrix<scalar_type, DIM, DIM> static_matrix_type;
    typedef static_matrix<scalar_type, 3, 3>     static_matrix_type3D;

  private:
    // law hpp at qp
    law_hpp_qp_type m_law_hpp_qp;

  public:
    LogarithmicStrain_qp(const point<scalar_type, DIM>& point, const scalar_type& weight)
    {
        m_law_hpp_qp = law_hpp_qp_type(point, weight);
    }

    point<scalar_type, DIM>
    point() const
    {
        return m_law_hpp_qp.point();
    }

    scalar_type
    weight() const
    {
        return m_law_hpp_qp.weight();
    }

    bool
    is_plastic() const
    {
        return m_law_hpp_qp.is_plastic();
    }

    static_matrix_type3D
    getElasticStrain() const
    {
        return m_law_hpp_qp.getElasticStrain();
    }

    static_matrix_type3D
    getPlasticStrain() const
    {
        return m_law_hpp_qp.getPlasticStrain();
    }

    static_matrix_type
    getTotalStrain() const
    {
        return m_law_hpp_qp.getTotalStrain();
    }

    static_matrix_type
    getTotalStrainPrev() const
    {
        return m_law_hpp_qp.getTotalStrainPrev();
    }

    scalar_type
    getAccumulatedPlasticStrain() const
    {
        return m_law_hpp_qp.getAccumulatedPlasticStrain();
    }

    void
    update()
    {
        m_law_hpp_qp.update();
    }

    static_matrix_type
    compute_stress(const data_type& data) const
    {
        return m_law_hpp_qp.compute_stress(data);
    }

    std::pair<static_matrix_type, static_tensor<scalar_type, DIM>>
    compute_whole(const static_matrix_type& incr_strain, const data_type& data, bool tangentmodulus = true)
    {
        return m_law_hpp_qp.compute_whole(incr_strain, data, tangentmodulus);
    }
};
}
}