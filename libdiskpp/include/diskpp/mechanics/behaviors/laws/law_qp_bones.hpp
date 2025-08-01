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
 * Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework.
 * M. Abbas, A. Ern, N. Pignet.
 * International Journal of Numerical Methods in Engineering (2019)
 * 120(3), 303-327
 * DOI: 10.1002/nme.6137
 */

#pragma once

#include "diskpp/common/eigen.hpp"
#include "diskpp/mechanics/behaviors/laws/materialData.hpp"
#include "diskpp/mechanics/behaviors/maths_tensor.hpp"
#include "diskpp/mechanics/behaviors/maths_utils.hpp"
#include "diskpp/mechanics/behaviors/tensor_conversion.hpp"
#include "diskpp/mesh/point.hpp"
#include "diskpp/quadratures/quadrature_point.hpp"

namespace disk
{

// Bones for the computation of a behavior law at a quadrature point

/* pqrst is there only to squelch -fpermissive
 * GCC warnings. sorry for this. */

template<typename T, size_t DIM>
using pqrst = point<T,DIM>;

template<typename T, size_t DIM>
using qprst = quadrature_point<T,DIM>;

template<typename T, int DIM>
class law_qp_bones
{
  public:
    typedef T                                    scalar_type;
    typedef static_matrix<scalar_type, DIM, DIM> static_matrix_type;
    typedef static_matrix<scalar_type, 3, 3>     static_matrix_type3D;
    const static size_t                          dimension = DIM;
    typedef MaterialData<scalar_type>            data_type;

  protected:
    // coordinat and weight of considered gauss point.
    pqrst<scalar_type, DIM> m_point;
    scalar_type             m_weight;

    // internal variables at current step
    static_matrix_type3D m_estrain_prev; // elastic strain

    // internal variables at current step
    static_matrix_type3D m_estrain_curr; // elastic strain

    static_tensor<scalar_type, DIM>
    elastic_modulus(const data_type& data) const
    {

        return 2 * data.getMu() * IdentitySymTensor4<scalar_type, DIM>() + data.getLambda() * IxI<scalar_type, DIM>();
    }

    static_tensor<scalar_type, 3>
    elastic_modulus3D(const data_type& data) const
    {

        return 2 * data.getMu() * IdentitySymTensor4<scalar_type, 3>() + data.getLambda() * IxI<scalar_type, 3>();
    }

    scalar_type
    sigmaeq(const static_matrix_type3D& dev) const
    {
        return sqrt(scalar_type(1.5) * dev.squaredNorm());
    }

  public:
    law_qp_bones() :
      m_weight(0), m_estrain_prev(static_matrix_type3D::Zero()), m_estrain_curr(static_matrix_type3D::Zero())
    {
    }

    law_qp_bones(const point<scalar_type, DIM>& point, const scalar_type& weight) :
      m_point(point), m_weight(weight), m_estrain_prev(static_matrix_type3D::Zero()),
      m_estrain_curr(static_matrix_type3D::Zero())
    {
    }

    qprst<scalar_type, DIM>
    quadrature_point() const
    {
        return make_qp(m_point, m_weight);
    }

    auto
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
    getEquivalentPlasticStrain() const
    {
        return scalar_type(0);
    }

    void
    update()
    {
        m_estrain_prev = m_estrain_curr;
    }
};
}
