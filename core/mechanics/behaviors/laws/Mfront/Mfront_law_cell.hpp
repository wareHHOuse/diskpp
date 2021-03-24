/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2018, 2019                 nicolas.pignet@enpc.fr
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

#ifdef HAVE_MGIS

#include "core/mechanics/behaviors/laws/Mfront/Mfront_qp.hpp"
#include "core/quadratures/quadratures.hpp"

#include "MGIS/Behaviour/Behaviour.hxx"

namespace disk
{

/// Law cell bones

template<typename MeshType>
class Mfront_law_cell
{
  public:
    typedef MeshType                            mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell            cell_type;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    const static size_t dimension = mesh_type::dimension;

    typedef Mfront_qp<typename MeshType::coordinate_type, MeshType::dimension> law_qp_type;

  private:
    typedef std::shared_ptr<mgis::behaviour::Behaviour> BehaviourPtr;

    BehaviourPtr m_behav;

    std::vector<law_qp_type> m_list_qp;

  public:
    Mfront_law_cell() : m_behav(nullptr) {}

    Mfront_law_cell(const mesh_type& msh, const cell_type& cl, const size_t degree, const BehaviourPtr& b):
    m_behav(b)
    {
        const auto qps = integrate(msh, cl, degree);

        m_list_qp.clear();
        m_list_qp.reserve(qps.size());

        for (auto& qp : qps)
        {
            m_list_qp.push_back(law_qp_type(qp.point(), qp.weight(), b));
        }
    }

    int
    getNumberOfQP() const
    {
        return m_list_qp.size();
    }

    void
    update()
    {
        for (auto& qp : m_list_qp)
        {
            qp.update();
        }
    }

    std::vector<law_qp_type>&
    getQPs()
    {
        return m_list_qp;
    }

    const std::vector<law_qp_type>&
    getQPs() const
    {
        return m_list_qp;
    }

    law_qp_type&
    getQP(const size_t& qp_id)
    {
        return m_list_qp[qp_id];
    }

    const law_qp_type&
    getQP(const size_t& qp_id) const
    {
        return m_list_qp[qp_id];
    }

    std::vector<law_qp_type>
    getIVs() const
    {
        return m_list_qp;
    }
};
}

#endif
