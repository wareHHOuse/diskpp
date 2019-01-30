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
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */

#pragma once

#include <vector>

#include "core/bases/bases.hpp"
#include "core/common/eigen.hpp"
#include "core/mechanics/behaviors/laws/behaviorlaws.hpp"
#include "core/mechanics/behaviors/maths_tensor.hpp"
#include "core/mechanics/behaviors/maths_utils.hpp"
#include "core/quadratures/quadratures.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

namespace disk
{

/// Law cell bones

template<typename MeshType, typename LawTypeQp, bool PlasticBehavior>
class LawTypeCellBones
{
  public:
    typedef MeshType                            mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell            cell_type;

    typedef LawTypeQp                     law_qp_type;
    typedef typename LawTypeQp::data_type data_type;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    const static size_t dimension = mesh_type::dimension;

  private:
    std::vector<law_qp_type> m_list_qp;

  public:
    LawTypeCellBones(const mesh_type& msh, const cell_type& cl, const size_t degree)
    {
        const auto qps = integrate(msh, cl, degree);

        m_list_qp.clear();
        m_list_qp.reserve(qps.size());

        for (auto& qp : qps)
        {
            m_list_qp.push_back(law_qp_type(qp.point(), qp.weight()));
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

    std::vector<law_qp_type>
    getIVs() const
    {
        return m_list_qp;
    }

    vector_type
    projectStressOnCell(const mesh_type&             msh,
                        const cell_type&             cl,
                        const hho_degree_info&       hdi,
                        const data_type&             material_data) const
    {
        const auto grad_degree     = hdi.grad_degree();
        const int  grad_basis_size = matrix_basis_size(grad_degree, dimension, dimension);
        const auto gb              = make_matrix_monomial_basis(msh, cl, grad_degree);

        const matrix_type mass = make_mass_matrix(msh, cl, gb);
        vector_type rhs  = vector_type::Zero(grad_basis_size);

        for (auto& qp : m_list_qp)
        {
            const auto stress = qp.compute_stress(material_data);
            const auto gphi   = gb.eval_functions(qp.point());
            assert(gphi.size() == grad_basis_size);

            const auto qp_stress = priv::inner_product(qp.weight(), stress);

            rhs += priv::outer_product(gphi, qp_stress);
        }

        return mass.llt().solve(rhs);
    }

    vector_type
    projectPOnCell(const mesh_type& msh, const cell_type& cl, const hho_degree_info& hdi) const
    {
        const auto grad_degree = hdi.grad_degree();
        const int  pbs         = scalar_basis_size(grad_degree, dimension);

        if (PlasticBehavior)
        {
            const auto pb = make_scalar_monomial_basis(msh, cl, grad_degree);

            const matrix_type mass = make_mass_matrix(msh, cl, pb);
            vector_type rhs  = vector_type::Zero(pbs);

            for (auto& qp : m_list_qp)
            {
                const auto pphi = pb.eval_functions(qp.point());
                rhs += qp.weight() * qp.getAccumulatedPlasticStrain() * pphi;
            }

            return mass.llt().solve(rhs);
        }

        return vector_type::Zero(pbs);
    }

    vector_type
    projectStateOnCell(const mesh_type& msh, const cell_type& cl, const hho_degree_info& hdi) const
    {
        const auto grad_degree = hdi.grad_degree();
        const int  pbs         = scalar_basis_size(grad_degree, dimension);

        if (PlasticBehavior)
        {
            const auto pb = make_scalar_monomial_basis(msh, cl, grad_degree);

            const matrix_type mass = make_mass_matrix(msh, cl, pb);
            vector_type rhs  = vector_type::Zero(pbs);

            for (auto& qp : m_list_qp)
            {
                const auto pphi = pb.eval_functions(qp.point());
                assert(pphi.size() == pbs);

                if (qp.is_plastic())
                {
                    rhs += qp.weight() * pphi;
                }
            }

            return mass.llt().solve(rhs);
        }

        return vector_type::Zero(pbs);
    }
};
}
