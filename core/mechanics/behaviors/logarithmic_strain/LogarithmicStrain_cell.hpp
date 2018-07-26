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

#include <vector>

#include "bases/bases.hpp"
#include "common/eigen.hpp"
#include "core/mechanics/behaviors/laws/behaviorlaws.hpp"
#include "mechanics/behaviors/logarithmic_strain/LogarithmicStrain_qp.hpp"
// #include "mechanics/behaviors/maths_tensor.hpp"
// #include "mechanics/behaviors/maths_utils.hpp"
#include "quadratures/quadratures.hpp"

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

template<typename LawTypeCell>
class LogarithmicStrain_cell
{
  public:
    typedef LawTypeCell                             law_hpp_cell_type;
    typedef typename law_hpp_cell_type::law_qp_type law_hpp_qp_type;
    typedef typename law_hpp_cell_type::mesh_type   mesh_type;
    typedef typename mesh_type::coordinate_type     scalar_type;
    typedef typename mesh_type::cell                cell_type;
    typedef typename law_hpp_qp_type::data_type     data_type;

    typedef typename disk::mechanics::LogarithmicStrain_qp<law_hpp_qp_type> law_qp_type;

  private:
    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    const static size_t dimension = mesh_type::dimension;

    std::vector<law_qp_type> m_list_qp;

  public:
    LogarithmicStrain_cell(const mesh_type& msh, const cell_type& cl, const size_t degree)
    {
        const auto qps = disk::integrate(msh, cl, degree);

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

    std::vector<law_qp_type>
    getIVs() const
    {
        return m_list_qp;
    }

    vector_type
    projectStressOnCell(const mesh_type&             msh,
                        const cell_type&             cl,
                        const disk::hho_degree_info& hdi,
                        const data_type&             material_data) const
    {
        const auto grad_degree     = hdi.grad_degree();
        const int  grad_basis_size = disk::sym_matrix_basis_size(grad_degree, dimension, dimension);
        auto       gb              = disk::make_sym_matrix_monomial_basis(msh, cl, grad_degree);

        matrix_type mass = matrix_type::Zero(grad_basis_size, grad_basis_size);
        vector_type rhs  = vector_type::Zero(grad_basis_size);

        for (auto& qp : m_list_qp)
        {
            const auto gphi = gb.eval_functions(qp.point());
            assert(gphi.size() == grad_basis_size);

            for (int j = 0; j < grad_basis_size; j++)
            {
                const auto qp_gphi_j = disk::priv::inner_product(qp.weight(), gphi[j]);
                for (int i = j; i < grad_basis_size; i++)
                {
                    mass(i, j) += disk::priv::inner_product(gphi[i], qp_gphi_j);
                }

                rhs(j) += disk::priv::inner_product(qp.compute_stress(material_data), qp_gphi_j);
            }
        }

        // lower part
        for (int j = 0; j < grad_basis_size; j++)
        {
            for (int i = 0; i < j; i++)
            {
                mass(i, j) = mass(j, i);
            }
        }

        return mass.llt().solve(rhs);
    }

    vector_type
    projectPOnCell(const mesh_type& msh, const cell_type& cl, const disk::hho_degree_info& hdi) const
    {
        const auto grad_degree = hdi.grad_degree();
        auto       pb          = disk::make_scalar_monomial_basis(msh, cl, grad_degree);
        const int  pbs         = disk::scalar_basis_size(grad_degree, dimension);

        matrix_type mass = matrix_type::Zero(pbs, pbs);
        vector_type rhs  = vector_type::Zero(pbs);

        for (auto& qp : m_list_qp)
        {
            const auto pphi = pb.eval_functions(qp.point());
            assert(pphi.size() == pbs);

            mass += qp.weight() * disk::priv::outer_product(pphi, pphi);
            rhs += qp.weight() * qp.getAccumulatedPlasticStrain() * pphi;
        }

        return mass.llt().solve(rhs);
    }

    vector_type
    projectStateOnCell(const mesh_type& msh, const cell_type& cl, const disk::hho_degree_info& hdi) const
    {
        const auto grad_degree = hdi.grad_degree();
        auto       pb          = disk::make_scalar_monomial_basis(msh, cl, grad_degree);
        const int  pbs         = disk::scalar_basis_size(grad_degree, dimension);

        matrix_type mass = matrix_type::Zero(pbs, pbs);
        vector_type rhs  = vector_type::Zero(pbs);

        for (auto& qp : m_list_qp)
        {
            const auto pphi = pb.eval_functions(qp.point());
            assert(pphi.size() == pbs);

            mass += qp.weight() * disk::priv::outer_product(pphi, pphi);
            if (qp.is_plastic())
            {
                rhs += qp.weight() * pphi;
            }
        }

        return mass.llt().solve(rhs);
    }
};
}
}