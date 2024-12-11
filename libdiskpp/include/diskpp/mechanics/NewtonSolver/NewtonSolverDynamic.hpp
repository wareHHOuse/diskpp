/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019, 2025               nicolas.pignet@enpc.fr
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

#include <cassert>

#include "diskpp/bases/bases.hpp"
#include "diskpp/common/eigen.hpp"
#include "diskpp/common/timecounter.hpp"
#include "diskpp/mechanics/NewtonSolver/Fields.hpp"
#include "diskpp/mechanics/NewtonSolver/NewtonSolverParameters.hpp"
#include "diskpp/mechanics/NewtonSolver/TimeManager.hpp"
#include "diskpp/methods/hho"
#include "diskpp/quadratures/quadratures.hpp"

namespace disk
{

    namespace mechanics
    {

        template <typename MeshType>
        class dynamic_computation
        {
            typedef MeshType mesh_type;
            typedef typename mesh_type::coordinate_type scalar_type;
            typedef typename mesh_type::cell cell_type;

            typedef dynamic_matrix<scalar_type> matrix_type;
            typedef dynamic_vector<scalar_type> vector_type;

            bool m_enable;

            std::map<std::string, scalar_type> m_param;
            DynamicType m_scheme;

            std::vector<vector_type> m_acce_pred;

        public:
            matrix_type K_iner;
            vector_type RTF;

            double time_dyna;

            dynamic_computation(const NewtonSolverParameter<scalar_type> &rp)
            {
                m_enable = rp.isUnsteady();
                m_param = rp.getUnsteadyParameters();
                m_scheme = rp.getUnsteadyScheme();
            }

            bool
            enable(void) const
            {
                return m_enable;
            }

            void
            prediction(const mesh_type &mesh,
                       const MultiTimeField<scalar_type> &fields,
                       const TimeStep<scalar_type> &time_step)
            {
                if (this->enable())
                {
                    switch (m_scheme)
                    {
                    case DynamicType::NEWMARK:
                    {
                        m_acce_pred.clear();
                        m_acce_pred.reserve(mesh.cells_size());

                        const auto depl_prev = fields.getField(-1, FieldName::DEPL_CELLS);
                        const auto velo_prev = fields.getField(-1, FieldName::VITE_CELLS);
                        const auto acce_prev = fields.getField(-1, FieldName::ACCE_CELLS);

                        auto beta = m_param.at("beta");
                        scalar_type dt = time_step.increment_time();
                        scalar_type denom = 1.0 / (beta * dt * dt);
                        scalar_type nume = dt * dt / 2.0 * (1.0 - 2.0 * beta);

                        for (auto &cl : mesh)
                        {
                            const auto cl_id = mesh.lookup(cl);

                            auto aT_pred = denom * (depl_prev[cl_id] + dt * velo_prev[cl_id] + nume * acce_prev[cl_id]);
                            m_acce_pred.push_back(aT_pred);
                        }
                        break;
                    }
                    default:
                        break;
                    }
                }
            }

            scalar_type
            postprocess(const mesh_type &msh,
                        MultiTimeField<scalar_type> &fields,
                        const TimeStep<scalar_type> &time_step) const
            {
                timecounter tc;
                tc.tic();

                if (this->enable())
                {
                    switch (m_scheme)
                    {
                    case DynamicType::NEWMARK:
                    {
                        std::vector<vector_type> vite, acce;
                        vite.reserve(msh.cells_size());
                        acce.reserve(msh.cells_size());
                        const auto vite_prev = fields.getField(-1, FieldName::VITE_CELLS);
                        const auto acce_prev = fields.getField(-1, FieldName::ACCE_CELLS);

                        const auto depl_cells = fields.getCurrentField(FieldName::DEPL_CELLS);

                        auto beta = m_param.at("beta");
                        auto gamma = m_param.at("gamma");
                        scalar_type dt = time_step.increment_time();

                        const scalar_type denom = 1.0 / (beta * dt * dt);
                        const scalar_type nume = dt * dt / 2.0 * (1.0 - 2.0 * beta);
                        const scalar_type g0 = (1.0 - gamma) * dt, g1 = gamma * dt;

                        for (auto &cl : msh)
                        {
                            const auto cell_i = msh.lookup(cl);

                            auto aT = denom * depl_cells.at(cell_i) - m_acce_pred.at(cell_i);
                            auto vT = vite_prev.at(cell_i) + g0 * acce_prev.at(cell_i) + g1 * aT;

                            vite.push_back(vT);
                            acce.push_back(aT);
                        }
                        fields.setCurrentField(FieldName::VITE_CELLS, vite);
                        fields.setCurrentField(FieldName::ACCE_CELLS, acce);
                    }
                    default:
                        break;
                    }
                }

                tc.toc();
                return tc.elapsed();
            }

            void compute(const mesh_type &msh,
                         const cell_type &cl,
                         const MeshDegreeInfo<mesh_type> &degree_infos,
                         const vector_type &uTF,
                         const TimeStep<scalar_type> &time_step)
            {
                // Unsteady Computation

                time_dyna = 0.0;
                timecounter tc;

                tc.tic();
                if (this->enable())
                {
                    const auto cell_i = msh.lookup(cl);

                    const auto cell_infos = degree_infos.cellDegreeInfo(msh, cl);
                    const auto cell_degree = cell_infos.cell_degree();

                    const auto cb = make_vector_monomial_basis(msh, cl, cell_degree);

                    const auto faces_infos = cell_infos.facesDegreeInfo();

                    const auto num_cell_dofs = vector_basis_size(cell_degree, mesh_type::dimension, mesh_type::dimension);
                    const auto num_faces_dofs = vector_faces_dofs(msh, faces_infos);
                    const auto num_total_dofs = num_cell_dofs + num_faces_dofs;

                    const vector_type uT = uTF.head(num_cell_dofs);

                    RTF = vector_type::Zero(num_total_dofs);
                    K_iner = matrix_type::Zero(num_total_dofs, num_total_dofs);

                    switch (m_scheme)
                    {
                    case DynamicType::NEWMARK:
                    {

                        auto dt = time_step.increment_time();
                        auto beta = m_param.at("beta");
                        auto rho = m_param.at("rho");

                        const matrix_type mass_mat = rho * make_mass_matrix(msh, cl, cb);

                        const auto coeff = 1.0 / (beta * dt * dt);

                        K_iner.topLeftCorner(num_cell_dofs, num_cell_dofs) = coeff * mass_mat;

                        auto F_iner = coeff * (mass_mat * uT);
                        RTF.head(num_cell_dofs) -= F_iner;

                        RTF.head(num_cell_dofs) += mass_mat * m_acce_pred.at(cell_i);
                    }
                    default:
                        break;
                    }
                }
                tc.toc();
                time_dyna += tc.elapsed();
            }
        };
    }

} // end namespace diskpp