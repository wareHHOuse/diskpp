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

// NewtonRaphson_step

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "../Informations.hpp"
#include "../Parameters.hpp"
#include "mechanics/BoundaryConditions.hpp"
#include "plasticity_elementary_computation.hpp"
#include "revolution/bases"
#include "revolution/methods/hho"
#include "revolution/quadratures"

#include "solvers/solver.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

namespace NLE
{

template<typename MeshType>
class NewtonRaphson_step_plasticity
{
    typedef MeshType                                       mesh_type;
    typedef typename mesh_type::scalar_type                scalar_type;
    typedef ParamRun<scalar_type>                          param_type;
    typedef typename revolution::hho_degree_info           hdi_type;
    typedef disk::mechanics::BoundaryConditions<mesh_type> bnd_type;

    const static int dimension = mesh_type::dimension;

    typedef dynamic_matrix<scalar_type> matrix_dynamic;
    typedef dynamic_vector<scalar_type> vector_dynamic;

    typedef revolution::assembler_mechanics<mesh_type> assembler_type;

    typedef plasticity<mesh_type> elem_type;

    vector_dynamic m_system_solution;

    const bnd_type&   m_bnd;
    const mesh_type&  m_msh;
    const hdi_type&   m_hdi;
    const param_type& m_rp;
    assembler_type    m_assembler;

    std::vector<vector_dynamic> m_bL;
    std::vector<matrix_dynamic> m_AL;

    std::vector<vector_dynamic> m_postprocess_data, m_solution_data;
    std::vector<vector_dynamic> m_solution_cells, m_solution_faces;

    scalar_type m_F_int;

    bool m_verbose;

  public:
    NewtonRaphson_step_plasticity(const mesh_type&  msh,
                                  const hdi_type&   hdi,
                                  const bnd_type&   bnd,
                                  const param_type& rp) :
      m_msh(msh),
      m_verbose(rp.m_verbose), m_rp(rp), m_hdi(hdi), m_bnd(bnd)
    {
        m_AL.clear();
        m_AL.resize(m_msh.cells_size());

        m_bL.clear();
        m_bL.resize(m_msh.cells_size());

        m_assembler = revolution::make_mechanics_assembler(m_msh, m_hdi, m_bnd);
    }

    bool
    verbose(void) const
    {
        return m_verbose;
    }
    void
    verbose(bool v)
    {
        m_verbose = v;
    }

    void
    initialize(const std::vector<vector_dynamic>& initial_solution_cells,
               const std::vector<vector_dynamic>& initial_solution_faces,
               const std::vector<vector_dynamic>& initial_solution)
    {
        m_solution_cells.clear();
        m_solution_cells = initial_solution_cells;
        assert(m_msh.cells_size() == m_solution_cells.size());

        m_solution_faces.clear();
        m_solution_faces = initial_solution_faces;
        assert(m_msh.faces_size() == m_solution_faces.size());

        m_solution_data.clear();
        m_solution_data = initial_solution;
        assert(m_msh.cells_size() == m_solution_data.size());
    }

    template<typename LoadFunction, typename Law>
    AssemblyInfo
    assemble(const LoadFunction&                lf,
             const std::vector<matrix_dynamic>& gradient_precomputed,
             const std::vector<matrix_dynamic>& stab_precomputed,
             Law&                               law,
             bool                               elastic_modulus = false)
    {
        typename revolution::static_condensation_vector<mesh_type> statcond;

        elem_type    elem(m_msh, m_hdi);
        AssemblyInfo ai;

        // set RHS to zero
        m_assembler.setZeroRhs();

        const auto material_data = law.getMaterialData();
        m_F_int                  = 0.0;

        timecounter tc, ttot;

        ttot.tic();
        int cell_i = 0;

        for (auto& cl : m_msh)
        {
            // Gradient Reconstruction
            matrix_dynamic GT;
            tc.tic();
            if (m_rp.m_precomputation)
            {
                GT = gradient_precomputed[cell_i];
            }
            else
            {
                const auto gradrec_full = make_hho_sym_gradrec_matrix(m_msh, cl, m_hdi);
                GT                      = gradrec_full.first;
            }
            tc.toc();
            ai.m_time_gradrec += tc.to_double();

            // Begin Assembly
            // Build rhs and lhs

            // Mechanical Computation

            auto& law_cell = law.getCellQPs(cell_i);
            tc.tic();
            elem.compute(cl, lf, GT, m_solution_data.at(cell_i), law_cell, material_data, elastic_modulus);

            matrix_dynamic lhs = elem.K_int;
            vector_dynamic rhs = elem.RTF;
            m_F_int += elem.F_int.squaredNorm();

            tc.toc();
            ai.m_time_elem += tc.to_double();
            ai.m_time_law += elem.time_law;

            // Stabilisation Contribution
            tc.tic();

            if (m_rp.m_stab)
            {
                if (m_rp.m_precomputation)
                {
                    const matrix_dynamic stab = stab_precomputed.at(cell_i);
                    assert(elem.K_int.rows() == stab.rows());
                    assert(elem.K_int.cols() == stab.cols());
                    assert(elem.RTF.rows() == (stab * m_solution_data.at(cell_i)).rows());
                    assert(elem.RTF.cols() == (stab * m_solution_data.at(cell_i)).cols());

                    lhs += m_rp.m_beta * stab;
                    rhs -= m_rp.m_beta * stab * m_solution_data.at(cell_i);
                }
                else
                {
                    switch (m_rp.m_stab_type)
                    {
                        case HHO:
                        {
                            const auto recons   = make_hho_vector_symmetric_laplacian(m_msh, cl, m_hdi);
                            const auto stab_HHO = make_hho_vector_stabilization(m_msh, cl, recons.first, m_hdi);
                            assert(elem.K_int.rows() == stab_HHO.rows());
                            assert(elem.K_int.cols() == stab_HHO.cols());
                            assert(elem.RTF.rows() == (stab_HHO * m_solution_data.at(cell_i)).rows());
                            assert(elem.RTF.cols() == (stab_HHO * m_solution_data.at(cell_i)).cols());

                            lhs += m_rp.m_beta * stab_HHO;
                            rhs -= m_rp.m_beta * stab_HHO * m_solution_data.at(cell_i);
                            break;
                        }
                        case HDG: {
                            const auto stab_HDG = make_hdg_vector_stabilization(m_msh, cl, m_hdi);
                            assert(elem.K_int.rows() == stab_HDG.rows());
                            assert(elem.K_int.cols() == stab_HDG.cols());
                            assert(elem.RTF.rows() == (stab_HDG * m_solution_data.at(cell_i)).rows());
                            assert(elem.RTF.cols() == (stab_HDG * m_solution_data.at(cell_i)).cols());

                            lhs += m_rp.m_beta * stab_HDG;
                            rhs -= m_rp.m_beta * stab_HDG * m_solution_data.at(cell_i);
                            break;
                        }
                        case NO: { break;
                        }
                        default: throw std::invalid_argument("Unknown stabilization");
                    }
                }
            }
            tc.toc();
            ai.m_time_stab += tc.to_double();

            // Static Condensation
            tc.tic();
            auto scnp = statcond.compute_rhsfull(m_msh, cl, lhs, rhs, m_hdi);

            m_AL[cell_i] = statcond.AL;
            m_bL[cell_i] = statcond.bL;

            tc.toc();
            ai.m_time_statcond += tc.to_double();

            m_assembler.assemble_nl(m_msh, cl, m_bnd, scnp, m_solution_faces);

            cell_i++;
        }

        m_F_int = sqrt(m_F_int);

        m_assembler.impose_neumann_boundary_conditions(m_msh, m_bnd);
        m_assembler.finalize();

        ttot.toc();
        ai.m_time_assembly      = ttot.to_double();
        ai.m_linear_system_size = m_assembler.LHS.rows();
        return ai;
    }

    SolveInfo
    solve(void)
    {
        timecounter tc;

        tc.tic();
        m_system_solution = vector_dynamic::Zero(m_assembler.LHS.rows());

        disk::solvers::pardiso_params<scalar_type> pparams;
        mkl_pardiso(pparams, m_assembler.LHS, m_assembler.RHS, m_system_solution);
        tc.toc();

        return SolveInfo(m_assembler.LHS.rows(), m_assembler.LHS.nonZeros(), tc.to_double());
    }

    scalar_type
    postprocess()
    {
        timecounter tc;
        tc.tic();

        const int fbs = revolution::vector_basis_size(m_hdi.face_degree(), dimension - 1, dimension);
        const int cbs = revolution::vector_basis_size(m_hdi.cell_degree(), dimension, dimension);

        const auto solF = m_assembler.expand_solution_nl(m_msh, m_bnd, m_system_solution, m_solution_faces);

        // Update  unknowns
        // Update face Uf^{i+1} = Uf^i + delta Uf^i
        for (int i = 0; i < m_solution_faces.size(); i++)
        {
            assert(m_solution_faces.at(i).size() == fbs);
            m_solution_faces.at(i) += solF.segment(i * fbs, fbs);
        }
        // Update cell
        int cell_i = 0;
        for (auto& cl : m_msh)
        {
            // Extract the solution
            const auto fcs             = faces(m_msh, cl);
            const auto num_faces       = fcs.size();
            const auto total_faces_dof = num_faces * fbs;

            vector_dynamic xFs = vector_dynamic::Zero(total_faces_dof);

            for (int face_i = 0; face_i < num_faces; face_i++)
            {
                const auto fc  = fcs[face_i];
                const auto eid = find_element_id(m_msh.faces_begin(), m_msh.faces_end(), fc);
                if (!eid.first)
                    throw std::invalid_argument("This is a bug: face not found");

                const auto face_id = eid.second;

                xFs.segment(face_i * fbs, fbs) = solF.segment(face_id * fbs, fbs);
            }

            const vector_dynamic xT = m_bL[cell_i] - m_AL[cell_i] * xFs;

            assert(xT.size() == cbs);
            assert(m_solution_data.at(cell_i).size() == xT.size() + xFs.size());
            // Update element U^{i+1} = U^i + delta U^i ///
            (m_solution_data.at(cell_i)).segment(0, cbs) += xT;
            (m_solution_data.at(cell_i)).segment(cbs, total_faces_dof) += xFs;

            // Update Cell Uc^{i+1} = Uc^i + delta Uc^i ///
            assert(m_solution_cells.at(cell_i).size() == cbs);
            m_solution_cells.at(cell_i) += xT;

            cell_i++;
        }

        tc.toc();
        return tc.to_double();
    }

    bool
    test_convergence(const int iter)
    {
        // norm of the solution
        scalar_type error_un = 0;
        for (int i = 0; i < m_solution_faces.size(); i++)
        {
            scalar_type norm = m_solution_faces[i].norm();
            error_un += norm * norm;
        }

        error_un = std::sqrt(error_un);

        if (error_un <= scalar_type(10E-15))
        {
            error_un = scalar_type(10E16);
        }

        // norm of the rhs
        const scalar_type residual  = m_assembler.RHS.norm();
        scalar_type       max_error = 0.0;
        for (int i = 0; i < m_assembler.RHS.size(); i++)
            max_error = std::max(max_error, std::abs(m_assembler.RHS(i)));

        // norm of the increment
        const scalar_type error_incr     = m_system_solution.norm();
        scalar_type       relative_displ = error_incr / error_un;
        scalar_type       relative_error = residual / m_F_int;

        if (iter == 0)
        {
            relative_displ = 1;
            relative_error = 1;
        }

        if (m_verbose)
        {
            std::string s_iter = "   " + std::to_string(iter) + "               ";
            s_iter.resize(9);

            if (iter == 0)
            {
                std::cout << "----------------------------------------------------------------------"
                             "------------------------"
                          << std::endl;
                std::cout << "| Iteration | Norme l2 incr | Relative incr |  Residual l2  | "
                             "Relative error | Maximum error |"
                          << std::endl;
                std::cout << "----------------------------------------------------------------------"
                             "------------------------"
                          << std::endl;
            }
            std::ios::fmtflags f(std::cout.flags());
            std::cout.precision(5);
            std::cout.setf(std::iostream::scientific, std::iostream::floatfield);
            std::cout << "| " << s_iter << " |   " << error_incr << " |   " << relative_displ << " |   " << residual
                      << " |   " << relative_error << "  |  " << max_error << "  |" << std::endl;
            std::cout << "-------------------------------------------------------------------------"
                         "---------------------"
                      << std::endl;
            std::cout.flags(f);
        }

        const scalar_type error = std::max(relative_displ, relative_error);

        if (error <= m_rp.m_epsilon)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    void
    save_solutions(std::vector<vector_dynamic>& solution_cells,
                   std::vector<vector_dynamic>& solution_faces,
                   std::vector<vector_dynamic>& solution)
    {
        solution_cells.clear();
        solution_cells = m_solution_cells;
        assert(m_solution_cells.size() == solution_cells.size());

        solution_faces.clear();
        solution_faces = m_solution_faces;
        assert(m_solution_faces.size() == solution_faces.size());

        solution.clear();
        solution = m_solution_data;
        assert(m_solution_data.size() == solution.size());
    }
};

} // end NLE