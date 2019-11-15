/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
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

// NewtonRaphson iteration

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "NewtonSolverComput.hpp"
#include "NewtonSolverInformations.hpp"
#include "NewtonSolverParameters.hpp"
#include "bases/bases.hpp"

#include "adaptivity/adaptivity.hpp"
#include "boundary_conditions/boundary_conditions.hpp"
#include "mechanics/behaviors/laws/behaviorlaws.hpp"
#include "methods/hho"

#include "solvers/solver.hpp"

#include "timecounter.h"

namespace disk
{

namespace mechanics
{

/**
 * @brief Newton-Raphson iteration for nonlinear solid mechanics
 *
 *  Specialized for HHO methods
 *
 *  Options :  - small and finite deformations
 *             - plasticity, hyperelasticity (various laws)
 *
 * @tparam MeshType type of the mesh
 */
template<typename MeshType>
class NewtonIteration
{
    typedef MeshType                            mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    typedef NewtonSolverParameter<scalar_type>    param_type;
    typedef vector_boundary_conditions<mesh_type> bnd_type;
    typedef hho_degree_info                       hdi_type;
    typedef Behavior<mesh_type>                   behavior_type;

    typedef vector_primal_hho_assembler<mesh_type> assembler_type;
    typedef mechanical_computation<mesh_type>      elem_type;

    vector_type m_system_solution;

    const mesh_type&  m_msh;
    const hdi_type&   m_hdi;
    const bnd_type&   m_bnd;
    const param_type& m_rp;
    assembler_type    m_assembler;

    std::vector<vector_type> m_bL;
    std::vector<matrix_type> m_AL;

    std::vector<vector_type> m_postprocess_data;
    std::vector<vector_type> m_solution, m_solution_faces;

    scalar_type m_F_int;

    bool m_verbose;

  public:
    NewtonIteration(const mesh_type& msh, const hdi_type& hdi, const bnd_type& bnd, const param_type& rp) :
      m_msh(msh), m_hdi(hdi), m_rp(rp), m_bnd(bnd), m_verbose(rp.m_verbose)
    {
        m_AL.clear();
        m_AL.resize(m_msh.cells_size());

        m_bL.clear();
        m_bL.resize(m_msh.cells_size());

        m_assembler = make_vector_primal_hho_assembler(m_msh, m_hdi, m_bnd);
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
    initialize(const std::vector<vector_type>& initial_solution, const std::vector<vector_type>& initial_solution_faces)
    {
        m_solution_faces.clear();
        m_solution_faces = initial_solution_faces;
        assert(m_msh.faces_size() == m_solution_faces.size());

        m_solution.clear();
        m_solution = initial_solution;
        assert(m_msh.cells_size() == m_solution.size());
    }

    template<typename LoadFunction>
    AssemblyInfo
    assemble(const LoadFunction&             lf,
             const std::vector<matrix_type>& gradient_precomputed,
             const std::vector<matrix_type>& stab_precomputed,
             const MeshDegree<mesh_type>&    degree_infos,
             behavior_type&                  behavior)
    {
        elem_type    elem(m_msh, m_hdi);
        AssemblyInfo ai;

        // set RHS to zero
        m_assembler.initialize();
        m_F_int = 0.0;

        // Get materail data
        const auto material_data = behavior.getMaterialData();
        // Get Law
        auto& law = behavior.law();

        const bool small_def = (behavior.getDeformation() == SMALL_DEF);

        timecounter tc, ttot;

        ttot.tic();
        size_t cell_i = 0;

        for (auto& cl : m_msh)
        {
            // Gradient Reconstruction
            matrix_type GT;
            tc.tic();
            if (m_rp.m_precomputation)
            {
                GT = gradient_precomputed[cell_i];
            }
            else
            {
                if (small_def)
                {
                    const auto gradrec_sym_full = make_matrix_symmetric_gradrec(m_msh, cl, m_hdi);
                    GT                          = gradrec_sym_full.first;
                }
                else
                {
                    const auto gradrec_full = make_marix_hho_gradrec(m_msh, cl, m_hdi);
                    GT                      = gradrec_full.first;
                }
            }
            tc.toc();
            ai.m_time_gradrec += tc.to_double();

            // Begin Assembly
            // Build rhs and lhs

            // Mechanical Computation

            auto& law_cell = law.getCellQPs(cell_i);
            //auto& law_qp   = behavior.getQPs(m_msh, cl);
            tc.tic();
            elem.compute(cl, lf, GT, m_solution.at(cell_i), law_cell, material_data, small_def);

            matrix_type lhs = elem.K_int;
            vector_type rhs = elem.RTF;
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
                    const matrix_type& stab = stab_precomputed.at(cell_i);
                    assert(elem.K_int.rows() == stab.rows());
                    assert(elem.K_int.cols() == stab.cols());
                    assert(elem.RTF.rows() == stab.rows());
                    assert(elem.RTF.cols() == m_solution.at(cell_i).cols());

                    lhs += m_rp.m_beta * stab;
                    rhs -= m_rp.m_beta * stab * m_solution.at(cell_i);
                }
                else
                {
                    switch (m_rp.m_stab_type)
                    {
                        case HHO:
                        {
                            matrix_type stab_HHO;
                            if (small_def)
                            {
                                const auto recons = make_vector_hho_symmetric_laplacian(m_msh, cl, m_hdi);
                                stab_HHO          = make_vector_hho_stabilization(m_msh, cl, recons.first, m_hdi);
                            }
                            else
                            {
                                const auto recons_scalar = make_scalar_hho_laplacian(m_msh, cl, m_hdi);
                                stab_HHO  = make_vector_hho_stabilization_optim(m_msh, cl, recons_scalar.first, m_hdi);
                            }

                            assert(elem.K_int.rows() == stab_oper.rows());
                            assert(elem.K_int.cols() == stab_oper.cols());
                            assert(elem.RTF.rows() == stab_oper.rows());
                            assert(elem.RTF.cols() == m_solution.at(cell_i).cols());

                            lhs += m_rp.m_beta * stab_HHO;
                            rhs -= m_rp.m_beta * stab_HHO * m_solution.at(cell_i);
                            break;
                        }
                        case HDG:
                        {
                            const auto stab_HDG = make_vector_hdg_stabilization(m_msh, cl, m_hdi);
                            assert(elem.K_int.rows() == stab_HDG.rows());
                            assert(elem.K_int.cols() == stab_HDG.cols());
                            assert(elem.RTF.rows() == stab_HDG.rows());
                            assert(elem.RTF.cols() == m_solution.at(cell_i).cols());

                            lhs += m_rp.m_beta * stab_HDG;
                            rhs -= m_rp.m_beta * stab_HDG * m_solution.at(cell_i);
                            break;
                        }
                        case DG:
                        {
                            const auto stab_DG = make_vector_dg_stabilization(m_msh, cl, m_hdi);
                            assert(elem.K_int.rows() == stab_DG.rows());
                            assert(elem.K_int.cols() == stab_DG.cols());
                            assert(elem.RTF.rows() == stab_DG.rows());
                            assert(elem.RTF.cols() == m_solution.at(cell_i).cols());

                            lhs += m_rp.m_beta * stab_DG;
                            rhs -= m_rp.m_beta * stab_DG * m_solution.at(cell_i);
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
            const auto scnp = make_vector_static_condensation_withMatrix(m_msh, cl, m_hdi, lhs, rhs);

            m_AL[cell_i] = std::get<1>(scnp);
            m_bL[cell_i] = std::get<2>(scnp);

            tc.toc();
            ai.m_time_statcond += tc.to_double();

            const auto& lc = std::get<0>(scnp);
            m_assembler.assemble_nonlinear(m_msh, cl, m_bnd, lc.first, lc.second, m_solution_faces);

            cell_i++;
        }

        m_F_int = sqrt(m_F_int);
        //std::cout << "F_int: " << m_F_int << std::endl;

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
        m_system_solution = vector_type::Zero(m_assembler.LHS.rows());

        solvers::pardiso_params<scalar_type> pparams;
        mkl_pardiso(pparams, m_assembler.LHS, m_assembler.RHS, m_system_solution);
        tc.toc();

        return SolveInfo(m_assembler.LHS.rows(), m_assembler.LHS.nonZeros(), tc.to_double());
    }

    scalar_type
    postprocess(const MeshDegree<mesh_type>&  degree_infos)
    {
        timecounter tc;
        tc.tic();

        // Update cell
        size_t cell_i = 0;
        for (auto& cl : m_msh)
        {
            const vector_type xdT = m_assembler.take_local_solution_nonlinear(m_msh, cl, m_bnd, m_system_solution, m_solution_faces);

            // static decondensation
            const vector_type xT = m_bL[cell_i] - m_AL[cell_i] * xdT;

            assert(m_solution.at(cell_i).size() == xT.size() + xdT.size());
            // Update element U^{i+1} = U^i + delta U^i ///
            m_solution.at(cell_i).head(xT.size()) += xT;
            m_solution.at(cell_i).segment(xT.size(), xdT.size()) += xdT;

            // std::cout << "KT_F " << m_AL[cell_i].norm() << std::endl;
            // std::cout << "sol_F" << std::endl;
            // std::cout << xdT.transpose() << std::endl;
            // std::cout << "ft" << std::endl;
            // std::cout << m_bL[cell_i].transpose() << std::endl;
            // std::cout << "sol_T" << std::endl;
            // std::cout << xT.transpose() << std::endl;
            // std::cout << (m_solution.at(cell_i)).segment(0, cbs).transpose() << std::endl;

            cell_i++;
        }

        const auto solF = m_assembler.expand_solution_nonlinear(m_msh, m_bnd, m_system_solution, m_solution_faces);

        // Update  unknowns
        // Update face Uf^{i+1} = Uf^i + delta Uf^i
        const auto fbs = vector_basis_size(m_hdi.face_degree(), mesh_type::dimension - 1, mesh_type::dimension);
        for (size_t i = 0; i < m_solution_faces.size(); i++)
        {
            assert(m_solution_faces.at(i).size() == fbs);
            m_solution_faces.at(i) += solF.segment(i * fbs, fbs);
        }

        tc.toc();
        return tc.to_double();
    }

    bool
    convergence(const size_t iter)
    {
        // norm of the solution
        scalar_type error_un = 0;
        for (size_t i = 0; i < m_solution_faces.size(); i++)
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
        for (size_t i = 0; i < m_assembler.RHS.size(); i++)
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
    save_solutions(std::vector<vector_type>& solution, std::vector<vector_type>& solution_faces)
    {
        solution_faces.clear();
        solution_faces = m_solution_faces;
        assert(m_solution_faces.size() == solution_faces.size());

        solution.clear();
        solution = m_solution;
        assert(m_solution.size() == solution.size());
    }
};
}
} // end diskpp