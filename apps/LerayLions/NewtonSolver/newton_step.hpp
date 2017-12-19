/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
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
#include "LerayLions_elementary_computation.hpp"
#include "hho/assembler.hpp"
#include "hho/gradient_reconstruction.hpp"
#include "hho/hho_bq.hpp"
#include "hho/stabilization.hpp"
#include "hho/static_condensation.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

template<typename BQData>
class NewtonRaphson_step_leraylions
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef ParamRun<scalar_type>           param_type;

   typedef dynamic_matrix<scalar_type> matrix_dynamic;
   typedef dynamic_vector<scalar_type> vector_dynamic;

   typedef disk::hho::gradient_reconstruction_bq<BQData>      gradrec_type;
   typedef disk::hho::gradient_reconstruction_full_bq<BQData> gradrec_full_type;
   typedef disk::hho::static_condensation_bq<BQData>          statcond_type;
   typedef disk::hho::assembler_by_elimination_bq<BQData>     assembler_type;

   typedef disk::hho::hho_stabilization_bq<BQData>  hho_stab_type;
   typedef disk::hho::hdg_stabilization_bq<BQData>  hdg_stab_type;
   typedef disk::hho::pikF_stabilization_bq<BQData> pikF_stab_type;

   typedef LerayLions::LerayLions<BQData> elem_type;

   typedef typename assembler_type::sparse_matrix_type sparse_matrix_type;
   typedef typename assembler_type::vector_type        vector_type;

#ifdef HAVE_INTEL_MKL
   typedef Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>> solver_type;
#else
   typedef Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>> solver_type;
#endif

   sparse_matrix_type m_system_matrix;
   vector_type        m_system_rhs, m_system_solution;

   solver_type solver;

   const BQData&     m_bqd;
   const param_type& m_rp;
   const mesh_type&  m_msh;

   std::vector<vector_dynamic> m_RT;
   std::vector<matrix_dynamic> m_KTT, m_KTF;

   std::vector<vector_dynamic> m_postprocess_data, m_solution_data;
   std::vector<vector_dynamic> m_solution_cells, m_solution_faces;

   bool m_verbose;

 public:
   NewtonRaphson_step_leraylions(const mesh_type& msh, const BQData& bqd, const param_type& rp) :
     m_msh(msh), m_verbose(rp.m_verbose), m_rp(rp), m_bqd(bqd)
   {
      m_RT.clear();
      m_RT.resize(m_msh.cells_size());

      m_KTF.clear();
      m_KTF.resize(m_msh.cells_size());

      m_KTT.clear();
      m_KTT.resize(m_msh.cells_size());
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

   template<typename LoadFunction, typename BoundaryConditionFunction>
   AssemblyInfo
   assemble(const LoadFunction&                lf,
            const BoundaryConditionFunction&   bcf,
            const std::vector<matrix_dynamic>& gradient_precomputed)
   {
      gradrec_type      gradrec(m_bqd);
      gradrec_full_type gradrec_full(m_bqd);
      elem_type         elem(m_bqd);
      statcond_type     statcond(m_bqd);
      assembler_type    assembler(m_msh, m_bqd);

      hho_stab_type  stab_HHO(m_bqd);
      hdg_stab_type  stab_HDG(m_bqd);
      pikF_stab_type stab_PIKF(m_bqd);

      AssemblyInfo ai;

      timecounter tc, ttot;

      ttot.tic();
      size_t cell_i = 0;

      for (auto& cl : m_msh) {
         // Gradient Reconstruction
         matrix_dynamic GT;
         tc.tic();
         if (m_rp.m_precomputation) {
            GT = gradient_precomputed[cell_i];
         } else {
            gradrec_full.compute(m_msh, cl);
            GT = gradrec_full.oper;
         }
         tc.toc();
         ai.m_time_gradrec += tc.to_double();

         // Stabilisation
         tc.tic();
         if (m_rp.m_stab) {
            switch (m_rp.m_stab_type) {
               case PIKF: {
                  stab_PIKF.compute(m_msh, cl);

                  break;
               }
               case HHO: {
                  gradrec.compute(m_msh, cl);
                  stab_HHO.compute(m_msh, cl, gradrec.oper);

                  break;
               }
               case HDG: {
                  stab_HDG.compute(m_msh, cl);

                  break;
               }
               case NO: {
                  break;
               }
               default: throw std::invalid_argument("Unknown stabilization");
            }
         }
         tc.toc();
         ai.m_time_stab += tc.to_double();

         // Begin Assembly
         // Build rhs and lhs

         // Mechanical Computation

         tc.tic();
         elem.compute(m_msh, cl, lf, GT, m_solution_data.at(cell_i), m_rp.m_leray_param);

         dynamic_matrix<scalar_type> lhs = elem.K_int;
         dynamic_vector<scalar_type> rhs = elem.RTF;

         tc.toc();
         ai.m_time_elem += tc.to_double();
         ai.m_time_law += elem.time_law;

         // Stabilisation Contribution
         tc.tic();
         if (m_rp.m_stab) {
            switch (m_rp.m_stab_type) {
               case PIKF: {
                  assert(elem.K_int.rows() == stab_PIKF.data.rows());
                  assert(elem.K_int.cols() == stab_PIKF.data.cols());
                  assert(elem.RTF.rows() == (stab_PIKF.data * m_solution_data.at(cell_i)).rows());
                  assert(elem.RTF.cols() == (stab_PIKF.data * m_solution_data.at(cell_i)).cols());

                  lhs += m_rp.m_beta * stab_PIKF.data;
                  rhs -= m_rp.m_beta * stab_PIKF.data * m_solution_data.at(cell_i);
                  break;
               }
               case HHO: {
                  assert(elem.K_int.rows() == stab_HHO.data.rows());
                  assert(elem.K_int.cols() == stab_HHO.data.cols());
                  assert(elem.RTF.rows() == (stab_HHO.data * m_solution_data.at(cell_i)).rows());
                  assert(elem.RTF.cols() == (stab_HHO.data * m_solution_data.at(cell_i)).cols());

                  lhs += m_rp.m_beta * stab_HHO.data;
                  rhs -= m_rp.m_beta * stab_HHO.data * m_solution_data.at(cell_i);
                  break;
               }
               case HDG: {
                  assert(elem.K_int.rows() == stab_HDG.data.rows());
                  assert(elem.K_int.cols() == stab_HDG.data.cols());
                  assert(elem.RTF.rows() == (stab_HDG.data * m_solution_data.at(cell_i)).rows());
                  assert(elem.RTF.cols() == (stab_HDG.data * m_solution_data.at(cell_i)).cols());

                  lhs += m_rp.m_beta * stab_HDG.data;
                  rhs -= m_rp.m_beta * stab_HDG.data * m_solution_data.at(cell_i);
                  break;
               }
               case NO: {
                  break;
               }
               default: throw std::invalid_argument("Unknown stabilization");
            }
         }
         tc.toc();
         ai.m_time_stab += tc.to_double();

         // Static Condensation
         tc.tic();
         auto scnp     = statcond.compute(m_msh, cl, lhs, rhs, true);
         m_KTT[cell_i] = statcond.KTT;
         m_KTF[cell_i] = statcond.KTF;
         m_RT[cell_i]  = statcond.RT;
         tc.toc();
         ai.m_time_statcond += tc.to_double();

         assembler.assemble_nl(m_msh, cl, scnp, bcf, m_solution_faces);

         cell_i++;
      }

      assembler.finalize(m_system_matrix, m_system_rhs);

      ttot.toc();
      ai.m_time_assembly      = ttot.to_double();
      ai.m_linear_system_size = m_system_matrix.rows();
      return ai;
   }

   SolveInfo
   solve(void)
   {
      solver_type solver;

      timecounter tc;

      tc.tic();
      solver.analyzePattern(m_system_matrix);
      solver.factorize(m_system_matrix);

      if (solver.info() != Eigen::Success) {
         std::cerr << "ERROR: Could not factorize the matrix" << std::endl;
      }

      m_system_solution = solver.solve(m_system_rhs);
      if (solver.info() != Eigen::Success) {
         std::cerr << "ERROR: Could not solve the linear system" << std::endl;
      }
      tc.toc();

      return SolveInfo(m_system_matrix.rows(), m_system_matrix.nonZeros(), tc.to_double());
   }

   template<typename BoundaryConditionFunction>
   scalar_type
   postprocess(const BoundaryConditionFunction& bcf)
   {
      timecounter tc;
      tc.tic();

      assembler_type assembler(m_msh, m_bqd);

      const size_t fbs = howmany_dofs(m_bqd.face_basis);
      const size_t cbs = howmany_dofs(m_bqd.cell_basis);

      const auto solF =
        assembler.expand_solution_nl(m_msh, m_system_solution, bcf, m_solution_faces);

      // Update  unknowns
      // Update face Uf^{i+1} = Uf^i + delta Uf^i
      for (size_t i = 0; i < m_solution_faces.size(); i++) {
         assert(m_solution_faces.at(i).size() == fbs);
         m_solution_faces.at(i) += solF.segment(i * fbs, fbs);
      }
      // Update cell
      size_t cell_i = 0;
      for (auto& cl : m_msh) {
         // Extract the solution
         const auto fcs             = faces(m_msh, cl);
         const auto num_faces       = fcs.size();
         const auto total_faces_dof = num_faces * fbs;

         dynamic_vector<scalar_type> xFs = dynamic_vector<scalar_type>::Zero(total_faces_dof);

         for (size_t face_i = 0; face_i < num_faces; face_i++) {
            const auto fc  = fcs[face_i];
            const auto eid = find_element_id(m_msh.faces_begin(), m_msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

            const auto face_id = eid.second;

            xFs.segment(face_i * fbs, fbs) = solF.segment(face_id * fbs, fbs);
         }

         auto                              K_TT_ldlt = m_KTT[cell_i].llt();
         const dynamic_vector<scalar_type> xT =
           -K_TT_ldlt.solve(-m_RT[cell_i] + m_KTF[cell_i] * xFs);

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
   test_convergence(const size_t iter)
   {
      // norm of the solution
      scalar_type error_un = 0;
      for (size_t i = 0; i < m_solution_faces.size(); i++) {
         scalar_type norm = m_solution_faces[i].norm();
         error_un += norm * norm;
      }

      error_un = std::sqrt(error_un);

      if (error_un <= scalar_type(10E-15)) {
         error_un = scalar_type(10E16);
      }

      // norm of the rhs
      const scalar_type residual = m_system_rhs.norm();
      scalar_type       max_error(0.0);
      for (size_t i = 0; i < m_system_rhs.size(); i++)
         max_error = std::max(max_error, std::abs(m_system_rhs(i)));

      // norm of the increment
      const scalar_type error_incr = m_system_solution.norm();

      scalar_type relative_displ = error_incr / error_un;
      scalar_type relative_error = residual / error_un;

      if (iter == 0) {
         relative_displ = 1;
         relative_error = 1;
      }

      const size_t nb_faces_dof = m_bqd.face_basis.size() * m_msh.faces_size();

      if (m_verbose) {
         std::string s_iter = "   " + std::to_string(iter) + "               ";
         s_iter.resize(9);

         if (iter == 0) {
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
         std::cout << "| " << s_iter << " |   " << error_incr << " |   " << relative_displ
                   << " |   " << residual << " |   " << relative_error << "  |  " << max_error
                   << "  |" << std::endl;
         std::cout << "-------------------------------------------------------------------------"
                      "---------------------"
                   << std::endl;
         std::cout.flags(f);
      }

      const scalar_type error = std::max(relative_displ, relative_error);

      if (error <= m_rp.m_epsilon) {
         return true;
      } else {
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
