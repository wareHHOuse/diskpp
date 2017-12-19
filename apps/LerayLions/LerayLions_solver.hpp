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

#include <iostream>
#include <sstream>

#include <list>
#include <vector>

#include "hho/gradient_reconstruction.hpp"
#include "hho/hho.hpp"
#include "hho/hho_bq.hpp"
#include "hho/projector.hpp"

#include "Informations.hpp"
#include "NewtonSolver/newton_solver.hpp"
#include "Parameters.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

struct time_step
{
   double time;
   size_t level;
};

template<template<typename, size_t, typename> class Mesh, typename T, size_t DIM, typename Storage>
class leraylions_solver
{
   typedef Mesh<T, DIM, Storage>           mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef ParamRun<T>                     param_type;

   typedef disk::hho::basis_quadrature_data_full<mesh_type,
                                                 disk::scaled_monomial_scalar_basis,
                                                 disk::scaled_monomial_vector_basis,
                                                 disk::quadrature>
     bqdata_type;

   typedef dynamic_matrix<scalar_type> matrix_dynamic;
   typedef dynamic_vector<scalar_type> vector_dynamic;

   typedef disk::hho::gradient_reconstruction_full_bq<bqdata_type> gradrec_full_type;
   typedef disk::hho::projector_bq<bqdata_type>                    projector_type;

   bqdata_type m_bqd;

   const mesh_type& m_msh;

   std::vector<vector_dynamic> m_solution_data;
   std::vector<vector_dynamic> m_solution_cells, m_solution_faces;
   std::vector<matrix_dynamic> m_gradient_precomputed;

   bool m_verbose, m_convergence;

   param_type m_rp;

   const static size_t dimension = mesh_type::dimension;

   size_t total_dof_depl_static;

   void
   init()
   {
      m_solution_data.clear();
      m_solution_cells.clear();
      m_solution_faces.clear();

      m_solution_data.reserve(m_msh.cells_size());
      m_solution_cells.reserve(m_msh.cells_size());
      m_solution_faces.reserve(m_msh.faces_size());

      const size_t num_cell_dofs = howmany_dofs(m_bqd.cell_basis);
      const size_t num_face_dofs = howmany_dofs(m_bqd.face_basis);
      const size_t total_dof =
        m_msh.cells_size() * num_cell_dofs + m_msh.faces_size() * num_face_dofs;

      for (auto& cl : m_msh) {
         auto         fcs       = faces(m_msh, cl);
         const size_t num_faces = fcs.size();
         m_solution_data.push_back(vector_dynamic::Zero(num_cell_dofs + num_faces * num_face_dofs));
         m_solution_cells.push_back(vector_dynamic::Zero(num_cell_dofs));
      }

      for (size_t i = 0; i < m_msh.faces_size(); i++) {
         m_solution_faces.push_back(vector_dynamic::Zero(num_face_dofs));
      }

      if (m_verbose) {
         std::cout << "** Numbers of cells: " << m_msh.cells_size() << std::endl;
         std::cout << "** Numbers of faces: " << m_msh.faces_size()
                   << " ( boundary faces: " << m_msh.boundary_faces_size() << " )" << std::endl;
         std::cout << "** Numbers of dofs: " << total_dof << std::endl;
         std::cout << "** After static condensation: " << std::endl;
         std::cout << "** Numbers of dofs: " << m_msh.faces_size() * num_face_dofs << std::endl;
      }

      total_dof_depl_static = m_msh.faces_size() * num_face_dofs;

      for (size_t i = 0; i < m_msh.cells_size(); i++) {
         m_solution_data[i].setConstant(1);
         m_solution_cells[i].setConstant(1);
      }

      for (size_t i = 0; i < m_msh.faces_size(); i++) {
         m_solution_faces[i].setConstant(1);
      }
   }

   void
   pre_computation()
   {
      gradrec_full_type gradrec_full(m_bqd);

      m_gradient_precomputed.clear();
      m_gradient_precomputed.reserve(m_msh.cells_size());

      for (auto& cl : m_msh) {
         /////// Gradient Reconstruction /////////
         gradrec_full.compute(m_msh, cl);
         m_gradient_precomputed.push_back(gradrec_full.oper);
      }
   }

   // compute solution for p = 2
   template<typename LoadFunction, typename BoundaryConditionFunction>
   void
   prediction(const LoadFunction& lf, const BoundaryConditionFunction& bcf, const scalar_type time)
   {
      param_type rp(m_rp);
      rp.m_leray_param = 2;

      // Newton solver
      NewtonRaphson_solver_leraylions<bqdata_type> newton_solver(m_msh, m_bqd, rp);

      newton_solver.initialize(m_solution_cells, m_solution_faces, m_solution_data);

      auto rlf = [&lf, &time ](const point<scalar_type, DIM>& p) -> auto
      {
         return disk::mm_prod(time, lf(p));
      };

      auto rbcf = [&bcf, &time ](const point<scalar_type, DIM>& p) -> auto
      {
         return disk::mm_prod(time, bcf(p));
      };

      if (m_verbose) {
         std::cout << "-------------------------------------------------------------------------"
                      "---------------------"
                   << std::endl;
         std::cout << "**************************************** Prediction "
                      "*****************************************|"
                   << std::endl;
      }

      NewtonSolverInfo newton_info = newton_solver.compute(rlf, rbcf, m_gradient_precomputed);

      newton_solver.save_solutions(m_solution_cells, m_solution_faces, m_solution_data);
   }

 public:
   leraylions_solver(const mesh_type& msh, const param_type& rp) :
     m_msh(msh), m_verbose(rp.m_verbose), m_convergence(false), m_rp(rp)
   {
      size_t face_degree = rp.m_face_degree;
      if (face_degree <= 0) {
         std::cout << "'face_degree' should be > 0. Reverting to 1." << std::endl;
         face_degree = 1;
      }

      m_rp.m_face_degree = face_degree;

      size_t cell_degree = rp.m_cell_degree;
      if (face_degree - 1 > cell_degree or cell_degree > face_degree + 1) {
         std::cout << "'cell_degree' should be 'face_degree + 1' => 'cell_degree' => 'face_degree "
                      "-1'. Reverting to 'face_degree'."
                   << std::endl;
         cell_degree = face_degree;
      }

      m_rp.m_cell_degree = cell_degree;

      size_t grad_degree = rp.m_grad_degree;
      if (grad_degree < cell_degree) {
         std::cout << "'grad_degree' should be > 'cell_degree'. Reverting to 'cell_degree'."
                   << std::endl;
         grad_degree = cell_degree;
      }

      m_rp.m_grad_degree = grad_degree;

      m_bqd = bqdata_type(face_degree, cell_degree, grad_degree);

      if (m_verbose) {
         m_bqd.info_degree();
         m_rp.infos();
      }

      init();
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

   template<typename LoadFunction, typename BoundaryConditionFunction>
   SolverInfo
   compute(const LoadFunction& lf, const BoundaryConditionFunction& bcf)
   {

      SolverInfo  si;
      timecounter ttot;
      ttot.tic();

      // time step
      std::list<time_step> list_step;

      scalar_type time1 = 0.0;
      scalar_type time2 = 0.0;
      for (size_t n = 0; n < m_rp.m_time_step.size(); n++) {
         auto time_info            = m_rp.m_time_step[n];
         time2                     = time_info.first;
         const scalar_type delta_t = (time2 - time1) / time_info.second;
         for (size_t i = 0; i < time_info.second; i++) {
            time_step step;
            step.time  = time1 + (i + 1) * delta_t;
            step.level = 1;
            list_step.push_back(step);
         }
         time1 = time2;
      }

      size_t current_step = 0;
      size_t total_step   = list_step.size();

      scalar_type old_time = 0.0;

      // time of saving
      bool time_saving(false);
      if (m_rp.m_n_time_save > 0) {
         time_saving = true;
      }

      // Precomputation
      if (m_rp.m_precomputation) {
         timecounter t1;
         t1.tic();
         this->pre_computation();
         t1.toc();
         if (m_verbose) std::cout << "-Precomputation: " << t1.to_double() << " sec" << std::endl;
      }

      // prediction (computation for p=2)
      if (m_rp.m_pred) {
         prediction(lf, bcf, list_step.front().time);
      }

      // Newton solver
      NewtonRaphson_solver_leraylions<bqdata_type> newton_solver(m_msh, m_bqd, m_rp);

      newton_solver.initialize(m_solution_cells, m_solution_faces, m_solution_data);

      // loading
      while (!list_step.empty()) {
         current_step += 1;
         time_step         step         = list_step.front();
         const scalar_type current_time = step.time;
         if (m_verbose) {
            std::cout << "------------------------------------------------------------------------"
                         "----------------------"
                      << std::endl;
            std::cout << "************************** Time : " << current_time
                      << " sec (step: " << current_step << "/" << total_step
                      << ", sublevel: " << step.level << " ) ***************************|"
                      << std::endl;
         }

         auto rlf = [&lf, &current_time ](const point<scalar_type, DIM>& p) -> auto
         {
            return disk::mm_prod(current_time, lf(p));
         };

         auto rbcf = [&bcf, &current_time ](const point<scalar_type, DIM>& p) -> auto
         {
            return disk::mm_prod(current_time, bcf(p));
         };

         // correction
         NewtonSolverInfo newton_info = newton_solver.compute(rlf, rbcf, m_gradient_precomputed);
         si.updateInfo(newton_info);

         if (m_verbose) {
            std::cout << "** Time in this step " << newton_info.m_time_newton << " sec"
                      << std::endl;
            std::cout << "**** Assembly time: " << newton_info.m_assembly_info.m_time_assembly
                      << " sec" << std::endl;
            std::cout << "****** Gradient reconstruction: "
                      << newton_info.m_assembly_info.m_time_gradrec << " sec" << std::endl;
            std::cout << "****** Stabilisation: " << newton_info.m_assembly_info.m_time_stab
                      << " sec" << std::endl;
            std::cout << "****** Mechanical computation: "
                      << newton_info.m_assembly_info.m_time_elem << " sec" << std::endl;
            std::cout << "       *** Behavior computation: "
                      << newton_info.m_assembly_info.m_time_law << " sec" << std::endl;
            std::cout << "****** Static condensation: "
                      << newton_info.m_assembly_info.m_time_statcond << " sec" << std::endl;
            std::cout << "****** Postprocess time: " << newton_info.m_assembly_info.m_time_postpro
                      << " sec" << std::endl;
            std::cout << "**** Solver time: " << newton_info.m_solve_info.m_time_solve << " sec"
                      << std::endl;
         }

         m_convergence = newton_solver.test_convergence();

         if (!m_convergence) {
            if (step.level > m_rp.m_sublevel) {
               std::cout << "***********************************************************"
                         << std::endl;
               std::cout << "***** PROBLEM OF CONVERGENCE: We stop the calcul here *****"
                         << std::endl;
               std::cout << "***********************************************************"
                         << std::endl;
               break;
            } else {
               if (m_verbose) {
                  std::cout << "***********************************************************"
                            << std::endl;
                  std::cout << "*****     NO CONVERGENCE: We split the time step     ******"
                            << std::endl;
                  std::cout << "***********************************************************"
                            << std::endl;
               }
               total_step += 1;
               current_step -= 1;
               time_step new_step;
               new_step.time  = old_time + (current_time - old_time) / 2.0;
               new_step.level = step.level + 1;
               list_step.push_front(new_step);
            }
         } else {
            old_time = current_time;
            list_step.pop_front();
            newton_solver.save_solutions(m_solution_cells, m_solution_faces, m_solution_data);

            if (time_saving) {
               if (m_rp.m_time_save.front() < old_time + 1E-5) {
                  m_rp.m_time_save.pop_front();
                  if (m_rp.m_time_save.empty()) time_saving = false;
               }
            }
         }
      }

      si.m_time_step = total_step;

      ttot.toc();
      si.m_time_solver = ttot.to_double();
      return si;
   }

   bool
   test_convergence() const
   {
      return m_convergence;
   }

   template<typename AnalyticalSolution>
   scalar_type
   compute_l2_error(const AnalyticalSolution& as)
   {
      scalar_type err_dof = scalar_type{0.0};

      projector_type projk(m_bqd);

      size_t i = 0;
      for (auto& cl : m_msh) {
         const auto                        x        = m_solution_cells.at(i++);
         const dynamic_vector<scalar_type> true_dof = projk.projectOnCell(m_msh, cl, as);
         const dynamic_vector<scalar_type> comp_dof = x.block(0, 0, true_dof.size(), 1);
         const dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
         err_dof += diff_dof.dot(projk.cell_mm * diff_dof);
      }

      return sqrt(err_dof);
   }

   size_t
   getDofs()
   {
      return total_dof_depl_static;
   }

   template<typename AnalyticalSolution>
   scalar_type
   compute_l2_gradient_error(const AnalyticalSolution& grad)
   {
      scalar_type err_dof = scalar_type{0.0};

      projector_type projk(m_bqd);

      gradrec_full_type gradrec_full(m_bqd);

      size_t i = 0;
      for (auto& cl : m_msh) {
         const auto x = m_solution_data.at(i++);
         gradrec_full.compute(m_msh, cl);
         const dynamic_vector<scalar_type> GTu      = gradrec_full.oper * x;
         const dynamic_vector<scalar_type> true_dof = projk.projectGradOnCell(m_msh, cl, grad);
         const dynamic_vector<scalar_type> comp_dof = GTu.block(0, 0, true_dof.size(), 1);
         const dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
         err_dof += diff_dof.dot(projk.grad_mm * diff_dof);
      }

      return sqrt(err_dof);
   }
};
