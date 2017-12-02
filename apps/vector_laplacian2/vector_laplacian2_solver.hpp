/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016, 2017 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#include <iostream>

#include <sstream>

#include "hho/assembler.hpp"
#include "hho/gradient_reconstruction.hpp"
#include "hho/hho_bq.hpp"
#include "hho/hho_utils.hpp"
#include "hho/projector.hpp"
#include "hho/static_condensation.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

struct assembly_info
{
   size_t linear_system_size;
   double time_gradrec, time_statcond, time_stab, time_assembly;
};

struct solver_info
{
   double time_solver;
};

struct postprocess_info
{
   double time_postprocess;
};

struct LaplacianParameters
{
   LaplacianParameters(double _lambda = 1.)
     : lambda(_lambda)
   {
      // Do nothing
   }

   double lambda;
};

template<typename Mesh>
class vector_laplacian_solver
{
   typedef Mesh                            mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef typename mesh_type::face        face_type;

   typedef dynamic_matrix<scalar_type> matrix_dynamic;
   typedef dynamic_vector<scalar_type> vector_dynamic;

   typedef disk::hho::basis_quadrature_data_full<mesh_type,
                                                 disk::scaled_monomial_vector_basis,
                                                 disk::Raviart_Thomas_matrix_basis,
                                                 disk::quadrature>
     bqdata_type;

   typedef disk::hho::gradient_reconstruction_bq<bqdata_type> gradrec_type;

   typedef disk::hho::static_condensation_bq<bqdata_type> statcond_type;
   typedef disk::hho::assembler_bq<bqdata_type>           assembler_type;
   typedef disk::hho::projector_bq<bqdata_type>           projector_type;

   typename assembler_type::sparse_matrix_type m_system_matrix;
   typename assembler_type::vector_type        m_system_rhs, m_system_solution;

   size_t m_cell_degree, m_face_degree;

   bqdata_type m_bqd;

   const mesh_type& m_msh;

   size_t m_dim;

   std::vector<vector_dynamic> m_solution_data;

   bool m_verbose;

   LaplacianParameters m_laplacian_parameters;

 public:
   vector_laplacian_solver(const mesh_type& msh, size_t degree, int l = 0)
     : m_msh(msh)
     , m_verbose(false)
   {
      if (l < -1 or l > 1) {
         std::cout << "'l' should be -1, 0 or 1. Reverting to 0." << std::endl;
         l = 0;
      }

      if (degree == 0 && l == -1) {
         std::cout << "'l' should be 0 or 1. Reverting to 0." << std::endl;
         l = 0;
      }

      m_cell_degree = degree + l;
      m_face_degree = degree;

      m_dim = msh.dimension;

      m_laplacian_parameters.lambda = 1.0;
      m_bqd                         = bqdata_type(m_face_degree, m_cell_degree, m_cell_degree + 1);
   }

   vector_laplacian_solver(const mesh_type&          msh,
                           size_t                    degree,
                           const LaplacianParameters data,
                           int                       l = 0)
     : m_msh(msh)
     , m_verbose(false)
   {
      if (l < -1 or l > 1) {
         std::cout << "'l' should be -1, 0 or 1. Reverting to 0." << std::endl;
         l = 0;
      }

      if (degree == 0 && l == -1) {
         std::cout << "'l' should be 0 or 1. Reverting to 0." << std::endl;
         l = 0;
      }

      m_cell_degree = degree + l;
      m_face_degree = degree;

      m_laplacian_parameters.lambda = data.lambda;
      m_bqd                         = bqdata_type(m_cell_degree, m_face_degree, m_cell_degree + 1);
   }

   size_t getDofs() { return m_msh.faces_size() * m_bqd.face_basis.size(); }

   void changeLaplacianParameters(const LaplacianParameters data)
   {
      m_laplacian_parameters.lambda = data.lambda;
   }

   bool verbose(void) const { return m_verbose; }
   void verbose(bool v) { m_verbose = v; }

   template<typename LoadFunction, typename BoundaryConditionFunction>
   assembly_info assemble(const LoadFunction& lf, const BoundaryConditionFunction& bcf)
   {
      auto gradrec   = gradrec_type(m_bqd);
      auto statcond  = statcond_type(m_bqd);
      auto assembler = assembler_type(m_msh, m_bqd);

      assembly_info ai;
      bzero(&ai, sizeof(ai));

      timecounter tc;

      for (auto& cl : m_msh) {
         tc.tic();
         gradrec.compute(m_msh, cl);
         tc.toc();
         ai.time_gradrec += tc.to_double();

         tc.tic();
         const auto cell_rhs = disk::hho::compute_rhs(m_msh, cl, lf, m_bqd);

         const dynamic_matrix<scalar_type> loc  = m_laplacian_parameters.lambda * gradrec.data;
         const auto                        scnp = statcond.compute(m_msh, cl, loc, cell_rhs);
         tc.toc();
         ai.time_statcond += tc.to_double();

         assembler.assemble(m_msh, cl, scnp);
      }

      assembler.impose_boundary_conditions(m_msh, bcf);
      assembler.finalize(m_system_matrix, m_system_rhs);

      ai.linear_system_size = m_system_matrix.rows();
      ai.time_assembly      = ai.time_gradrec + ai.time_stab + ai.time_statcond;
      return ai;
   }

   solver_info solve(void)
   {
#ifdef HAVE_INTEL_MKL
      Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>> solver;
#else
      Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>> solver;
#endif

      solver_info si;

      const size_t systsz = m_system_matrix.rows();
      const size_t nnz    = m_system_matrix.nonZeros();

      if (verbose()) {
         std::cout << "Starting linear solver..." << std::endl;
         std::cout << " * Solving for " << systsz << " unknowns." << std::endl;
         std::cout << " * Matrix fill: " << 100.0 * double(nnz) / (systsz * systsz) << "%"
                   << std::endl;
      }

      timecounter tc;

      tc.tic();
      solver.analyzePattern(m_system_matrix);
      solver.factorize(m_system_matrix);
      m_system_solution = solver.solve(m_system_rhs);
      tc.toc();
      si.time_solver = tc.to_double();

      return si;
   }

   template<typename LoadFunction>
   postprocess_info postprocess(const LoadFunction& lf)
   {
      auto gradrec  = gradrec_type(m_bqd);
      auto statcond = statcond_type(m_bqd);

      const size_t fbs = m_bqd.face_basis.size();

      postprocess_info pi;

      m_solution_data.reserve(m_msh.cells_size());

      timecounter tc;
      tc.tic();
      for (auto& cl : m_msh) {
         const auto fcs       = faces(m_msh, cl);
         const auto num_faces = fcs.size();

         dynamic_vector<scalar_type> xFs = dynamic_vector<scalar_type>::Zero(num_faces * fbs);

         for (size_t face_i = 0; face_i < num_faces; face_i++) {
            const auto fc  = fcs[face_i];
            const auto eid = find_element_id(m_msh.faces_begin(), m_msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

            const auto face_id = eid.second;

            dynamic_vector<scalar_type> xF     = dynamic_vector<scalar_type>::Zero(fbs);
            xF                                 = m_system_solution.block(face_id * fbs, 0, fbs, 1);
            xFs.block(face_i * fbs, 0, fbs, 1) = xF;
         }

         gradrec.compute(m_msh, cl);
         const dynamic_matrix<scalar_type> loc      = m_laplacian_parameters.lambda * gradrec.data;
         const auto                        cell_rhs = disk::hho::compute_rhs(m_msh, cl, lf, m_bqd);
         const dynamic_vector<scalar_type> x = statcond.recover(m_msh, cl, loc, cell_rhs, xFs);
         m_solution_data.push_back(x);
      }
      tc.toc();

      pi.time_postprocess = tc.to_double();

      return pi;
   }

   template<typename AnalyticalSolution>
   scalar_type compute_l2_error(const AnalyticalSolution& as)
   {
      scalar_type err_dof = scalar_type{0.0};

      projector_type projk(m_bqd);

      size_t i = 0;

      for (auto& cl : m_msh) {
         const auto                        x        = m_solution_data.at(i++);
         const dynamic_vector<scalar_type> true_dof = projk.projectOnCell(m_msh, cl, as);
         const dynamic_vector<scalar_type> comp_dof = x.block(0, 0, true_dof.size(), 1);
         const dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
         err_dof += diff_dof.dot(projk.cell_mm * diff_dof);
      }

      return sqrt(err_dof);
   }

   template<typename AnalyticalSolution>
   scalar_type compute_l2_gradient_error(const AnalyticalSolution& grad)
   {
      scalar_type err_dof = scalar_type{0.0};

      projector_type projk(m_bqd);

      gradrec_type gradrec(m_bqd);

      size_t i = 0;

      for (auto& cl : m_msh) {
         const auto x = m_solution_data.at(i++);
         gradrec.compute(m_msh, cl);
         const dynamic_vector<scalar_type> RTu = gradrec.oper * x;

         const dynamic_vector<scalar_type> true_dof = projk.projectGradOnCell(m_msh, cl, grad);
         const dynamic_vector<scalar_type> comp_dof = RTu.block(0, 0, true_dof.size(), 1);
         const dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
         err_dof += diff_dof.dot(projk.grad_mm * diff_dof);
      }

      return sqrt(err_dof);
   }
};
