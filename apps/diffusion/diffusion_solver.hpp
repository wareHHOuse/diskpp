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

#include "config.h"

#include "hho/assembler.hpp"
#include "hho/gradient_reconstruction.hpp"
#include "hho/hho.hpp"
#include "hho/hho_bq.hpp"
#include "hho/hho_utils.hpp"
#include "hho/projector.hpp"
#include "hho/stabilization.hpp"
#include "hho/static_condensation.hpp"

#include <Eigen/Eigenvalues>

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

template<typename T>
bool
conjugated_gradient(const Eigen::SparseMatrix<T>&              A,
                    const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
                    Eigen::Matrix<T, Eigen::Dynamic, 1>&       x)
{
   size_t N    = A.cols();
   size_t iter = 0;
   T      nr, nr0;
   T      alpha, beta, rho;

   Eigen::Matrix<T, Eigen::Dynamic, 1> d(N), r(N), r0(N), y(N);
   x = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(N);

   r0 = d = r = b - A * x;
   nr = nr0 = r.norm();

   std::ofstream ofs("cocg_nopre_convergence.txt");

   while (nr / nr0 > 1e-8 && iter < 40000 && nr / nr0 < 10000) {
      std::cout << "                                                 \r";
      std::cout << " -> Iteration " << iter << ", rr = ";
      std::cout << nr / nr0 << "\b\r";
      std::cout.flush();

      ofs << nr / nr0 << std::endl;
      y     = A * d;
      rho   = r.dot(r);
      alpha = rho / d.dot(y);
      x     = x + alpha * d;
      r     = r - alpha * y;
      beta  = r.dot(r) / rho;
      d     = r + beta * d;

      nr = r.norm();
      iter++;
   }

   ofs << nr / nr0 << std::endl;
   ofs.close();

   std::cout << " -> Iteration " << iter << ", rr = " << nr / nr0 << std::endl;

   return true;
}

struct assembly_info
{
   size_t linear_system_size;
   double time_gradrec, time_statcond, time_stab;
};

struct solver_info
{
   double time_solver;
};

struct postprocess_info
{
   double time_postprocess;
};

template<typename Mesh>
class diffusion_solver
{
   typedef Mesh                            mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef typename mesh_type::face        face_type;

   typedef disk::hho::
     basis_quadrature_data<mesh_type, disk::scaled_monomial_scalar_basis, disk::quadrature>
       bqdata_type;

   typedef disk::hho::gradient_reconstruction_bq<bqdata_type>                     gradrec_type;
   typedef disk::hho::hho_stabilization_bq<bqdata_type>                           stab_type;
   typedef disk::hho::static_condensation_bq<bqdata_type>                         statcond_type;
   typedef disk::hho::assembler_bq<bqdata_type>                                   assembler_type;
   typedef static_matrix<scalar_type, mesh_type::dimension, mesh_type::dimension> tensor_type;

   size_t m_cell_degree, m_face_degree;

   bqdata_type m_bqd;

   typename assembler_type::sparse_matrix_type m_system_matrix;
   typename assembler_type::vector_type        m_system_rhs, m_system_solution;

   const mesh_type& m_msh;

   std::vector<dynamic_vector<scalar_type>> m_postprocess_data;

   bool m_verbose;

 public:
   diffusion_solver(const mesh_type& msh, size_t degree, int l = 0)
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

      // std::cout << "HHO(" << m_cell_degree << "," << m_face_degree << ")" << std::endl;

      m_bqd = bqdata_type(m_cell_degree, m_face_degree);
   }

   bool verbose(void) const { return m_verbose; }
   void verbose(bool v) { m_verbose = v; }

   template<typename LoadFunction, typename BoundaryConditionFunction>
   assembly_info assemble(const LoadFunction& lf, const BoundaryConditionFunction& bcf)
   {
      auto gradrec  = gradrec_type(m_bqd);
      auto stab     = stab_type(m_bqd);
      auto statcond = statcond_type(m_bqd);
      // auto assembler  = assembler_type(m_msh, m_face_degree);
      auto assembler = assembler_type(m_msh, m_bqd);

      assembly_info ai;
      bzero(&ai, sizeof(ai));

      timecounter tc;

      auto tf = [](const typename mesh_type::point_type& pt) -> tensor_type {
         tensor_type ret = tensor_type::Identity();
         return ret;
         // return 6.72071 * ret;
         auto c = cos(M_PI * pt.x() / 0.004);
         auto s = sin(M_PI * pt.y() / 0.004);
         return ret * (1 + 100 * c * c * s * s);
      };

      size_t elem_i = 0;
      for (auto& cl : m_msh) {

         // if (elem_i%10000 == 0)
         //    std::cout << elem_i << std::endl;

         elem_i++;

         tc.tic();
         gradrec.compute(m_msh, cl, tf);
         tc.toc();
         ai.time_gradrec += tc.to_double();

         tc.tic();
         stab.compute(m_msh, cl, gradrec.oper);
         tc.toc();
         ai.time_stab += tc.to_double();

         tc.tic();
         const auto cell_rhs = disk::hho::compute_rhs(m_msh, cl, lf, m_bqd, m_bqd.cell_degree());
         dynamic_matrix<scalar_type> loc  = gradrec.data + stab.data;
         const auto                  scnp = statcond.compute(m_msh, cl, loc, cell_rhs);
         tc.toc();
         ai.time_statcond += tc.to_double();

         assembler.assemble(m_msh, cl, scnp);
      }

      assembler.impose_boundary_conditions(m_msh, bcf);
      assembler.finalize(m_system_matrix, m_system_rhs);

      ai.linear_system_size = m_system_matrix.rows();
      return ai;
   }

   solver_info solve(void)
   {
#ifdef HAVE_INTEL_MKL
      Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>> solver;
      // solver.pardisoParameterArray()[59] = 0; //out-of-core
#else
      Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>> solver;
#endif

      solver_info si;

      size_t systsz = m_system_matrix.rows();
      size_t nnz    = m_system_matrix.nonZeros();

      if (verbose()) {
         std::cout << "Starting linear solver..." << std::endl;
         std::cout << " * Solving for " << systsz << " unknowns." << std::endl;
         std::cout << " * Matrix fill: " << 100.0 * double(nnz) / (systsz * systsz) << "%"
                   << std::endl;
      }

      timecounter_new tc;
      tc.tic();

      solver.analyzePattern(m_system_matrix);
      solver.factorize(m_system_matrix);
      m_system_solution = solver.solve(m_system_rhs);

      tc.toc();
      si.time_solver = tc.to_double();

      return si;
   }

   template<typename LoadFunction, typename BoundaryConditionFunction>
   postprocess_info postprocess(const LoadFunction& lf, const BoundaryConditionFunction& bcf)
   {
      // expand solution
      auto assembler = assembler_type(m_msh, m_bqd);

      auto gradrec  = gradrec_type(m_bqd);
      auto stab     = stab_type(m_bqd);
      auto statcond = statcond_type(m_bqd);

      const size_t fbs = howmany_dofs(m_bqd.face_basis);

      postprocess_info pi;

      m_postprocess_data.reserve(m_msh.cells_size());

      auto tf = [](const typename mesh_type::point_type& pt) -> tensor_type {
         tensor_type ret = tensor_type::Identity();
         return ret;
         // return 6.72071 * ret;
         auto c = cos(M_PI * pt.x() / 0.004);
         auto s = sin(M_PI * pt.y() / 0.004);
         return ret * (1 + 100 * c * c * s * s);
      };

      timecounter tc;
      tc.tic();

#define DUMP_SOLUTION_DATA

#ifdef DUMP_SOLUTION_DATA
      auto          elem_num = m_msh.cells_size();
      auto          cb_deg   = m_bqd.cell_basis.degree();
      auto          fb_deg   = m_bqd.face_basis.degree();
      std::ofstream sol_data("solution.bin", std::ios::binary);
      sol_data.write(reinterpret_cast<char*>(&elem_num), sizeof(elem_num));
      sol_data.write(reinterpret_cast<char*>(&cb_deg), sizeof(elem_num));
      sol_data.write(reinterpret_cast<char*>(&fb_deg), sizeof(elem_num));
#endif
      for (auto& cl : m_msh) {
         auto   fcs       = faces(m_msh, cl);
         size_t num_faces = fcs.size();

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

         gradrec.compute(m_msh, cl, tf);
         stab.compute(m_msh, cl, gradrec.oper);
         const dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
         const auto cell_rhs = disk::hho::compute_rhs(m_msh, cl, lf, m_bqd, m_bqd.cell_degree());
         const dynamic_vector<scalar_type> x = statcond.recover(m_msh, cl, loc, cell_rhs, xFs);
         m_postprocess_data.push_back(x);
         dynamic_vector<scalar_type> Rx = gradrec.oper * x;
#ifdef DUMP_SOLUTION_DATA
         // dump number of faces
         // sol_data.write( reinterpret_cast<char *>(&num_faces), sizeof(size_t) );

         // dump vT and vF coefficients
         // for (size_t i = 0; i < x.size(); i++)
         //    sol_data.write( reinterpret_cast<char *>(&x(i)), sizeof(x(i)) );

         // dump R(v)
         dynamic_vector<scalar_type> Rxd;
         Rxd                 = dynamic_vector<scalar_type>::Zero(Rx.size() + 1);
         Rxd(0)              = x(0);
         Rxd.tail(Rx.size()) = Rx;

         for (size_t i = 0; i < Rxd.size(); i++)
            sol_data.write(reinterpret_cast<char*>(&Rxd(i)), sizeof(Rxd(i)));
#endif
      }
#ifdef DUMP_SOLUTION_DATA
      sol_data.close();
#endif
      tc.toc();

      pi.time_postprocess = tc.to_double();

      return pi;
   }

   template<typename AnalyticalSolution>
   scalar_type compute_l2_error(const AnalyticalSolution& as)
   {
      scalar_type err_dof = 0.0;

      disk::hho::projector_bq<bqdata_type> projk(m_bqd);

      size_t i = 0;
      for (auto& cl : m_msh) {
         auto                        x        = m_postprocess_data.at(i++);
         dynamic_vector<scalar_type> true_dof = projk.projectOnCell(m_msh, cl, as);
         dynamic_vector<scalar_type> comp_dof = x.block(0, 0, true_dof.size(), 1);
         dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
         err_dof += diff_dof.dot(projk.cell_mm * diff_dof);
      }

      return sqrt(err_dof);
   }

   void plot_solution(const std::string& filename)
   {
      std::ofstream ofs(filename);

      size_t cell_i = 0;
      for (auto& cl : m_msh) {
         auto x = m_postprocess_data.at(cell_i++);
         // auto qps = m_bqd.cell_quadrature.integrate(m_msh, cl);
         // for (auto& qp : qps)
         auto tps = make_test_points(m_msh, cl, 10);
         for (auto& tp : tps) {
            auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, tp);

            scalar_type pot = 0.0;
            for (size_t i = 0; i < howmany_dofs(m_bqd.cell_basis); i++)
               pot += phi[i] * x(i);

            // auto tp = qp.point();
            for (size_t i = 0; i < mesh_type::dimension; i++)
               ofs << tp[i] << " ";
            ofs << pot << std::endl;
         }
      }

      ofs.close();
   }
};
