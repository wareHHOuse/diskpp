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

/* Solve the linear elasticity problem with HHO method
 * A hybrid high-order locking-free method for linear elasticity on general meshes
 * Di Pietro, Daniele A. and Ern, Alexandre.
 * Comput. Methods Appl. Mech. Engrg. (2015)
 *  DOI: 10.1016/j.cma.2014.09.009
 */

#include <iostream>
#include <sstream>

#include "hho/assembler.hpp"
#include "hho/divergence.hpp"
#include "hho/hho_bq.hpp"
#include "hho/hho_utils.hpp"
#include "hho/projector.hpp"
#include "hho/stabilization.hpp"
#include "hho/static_condensation.hpp"
#include "hho/sym_gradient_reconstruction.hpp"
#include "mechanics/BoundaryConditions.hpp"

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

struct assembly_info
{
   size_t linear_system_size;
   double time_gradrec, time_statcond, time_stab, time_assembly, time_divrec;
};

struct solver_info
{
   double time_solver;
};

struct postprocess_info
{
   double time_postprocess;
};

struct ElasticityParameters
{
   double lambda;
   double mu;
};

template<typename Mesh>
class linear_elasticity_solver
{
   typedef Mesh                                           mesh_type;
   typedef typename mesh_type::scalar_type                scalar_type;
   typedef typename mesh_type::cell                       cell_type;
   typedef typename mesh_type::face                       face_type;
   typedef disk::mechanics::BoundaryConditions<mesh_type> bnd_type;

   typedef disk::hho::basis_quadrature_data_linear_elasticity<mesh_type,
                                                              disk::scaled_monomial_vector_basis,
                                                              disk::scaled_monomial_scalar_basis,
                                                              disk::quadrature>
     bqdata_type;

   typedef dynamic_matrix<scalar_type> matrix_dynamic;
   typedef dynamic_vector<scalar_type> vector_dynamic;

   typedef disk::hho::sym_gradient_reconstruction_bq<bqdata_type> sgradrec_type;
   typedef disk::hho::divergence_reconstruction_bq<bqdata_type>   divrec_type;

   typedef disk::hho::hho_stabilization_bq<bqdata_type> stab_type;

   typedef disk::hho::static_condensation_bq<bqdata_type> statcond_type;

   typedef disk::hho::assembler_by_elimination_mechanics_bq<bqdata_type, bnd_type> assembler_type;

   typedef disk::hho::projector_bq<bqdata_type> projector_type;

   typename assembler_type::sparse_matrix_type m_system_matrix;
   typename assembler_type::vector_type        m_system_rhs, m_system_solution;

   size_t m_cell_degree, m_face_degree;

   const static size_t dimension = mesh_type::dimension;

   const bnd_type&  m_bnd;
   const mesh_type& m_msh;
   bqdata_type      m_bqd;

   std::vector<vector_dynamic> m_solution_data;

   bool m_verbose;

   ElasticityParameters m_elas_parameters;

 public:
   linear_elasticity_solver(const mesh_type&           msh,
                            const bnd_type&            bnd,
                            const ElasticityParameters data,
                            size_t                     degree,
                            int                        l = 0) :
     m_msh(msh),
     m_verbose(false), m_bnd(bnd)
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

      m_elas_parameters.mu     = data.mu;
      m_elas_parameters.lambda = data.lambda;

      m_bqd = bqdata_type(m_face_degree, m_cell_degree);
   }

   void
   changeElasticityParameters(const ElasticityParameters data)
   {
      m_elas_parameters.mu     = data.mu;
      m_elas_parameters.lambda = data.lambda;
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

   size_t
   getDofs()
   {
      return m_msh.faces_size() * disk::howmany_dofs(m_bqd.face_basis);
   }

   template<typename LoadFunction>
   assembly_info
   assemble(const LoadFunction& lf)
   {
      auto sgradrec  = sgradrec_type(m_bqd);
      auto divrec    = divrec_type(m_bqd);
      auto stab      = stab_type(m_bqd);
      auto statcond  = statcond_type(m_bqd);
      auto assembler = assembler_type(m_msh, m_bqd, m_bnd);

      assembly_info ai;
      bzero(&ai, sizeof(ai));

      timecounter tc;

      for (auto& cl : m_msh) {
         tc.tic();
         sgradrec.compute(m_msh, cl);
         tc.toc();
         ai.time_gradrec += tc.to_double();

         tc.tic();
         divrec.compute(m_msh, cl);
         tc.toc();
         ai.time_divrec += tc.to_double();

         tc.tic();
         stab.compute(m_msh, cl, sgradrec.oper);
         tc.toc();
         ai.time_stab += tc.to_double();

         tc.tic();
         const auto cell_rhs = disk::hho::compute_rhs(m_msh, cl, lf, m_bqd, m_cell_degree);
         const dynamic_matrix<scalar_type> loc =
           2.0 * m_elas_parameters.mu * (sgradrec.data + stab.data) +
           m_elas_parameters.lambda * divrec.data;
         const auto scnp = statcond.compute(m_msh, cl, loc, cell_rhs);
         tc.toc();
         ai.time_statcond += tc.to_double();

         assembler.assemble(m_msh, cl, scnp);
      }

      assembler.impose_neumann_boundary_conditions(m_msh);
      assembler.finalize(m_system_matrix, m_system_rhs);

      ai.linear_system_size = m_system_matrix.rows();
      ai.time_assembly      = ai.time_gradrec + ai.time_divrec + ai.time_stab + ai.time_statcond;
      return ai;
   }

   solver_info
   solve(void)
   {
#ifdef HAVE_INTEL_MKL
      Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>> solver;
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
   postprocess_info
   postprocess(const LoadFunction& lf)
   {
      auto sgradrec  = sgradrec_type(m_bqd);
      auto divrec    = divrec_type(m_bqd);
      auto stab      = stab_type(m_bqd);
      auto statcond  = statcond_type(m_bqd);
      auto assembler = assembler_type(m_msh, m_bqd, m_bnd);

      const size_t fbs = disk::howmany_dofs(m_bqd.face_basis);

      postprocess_info pi;

      m_solution_data.reserve(m_msh.cells_size());

      const auto solF = assembler.expand_solution(m_msh, m_system_solution);

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
            xF                                 = solF.block(face_id * fbs, 0, fbs, 1);
            xFs.block(face_i * fbs, 0, fbs, 1) = xF;
         }

         sgradrec.compute(m_msh, cl);
         divrec.compute(m_msh, cl);
         stab.compute(m_msh, cl, sgradrec.oper);
         const dynamic_matrix<scalar_type> loc =
           2.0 * m_elas_parameters.mu * (sgradrec.data + stab.data) +
           m_elas_parameters.lambda * divrec.data;
         const auto cell_rhs = disk::hho::compute_rhs(m_msh, cl, lf, m_bqd, m_cell_degree);

         const dynamic_vector<scalar_type> x = statcond.recover(m_msh, cl, loc, cell_rhs, xFs);
         m_solution_data.push_back(x);
      }
      tc.toc();

      pi.time_postprocess = tc.to_double();

      return pi;
   }

   template<typename AnalyticalSolution>
   scalar_type
   compute_l2_displacement_error(const AnalyticalSolution& as)
   {
      scalar_type err_dof = 0;

      projector_type projk(m_bqd);

      size_t cell_i = 0;

      for (auto& cl : m_msh) {
         const auto                        x        = m_solution_data.at(cell_i++);
         const dynamic_vector<scalar_type> true_dof = projk.projectOnCell(m_msh, cl, as);
         const dynamic_vector<scalar_type> comp_dof = x.block(0, 0, true_dof.size(), 1);
         const dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
         err_dof += diff_dof.dot(projk.cell_mm * diff_dof);
      }

      return sqrt(err_dof);
   }

   template<typename AnalyticalSolution>
   scalar_type
   compute_l2_gradient_error(const AnalyticalSolution& grad)
   {
      scalar_type err_dof = scalar_type{0.0};

      projector_type projk(m_bqd);

      sgradrec_type gradrec(m_bqd);

      size_t i = 0;

      for (auto& cl : m_msh) {
         const auto x = m_solution_data.at(i++);
         gradrec.compute(m_msh, cl);
         const dynamic_vector<scalar_type> RTu = gradrec.oper * x;

         const dynamic_vector<scalar_type> true_dof =
           projk.projectOnSymStiffnessSpace(m_msh, cl, grad);
         const dynamic_vector<scalar_type> comp_dof = RTu;
         assert(true_dof.size() == comp_dof.size());
         const dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
         err_dof += diff_dof.dot(projk.grad_mm * diff_dof);
      }

      return sqrt(err_dof);
   }

   template<typename AnalyticalSolution>
   scalar_type
   compute_l2_stress_error(const AnalyticalSolution& stress) const
   {
      auto sgradrec = sgradrec_type(m_bqd);
      auto divrec   = divrec_type(m_bqd);

      projector_type projk(m_bqd);

      const auto cell_degree = m_bqd.cell_degree();

      size_t      cell_i(0);
      scalar_type error_stress(0.0);

      for (auto& cl : m_msh) {
         const auto x = m_solution_data.at(cell_i++);
         sgradrec.compute(m_msh, cl);

         const dynamic_vector<scalar_type> GTu = sgradrec.oper * x;

         divrec.compute(m_msh, cl);

         const dynamic_vector<scalar_type> divu = divrec.oper * x;

         const auto cell_quadpoints = m_bqd.cell_quadrature.integrate(m_msh, cl);

         for (auto& qp : cell_quadpoints) {
            const auto gphi =
              m_bqd.cell_basis.eval_sgradients(m_msh, cl, qp.point(), 1, cell_degree + 1);
            const auto GT_iqn = disk::hho::eval(GTu, gphi);

            const auto divphi =
              m_bqd.div_cell_basis.eval_functions(m_msh, cl, qp.point(), 0, cell_degree);
            const auto divu_iqn = disk::hho::eval(divu, divphi);

            const auto sigma = 2.0 * m_elas_parameters.mu * GT_iqn +
                               m_elas_parameters.lambda * divu_iqn *
                                 static_matrix<scalar_type, dimension, dimension>::Identity();

            const auto stress_diff = (stress(qp.point()) - sigma).eval();

            error_stress += qp.weight() * disk::mm_prod(stress_diff, stress_diff);
         }
      }

      return sqrt(error_stress);
   }

   // post-processing
   void
   compute_discontinuous_displacement(const std::string& filename) const
   {
      const size_t cell_degree = m_bqd.cell_degree();
      const size_t cbs         = disk::howmany_dofs(m_bqd.cell_basis);

      gmsh::Gmesh gmsh(dimension);
      auto        storage = m_msh.backend_storage();

      std::vector<gmsh::Data>          data;    // create data (not used)
      const std::vector<gmsh::SubData> subdata; // create subdata to save soution at gauss point

      size_t cell_i   = 0;
      size_t nb_nodes = 0;
      for (auto& cl : m_msh) {
         const vector_dynamic    x          = m_solution_data.at(cell_i++).head(cbs);
         const auto              cell_nodes = disk::cell_nodes(m_msh, cl);
         std::vector<gmsh::Node> new_nodes;

         // loop on the nodes of the cell
         for (size_t i = 0; i < cell_nodes.size(); i++) {
            nb_nodes++;
            const auto point_ids = cell_nodes[i];
            const auto pt        = storage->points[point_ids];

            const auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt, 0, cell_degree);

            const auto depl = disk::hho::eval(x, phi);

            const std::vector<double>   deplv = disk::convertToVectorGmsh(depl);
            const std::array<double, 3> coor  = disk::init_coor(pt);

            // Add a node
            const gmsh::Node tmp_node(coor, nb_nodes, 0);
            new_nodes.push_back(tmp_node);
            gmsh.addNode(tmp_node);

            const gmsh::Data datatmp(nb_nodes, deplv);
            data.push_back(datatmp);
         }
         // Add new element
         disk::add_element(gmsh, new_nodes);
      }

      // Create and init a nodedata view
      gmsh::NodeData nodedata(3, 0.0, "depl_node", data, subdata);

      // Save the view
      nodedata.saveNodeData(filename, gmsh);
   }

   void
   compute_discontinuous_stress(const std::string& filename) const
   {
      auto sgradrec = sgradrec_type(m_bqd);
      auto divrec   = divrec_type(m_bqd);

      const size_t cbs         = disk::howmany_dofs(m_bqd.cell_basis);
      const size_t cell_degree = m_bqd.cell_degree();

      gmsh::Gmesh gmsh(dimension);
      auto        storage = m_msh.backend_storage();

      std::vector<gmsh::Data>          data;    // create data (not used)
      const std::vector<gmsh::SubData> subdata; // create subdata to save soution at gauss point

      size_t cell_i   = 0;
      size_t nb_nodes = 0;
      for (auto& cl : m_msh) {
         const auto x = m_solution_data.at(cell_i++).head(cbs);
         sgradrec.compute(m_msh, cl);

         const dynamic_vector<scalar_type> GsTu = sgradrec.oper * x;

         divrec.compute(m_msh, cl);

         const dynamic_vector<scalar_type> Divu = divrec.oper * x;

         const auto              cell_nodes = disk::cell_nodes(m_msh, cl);
         std::vector<gmsh::Node> new_nodes;

         // loop on the nodes of the cell
         for (size_t i = 0; i < cell_nodes.size(); i++) {
            nb_nodes++;
            const auto point_ids = cell_nodes[i];
            const auto pt        = storage->points[point_ids];

            const auto sdphi = m_bqd.cell_basis.eval_sgradients(m_msh, cl, pt, 0, cell_degree + 1);
            const auto gsT   = disk::hho::eval(GsTu, sdphi);

            const auto divphi = m_bqd.div_cell_basis.eval_functions(m_msh, cl, pt, 0, cell_degree);
            const auto divu   = disk::hho::eval(Divu, divphi);

            const static_matrix<scalar_type, dimension, dimension> stress =
              2 * m_elas_parameters.mu * gsT +
              m_elas_parameters.lambda * divu *
                static_matrix<scalar_type, dimension, dimension>::Identity();

            const std::vector<double>   tens = disk::convertToVectorGmsh(stress);
            const std::array<double, 3> coor = disk::init_coor(pt);

            // Add a node
            const gmsh::Node tmp_node(coor, nb_nodes, 0);
            new_nodes.push_back(tmp_node);
            gmsh.addNode(tmp_node);

            const gmsh::Data datatmp(nb_nodes, tens);
            data.push_back(datatmp);
         }

         // Add new element
         disk::add_element(gmsh, new_nodes);
         cell_i++;
      }

      // Create and init a nodedata view
      gmsh::NodeData nodedata(9, 0.0, "stress_node", data, subdata);

      // Save the view
      nodedata.saveNodeData(filename, gmsh);
   }
};
