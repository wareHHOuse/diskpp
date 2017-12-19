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

#include "hho/assembler.hpp"
#include "hho/gradient_reconstruction.hpp"
#include "hho/hho_bq.hpp"
#include "hho/hho_utils.hpp"
#include "hho/projector.hpp"
#include "hho/stabilization.hpp"
#include "hho/static_condensation.hpp"

#include "timecounter.h"

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

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

   typedef disk::hho::gradient_reconstruction_bq<bqdata_type> gradrec_type;
   typedef disk::hho::hho_stabilization_bq<bqdata_type>       stab_type;
   typedef disk::hho::static_condensation_bq<bqdata_type>     statcond_type;
   typedef disk::hho::assembler_bq<bqdata_type>               assembler_type;

   size_t m_cell_degree, m_face_degree;

   bqdata_type m_bqd;

   typename assembler_type::sparse_matrix_type m_system_matrix;
   typename assembler_type::vector_type        m_system_rhs, m_system_solution;

   const mesh_type& m_msh;

   std::vector<dynamic_vector<scalar_type>> m_postprocess_data;

   bool m_verbose;

   typedef dynamic_vector<scalar_type> vector_dynamic;

 public:
   diffusion_solver(const mesh_type& msh, size_t degree, int l = 0) : m_msh(msh), m_verbose(false)
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

      m_bqd = bqdata_type(m_cell_degree, m_face_degree);
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
   assembly_info
   assemble(const LoadFunction& lf, const BoundaryConditionFunction& bcf)
   {
      auto gradrec   = gradrec_type(m_bqd);
      auto stab      = stab_type(m_bqd);
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
         stab.compute(m_msh, cl, gradrec.oper);
         tc.toc();
         ai.time_stab += tc.to_double();

         tc.tic();
         const auto cell_rhs = disk::hho::compute_rhs(m_msh, cl, lf, m_bqd, m_bqd.cell_degree());
         const dynamic_matrix<scalar_type> loc  = gradrec.data + stab.data;
         const auto                        scnp = statcond.compute(m_msh, cl, loc, cell_rhs);
         tc.toc();
         ai.time_statcond += tc.to_double();

         assembler.assemble(m_msh, cl, scnp);
      }

      assembler.impose_boundary_conditions(m_msh, bcf);
      assembler.finalize(m_system_matrix, m_system_rhs);

      ai.linear_system_size = m_system_matrix.rows();
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
   postprocess_info
   postprocess(const LoadFunction& lf)
   {
      auto gradrec  = gradrec_type(m_bqd);
      auto stab     = stab_type(m_bqd);
      auto statcond = statcond_type(m_bqd);

      const size_t fbs = howmany_dofs(m_bqd.face_basis);

      postprocess_info pi;

      m_postprocess_data.reserve(m_msh.cells_size());

      timecounter tc;
      tc.tic();
      for (auto& cl : m_msh) {
         const auto fcs       = faces(m_msh, cl);
         const auto num_faces = fcs.size();

         dynamic_vector<scalar_type> xFs = dynamic_vector<scalar_type>::Zero(num_faces * fbs);

         for (size_t face_i = 0; face_i < num_faces; face_i++) {
            auto       fc  = fcs[face_i];
            const auto eid = find_element_id(m_msh.faces_begin(), m_msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

            const auto face_id = eid.second;

            dynamic_vector<scalar_type> xF     = dynamic_vector<scalar_type>::Zero(fbs);
            xF                                 = m_system_solution.block(face_id * fbs, 0, fbs, 1);
            xFs.block(face_i * fbs, 0, fbs, 1) = xF;
         }

         gradrec.compute(m_msh, cl);
         stab.compute(m_msh, cl, gradrec.oper);
         const dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
         const auto cell_rhs = disk::hho::compute_rhs(m_msh, cl, lf, m_bqd, m_bqd.cell_degree());
         const dynamic_vector<scalar_type> x = statcond.recover(m_msh, cl, loc, cell_rhs, xFs);
         m_postprocess_data.push_back(x);
      }
      tc.toc();

      pi.time_postprocess = tc.to_double();

      return pi;
   }

   template<typename AnalyticalSolution>
   scalar_type
   compute_l2_error(const AnalyticalSolution& as)
   {
      scalar_type err_dof = 0.0;

      disk::hho::projector_bq<bqdata_type> projk(m_bqd);

      size_t i = 0;
      for (auto& cl : m_msh) {
         const auto                        x        = m_postprocess_data.at(i++);
         const dynamic_vector<scalar_type> true_dof = projk.projectOnCell(m_msh, cl, as);
         const dynamic_vector<scalar_type> comp_dof = x.block(0, 0, true_dof.size(), 1);
         const dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
         err_dof += diff_dof.dot(projk.cell_mm * diff_dof);
      }

      return sqrt(err_dof);
   }

   template<typename AnalyticalSolution>
   void
   plot_l2error_at_gausspoint(const std::string& filename, const AnalyticalSolution& as)
   {
      const size_t cbs         = disk::howmany_dofs(m_bqd.cell_basis);
      const size_t cell_degree = m_bqd.cell_degree();

      std::cout << "Compute L2 error at Gauss points" << std::endl;
      gmsh::Gmesh msh(m_msh.dimension); // creta a mesh

      std::vector<gmsh::Data>    data;    // create data (not used)
      std::vector<gmsh::SubData> subdata; // create subdata to save soution at gauss point
      size_t                     nb_node = msh.getNumberofNodes();

      size_t cell_i = 0;
      for (auto& cl : m_msh) {
         const vector_dynamic x   = m_postprocess_data.at(cell_i++).head(cbs);
         const auto           qps = m_bqd.cell_quadrature.integrate(m_msh, cl);
         for (auto& qp : qps) {
            const auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, qp.point(), 0, cell_degree);

            const scalar_type pot = disk::hho::eval(x, phi);

            const scalar_type true_pot = as(qp.point()); // a voir et projeté

            nb_node += 1;
            const gmsh::Node snode =
              disk::convertPoint(qp.point(), nb_node); // create a node at gauss point
            const std::vector<double> value = {
              std::abs(pot - true_pot)}; // save the solution at gauss point
            const gmsh::SubData sdata(value, snode);
            subdata.push_back(sdata); // add subdata
         }
      }
      gmsh::NodeData nodedata(
        1, 0.0, "sol_scalar", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, msh); // save the view
   }

   void
   plot_solution_at_gausspoint(const std::string& filename)
   {
      const size_t cbs         = disk::howmany_dofs(m_bqd.cell_basis);
      const size_t cell_degree = m_bqd.cell_degree();

      std::cout << "Compute solution at Gauss points" << std::endl;
      gmsh::Gmesh msh(m_msh.dimension); // creta a mesh

      std::vector<gmsh::Data>    data;    // create data (not used)
      std::vector<gmsh::SubData> subdata; // create subdata to save soution at gauss point
      size_t                     nb_node = msh.getNumberofNodes();

      size_t cell_i = 0;
      for (auto& cl : m_msh) {
         const vector_dynamic x   = m_postprocess_data.at(cell_i++).head(cbs);
         const auto           qps = m_bqd.cell_quadrature.integrate(m_msh, cl);
         for (auto& qp : qps) {
            const auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, qp.point(), 0, cell_degree);

            const scalar_type pot = disk::hho::eval(x, phi);

            nb_node += 1;
            gmsh::Node snode =
              disk::convertPoint(qp.point(), nb_node); // create a node at gauss point
            std::vector<double> value = {double(pot)}; // save the solution at gauss point
            gmsh::SubData       sdata(value, snode);
            subdata.push_back(sdata); // add subdata
         }
      }

      gmsh::NodeData nodedata(
        1, 0.0, "sol_scalar", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, msh); // save the view
   }

   void
   plot_conforme_solution(const std::string& filename)
   {
      const size_t cbs         = disk::howmany_dofs(m_bqd.cell_basis);
      const size_t cell_degree = m_bqd.cell_degree();

      gmsh::Gmesh gmsh    = disk::convertMesh(m_msh);
      auto        storage = m_msh.backend_storage();
      size_t      nb_nodes(gmsh.getNumberofNodes());

      // first(number of data at this node), second(cumulated value)
      std::vector<std::pair<size_t, scalar_type>> value(nb_nodes, std::make_pair(0, double(0.0)));

      size_t cell_i(0);
      for (auto& cl : m_msh) {
         const vector_dynamic    x          = m_postprocess_data.at(cell_i++).head(cbs);
         const auto              cell_nodes = disk::cell_nodes(m_msh, cl);
         std::vector<gmsh::Node> new_nodes;
         for (size_t i = 0; i < cell_nodes.size(); i++) {
            nb_nodes++;
            const auto point_ids = cell_nodes[i];
            const auto pt        = storage->points[point_ids];

            const auto        phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt, 0, cell_degree);
            const scalar_type pot = disk::hho::eval(x, phi);

            value[point_ids].first += 1;
            value[point_ids].second += pot;
         }
      }

      std::vector<gmsh::Data>    data;    // create data
      std::vector<gmsh::SubData> subdata; // create subdata
      data.reserve(nb_nodes);             // data has a size of nb_node

      for (size_t i_node = 0; i_node < value.size(); i_node++) {
         const std::vector<double> tmp_value(1, value[i_node].second / double(value[i_node].first));
         const gmsh::Data          tmp_data(i_node + 1, tmp_value);
         data.push_back(tmp_data); // add data
      }

      gmsh::NodeData nodedata(1, 0.0, "sol_node", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }

   void
   plot_discontinuous_solution(const std::string& filename)
   {
      const size_t cbs         = disk::howmany_dofs(m_bqd.cell_basis);
      const size_t cell_degree = m_bqd.cell_degree();

      const size_t dim = m_msh.dimension;
      gmsh::Gmesh  gmsh(dim);
      auto         storage = m_msh.backend_storage();

      std::vector<gmsh::Data>    data;    // create data (not used)
      std::vector<gmsh::SubData> subdata; // create subdata to save soution at gauss point

      size_t cell_i(0);
      size_t nb_nodes(0);
      for (auto& cl : m_msh) {
         const vector_dynamic    x          = m_postprocess_data.at(cell_i++).head(cbs);
         const auto              cell_nodes = disk::cell_nodes(m_msh, cl);
         std::vector<gmsh::Node> new_nodes;
         for (size_t i = 0; i < cell_nodes.size(); i++) {
            nb_nodes++;
            const auto point_ids = cell_nodes[i];
            const auto pt        = storage->points[point_ids];

            const auto        phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt, 0, cell_degree);
            const scalar_type pot = disk::hho::eval(x, phi);

            const std::array<double, 3> coor = disk::init_coor(pt);

            const gmsh::Node tmp_node(coor, nb_nodes, 0);
            new_nodes.push_back(tmp_node);
            gmsh.addNode(tmp_node);

            const std::vector<double> value = {double(pot)};
            const gmsh::Data          datatmp(nb_nodes, value);
            data.push_back(datatmp);
         }
         // add new element
         disk::add_element(gmsh, new_nodes);
      }

      gmsh::NodeData nodedata(1, 0.0, "sol_node", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }

   void
   plot_deformed_conforme(const std::string& filename)
   {
      const size_t cbs         = disk::howmany_dofs(m_bqd.cell_basis);
      const size_t cell_degree = m_bqd.cell_degree();

      static_assert(mesh_type::dimension < 3, "dim<=3");

      gmsh::Gmesh gmsh    = disk::convertMesh(m_msh);
      auto        storage = m_msh.backend_storage();
      size_t      nb_nodes(gmsh.getNumberofNodes());

      // first(number of data at this node), second(cumulated value)
      std::vector<std::pair<size_t, scalar_type>> value(nb_nodes, std::make_pair(0, double(0.0)));

      size_t cell_i(0);
      for (auto& cl : m_msh) {
         const vector_dynamic    x          = m_postprocess_data.at(cell_i++).head(cbs);
         const auto              cell_nodes = disk::cell_nodes(m_msh, cl);
         std::vector<gmsh::Node> new_nodes;
         for (size_t i = 0; i < cell_nodes.size(); i++) {
            nb_nodes++;
            const auto point_ids = cell_nodes[i];
            const auto pt        = storage->points[point_ids];

            const auto        phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt, 0, cell_degree);
            const scalar_type pot = disk::hho::eval(x, phi);

            value[point_ids].first += 1;
            value[point_ids].second += pot;
         }
      }

      std::vector<gmsh::Data>    data;    // create data
      std::vector<gmsh::SubData> subdata; // create subdata
      data.reserve(nb_nodes);             // data has a size of nb_node

      for (size_t i_node = 0; i_node < value.size(); i_node++) {
         std::vector<double> tmp_value(3, 0.0);
         tmp_value[m_msh.dimension] = value[i_node].second / double(value[i_node].first);
         const gmsh::Data tmp_data(i_node + 1, tmp_value);
         data.push_back(tmp_data); // add data
      }

      gmsh::NodeData nodedata(3, 0.0, "sol_node", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }

   void
   plot_deformed_discontinuous(const std::string& filename)
   {
      const size_t cbs         = disk::howmany_dofs(m_bqd.cell_basis);
      const size_t cell_degree = m_bqd.cell_degree();

      static_assert(mesh_type::dimension < 3, "dim<=3");

      gmsh::Gmesh gmsh(m_msh.dimension);
      auto        storage = m_msh.backend_storage();

      std::vector<gmsh::Data>    data;    // create data (not used)
      std::vector<gmsh::SubData> subdata; // create subdata to save soution at gauss point

      size_t cell_i(0);
      size_t nb_nodes(0);
      for (auto& cl : m_msh) {
         const vector_dynamic    x          = m_postprocess_data.at(cell_i++).head(cbs);
         const auto              cell_nodes = disk::cell_nodes(m_msh, cl);
         std::vector<gmsh::Node> new_nodes;
         for (size_t i = 0; i < cell_nodes.size(); i++) {
            nb_nodes++;
            const auto point_ids = cell_nodes[i];
            const auto pt        = storage->points[point_ids];

            const auto        phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt, 0, cell_degree);
            const scalar_type pot = disk::hho::eval(x, phi);

            const std::array<double, 3> coor = disk::init_coor(pt);

            const gmsh::Node tmp_node(coor, nb_nodes, 0);
            new_nodes.push_back(tmp_node);
            gmsh.addNode(tmp_node);

            std::vector<double> value(3, 0.0);
            value[m_msh.dimension + 1 - 1] = {double(pot)};
            gmsh::Data datatmp(nb_nodes, value);

            data.push_back(datatmp);
         }
         // add new element
         disk::add_element(gmsh, new_nodes);
      }

      gmsh::NodeData nodedata(
        3, 0.0, "sol_node_deformed", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }

   void
   saveMesh(const std::string& filename)
   {
      gmsh::Gmesh gmsh = disk::convertMesh(m_msh);
      gmsh.writeGmesh(filename, 2);
   }
};
