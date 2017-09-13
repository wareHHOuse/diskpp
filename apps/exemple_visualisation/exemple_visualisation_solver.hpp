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


#ifdef HAVE_SOLVER_WRAPPERS
#include "agmg/agmg.hpp"
#endif

#include "hho/hho.hpp"

#include "timecounter.h"

#include "visualisation/gmshConvertMesh.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

struct assembly_info
{
    size_t  linear_system_size;
    double  time_gradrec, time_statcond, time_stab;
};

struct solver_info
{
    double  time_solver;
};

struct postprocess_info
{
    double  time_postprocess;
};

template<typename Mesh>
class diffusion_solver
{
    typedef Mesh                                       mesh_type;
    typedef typename mesh_type::scalar_type            scalar_type;
    typedef typename mesh_type::cell                   cell_type;
    typedef typename mesh_type::face                   face_type;

    typedef disk::quadrature<mesh_type, cell_type>      cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>      face_quadrature_type;

    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

    typedef
    disk::basis_quadrature_data<mesh_type,
                                disk::scaled_monomial_scalar_basis,
                                disk::quadrature> bqdata_type;

    typedef disk::gradient_reconstruction_bq<bqdata_type>               gradrec_type;
    typedef disk::diffusion_like_stabilization_bq<bqdata_type>          stab_type;
    typedef disk::diffusion_like_static_condensation_bq<bqdata_type>    statcond_type;
    typedef disk::assembler<mesh_type, face_basis_type, face_quadrature_type> assembler_type;

    size_t m_cell_degree, m_face_degree;

    bqdata_type     m_bqd;

    typename assembler_type::sparse_matrix_type     m_system_matrix;
    typename assembler_type::vector_type            m_system_rhs, m_system_solution;

    const mesh_type& m_msh;

    std::vector<dynamic_vector<scalar_type>>        m_postprocess_data;

    bool m_verbose;

    typedef dynamic_vector<scalar_type>         vector_dynamic;

public:
   diffusion_solver(const mesh_type& msh, size_t degree, int l = 0)
   : m_msh(msh), m_verbose(false)
   {
      if ( l < -1 or l > 1)
      {
         std::cout << "'l' should be -1, 0 or 1. Reverting to 0." << std::endl;
         l = 0;
      }

      if (degree == 0 && l == -1)
      {
         std::cout << "'l' should be 0 or 1. Reverting to 0." << std::endl;
         l = 0;
      }

      m_cell_degree = degree + l;
      m_face_degree = degree;

      m_bqd = bqdata_type(m_cell_degree, m_face_degree);
   }

   bool    verbose(void) const     { return m_verbose; }
   void    verbose(bool v)         { m_verbose = v; }

   template<typename LoadFunction, typename BoundaryConditionFunction>
   assembly_info
   assemble(const LoadFunction& lf, const BoundaryConditionFunction& bcf)
   {
      auto gradrec    = gradrec_type(m_bqd);
      auto stab       = stab_type(m_bqd);
      auto statcond   = statcond_type(m_bqd);
      auto assembler  = assembler_type(m_msh, m_face_degree);

      assembly_info ai;
      bzero(&ai, sizeof(ai));

      timecounter tc;

      for (auto& cl : m_msh)
      {
         tc.tic();
         gradrec.compute(m_msh, cl);
         tc.toc();
         ai.time_gradrec += tc.to_double();

         tc.tic();
         stab.compute(m_msh, cl, gradrec.oper);
         tc.toc();
         ai.time_stab += tc.to_double();

         tc.tic();
         auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(m_msh, cl, lf, m_cell_degree);
         dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
         auto scnp = statcond.compute(m_msh, cl, loc, cell_rhs);
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
      Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
      #else
      Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
      #endif

      solver_info si;

      size_t systsz = m_system_matrix.rows();
      size_t nnz = m_system_matrix.nonZeros();

      if (verbose())
      {
         std::cout << "Starting linear solver..." << std::endl;
         std::cout << " * Solving for " << systsz << " unknowns." << std::endl;
         std::cout << " * Matrix fill: " << 100.0*double(nnz)/(systsz*systsz) << "%" << std::endl;
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
      auto gradrec    = gradrec_type(m_bqd);
      auto stab       = stab_type(m_bqd);
      auto statcond   = statcond_type(m_bqd);

      size_t fbs = m_bqd.face_basis.size();

      postprocess_info pi;

      m_postprocess_data.reserve(m_msh.cells_size());

      timecounter tc;
      tc.tic();
      for (auto& cl : m_msh)
      {
         auto fcs = faces(m_msh, cl);
         auto num_faces = fcs.size();

         dynamic_vector<scalar_type> xFs = dynamic_vector<scalar_type>::Zero(num_faces*fbs);

         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            auto fc = fcs[face_i];
            auto eid = find_element_id(m_msh.faces_begin(), m_msh.faces_end(), fc);
            if (!eid.first)
            throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            dynamic_vector<scalar_type> xF = dynamic_vector<scalar_type>::Zero(fbs);
            xF = m_system_solution.block(face_id * fbs, 0, fbs, 1);
            xFs.block(face_i * fbs, 0, fbs, 1) = xF;
         }

         gradrec.compute(m_msh, cl);
         stab.compute(m_msh, cl, gradrec.oper);
         dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
         auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(m_msh, cl, lf, m_cell_degree);
         dynamic_vector<scalar_type> x = statcond.recover(m_msh, cl, loc, cell_rhs, xFs);
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

      disk::projector_bq<bqdata_type> projk(m_bqd);

      size_t i = 0;
      for (auto& cl : m_msh)
      {
         auto x = m_postprocess_data.at(i++);
         dynamic_vector<scalar_type> true_dof = projk.compute_cell(m_msh, cl, as);
         dynamic_vector<scalar_type> comp_dof = x.block(0,0,true_dof.size(), 1);
         dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
         err_dof += diff_dof.dot(projk.cell_mm * diff_dof);
      }

      return sqrt(err_dof);
   }

   template<typename AnalyticalSolution>
   void
   plot_l2error_at_gausspoint(const std::string& filename, const AnalyticalSolution& as)
   {
      std::cout << "Compute L2 error at Gauss points" << std::endl;
      visu::Gmesh msh(m_msh.dimension); //creta a mesh

      std::vector<visu::Data> data; //create data (not used)
      std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point
      size_t nb_node =  msh.getNumberofNodes();

      size_t cell_i = 0;
      for (auto& cl : m_msh)
      {
         auto x = m_postprocess_data.at(cell_i++);
         auto qps = m_bqd.cell_quadrature.integrate(m_msh, cl);
         for (auto& qp : qps)
         {
            auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, qp.point());

            scalar_type pot = 0.0;
            for (size_t i = 0; i < m_bqd.cell_basis.range(0, m_cell_degree).size(); i++){
               pot += phi[i] * x(i);
            }

            scalar_type true_pot = as(qp.point()); // a voir et projeté

            nb_node += 1;
            visu::Node snode = visu::convertPoint(qp.point(), nb_node); //create a node at gauss point
            std::vector<double> value = {std::fabs(pot - true_pot)}; // save the solution at gauss point
            visu::SubData sdata(value, snode);
            subdata.push_back(sdata); // add subdata
         }
      }
      visu::NodeData nodedata(1, 0.0, "sol_scalar", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, msh); // save the view
   }

   void
   plot_solution_at_gausspoint(const std::string& filename)
   {
      std::cout << "Compute solution at Gauss points" << std::endl;
      visu::Gmesh msh(m_msh.dimension); //creta a mesh

      std::vector<visu::Data> data; //create data (not used)
      std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point
      size_t nb_node =  msh.getNumberofNodes();

      size_t cell_i = 0;
      for (auto& cl : m_msh)
      {
         auto x = m_postprocess_data.at(cell_i++);
         auto qps = m_bqd.cell_quadrature.integrate(m_msh, cl);
         for (auto& qp : qps)
         {
            auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, qp.point());

            scalar_type pot = 0.0;
            for (size_t i = 0; i < m_bqd.cell_basis.range(0, m_cell_degree).size(); i++)
               pot += phi[i] * x(i);

            visu::Node snode = visu::convertPoint(qp.point(), nb_node++); //create a node at gauss point
            std::vector<double> value = {double(pot)}; // save the solution at gauss point
            visu::SubData sdata(value, snode);
            subdata.push_back(sdata); // add subdata
         }
      }

      visu::NodeData nodedata(1, 0.0, "sol_scalar", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, msh); // save the view
   }

   void
   plot_conforme_solution(const std::string& filename)
   {
      visu::Gmesh gmsh = visu::convertMesh(m_msh);
      auto storage = m_msh.backend_storage();
      size_t nb_nodes(gmsh.getNumberofNodes());

      //first(number of data at this node), second(cumulated value)
      std::vector<std::pair<size_t, scalar_type > > value(nb_nodes, std::make_pair(0, double(0.0)));

      size_t cell_i(0);
      for (auto& cl : m_msh)
      {
         vector_dynamic x = m_postprocess_data.at(cell_i++);
         auto cell_nodes = cl.point_ids();
         std::vector<visu::Node> new_nodes;
         for (size_t i = 0; i < cell_nodes.size(); i++)
         {
            nb_nodes++;
            auto point_ids = cell_nodes[i];
            auto pt = storage->points[point_ids];

            auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt);

            scalar_type pot(0.0);

            //compute solution at the node
            for (size_t j = 0; j < m_bqd.cell_basis.range(0, m_cell_degree).size(); j++)
               pot += phi[j] * x(j);

            value[point_ids].first += 1;
            value[point_ids].second += pot;
         }
      }

      std::vector<visu::Data> data; //create data
      std::vector<visu::SubData> subdata; //create subdata
      data.reserve(nb_nodes); // data has a size of nb_node

      for(size_t  i_node = 0; i_node < value.size(); i_node++){
         std::vector<double> tmp_value(1, value[i_node].second/ double(value[i_node].first));
         visu::Data tmp_data(i_node + 1, tmp_value);
         data.push_back(tmp_data); //add data
      }

      visu::NodeData nodedata(1, 0.0, "sol_node", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }

   void
   plot_discontinuous_solution(const std::string& filename)
   {
      visu::Gmesh gmsh(m_msh.dimension);
      std::vector<visu::Data> data; //create data

      size_t cell_i(0);
      size_t nb_nodes(0);
      for (auto& cl : m_msh)
      {
         vector_dynamic x = m_postprocess_data.at(cell_i++);
         auto cell_nodes = points(m_msh, cl);
         std::vector<visu::Node> new_nodes;
         for (auto& pt : cell_nodes)
         {
            nb_nodes++;

            auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt);

            std::array<double, 3> coor = {double{0.0}, double{0.0}, double{0.0}};

            visu::init_coordinate(pt, coor);
            visu::Node tmp_node(coor, nb_nodes, 0);
            new_nodes.push_back(tmp_node);
            gmsh.addNode(tmp_node);

            // plot magnitude at node
            scalar_type pot(0.0);
            for (size_t j = 0; j < m_bqd.cell_basis.range(0, m_cell_degree).size(); j++) //compute solution at the node
               pot += phi[j] * x(j);

            std::vector<double> value = {double(pot)};
            visu::Data datatmp(nb_nodes, value);
            data.push_back(datatmp);
         }
         // add new element
         visu::add_element(gmsh, new_nodes);
      }
      std::vector<visu::SubData> subdata; //(not used)
      visu::NodeData nodedata(1, 0.0, "sol_node", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }

   void
   plot_deformed_conforme(const std::string& filename)
   {
      const size_t DIM = m_msh.dimension;
      if(DIM >= 3)
      std::cout << "Compute deformed only in 1D or 2D" << '\n';
      else {
         visu::Gmesh gmsh = visu::convertMesh(m_msh);
         auto storage = m_msh.backend_storage();
         size_t nb_nodes(gmsh.getNumberofNodes());

         //first(number of data at this node), second(cumulated value)
         std::vector<std::pair<size_t, scalar_type > > value(nb_nodes, std::make_pair(0, double(0.0)));

         size_t cell_i(0);
         for (auto& cl : m_msh)
         {
            vector_dynamic x = m_postprocess_data.at(cell_i++);
            auto cell_nodes = cl.point_ids();
            std::vector<visu::Node> new_nodes;
            for (size_t i = 0; i < cell_nodes.size(); i++)
            {
               nb_nodes++;
               auto point_ids = cell_nodes[i];
               auto pt = storage->points[point_ids];

               auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt);

               scalar_type pot(0.0);

               //compute solution at the node
               for (size_t j = 0; j < m_bqd.cell_basis.range(0, m_cell_degree).size(); j++)
                  pot += phi[j] * x(j);

               value[point_ids].first +=1;
               value[point_ids].second += pot;
            }
         }

         std::vector<visu::Data> data; //create data
         std::vector<visu::SubData> subdata; //create subdata
         data.reserve(nb_nodes); // data has a size of nb_node

         for(size_t  i_node = 0; i_node < value.size(); i_node++){
            std::vector<double> tmp_value(3, 0.0);
            tmp_value[DIM] = value[i_node].second/ double(value[i_node].first);
            visu::Data tmp_data(i_node + 1, tmp_value);
            data.push_back(tmp_data); //add data
         }

         visu::NodeData nodedata(3, 0.0, "sol_node", data, subdata); // create and init a nodedata view

         nodedata.saveNodeData(filename, gmsh); // save the view
      }
   }

   void
   plot_deformed_discontinuous(const std::string& filename)
   {
      const size_t DIM = m_msh.dimension;
      if(DIM >= 3)
         std::cout << "Compute deformed only in 1D or 2D" << '\n';
      else {
         visu::Gmesh gmsh(DIM);
         std::vector<visu::Data> data; //create data

         size_t cell_i(0);
         size_t nb_nodes(0);
         for (auto& cl : m_msh)
         {
            vector_dynamic x = m_postprocess_data.at(cell_i++);
            auto cell_nodes = points(m_msh, cl);
            std::vector<visu::Node> new_nodes;
            for (auto& pt : cell_nodes)
            {
               nb_nodes++;
               auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt);

               std::array<double, 3> coor = {double{0.0}, double{0.0}, double{0.0}};

               visu::init_coordinate(pt, coor);
               visu::Node tmp_node(coor, nb_nodes, 0);
               new_nodes.push_back(tmp_node);
               gmsh.addNode(tmp_node);

               // plot magnitude at node
               scalar_type pot(0.0);
               for (size_t j = 0; j < m_bqd.cell_basis.range(0, m_cell_degree).size(); j++) //compute solution at the node
                  pot += phi[j] * x(j);

               std::vector<double> value(3, 0.0);
               value[DIM+1 -1] = {double(pot)};
               visu::Data datatmp(nb_nodes, value);

               data.push_back(datatmp);
            }
            // add new element
            visu::add_element(gmsh, new_nodes);
         }
         std::vector<visu::SubData> subdata; //not used
         visu::NodeData nodedata(3, 0.0, "sol_node_deformed", data, subdata); // create and init a nodedata view

         nodedata.saveNodeData(filename, gmsh); // save the view
      }
   }

   void
   saveMesh(const std::string& filename)
   {
      visu::Gmesh gmsh = visu::convertMesh(m_msh);
      gmsh.writeGmesh(filename, 2);
   }
};
