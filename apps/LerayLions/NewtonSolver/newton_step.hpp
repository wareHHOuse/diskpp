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

 // NewtonRaphson_step

#pragma once

#include <iostream>

#include <sstream>
#include <string>
#include <fstream>


#ifdef HAVE_SOLVER_WRAPPERS
#include "agmg/agmg.hpp"
#endif

#include "hho/hho.hpp"
#include "hho/hho_bq.hpp"
#include "../Informations.hpp"
#include "../LerayLions_elementary_computation.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>



template<typename BQData>
class NewtonRaphson_step_leraylions
{
   typedef typename BQData::mesh_type          mesh_type;
   typedef typename mesh_type::scalar_type     scalar_type;

   typedef dynamic_matrix<scalar_type>         matrix_dynamic;
   typedef dynamic_vector<scalar_type>         vector_dynamic;


   typedef disk::gradient_reconstruction_full_bq<BQData>          gradrec_type;
   typedef LerayLions::LerayLions<BQData>                         elem_type;
   typedef disk::diffusion_like_static_condensation_bq<BQData>    statcond_type;
   typedef disk::assembler_nl_bq<BQData>                          assembler_type;


   typedef typename assembler_type::sparse_matrix_type       sparse_matrix_type;
   typedef typename assembler_type::vector_type              vector_type;

    #ifdef HAVE_INTEL_MKL
        typedef Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver_type;
    #else
        typedef Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver_type;
    #endif

    sparse_matrix_type     m_system_matrix;
    vector_type            m_system_rhs, m_system_solution;

    solver_type             solver;

    const BQData&                               m_bqd;

    const mesh_type& m_msh;

    std::vector<vector_dynamic>        m_postprocess_data, m_solution_data;
    std::vector<vector_dynamic>         m_solution_cells, m_solution_faces, m_solution_lagr;

    bool m_verbose;

    scalar_type m_leray_param;

    scalar_type initial_residual;


public:
   NewtonRaphson_step_leraylions(const mesh_type& msh, const BQData& bqd, const scalar_type leray_param)
   : m_msh(msh), m_verbose(false), m_leray_param(leray_param), m_bqd(bqd)
    {}

    bool    verbose(void) const     { return m_verbose; }
    void    verbose(bool v)         { m_verbose = v; }



    void
    initialize( const std::vector<vector_dynamic>& initial_solution_cells,
                const std::vector<vector_dynamic>& initial_solution_faces,
                const std::vector<vector_dynamic>& initial_solution_lagr,
                const std::vector<vector_dynamic>& initial_solution)
    {
      m_solution_cells.clear();
      m_solution_cells = initial_solution_cells;
      assert(m_msh.cells_size() == m_solution_cells.size());

      m_solution_faces.clear();
      m_solution_faces = initial_solution_faces;
      assert(m_msh.faces_size() == m_solution_faces.size());

      m_solution_lagr.clear();
      m_solution_lagr = initial_solution_lagr;

      m_solution_data.clear();
      m_solution_data = initial_solution;
      assert(m_msh.cells_size() == m_solution_data.size());
    }

    template<typename LoadFunction, typename BoundaryConditionFunction>
    AssemblyInfo
    assemble(const LoadFunction& lf, const BoundaryConditionFunction& bcf)
    {
        gradrec_type    gradrec(m_bqd);
        elem_type       elem(m_bqd);
        statcond_type   statcond(m_bqd);
        assembler_type  assembler(m_msh, m_bqd);

        AssemblyInfo ai;

        timecounter tc, ttot;

        ttot.tic();
        size_t i = 0;

        for (auto& cl : m_msh)
        {
            tc.tic();
            gradrec.compute_optim(m_msh, cl, false);
            tc.toc();
            ai.m_time_gradrec += tc.to_double();

            tc.tic();
            elem.compute(m_msh, cl, lf, gradrec.oper(), m_solution_data.at(i), m_leray_param);


            /////// NON LINEAIRE /////////
            dynamic_matrix<scalar_type> lhs = elem.K_int;
            dynamic_vector<scalar_type> rhs = elem.RTF;

            tc.toc();
            ai.m_time_elem += tc.to_double();
            ai.m_time_law += elem.time_law;
            tc.tic();
            auto scnp = statcond.compute(m_msh, cl, lhs, rhs, true);
            tc.toc();
            ai.m_time_statcond += tc.to_double();

            assembler.assemble(m_msh, cl, scnp);

            i++;
        }

         assembler.impose_boundary_conditions(m_msh, bcf, m_solution_faces, m_solution_lagr);
         assembler.finalize(m_system_matrix, m_system_rhs);

         ttot.toc();
         ai.m_time_assembly = ttot.to_double();
         ai.m_linear_system_size = m_system_matrix.rows();
        return ai;
    }


    void
    saveMatrix(const std::string& filename)
    {
       std::ofstream fichier(filename, std::ios::out | std::ios::trunc);

       for (int k=0; k<m_system_matrix.outerSize(); ++k)
          for (typename Eigen::SparseMatrix<scalar_type>::InnerIterator it(m_system_matrix,k); it; ++it)
         {
            fichier << it.row() << " ; " << it.col() << " ; " << it.value() << std::endl;
         }

       fichier.close();
    }

    size_t
    test_coercivity(void)
    {

      Eigen::SelfAdjointEigenSolver<sparse_matrix_type> es;
      es.compute(m_system_matrix);

      if(es.info() != Eigen::Success) {
          std::cerr << "ERROR: Could not compute eigenvalues of the matrix" << std::endl;
      }

      auto ev = es.eigenvalues();

      size_t nb_negative_eigenvalue(0);

      size_t i_ev = 0;
      while(ev(i_ev) <= scalar_type(0.0))
      {
         nb_negative_eigenvalue++;
         i_ev++;
      }

      const scalar_type min_ev = ev.minCoeff();
      const scalar_type max_ev = ev.maxCoeff();

      if(m_verbose){
         std::cout << "******* Eigenvalues test ********" << std::endl;
         std::cout << "Number of eigenvalues: " << ev.size()  << std::endl;
         std::cout << "Number of negative eigenvalues: " << nb_negative_eigenvalue  << std::endl;
         std::cout << "Maximum eigenvalue: " << max_ev << std::endl;
         std::cout << "Minimum eigenvalue: " << min_ev << std::endl;
         std::cout << "Conditionning number: " << std::abs(max_ev/min_ev) << std::endl;
      }

      return nb_negative_eigenvalue;
    }


    SolveInfo
    solve(void)
    {
       #ifdef HAVE_INTEL_MKL
       Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
       #else
       Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
       #endif

       timecounter tc;

       tc.tic();
       solver.analyzePattern(m_system_matrix);
       solver.factorize(m_system_matrix);

       if(solver.info() != Eigen::Success) {
          std::cerr << "ERROR: Could not factorize the matrix" << std::endl;
       }

       m_system_solution = solver.solve(m_system_rhs);
       if(solver.info() != Eigen::Success) {
          std::cerr << "ERROR: Could not solve the linear system" << std::endl;
       }
       tc.toc();

       return SolveInfo(m_system_matrix.rows(), m_system_matrix.nonZeros(), tc.to_double());
    }

    //a améliorer supprimer condensation statique
    template<typename LoadFunction, typename BoundaryConditionFunction>
    bool
    line_search(const LoadFunction& lf, const BoundaryConditionFunction& bcf,
                scalar_type& gamma, const size_t max_iter)
    {

       gamma = scalar_type(1.0);
       const scalar_type r0 = m_system_rhs.norm();



       for(size_t iter = 0; iter < max_iter; iter++){


         std::vector<vector_dynamic>        solution_data(m_solution_data);
         std::vector<vector_dynamic>        solution_faces(m_solution_faces);
         std::vector<vector_dynamic>        solution_lagr(m_solution_lagr);

         for(size_t i=0; i < m_postprocess_data.size(); i++){
            solution_data.at(i) += gamma * m_postprocess_data.at(i);
         }
         size_t fbs = m_bqd.face_basis.size();

         for(size_t i=0; i < m_solution_faces.size(); i++){
            solution_faces.at(i) += gamma * m_system_solution.block(i * fbs, 0, fbs, 1);
         }
         const size_t lagrange_offset = m_solution_faces.size() * fbs;
         for(size_t i=0; i < m_solution_lagr.size(); i++){
            solution_lagr.at(i) += gamma * m_system_solution.block(lagrange_offset + i * fbs, 0, fbs, 1);
         }


         //compute rhs
         gradrec_type    gradrec(m_bqd);
         elem_type       elem(m_bqd);
         statcond_type   statcond(m_bqd);
         assembler_type  assembler(m_msh, m_bqd);
         sparse_matrix_type     system_matrix;
         vector_type            system_rhs;

         AssemblyInfo ai;

         timecounter tc, ttot;

         ttot.tic();
         size_t i = 0;




         for (auto& cl : m_msh)
         {
            tc.tic();
            gradrec.compute_optim(m_msh, cl);
            tc.toc();
            ai.m_time_gradrec += tc.to_double();

            tc.tic();
            elem.compute(m_msh, cl, lf, gradrec.oper, solution_data.at(i), m_leray_param);


            /////// NON LINEAIRE /////////
            dynamic_matrix<scalar_type> lhs = elem.K_int;
            dynamic_vector<scalar_type> rhs = elem.RTF;

            tc.toc();
            ai.m_time_elem += tc.to_double();
            ai.m_time_law += elem.time_law;
            tc.tic();
            auto scnp = statcond.compute(m_msh, cl, lhs, rhs, true);
            tc.toc();
            ai.m_time_statcond += tc.to_double();

            assembler.assemble(m_msh, cl, scnp);

            i++;
         }

         assembler.impose_boundary_conditions(m_msh, bcf, solution_faces, solution_lagr);
         assembler.finalize(system_matrix, system_rhs);

         ttot.toc();
         ai.m_time_assembly = ttot.to_double();
         ai.m_linear_system_size = system_matrix.rows();

         const scalar_type r = system_rhs.norm();

         std::cout << "r0= " << r0 << " , r= " << r << std::endl;

         if(std::abs(r) < 0.1*std::abs(r0)){
            m_system_matrix = system_matrix;
            m_system_rhs = system_rhs;
            return true;
         }
         else
            gamma /=2;
       }//fin for
       gamma = 1.0;
       return false;
    }



    template<typename LoadFunction>
    PostprocessInfo
    postprocess(const LoadFunction& lf)
    {
        gradrec_type gradrec(m_bqd);
        elem_type elem(m_bqd);
        statcond_type statcond(m_bqd);

        const size_t fbs = m_bqd.face_basis.size();
        const size_t cbs = (m_bqd.cell_basis.range(0, m_bqd.cell_degree())).size();

        PostprocessInfo pi;

        m_postprocess_data.clear();

        m_postprocess_data.reserve(m_msh.cells_size());

        timecounter tc, ttot;
        tc.tic(); ttot.tic();

        size_t i = 0;
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

            tc.tic();
            gradrec.compute_optim(m_msh, cl, false);
            tc.toc();
            pi.m_time_gradrec += tc.to_double();

            tc.tic();
            elem.compute(m_msh, cl, lf, gradrec.oper(), m_solution_data.at(i), m_leray_param);

            /////// NON LINEAIRE /////////
            dynamic_matrix<scalar_type> lhs = elem.K_int;
            dynamic_vector<scalar_type> rhs = elem.RTF;
            dynamic_vector<scalar_type> rhs_cell = rhs.block(0,0, cbs, 1);

            tc.toc();
            pi.m_time_elem += tc.to_double();
            pi.m_time_law += elem.time_law;

            tc.tic();
            dynamic_vector<scalar_type> x = statcond.recover(m_msh, cl, lhs, rhs_cell, xFs);
            tc.toc();
            pi.m_time_statcond += tc.to_double();

            m_postprocess_data.push_back(x);

            i++;
        }
        ttot.toc();

        pi.m_time_post = ttot.to_double();

        return pi;
    }



    void
    update_solution(const scalar_type gamma = 1.0)
    {
        assert(m_postprocess_data.size() == m_solution_data.size());

        for(size_t i=0; i < m_postprocess_data.size(); i++){
            assert(m_solution_data.at(i).size() == m_postprocess_data.at(i).size());
            m_solution_data.at(i) += gamma * m_postprocess_data.at(i);
        }

        size_t cbs = (m_bqd.cell_basis.range(0, m_bqd.cell_degree())).size();

        assert(m_postprocess_data.size() == m_solution_cells.size());

        for(size_t i=0; i < m_postprocess_data.size(); i++){
            assert(m_solution_cells.at(i).size() == cbs);
            m_solution_cells.at(i) += gamma * (m_postprocess_data.at(i)).block(0,0,cbs,1);
        }

        size_t fbs = m_bqd.face_basis.size();


        for(size_t i=0; i < m_solution_faces.size(); i++){
            assert(m_solution_faces.at(i).size() == fbs);
            m_solution_faces.at(i) += gamma * m_system_solution.block(i * fbs, 0, fbs, 1);
        }

        const size_t lagrange_offset = m_solution_faces.size() * fbs;
        for(size_t i=0; i < m_solution_lagr.size(); i++){
            assert(m_solution_lagr.at(i).size() == fbs);
            m_solution_lagr.at(i) += gamma * m_system_solution.block(lagrange_offset + i * fbs, 0, fbs, 1);
        }
    }



    bool test_convergence(const scalar_type epsilon, const size_t iter, scalar_type& error)
    {
      if(iter == 0)
         initial_residual = m_system_rhs.norm();

      scalar_type residual = m_system_rhs.norm();
      scalar_type max_error(0.0);
      for(size_t i = 0; i < m_system_rhs.size(); i++)
         max_error = std::max( max_error, std::abs(m_system_rhs(i)));

      scalar_type relative_error = residual / initial_residual;

      if(initial_residual == scalar_type{0.0})
      {
         relative_error = 0.0;
         max_error = 0.0;
      }

      const scalar_type error_incr = m_system_solution.norm();

      scalar_type error_un(0.0);

      for (size_t i = 0; i < m_solution_faces.size(); i++) {
         scalar_type norm = m_solution_faces[i].norm();
         error_un += norm * norm;
      }

      for (size_t i = 0; i < m_solution_lagr.size(); i++) {
         scalar_type norm = m_solution_lagr[i].norm();
         error_un += norm * norm;
      }

      error_un = sqrt(error_un);

      if(error_un == scalar_type(0.0))
         error_un = 10E-6;

      scalar_type relative_displ = error_incr/error_un;

      if(iter == 0)
         relative_displ = 1.0;

      size_t nb_faces_dof = m_bqd.face_basis.size() * m_msh.faces_size();

      scalar_type norm_depl = (m_system_rhs.head(nb_faces_dof)).norm();
      scalar_type norm_lag = (m_system_rhs.tail(m_system_rhs.size() - nb_faces_dof)).norm();

      if(m_verbose){
         std::string s_iter = "   " + std::to_string(iter) + "               ";
         s_iter.resize(9);

         if(iter == 0){
           std::cout << "--------------------------------------------------------------------------------------------------------------------------------" << std::endl;
           std::cout << "| Iteration | Norme l2 incr | Relative depl |  Residual l2  | Relative error | Maximum error | Residual face  |  Residual BC   |" << std::endl;
           std::cout << "--------------------------------------------------------------------------------------------------------------------------------" << std::endl;

         }
         std::ios::fmtflags f( std::cout.flags() );
         std::cout.precision(5);
         std::cout.setf(std::iostream::scientific, std::iostream::floatfield);
         std::cout << "| " << s_iter << " |   " << error_incr << " |   " << relative_displ << " |   " << residual << " |   " << relative_error << "  |  " << max_error << "  |   "
         << norm_depl << "  |   " << norm_lag << "  |" << std::endl;
         std::cout << "--------------------------------------------------------------------------------------------------------------------------------" << std::endl;
         std::cout.flags( f );
      }

      error = std::min(relative_displ, relative_error);

      if(error <= epsilon ){
         return true;
      }
      else {
         return false;
      }

   }

   void
    save_solutions( std::vector<vector_dynamic>& solution_cells,
                    std::vector<vector_dynamic>& solution_faces,
                    std::vector<vector_dynamic>& solution_lagr,
                    std::vector<vector_dynamic>& solution)
    {
      solution_cells.clear();
      solution_cells = m_solution_cells;
      assert(m_solution_cells.size() == solution_cells.size());

      solution_faces.clear();
      solution_faces = m_solution_faces;
      assert(m_solution_faces.size() == solution_faces.size());

      solution_lagr.clear();
      solution_lagr = m_solution_lagr;
      assert(m_solution_lagr.size() == solution_lagr.size());

      solution.clear();
      solution = m_solution_data;
      assert(m_solution_data.size() == solution.size());
    }
};
