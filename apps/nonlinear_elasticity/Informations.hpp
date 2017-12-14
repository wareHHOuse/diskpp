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

#pragma once

class AssemblyInfo
{
 public:
   size_t m_linear_system_size;
   double m_time_assembly;
   double m_time_gradrec, m_time_statcond, m_time_stab, m_time_elem, m_time_law, m_time_postpro;

   AssemblyInfo() :
     m_linear_system_size(0), m_time_assembly(0.0), m_time_gradrec(0.0), m_time_statcond(0.0),
     m_time_stab(0.0), m_time_elem(0.0), m_time_law(0.0), m_time_postpro(0.0)
   {}
};

class SolveInfo
{
 public:
   size_t m_linear_system_size, m_nonzeros;
   double m_time_solve;

   SolveInfo() : m_linear_system_size(0), m_nonzeros(0), m_time_solve(0.0) {}
   SolveInfo(const size_t linear_system_size, const size_t nonzeros, const double time_solve) :
     m_linear_system_size(linear_system_size), m_nonzeros(nonzeros), m_time_solve(time_solve)
   {}
};

class NewtonSolverInfo
{
 public:
   AssemblyInfo m_assembly_info;
   SolveInfo    m_solve_info;
   double       m_time_newton;
   size_t       m_iter;

   NewtonSolverInfo() : m_assembly_info(), m_solve_info(), m_time_newton(0.0), m_iter(0) {}

   void
   updateAssemblyInfo(const AssemblyInfo& assembly_info)
   {
      m_assembly_info.m_linear_system_size = assembly_info.m_linear_system_size;
      m_assembly_info.m_time_assembly += assembly_info.m_time_assembly;
      m_assembly_info.m_time_gradrec += assembly_info.m_time_gradrec;
      m_assembly_info.m_time_statcond += assembly_info.m_time_statcond;
      m_assembly_info.m_time_stab += assembly_info.m_time_stab;
      m_assembly_info.m_time_elem += assembly_info.m_time_elem;
      m_assembly_info.m_time_law += assembly_info.m_time_law;
      m_assembly_info.m_time_postpro += assembly_info.m_time_postpro;
   }

   void
   updateSolveInfo(const SolveInfo& solve_info)
   {
      m_solve_info.m_linear_system_size = solve_info.m_linear_system_size;
      m_solve_info.m_nonzeros           = solve_info.m_nonzeros;
      m_solve_info.m_time_solve += solve_info.m_time_solve;
   }
};

class SolverInfo
{
 public:
   NewtonSolverInfo m_newton_info;
   double           m_time_solver;
   size_t           m_iter, m_time_step;

   SolverInfo() : m_newton_info(), m_time_solver(0.0), m_iter(0), m_time_step(0) {}

   void
   updateInfo(const NewtonSolverInfo& newton_solver_info)
   {
      m_iter += newton_solver_info.m_iter;
      m_newton_info.updateAssemblyInfo(newton_solver_info.m_assembly_info);
      m_newton_info.updateSolveInfo(newton_solver_info.m_solve_info);
      m_newton_info.m_time_newton += newton_solver_info.m_time_newton;
   }
};