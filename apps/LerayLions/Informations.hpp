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
 * cite
*/

#pragma once

class AssemblyInfo
{
public:
   size_t  m_linear_system_size;
   double  m_time_assembly;
   double  m_time_gradrec, m_time_statcond, m_time_elem, m_time_law;


   AssemblyInfo() :  m_linear_system_size(0), m_time_assembly(0.0),
   m_time_gradrec(0.0), m_time_statcond(0.0), m_time_elem(0.0), m_time_law(0.0) {}
};


class SolveInfo
{
public:
   size_t  m_linear_system_size, m_nonzeros;
   double  m_time_solve;


   SolveInfo() :  m_linear_system_size(0), m_nonzeros(0), m_time_solve(0.0) {}
   SolveInfo(const size_t linear_system_size, const size_t nonzeros, const double time_solve) :
   m_linear_system_size(linear_system_size), m_nonzeros(nonzeros), m_time_solve(time_solve) {}

};


class PostprocessInfo
{
public:
   double  m_time_post;
   double  m_time_gradrec, m_time_statcond, m_time_elem, m_time_law;


   PostprocessInfo() :  m_time_post(0.0),
   m_time_gradrec(0.0), m_time_statcond(0.0), m_time_elem(0.0), m_time_law(0.0){}
};



class NewtonSolverInfo
{
public:
   AssemblyInfo     m_assembly_info;
   SolveInfo        m_solve_info;
   PostprocessInfo  m_postprocess_info;
   double           m_time_newton;


   NewtonSolverInfo() :  m_assembly_info(), m_solve_info(), m_postprocess_info(), m_time_newton(0.0) {}

   void updateAssemblyInfo(const AssemblyInfo&  assembly_info)
   {
      m_assembly_info.m_linear_system_size = assembly_info.m_linear_system_size;
      m_assembly_info.m_time_assembly += assembly_info.m_time_assembly;
      m_assembly_info.m_time_gradrec += assembly_info.m_time_gradrec;
      m_assembly_info.m_time_statcond += assembly_info.m_time_statcond;
      m_assembly_info.m_time_elem += assembly_info.m_time_elem;
      m_assembly_info.m_time_law += assembly_info.m_time_law;
   }

   void updateSolveInfo(const SolveInfo&  solve_info)
   {
      m_solve_info.m_linear_system_size = solve_info.m_linear_system_size;
      m_solve_info.m_nonzeros = solve_info.m_nonzeros;
      m_solve_info.m_time_solve += solve_info.m_time_solve;
   }

   void updatePostProcessInfo(const PostprocessInfo&  postprocess_info)
   {
      m_postprocess_info.m_time_post += postprocess_info.m_time_post;
      m_postprocess_info.m_time_gradrec += postprocess_info.m_time_gradrec;
      m_postprocess_info.m_time_statcond += postprocess_info.m_time_statcond;
      m_postprocess_info.m_time_elem += postprocess_info.m_time_elem;
      m_postprocess_info.m_time_law += postprocess_info.m_time_law;
   }
};


class SolverInfo
{
public:
   NewtonSolverInfo  m_newton_info;
   double  m_time_solver;


   SolverInfo() :  m_newton_info(), m_time_solver(0.0) {}

   void updateInfo(const NewtonSolverInfo& newton_solver_info)
   {
      m_newton_info.updateAssemblyInfo(newton_solver_info.m_assembly_info);
      m_newton_info.updateSolveInfo(newton_solver_info.m_solve_info);
      m_newton_info.updatePostProcessInfo(newton_solver_info.m_postprocess_info);
      m_newton_info.m_time_newton += newton_solver_info.m_time_newton;
   }
};
