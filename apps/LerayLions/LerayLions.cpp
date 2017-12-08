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

#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <unistd.h>

#include <map>

#include "Informations.hpp"
#include "Parameters.hpp"
#include "loaders/loader.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "LerayLions_solver.hpp"
/*
template<template<typename, size_t , typename> class Mesh,
         typename T, typename Storage>
void
run_leraylions_solver(const Mesh<T, 1, Storage>& msh, ParamRun<T>& rp, const T leray_param)
{
   typedef T result_type;

   auto load = [leray_param](const point<T,1>& pt) -> result_type {
      const T p = leray_param;
      result_type norm_G = M_PI *  std::abs(cos(pt.x() * M_PI));

      T fx = 0;

      return result_type{fx};
   };

   auto solution = [leray_param](const point<T,1>& pt) -> result_type {
      return sin(pt.x() * M_PI);
   };


   leraylions_solver<Mesh, T, 1, Storage,  point<T, 1> > nl(msh, rp, leray_param);


   nl.compute_initial_state();

   if(nl.verbose()){
      std::cout << "Solving the problem ..."  << '\n';
   }

   SolverInfo solve_info = nl.compute(load, solution);

   if(nl.verbose()){
      std::cout << " " << std::endl;
      std::cout << "------------------------------------------------------- " << std::endl;
      std::cout << "Summaring: " << std::endl;
      std::cout << "Total time to solve the problem: " << solve_info.m_time_solver << " sec" <<
std::endl; std::cout << "**** Assembly time: " <<
solve_info.m_newton_info.m_assembly_info.m_time_assembly << " sec" << std::endl; std::cout <<
"****** Gradient reconstruction: " << solve_info.m_newton_info.m_assembly_info.m_time_gradrec << "
sec" << std::endl; std::cout << "****** Elementary computation: " <<
solve_info.m_newton_info.m_assembly_info.m_time_elem << " sec" << std::endl; std::cout << "
*** Behavior computation: " << solve_info.m_newton_info.m_assembly_info.m_time_law << " sec" <<
std::endl; std::cout << "****** Static condensation: " <<
solve_info.m_newton_info.m_assembly_info.m_time_statcond << " sec" << std::endl; std::cout << "****
Solver time: " << solve_info.m_newton_info.m_solve_info.m_time_solve << " sec" << std::endl;
      std::cout << "**** Postprocess time: " <<
solve_info.m_newton_info.m_postprocess_info.m_time_post << " sec" << std::endl; std::cout <<
"------------------------------------------------------- " << std::endl; std::cout << " " <<
std::endl;
   }


   if(nl.test_convergence()){
      std::cout << "average diameter h: " << average_diameter(msh) << std::endl;
      std::cout << "l2 error displacement: " << nl.compute_l2_error(solution) << std::endl;

        std::cout << "Post-processing: " << std::endl;
         nl.compute_discontinuous_solution("sol_disc1D.msh");
//         nl.compute_conforme_solution("sol_conf2D.msh");
         nl.plot_solution_at_gausspoint("depl_gp1D.msh");
         nl.compute_deformed("def2d.msh");
   }
}*/

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_leraylions_solver(const Mesh<T, 2, Storage>& msh, const ParamRun<T>& rp)
{
   typedef T result_type;

   auto load = [rp](const point<T, 2>& pt) -> result_type {
      const T     p      = rp.m_leray_param;
      result_type norm_G = M_PI * sqrt(std::pow(cos(pt.x() * M_PI) * sin(pt.y() * M_PI), 2) +
                                       std::pow(sin(pt.x() * M_PI) * cos(pt.y() * M_PI), 2));

      T fx = sin(M_PI * pt.x()) * sin(M_PI * pt.y()) *
             (2. * std::pow(M_PI, 2) * std::pow(norm_G, p - 2.0) -
              (p - 2) * std::pow(M_PI, 4) * std::pow(norm_G, p - 4.0) *
                (std::pow(cos(M_PI * pt.x()), 2) * cos(2. * M_PI * pt.y()) +
                 cos(2. * M_PI * pt.x()) * std::pow(cos(M_PI * pt.y()), 2)));

      return result_type{fx};
   };

   auto solution = [rp](const point<T, 2>& pt) -> result_type {
      return sin(pt.x() * M_PI) * sin(pt.y() * M_PI);
   };

   leraylions_solver<Mesh, T, 2, Storage> nl(msh, rp);

   if (nl.verbose()) {
      std::cout << "Solving the problem ..." << '\n';
   }

   SolverInfo solve_info = nl.compute(load, solution);

   if (nl.verbose()) {
      std::cout << " " << std::endl;
      std::cout << "------------------------------------------------------- " << std::endl;
      std::cout << "Summaring: " << std::endl;
      std::cout << "Total Newton's iterations: " << solve_info.m_iter << " in "
                << solve_info.m_time_step << " load increments" << std::endl;
      std::cout << "Total time to solve the problem: " << solve_info.m_time_solver << " sec"
                << std::endl;
      std::cout << "**** Assembly time: "
                << solve_info.m_newton_info.m_assembly_info.m_time_assembly << " sec" << std::endl;
      std::cout << "****** Gradient reconstruction: "
                << solve_info.m_newton_info.m_assembly_info.m_time_gradrec << " sec" << std::endl;
      std::cout << "****** Stabilisation: " << solve_info.m_newton_info.m_assembly_info.m_time_stab
                << " sec" << std::endl;
      std::cout << "****** Elementary computation: "
                << solve_info.m_newton_info.m_assembly_info.m_time_elem << " sec" << std::endl;
      std::cout << "       *** Behavior computation: "
                << solve_info.m_newton_info.m_assembly_info.m_time_law << " sec" << std::endl;
      std::cout << "****** Static condensation: "
                << solve_info.m_newton_info.m_assembly_info.m_time_statcond << " sec" << std::endl;
      std::cout << "**** Postprocess time: "
                << solve_info.m_newton_info.m_assembly_info.m_time_postpro << " sec" << std::endl;
      std::cout << "**** Solver time: " << solve_info.m_newton_info.m_solve_info.m_time_solve
                << " sec" << std::endl;
      std::cout << "------------------------------------------------------- " << std::endl;
      std::cout << " " << std::endl;
   }

   if (nl.test_convergence()) {
      std::cout << "average diameter h: " << average_diameter(msh) << std::endl;
      std::cout << "l2 error displacement: " << nl.compute_l2_error(solution) << std::endl;

      std::cout << "Post-processing: " << std::endl;
      //   nl.compute_discontinuous_solution("sol_disc2D.msh");
      //   //         nl.compute_conforme_solution("sol_conf2D.msh");
      //   nl.plot_solution_at_gausspoint("depl_gp2D.msh");
      //   nl.compute_deformed("def2d.msh");
   }
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_leraylions_solver(const Mesh<T, 3, Storage>& msh, const ParamRun<T>& rp)
{
   typedef T result_type;

   auto load = [rp](const point<T, 3>& pt) -> auto
   {
      const T p  = rp.m_leray_param;
      T       fx = -(std::pow(3.0, p / 2.0) * (p - 1) * exp((p - 1) * (pt.x() + pt.y() + pt.z())));
      return result_type{fx};
   };

   auto solution = [](const point<T, 3>& pt) -> auto
   {
      T fx = exp(pt.x() + pt.y() + pt.z());

      return result_type{fx};
   };

   leraylions_solver<Mesh, T, 3, Storage> nl(msh, rp);

   if (nl.verbose()) {
      std::cout << "Solving the problem ..." << '\n';
   }

   SolverInfo solve_info = nl.compute(load, solution);

   if (nl.verbose()) {
      std::cout << " " << std::endl;
      std::cout << "------------------------------------------------------- " << std::endl;
      std::cout << "Summaring: " << std::endl;
      std::cout << "Total Newton's iterations: " << solve_info.m_iter << " in "
                << solve_info.m_time_step << " load increments" << std::endl;
      std::cout << "Total time to solve the problem: " << solve_info.m_time_solver << " sec"
                << std::endl;
      std::cout << "**** Assembly time: "
                << solve_info.m_newton_info.m_assembly_info.m_time_assembly << " sec" << std::endl;
      std::cout << "****** Gradient reconstruction: "
                << solve_info.m_newton_info.m_assembly_info.m_time_gradrec << " sec" << std::endl;
      std::cout << "****** Stabilisation: " << solve_info.m_newton_info.m_assembly_info.m_time_stab
                << " sec" << std::endl;
      std::cout << "****** Elementary computation: "
                << solve_info.m_newton_info.m_assembly_info.m_time_elem << " sec" << std::endl;
      std::cout << "       *** Behavior computation: "
                << solve_info.m_newton_info.m_assembly_info.m_time_law << " sec" << std::endl;
      std::cout << "****** Static condensation: "
                << solve_info.m_newton_info.m_assembly_info.m_time_statcond << " sec" << std::endl;
      std::cout << "**** Postprocess time: "
                << solve_info.m_newton_info.m_assembly_info.m_time_postpro << " sec" << std::endl;
      std::cout << "**** Solver time: " << solve_info.m_newton_info.m_solve_info.m_time_solve
                << " sec" << std::endl;
      std::cout << "------------------------------------------------------- " << std::endl;
      std::cout << " " << std::endl;
   }

   if (nl.test_convergence()) {
      std::cout << "avrage diameter h: " << average_diameter(msh) << std::endl;
      std::cout << "l2 error: " << nl.compute_l2_error(solution) << std::endl;

      std::cout << "Post-processing: " << std::endl;
      //         nl.compute_discontinuous_solution("sol_disc3D.msh");
      //         nl.plot_solution_at_gausspoint("sol_gp3D.msh");
   }
}

int
main(int argc, char** argv)
{
   using RealType = double;

   char* mesh_filename = nullptr;
   int   elems_1d      = 8;

   ParamRun<RealType> rp;

   int ch;

   while ((ch = getopt(argc, argv, "r:")) != -1) {
      switch (ch) {
         case 'r':
            if (!rp.readParameters(optarg)) exit(1);
            break;
         default: std::cout << "wrong arguments" << std::endl; exit(1);
      }
   }

   argc -= optind;
   argv += optind;

   if (argc == 0) {
      std::cout << "Mesh format: 1D uniform" << std::endl;
      // auto msh = disk::load_uniform_1d_mesh<RealType>(0, 1, elems_1d);
      // run_leraylions_solver(msh, rp, leray_param);
      return 0;
   }

   mesh_filename = argv[0];

   /* FVCA5 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$"))) {
      std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
      auto msh = disk::load_fvca5_2d_mesh<RealType>(mesh_filename);
      run_leraylions_solver(msh, rp);
      return 0;
   }

   /* Netgen 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$"))) {
      std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
      auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);
      run_leraylions_solver(msh, rp);
      return 0;
   }

   /* DiSk++ cartesian 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.quad$"))) {
      std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
      auto msh = disk::load_cartesian_2d_mesh<RealType>(mesh_filename);
      run_leraylions_solver(msh, rp);
      return 0;
   }

   /* DiSk++ cartesian 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.quad2$"))) {
      std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
      auto msh = disk::load_cartesian_2d_mesh2<RealType>(mesh_filename);
      run_leraylions_solver(msh, rp);
      return 0;
   }

   /* Medit 2d*/
   if (std::regex_match(mesh_filename, std::regex(".*\\.medit2d$"))) {
      std::cout << "Guessed mesh format: Medit format" << std::endl;
      auto msh = disk::load_medit_2d_mesh<RealType>(mesh_filename);
      return 0;
   }

   /* Netgen 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$"))) {
      std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
      auto msh = disk::load_netgen_3d_mesh<RealType>(mesh_filename);
      run_leraylions_solver(msh, rp);
      return 0;
   }

   /* DiSk++ cartesian 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.hex$"))) {
      std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
      auto msh = disk::load_cartesian_3d_mesh<RealType>(mesh_filename);
      run_leraylions_solver(msh, rp);
      return 0;
   }

   /* FVCA6 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.msh$"))) {
      std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
      auto msh = disk::load_fvca6_3d_mesh<RealType>(mesh_filename);
      run_leraylions_solver(msh, rp);
      return 0;
   }
}
