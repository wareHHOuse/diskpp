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
#include "mechanics/BoundaryConditions.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "nonlinear_elasticity_solver.hpp"

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_nl_elasticity_solver(const Mesh<T, 2, Storage>&        msh,
                         const ParamRun<T>&                rp,
                         const NLE::MaterialParameters<T>& material_data)
{
   typedef Mesh<T, 2, Storage>                            mesh_type;
   typedef static_vector<T, 2>                            result_type;
   typedef static_matrix<T, 2, 2>                         result_grad_type;
   typedef disk::mechanics::BoundaryConditions<mesh_type> Bnd_type;

   auto load = [material_data](const point<T, 2>& p) -> result_type {
      const T lambda = material_data.lambda;
      const T mu     = material_data.mu;

      T fx = lambda * cos(M_PI * (p.x() + p.y())) -
             2.0 * mu *
               ((4 * lambda + 4) * sin(2 * M_PI * p.y()) * cos(2 * M_PI * p.x()) +
                sin(M_PI * p.x()) * sin(M_PI * p.y())) +
             2.0 * mu *
               (2.0 * lambda * sin(2 * M_PI * p.y()) + 2.0 * sin(2 * M_PI * p.y()) +
                0.5 * cos(M_PI * (p.x() + p.y())));
      T fy = lambda * cos(M_PI * (p.x() + p.y())) +
             2.0 * mu *
               ((4 * lambda + 4) * sin(2 * M_PI * p.x()) * cos(2 * M_PI * p.y()) -
                sin(M_PI * p.x()) * sin(M_PI * p.y())) -
             2.0 * mu *
               (2.0 * lambda * sin(2 * M_PI * p.x()) + 2.0 * sin(2 * M_PI * p.x()) -
                0.5 * cos(M_PI * (p.x() + p.y())));

      return -M_PI * M_PI / (lambda + 1) * result_type{0, 0};
   };

   auto solution = [material_data](const point<T, 2>& p) -> result_type {
      T fx = sin(2 * M_PI * p.y()) * (cos(2 * M_PI * p.x()) - 1) +
             1.0 / (1 + material_data.lambda) * sin(M_PI * p.x()) * sin(M_PI * p.y());
      T fy = -sin(2 * M_PI * p.x()) * (cos(2 * M_PI * p.y()) - 1) +
             1.0 / (1 + material_data.lambda) * sin(M_PI * p.x()) * sin(M_PI * p.y());

      return result_type{fx, fy};
   };

   auto gradient = [material_data](const point<T, 2>& p) -> result_grad_type {

      const T lambda = material_data.lambda;

      T g11 = -(2 * (lambda + 1) * sin(2 * M_PI * p.x()) * sin(2 * M_PI * p.y()) -
                sin(M_PI * p.y()) * cos(M_PI * p.x()));
      T g12 = (2 * lambda + 2) * (cos(2 * M_PI * p.x()) - 1) * cos(2 * M_PI * p.y()) +
              sin(M_PI * p.x()) * cos(M_PI * p.y());

      T g21 = (-2 * lambda + 2) * (cos(2 * M_PI * p.y()) - 1) * cos(2 * M_PI * p.x()) +
              sin(M_PI * p.y()) * cos(M_PI * p.x());

      T g22 = 2 * (lambda + 1) * sin(2 * M_PI * p.x()) * sin(2 * M_PI * p.y()) +
              sin(M_PI * p.x()) * cos(M_PI * p.y());

      result_grad_type g = result_grad_type::Zero();

      g(0, 0) = g11;
      g(0, 1) = g12;
      g(1, 0) = g21;
      g(1, 1) = g22;

      return M_PI / (lambda + 1) * g;
   };

   Bnd_type bnd(msh);
   // bnd.addDirichletEverywhere(solution);

   bnd.addDirichletBC(disk::mechanics::CLAMPED, 10, solution);

   auto neu = [material_data](const point<T, 2>& p) -> result_type {
      return result_type{0, 1.0 / 16.0};
   };

   bnd.addNeumannBC(disk::mechanics::NEUMANN, 6, neu);

   nl_elasticity_solver<mesh_type> nl(msh, bnd, rp, material_data);

   if (nl.verbose()) {
      std::cout << "Solving the problem ..." << '\n';
   }

   SolverInfo solve_info = nl.compute(load);

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
   }

   nl.compute_discontinuous_displacement("sol2D.msh");
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_nl_elasticity_solver(const Mesh<T, 3, Storage>&        msh,
                         const ParamRun<T>&                rp,
                         const NLE::MaterialParameters<T>& material_data)
{
   typedef Mesh<T, 3, Storage>                            mesh_type;
   typedef static_vector<T, 3>                            result_type;
   typedef static_matrix<T, 3, 3>                         result_grad_type;
   typedef disk::mechanics::BoundaryConditions<mesh_type> Bnd_type;

   auto load = [material_data](const point<T, 3>& p) -> result_type {
      const T lambda = material_data.lambda;
      const T mu     = material_data.mu;

      T fx = lambda * cos(M_PI * (p.x() + p.y())) -
             2.0 * mu *
               ((4 * lambda + 4) * sin(2 * M_PI * p.y()) * cos(2 * M_PI * p.x()) +
                sin(M_PI * p.x()) * sin(M_PI * p.y())) +
             2.0 * mu *
               (2.0 * lambda * sin(2 * M_PI * p.y()) + 2.0 * sin(2 * M_PI * p.y()) +
                0.5 * cos(M_PI * (p.x() + p.y())));
      T fy = lambda * cos(M_PI * (p.x() + p.y())) +
             2.0 * mu *
               ((4 * lambda + 4) * sin(2 * M_PI * p.x()) * cos(2 * M_PI * p.y()) -
                sin(M_PI * p.x()) * sin(M_PI * p.y())) -
             2.0 * mu *
               (2.0 * lambda * sin(2 * M_PI * p.x()) + 2.0 * sin(2 * M_PI * p.x()) -
                0.5 * cos(M_PI * (p.x() + p.y())));

      return -M_PI * M_PI / (lambda + 1) * result_type{fx, fy, 0};
   };

   auto solution = [material_data](const point<T, 3>& p) -> result_type {
      T fx = sin(2 * M_PI * p.y()) * (cos(2 * M_PI * p.x()) - 1) +
             1.0 / (1 + material_data.lambda) * sin(M_PI * p.x()) * sin(M_PI * p.y());
      T fy = -sin(2 * M_PI * p.x()) * (cos(2 * M_PI * p.y()) - 1) +
             1.0 / (1 + material_data.lambda) * sin(M_PI * p.x()) * sin(M_PI * p.y());

      return result_type{fx, fy, 0};
   };

   auto gradient = [material_data](const point<T, 3>& p) -> result_grad_type {

      const T lambda = material_data.lambda;

      T g11 = -(2 * (lambda + 1) * sin(2 * M_PI * p.x()) * sin(2 * M_PI * p.y()) -
                sin(M_PI * p.y()) * cos(M_PI * p.x()));
      T g12 = (2 * lambda + 2) * (cos(2 * M_PI * p.x()) - 1) * cos(2 * M_PI * p.y()) +
              sin(M_PI * p.x()) * cos(M_PI * p.y());

      T g21 = (-2 * lambda + 2) * (cos(2 * M_PI * p.y()) - 1) * cos(2 * M_PI * p.x()) +
              sin(M_PI * p.y()) * cos(M_PI * p.x());

      T g22 = 2 * (lambda + 1) * sin(2 * M_PI * p.x()) * sin(2 * M_PI * p.y()) +
              sin(M_PI * p.x()) * cos(M_PI * p.y());

      result_grad_type g = result_grad_type::Zero();

      g(0, 0) = g11;
      g(0, 1) = g12;
      g(1, 0) = g21;
      g(1, 1) = g22;

      return M_PI / (lambda + 1) * g;
   };

   Bnd_type bnd(msh);
   bnd.addDirichletEverywhere(solution);

   nl_elasticity_solver<mesh_type> nl(msh, bnd, rp, material_data);

   if (nl.verbose()) {
      std::cout << "Solving the problem ..." << '\n';
   }

   SolverInfo solve_info = nl.compute(load);

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

      nl.compute_discontinuous_displacement("sol3D.msh");
   }
}

int
main(int argc, char** argv)
{
   using RealType = double;

   char* mesh_filename = nullptr;

   ParamRun<RealType> rp;

   // Elasticity Parameters
   NLE::MaterialParameters<RealType> material_data;

   material_data.mu     = 0.375;
   material_data.lambda = 7.5 * 10E6;

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
      std::cout << "Error" << std::endl;
      return 0;
   }

   mesh_filename = argv[0];

   /* FVCA5 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$"))) {
      std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
      auto msh = disk::load_fvca5_2d_mesh<RealType>(mesh_filename);
      run_nl_elasticity_solver(msh, rp, material_data);
      return 0;
   }

   /* Netgen 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$"))) {
      std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
      auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);
      run_nl_elasticity_solver(msh, rp, material_data);
      return 0;
   }

   /* DiSk++ cartesian 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.quad$"))) {
      std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
      auto msh = disk::load_cartesian_2d_mesh<RealType>(mesh_filename);
      run_nl_elasticity_solver(msh, rp, material_data);
      return 0;
   }

   /* DiSk++ cartesian 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.quad2$"))) {
      std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
      auto msh = disk::load_cartesian_2d_mesh2<RealType>(mesh_filename);
      run_nl_elasticity_solver(msh, rp, material_data);
      return 0;
   }

   /* Medit 2d*/
   if (std::regex_match(mesh_filename, std::regex(".*\\.medit2d$"))) {
      std::cout << "Guessed mesh format: Medit format" << std::endl;
      auto msh = disk::load_medit_2d_mesh<RealType>(mesh_filename);
      run_nl_elasticity_solver(msh, rp, material_data);
      return 0;
   }

   /* Netgen 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$"))) {
      std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
      auto msh = disk::load_netgen_3d_mesh<RealType>(mesh_filename);
      run_nl_elasticity_solver(msh, rp, material_data);
      return 0;
   }

   /* DiSk++ cartesian 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.hex$"))) {
      std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
      auto msh = disk::load_cartesian_3d_mesh<RealType>(mesh_filename);
      run_nl_elasticity_solver(msh, rp, material_data);
      return 0;
   }

   /* FVCA6 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.msh$"))) {
      std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
      auto msh = disk::load_fvca6_3d_mesh<RealType>(mesh_filename);
      run_nl_elasticity_solver(msh, rp, material_data);
      return 0;
   }
}
