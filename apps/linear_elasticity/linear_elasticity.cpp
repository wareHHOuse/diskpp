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

#include "geometry/geometry.hpp"
#include "loaders/loader.hpp"
#include "mechanics/BoundaryConditions.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "linear_elasticity_solver.hpp"

struct run_params
{
   size_t degree;
   int    l;
   bool   verbose;
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_linear_elasticity_solver(const Mesh<T, 2, Storage>& msh,
                             const run_params&          rp,
                             const ElasticityParameters material_data)
{
   typedef Mesh<T, 2, Storage> mesh_type;
   typedef static_vector<T, 2> result_type;

   timecounter tc;
   tc.tic();

   auto load = [material_data](const point<T, 2>& p) -> result_type {
      T fx = 2. * material_data.mu * M_PI * M_PI * sin(M_PI * p.x()) * sin(M_PI * p.y());
      T fy = 2. * material_data.mu * M_PI * M_PI * cos(M_PI * p.x()) * cos(M_PI * p.y());

      return result_type{fx, fy};
   };

   auto solution = [material_data](const point<T, 2>& p) -> result_type {
      T fx = sin(M_PI * p.x()) * sin(M_PI * p.y()) + 1.0 / (2.0 * material_data.lambda) * p.x();
      T fy = cos(M_PI * p.x()) * cos(M_PI * p.y()) + 1.0 / (2.0 * material_data.lambda) * p.y();

      return result_type{fx, fy};
   };

   typedef BoundaryConditions<mesh_type, decltype(solution)> Bnd_type;
   Bnd_type                                                  bnd(msh);

   bnd.addDirichletEverywhere(solution);

   //    auto b1 = [material_data](const point<T, 2>& p) -> result_type { return result_type{0, 0};
   //    };

   //    auto b2 = [material_data](const point<T, 2>& p) -> result_type { return result_type{0, 0};
   //    };

   //    auto b3 = [material_data](const point<T, 2>& p) -> result_type { return result_type{1, 1};
   //    };

   //    bnd.addDirichletBC(DX, 1, solution);
   //    bnd.addDirichletBC(DY, 2, solution);
   //    bnd.addDirichletBC(DXDY, 3, solution);

   linear_elasticity_solver<mesh_type, Bnd_type> le(msh, bnd, material_data, rp.degree);
   le.verbose(rp.verbose);

   le.changeElasticityParameters(material_data);

   assembly_info assembling_info = le.assemble(load);

   solver_info solve_info = le.solve();

   postprocess_info post_info = le.postprocess(load);

   tc.toc();

   if (le.verbose()) {
      std::cout << std::endl;
      std::cout << "************************************************************" << std::endl;
      std::cout << "** Time to solve the problem " << tc.to_double() << " sec" << std::endl;
      std::cout << "**** Assembly time: " << assembling_info.time_assembly << " sec" << std::endl;
      std::cout << "****** Gradient reconstruction: " << assembling_info.time_gradrec << " sec"
                << std::endl;
      std::cout << "****** Divergence reconstruction: " << assembling_info.time_divrec << " sec"
                << std::endl;
      std::cout << "****** Stabilisation: " << assembling_info.time_stab << " sec" << std::endl;
      std::cout << "****** Static condensation: " << assembling_info.time_statcond << " sec"
                << std::endl;
      std::cout << "**** Solver time: " << solve_info.time_solver << " sec" << std::endl;
      std::cout << "**** Postprocess time: " << post_info.time_postprocess << " sec" << std::endl;
      std::cout << "***********************************************************" << std::endl;
   }

   std::cout << "Discetisation h: " << disk::mesh_h(msh) << std::endl;
   std::cout << "L2 error: " << le.compute_l2_error(solution) << std::endl;
   //    le.plot_solution_at_gausspoint("sol_elas_2d.msh");
   //    le.plot_l2error_at_gausspoint("error_gp_2d_.msh", solution_lin);
   //    le.compute_deformed("deforme2d.msh");
   le.compute_discontinuous_displacement("depl2d.msh");
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_linear_elasticity_solver(const Mesh<T, 3, Storage>& msh,
                             run_params&                rp,
                             ElasticityParameters       material_data)
{
   typedef Mesh<T, 3, Storage> mesh_type;
   typedef static_vector<T, 3> result_type;

   timecounter tc;
   tc.tic();

   auto load = [material_data](const point<T, 3>& p) -> auto
   {
      const T mu     = material_data.mu;
      const T lambda = material_data.lambda;

      T fx = M_PI * M_PI *
             (12 * lambda * cos(2 * M_PI * p.x()) * sin(2 * M_PI * p.y()) * sin(2 * M_PI * p.z()) +
              8 * (3 * cos(2 * M_PI * p.x()) - 1) * sin(2 * M_PI * p.y()) * sin(2 * M_PI * p.z()) -
              cos(M_PI * p.x()) * sin(M_PI * (p.y() + p.z())) +
              (1 + 3. / (1 + lambda)) * sin(M_PI * p.x()) * sin(M_PI * p.y()) * sin(M_PI * p.z()));

      T fy = M_PI * M_PI *
             (12 * lambda * cos(2 * M_PI * p.y()) * sin(2 * M_PI * p.x()) * sin(2 * M_PI * p.z()) +
              8 * (3 * cos(2 * M_PI * p.y()) - 1) * sin(2 * M_PI * p.x()) * sin(2 * M_PI * p.z()) -
              cos(M_PI * p.y()) * sin(M_PI * (p.x() + p.z())) +
              (1 + 3. / (1 + lambda)) * sin(M_PI * p.x()) * sin(M_PI * p.y()) * sin(M_PI * p.z()));

      T fz = M_PI * M_PI *
             (12 * lambda * cos(2 * M_PI * p.z()) * sin(2 * M_PI * p.y()) * sin(2 * M_PI * p.x()) +
              8 * (3 * cos(2 * M_PI * p.z()) - 1) * sin(2 * M_PI * p.y()) * sin(2 * M_PI * p.x()) -
              cos(M_PI * p.z()) * sin(M_PI * (p.y() + p.x())) +
              (1 + 3. / (1 + lambda)) * sin(M_PI * p.x()) * sin(M_PI * p.y()) * sin(M_PI * p.z()));

      return result_type{fx, fy, fz};
   };

   auto solution = [material_data](const point<T, 3>& p) -> auto
   {
      const T mu     = material_data.mu;
      const T lambda = material_data.lambda;
      T       fx = sin(2 * M_PI * p.y()) * sin(2 * M_PI * p.z()) * (-1 + cos(2 * M_PI * p.x())) +
             (1. / (1 + lambda)) * sin(M_PI * p.x()) * sin(M_PI * p.y()) * sin(M_PI * p.z());
      T fy = sin(2 * M_PI * p.z()) * sin(2 * M_PI * p.x()) * (-1 + cos(2 * M_PI * p.y())) +
             (1. / (1 + lambda)) * sin(M_PI * p.x()) * sin(M_PI * p.y()) * sin(M_PI * p.z());
      T fz = sin(2 * M_PI * p.x()) * sin(2 * M_PI * p.y()) * (-1 + cos(2 * M_PI * p.z())) +
             (1. / (1 + lambda)) * sin(M_PI * p.x()) * sin(M_PI * p.y()) * sin(M_PI * p.z());
      return result_type{fx, fy, fz};
   };

   typedef BoundaryConditions<mesh_type, decltype(solution)> Bnd_type;
   Bnd_type                                                  bnd(msh);

   bnd.addDirichletEverywhere(solution);

   linear_elasticity_solver<mesh_type, Bnd_type> le(msh, bnd, material_data, rp.degree);
   le.verbose(rp.verbose);

   le.changeElasticityParameters(material_data);

   assembly_info assembling_info = le.assemble(load);

   solver_info solve_info = le.solve();

   postprocess_info post_info = le.postprocess(load);

   tc.toc();

   if (le.verbose()) {
      std::cout << std::endl;
      std::cout << "***********************************************************" << std::endl;
      std::cout << "** Time to solve the problem " << tc.to_double() << " sec" << std::endl;
      std::cout << "**** Assembly time: " << assembling_info.time_assembly << " sec" << std::endl;
      std::cout << "****** Gradient reconstruction: " << assembling_info.time_gradrec << " sec"
                << std::endl;
      std::cout << "****** Divergence reconstruction: " << assembling_info.time_divrec << " sec"
                << std::endl;
      std::cout << "****** Stabilisation: " << assembling_info.time_stab << " sec" << std::endl;
      std::cout << "****** Static condensation: " << assembling_info.time_statcond << " sec"
                << std::endl;
      std::cout << "**** Solver time: " << solve_info.time_solver << " sec" << std::endl;
      std::cout << "**** Postprocess time: " << post_info.time_postprocess << " sec" << std::endl;
      std::cout << "***********************************************************" << std::endl;
   }

   std::cout << "Discetisation h: " << disk::mesh_h(msh) << std::endl;
   std::cout << "L2 error: " << le.compute_l2_error(solution) << std::endl;
   //    le.plot_solution_at_gausspoint("sol_elas_3d.msh");
   //    le.plot_l2error_at_gausspoint("error_gp_3d_.msh", solution);
   //    le.compute_deformed("deforme3d.msh");
   //    le.compute_discontinuous_solution("depl3d.msh");
}

int
main(int argc, char** argv)
{
   using RealType = double;

   char* mesh_filename = nullptr;
   int   degree        = 1;
   int   l             = 0;

   run_params rp;
   rp.degree  = 1;
   rp.l       = 0;
   rp.verbose = true;

   // Elasticity Parameters
   ElasticityParameters material_data;

   material_data.mu     = 1.0;
   material_data.lambda = 10E5;

   int ch;

   while ((ch = getopt(argc, argv, "k:l:v:")) != -1) {
      switch (ch) {
         case 'k':
            degree = atoi(optarg);
            if (degree < 0) {
               std::cout << "Degree must be positive. Falling back to 1." << std::endl;
               degree = 1;
            }
            rp.degree = degree;
            break;

         case 'l':
            rp.l = atoi(optarg);
            if (rp.l < -1 or rp.l > 1) {
               std::cout << "l can be -1, 0 or 1. Falling back to 0." << std::endl;
               rp.l = 0;
            }
            break;

         case 'v': rp.verbose = true; break;

         case '?':
         default: std::cout << "wrong arguments" << std::endl; exit(1);
      }
   }

   argc -= optind;
   argv += optind;

   mesh_filename = argv[0];

   /* FVCA5 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$"))) {
      std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
      auto msh = disk::load_fvca5_2d_mesh<RealType>(mesh_filename);
      run_linear_elasticity_solver(msh, rp, material_data);
      return 0;
   }

   /* Netgen 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$"))) {
      std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
      auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);
      run_linear_elasticity_solver(msh, rp, material_data);
      return 0;
   }

   /* DiSk++ cartesian 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.quad$"))) {
      std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
      auto msh = disk::load_cartesian_2d_mesh<RealType>(mesh_filename);
      run_linear_elasticity_solver(msh, rp, material_data);
      return 0;
   }

   /* Netgen 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$"))) {
      std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
      auto msh = disk::load_netgen_3d_mesh<RealType>(mesh_filename);
      run_linear_elasticity_solver(msh, rp, material_data);
      return 0;
   }

   /* DiSk++ cartesian 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.hex$"))) {
      std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
      auto msh = disk::load_cartesian_3d_mesh<RealType>(mesh_filename);
      run_linear_elasticity_solver(msh, rp, material_data);
      return 0;
   }

   std::cout << "Unkwnon mesh format" << std::endl;
   return 0;
}
