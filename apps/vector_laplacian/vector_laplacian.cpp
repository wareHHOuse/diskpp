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

#include "colormanip.h"

#include "geometry/geometry.hpp"
#include "loaders/loader.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "vector_laplacian_solver.hpp"

void
usage(const char* progname)
{
   printf("Usage: %s <options> <filename>\n\n", progname);
   printf("    -n: number of elements for 1D mesh (default)\n");
   printf("    -k: face degree (>=0)\n");
   printf("    -l: difference beetween cell and face degree (-1 <= l <= 1) \n");
   printf("    -v: verbose\n");
}

struct run_params
{
   size_t degree;
   int    l;
   bool   verbose;
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_vector_laplacian_solver(const Mesh<T, 1, Storage>& msh,
                            const run_params&          rp,
                            LaplacianParameters        material_data)
{
   typedef Mesh<T, 1, Storage> mesh_type;

   timecounter tc;
   tc.tic();

   auto load = [material_data](const point<T, 1>& p) -> T {
      return material_data.lambda * M_PI * M_PI * sin(M_PI * p.x());
   };

   auto solution = [material_data](const point<T, 1>& p) -> T { return sin(M_PI * p.x()); };

   auto gradient = [material_data](const point<T, 1>& p) -> T { return M_PI * cos(M_PI * p.x()); };

   vector_laplacian_solver<mesh_type> vl(msh, rp.degree, rp.l);
   vl.verbose(rp.verbose);

   vl.changeLaplacianParameters(material_data);

   assembly_info assembling_info = vl.assemble(load, solution);

   if (vl.verbose()) {
      std::cout << "Assembling: " << assembling_info.time_assembly << " sec" << '\n';
   }

   solver_info solve_info = vl.solve();

   if (vl.verbose()) {
      std::cout << "Total time to solve the problem: " << solve_info.time_solver << " sec" << '\n';
   }

   postprocess_info post_info = vl.postprocess(load);

   if (vl.verbose()) {
      std::cout << "Post-Processing: " << post_info.time_postprocess << " sec" << '\n';
   }

   std::cout << "Discetisation h: " << disk::mesh_h(msh) << std::endl;
   std::cout << "L2 error: " << vl.compute_l2_error(solution) << std::endl;
   std::cout << "L2 gradient error: " << vl.compute_l2_gradient_error(gradient) << std::endl;

   tc.toc();

   if (vl.verbose()) {
      std::cout << std::endl;
      std::cout << "************************************************************" << std::endl;
      std::cout << "** Time to solve the problem " << tc.to_double() << " sec" << std::endl;
      std::cout << "**** Assembly time: " << assembling_info.time_assembly << " sec" << std::endl;
      std::cout << "****** Gradient reconstruction: " << assembling_info.time_gradrec << " sec"
                << std::endl;
      std::cout << "****** Stabilisation: " << assembling_info.time_stab << " sec" << std::endl;
      std::cout << "****** Static condensation: " << assembling_info.time_statcond << " sec"
                << std::endl;
      std::cout << "**** Solver time: " << solve_info.time_solver << " sec" << std::endl;
      std::cout << "**** Postprocess time: " << post_info.time_postprocess << " sec" << std::endl;
      std::cout << "***********************************************************" << std::endl;
   }
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_vector_laplacian_solver(const Mesh<T, 2, Storage>& msh,
                            run_params&                rp,
                            LaplacianParameters        material_data)
{
   typedef Mesh<T, 2, Storage> mesh_type;
   typedef static_vector<T, 2> result_type;

   timecounter tc;
   tc.tic();

   auto load = [material_data](const point<T, 2>& p) -> result_type {
      T fx = 2. * material_data.lambda * M_PI * M_PI * sin(M_PI * p.x()) * sin(M_PI * p.y());
      T fy = 2. * material_data.lambda * M_PI * M_PI * cos(M_PI * p.x()) * cos(M_PI * p.y());

      return result_type{fx, fy};
   };

   auto solution = [material_data](const point<T, 2>& p) -> result_type {
      T fx = material_data.lambda * sin(M_PI * p.x()) * sin(M_PI * p.y());
      T fy = material_data.lambda * cos(M_PI * p.x()) * cos(M_PI * p.y());

      return result_type{fx, fy};
   };

   vector_laplacian_solver<mesh_type> vl(msh, rp.degree);
   vl.verbose(rp.verbose);

   vl.changeLaplacianParameters(material_data);

   assembly_info assembling_info = vl.assemble(load, solution);

   if (vl.verbose()) {
      std::cout << "Assembling: " << assembling_info.time_assembly << " sec" << '\n';
   }

   solver_info solve_info = vl.solve();

   if (vl.verbose()) {
      std::cout << "Total time to solve the problem: " << solve_info.time_solver << " sec" << '\n';
   }

   postprocess_info post_info = vl.postprocess(load);

   if (vl.verbose()) {
      std::cout << "Post-Processing: " << post_info.time_postprocess << " sec" << '\n';
   }

   std::cout << "Discetisation h: " << disk::mesh_h(msh) << std::endl;
   std::cout << "L2 error: " << vl.compute_l2_error(solution) << std::endl;

   tc.toc();

   if (vl.verbose()) {
      std::cout << std::endl;
      std::cout << "************************************************************" << std::endl;
      std::cout << "** Time to solve the problem " << tc.to_double() << " sec" << std::endl;
      std::cout << "**** Assembly time: " << assembling_info.time_assembly << " sec" << std::endl;
      std::cout << "****** Gradient reconstruction: " << assembling_info.time_gradrec << " sec"
                << std::endl;
      std::cout << "****** Stabilisation: " << assembling_info.time_stab << " sec" << std::endl;
      std::cout << "****** Static condensation: " << assembling_info.time_statcond << " sec"
                << std::endl;
      std::cout << "**** Solver time: " << solve_info.time_solver << " sec" << std::endl;
      std::cout << "**** Postprocess time: " << post_info.time_postprocess << " sec" << std::endl;
      std::cout << "***********************************************************" << std::endl;
   }
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_vector_laplacian_solver(const Mesh<T, 3, Storage>& msh,
                            run_params&                rp,
                            LaplacianParameters        material_data)
{
   typedef Mesh<T, 3, Storage> mesh_type;
   typedef static_vector<T, 3> result_type;

   timecounter tc;
   tc.tic();

   auto load = [material_data](const point<T, 3>& p) -> auto
   {
      T fx = 2. * material_data.lambda * M_PI * M_PI * cos(M_PI * p.x()) * sin(M_PI * p.y());

      T fy = 2. * material_data.lambda * M_PI * M_PI * cos(M_PI * p.y()) * sin(M_PI * p.z());

      T fz = 2. * material_data.lambda * M_PI * M_PI * cos(M_PI * p.z()) * sin(M_PI * p.x());

      return result_type{fx, fy, fz};
   };

   auto solution = [material_data](const point<T, 3>& p) -> auto
   {
      T fx = cos(M_PI * p.x()) * sin(M_PI * p.y());
      T fy = cos(M_PI * p.y()) * sin(M_PI * p.z());
      T fz = cos(M_PI * p.z()) * sin(M_PI * p.x());
      return result_type{fx, fy, fz};
   };

   vector_laplacian_solver<mesh_type> vl(msh, rp.degree);
   vl.verbose(rp.verbose);

   vl.changeLaplacianParameters(material_data);

   assembly_info assembling_info = vl.assemble(load, solution);

   if (vl.verbose()) {
      std::cout << "Assembling: " << assembling_info.time_assembly << " sec" << '\n';
   }

   solver_info solve_info = vl.solve();

   if (vl.verbose()) {
      std::cout << "Total time to solve the problem: " << solve_info.time_solver << " sec" << '\n';
   }

   postprocess_info post_info = vl.postprocess(load);

   if (vl.verbose()) {
      std::cout << "Post-Processing: " << post_info.time_postprocess << " sec" << '\n';
   }

   std::cout << "Discetisation h: " << disk::mesh_h(msh) << std::endl;
   std::cout << "L2 error: " << vl.compute_l2_error(solution) << std::endl;

   tc.toc();

   if (vl.verbose()) {
      std::cout << std::endl;
      std::cout << "***********************************************************" << std::endl;
      std::cout << "** Time to solve the problem " << tc.to_double() << " sec" << std::endl;
      std::cout << "**** Assembly time: " << assembling_info.time_assembly << " sec" << std::endl;
      std::cout << "****** Gradient reconstruction: " << assembling_info.time_gradrec << " sec"
                << std::endl;
      std::cout << "****** Stabilisation: " << assembling_info.time_stab << " sec" << std::endl;
      std::cout << "****** Static condensation: " << assembling_info.time_statcond << " sec"
                << std::endl;
      std::cout << "**** Solver time: " << solve_info.time_solver << " sec" << std::endl;
      std::cout << "**** Postprocess time: " << post_info.time_postprocess << " sec" << std::endl;
      std::cout << "***********************************************************" << std::endl;
   }
}

int
main(int argc, char** argv)
{
   using RealType = double;

   char* mesh_filename = nullptr;
   int   degree        = 1;
   int   l             = 0;
   int   elems_1d      = 8;

   run_params rp;
   rp.degree  = 1;
   rp.l       = 0;
   rp.verbose = true;

   // Elasticity Parameters
   LaplacianParameters material_data;

   material_data.lambda = 1.0;

   int ch;

   while ((ch = getopt(argc, argv, "k:l:n:v")) != -1) {
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

         case 'n':
            elems_1d = atoi(optarg);
            if (elems_1d < 0) {
               std::cout << "Num of elems must be positive. Falling back to 8." << std::endl;
               elems_1d = 8;
            }
            break;

         case 'v': rp.verbose = true; break;

         case 'h':
         case '?':
         default:
            std::cout << "wrong arguments" << std::endl;
            usage(argv[0]);
            exit(1);
      }
   }

   argc -= optind;
   argv += optind;

   if (argc == 0) {
      std::cout << "Mesh format: 1D uniform" << std::endl;
      auto msh = disk::load_uniform_1d_mesh<RealType>(0, 1, elems_1d);
      run_vector_laplacian_solver(msh, rp, material_data);
      return 0;
   }

   mesh_filename = argv[0];

   /* FVCA5 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$"))) {
      std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
      auto msh = disk::load_fvca5_2d_mesh<RealType>(mesh_filename);
      run_vector_laplacian_solver(msh, rp, material_data);
      return 0;
   }

   /* Netgen 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$"))) {
      std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
      auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);
      run_vector_laplacian_solver(msh, rp, material_data);
      return 0;
   }

   /* DiSk++ cartesian 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.quad$"))) {
      std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
      auto msh = disk::load_cartesian_2d_mesh<RealType>(mesh_filename);
      run_vector_laplacian_solver(msh, rp, material_data);
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
      run_vector_laplacian_solver(msh, rp, material_data);
      return 0;
   }

   /* DiSk++ cartesian 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.hex$"))) {
      std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
      auto msh = disk::load_cartesian_3d_mesh<RealType>(mesh_filename);
      run_vector_laplacian_solver(msh, rp, material_data);
      return 0;
   }

   std::cout << "Unkwnon mesh format" << std::endl;
   return 0;
}
