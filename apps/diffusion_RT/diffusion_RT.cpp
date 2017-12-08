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

#include "loaders/loader.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "diffusion_RT_solver.hpp"

void
usage(const char* progname)
{
   printf("Usage: %s <options> <filename>\n\n", progname);
   printf("    -n: number of elemnt in 1D (>0)\n");
   printf("    -k: degree (>0)\n");
   printf("    -v: verbose \n");
}

struct run_params
{
   size_t degree;
   int    l;
   bool   verbose;
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_diffusion_solver(const Mesh<T, 1, Storage>& msh, run_params& rp)
{
   typedef Mesh<T, 1, Storage> mesh_type;

   auto load = [](const point<T, 1>& p) -> auto { return M_PI * M_PI * sin(p.x() * M_PI); };

   auto solution = [](const point<T, 1>& p) -> auto { return sin(p.x() * M_PI); };

   diffusion_solver<mesh_type> dp(msh, rp.degree, rp.l);
   dp.verbose(rp.verbose);

   dp.assemble(load, solution);
   dp.solve();
   dp.postprocess(load);
   std::cout << dp.compute_l2_error(solution) << std::endl;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_diffusion_solver(const Mesh<T, 2, Storage>& msh, run_params& rp)
{
   typedef Mesh<T, 2, Storage> mesh_type;

   auto load = [](const point<T, 2>& p) -> auto
   {

      return 2.0 * M_PI * M_PI * sin(p.x() * M_PI) * sin(p.y() * M_PI);
   };

   auto solution = [](const point<T, 2>& p) -> auto
   {
      return sin(p.x() * M_PI) * sin(p.y() * M_PI);
   };

   diffusion_solver<mesh_type> dp(msh, rp.degree, rp.l);
   dp.verbose(rp.verbose);

   dp.assemble(load, solution);
   dp.solve();
   dp.postprocess(load);
   std::cout << dp.compute_l2_error(solution) << std::endl;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_diffusion_solver(const Mesh<T, 3, Storage>& msh, run_params& rp)
{
   typedef Mesh<T, 3, Storage> mesh_type;

   auto load = [](const point<T, 3>& p) -> auto
   {
      return 3.0 * M_PI * M_PI * sin(p.x() * M_PI) * sin(p.y() * M_PI) * sin(p.z() * M_PI);
   };

   auto solution = [](const point<T, 3>& p) -> auto
   {
      return sin(p.x() * M_PI) * sin(p.y() * M_PI) * sin(p.z() * M_PI);
   };

   diffusion_solver<mesh_type> dp(msh, rp.degree, rp.l);
   dp.verbose(rp.verbose);

   dp.assemble(load, solution);
   dp.solve();
   dp.postprocess(load);
   std::cout << dp.compute_l2_error(solution) << std::endl;
}

int
main(int argc, char** argv)
{
   using RealType = double;

   char* mesh_filename = nullptr;
   char* plot_filename = nullptr;
   int   degree        = 1;
   int   l             = 0;
   int   elems_1d      = 8;

   run_params rp;
   rp.degree  = 1;
   rp.l       = 0;
   rp.verbose = false;

   int ch;

   while ((ch = getopt(argc, argv, "k:n:v")) != -1) {
      switch (ch) {
         case 'k':
            degree = atoi(optarg);
            if (degree < 0) {
               std::cout << "Degree must be positive. Falling back to 1." << std::endl;
               degree = 1;
            }
            rp.degree = degree;
            break;

         case 'n':
            elems_1d = atoi(optarg);
            if (elems_1d < 0) {
               std::cout << "Num of elems must be positive. Falling back to 8." << std::endl;
               elems_1d = 8;
            }
            break;

         case 'v': rp.verbose = true; break;

         case '?':
         default:
            std::cout << "wrong arguments" << std::endl;
            usage(argv[0]);
            exit(1);
      }
   }

   argc -= optind;
   argv += optind;

   mesh_filename = argv[0];

   if (argc == 0) {
      std::cout << "Mesh format: 1D uniform" << std::endl;
      const auto msh = disk::load_uniform_1d_mesh<RealType>(0, 1, elems_1d);
      run_diffusion_solver(msh, rp);
      return 0;
   }

   /* FVCA5 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$"))) {
      std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
      const auto msh = disk::load_fvca5_2d_mesh<RealType>(mesh_filename);
      run_diffusion_solver(msh, rp);
      return 0;
   }
   /* Netgen 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$"))) {
      std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
      const auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);

      run_diffusion_solver(msh, rp);
      return 0;
   }

   /* FVCA6 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.msh$"))) {
      std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
      const auto msh = disk::load_fvca6_3d_mesh<RealType>(mesh_filename);
      run_diffusion_solver(msh, rp);
      return 0;
   }

   /* Netgen 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$"))) {
      std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
      const auto msh = disk::load_netgen_3d_mesh<RealType>(mesh_filename);
      run_diffusion_solver(msh, rp);
      return 0;
   }
}
