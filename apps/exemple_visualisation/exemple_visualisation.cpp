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

#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <unistd.h>

#include <map>

#include "colormanip.h"

#include "config.h"

#ifdef HAVE_SOLVER_WRAPPERS
#include "agmg/agmg.hpp"
#endif

#include "hho/hho.hpp"
#include "loaders/loader.hpp"
#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "exemple_visualisation_solver.hpp"

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

   diffusion_solver<mesh_type> dp(msh, rp.degree);
   dp.verbose(rp.verbose);

   dp.assemble(load, solution);
   dp.solve();
   dp.postprocess(load);
   dp.plot_solution_at_gausspoint("solgp1d.msh");
   dp.plot_l2error_at_gausspoint("errorgp1d.msh", solution);
   dp.plot_conforme_solution("visu1d_conforme.msh");
   dp.plot_discontinuous_solution("visu1d.msh");
   dp.plot_deformed_conforme("visu1d_deformed_conforme.msh");
   dp.plot_deformed_discontinuous("visu1d_deformed_discontinuous.msh");
   dp.saveMesh("saveMesh1D.msh");
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

   diffusion_solver<mesh_type> dp(msh, rp.degree);
   dp.verbose(rp.verbose);

   dp.assemble(load, solution);
   dp.solve();
   dp.postprocess(load);
   dp.plot_solution_at_gausspoint("solgp2d.msh");
   dp.plot_l2error_at_gausspoint("errorgp2d.msh", solution);
   dp.plot_conforme_solution("visu2d_conforme.msh");
   dp.plot_discontinuous_solution("visu2d.msh");
   dp.plot_deformed_conforme("visu2d_deformed_conforme.msh");
   dp.plot_deformed_discontinuous("visu2d_deformed_discontinuous.msh");
   dp.saveMesh("saveMesh2D.msh");
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

   diffusion_solver<mesh_type> dp(msh, rp.degree);
   dp.verbose(rp.verbose);

   dp.assemble(load, solution);
   dp.solve();
   dp.postprocess(load);
   dp.plot_solution_at_gausspoint("solgp3d.msh");
   dp.plot_l2error_at_gausspoint("errorgp3d.msh", solution);
   dp.plot_conforme_solution("visu3d_conforme.msh");
   dp.plot_discontinuous_solution("visu3d.msh");
   dp.saveMesh("saveMesh3D.msh");
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

   while ((ch = getopt(argc, argv, "k:l:n:p:v")) != -1) {
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
            if (l < -1 or l > 1) {
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

         case 'p': plot_filename = optarg; break;

         case 'v': rp.verbose = true; break;

         case 'h':
         case '?':
         default: std::cout << "wrong arguments" << std::endl; exit(1);
      }
   }

   argc -= optind;
   argv += optind;

   if (argc == 0) {
      std::cout << "Mesh format: 1D uniform" << std::endl;
      auto msh = disk::load_uniform_1d_mesh<RealType>(0, 1, elems_1d);
      run_diffusion_solver(msh, rp);
      return 0;
   }

   mesh_filename = argv[0];

   /* FVCA5 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$"))) {
      std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
      auto msh = disk::load_fvca5_2d_mesh<RealType>(mesh_filename);
      run_diffusion_solver(msh, rp);
      return 0;
   }

   /* Netgen 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$"))) {
      std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
      auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);
      run_diffusion_solver(msh, rp);
      return 0;
   }

   /* DiSk++ cartesian 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.quad$"))) {
      std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
      auto msh = disk::load_cartesian_2d_mesh<RealType>(mesh_filename);
      run_diffusion_solver(msh, rp);
      return 0;
   }

   /* Netgen 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$"))) {
      std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
      auto msh = disk::load_netgen_3d_mesh<RealType>(mesh_filename);
      run_diffusion_solver(msh, rp);
      return 0;
   }

   /* DiSk++ cartesian 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.hex$"))) {
      std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
      auto msh = disk::load_cartesian_3d_mesh<RealType>(mesh_filename);
      run_diffusion_solver(msh, rp);
      return 0;
   }

   /* Medit 2d*/
   if (std::regex_match(mesh_filename, std::regex(".*\\.medit2d$"))) {
      std::cout << "Guessed mesh format: Medit format" << std::endl;
      auto msh = disk::load_medit_2d_mesh<RealType>(mesh_filename);
      run_diffusion_solver(msh, rp);
      visu::Gmesh gmsh = visu::convertMesh(msh);
      gmsh.writeGmesh("test_M2.msh", 2);
      return 0;
   }

   /* Medit 3d*/
   //  if (std::regex_match(mesh_filename, std::regex(".*\\.medit3d$") ))
   //  {
   //     std::cout << "Guessed mesh format: Medit format" << std::endl;
   //     auto msh = disk::load_medit_3d_mesh<RealType>(mesh_filename);
   //     ///run_diffusion_solver(msh, rp);
   //     visu::Gmesh gmsh = visu::convertMesh(msh);
   //     gmsh.writeGmesh("test_M3.msh",2);
   //     return 0;
   //  }
}
