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

#include "loaders/loader.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "diffusion_RT_solver.hpp"

struct error_type
{
   size_t degree;
   size_t nb_dof;
   double h;
   double error_depl;
   double error_grad;
};

struct run_params
{
   size_t degree;
   int    l;
   bool   verbose;
};

void
usage(const char* progname)
{
   printf("Usage: %s <options> <filename>\n\n", progname);
   printf("    -1: test 1D mesh (default)\n");
   printf("    -2: test 2D mesh\n");
   printf("    -3: test 3D mesh\n");
   printf("    -k: degree (>0)\n");
   printf("    -v: verbose \n");
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
error_type
run_diffusion_solver(const Mesh<T, 1, Storage>& msh, const run_params& rp)
{
   typedef Mesh<T, 1, Storage> mesh_type;

   timecounter tc;
   tc.tic();

   auto load = [](const point<T, 1>& p) -> auto { return M_PI * M_PI * sin(p.x() * M_PI); };

   auto solution = [](const point<T, 1>& p) -> auto { return sin(p.x() * M_PI); };

   auto gradient = [](const point<T, 1>& p) -> auto { return M_PI * cos(M_PI * p.x()); };

   diffusion_solver<mesh_type> dp(msh, rp.degree, rp.l);
   dp.verbose(rp.verbose);

   auto assembly_info = dp.assemble(load, solution);
   auto solve_info    = dp.solve();
   auto postpro_info  = dp.postprocess(load);

   error_type error;
   error.h          = average_diameter(msh);
   error.degree     = rp.degree;
   error.nb_dof     = dp.getDofs();
   error.error_depl = dp.compute_l2_error(solution);
   error.error_grad = dp.compute_l2_gradient_error(gradient);

   return error;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
error_type
run_diffusion_solver(const Mesh<T, 2, Storage>& msh, const run_params& rp)
{
   typedef Mesh<T, 2, Storage> mesh_type;
   typedef static_vector<T, 2> result_type;

   auto load = [](const point<T, 2>& p) -> auto
   {
      return 2.0 * M_PI * M_PI * sin(p.x() * M_PI) * sin(p.y() * M_PI);
   };

   auto solution = [](const point<T, 2>& p) -> auto
   {
      return sin(p.x() * M_PI) * sin(p.y() * M_PI);
   };

   auto gradient = [](const point<T, 2>& p) -> result_type {
      result_type grad = result_type::Zero();

      grad(0, 0) = M_PI * cos(p.x() * M_PI) * sin(p.y() * M_PI);
      grad(1, 0) = M_PI * sin(p.x() * M_PI) * cos(p.y() * M_PI);

      return grad;
   };

   diffusion_solver<mesh_type> dp(msh, rp.degree, rp.l);
   dp.verbose(rp.verbose);

   auto assembly_info = dp.assemble(load, solution);
   auto solve_info    = dp.solve();
   auto postpro_info  = dp.postprocess(load);

   error_type error;
   error.h          = average_diameter(msh);
   error.degree     = rp.degree;
   error.nb_dof     = dp.getDofs();
   error.error_depl = dp.compute_l2_error(solution);
   error.error_grad = dp.compute_l2_gradient_error(gradient);

   return error;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
error_type
run_diffusion_solver(const Mesh<T, 3, Storage>& msh, const run_params& rp)
{
   typedef Mesh<T, 3, Storage> mesh_type;
   typedef static_vector<T, 3> result_type;

   auto load = [](const point<T, 3>& p) -> auto
   {
      return 3.0 * M_PI * M_PI * sin(p.x() * M_PI) * sin(p.y() * M_PI) * sin(p.z() * M_PI);
   };

   auto solution = [](const point<T, 3>& p) -> auto
   {
      return sin(p.x() * M_PI) * sin(p.y() * M_PI) * sin(p.z() * M_PI);
   };

   auto gradient = [](const point<T, 3>& p) -> result_type {
      result_type grad = result_type::Zero();

      grad(0, 0) = M_PI * cos(p.x() * M_PI) * sin(p.y() * M_PI) * sin(p.z() * M_PI);
      grad(1, 0) = M_PI * sin(p.x() * M_PI) * cos(p.y() * M_PI) * sin(p.z() * M_PI);
      grad(2, 0) = M_PI * sin(p.x() * M_PI) * sin(p.y() * M_PI) * cos(p.z() * M_PI);

      return grad;
   };

   diffusion_solver<mesh_type> dp(msh, rp.degree, rp.l);
   dp.verbose(rp.verbose);

   auto assembly_info = dp.assemble(load, solution);
   auto solve_info    = dp.solve();
   auto postpro_info  = dp.postprocess(load);

   error_type error;
   error.h          = average_diameter(msh);
   error.degree     = rp.degree;
   error.nb_dof     = dp.getDofs();
   error.error_depl = dp.compute_l2_error(solution);
   error.error_grad = dp.compute_l2_gradient_error(gradient);

   return error;
}

void
printResults(const std::vector<error_type>& error)
{
   if (error.size() > 0) {
      std::ios::fmtflags f(std::cout.flags());
      std::cout.precision(4);
      std::cout.setf(std::iostream::scientific, std::iostream::floatfield);

      std::cout << "Convergence test for k = " << error[0].degree << std::endl;
      std::cout
        << "-----------------------------------------------------------------------------------"
        << std::endl;
      std::cout
        << "| Size mesh  | Displacement | Convergence |  Gradient  | Convergence |    Total   |"
        << std::endl;
      std::cout
        << "|    h       |   L2 error   |     rate    |  L2 error  |     rate    | faces DOF  |"
        << std::endl;
      std::cout
        << "-----------------------------------------------------------------------------------"
        << std::endl;

      std::string s_dof = " " + std::to_string(error[0].nb_dof) + "                  ";
      s_dof.resize(10);

      std::cout << "| " << error[0].h << " |  " << error[0].error_depl << "  | "
                << "     -     "
                << " | " << error[0].error_grad << " | "
                << "     -     "
                << " | " << s_dof << " |" << std::endl;

      for (size_t i = 1; i < error.size(); i++) {
         s_dof = " " + std::to_string(error[i].nb_dof) + "                  ";
         s_dof.resize(10);
         double rate_depl = (log10(error[i - 1].error_depl) - log10(error[i].error_depl)) /
                            (log10(error[i - 1].h) - log10(error[i].h));
         double rate_grad = (log10(error[i - 1].error_grad) - log10(error[i].error_grad)) /
                            (log10(error[i - 1].h) - log10(error[i].h));

         std::cout << "| " << error[i].h << " |  " << error[i].error_depl << "  |  " << rate_depl
                   << " | " << error[i].error_grad << " |  " << rate_grad << " | " << s_dof << " |"
                   << std::endl;
      }

      std::cout
        << "-----------------------------------------------------------------------------------"
        << std::endl;
      std::cout << "  " << std::endl;
      std::cout.flags(f);
   } else
      std::cout << "The file error is empty" << std::endl;
}

template<typename T>
void
test_1D(const run_params& rp)
{
   size_t runs = 4;

   std::vector<size_t> n_elems;
   n_elems.push_back(12);
   n_elems.push_back(24);
   n_elems.push_back(48);
   n_elems.push_back(96);

   std::vector<error_type> error_sumup;

   for (size_t i = 0; i < runs; i++) {
      const auto msh = disk::load_uniform_1d_mesh<T>(0, 1, n_elems[i]);
      error_sumup.push_back(run_diffusion_solver(msh, rp));
   }

   printResults(error_sumup);
}

template<typename T>
void
test_triangles_fvca5(const run_params& rp)
{
   size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_1.typ1");
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_2.typ1");
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_3.typ1");
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_4.typ1");
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_5.typ1");

   std::vector<error_type> error_sumup;

   for (size_t i = 0; i < runs; i++) {
      auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_diffusion_solver(msh, rp));
   }
   printResults(error_sumup);
}

template<typename T>
void
test_triangles_netgen(const run_params& rp)
{
   size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri01.mesh2d");
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri02.mesh2d");
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri03.mesh2d");
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri04.mesh2d");
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri05.mesh2d");

   std::vector<error_type> error_sumup;

   for (size_t i = 0; i < runs; i++) {
      auto msh = disk::load_netgen_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_diffusion_solver(msh, rp));
   }
   printResults(error_sumup);
}

template<typename T>
void
test_hexagons(const run_params& rp)
{
   size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_1.typ1");
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_2.typ1");
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_3.typ1");
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_4.typ1");
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_5.typ1");

   std::vector<error_type> error_sumup;

   for (size_t i = 0; i < runs; i++) {
      auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_diffusion_solver(msh, rp));
   }
   printResults(error_sumup);
}

template<typename T>
void
test_kershaws(const run_params& rp)
{
   size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_1.typ1");
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_2.typ1");
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_3.typ1");
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_4.typ1");
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_5.typ1");

   std::vector<error_type> error_sumup;

   for (size_t i = 0; i < runs; i++) {
      auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_diffusion_solver(msh, rp));
   }
   printResults(error_sumup);
}

template<typename T>
void
test_quads_fvca5(const run_params& rp)
{
   size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_1.typ1");
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_2.typ1");
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_3.typ1");
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_4.typ1");
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_5.typ1");

   std::vector<error_type> error_sumup;

   for (size_t i = 0; i < runs; i++) {
      auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_diffusion_solver(msh, rp));
   }
   printResults(error_sumup);
}

template<typename T>
void
test_quads_diskpp(const run_params& rp)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-4-4.quad");
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-8-8.quad");
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-16-16.quad");
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-32-32.quad");
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-256-256.quad");

   std::vector<error_type> error_sumup;

   for (size_t i = 0; i < runs; i++) {
      auto msh = disk::load_cartesian_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_diffusion_solver(msh, rp));
   }
   printResults(error_sumup);
}

template<typename T>
void
test_hexahedra_diskpp(const run_params& rp)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-2-2-2.hex");
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-4-4-4.hex");
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-8-8-8.hex");
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-16-16-16.hex");
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-32-32-32.hex");

   std::vector<error_type> error_sumup;

   for (int i = 0; i < runs; i++) {
      auto msh = disk::load_cartesian_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_diffusion_solver(msh, rp));
   }
   printResults(error_sumup);
}

template<typename T>
void
test_hexahedra_fvca6(const run_params& rp)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_2x2x2.msh");
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_4x4x4.msh");
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_8x8x8.msh");
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_16x16x16.msh");
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_32x32x32.hex");

   std::vector<error_type> error_sumup;

   for (int i = 0; i < runs; i++) {
      auto msh = disk::load_fvca6_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_diffusion_solver(msh, rp));
   }
   printResults(error_sumup);
}

template<typename T>
void
test_tetrahedra_netgen(const run_params& rp)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet0.mesh");
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet1.mesh");
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet2.mesh");
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet3.mesh");
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet4.mesh");

   std::vector<error_type> error_sumup;

   for (int i = 0; i < runs; i++) {
      auto msh = disk::load_netgen_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_diffusion_solver(msh, rp));
   }
   printResults(error_sumup);
}

template<typename T>
void
test_polyhedra_fvca6(const run_params& rp)
{
   size_t runs = 3;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_10.msh");
   paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_20.msh");
   paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_30.msh");
   paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_40.msh");

   std::vector<error_type> error_sumup;

   for (int i = 0; i < runs; i++) {
      auto msh = disk::load_fvca6_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_diffusion_solver(msh, rp));
   }
   printResults(error_sumup);
}

template<typename T>
void
test_tetrahedra_fvca6(const run_params& rp)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.0.msh");
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.1.msh");
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.2.msh");
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.3.msh");
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.4.msh");

   std::vector<error_type> error_sumup;

   for (int i = 0; i < runs; i++) {
      auto msh = disk::load_fvca6_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_diffusion_solver(msh, rp));
   }
   printResults(error_sumup);
}

int
main(int argc, char** argv)
{
   using RealType = double;

   char*  mesh_filename = nullptr;
   char*  plot_filename = nullptr;
   int    degree        = 1;
   int    l             = 0;
   size_t dim           = 1;

   run_params rp;
   rp.degree  = 1;
   rp.l       = 0;
   rp.verbose = false;

   int ch;

   while ((ch = getopt(argc, argv, "123k:v")) != -1) {
      switch (ch) {
         case '1': dim = 1; break;
         case '2': dim = 2; break;
         case '3': dim = 3; break;

         case 'k':
            degree = atoi(optarg);
            if (degree < 0) {
               std::cout << "Degree must be positive. Falling back to 0." << std::endl;
               degree = 0;
            }
            rp.degree = degree;
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

   timecounter tc;

   std::cout << " Test convergence rates for: " << std::endl;
   std::cout << " ** Face_Degree = " << rp.degree << std::endl;
   std::cout << " ** Cell_Degree  = " << rp.degree + rp.l << std::endl;
   std::cout << " " << std::endl;

   if (dim == 3) {
      tc.tic();
      std::cout << "-Tetrahedras fvca6:" << std::endl;
      test_tetrahedra_fvca6<RealType>(rp);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " " << std::endl;

      tc.tic();
      std::cout << "-Tetrahedras netgen:" << std::endl;
      test_tetrahedra_netgen<RealType>(rp);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " " << std::endl;

      // tc.tic();
      // std::cout << "-Hexahedras fvca6:" << std::endl;
      // test_hexahedra_fvca6<RealType>(rp);
      // tc.toc();
      // std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      // std::cout << " " << std::endl;

      // tc.tic();
      // std::cout << "-Hexahedras diskpp:" << std::endl;
      // test_hexahedra_diskpp<RealType>(rp);
      // tc.toc();
      // std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      // std::cout << " " << std::endl;

      // tc.tic();
      // std::cout << "-Polyhedra:" << std::endl;
      // test_polyhedra_fvca6<RealType>(rp);
      // tc.toc();
      // std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      // std::cout << " " << std::endl;
   } else if (dim == 2) {

      tc.tic();
      std::cout << "-Triangles fvca5:" << std::endl;
      test_triangles_fvca5<RealType>(rp);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " " << std::endl;

      tc.tic();
      std::cout << "-Triangles netgen:" << std::endl;
      test_triangles_netgen<RealType>(rp);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " " << std::endl;

      // tc.tic();
      // std::cout << "-Quadrangles fvca5:" << std::endl;
      // test_quads_fvca5<RealType>(rp);
      // tc.toc();
      // std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      // std::cout << " " << std::endl;

      // tc.tic();
      // std::cout << "-Quadrangles diskpp:" << std::endl;
      // test_quads_diskpp<RealType>(rp);
      // tc.toc();
      // std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      // std::cout << " " << std::endl;

      // tc.tic();
      // std::cout << "-Hexagons:" << std::endl;
      // test_hexagons<RealType>(rp);
      // tc.toc();
      // std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      // std::cout << " " << std::endl;

      // tc.tic();
      // std::cout << "-Kershaws:" << std::endl;
      // test_kershaws<RealType>(rp);
      // tc.toc();
      // std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      // std::cout << " " << std::endl;
   } else {
      tc.tic();
      std::cout << "-Test 1D:" << std::endl;
      test_1D<RealType>(rp);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
   }
}
