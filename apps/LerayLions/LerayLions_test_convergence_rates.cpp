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

#include <iostream>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include <map>

#include "colormanip.h"

#include "../../config.h"

#ifdef HAVE_SOLVER_WRAPPERS
#include "agmg/agmg.hpp"
#endif

#include "loaders/loader.hpp"
#include "Parameters.hpp"
#include "Informations.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "LerayLions_solver.hpp"

struct error_type
{
   size_t  degree;
   size_t nb_dof;
   double h;
   double error_depl;
   double error_grad;
};


void
usage(const char *progname)
{
   printf("Usage: %s <options> <filename>\n\n", progname);
   printf("    -2: test 2D mesh (default)\n");
   printf("    -3: test 3D mesh\n");
   printf("    -r: reads parameters from an file\n");
   printf("    -p: Leray-Lions parameter (>=2)\n");
   printf("    -v: verbose\n");
}


template<template<typename, size_t , typename> class Mesh,
typename T, typename Storage>
error_type
run_leraylions_solver(const Mesh<T, 2, Storage>& msh, const ParamRun<T>& rp, const T leray_param)
{
   typedef T result_type;
   typedef static_vector<T, 2> result_grad_type;

   auto load = [leray_param](const point<T,2>& pt) -> result_type {
      const T p = leray_param;
      result_type norm_G = M_PI * sqrt( std::pow(cos(pt.x() * M_PI) * sin(pt.y() * M_PI),2) +
      std::pow(sin(pt.x() * M_PI) * cos(pt.y() * M_PI),2));

      T fx = sin(M_PI*pt.x())*sin(M_PI*pt.y())*
      (
         2.*std::pow(M_PI,2)*std::pow(norm_G, p-2.0) - (p-2)*std::pow(M_PI,4)*std::pow(norm_G, p-4.0)*(
            std::pow(cos(M_PI*pt.x()),2)*cos(2.*M_PI*pt.y()) + cos(2.*M_PI*pt.x())*std::pow(cos(M_PI*pt.y()),2))
      );

      return result_type{fx};
   };

   auto solution = [](const point<T,2>& pt) -> result_type {
      return sin(pt.x() * M_PI) * sin(pt.y() * M_PI);
   };

   auto gradient = [](const point<T,2>& pt) -> auto {
      result_grad_type grad = result_grad_type::Zero();

      grad(0) = M_PI * cos(pt.x() * M_PI) * sin(pt.y() * M_PI);
      grad(1) = M_PI * sin(pt.x() * M_PI) * cos(pt.y() * M_PI);

      return grad;
   };


   leraylions_solver<Mesh, T, 2, Storage,  point<T, 2> > nl(msh, rp, leray_param);


   nl.compute_initial_state();

   if(nl.verbose()){
      std::cout << "Solving the problem ..."  << '\n';
   }

   SolverInfo solve_info = nl.compute(load, solution);

   error_type error;
   error.h = average_diameter(msh);
   error.degree = rp.m_grad_degree;
   error.nb_dof = nl.getDofs();
   error.error_depl = nl.compute_l2_error(solution);
   error.error_grad = nl.compute_l2_gradient_error(gradient);

   return error;
}


template<template<typename, size_t , typename> class Mesh,
typename T, typename Storage>
error_type
run_leraylions_solver(const Mesh<T, 3, Storage>& msh, const ParamRun<T>& rp, const T leray_param)
{
   typedef T result_type;
   typedef static_vector<T, 3> result_grad_type;

   auto load = [leray_param](const point<T,3>& pt) -> auto {
      const T p = leray_param;
      T fx = -(std::pow(3.0, p/2.0) *(p-1)*exp((p-1)*(pt.x() + pt.y() + pt.z())));
      return result_type{fx};
   };

   auto solution = [](const point<T,3>& pt) -> auto {
      T fx = exp(pt.x() + pt.y() + pt.z());

      return result_type{fx};
   };


   auto gradient = [](const point<T,3>& pt) -> result_grad_type {
      result_grad_type grad = result_grad_type::Zero();

      grad(0) = exp(pt.x() + pt.y() + pt.z());
      grad(1) = exp(pt.x() + pt.y() + pt.z());
      grad(2) = exp(pt.x() + pt.y() + pt.z());

      return grad;
   };

   leraylions_solver<Mesh, T, 3, Storage,  point<T, 3> > nl(msh, rp, leray_param);

   nl.compute_initial_state();

   if(nl.verbose()){
      std::cout << "Solving the problem ..." << '\n';
   }

   SolverInfo solve_info = nl.compute(load, solution);

   error_type error;
   error.h = average_diameter(msh);
   error.degree = rp.m_grad_degree;
   error.nb_dof = nl.getDofs();
   error.error_depl = nl.compute_l2_error(solution);
   error.error_grad = nl.compute_l2_gradient_error(gradient);

   return error;
}

void
printResults(const std::vector<error_type>& error)
{
   if(error.size() > 0){
      std::ios::fmtflags f( std::cout.flags() );
      std::cout.precision(4);
      std::cout.setf(std::iostream::scientific, std::iostream::floatfield);

      std::cout << "Convergence test for k = " << error[0].degree << std::endl;
      std::cout << "-----------------------------------------------------------------------------------" << std::endl;
      std::cout << "| Size mesh  | Displacement | Convergence |  Gradient  | Convergence |    Total   |" << std::endl;
      std::cout << "|    h       |   L2 error   |     rate    |  L2 error  |     rate    | faces DOF  |" << std::endl;
      std::cout << "-----------------------------------------------------------------------------------" << std::endl;


      std::string s_dof = " " + std::to_string(error[0].nb_dof) + "                  ";
      s_dof.resize(10);

      std::cout << "| " <<  error[0].h << " |  " << error[0].error_depl << "  | " << "     -     " << " | " <<
      error[0].error_grad <<  " | " << "     -     "  <<  " | " << s_dof  <<  " |" << std::endl;

      for(size_t i = 1; i < error.size(); i++){
         s_dof = " " + std::to_string(error[i].nb_dof) + "                  ";
         s_dof.resize(10);
         double rate_depl = (log10(error[i-1].error_depl) - log10(error[i].error_depl))/(log10(error[i-1].h) - log10(error[i].h));
         double rate_grad = (log10(error[i-1].error_grad) - log10(error[i].error_grad))/(log10(error[i-1].h) - log10(error[i].h));

         std::cout << "| " <<  error[i].h << " |  " << error[i].error_depl << "  |  " << rate_depl << " | " <<
         error[i].error_grad <<  " |  " << rate_grad  <<  " | " << s_dof <<  " |" << std::endl;
      }

      std::cout << "-----------------------------------------------------------------------------------" << std::endl;
      std::cout << "  " <<std::endl;
      std::cout.flags( f );
   }
   else
      std::cout << "The file error is empty" << std::endl;
}


template< typename T>
void test_triangles_fvca5(const ParamRun<T>& rp, const T leray_param)
{
   size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_1.typ1");
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_2.typ1");
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_3.typ1");
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_4.typ1");
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_5.typ1");

   std::vector<error_type> error_sumup;

   for(size_t i = 0; i < runs; i++){
      auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_leraylions_solver(msh, rp, leray_param));
   }
   printResults(error_sumup);
}


template< typename T>
void test_triangles_netgen(const ParamRun<T>& rp, const T leray_param)
{
   size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri01.mesh2d");
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri02.mesh2d");
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri03.mesh2d");
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri04.mesh2d");
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri05.mesh2d");

   std::vector<error_type> error_sumup;

   for(size_t i = 0; i < runs; i++){
      auto msh = disk::load_netgen_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_leraylions_solver(msh, rp, leray_param));
   }
   printResults(error_sumup);
}



template< typename T>
void test_hexagons(const ParamRun<T>& rp, const T leray_param)
{
   size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_1.typ1");
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_2.typ1");
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_3.typ1");
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_4.typ1");
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_5.typ1");

   std::vector<error_type> error_sumup;

   for(size_t i = 0; i < runs; i++){
      auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_leraylions_solver(msh, rp, leray_param));
   }
   printResults(error_sumup);
}


template< typename T>
void test_kershaws(const ParamRun<T>& rp, const T leray_param)
{
   size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_1.typ1");
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_2.typ1");
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_3.typ1");
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_4.typ1");
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_5.typ1");

   std::vector<error_type> error_sumup;

   for(size_t i = 0; i < runs; i++){
      auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_leraylions_solver(msh, rp, leray_param));
   }
   printResults(error_sumup);
}


template< typename T>
void test_quads_fvca5(const ParamRun<T>& rp, const T leray_param)
{
   size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_1.typ1");
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_2.typ1");
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_3.typ1");
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_4.typ1");
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_5.typ1");

   std::vector<error_type> error_sumup;

   for(size_t i = 0; i < runs; i++){
      auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_leraylions_solver(msh, rp, leray_param));
   }
   printResults(error_sumup);
}


template< typename T>
void test_quads_diskpp(const ParamRun<T>& rp, const T leray_param)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-4-4.quad");
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-8-8.quad");
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-16-16.quad");
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-32-32.quad");
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-256-256.quad");

   std::vector<error_type> error_sumup;

   for(size_t i = 0; i < runs; i++){
      auto msh = disk::load_cartesian_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_leraylions_solver(msh, rp, leray_param));
   }
   printResults(error_sumup);
}

template< typename T>
void test_hexahedra_diskpp(const ParamRun<T>& rp, const T leray_param)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-2-2-2.hex");
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-4-4-4.hex");
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-8-8-8.hex");
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-16-16-16.hex");
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-32-32-32.hex");

   std::vector<error_type> error_sumup;

   for(int i = 0; i < runs; i++){
      auto msh = disk::load_cartesian_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_leraylions_solver(msh, rp, leray_param));
   }
   printResults(error_sumup);
}


template< typename T>
void test_hexahedra_fvca6(const ParamRun<T>& rp, const T leray_param)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_2x2x2.msh");
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_4x4x4.msh");
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_8x8x8.msh");
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_16x16x16.msh");
   //paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_32x32x32.hex");

   std::vector<error_type> error_sumup;

   for(int i = 0; i < runs; i++){
      auto msh = disk::load_fvca6_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_leraylions_solver(msh, rp, leray_param));
   }
   printResults(error_sumup);
}

template< typename T>
void test_tetrahedra_netgen(const ParamRun<T>& rp, const T leray_param)
{
   size_t runs = 3;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/tetra01.mesh");
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/tetra02.mesh");
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/tetra03.mesh");
   //paths.push_back("../diskpp/meshes/3D_tetras/netgen/tetra04.mesh");

   std::vector<error_type> error_sumup;

   for(int i = 0; i < runs; i++){
      auto msh = disk::load_netgen_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_leraylions_solver(msh, rp, leray_param));
   }
   printResults(error_sumup);
}


template< typename T>
void test_polyhedra_fvca6(const ParamRun<T>& rp, const T leray_param)
{
   size_t runs = 3;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_10.msh");
   paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_20.msh");
   paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_30.msh");
   //paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_40.msh");

   std::vector<error_type> error_sumup;

   for(int i = 0; i < runs; i++){
      auto msh = disk::load_fvca6_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_leraylions_solver(msh, rp, leray_param));
   }
   printResults(error_sumup);
}


template< typename T>
void test_tetrahedra_fvca6(const ParamRun<T>& rp, const T leray_param)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.1.msh");
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.2.msh");
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.3.msh");
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.4.msh");

   std::vector<error_type> error_sumup;

   for(int i = 0; i < runs; i++){
      auto msh = disk::load_fvca6_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_leraylions_solver(msh, rp, leray_param));
   }
   printResults(error_sumup);
}


int main(int argc, char **argv)
{
   using RealType = double;

   char    *mesh_filename  = nullptr;
   char    *plot_filename  = nullptr;
   int     degree          = 1;
   int     l               = 0;
   int     n_time_step     = 1;
   int     sublevel        = 1;
   bool three_dimensions = false;
   RealType leray_param = 2;

   ParamRun<RealType> rp;
   rp.m_sublevel = 4;
   rp.m_epsilon = 1E-11;



   int ch;

   while ( (ch = getopt(argc, argv, "23p:r:v")) != -1 )
   {
      switch(ch)
      {
         case '2': three_dimensions = false; break;

         case '3': three_dimensions = true; break;

         case 'p':
            leray_param = atof(optarg);
            break;

         case 'r':
            if(!rp.readParameters(optarg))
            exit(1);

            break;

         case 'v':
             rp.m_verbose = true;
             break;
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

   std::cout << " Test convergence rates for: "<< std::endl;
   std::cout << " ** Face_Degree = " << rp.m_face_degree << std::endl;
   std::cout << " ** Cell_Degree  = " << rp.m_cell_degree << std::endl;
   std::cout << " ** Grad_Degree  = " << rp.m_grad_degree << std::endl;
   std::cout << " ** Leray Parameters = " << leray_param << std::endl;
   std::cout << " "<< std::endl;

   if(three_dimensions){
      tc.tic();
      std::cout << "-Tetrahedras fvca6:" << std::endl;
      test_tetrahedra_fvca6<RealType>(rp, leray_param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout <<  "-Tetrahedras netgen:" << std::endl;
      test_tetrahedra_netgen<RealType>(rp, leray_param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout << "-Hexahedras fvca6:"  << std::endl;
      test_hexahedra_fvca6<RealType>(rp, leray_param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout << "-Hexahedras diskpp:"  << std::endl;
      test_hexahedra_diskpp<RealType>(rp, leray_param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;


      tc.tic();
      std::cout << "-Polyhedra:"  << std::endl;
      test_polyhedra_fvca6<RealType>(rp, leray_param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;
   }
   else{

      tc.tic();
      std::cout << "-Triangles fvca5:" << std::endl;
      test_triangles_fvca5<RealType>(rp, leray_param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout <<  "-Triangles netgen:" << std::endl;
      test_triangles_netgen<RealType>(rp, leray_param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout << "-Quadrangles fvca5:"  << std::endl;
      test_quads_fvca5<RealType>(rp, leray_param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout << "-Quadrangles diskpp:"  << std::endl;
      test_quads_diskpp<RealType>(rp, leray_param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;


      tc.tic();
      std::cout << "-Hexagons:"  << std::endl;
      test_hexagons<RealType>(rp, leray_param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;


      tc.tic();
      std::cout << "-Kershaws:"  << std::endl;
      test_kershaws<RealType>(rp, leray_param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

   }

}
