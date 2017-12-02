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

#include "config.h"

#ifdef HAVE_SOLVER_WRAPPERS
#include "agmg/agmg.hpp"
#endif

#include "loaders/loader.hpp"
#include "hho/hho.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "vector_laplacian3_solver.hpp"

template< typename T>
class conditionning_type
{
public:
   size_t   degree;
   size_t   nb_dof;
   T        h;
   T        num_cond;
   T        num_cond_full;

   conditionning_type() : degree(0), nb_dof(0), h(static_cast<T>(0)),
   num_cond(static_cast<T>(0)), num_cond_full(static_cast<T>(0))
   {}
};

struct run_params
{
   size_t  degree;
   int     l;
   bool    verbose;
};


void
usage(const char *progname)
{
   printf("Usage: %s <options> <filename>\n\n", progname);
   printf("    -2: test 2D mesh (default)\n");
   printf("    -3: test 3D mesh\n");
   printf("    -k: face degree (>=0)\n");
   printf("    -l: difference beetween cell and face degree (-1 <= l <= 1) \n");
}


template<template<typename, size_t , typename> class Mesh,
typename T, typename Storage>
conditionning_type<T>
run_vector_laplacian_solver(const Mesh<T, 2, Storage>& msh, const run_params& rp, LaplacianParameters material_data)
{
   typedef Mesh<T, 2, Storage> mesh_type;
   typedef static_vector<T, 2> result_type;
   typedef static_matrix<T, 2, 2> result_grad_type;

   timecounter tc;
   tc.tic();

   auto load = [material_data](const point<T,2>& p) -> result_type {
      T fx = 2.*material_data.lambda*M_PI*M_PI*sin(M_PI*p.x())*sin(M_PI*p.y());
      T fy = 2.*material_data.lambda*M_PI*M_PI*cos(M_PI*p.x())*cos(M_PI*p.y());

      return result_type{fx,fy};
   };

   auto solution = [material_data](const point<T,2>& p) -> result_type {
      T fx = sin(M_PI*p.x())*sin(M_PI*p.y());
      T fy = cos(M_PI*p.x())*cos(M_PI*p.y());

      return result_type{fx,fy};
   };

   auto gradient = [material_data](const point<T,2>& p) -> result_grad_type {
      result_grad_type grad = result_grad_type::Zero();

      grad(0,0) = cos(M_PI*p.x())*sin(M_PI*p.y());
      grad(0,1) = sin(M_PI*p.x())*cos(M_PI*p.y());
      grad(1,0) = -sin(M_PI*p.x())*cos(M_PI*p.y());
      grad(1,1) = -cos(M_PI*p.x())*sin(M_PI*p.y());
      return M_PI*grad;
   };


   vector_laplacian_solver<mesh_type> vl(msh, rp.degree, rp.l);
   vl.verbose(false);

   vl.changeLaplacianParameters(material_data);

   const assembly_info assembling_info = vl.assemble(load, solution);
   const auto cond_numbers = vl.conditioning();

   conditionning_type<T> cond;
   cond.h = average_diameter(msh);
   cond.degree = rp.degree;
   cond.nb_dof = vl.getDofs();
   cond.num_cond = cond_numbers.first;
   cond.num_cond_full = cond_numbers.second;

   return cond;
}


template<template<typename, size_t , typename> class Mesh,
typename T, typename Storage>
conditionning_type<T>
run_vector_laplacian_solver(const Mesh<T, 3, Storage>& msh, const run_params& rp, LaplacianParameters material_data)
{
   typedef Mesh<T, 3, Storage> mesh_type;
   typedef static_vector<T, 3> result_type;
   typedef static_matrix<T, 3, 3> result_grad_type;

   timecounter tc;
   tc.tic();

   auto load = [material_data](const point<T,3>& p) -> auto {
      T fx = 2.*material_data.lambda*M_PI*M_PI *cos(M_PI*p.x())*sin(M_PI*p.y());

      T fy = 2.*material_data.lambda*M_PI*M_PI *cos(M_PI*p.y())*sin(M_PI*p.z());

      T fz = 2.*material_data.lambda*M_PI*M_PI *cos(M_PI*p.z())*sin(M_PI*p.x());

      return result_type{fx,fy,fz};
   };

   auto solution = [material_data](const point<T,3>& p) -> auto {
      T fx = cos(M_PI*p.x())*sin(M_PI*p.y());
      T fy = cos(M_PI*p.y())*sin(M_PI*p.z());
      T fz = cos(M_PI*p.z())*sin(M_PI*p.x());
      return result_type{fx,fy,fz};
   };


   auto gradient = [material_data](const point<T,3  >& p) -> result_grad_type {
      result_grad_type grad = result_grad_type::Zero();

      grad(0,0) = -sin(M_PI*p.x())*sin(M_PI*p.y()); grad(0,1) = cos(M_PI*p.x())*cos(M_PI*p.y());
      grad(1,1) = -sin(M_PI*p.y())*sin(M_PI*p.z()); grad(1,2) = cos(M_PI*p.y())*cos(M_PI*p.z());
      grad(2,2) = -sin(M_PI*p.z())*sin(M_PI*p.x()); grad(2,0) = cos(M_PI*p.z())*cos(M_PI*p.x());
      return M_PI*grad;
   };

   vector_laplacian_solver<mesh_type> vl(msh, rp.degree, rp.l);
   vl.verbose(false);

   vl.changeLaplacianParameters(material_data);

   const assembly_info assembling_info = vl.assemble(load, solution);
   const auto cond_numbers = vl.conditioning();

   conditionning_type<T> cond;
   cond.h = average_diameter(msh);
   cond.degree = rp.degree;
   cond.nb_dof = vl.getDofs();
   cond.num_cond = cond_numbers.first;
   cond.num_cond_full = cond_numbers.second;

   return cond;
}

template< typename T>
void
printResults(const std::vector<conditionning_type<T> >& cond)
{
   if(cond.size() > 0){
      std::ios::fmtflags f( std::cout.flags() );
      std::cout.precision(4);
      std::cout.setf(std::iostream::scientific, std::iostream::floatfield);

      std::cout << "Conditioning test for k = " << cond[0].degree << std::endl;
      std::cout << "---------------------------------------------------------------------" << std::endl;
      std::cout << "| Size mesh  |    Conditioning   |    Conditioning     |    Total   |" << std::endl;
      std::cout << "|    h       | with condensation | without condensation| faces DOF  |" << std::endl;
      std::cout << "---------------------------------------------------------------------" << std::endl;

      std::string s_dof;
      for(size_t i = 0; i < cond.size(); i++){
         s_dof = " " + std::to_string(cond[i].nb_dof) + "                  ";
         s_dof.resize(10);
         std::cout << "| " <<  cond[i].h << " |  " << cond[i].num_cond << "  |  " << cond[i].num_cond_full
          << " | " << s_dof <<  " |" << std::endl;
      }

      std::cout << "---------------------------------------------------------------------" << std::endl;
      std::cout << "  " <<std::endl;
      std::cout.flags( f );
   }
   else
      std::cout << "The file is empty" << std::endl;
}


template< typename T>
void test_triangles_fvca5(const run_params& rp, LaplacianParameters material_data)
{
   const size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_1.typ1");
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_2.typ1");
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_3.typ1");
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_4.typ1");
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_5.typ1");

   std::vector<conditionning_type<T> > cond_sumup;

   for(size_t i = 0; i < runs; i++){
      const auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      cond_sumup.push_back(run_vector_laplacian_solver(msh, rp, material_data));
   }
   printResults(cond_sumup);
}


template< typename T>
void test_triangles_netgen(const run_params& rp, LaplacianParameters material_data)
{
   const size_t runs = 3;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri01.mesh2d");
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri02.mesh2d");
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri03.mesh2d");
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri04.mesh2d");
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri05.mesh2d");

   std::vector<conditionning_type<T> > cond_sumup;

   for(size_t i = 0; i < runs; i++){
      const auto msh = disk::load_netgen_2d_mesh<T>(paths[i].c_str());
      cond_sumup.push_back(run_vector_laplacian_solver(msh, rp, material_data));
   }
   printResults(cond_sumup);
}



template< typename T>
void test_hexagons(const run_params& rp, LaplacianParameters material_data)
{
   const size_t runs = 3;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_1.typ1");
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_2.typ1");
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_3.typ1");
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_4.typ1");
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_5.typ1");

   std::vector<conditionning_type<T> > cond_sumup;

   for(size_t i = 0; i < runs; i++){
      const auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      cond_sumup.push_back(run_vector_laplacian_solver(msh, rp, material_data));
   }
   printResults(cond_sumup);
}


template< typename T>
void test_kershaws(const run_params& rp, LaplacianParameters material_data)
{
   const size_t runs = 3;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_1.typ1");
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_2.typ1");
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_3.typ1");
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_4.typ1");
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_5.typ1");

   std::vector<conditionning_type<T> > cond_sumup;

   for(size_t i = 0; i < runs; i++){
      const auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      cond_sumup.push_back(run_vector_laplacian_solver(msh, rp, material_data));
   }
   printResults(cond_sumup);
}


template< typename T>
void test_quads_fvca5(const run_params& rp, LaplacianParameters material_data)
{
   const size_t runs = 3;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_1.typ1");
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_2.typ1");
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_3.typ1");
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_4.typ1");
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_5.typ1");

   std::vector<conditionning_type<T> > cond_sumup;

   for(size_t i = 0; i < runs; i++){
      const auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      cond_sumup.push_back(run_vector_laplacian_solver(msh, rp, material_data));
   }
   printResults(cond_sumup);
}


template< typename T>
void test_quads_diskpp(const run_params& rp, LaplacianParameters material_data)
{
   const size_t runs = 3;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-4-4.quad");
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-8-8.quad");
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-16-16.quad");
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-32-32.quad");
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-256-256.quad");

   std::vector<conditionning_type<T> > cond_sumup;

   for(size_t i = 0; i < runs; i++){
      auto msh = disk::load_cartesian_2d_mesh<T>(paths[i].c_str());
      cond_sumup.push_back(run_vector_laplacian_solver(msh, rp, material_data));
   }
   printResults(cond_sumup);
}

template< typename T>
void test_hexahedra_diskpp(const run_params& rp, LaplacianParameters material_data)
{
   size_t runs = 2;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-2-2-2.hex");
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-4-4-4.hex");
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-8-8-8.hex");
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-16-16-16.hex");
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-32-32-32.hex");

   std::vector<conditionning_type<T> > cond_sumup;

   for(int i = 0; i < runs; i++){
      const auto msh = disk::load_cartesian_3d_mesh<T>(paths[i].c_str());
      cond_sumup.push_back(run_vector_laplacian_solver(msh, rp, material_data));
   }
   printResults(cond_sumup);
}


template< typename T>
void test_hexahedra_fvca6(const run_params& rp, LaplacianParameters material_data)
{
   const size_t runs = 2;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_2x2x2.msh");
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_4x4x4.msh");
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_8x8x8.msh");
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_16x16x16.msh");
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_32x32x32.hex");

   std::vector<conditionning_type<T> > cond_sumup;

   for(int i = 0; i < runs; i++){
      const auto msh = disk::load_fvca6_3d_mesh<T>(paths[i].c_str());
      cond_sumup.push_back(run_vector_laplacian_solver(msh, rp, material_data));
   }
   printResults(cond_sumup);
}

template< typename T>
void test_tetrahedra_netgen(const run_params& rp, LaplacianParameters material_data)
{
   const size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet0.mesh");
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet1.mesh");
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet2.mesh");
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet3.mesh");
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet4.mesh");

   std::vector<conditionning_type<T> > cond_sumup;

   for(int i = 0; i < runs; i++){
      const auto msh = disk::load_netgen_3d_mesh<T>(paths[i].c_str());
      cond_sumup.push_back(run_vector_laplacian_solver(msh, rp, material_data));
   }
   printResults(cond_sumup);
}


template< typename T>
void test_polyhedra_fvca6(const run_params& rp, LaplacianParameters material_data)
{
   const size_t runs = 2;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_10.msh");
   paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_20.msh");
   paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_30.msh");
   paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_40.msh");

   std::vector<conditionning_type<T> > cond_sumup;

   for(int i = 0; i < runs; i++){
      const auto msh = disk::load_fvca6_3d_mesh<T>(paths[i].c_str());
      cond_sumup.push_back(run_vector_laplacian_solver(msh, rp, material_data));
   }
   printResults(cond_sumup);
}


template< typename T>
void test_tetrahedra_fvca6(const run_params& rp, LaplacianParameters material_data)
{
   const size_t runs = 2;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.1.msh");
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.2.msh");
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.3.msh");
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.4.msh");

   std::vector<conditionning_type<T> > cond_sumup;

   for(int i = 0; i < runs; i++){
      const auto msh = disk::load_fvca6_3d_mesh<T>(paths[i].c_str());
      cond_sumup.push_back(run_vector_laplacian_solver(msh, rp, material_data));
   }
   printResults(cond_sumup);
}


int main(int argc, char **argv)
{
   using RealType = double;

   int     degree          = 1;
   int     l               = 0;
   bool three_dimensions = false;

   run_params rp;
   rp.degree   = 1;
   rp.l        = 0;
   rp.verbose  = false;

   LaplacianParameters material_data;

   material_data.lambda = 1.0;

   int ch;

   while ( (ch = getopt(argc, argv, "23k:l:v")) != -1 )
   {
      switch(ch)
      {
         case '2': three_dimensions = false; break;

         case '3': three_dimensions = true; break;

         case 'k':
            degree = atoi(optarg);
            if (degree < 0)
            {
               std::cout << "Degree must be positive. Falling back to 1." << std::endl;
               degree = 1;
            }
            rp.degree = degree;
            break;

         case 'l':
            l = atoi(optarg);
            if (l < -1 or l > 1)
            {
               std::cout << "l can be -1, 0 or 1. Falling back to 0." << std::endl;
               l = 0;
            }
            rp.l = l;
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

   std::cout << " Test conditionning for: "<< std::endl;
   std::cout << " ** Face_Degree = " << rp.degree << std::endl;
   std::cout << " ** Cell_Degree  = " << rp.degree + rp.l << std::endl;
   std::cout << " "<< std::endl;

   if(three_dimensions){
      tc.tic();
      std::cout << "-Tetrahedras fvca6:" << std::endl;
      //test_tetrahedra_fvca6<RealType>(rp, material_data);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout <<  "-Tetrahedras netgen:" << std::endl;
      test_tetrahedra_netgen<RealType>(rp, material_data);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout << "-Hexahedras fvca6:"  << std::endl;
      test_hexahedra_fvca6<RealType>(rp, material_data);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout << "-Hexahedras diskpp:"  << std::endl;
      test_hexahedra_diskpp<RealType>(rp, material_data);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout << "-Polyhedra:"  << std::endl;
      test_polyhedra_fvca6<RealType>(rp, material_data);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;
   }
   else{
      tc.tic();
      std::cout << "-Triangles fvca5:" << std::endl;
      test_triangles_fvca5<RealType>(rp, material_data);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout <<  "-Triangles netgen:" << std::endl;
      test_triangles_netgen<RealType>(rp, material_data);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout << "-Quadrangles fvca5:"  << std::endl;
      test_quads_fvca5<RealType>(rp, material_data);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout << "-Quadrangles diskpp:"  << std::endl;
      test_quads_diskpp<RealType>(rp, material_data);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout << "-Hexagons:"  << std::endl;
      test_hexagons<RealType>(rp, material_data);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout << "-Kershaws:"  << std::endl;
      test_kershaws<RealType>(rp, material_data);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;
   }
}
