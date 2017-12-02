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

#include "loaders/loader.hpp"
#include "hho/hho.hpp"
#include "geometry/geometry.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "vector_laplacian3_solver.hpp"

struct run_params
{
   size_t  degree;
   int     l;
   bool    verbose;
};


template<template<typename, size_t , typename> class Mesh,
typename T, typename Storage>
void
run_vector_laplacian_solver(const Mesh<T, 2, Storage>& msh, run_params& rp, LaplacianParameters material_data)
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
      T fx = material_data.lambda*sin(M_PI*p.x())*sin(M_PI*p.y());
      T fy = material_data.lambda*cos(M_PI*p.x())*cos(M_PI*p.y());

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
   vl.verbose(rp.verbose);

   vl.changeLaplacianParameters(material_data);

   assembly_info assembling_info = vl.assemble(load, solution);

   if(vl.verbose()){
      std::cout << "Assembling: " << assembling_info.time_assembly << " sec"  << '\n';
      std::cout << "Conditionning: " << std::endl;
   }

   const auto cond_numbers = vl.conditioning();

   if(vl.verbose()){
      std::cout << "*********************************************" << std::endl;
      std::cout << "* conditioning_number: static condensation *" << std::endl;
      std::cout << "*         with        ***      without      *" << std::endl;
      std::cout << "*********************************************" << std::endl;
      std::cout << "*    " << cond_numbers.first << "    ***    " << cond_numbers.second
      << "    *" << std::endl;
      std::cout << "*********************************************" << std::endl;
   }

   if(vl.verbose()){
      solver_info solve_info = vl.solve(true);
      //postprocess_info post_info = vl.postprocess(load);
      postprocess_info post_info = vl.postprocess_full();
      std::cout << "Discetisation h: " << disk::mesh_h(msh) << std::endl;
      std::cout << "L2 error: " << vl.compute_l2_error(solution) << std::endl;
      std::cout << "L2 gradient error: " << vl.compute_l2_gradient_error(gradient) << std::endl;
   }

}

template<template<typename, size_t, typename> class Mesh,
typename T, typename Storage>
void
run_vector_laplacian_solver(const Mesh<T, 3, Storage>& msh, run_params& rp, LaplacianParameters material_data)
{
   typedef Mesh<T, 3, Storage> mesh_type;
   typedef static_vector<T, 3> result_type;
   typedef static_matrix<T, 3, 3> result_grad_type;

   timecounter tc;
   tc.tic();

   auto load = [material_data](const point<T,3>& p) -> auto {
      const T lambda = material_data.lambda;

      T fx = -2.0 * lambda * p.y() * p.z();

      T fy = -2.0 * lambda * p.z() * p.x();

      T fz = -2.0 * lambda * p.x() * p.y();

      return result_type{fx,fy,fz};
   };

   auto solution = [material_data](const point<T,3>& p) -> auto {
      T fx = p.x()*p.x()*p.y()*p.z();
      T fy = p.y()*p.x()*p.y()*p.z();
      T fz = p.z()*p.x()*p.y()*p.z();
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
   vl.verbose(rp.verbose);

   vl.changeLaplacianParameters(material_data);

   assembly_info assembling_info = vl.assemble(load, solution);

   if(vl.verbose()){
      std::cout << "Assembling: " << assembling_info.time_assembly << " sec"  << '\n';
   }

   const auto cond_numbers = vl.conditioning();

   if(vl.verbose()){
      std::cout << "*********************************************" << std::endl;
      std::cout << "* conditioning_number: static condensation *" << std::endl;
      std::cout << "*         with        ***      without      *" << std::endl;
      std::cout << "*********************************************" << std::endl;
      std::cout << "*    " << cond_numbers.first << "    ***    " << cond_numbers.second
      << "    *" << std::endl;
      std::cout << "*********************************************" << std::endl;
   }

   if(vl.verbose()){
      solver_info solve_info = vl.solve(true);
      //postprocess_info post_info = vl.postprocess(load);
      postprocess_info post_info = vl.postprocess_full();
      std::cout << "Discetisation h: " << disk::mesh_h(msh) << std::endl;
      std::cout << "L2 error: " << vl.compute_l2_error(solution) << std::endl;
      std::cout << "L2 gradient error: " << vl.compute_l2_gradient_error(gradient) << std::endl;
   }
}


int main(int argc, char **argv)
{
   using RealType = double;

   char    *mesh_filename  = nullptr;
   int     degree          = 1;

   run_params rp;
   rp.degree   = 1;
   rp.l        = 0;
   rp.verbose  = true;

   //Elasticity Parameters
   LaplacianParameters material_data;

   material_data.lambda = 1.0;

   int ch;

   while ( (ch = getopt(argc, argv, "k:l:n:p:v")) != -1 )
   {
      switch(ch)
      {
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
            rp.l = atoi(optarg);
            if (rp.l < -1 or rp.l > 1)
            {
               std::cout << "l can be -1, 0 or 1. Falling back to 0." << std::endl;
               rp.l = 0;
            }
            break;

         case 'v':
            rp.verbose = true;
            break;

         case 'h':
         case '?':
         default:
            std::cout << "wrong arguments" << std::endl;
            exit(1);
      }
   }

   argc -= optind;
   argv += optind;

   if (argc == 0)
   {
      std::cout << "Mesh format: 1D uniform (Not avaible)" << std::endl;
      //auto msh = disk::load_uniform_1d_mesh<RealType>(0, 1, elems_1d);
      //run_NL_elasticity_solver(msh, rp);
      return 0;
   }

   mesh_filename = argv[0];

   /* FVCA5 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
   {
      std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
      auto msh = disk::load_fvca5_2d_mesh<RealType>(mesh_filename);
      run_vector_laplacian_solver(msh, rp, material_data);
      return 0;
   }

   /* Netgen 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$") ))
   {
      std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
      auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);
      run_vector_laplacian_solver(msh, rp, material_data);
      return 0;
   }

   /* DiSk++ cartesian 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.quad$") ))
   {
      std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
      auto msh = disk::load_cartesian_2d_mesh<RealType>(mesh_filename);
      run_vector_laplacian_solver(msh, rp, material_data);
      return 0;
   }

   /* Medit 2d*/
   if (std::regex_match(mesh_filename, std::regex(".*\\.medit2d$") ))
   {
      std::cout << "Guessed mesh format: Medit format" << std::endl;
      auto msh = disk::load_medit_2d_mesh<RealType>(mesh_filename);
      return 0;
   }

   /* Netgen 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
   {
      std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
      auto msh = disk::load_netgen_3d_mesh<RealType>(mesh_filename);
      run_vector_laplacian_solver(msh, rp, material_data);
      return 0;
   }

   /* DiSk++ cartesian 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
   {
      std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
      auto msh = disk::load_cartesian_3d_mesh<RealType>(mesh_filename);
      run_vector_laplacian_solver(msh, rp, material_data);
      return 0;
   }

   /* DiSk++ cartesian 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.msh$") ))
   {
      std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
      auto msh = disk::load_fvca6_3d_mesh<RealType>(mesh_filename);
      run_vector_laplacian_solver(msh, rp, material_data);
      return 0;
   }

   std::cout << "Unkwnon mesh format" << std::endl;
   return 0;

}
