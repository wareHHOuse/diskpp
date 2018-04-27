/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
 * Nicolas Pignet  (C) 2018                     nicolas.pignet@enpc.fr
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

#include "colormanip.h"

#include "loaders/loader.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "linear_elasticity_solver.hpp"

struct error_type
{
   size_t degree;
   size_t nb_dof;
   double h;
   double error_depl;
   double error_stress;
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
   printf("    -2: test 2D mesh (default)\n");
   printf("    -3: test 3D mesh\n");
   printf("    -k: face degree (>=0)\n");
   printf("    -l: difference beetween cell and face degree (-1 <= l <= 1) \n");
   printf("    -v: verbose\n");
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
error_type
run_linear_elasticity_solver(const Mesh<T, 2, Storage>& msh,
                             const run_params&          rp,
                             const ElasticityParameters material_data)
{
   typedef Mesh<T, 2, Storage>                            mesh_type;
   typedef static_vector<T, 2>                            result_type;
   typedef static_matrix<T, 2, 2>                         grad_type;
   typedef disk::mechanics::BoundaryConditions<mesh_type> Bnd_type;

   timecounter tc;
   tc.tic();

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

      return -M_PI * M_PI / (lambda + 1) * result_type{fx, fy};
   };

   auto solution = [material_data](const point<T, 2>& p) -> result_type {
      T fx = sin(2 * M_PI * p.y()) * (cos(2 * M_PI * p.x()) - 1) +
             1.0 / (1 + material_data.lambda) * sin(M_PI * p.x()) * sin(M_PI * p.y());
      T fy = -sin(2 * M_PI * p.x()) * (cos(2 * M_PI * p.y()) - 1) +
             1.0 / (1 + material_data.lambda) * sin(M_PI * p.x()) * sin(M_PI * p.y());

      return result_type{fx, fy};
   };

   auto sigma = [material_data](const point<T, 2>& p) -> grad_type {

      const T lambda = material_data.lambda;
      const T mu     = material_data.mu;

      T g11 = -(2 * (lambda + 1) * sin(2 * M_PI * p.x()) * sin(2 * M_PI * p.y()) -
                sin(M_PI * p.y()) * cos(M_PI * p.x()));
      T g12 = (2 * lambda + 2) * (cos(2 * M_PI * p.x()) - 1) * cos(2 * M_PI * p.y()) +
              sin(M_PI * p.x()) * cos(M_PI * p.y());

      T g21 = (-2 * lambda + 2) * (cos(2 * M_PI * p.y()) - 1) * cos(2 * M_PI * p.x()) +
              sin(M_PI * p.y()) * cos(M_PI * p.x());

      T g22 = 2 * (lambda + 1) * sin(2 * M_PI * p.x()) * sin(2 * M_PI * p.y()) +
              sin(M_PI * p.x()) * cos(M_PI * p.y());

      grad_type g = grad_type::Zero();

      g(0, 0) = g11;
      g(0, 1) = g12;
      g(1, 0) = g21;
      g(1, 1) = g22;

      g *= M_PI / (lambda + 1);

      const grad_type gs = 0.5 * (g + g.transpose());

      const T divu = gs.trace();

      return 2 * mu * gs + lambda * divu * static_matrix<T, 2, 2>::Identity();
   };

   Bnd_type bnd(msh);
   bnd.addDirichletEverywhere(solution);

   linear_elasticity_solver<mesh_type> le(msh, bnd, material_data, rp.degree, rp.l);
   le.verbose(rp.verbose);

   le.changeElasticityParameters(material_data);

   assembly_info    assembling_info = le.assemble(load);
   solver_info      solve_info      = le.solve();
   postprocess_info post_info       = le.postprocess(load);

   error_type error;
   error.h            = average_diameter(msh);
   error.degree       = rp.degree;
   error.nb_dof       = le.getDofs();
   error.error_depl   = le.compute_l2_displacement_error(solution);
   error.error_stress = le.compute_l2_stress_error(sigma);

   return error;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
error_type
run_linear_elasticity_solver(const Mesh<T, 3, Storage>&  msh,
                             const run_params&           rp,
                             const ElasticityParameters& material_data)
{
   typedef Mesh<T, 3, Storage>                            mesh_type;
   typedef static_vector<T, 3>                            result_type;
   typedef static_matrix<T, 3, 3>                         grad_type;
   typedef disk::mechanics::BoundaryConditions<mesh_type> Bnd_type;

   timecounter tc;
   tc.tic();

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

   auto sigma = [material_data](const point<T, 3>& p) -> grad_type {

      const T lambda = material_data.lambda;
      const T mu     = material_data.mu;

      T g11 = -(2 * (lambda + 1) * sin(2 * M_PI * p.x()) * sin(2 * M_PI * p.y()) -
                sin(M_PI * p.y()) * cos(M_PI * p.x()));
      T g12 = (2 * lambda + 2) * (cos(2 * M_PI * p.x()) - 1) * cos(2 * M_PI * p.y()) +
              sin(M_PI * p.x()) * cos(M_PI * p.y());

      T g21 = (-2 * lambda + 2) * (cos(2 * M_PI * p.y()) - 1) * cos(2 * M_PI * p.x()) +
              sin(M_PI * p.y()) * cos(M_PI * p.x());

      T g22 = 2 * (lambda + 1) * sin(2 * M_PI * p.x()) * sin(2 * M_PI * p.y()) +
              sin(M_PI * p.x()) * cos(M_PI * p.y());

      grad_type g = grad_type::Zero();

      g(0, 0) = g11;
      g(0, 1) = g12;
      g(1, 0) = g21;
      g(1, 1) = g22;

      g *= M_PI / (lambda + 1);

      const grad_type gs = 0.5 * (g + g.transpose());

      const T divu = gs.trace();

      return 2 * mu * gs + lambda * divu * static_matrix<T, 3, 3>::Identity();
   };

   Bnd_type bnd(msh);
   bnd.addDirichletEverywhere(solution);

   linear_elasticity_solver<mesh_type> le(msh, bnd, material_data, rp.degree, rp.l);
   le.verbose(rp.verbose);

   le.changeElasticityParameters(material_data);

   assembly_info    assembling_info = le.assemble(load);
   solver_info      solve_info      = le.solve();
   postprocess_info post_info       = le.postprocess(load);

   error_type error;
   error.h            = average_diameter(msh);
   error.degree       = rp.degree;
   error.nb_dof       = le.getDofs();
   error.error_depl   = le.compute_l2_displacement_error(solution);
   error.error_stress = le.compute_l2_stress_error(sigma);

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
        << "| Size mesh  | Displacement | Convergence |  Stress    | Convergence |    Total   |"
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
                << " | " << error[0].error_stress << " | "
                << "     -     "
                << " | " << s_dof << " |" << std::endl;

      for (size_t i = 1; i < error.size(); i++) {
         s_dof = " " + std::to_string(error[i].nb_dof) + "                  ";
         s_dof.resize(10);
         double rate_depl = (log10(error[i - 1].error_depl) - log10(error[i].error_depl)) /
                            (log10(error[i - 1].h) - log10(error[i].h));
         double rate_stress = (log10(error[i - 1].error_stress) - log10(error[i].error_stress)) /
                              (log10(error[i - 1].h) - log10(error[i].h));

         std::cout << "| " << error[i].h << " |  " << error[i].error_depl << "  |  " << rate_depl
                   << " | " << error[i].error_stress << " |  " << rate_stress << " | " << s_dof
                   << " |" << std::endl;
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
test_triangles_fvca5(const run_params& rp, const ElasticityParameters material_data)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_1.typ1");
   paths.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_2.typ1");
   paths.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_3.typ1");
   paths.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_4.typ1");
   paths.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_5.typ1");

   std::vector<error_type> error_sumup;

   for (size_t i = 0; i < runs; i++) {
      auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_linear_elasticity_solver(msh, rp, material_data));
   }
   printResults(error_sumup);
}

template<typename T>
void
test_triangles_netgen(const run_params& rp, const ElasticityParameters material_data)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri01.mesh2d");
   paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri02.mesh2d");
   paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri03.mesh2d");
   paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri04.mesh2d");
   paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri05.mesh2d");

   std::vector<error_type> error_sumup;

   for (size_t i = 0; i < runs; i++) {
      auto msh = disk::load_netgen_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_linear_elasticity_solver(msh, rp, material_data));
   }
   printResults(error_sumup);
}

template<typename T>
void
test_hexagons(const run_params& rp, const ElasticityParameters material_data)
{
   size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_1.typ1");
   paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_2.typ1");
   paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_3.typ1");
   paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_4.typ1");
   paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_5.typ1");

   std::vector<error_type> error_sumup;

   for (size_t i = 0; i < runs; i++) {
      auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_linear_elasticity_solver(msh, rp, material_data));
   }
   printResults(error_sumup);
}

template<typename T>
void
test_kershaws(const run_params& rp, const ElasticityParameters material_data)
{
   size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../../../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_1.typ1");
   paths.push_back("../../../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_2.typ1");
   paths.push_back("../../../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_3.typ1");
   paths.push_back("../../../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_4.typ1");
   paths.push_back("../../../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_5.typ1");

   std::vector<error_type> error_sumup;

   for (size_t i = 0; i < runs; i++) {
      auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_linear_elasticity_solver(msh, rp, material_data));
   }
   printResults(error_sumup);
}

template<typename T>
void
test_quads_fvca5(const run_params& rp, const ElasticityParameters material_data)
{
   size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_1.typ1");
   paths.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_2.typ1");
   paths.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_3.typ1");
   paths.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_4.typ1");
   paths.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_5.typ1");

   std::vector<error_type> error_sumup;

   for (size_t i = 0; i < runs; i++) {
      auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_linear_elasticity_solver(msh, rp, material_data));
   }
   printResults(error_sumup);
}

template<typename T>
void
test_quads_diskpp(const run_params& rp, const ElasticityParameters material_data)
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
      error_sumup.push_back(run_linear_elasticity_solver(msh, rp, material_data));
   }
   printResults(error_sumup);
}

// template<typename T>
// void
// test_hexahedra_diskpp(const run_params& rp, const ElasticityParameters material_data)
// {
//    size_t runs = 4;

//    std::vector<std::string> paths;
//    paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-2-2-2.hex");
//    paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-4-4-4.hex");
//    paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-8-8-8.hex");
//    paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-16-16-16.hex");
//    paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-32-32-32.hex");

//    std::vector<error_type> error_sumup;

//    for (int i = 0; i < runs; i++) {
//       auto msh = disk::load_cartesian_3d_mesh<T>(paths[i].c_str());
//       error_sumup.push_back(run_linear_elasticity_solver(msh, rp, material_data));
//    }
//    printResults(error_sumup);
// }

// template<typename T>
// void
// test_hexahedra_fvca6(const run_params& rp, const ElasticityParameters material_data)
// {
//    size_t runs = 4;

//    std::vector<std::string> paths;
//    paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_2x2x2.msh");
//    paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_4x4x4.msh");
//    paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_8x8x8.msh");
//    paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_16x16x16.msh");
//    paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_32x32x32.hex");

//    std::vector<error_type> error_sumup;

//    for (int i = 0; i < runs; i++) {
//       auto msh = disk::load_fvca6_3d_mesh<T>(paths[i].c_str());
//       error_sumup.push_back(run_linear_elasticity_solver(msh, rp, material_data));
//    }
//    printResults(error_sumup);
// }

// template<typename T>
// void
// test_tetrahedra_netgen(const run_params& rp, const ElasticityParameters material_data)
// {
//    size_t runs = 4;

//    std::vector<std::string> paths;
//    paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet0.mesh");
//    paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet1.mesh");
//    paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet2.mesh");
//    paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet3.mesh");
//    paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet4.mesh");

//    std::vector<error_type> error_sumup;

//    for (int i = 0; i < runs; i++) {
//       auto msh = disk::load_netgen_3d_mesh<T>(paths[i].c_str());
//       error_sumup.push_back(run_linear_elasticity_solver(msh, rp, material_data));
//    }
//    printResults(error_sumup);
// }

// template<typename T>
// void
// test_polyhedra_fvca6(const run_params& rp, const ElasticityParameters material_data)
// {
//    size_t runs = 3;

//    std::vector<std::string> paths;
//    paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_10.msh");
//    paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_20.msh");
//    paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_30.msh");
//    paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_40.msh");

//    std::vector<error_type> error_sumup;

//    for (int i = 0; i < runs; i++) {
//       auto msh = disk::load_fvca6_3d_mesh<T>(paths[i].c_str());
//       error_sumup.push_back(run_linear_elasticity_solver(msh, rp, material_data));
//    }
//    printResults(error_sumup);
// }

// template<typename T>
// void
// test_tetrahedra_fvca6(const run_params& rp, const ElasticityParameters material_data)
// {
//    size_t runs = 4;

//    std::vector<std::string> paths;
//    paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.0.msh");
//    paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.1.msh");
//    paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.2.msh");
//    paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.3.msh");
//    paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.4.msh");

//    std::vector<error_type> error_sumup;

//    for (int i = 0; i < runs; i++) {
//       auto msh = disk::load_fvca6_3d_mesh<T>(paths[i].c_str());
//       error_sumup.push_back(run_linear_elasticity_solver(msh, rp, material_data));
//    }
//    printResults(error_sumup);
// }

int
main(int argc, char** argv)
{
   using RealType = double;

   int    degree = 1;
   int    l      = 0;
   size_t dim    = 2;

   run_params rp;
   rp.degree  = 1;
   rp.l       = 0;
   rp.verbose = false;

   // Elasticity Parameters
   ElasticityParameters material_data;

   material_data.mu     = 1.0;
   material_data.lambda = 10E5;

   int ch;

   while ((ch = getopt(argc, argv, "23k:l:v")) != -1) {
      switch (ch) {
         case '2': dim = 2; break;
         case '3': dim = 3; break;

         case 'k':
            degree = atoi(optarg);
            if (degree < 0) {
               std::cout << "Degree must be positive. Falling back to 1." << std::endl;
               degree = 1;
            }
            rp.degree = degree;
            break;

         case 'l':
            l = atoi(optarg);
            if (l < -1 or l > 1) {
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

   std::cout << " Test convergence rates for: " << std::endl;
   std::cout << " ** Face_Degree = " << rp.degree << std::endl;
   std::cout << " ** Cell_Degree  = " << rp.degree + rp.l << std::endl;
   std::cout << " " << std::endl;

   if (dim == 3) {
      //       tc.tic();
      //       std::cout << "-Tetrahedras fvca6:" << std::endl;
      //       test_tetrahedra_fvca6<RealType>(rp, material_data);
      //       tc.toc();
      //       std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      //       std::cout << " " << std::endl;

      //       tc.tic();
      //       std::cout << "-Tetrahedras netgen:" << std::endl;
      //       test_tetrahedra_netgen<RealType>(rp, material_data);
      //       tc.toc();
      //       std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      //       std::cout << " " << std::endl;

      //       tc.tic();
      //       std::cout << "-Hexahedras fvca6:" << std::endl;
      //       test_hexahedra_fvca6<RealType>(rp, material_data);
      //       tc.toc();
      //       std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      //       std::cout << " " << std::endl;

      //       tc.tic();
      //       std::cout << "-Hexahedras diskpp:" << std::endl;
      //       test_hexahedra_diskpp<RealType>(rp, material_data);
      //       tc.toc();
      //       std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      //       std::cout << " " << std::endl;

      //       tc.tic();
      //       std::cout << "-Polyhedra:" << std::endl;
      //       test_polyhedra_fvca6<RealType>(rp, material_data);
      //       tc.toc();
      //       std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      //       std::cout << " " << std::endl;
   } else if (dim == 2) {

      tc.tic();
      std::cout << "-Triangles fvca5:" << std::endl;
      test_triangles_fvca5<RealType>(rp, material_data);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " " << std::endl;

      tc.tic();
      std::cout << "-Triangles netgen:" << std::endl;
      test_triangles_netgen<RealType>(rp, material_data);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " " << std::endl;

      tc.tic();
      std::cout << "-Quadrangles fvca5:" << std::endl;
      test_quads_fvca5<RealType>(rp, material_data);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " " << std::endl;

      tc.tic();
      std::cout << "-Quadrangles diskpp:" << std::endl;
      test_quads_diskpp<RealType>(rp, material_data);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " " << std::endl;

      tc.tic();
      std::cout << "-Hexagons:" << std::endl;
      test_hexagons<RealType>(rp, material_data);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " " << std::endl;

      tc.tic();
      std::cout << "-Kershaws:" << std::endl;
      test_kershaws<RealType>(rp, material_data);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " " << std::endl;
   }
}
