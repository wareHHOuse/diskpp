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

#include <iostream>
#include <regex>
#include <unistd.h>

#include <map>

#include "hho/hho.hpp"
#include "hho/hho_stiffness_matrix.hpp"
#include "loaders/loader.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

template<typename MeshType>
void
test_quadrature_loop(MeshType& msh, size_t degree)
{
   typedef MeshType                        mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef typename mesh_type::face        face_type;

   typedef disk::hho::
     basis_quadrature_data<mesh_type, disk::scaled_monomial_scalar_basis, disk::quadrature>
       bqdata_type;

   const bqdata_type bqd = bqdata_type(degree, degree);

   std::cout << "Cell basis size: " << bqd.cell_basis.size() << std::endl;

   timecounter tc;

   tc.tic();
   for (auto& cl : msh) {
      const auto stiff_mat = stiffness_matrix(msh, cl, bqd, degree + 1);
   }
   tc.toc();

   std::cout << tc << " seconds" << std::endl;
}

int
main(int argc, char** argv)
{
   using RealType = double;

   char* filename     = nullptr;
   int   degree       = 1;
   int   l            = 0;
   int   elems_1d     = 8;
   bool  submesh_flag = false;
   int   ch;

   while ((ch = getopt(argc, argv, "k:")) != -1) {
      switch (ch) {
         case 'k':
            degree = atoi(optarg);
            if (degree < 0) {
               std::cout << "Degree must be positive. Falling back to 1." << std::endl;
               degree = 1;
            }
            break;

         case 'h':
         case '?':
         default: std::cout << "wrong arguments" << std::endl; exit(1);
      }
   }

   argc -= optind;
   argv += optind;

   if (argc == 0) {
      std::cout << "Running 1D test simulation" << std::endl;

      typedef disk::generic_mesh<RealType, 1> mesh_type;

      mesh_type                              msh;
      disk::uniform_mesh_loader<RealType, 1> loader(0, 1, elems_1d);
      loader.populate_mesh(msh);

      auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto
      {
         // return 2.;
         return M_PI * M_PI * sin(p.x() * M_PI);
      };

      auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto
      {
         // return - p.x() * p.x();
         return sin(p.x() * M_PI);
      };

      test_quadrature_loop(msh, degree);

      return 0;
   }

   filename = argv[0];

   if (std::regex_match(filename, std::regex(".*\\.typ1$"))) {
      std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;

      typedef disk::generic_mesh<RealType, 2> mesh_type;

      mesh_type                            msh;
      disk::fvca5_mesh_loader<RealType, 2> loader;
      if (!loader.read_mesh(filename)) {
         std::cout << "Problem loading mesh." << std::endl;
         return 1;
      }
      loader.populate_mesh(msh);

      auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto
      {
         return M_PI * M_PI * sin(p.x() * M_PI);
         // return 6. * p.x();
      };

      auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto
      {
         return sin(p.x() * M_PI);
         // return - p.x() * p.x() * p.x();
      };

      test_quadrature_loop(msh, degree);
      // test_gradrec(msh, degree);
   }

   if (std::regex_match(filename, std::regex(".*\\.mesh2d$"))) {
      std::cout << "Guessed mesh format: Netgen 2D" << std::endl;

      typedef disk::simplicial_mesh<RealType, 2> mesh_type;

      mesh_type                             msh;
      disk::netgen_mesh_loader<RealType, 2> loader;
      if (!loader.read_mesh(filename)) {
         std::cout << "Problem loading mesh." << std::endl;
         return 1;
      }
      loader.populate_mesh(msh);

      auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto
      {
         return M_PI * M_PI * sin(p.x() * M_PI);
         // return 6. * p.x();
      };

      auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto
      {
         return sin(p.x() * M_PI);
         // return - p.x() * p.x() * p.x();
      };

      test_quadrature_loop(msh, degree);
   }

   if (std::regex_match(filename, std::regex(".*\\.msh$"))) {
      std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;

      typedef disk::generic_mesh<RealType, 3> mesh_type;

      mesh_type                            msh;
      disk::fvca6_mesh_loader<RealType, 3> loader;

      if (!loader.read_mesh(filename)) {
         std::cout << "Problem loading mesh." << std::endl;
         return 1;
      }
      /*
      loader.populate_mesh(msh);

      auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto {
          return M_PI * M_PI * sin(p.x() * M_PI);
          //return 1.0;
      };

      auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
          return sin(p.x() * M_PI);
          //return -p.x() * p.x() * 0.5;
      };
      */

      // test_diffusion(msh, f, sf, degree, "plot.dat");
      // test_gradrec(msh, degree);
   }

   if (std::regex_match(filename, std::regex(".*\\.mesh$"))) {
      std::cout << "Guessed mesh format: Netgen 3D" << std::endl;

      typedef disk::simplicial_mesh<RealType, 3> mesh_type;

      mesh_type                             msh;
      disk::netgen_mesh_loader<RealType, 3> loader;
      if (!loader.read_mesh(filename)) {
         std::cout << "Problem loading mesh." << std::endl;
         return 1;
      }
      loader.populate_mesh(msh);

      auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto
      {
         return M_PI * M_PI * sin(p.x() * M_PI);
         // return 1.0;
      };

      auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto
      {
         return sin(p.x() * M_PI);
         // return -p.x() * p.x() * 0.5;
      };

      test_quadrature_loop(msh, degree);
   }

   if (std::regex_match(filename, std::regex(".*\\.quad$"))) {
      std::cout << "Guessed mesh format: Cartesian 2D" << std::endl;

      typedef disk::cartesian_mesh<RealType, 2> mesh_type;

      mesh_type                                msh;
      disk::cartesian_mesh_loader<RealType, 2> loader;
      if (!loader.read_mesh(filename)) {
         std::cout << "Problem loading mesh." << std::endl;
         return 1;
      }
      loader.populate_mesh(msh);

      auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto
      {
         return M_PI * M_PI * sin(p.x() * M_PI);
         // return 1.0;
      };

      auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto
      {
         return sin(p.x() * M_PI);
         // return -p.x() * p.x() * 0.5;
      };

      test_quadrature_loop(msh, degree);
      // test_gradrec(msh, degree);
   }

   if (std::regex_match(filename, std::regex(".*\\.hex$"))) {
      std::cout << "Guessed mesh format: Cartesian 3D" << std::endl;

      typedef disk::cartesian_mesh<RealType, 3> mesh_type;

      mesh_type                                msh;
      disk::cartesian_mesh_loader<RealType, 3> loader;
      if (!loader.read_mesh(filename)) {
         std::cout << "Problem loading mesh." << std::endl;
         return 1;
      }
      loader.populate_mesh(msh);

      auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto
      {
         return M_PI * M_PI * sin(p.x() * M_PI);
         // return 1.0;
      };

      auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto
      {
         return sin(p.x() * M_PI);
         // return -p.x() * p.x() * 0.5;
      };

      test_quadrature_loop(msh, degree);
      ;
      // test_gradrec(msh, degree);
   }

   return 0;
}
