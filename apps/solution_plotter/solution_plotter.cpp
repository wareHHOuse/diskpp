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

#include "common/eigen.hpp"
#include "hho/hho.hpp"
#include "loaders/loader.hpp"

int
main(int argc, char** argv)
{
   using RealType = double;

   // const char *msh_fn = "/Users/matteo/Desktop/MsHHO meshes/mesh6.mesh2d";
   const char* msh_fn = "../multiscale/newsol/L11-2.mesh2d";
   const char* sol_fn = "../multiscale/newsol/solution-L11-2.bin";

   typedef disk::simplicial_mesh<RealType, 2> mesh_type;
   typedef typename mesh_type::cell           cell_type;
   typedef typename mesh_type::scalar_type    scalar_type;

   mesh_type msh;

   disk::netgen_mesh_loader<RealType, 2> loader;
   loader.verbose(true);
   if (!loader.read_mesh(msh_fn)) {
      std::cout << "Problem loading mesh." << std::endl;
      return 1;
   }
   loader.populate_mesh(msh);

   std::ifstream ifs(sol_fn, std::ifstream::binary);

   size_t sol_num_elements, cell_basis_deg, face_basis_deg;

   ifs.read(reinterpret_cast<char*>(&sol_num_elements), sizeof(size_t));
   ifs.read(reinterpret_cast<char*>(&cell_basis_deg), sizeof(size_t));
   ifs.read(reinterpret_cast<char*>(&face_basis_deg), sizeof(size_t));

   if (sol_num_elements != msh.cells_size()) {
      std::cout << "Solution has a different number of elements than the mesh (";
      std::cout << sol_num_elements << " vs. " << msh.cells_size() << ")" << std::endl;
      return 1;
   }

   // std::cout << "Solution has " << sol_num_elements << " elements" << std::endl;
   // std::cout << "Cell basis degree: " << cell_basis_deg << std::endl;
   // std::cout << "Face basis degree: " << face_basis_deg << std::endl;

   typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type> cell_basis_type;
   typedef disk::quadrature<mesh_type, cell_type>                   cell_quadrature_type;

   cell_basis_type      cell_basis_k1(cell_basis_deg + 1);
   cell_quadrature_type cell_quadrature(2 * cell_basis_deg + 2);

   size_t                      local_dofs_size = cell_basis_k1.size();
   size_t                      dofs_vec_size   = msh.cells_size() * local_dofs_size;
   dynamic_vector<scalar_type> cell_dofs       = dynamic_vector<scalar_type>::Zero(dofs_vec_size);

   for (size_t i = 0; i < dofs_vec_size; i++)
      ifs.read(reinterpret_cast<char*>(&cell_dofs(i)), sizeof(scalar_type));

   ifs.close();

   std::ofstream ofs("pippo.dat");

   size_t elemnum = 0;
   for (auto& cl : msh) {
      auto                        offset_start = elemnum * local_dofs_size;
      dynamic_vector<scalar_type> local_dofs;
      local_dofs = cell_dofs.block(offset_start, 0, local_dofs_size, 1);

      auto bar = barycenter(msh, cl);

      auto phi = cell_basis_k1.eval_functions(msh, cl, bar);
      auto val = local_dofs.dot(phi);

      elemnum++;

      ofs << bar.x() << " " << bar.y() << " " << val << std::endl;
   }

   ofs.close();

   /*
   double norm = 0.0;
   double norm_chk = 0.0;
   size_t elemnum = 0;
   for (auto& cl : msh)
   {
       auto offset_start = elemnum * local_dofs_size;
       dynamic_vector<scalar_type> local_dofs;
       local_dofs = cell_dofs.block(offset_start, 0, local_dofs_size, 1);


       auto qps = cell_quadrature.integrate(msh, cl);
       for (auto& qp : qps)
       {
           auto pt = qp.point();

           dynamic_vector<scalar_type> grad_u = dynamic_vector<scalar_type>::Zero(2);
           auto dphi = cell_basis_k1.eval_gradients(msh, cl, pt);
           for (size_t i = 0; i < local_dofs.size(); i++)
               grad_u += (dphi.block(i, 0, 1, 2) * local_dofs(i)).transpose();

           static_matrix<scalar_type,2,2> ret = static_matrix<scalar_type,2,2>::Identity();
           auto c = std::cos(M_PI * pt.x()/0.02);
           auto s = std::sin(M_PI * pt.y()/0.02);
           ret = ret * sqrt((1 + 100*c*c*s*s));

           auto vec = ret * grad_u;

           norm += qp.weight() * vec.dot(vec);



           auto phi = cell_basis_k1.eval_functions(msh, cl, pt);
           norm_chk += qp.weight() * sin(pt.x()) * sin(pt.y()) * local_dofs.dot(phi);
       }

       elemnum++;
   }

   std::cout << norm << std::endl;
   std::cout << norm_chk << std::endl;
   */

   return 0;
}
