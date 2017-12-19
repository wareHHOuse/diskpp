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

#include "geometry/geometry.hpp"
#include "hho/gradient_reconstruction.hpp"
#include "hho/hho_bq.hpp"
#include "hho/hho_utils.hpp"
#include "hho/stabilization.hpp"
#include "hho/static_condensation.hpp"
#include "loaders/loader.hpp"
#include <common/condition_number.hpp>

#include "contrib/sol2/sol.hpp"
#include "contrib/timecounter.h"

template<typename Mesh>
bool
estimate_element_cond(sol::state& lua, const Mesh& msh)
{
   typedef Mesh                            mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef typename mesh_type::face        face_type;

   typedef disk::quadrature<mesh_type, cell_type> cell_quadrature_type;
   typedef disk::quadrature<mesh_type, face_type> face_quadrature_type;

   typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type> cell_basis_type;
   typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type> face_basis_type;

   typedef disk::hho::
     basis_quadrature_data<mesh_type, disk::scaled_monomial_scalar_basis, disk::quadrature>
       bqdata_type;

   typedef disk::hho::gradient_reconstruction_bq<bqdata_type> gradrec_type;
   typedef disk::hho::hho_stabilization_bq<bqdata_type>       stab_type;
   typedef disk::hho::static_condensation_bq<bqdata_type>     statcond_type;

   size_t degree = lua["config"]["degree"].get_or(1);

   auto bqd = bqdata_type(degree, degree);

   auto gradrec  = gradrec_type(bqd);
   auto stab     = stab_type(bqd);
   auto statcond = statcond_type(bqd);

   auto lua_rhs_fun = lua["right_hand_side"];
   if (!lua_rhs_fun.valid()) {
      std::cout << "[config] right_hand_side() not defined" << std::endl;
      return false;
   }

   auto f = [&](const typename mesh_type::point_type& pt) -> scalar_type {
      return lua_rhs_fun(pt.x(), pt.y());
   };

   for (auto& cl : msh) {
      gradrec.compute(msh, cl);
      stab.compute(msh, cl, gradrec.oper);
      const auto                        cell_rhs = disk::hho::compute_rhs(msh, cl, f, bqd, degree);
      const dynamic_matrix<scalar_type> loc      = gradrec.data + stab.data;
      const auto                        scnp     = statcond.compute(msh, cl, loc, cell_rhs);

      auto cbs = bqd.cell_basis.size();

      const dynamic_matrix<scalar_type> cloc        = loc.block(0, 0, cbs, cbs);
      const scalar_type                 largest_ev  = disk::largest_eigenvalue(cloc);
      const scalar_type                 smallest_ev = disk::smallest_eigenvalue(cloc);
      const scalar_type                 cond        = largest_ev / smallest_ev;
      std::cout << "Condition number: " << cond << std::endl;
      std::cout << "Sigma max: " << largest_ev << std::endl;
      std::cout << "Sigma min: " << smallest_ev << std::endl;
   }

   return true;
}

int
main(int argc, char** argv)
{
   using scalar_type = double;

   if (argc != 2) {
      std::cout << "Please specify configuration file" << std::endl;
      return 1;
   }

   sol::state lua;
   lua.open_libraries(sol::lib::math);
   lua["config"] = lua.create_table();

   auto r = lua.do_file(argv[1]);

   if (!r.valid()) {
      std::cout << "Problems opening configuration file" << std::endl;
      return 1;
   }

   std::string input_mesh = lua["config"]["input_mesh"];

   if (std::regex_match(input_mesh, std::regex(".*\\.typ1$"))) {
      std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;

      typedef disk::generic_mesh<scalar_type, 2> mesh_type;

      mesh_type                               msh;
      disk::fvca5_mesh_loader<scalar_type, 2> loader;

      if (!loader.read_mesh(input_mesh)) {
         std::cout << "Problem loading mesh." << std::endl;
         return 1;
      }

      loader.populate_mesh(msh);

      estimate_element_cond(lua, msh);
   }
}
