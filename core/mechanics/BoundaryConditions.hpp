/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet (C) 2018                      nicolas.pignet@enpc.fr
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

#pragma once

#include <tuple>
#include <vector>

#include "bases/bases_utils.hpp"
#include "mesh/point.hpp"
#include "revolution/bases"

namespace disk {
namespace mechanics {

enum DirichletType : size_t
{
   DIRICHLET = 0,
   CLAMPED   = 1,
   DX        = 2,
   DY        = 3,
   DZ        = 4,
   DXDY      = 5,
   DXDZ      = 6,
   DYDZ      = 7,
   NOTHING   = 8
};

enum NeumannType : size_t
{
   NEUMANN = 9,
   FREE    = 10,
};

template<typename MeshType>
class BoundaryConditions
{
 public:
   typedef MeshType                                         mesh_type;
   typedef typename mesh_type::scalar_type                  scalar_type;
   typedef static_vector<scalar_type, mesh_type::dimension> function_type;
   typedef point<scalar_type, mesh_type::dimension>         point_type;

 private:
   const mesh_type& m_msh;

   std::vector<std::function<function_type(point_type)>> m_dirichlet_func;
   std::vector<std::function<function_type(point_type)>> m_neumann_func;

   std::vector<std::tuple<bool, size_t, size_t, size_t>> m_faces_is_dirichlet, m_faces_is_neumann;

   size_t m_dirichlet_faces, m_neumann_faces;

   scalar_type m_factor;

 public:
   BoundaryConditions() = delete;

   BoundaryConditions(const mesh_type& msh) :
     m_msh(msh), m_dirichlet_faces(0), m_neumann_faces(0), m_factor(1)
   {
      m_faces_is_dirichlet.assign(m_msh.faces_size(), std::make_tuple(false, NOTHING, 0, 0));
      m_faces_is_neumann.assign(m_msh.faces_size(), std::make_tuple(false, FREE, 0, 0));
   }

   template<typename Function>
   void
   addDirichletEverywhere(const Function& bcf)
   {
      const size_t bcf_id = m_dirichlet_func.size();
      m_dirichlet_func.push_back(bcf);

      for (auto itor = m_msh.boundary_faces_begin(); itor != m_msh.boundary_faces_end(); itor++) {
         const auto bfc = *itor;

         const auto eid = find_element_id(m_msh.faces_begin(), m_msh.faces_end(), bfc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const auto face_id = eid.second;

         m_faces_is_dirichlet.at(face_id) = std::make_tuple(true, DIRICHLET, 0, bcf_id);
         m_dirichlet_faces++;
      }
   }

   template<typename Function>
   void
   addDirichletBC(const size_t& btype, const size_t& b_id, const Function& bcf)
   {
      const size_t bcf_id = m_dirichlet_func.size();
      m_dirichlet_func.push_back(bcf);

      for (auto itor = m_msh.boundary_faces_begin(); itor != m_msh.boundary_faces_end(); itor++) {
         const auto bfc = *itor;

         const auto eid = find_element_id(m_msh.faces_begin(), m_msh.faces_end(), bfc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const auto face_id = eid.second;

         if (m_msh.boundary_id(face_id) == b_id) {
            m_faces_is_dirichlet.at(face_id) = std::make_tuple(true, btype, b_id, bcf_id);
            m_dirichlet_faces++;
         }
      }
   }

   template<typename Function>
   void
   addNeumannBC(const size_t& btype, const size_t& b_id, const Function& bcf)
   {
      const size_t bcf_id = m_neumann_func.size();
      m_neumann_func.push_back(bcf);

      for (auto itor = m_msh.boundary_faces_begin(); itor != m_msh.boundary_faces_end(); itor++) {
         const auto bfc = *itor;

         const auto eid = find_element_id(m_msh.faces_begin(), m_msh.faces_end(), bfc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const auto face_id = eid.second;

         if (m_msh.boundary_id(face_id) == b_id) {
            m_faces_is_neumann.at(face_id) = std::make_tuple(true, btype, b_id, bcf_id);
            m_neumann_faces++;
         }
      }
   }

   void
   multiplyAllFunctionsOfAFactor(const scalar_type& factor)
   {
      m_factor = factor;
   }

   size_t
   nb_faces_boundary() const
   {
      return m_dirichlet_faces + m_neumann_faces;
   }
   size_t
   nb_faces_dirichlet() const
   {
      return m_dirichlet_faces;
   }
   size_t
   nb_faces_neumann() const
   {
      return m_neumann_faces;
   }

   bool
   is_dirichlet_face(const size_t face_i) const
   {
      return std::get<0>(m_faces_is_dirichlet.at(face_i));
   }

   bool
   is_neumann_face(const size_t face_i) const
   {
      return std::get<0>(m_faces_is_neumann.at(face_i));
   }

   size_t
   dirichlet_boundary_type(const size_t face_i) const
   {
      return std::get<1>(m_faces_is_dirichlet.at(face_i));
   }

   size_t
   neumann_boundary_type(const size_t face_i) const
   {
      return std::get<1>(m_faces_is_neumann.at(face_i));
   }

   size_t
   dirichlet_boundary_id(const size_t face_i) const
   {
      return std::get<2>(m_faces_is_dirichlet.at(face_i));
   }

   size_t
   neumann_boundary_id(const size_t face_i) const
   {
      return std::get<2>(m_faces_is_neumann.at(face_i));
   }

   auto
   dirichlet_boundary_func(const size_t face_i) const
   {
      if (!is_dirichlet_face(face_i)) {
         throw std::logic_error(
           "You want the Dirichlet function of face which is not a Dirichlet face");
      }

      const auto        func   = m_dirichlet_func.at(std::get<3>(m_faces_is_dirichlet.at(face_i)));
      const scalar_type factor = m_factor;

      auto rfunc = [ func, factor ](const point_type& p) -> auto
      {
         return disk::mm_prod(factor, func(p));
      };

      return rfunc;
   }

   auto
   neumann_boundary_func(const size_t face_i) const
   {
      if (!is_neumann_face(face_i)) {
         throw std::logic_error(
           "You want the Neumann function of face which is not a Neumann face");
      }

      const auto        func   = m_neumann_func.at(std::get<3>(m_faces_is_neumann.at(face_i)));
      const scalar_type factor = m_factor;

      auto rfunc = [ func, factor ](const point_type& p) -> auto
      {
         return disk::mm_prod(factor, func(p));
      };

      return rfunc;
   }

   void
   boundary_info() const
   {
      std::cout << "Number of boundary faces: " << nb_faces_boundary() << std::endl;
      std::cout << "including: " << std::endl;
      std::cout << " - Number of Dirichlet faces: " << nb_faces_dirichlet() << std::endl;
      std::cout << " - Number of Neumann faces: " << nb_faces_neumann() << std::endl;
   }

   size_t
   dirichlet_imposed_dofs(const size_t& face_id, const size_t face_degree) const
   {
      const size_t dimension = mesh_type::dimension;
      const size_t num_face_dofs =
        revolution::vector_basis_size(face_degree, dimension - 1, dimension);
      const size_t num_dim_dofs = num_face_dofs / dimension;

      if (is_dirichlet_face(face_id)) {
         const size_t btype = dirichlet_boundary_type(face_id);

         switch (btype) {
            case DIRICHLET: {
               return num_face_dofs;
               break;
            }
            case CLAMPED: {
               return num_face_dofs;
               break;
            }
            case DX: {
               return num_dim_dofs;
               break;
            }
            case DY: {
               return num_dim_dofs;
               break;
            }
            case DZ: {
               if (dimension != 3) throw std::invalid_argument("You are not in 3D");
               return num_dim_dofs;
               break;
            }
            case DXDY: {
               return 2 * num_dim_dofs;
               break;
            }
            case DXDZ: {
               if (dimension != 3) throw std::invalid_argument("You are not in 3D");
               return 2 * num_dim_dofs;
               break;
            }
            case DYDZ: {
               if (dimension != 3) throw std::invalid_argument("You are not in 3D");
               return 2 * num_dim_dofs;
               break;
            }
            case NOTHING: {
               return 0;
               break;
            }
            default: {
               throw std::logic_error("Unknown Boundary condition");
               break;
            }
         }
      }

      return 0;
   }
};

} // end mechanics
} // end disk
