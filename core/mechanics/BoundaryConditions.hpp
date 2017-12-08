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

#pragma once

#include "common/eigen.hpp"
#include <vector>

enum DirichletType : size_t
{
   DIRICHLET = 0,
   CLAMPED   = 1,
   DX        = 2,
   DY        = 3,
   DZ        = 4,
   DXDY      = 5,
   DXDZ      = 6,
   DYDZ      = 7
};

enum NeumannType : size_t
{
   NEUMANN = 0,
   FREE    = 1,
   NDX     = 2,
   NDY     = 3,
   NDZ     = 4,
   NDXDY   = 5,
   NDXDZ   = 6,
   NDYDZ   = 7
};

struct BoundaryType
{
   size_t id;
   size_t boundary_type;
};

class BoundaryConditions
{
 private:
   std::vector<BoundaryType>            m_neumann_conditions;
   std::vector<BoundaryType>            m_dirichlet_conditions;
   std::vector<std::pair<bool, size_t>> m_faces_dirichlet;
   std::vector<std::array<size_t, 2>>   m_lagranges_info;

   dynamic_vector<bool>   m_faces_is_dirichlet, m_faces_is_neumann;
   dynamic_vector<size_t> m_type_boundary_dirichlet;

   size_t m_nb_faces_boundary;
   size_t m_nb_faces_dirichlet;
   size_t m_nb_faces_neumann;
   size_t m_nb_lag;

   template<typename TypeMesh>
   void
   find_dirichlet_faces(const TypeMesh& msh)
   {
      for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++) {
         auto bfc = *itor;

         const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const auto face_id                 = eid.second;
         m_faces_is_dirichlet(face_id)      = true;
         m_type_boundary_dirichlet(face_id) = DirichletType::DIRICHLET;
      }
   }

   template<typename TypeMesh>
   void
   find_neumann_faces(const TypeMesh& msh)
   {
      m_faces_dirichlet.reserve(m_nb_faces_boundary);
      // m_faces_dirichlet.assign(m_nb_faces_boundary, std::make_pair(true, OTHER));
      m_nb_faces_dirichlet = m_nb_faces_boundary;
      m_nb_faces_neumann   = 0;

      if (!m_neumann_conditions.empty()) {
         size_t face(0);
         for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++) {
            auto bfc = *itor;

            const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

            const auto   face_id = eid.second;
            const size_t b_id    = msh.boundary_id(face_id);

            // Find if this face is a boundary face with Neumann Condition
            for (auto& elem : m_neumann_conditions) {
               if (b_id == elem.id) {
                  m_faces_dirichlet[face] = std::make_pair(false, elem.boundary_type);
                  m_nb_faces_dirichlet--;
                  m_nb_faces_neumann++;
                  m_faces_is_dirichlet(face_id) = false;
                  m_faces_is_neumann(face_id)   = true;
                  break;
               }
            }
            face++;
         }
      }
   }

   template<typename TypeMesh>
   void
   number_of_lag_conditions(const TypeMesh& msh)
   {
      // By default all boundary faces are dirichlet faces
      const size_t DIM = msh.dimension;
      m_nb_lag         = 0;

      size_t                face(0);
      std::array<size_t, 2> info = {0, 0};
      for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++) {
         if (m_faces_dirichlet[face].first) {
            auto bfc = *itor;

            const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

            const auto face_id            = eid.second;
            bool       dirichlet_standart = true;

            if (!m_dirichlet_conditions.empty()) {
               for (auto& elem : m_dirichlet_conditions) {
                  if (msh.boundary_id(face_id) == elem.id) {

                     switch (elem.boundary_type) {
                        case CLAMPED:
                           dirichlet_standart = false;
                           info[0]            = m_nb_lag;
                           info[1]            = DIM;
                           m_lagranges_info.push_back(info);
                           m_nb_lag += DIM;
                           m_faces_dirichlet[face].second     = CLAMPED;
                           m_type_boundary_dirichlet(face_id) = CLAMPED;
                           break;
                        // case DX:
                        // dirichlet_standart = false;
                        // info[0] = m_nb_lag; info[1] = 1;
                        // m_lagranges_info.push_back(info);
                        // m_nb_lag += 1;
                        // m_faces_dirichlet[face].second = DX;
                        // break;
                        // case DY:
                        // dirichlet_standart = false;
                        // info[0] = m_nb_lag; info[1] = 1;
                        // m_lagranges_info.push_back(info);
                        // m_nb_lag += 1;
                        // m_faces_dirichlet[face].second = DY;
                        // break;
                        // case DZ:
                        // dirichlet_standart = false;
                        // if( DIM != 3){
                        //    std::cout << "Invalid condition for face:" << face_id << std::endl;
                        //    throw std::invalid_argument(" ONLY DIM = 3 for this Dirichlet
                        //    Conditions");
                        // }
                        // else{
                        //    info[0] = m_nb_lag; info[1] = 1;
                        //    m_lagranges_info.push_back(info);
                        //    m_nb_lag += 1;
                        //    m_faces_dirichlet[face].second = DZ;
                        // }
                        // break;
                        // case DXDY:
                        // dirichlet_standart = false;
                        // if( DIM != 3){
                        //    std::cout << "Invalid condition for face:" << face_id << std::endl;
                        //    throw std::invalid_argument(" ONLY DIM = 3 for this Dirichlet
                        //    Conditions");
                        // }
                        // else{
                        //    info[0] = m_nb_lag; info[1] = 2;
                        //    m_lagranges_info.push_back(info);
                        //    m_nb_lag += 2;
                        //    m_faces_dirichlet[face].second = DXDY;
                        // }
                        // break;
                        // case DXDZ:
                        // dirichlet_standart = false;
                        // if( DIM != 3){
                        //    std::cout << "Invalid condition for face:" << face_id << std::endl;
                        //    throw std::invalid_argument(" ONLY DIM = 3 for this Dirichlet
                        //    Conditions");
                        // }
                        // else{
                        //    info[0] = m_nb_lag; info[1] = 2;
                        //    m_lagranges_info.push_back(info);
                        //    m_nb_lag += 2;
                        //    m_faces_dirichlet[face].second = DXDZ;
                        // }
                        // break;
                        // case DYDZ:
                        // dirichlet_standart = false;
                        // if( DIM != 3){
                        //    std::cout << "Invalid condition for face:" << face_id << std::endl;
                        //    throw std::invalid_argument(" ONLY DIM = 3 for this Dirichlet
                        //    Conditions");
                        // }
                        // else{
                        //    info[0] = m_nb_lag; info[1] = 2;
                        //    m_lagranges_info.push_back(info);
                        //    m_nb_lag += 2;
                        //    m_faces_dirichlet[face].second = DYDZ;
                        // }
                        // break;
                        default:
                           std::cout << "Unknown Dirichlet Conditions" << std::endl;
                           throw std::invalid_argument(" Unknown Dirichlet Conditions");
                           break;
                     }
                     break;
                  }
               }
            }

            if (dirichlet_standart) {
               info[0] = m_nb_lag;
               info[1] = DIM;
               m_lagranges_info.push_back(info);
               m_nb_lag += DIM;
            }
         }
         face++;
      }
   }

 public:
   BoundaryConditions() :
     m_nb_faces_boundary(0), m_nb_faces_dirichlet(0), m_nb_faces_neumann(0), m_nb_lag(0)
   {}

   template<typename MeshType>
   BoundaryConditions(const MeshType&                  msh,
                      const std::vector<BoundaryType>& neumann_conditions,
                      const std::vector<BoundaryType>& dirichlet_conditions) :
     m_nb_faces_boundary(msh.boundary_faces_size()),
     m_neumann_conditions(neumann_conditions), m_dirichlet_conditions(dirichlet_conditions)
   {
      m_faces_is_dirichlet.resize(msh.faces_size());
      m_faces_is_dirichlet.setConstant(false);
      m_type_boundary_dirichlet.resize(msh.faces_size());
      // m_type_boundary_dirichlet.setConstant(DirichletType::OTHER);
      find_dirichlet_faces(msh);

      m_faces_is_neumann.resize(msh.faces_size());
      m_faces_is_neumann.setConstant(false);
      find_neumann_faces(msh);
      number_of_lag_conditions(msh);
   }

   size_t
   nb_faces_boundary() const
   {
      return m_nb_faces_boundary;
   }
   size_t
   nb_faces_dirichlet() const
   {
      return m_nb_faces_dirichlet;
   }
   size_t
   nb_faces_neumann() const
   {
      return m_nb_faces_neumann;
   }
   size_t
   nb_lags() const
   {
      return m_nb_lag;
   }

   std::vector<BoundaryType>
   boundary_neumann() const
   {
      return m_neumann_conditions;
   }

   bool
   is_boundary_dirichlet(const size_t face_i) const
   {
      return m_faces_dirichlet[face_i].first;
   }

   bool
   is_dirichlet_face(const size_t face_i) const
   {
      return m_faces_is_dirichlet(face_i);
   }

   bool
   is_neumann_face(const size_t face_i) const
   {
      return m_faces_is_neumann(face_i);
   }

   size_t
   boundary_type(const size_t face_i) const
   {
      return m_faces_dirichlet[face_i].second;
   }

   size_t
   dirichlet_boundary_type(const size_t face_i) const
   {
      return m_type_boundary_dirichlet(face_i);
   }

   size_t
   nb_lag_conditions_faceI(const size_t face_i) const
   {
      return m_lagranges_info[face_i][1];
   }

   size_t
   begin_lag_conditions_faceI(const size_t face_i) const
   {
      return m_lagranges_info[face_i][0];
   }

   void
   boundary_info() const
   {
      std::cout << "Number of boundary faces: " << m_nb_faces_boundary << std::endl;
      std::cout << "including: " << std::endl;
      std::cout << " - Number of Dirichlet faces: " << m_nb_faces_dirichlet << std::endl;
      std::cout << " - Number of Neumann faces: " << m_nb_faces_neumann << std::endl;
      std::cout << " - Number of Lagrangian conditions: " << m_nb_lag << std::endl;
   }
};
