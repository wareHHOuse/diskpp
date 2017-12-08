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

#include <cassert>

#include "common/eigen.hpp"
#include "hho/hho_bq.hpp"
#include "hho/hho_mass_matrix.hpp"
#include "hho/hho_utils.hpp"
#include "hho/projector.hpp"
#include "mechanics/BoundaryConditions.hpp"

namespace disk {

namespace hho {

template<typename BQData>
class assembler_bq
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef typename mesh_type::face        face_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   typedef Eigen::Triplet<scalar_type> triplet_type;

   std::vector<triplet_type> m_triplets;
   size_t                    m_num_unknowns;

   const BQData& m_bqd;

 public:
   typedef Eigen::SparseMatrix<scalar_type> sparse_matrix_type;
   typedef dynamic_vector<scalar_type>      vector_type;

   sparse_matrix_type matrix;
   vector_type        rhs;

   assembler_bq() = delete;

   assembler_bq(const mesh_type& msh, const BQData& bqd) : m_bqd(bqd)
   {
      m_num_unknowns =
        howmany_dofs(m_bqd.face_basis) * (msh.faces_size() + msh.boundary_faces_size());
      matrix = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
      rhs    = vector_type::Zero(m_num_unknowns);
   }

   template<typename LocalContrib>
   void
   assemble(const mesh_type& msh, const cell_type& cl, const LocalContrib& lc)
   {
      const size_t num_face_dofs = howmany_dofs(m_bqd.face_basis);

      const auto          fcs = faces(msh, cl);
      std::vector<size_t> l2g(fcs.size() * num_face_dofs);

      for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
         const auto fc  = fcs[face_i];
         const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const auto face_id     = eid.second;
         const auto face_offset = face_id * num_face_dofs;
         const auto pos         = face_i * num_face_dofs;

         for (size_t i = 0; i < num_face_dofs; i++) {
            l2g[pos + i] = face_offset + i;
         }
      }

      assert(lc.first.rows() == fcs.size() * num_face_dofs);
      assert(lc.first.rows() == lc.first.cols());
      assert(lc.first.rows() == lc.second.size());
      assert(lc.second.size() == l2g.size());

#ifdef FILL_COLMAJOR
      for (size_t j = 0; j < lc.first.cols(); j++) {
         for (size_t i = 0; i < lc.first.rows(); i++)
            m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));

         rhs(l2g.at(j)) += lc.second(j);
      }
#else
      for (size_t i = 0; i < lc.first.rows(); i++) {
         for (size_t j = 0; j < lc.first.cols(); j++)
            m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));

         rhs(l2g.at(i)) += lc.second(i);
      }
#endif
   }

   template<typename Function>
   void
   impose_boundary_conditions(const mesh_type& msh, const Function& bc)
   {
      const size_t face_degree   = m_bqd.face_degree();
      const size_t num_face_dofs = howmany_dofs(m_bqd.face_basis);
      size_t       face_i        = 0;
      for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++) {
         auto bfc = *itor;

         const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const auto face_id              = eid.second;
         const auto face_offset          = face_id * num_face_dofs;
         const auto face_offset_lagrange = (msh.faces_size() + face_i) * num_face_dofs;

         const auto        MFF   = mass_matrix(msh, bfc, m_bqd, 0, face_degree);
         const vector_type rhs_f = compute_rhs(msh, bfc, bc, m_bqd, face_degree);

         assert(MFF.rows() == num_face_dofs && MFF.cols() == num_face_dofs);
         assert(rhs_f.rows() == num_face_dofs);

#ifdef FILL_COLMAJOR
         for (size_t j = 0; j < MFF.cols(); j++) {
            for (size_t i = 0; i < MFF.rows(); i++) {
               m_triplets.push_back(
                 triplet_type(face_offset + i, face_offset_lagrange + j, MFF(i, j)));
               m_triplets.push_back(
                 triplet_type(face_offset_lagrange + j, face_offset + i, MFF(i, j)));
            }
            rhs(face_offset_lagrange + j) = rhs_f(j);
         }
#else
         for (size_t i = 0; i < MFF.rows(); i++) {
            for (size_t j = 0; j < MFF.cols(); j++) {
               m_triplets.push_back(
                 triplet_type(face_offset + i, face_offset_lagrange + j, MFF(i, j)));
               m_triplets.push_back(
                 triplet_type(face_offset_lagrange + j, face_offset + i, MFF(i, j)));
            }
            rhs(face_offset_lagrange + i) = rhs_f(i);
         }
#endif

         face_i++;
      }
   }

   // impose boundary condition for nonlinear problems
   template<typename Function>
   void
   impose_boundary_conditions_nl(const mesh_type&                msh,
                                 const Function&                 bc,
                                 const std::vector<vector_type>& sol_faces,
                                 const std::vector<vector_type>& sol_lagr)
   {
      const size_t face_degree   = m_bqd.face_degree();
      const size_t num_face_dofs = howmany_dofs(m_bqd.face_basis);
      size_t       face_i        = 0;
      for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++) {
         auto bfc = *itor;

         const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const auto face_id              = eid.second;
         const auto face_offset          = face_id * num_face_dofs;
         const auto face_offset_lagrange = (msh.faces_size() + face_i) * num_face_dofs;

         const auto        MFF = mass_matrix(msh, bfc, m_bqd, 0, face_degree);
         const vector_type rhs_f =
           compute_rhs(msh, bfc, bc, m_bqd, face_degree) - MFF * sol_faces.at(face_id);
         const vector_type rhs_l = MFF * sol_lagr.at(face_i);

         assert(MFF.rows() == num_face_dofs && MFF.cols() == num_face_dofs);
         assert(rhs_f.rows() == num_face_dofs && rhs_l.rows() == num_face_dofs);

#ifdef FILL_COLMAJOR
         for (size_t j = 0; j < MFF.cols(); j++) {
            for (size_t i = 0; i < MFF.rows(); i++) {
               m_triplets.push_back(
                 triplet_type(face_offset + i, face_offset_lagrange + j, MFF(i, j)));
               m_triplets.push_back(
                 triplet_type(face_offset_lagrange + j, face_offset + i, MFF(i, j)));
            }
            rhs(face_offset_lagrange + j) = rhs_f(j);
            rhs(face_offset + j) -= rhs_l(j);
         }
#else
         for (size_t i = 0; i < MFF.rows(); i++) {
            for (size_t j = 0; j < MFF.cols(); j++) {
               m_triplets.push_back(
                 triplet_type(face_offset + i, face_offset_lagrange + j, MFF(i, j)));
               m_triplets.push_back(
                 triplet_type(face_offset_lagrange + j, face_offset + i, MFF(i, j)));
            }
            rhs(face_offset_lagrange + i) = rhs_f(i);
            rhs(face_offset + i) -= rhs_l(i);
         }
#endif

         face_i++;
      }
   }

   void
   finalize()
   {
      matrix.setFromTriplets(m_triplets.begin(), m_triplets.end());
      m_triplets.clear();
   }

   void
   finalize(sparse_matrix_type& mat, vector_type& vec)
   {
      mat = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
      mat.setFromTriplets(m_triplets.begin(), m_triplets.end());
      m_triplets.clear();
      vec = rhs;
   }
};

template<typename BQData>
class assembler_by_elimination_bq
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef typename mesh_type::face        face_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   typedef Eigen::Triplet<scalar_type> triplet_type;

   std::vector<triplet_type> m_triplets;
   size_t                    m_num_unknowns;

   const BQData& m_bqd;

   typedef projector_bq<BQData> projector_type;

 public:
   typedef Eigen::SparseMatrix<scalar_type> sparse_matrix_type;
   typedef dynamic_vector<scalar_type>      vector_type;

   sparse_matrix_type  matrix;
   vector_type         rhs;
   std::vector<size_t> face_compress_map, face_expand_map;

   assembler_by_elimination_bq() = delete;

   assembler_by_elimination_bq(const mesh_type& msh, const BQData& bqd) : m_bqd(bqd)
   {
      m_num_unknowns = howmany_dofs(m_bqd.face_basis) * msh.internal_faces_size();
      matrix         = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
      rhs            = vector_type::Zero(m_num_unknowns);

      face_compress_map.resize(msh.faces_size());
      face_expand_map.resize(msh.internal_faces_size());
      size_t fn = 0, fi = 0;
      for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++, fi++) {
         if (msh.is_boundary(*itor)) continue;

         face_compress_map.at(fi) = fn;
         face_expand_map.at(fn)   = fi;
         fn++;
      }
   }

   template<typename LocalContrib, typename Function>
   void
   assemble(const mesh_type& msh, const cell_type& cl, const LocalContrib& lc, const Function& bc)
   {
      const auto          face_degree   = m_bqd.face_degree();
      const auto          num_face_dofs = howmany_dofs(m_bqd.face_basis);
      const auto          fcs           = faces(msh, cl);
      std::vector<size_t> l2g(fcs.size() * num_face_dofs);
      vector_type         rhs_bc = vector_type::Zero(fcs.size() * num_face_dofs);

      projector_type projector(m_bqd);

      for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
         const auto fc             = fcs[face_i];
         const bool fc_is_boundary = msh.is_boundary(fc);
         const auto eid            = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const auto face_id     = eid.second;
         const auto face_offset = face_compress_map.at(face_id) * num_face_dofs;
         const auto pos         = face_i * num_face_dofs;

         for (size_t i = 0; i < num_face_dofs; i++) {
            if (fc_is_boundary)
               l2g.at(pos + i) = 0xDEADBEEF;
            else
               l2g.at(pos + i) = face_offset + i;
         }

         if (fc_is_boundary) {
            const vector_type sol_F = projector.projectOnFace(msh, fc, bc, face_degree);
            assert(sol_F.size() == num_face_dofs);

            for (size_t face_j = 0; face_j < fcs.size(); face_j++) {
               const auto fcj             = fcs[face_j];
               const bool fcj_is_boundary = msh.is_boundary(fcj);

               if (!fcj_is_boundary) {
                  rhs_bc.block(face_j * num_face_dofs, 0, num_face_dofs, 1) +=
                    lc.first.block(face_j * num_face_dofs, pos, num_face_dofs, num_face_dofs) *
                    sol_F;
               }
            }
         }
      }

      assert(lc.first.rows() == lc.first.cols());
      assert(lc.first.rows() == lc.second.size());
      assert(lc.second.size() == l2g.size());

#ifdef FILL_COLMAJOR
      for (size_t j = 0; j < lc.first.cols(); j++) {
         if (l2g[j] == 0xDEADBEEF) continue;

         for (size_t i = 0; i < lc.first.rows(); i++) {
            if (l2g[i] == 0xDEADBEEF) continue;

            m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));
         }
         rhs(l2g[j]) += lc.second(j) - rhs_bc(j);
      }
#else
      for (size_t i = 0; i < lc.first.rows(); i++) {
         if (l2g[i] == 0xDEADBEEF) continue;

         for (size_t j = 0; j < lc.first.cols(); j++) {
            if (l2g[j] == 0xDEADBEEF) continue;

            m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));
         }
         rhs(l2g.at(i)) += lc.second(i) - rhs_bc(i);
      }
#endif
   }

   // assembler for non-linear problems
   template<typename LocalContrib, typename Function>
   void
   assemble_nl(const mesh_type&                msh,
               const cell_type&                cl,
               const LocalContrib&             lc,
               const Function&                 bc,
               const std::vector<vector_type>& sol_F)
   {
      assert(sol_F.size() == msh.faces_size());

      const auto face_degree   = m_bqd.face_degree();
      const auto num_face_dofs = howmany_dofs(m_bqd.face_basis);
      const auto fcs           = faces(msh, cl);

      std::vector<size_t> l2g(fcs.size() * num_face_dofs);
      vector_type         rhs_bc = vector_type::Zero(fcs.size() * num_face_dofs);

      projector_type projector(m_bqd);

      for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
         const auto fc             = fcs[face_i];
         const bool fc_is_boundary = msh.is_boundary(fc);
         const auto eid            = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const auto face_id     = eid.second;
         const auto face_offset = face_compress_map.at(face_id) * num_face_dofs;
         const auto pos         = face_i * num_face_dofs;

         for (size_t i = 0; i < num_face_dofs; i++) {
            if (fc_is_boundary)
               l2g.at(pos + i) = 0xDEADBEEF;
            else
               l2g.at(pos + i) = face_offset + i;
         }

         if (fc_is_boundary) {
            const vector_type proj_bcf = projector.projectOnFace(msh, fc, bc, face_degree);
            assert(proj_bcf.size() == num_face_dofs);

            for (size_t face_j = 0; face_j < fcs.size(); face_j++) {
               const auto fcj             = fcs[face_j];
               const bool fcj_is_boundary = msh.is_boundary(fcj);

               if (!fcj_is_boundary) {
                  rhs_bc.block(face_j * num_face_dofs, 0, num_face_dofs, 1) +=
                    lc.first.block(face_j * num_face_dofs, pos, num_face_dofs, num_face_dofs) *
                    (proj_bcf - sol_F[face_id]);
               }
            }
         }
      }

      assert(lc.first.rows() == lc.first.cols());
      assert(lc.first.rows() == lc.second.size());
      assert(lc.second.size() == l2g.size());

#ifdef FILL_COLMAJOR
      for (size_t j = 0; j < lc.first.cols(); j++) {
         if (l2g[j] == 0xDEADBEEF) continue;

         for (size_t i = 0; i < lc.first.rows(); i++) {
            if (l2g[i] == 0xDEADBEEF) continue;

            m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));
         }
         rhs(l2g[j]) += lc.second(j) - rhs_bc(j);
      }
#else
      for (size_t i = 0; i < lc.first.rows(); i++) {
         if (l2g[i] == 0xDEADBEEF) continue;

         for (size_t j = 0; j < lc.first.cols(); j++) {
            if (l2g[j] == 0xDEADBEEF) continue;

            m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));
         }
         rhs(l2g.at(i)) += lc.second(i) - rhs_bc(i);
      }
#endif
   }

   template<typename Function>
   vector_type
   expand_solution(const mesh_type& msh, const vector_type& solution, const Function& bc)
   {
      const auto face_degree        = m_bqd.face_degree();
      const auto num_face_dofs      = howmany_dofs(m_bqd.face_basis);
      const auto num_internal_faces = msh.internal_faces_size();

      vector_type ret = vector_type::Zero(num_face_dofs * msh.faces_size());
      size_t      cfacenum(0);

      projector_type projector(m_bqd);

      for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++) {
         if (msh.is_boundary(*itor)) {
            const auto bfc = *itor;
            const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

            const auto face_id     = eid.second;
            const auto face_offset = face_id * num_face_dofs;
            ret.block(face_offset, 0, num_face_dofs, 1) =
              projector.projectOnFace(msh, bfc, bc, face_degree);
         } else {
            const size_t src_block_offset = cfacenum * num_face_dofs;
            const size_t dst_block_offset = face_expand_map.at(cfacenum) * num_face_dofs;
            const size_t block_size       = num_face_dofs;
            ret.block(dst_block_offset, 0, block_size, 1) =
              solution.block(src_block_offset, 0, block_size, 1);

            cfacenum++;
         }
      }

      return ret;
   }

   template<typename Function>
   vector_type
   expand_solution_nl(const mesh_type&                msh,
                      const vector_type&              solution,
                      const Function&                 bc,
                      const std::vector<vector_type>& sol_F)
   {
      assert(sol_F.size() == msh.faces_size());

      const auto face_degree        = m_bqd.face_degree();
      const auto num_face_dofs      = howmany_dofs(m_bqd.face_basis);
      const auto num_internal_faces = msh.internal_faces_size();

      vector_type ret = vector_type::Zero(num_face_dofs * msh.faces_size());
      size_t      cfacenum(0);

      projector_type projector(m_bqd);

      for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++) {
         if (msh.is_boundary(*itor)) {
            const auto bfc = *itor;
            const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

            const auto face_id     = eid.second;
            const auto face_offset = face_id * num_face_dofs;

            const auto proj_bcf = projector.projectOnFace(msh, bfc, bc, face_degree);
            ret.block(face_offset, 0, num_face_dofs, 1) = (proj_bcf - sol_F[face_id]);

         } else {
            const size_t src_block_offset = cfacenum * num_face_dofs;
            const size_t dst_block_offset = face_expand_map.at(cfacenum) * num_face_dofs;
            const size_t block_size       = num_face_dofs;
            ret.block(dst_block_offset, 0, block_size, 1) =
              solution.block(src_block_offset, 0, block_size, 1);

            cfacenum++;
         }
      }

      return ret;
   }

   void
   finalize()
   {
      matrix.setFromTriplets(m_triplets.begin(), m_triplets.end());
      m_triplets.clear();
   }

   void
   finalize(sparse_matrix_type& mat, vector_type& vec)
   {
      mat = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
      mat.setFromTriplets(m_triplets.begin(), m_triplets.end());
      m_triplets.clear();
      vec = rhs;
   }
};

// template<typename BQData>
// class assembler_by_elimination_mechanics_bq
// {
//    typedef typename BQData::mesh_type      mesh_type;
//    typedef typename mesh_type::scalar_type scalar_type;
//    typedef typename mesh_type::cell        cell_type;
//    typedef typename mesh_type::face        face_type;

//    typedef dynamic_matrix<scalar_type> matrix_type;

//    typedef Eigen::Triplet<scalar_type> triplet_type;

//    const static size_t dimension = mesh_type::dimension;

//    std::vector<triplet_type> m_triplets;
//    size_t                    m_num_unknowns;

//    const BQData& m_bqd;

//    typedef projector_bq<BQData> projector_type;

//  public:
//    typedef Eigen::SparseMatrix<scalar_type> sparse_matrix_type;
//    typedef dynamic_vector<scalar_type>      vector_type;

//    sparse_matrix_type  matrix;
//    vector_type         rhs;
//    std::vector<size_t> face_compress_map;

//    assembler_by_elimination_mechanics_bq() = delete;

//    assembler_by_elimination_mechanics_bq(const mesh_type&          msh,
//                                          const BQData&             bqd,
//                                          const BoundaryConditions& boundary_conditions) :
//      m_bqd(bqd)
//    {
//       const auto num_face_dofs = howmany_dofs(m_bqd.face_basis);

//       face_compress_map.resize(msh.faces_size());

//       size_t total_dofs = 0;
//       for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++) {
//          const auto bfc = *itor;
//          const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
//          if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

//          const auto face_id = eid.second;

//          if (boundary_conditions.is_dirichlet_face(face_id)) {
//             face_compress_map.at(face_id) = total_dofs;
//             switch (boundary_conditions.dirichlet_boundary_type(face_id)) {
//                case DIRICHLET: {
//                   break;
//                }
//                case CLAMPED: {
//                   break;
//                }

//                case default: {
//                   const auto free_dofs =
//                     num_face_dofs - boundary_conditions.dirichlet_imposed_dofs(face_id);
//                   total_dofs += free_dofs;
//                   break;
//                }
//             }
//          } else {
//             face_compress_map.at(face_id) = total_dofs;
//             total_dofs += num_face_dofs;
//          }
//       }

//       m_num_unknowns = total_dofs;
//       matrix         = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
//       rhs            = vector_type::Zero(m_num_unknowns);
//    }

//    template<typename LocalContrib, typename Function>
//    void
//    assemble(const mesh_type&          msh,
//             const cell_type&          cl,
//             const LocalContrib&       lc,
//             const Function&           bc,
//             const BoundaryConditions& boundary_conditions)
//    {
//       const auto num_face_dofs = howmany_dofs(m_bqd.face_basis);
//       const auto num_dim_dofs  = num_face_dofs / dimension;

//       const auto          fcs = faces(msh, cl);
//       std::vector<size_t> l2g(fcs.size() * num_face_dofs);
//       vector_type         rhs_bc = vector_type::Zero(fcs.size() * num_face_dofs);

//       projector_type projector(m_bqd);

//       for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
//          const auto fc = fcs[face_i];

//          auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
//          if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

//          const auto face_id                  = eid.second;
//          const bool fc_is_dirichlet_boundary = boundary_conditions.is_dirichlet_face(face_id);

//          const auto face_offset = face_compress_map.at(face_id);

//          const auto pos = face_i * num_face_dofs;

//          for (size_t i = 0; i < face_basis_size; i++) {
//             if (fc_is_dirichlet_boundary)
//                l2g.at(pos + i) = 0xDEADBEEF;
//             else
//                l2g.at(pos + i) = face_offset + i;
//          }

//          if (fc_is_dirichlet_boundary) {
//             const vector_type proj_bcf = projector.projectOnFace(msh, fc, bc);
//             for (size_t face_j = 0; face_j < fcs.size(); face_j++) {
//                const auto fcj  = fcs[face_j];
//                auto       eidj = find_element_id(msh.faces_begin(), msh.faces_end(), fcj);
//                if (!eidj.first) throw std::invalid_argument("This is a bug: face not found");

//                const auto face_idj = eidj.second;
//                const bool fcj_is_dirichlet_boundary =
//                  boundary_conditions.is_dirichlet_face(face_idj);

//                if (!fcj_is_dirichlet_boundary) {
//                   const auto rhs_bcf =
//                     lc.first.block(face_j * num_face_dofs, pos, num_face_dofs, num_face_dofs) *
//                     (proj_bcf - sol_F[face_id]);
//                   switch (boundary_conditions.dirichlet_boundary_type(face_id)) {
//                      case DIRICHLET: {
//                         rhs_bc.block(face_j * num_face_dofs, 0, num_face_dofs, 1) += bcf;
//                         break;
//                      }
//                      case CLAMPED: {
//                         rhs_bc.block(face_j * num_face_dofs, 0, num_face_dofs, 1) -=
//                           lc.first.block(
//                             face_j * face_basis_size, pos, face_basis_size, face_basis_size) *
//                           sol_F[face_id];
//                         break;
//                      }
//                      case DX: {

//                         break;
//                      }
//                      case DY: {

//                         break;
//                      }
//                      case DZ: {

//                         break;
//                      }
//                      case DXDY: {

//                         break;
//                      }
//                      case DXDZ: {

//                         break;
//                      }
//                      case DYDZ: {

//                         break;
//                      }
//                      default: {
//                         std::cout << "Unknown Dirichlet Conditions: we do nothing (assembly)"
//                                   << std::endl;
//                         break;
//                      }
//                   }
//                }
//             }
//          }
//       }

//       assert(lc.first.rows() == lc.first.cols());
//       assert(lc.first.rows() == lc.second.size());
//       assert(lc.second.size() == l2g.size());

//       // std::cout << lc.second.size() << " " << l2g.size() << std::endl;

//       for (size_t i = 0; i < lc.first.rows(); i++) {
//          if (l2g[i] == 0xDEADBEEF) continue;

//          for (size_t j = 0; j < lc.first.cols(); j++) {
//             if (l2g[j] == 0xDEADBEEF) continue;

//             m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));
//          }
//          rhs(l2g.at(i)) += lc.second(i) - rhs_bc(i);
//       }
//    }

//    template<typename LocalContrib, typename Function>
//    void
//    assemble_nl(const mesh_type&                msh,
//                const cell_type&                cl,
//                const LocalContrib&             lc,
//                const Function&                 bc,
//                const BoundaryConditions&       boundary_conditions,
//                const std::vector<vector_type>& sol_F)
//    {
//       const auto          face_basis_size = howmany_dofs(m_bqd.face_basis);
//       const auto          fcs             = faces(msh, cl);
//       std::vector<size_t> l2g(fcs.size() * face_basis_size);
//       vector_type         rhs_bc = vector_type::Zero(fcs.size() * face_basis_size);

//       projector_type projector(m_bqd);

//       for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
//          const auto fc = fcs[face_i];

//          auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
//          if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

//          const auto face_id                  = eid.second;
//          const bool fc_is_dirichlet_boundary = boundary_conditions.is_dirichlet_face(face_id);

//          const auto face_offset = face_compress_map.at(face_id) * face_basis_size;

//          const auto pos = face_i * face_basis_size;

//          for (size_t i = 0; i < face_basis_size; i++) {
//             if (fc_is_dirichlet_boundary)
//                l2g.at(pos + i) = 0xDEADBEEF;
//             else
//                l2g.at(pos + i) = face_offset + i;
//          }

//          if (fc_is_dirichlet_boundary) {
//             const vector_type proj_bcf = projector.projectOnFace(msh, fc, bc);
//             for (size_t face_j = 0; face_j < fcs.size(); face_j++) {
//                const auto fcj  = fcs[face_j];
//                auto       eidj = find_element_id(msh.faces_begin(), msh.faces_end(), fcj);
//                if (!eidj.first) throw std::invalid_argument("This is a bug: face not found");

//                const auto face_idj = eidj.second;
//                const bool fcj_is_dirichlet_boundary =
//                  boundary_conditions.is_dirichlet_face(face_idj);

//                if (!fcj_is_dirichlet_boundary) {
//                   switch (boundary_conditions.dirichlet_boundary_type(face_id)) {
//                      case DIRICHLET: {
//                         rhs_bc.block(face_j * face_basis_size, 0, face_basis_size, 1) +=
//                           lc.first.block(
//                             face_j * face_basis_size, pos, face_basis_size, face_basis_size) *
//                           (proj_bcf - sol_F[face_id]);
//                         break;
//                      }
//                      case CLAMPED: {
//                         rhs_bc.block(face_j * face_basis_size, 0, face_basis_size, 1) -=
//                           lc.first.block(
//                             face_j * face_basis_size, pos, face_basis_size, face_basis_size) *
//                           sol_F[face_id];
//                         break;
//                      }
//                      default: {
//                         std::cout << "Unknown Dirichlet Conditions: we do nothing (assembly)"
//                                   << std::endl;
//                         break;
//                      }
//                   }
//                }
//             }
//          }
//       }

//       assert(lc.first.rows() == lc.first.cols());
//       assert(lc.first.rows() == lc.second.size());
//       assert(lc.second.size() == l2g.size());

//       // std::cout << lc.second.size() << " " << l2g.size() << std::endl;

//       for (size_t i = 0; i < lc.first.rows(); i++) {
//          if (l2g[i] == 0xDEADBEEF) continue;

//          for (size_t j = 0; j < lc.first.cols(); j++) {
//             if (l2g[j] == 0xDEADBEEF) continue;

//             m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));
//          }
//          rhs(l2g.at(i)) += lc.second(i) - rhs_bc(i);
//       }
//    }

//    template<typename Function>
//    vector_type
//    expand_solution_nl(const mesh_type&                msh,
//                       const vector_type&              solution,
//                       const Function&                 bc,
//                       const BoundaryConditions&       boundary_conditions,
//                       const std::vector<vector_type>& sol_F)
//    {
//       const auto face_basis_size    = howmany_dofs(m_bqd.face_basis);
//       const auto num_internal_faces = msh.internal_faces_size();

//       vector_type ret = vector_type::Zero(face_basis_size * msh.faces_size());
//       size_t      cfacenum(0);

//       projector_type projector(m_bqd);

//       for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++) {
//          auto bfc = *itor;
//          auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
//          if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

//          const auto face_id = eid.second;

//          if (boundary_conditions.is_dirichlet_face(face_id)) {
//             const vector_type proj_bcf    = projector.projectOnFace(msh, bfc, bc);
//             const auto        face_offset = face_id * face_basis_size;
//             switch (boundary_conditions.dirichlet_boundary_type(face_id)) {
//                case DIRICHLET: {
//                   ret.block(face_offset, 0, face_basis_size, 1) = (proj_bcf - sol_F[face_id]);
//                   break;
//                }
//                case CLAMPED: {
//                   ret.block(face_offset, 0, face_basis_size, 1) = -sol_F[face_id];
//                   break;
//                }
//                default: {
//                   std::cout << "Unknown Dirichlet Conditions: we do nothing (expand)" <<
//                   std::endl; break;
//                }
//             }
//          } else {
//             const size_t src_block_offset = cfacenum * face_basis_size;
//             const size_t dst_block_offset = face_expand_map.at(cfacenum) * face_basis_size;
//             const size_t block_size       = face_basis_size;
//             ret.block(dst_block_offset, 0, block_size, 1) =
//               solution.block(src_block_offset, 0, block_size, 1);

//             cfacenum++;
//          }
//       }

//       return ret;
//    }

//    // ne marche pas
//    template<typename NeumannFunction>
//    void
//    impose_neumann_boundary_conditions(const mesh_type&          msh,
//                                       const NeumannFunction&    g,
//                                       const BoundaryConditions& boundary_conditions)
//    {
//       const size_t face_basis_size = howmany_dofs(m_bqd.face_basis);

//       for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++) {
//          auto bfc = *itor;

//          const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
//          if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

//          const auto face_id = eid.second;

//          if (boundary_conditions.is_neumann_face(face_id)) {
//             const size_t face_offset =
//               face_compress_map.at(face_id) * face_basis_size; // trouver bon decalage
//             vector_type TF = vector_type::Zero(face_basis_size);

//             const auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, bfc);
//             for (auto& qp : face_quadpoints) {
//                auto fphi = m_bqd.face_basis.eval_functions(msh, bfc, qp.point());
//                assert(fphi.size() == face_basis_size);

//                for (size_t i = 0; i < face_basis_size; i++)
//                   TF(i) += qp.weight() * mm_prod(g(qp.point()), fphi[i]);
//             }

//             for (size_t i = 0; i < TF.rows(); i++) {
//                rhs(face_offset + i) += TF(i);
//             }
//          }
//       }
//    }

//    void
//    finalize()
//    {
//       matrix.setFromTriplets(m_triplets.begin(), m_triplets.end());
//       m_triplets.clear();
//    }

//    void
//    finalize(sparse_matrix_type& mat, vector_type& vec)
//    {
//       mat = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
//       mat.setFromTriplets(m_triplets.begin(), m_triplets.end());
//       m_triplets.clear();
//       vec = rhs;
//    }
// };

template<typename BQData>
class assembler_by_elimination_without_static_condensation_bq
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef typename mesh_type::face        face_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   typedef Eigen::Triplet<scalar_type> triplet_type;

   std::vector<triplet_type> m_triplets;
   size_t                    m_num_unknowns;

   const BQData& m_bqd;

   typedef projector_bq<BQData> projector_type;

 public:
   typedef Eigen::SparseMatrix<scalar_type> sparse_matrix_type;
   typedef dynamic_vector<scalar_type>      vector_type;

   sparse_matrix_type  matrix;
   vector_type         rhs;
   std::vector<size_t> face_compress_map, face_expand_map;

   assembler_by_elimination_without_static_condensation_bq() = delete;

   assembler_by_elimination_without_static_condensation_bq(const mesh_type& msh,
                                                           const BQData&    bqd) :
     m_bqd(bqd)
   {
      const size_t cells_dofs = howmany_dofs(m_bqd.cell_basis) * msh.cells_size();
      const size_t faces_dofs = howmany_dofs(m_bqd.face_basis) * msh.internal_faces_size();
      m_num_unknowns          = cells_dofs + faces_dofs;
      matrix                  = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
      rhs                     = vector_type::Zero(m_num_unknowns);

      face_compress_map.resize(msh.faces_size());
      face_expand_map.resize(msh.internal_faces_size());
      size_t fn = 0, fi = 0;
      for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++, fi++) {
         if (msh.is_boundary(*itor)) continue;

         face_compress_map.at(fi) = fn;
         face_expand_map.at(fn)   = fi;
         fn++;
      }
   }

   template<typename LocalContrib, typename Function>
   void
   assemble(const mesh_type& msh, const cell_type& cl, const LocalContrib& lc, const Function& bc)
   {
      const auto   cell_basis_size = howmany_dofs(m_bqd.cell_basis);
      const auto   face_basis_size = howmany_dofs(m_bqd.face_basis);
      const auto   fcs             = faces(msh, cl);
      const size_t total_dof       = cell_basis_size + fcs.size() * face_basis_size;

      std::vector<size_t> l2g(total_dof);

      const auto cid = find_element_id(msh.cells_begin(), msh.cells_end(), cl);
      if (!cid.first) throw std::invalid_argument("This is a bug: cell not found");

      const auto cell_id     = cid.second;
      const auto cell_offset = cell_id * cell_basis_size;

      for (size_t i = 0; i < cell_basis_size; i++) {
         l2g[i] = cell_offset + i;
      }

      const size_t cells_offset = msh.cells_size() * cell_basis_size;

      vector_type rhs_bc = vector_type::Zero(total_dof);

      projector_type projector(m_bqd);

      for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
         const auto fc             = fcs[face_i];
         const bool fc_is_boundary = msh.is_boundary(fc);
         auto       eid            = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const auto face_id = eid.second;

         const auto face_offset = cells_offset + face_compress_map.at(face_id) * face_basis_size;

         const auto pos = cell_basis_size + face_i * face_basis_size;

         for (size_t i = 0; i < face_basis_size; i++) {
            if (fc_is_boundary)
               l2g.at(pos + i) = 0xDEADBEEF;
            else
               l2g.at(pos + i) = face_offset + i;
         }

         if (fc_is_boundary) {
            const vector_type sol_F = projector.projectOnFace(msh, fc, bc);
            for (size_t face_j = 0; face_j < fcs.size(); face_j++) {
               const auto fcj             = fcs[face_j];
               const bool fcj_is_boundary = msh.is_boundary(fcj);

               if (!fcj_is_boundary) {
                  rhs_bc.block(cell_basis_size + face_j * face_basis_size, 0, face_basis_size, 1) +=
                    lc.first.block(cell_basis_size + face_j * face_basis_size,
                                   pos,
                                   face_basis_size,
                                   face_basis_size) *
                    sol_F;
               }
            }
         }
      }

      assert(lc.first.rows() == lc.first.cols());
      assert(lc.first.rows() == lc.second.size());
      assert(lc.second.size() == l2g.size());

      // std::cout << lc.second.size() << " " << l2g.size() << std::endl;

      for (size_t i = 0; i < lc.first.rows(); i++) {
         if (l2g[i] == 0xDEADBEEF) continue;

         for (size_t j = 0; j < lc.first.cols(); j++) {
            if (l2g[j] == 0xDEADBEEF) continue;

            m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));
         }
         rhs(l2g.at(i)) += lc.second(i) - rhs_bc(i);
      }
   }

   void
   finalize()
   {
      matrix.setFromTriplets(m_triplets.begin(), m_triplets.end());
      m_triplets.clear();
   }

   void
   finalize(sparse_matrix_type& mat, vector_type& vec)
   {
      mat = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
      mat.setFromTriplets(m_triplets.begin(), m_triplets.end());
      m_triplets.clear();
      vec = rhs;
   }
};

} // hho
} // disk
