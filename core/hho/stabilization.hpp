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
#include "hho/hho_bq.hpp"
#include "hho/hho_mass_matrix.hpp"
#include "hho/hho_trace_matrix.hpp"
#include "hho/hho_utils.hpp"

namespace disk {

namespace hho {

template<typename BQData>
class hho_stabilization_bq
{
 private:
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef dynamic_matrix<scalar_type> matrix_type;
   typedef dynamic_vector<scalar_type> vector_type;

   const BQData& m_bqd;

 public:
   matrix_type data;

   hho_stabilization_bq(const BQData& bqd) : m_bqd(bqd) {}

   void
   compute(const mesh_type& msh, const cell_type& cl, const matrix_type& gradrec_oper)
   {
      const auto cell_degree     = m_bqd.cell_degree();
      const auto face_degree     = m_bqd.face_degree();
      const auto cell_basis_size = m_bqd.cell_basis.computed_size();
      const auto num_face_dofs   = howmany_dofs(m_bqd.face_basis);
      const auto num_cell_dofs   = howmany_dofs(m_bqd.cell_basis);

      const matrix_type mass_mat = mass_matrix(msh, cl, m_bqd, 0, cell_degree + 1);
      assert(mass_mat.rows() == cell_basis_size && mass_mat.cols() == cell_basis_size);

      const auto zero_range = m_bqd.cell_basis.range(0, cell_degree);
      const auto one_range  = m_bqd.cell_basis.range(1, cell_degree + 1);

      const auto            fcs       = faces(msh, cl);
      const auto            num_faces = fcs.size();
      const dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

      assert(gradrec_oper.rows() == one_range.size() && gradrec_oper.cols() == dsr.total_size());

      // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

      // Step 1: compute \pi_T^k p_T^k v (third term).
      const matrix_type M1 = take(mass_mat, zero_range, zero_range);
      const matrix_type M2 = take(mass_mat, zero_range, one_range);

      assert(M2.cols() == gradrec_oper.rows());

      matrix_type proj1 = -M1.llt().solve(M2 * gradrec_oper);

      // Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
      const matrix_type I_T = matrix_type::Identity(zero_range.size(), zero_range.size());
      proj1.block(0, 0, zero_range.size(), zero_range.size()) += I_T;

      data = matrix_type::Zero(dsr.total_size(), dsr.total_size());

      // Step 3: project on faces (eqn. 21)
      for (size_t face_i = 0; face_i < num_faces; face_i++) {
         const auto current_face_range = dsr.face_range(face_i);
         const auto h                  = diameter(msh, /*fcs[face_i]*/ cl);
         const auto fc                 = fcs[face_i];

         const matrix_type face_mass_matrix = mass_matrix(msh, fc, m_bqd, 0, face_degree);
         assert(face_mass_matrix.rows() == num_face_dofs &&
                face_mass_matrix.cols() == num_face_dofs);

         const matrix_type face_trace_matrix =
           trace_matrix(msh, cl, fc, m_bqd, cell_degree + 1, face_degree);
         assert(face_trace_matrix.rows() == num_face_dofs &&
                face_trace_matrix.cols() == cell_basis_size);

         Eigen::LLT<matrix_type> piKF;
         piKF.compute(face_mass_matrix);

         // Step 3a: \pi_F^k( v_F - p_T^k v )
         const auto        face_range = current_face_range.remove_offset();
         const matrix_type MR1        = take(face_trace_matrix, face_range, one_range);
         assert(MR1.rows() == num_face_dofs && MR1.cols() == gradrec_oper.rows());

         matrix_type       proj2        = piKF.solve(MR1 * gradrec_oper);
         const matrix_type I_F          = matrix_type::Identity(num_face_dofs, num_face_dofs);
         const auto        block_offset = current_face_range.min();
         proj2.block(0, block_offset, num_face_dofs, num_face_dofs) -= I_F;

         // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
         const matrix_type MR2 = take(face_trace_matrix, face_range, zero_range);
         assert(MR2.rows() == num_face_dofs && MR2.cols() == proj1.rows());

         const matrix_type proj3 = piKF.solve(MR2 * proj1);

         assert(proj2.rows() == proj3.rows() && proj2.cols() == proj3.cols());
         const matrix_type BRF = proj2 + proj3;

         data += BRF.transpose() * face_mass_matrix * BRF / h;
      }
   }
};

// 1/h_F (v_F - v_T, v_F-v_T)_F
template<typename BQData>
class hdg_stabilization_bq
{
 private:
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef dynamic_matrix<scalar_type> matrix_type;
   typedef dynamic_vector<scalar_type> vector_type;

   const BQData& m_bqd;

 public:
   matrix_type data;

   hdg_stabilization_bq(const BQData& bqd) : m_bqd(bqd) {}

   void
   compute(const mesh_type& msh, const cell_type& cl)
   {
      const auto cell_degree = m_bqd.cell_degree();
      const auto face_degree = m_bqd.face_degree();

      const auto num_face_dofs = howmany_dofs(m_bqd.face_basis);
      const auto num_cell_dofs = howmany_dofs(m_bqd.cell_basis);

      const auto fcs       = faces(msh, cl);
      const auto num_faces = fcs.size();

      const dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

      // Build v_F - v_T

      data = matrix_type::Zero(dsr.total_size(), dsr.total_size());

      // Step 3: project on faces (eqn. 21)
      for (size_t face_i = 0; face_i < num_faces; face_i++) {
         const auto current_face_range = dsr.face_range(face_i);
         const auto h                  = diameter(msh, /*fcs[face_i]*/ cl);
         const auto fc                 = fcs[face_i];

         const matrix_type face_mass_matrix = mass_matrix(msh, fc, m_bqd, 0, face_degree);
         assert(face_mass_matrix.rows() == num_face_dofs &&
                face_mass_matrix.cols() == num_face_dofs);

         const matrix_type face_trace_matrix =
           trace_matrix(msh, cl, fc, m_bqd, cell_degree, face_degree);
         assert(face_trace_matrix.rows() == num_face_dofs &&
                face_trace_matrix.cols() == num_cell_dofs);

         // Step 3a:  v_F - vT
         matrix_type       BRF          = matrix_type::Zero(num_face_dofs, dsr.total_size());
         const matrix_type I_F          = matrix_type::Identity(num_face_dofs, num_face_dofs);
         const auto        block_offset = current_face_range.min();
         BRF.block(0, block_offset, num_face_dofs, num_face_dofs) = I_F;
         BRF.block(0, 0, num_face_dofs, num_cell_dofs) -= face_trace_matrix;

         data += BRF.transpose() * face_mass_matrix * BRF / h;
      }
   }
};

// 1/h_F (v_F - pi^k_F(v_T), v_F-pi^k_F(v_T))_F
template<typename BQData>
class pikF_stabilization_bq
{
 private:
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef dynamic_matrix<scalar_type> matrix_type;
   typedef dynamic_vector<scalar_type> vector_type;

   const BQData& m_bqd;

 public:
   matrix_type data;

   pikF_stabilization_bq(const BQData& bqd) : m_bqd(bqd) {}

   void
   compute(const mesh_type& msh, const cell_type& cl)
   {
      const auto cell_degree = m_bqd.cell_degree();
      const auto face_degree = m_bqd.face_degree();

      const auto num_face_dofs = howmany_dofs(m_bqd.face_basis);
      const auto num_cell_dofs = howmany_dofs(m_bqd.cell_basis);

      const auto fcs       = faces(msh, cl);
      const auto num_faces = fcs.size();

      const auto            zero_range = m_bqd.cell_basis.range(0, cell_degree);
      const dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

      // Build \pi_F^k (v_F - v_T)

      // Step 2: v_T
      const matrix_type I_T   = matrix_type::Identity(num_cell_dofs, num_cell_dofs);
      matrix_type       proj1 = matrix_type::Zero(num_cell_dofs, dsr.total_size());
      proj1.block(0, 0, num_cell_dofs, num_cell_dofs) += I_T;

      data = matrix_type::Zero(dsr.total_size(), dsr.total_size());

      // Step 3: project on faces (eqn. 21)
      for (size_t face_i = 0; face_i < num_faces; face_i++) {
         const auto current_face_range = dsr.face_range(face_i);
         const auto h                  = diameter(msh, /*fcs[face_i]*/ cl);
         const auto fc                 = fcs[face_i];

         const matrix_type face_mass_matrix = mass_matrix(msh, fc, m_bqd, 0, face_degree);
         assert(face_mass_matrix.rows() == num_face_dofs &&
                face_mass_matrix.cols() == num_face_dofs);

         const matrix_type face_trace_matrix =
           trace_matrix(msh, cl, fc, m_bqd, cell_degree, face_degree);
         assert(face_trace_matrix.rows() == num_face_dofs &&
                face_trace_matrix.cols() == num_cell_dofs);

         Eigen::LLT<matrix_type> piKF;
         piKF.compute(face_mass_matrix);

         // Step 3a: -v_F
         const auto        face_range   = current_face_range.remove_offset();
         matrix_type       proj2        = matrix_type::Zero(num_face_dofs, dsr.total_size());
         const matrix_type I_F          = matrix_type::Identity(num_face_dofs, num_face_dofs);
         const auto        block_offset = current_face_range.min();
         proj2.block(0, block_offset, num_face_dofs, num_face_dofs) -= I_F;

         // Step 3b: \pi_F^k( v_T )
         const matrix_type proj3 = piKF.solve(face_trace_matrix * proj1);

         assert(proj2.rows() == proj3.rows() && proj2.cols() == proj3.cols());

         const matrix_type BRF = proj2 + proj3;

         data += BRF.transpose() * face_mass_matrix * BRF / h;
      }
   }
};

} // hho
} // disk
