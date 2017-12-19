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

namespace disk {

namespace hho {

template<typename BQData>
class static_condensation_bq
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef dynamic_matrix<scalar_type> matrix_type;
   typedef dynamic_vector<scalar_type> vector_type;

   const BQData& m_bqd;

 public:
   matrix_type KTT, KTF;
   vector_type RT;
   static_condensation_bq(const BQData& bqd)
     : m_bqd(bqd)
   {}

   std::pair<matrix_type, vector_type> compute(const mesh_type&   msh,
                                               const cell_type&   cl,
                                               const matrix_type& local_mat,
                                               const vector_type& cell_rhs)
   {
      const size_t num_cell_dofs = howmany_dofs(m_bqd.cell_basis);
      const size_t num_face_dofs = howmany_dofs(m_bqd.face_basis);

      const auto   fcs       = faces(msh, cl);
      const size_t num_faces = fcs.size();

      dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

      const size_t cell_size = dsr.cell_range().size();
      const size_t face_size = dsr.all_faces_range().size();

      assert(cell_size == cell_rhs.rows() && "wrong rhs dimension");
      assert((cell_size + face_size) == local_mat.rows() && "wrong lhs rows dimension");
      assert((cell_size + face_size) == local_mat.cols() && "wrong lhs cols dimension");

      const matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
      const matrix_type K_TF = local_mat.topRightCorner(cell_size, face_size);
      const matrix_type K_FT = local_mat.bottomLeftCorner(face_size, cell_size);
      const matrix_type K_FF = local_mat.bottomRightCorner(face_size, face_size);

      assert(K_TT.cols() == cell_size && "wrong K_TT dimension");
      assert(K_TT.cols() + K_TF.cols() == local_mat.cols());
      assert(K_TT.rows() + K_FT.rows() == local_mat.rows());
      assert(K_TF.rows() + K_FF.rows() == local_mat.rows());
      assert(K_FT.cols() + K_FF.cols() == local_mat.cols());

      const auto        K_TT_ldlt = K_TT.llt();
      const matrix_type AL        = K_TT_ldlt.solve(K_TF);
      const vector_type bL        = K_TT_ldlt.solve(cell_rhs);

      const matrix_type AC = K_FF - K_FT * AL;
      const vector_type bC = /* no projection on faces, eqn. 26*/ -K_FT * bL;

      return std::make_pair(AC, bC);
   }

   std::pair<matrix_type, vector_type> compute(const mesh_type&   msh,
                                               const cell_type&   cl,
                                               const matrix_type& local_mat,
                                               const vector_type& rhs,
                                               bool               l_face)
   {
      const size_t num_cell_dofs = howmany_dofs(m_bqd.cell_basis);
      const size_t num_face_dofs = howmany_dofs(m_bqd.face_basis);

      const auto   fcs       = faces(msh, cl);
      const size_t num_faces = fcs.size();

      dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

      const size_t cell_size  = dsr.cell_range().size();
      const size_t faces_size = dsr.all_faces_range().size();

      assert((cell_size + faces_size) == rhs.rows() && "wrong rhs dimension");
      assert((cell_size + faces_size) == local_mat.rows() && "wrong lhs rows dimension");
      assert((cell_size + faces_size) == local_mat.cols() && "wrong lhs cols dimension");

      const matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
      const matrix_type K_TF = local_mat.topRightCorner(cell_size, faces_size);
      const matrix_type K_FT = local_mat.bottomLeftCorner(faces_size, cell_size);
      const matrix_type K_FF = local_mat.bottomRightCorner(faces_size, faces_size);

      assert(K_TT.cols() == cell_size && "wrong K_TT dimension");
      assert(K_TT.cols() + K_TF.cols() == local_mat.cols());
      assert(K_TT.rows() + K_FT.rows() == local_mat.rows());
      assert(K_TF.rows() + K_FF.rows() == local_mat.rows());
      assert(K_FT.cols() + K_FF.cols() == local_mat.cols());

      const vector_type rhs_cell  = rhs.block(0, 0, cell_size, 1);
      const vector_type rhs_faces = rhs.block(cell_size, 0, faces_size, 1);

      assert((rhs_cell.rows() + rhs_faces.rows()) == rhs.rows() && "wrong rhs decomposition");

      const auto        K_TT_ldlt = K_TT.llt();
      const matrix_type AL        = K_TT_ldlt.solve(K_TF);
      const vector_type bL        = K_TT_ldlt.solve(rhs_cell);

      const matrix_type AC = K_FF - K_FT * AL;
      const vector_type bC = rhs_faces - K_FT * bL;

      KTT = K_TT;
      KTF = K_TF;
      RT  = rhs_cell;

      return std::make_pair(AC, bC);
   }

   vector_type recover(const mesh_type&   msh,
                       const cell_type&   cl,
                       const matrix_type& local_mat,
                       const vector_type& cell_rhs,
                       const vector_type& solF)
   {
      const size_t num_cell_dofs = howmany_dofs(m_bqd.cell_basis);
      const size_t num_face_dofs = howmany_dofs(m_bqd.face_basis);

      const auto   fcs       = faces(msh, cl);
      const size_t num_faces = fcs.size();

      dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

      const size_t cell_size      = dsr.cell_range().size();
      const size_t all_faces_size = dsr.all_faces_range().size();

      assert(cell_rhs.rows() == cell_size && "wrong rhs_cell dimension");
      assert(solF.rows() == all_faces_size && "wrong solF dimension");
      assert((cell_size + all_faces_size) == local_mat.rows() && "wrong lhs rows dimension");
      assert((cell_size + all_faces_size) == local_mat.cols() && "wrong lhs cols dimension");

      vector_type ret(dsr.total_size());

      const matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
      const matrix_type K_TF = local_mat.topRightCorner(cell_size, all_faces_size);

      assert(K_TT.cols() == cell_size && "wrong K_TT dimension");
      assert(K_TT.cols() + K_TF.cols() == local_mat.cols());

      assert(num_cell_dofs == cell_rhs.size());
      assert((num_face_dofs * num_faces) == solF.size());

      const vector_type solT = K_TT.llt().solve(cell_rhs - K_TF * solF);

      ret.head(cell_size)      = solT;
      ret.tail(all_faces_size) = solF;

      return ret;
   }
};

} // end hho
} // end disk