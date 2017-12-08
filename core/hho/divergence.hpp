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

#include "bases/bases_ranges.hpp"
#include "bases/bases_utils.hpp"
#include "common/eigen.hpp"
#include "hho/hho_bq.hpp"
#include "hho/hho_utils.hpp"
#include "quadratures/quadratures.hpp"

namespace disk {
namespace hho {

template<typename BQData>
class divergence_reconstruction_bq
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef dynamic_matrix<scalar_type> matrix_type;
   typedef dynamic_vector<scalar_type> vector_type;

   const BQData& m_bqd;

 public:
   matrix_type oper;
   matrix_type data;

   divergence_reconstruction_bq() = delete;

   divergence_reconstruction_bq(const BQData& bqd) : m_bqd(bqd) {}

   void
   compute(const mesh_type& msh, const cell_type& cl)
   {
      const auto cell_degree = m_bqd.cell_degree();
      const auto face_degree = m_bqd.face_degree();

      const auto num_face_dofs  = howmany_dofs(m_bqd.face_basis);
      const auto num_cell_dofs  = howmany_dofs(m_bqd.cell_basis);
      const auto div_basis_size = howmany_dofs(m_bqd.div_cell_basis);

      // LHS
      matrix_type MD = matrix_type::Zero(div_basis_size, div_basis_size);

      const auto cell_quadpoints = m_bqd.div_cell_quadrature.integrate(msh, cl);
      assert(2 * cell_degree <= m_bqd.div_cell_quadrature.order());

      for (auto& qp : cell_quadpoints) {
         const auto phi_d =
           m_bqd.div_cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);
         assert(phi_d.rows() == div_basis_size);

         MD += qp.weight() * phi_d * phi_d.transpose();
      }

      /* RHS, volumetric part. */
      const auto fcs       = faces(msh, cl);
      const auto num_faces = fcs.size();

      dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

      matrix_type BD = matrix_type::Zero(div_basis_size, dsr.total_size());

      // -(v_F, div(q))_T
      for (auto& qp : cell_quadpoints) {
         const auto phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);
         const auto dphi_d =
           m_bqd.div_cell_basis.eval_gradients(msh, cl, qp.point(), 0, cell_degree);

         assert(dphi_d.rows() == div_basis_size);
         assert(phi.size() == num_cell_dofs);

         const auto dphi_d_i = convert_to_vector(dphi_d);
         for (size_t j = 0; j < num_cell_dofs; j++) {
            for (size_t i = 0; i < div_basis_size; i++) {
               BD(i, j) -= qp.weight() * mm_prod(dphi_d_i[i], phi[j]);
            }
         }
      }

      size_t face_offset = num_cell_dofs;

      //  (v_F.nTF, q)_F
      for (size_t face_i = 0; face_i < num_faces; face_i++) {
         const auto fc              = fcs[face_i];
         const auto n               = normal(msh, cl, fc);
         const auto face_quadpoints = m_bqd.div_face_quadrature.integrate(msh, fc);
         assert(2 * face_degree <= m_bqd.div_face_quadrature.order());

         for (auto& qp : face_quadpoints) {
            const auto fphi_d =
              m_bqd.div_cell_basis.eval_functions(msh, cl, qp.point(), 0, face_degree);
            const auto fphi = m_bqd.face_basis.eval_functions(msh, fc, qp.point(), 0, face_degree);

            assert(fphi_d.rows() == div_basis_size);
            assert(fphi.size() == num_face_dofs);

            for (size_t j = 0; j < num_face_dofs; j++) {
               const scalar_type qp_fphi_n = qp.weight() * mm_prod(fphi[j], n);
               for (size_t i = 0; i < div_basis_size; i++) {
                  BD(i, face_offset + j) += mm_prod(fphi_d[i], qp_fphi_n);
               }
            }
         }

         face_offset += num_face_dofs;
      }

      oper = MD.partialPivLu().solve(BD);
      data = BD.transpose() * oper;

      assert(oper.rows() == div_basis_size);
      assert(oper.cols() == dsr.total_size());

      assert(data.rows() == dsr.total_size());
      assert(data.cols() == dsr.total_size());
   }
};

} // end hho
} // end disk