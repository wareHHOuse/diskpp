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

#include "LerayLions_law.hpp"
#include "bases/bases_utils.hpp"
#include "common/eigen.hpp"
#include "hho/hho_bq.hpp"
#include "hho/hho_utils.hpp"
#include "timecounter.h"

namespace LerayLions {

template<typename BQData>
class LerayLions
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef typename BQData::grad_basis_type::function_value_type gvt;

   typedef dynamic_matrix<scalar_type> matrix_type;
   typedef dynamic_vector<scalar_type> vector_type;

   const BQData& m_bqd;

 public:
   matrix_type K_int;
   vector_type RTF;
   double      time_law;

   LerayLions(const BQData& bqd) : m_bqd(bqd) {}

   template<typename Function>
   void
   compute(const mesh_type&   msh,
           const cell_type&   cl,
           const Function&    load,
           const matrix_type& GT,
           const vector_type& uTF,
           const scalar_type  leray_param)
   {
      const size_t DIM         = msh.dimension;
      const size_t cell_degree = m_bqd.cell_degree();
      const size_t grad_degree = m_bqd.grad_degree();

      const size_t cell_basis_size = howmany_dofs(m_bqd.cell_basis);
      const size_t grad_basis_size = (m_bqd.grad_basis.range(0, grad_degree)).size();
      const size_t face_basis_size = howmany_dofs(m_bqd.face_basis);

      time_law = 0.0;
      timecounter tc;

      const auto   fcs            = faces(msh, cl);
      const size_t num_faces      = fcs.size();
      const size_t num_total_dofs = cell_basis_size + num_faces * face_basis_size;

      matrix_type AT = matrix_type::Zero(grad_basis_size, grad_basis_size);
      vector_type aT = vector_type::Zero(grad_basis_size);

      RTF = vector_type::Zero(num_total_dofs);

      const vector_type GT_uTF = GT * uTF;

      pLaplaceLaw<scalar_type> law(leray_param);

      const auto grad_quadpoints = m_bqd.grad_quadrature.integrate(msh, cl);
      assert(2 * grad_degree <= m_bqd.grad_quadrature.order());

      for (auto& qp : grad_quadpoints) {
         const auto c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);
         const auto gphi  = m_bqd.grad_basis.eval_functions(msh, cl, qp.point(), 0, grad_degree);

         assert(c_phi.size() == cell_basis_size && gphi.size() == grad_basis_size);

         // Compute local gradient and norm
         const auto GT_iqn = disk::hho::eval(GT_uTF, gphi);

         // Compute bahavior
         tc.tic();
         const auto tensor_behavior = law.compute_whole(GT_iqn);
         tc.toc();
         time_law += tc.to_double();

         for (size_t j = 0; j < grad_basis_size; j++) {
            const gvt Agphi_j = qp.weight() * disk::mm_prod(tensor_behavior.second, gphi[j]);
            for (size_t i = j; i < grad_basis_size; i++) {
               // compute (Gkt v, A(u) : Gkt du)
               AT(i, j) += disk::mm_prod(gphi[i], Agphi_j);
            }
         }

         // compute (PK1(u), G^k_T v)_T
         for (size_t i = 0; i < grad_basis_size; i++) {
            aT(i) += qp.weight() * disk::mm_prod(tensor_behavior.first, gphi[i]);
         }
      }

      // compute (f,v)_T
      RTF.block(0, 0, cell_basis_size, 1) =
        disk::hho::compute_rhs(msh, cl, load, m_bqd, cell_degree);

      // lower part AT
      for (size_t j = 0; j < grad_basis_size; j++)
         for (size_t i = 0; i <= j; i++)
            AT(i, j) = AT(j, i);

      K_int = GT.transpose() * AT * GT;
      RTF -= GT.transpose() * aT;

      assert(K_int.rows() == num_total_dofs);
      assert(K_int.cols() == num_total_dofs);
      assert(RTF.rows() == num_total_dofs);
   }
};

} // end namespace LerayLions