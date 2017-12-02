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
#include "hho/hho_utils.hpp"

namespace disk {

namespace hho {

template<typename BQData>
class projector_bq
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef typename mesh_type::face        face_type;

   typedef dynamic_matrix<scalar_type> matrix_type;
   typedef dynamic_vector<scalar_type> vector_type;

   const BQData& m_bqd;

 public:
   projector_bq(const BQData& bqd)
     : m_bqd(bqd)
   {}

   matrix_type cell_mm;
   matrix_type face_mm;
   matrix_type grad_mm;

   template<typename Function>
   vector_type projectOnCell(const mesh_type& msh,
                             const cell_type& cl,
                             const Function&  func,
                             int              degree = -1)
   {
      if (degree < 0) degree = m_bqd.cell_degree();

      cell_mm               = mass_matrix(msh, cl, m_bqd, degree);
      const vector_type rhs = compute_rhs(msh, cl, func, m_bqd, degree);

      return cell_mm.llt().solve(rhs);
   }

   template<typename Function>
   vector_type projectOnFace(const mesh_type& msh,
                             const face_type& fc,
                             const Function&  func,
                             int              degree = -1)
   {
      if (degree < 0) degree = m_bqd.cell_degree();

      face_mm               = mass_matrix(msh, fc, m_bqd, degree);
      const vector_type rhs = compute_rhs(msh, fc, func, m_bqd, degree);

      return face_mm.llt().solve(rhs);
   }

   template<typename Function>
   vector_type projectGradOnCell(const mesh_type& msh,
                                 const cell_type& cl,
                                 const Function&  f,
                                 int              degree = -1)
   {
      const size_t grad_basis_size = m_bqd.grad_basis.size();

      // on doit enore creer un compute rhs et mass du grad  !!!
      matrix_type mm  = matrix_type::Zero(grad_basis_size, grad_basis_size);
      vector_type rhs = vector_type::Zero(grad_basis_size);

      const auto grad_quadpoints = m_bqd.grad_quadrature.integrate(msh, cl);
      for (auto& qp : grad_quadpoints) {
         auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());

         for (size_t j = 0; j < grad_basis_size; j++) {
            for (size_t i = j; i < grad_basis_size; i++) {
               mm(i, j) += qp.weight() * mm_prod(gphi[i], gphi[j]);
            }
         }

         for (size_t i = 0; i < grad_basis_size; i++) {
            rhs(i) += qp.weight() * mm_prod(f(qp.point()), gphi[i]);
         }
      }

      // lower part
      for (size_t i = 0; i < grad_basis_size; i++)
         for (size_t j = i; j < grad_basis_size; j++)
            mm(i, j) = mm(j, i);

      grad_mm = mm;
      return grad_mm.llt().solve(rhs);
   }
};

} // hho
} // disk
