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
#include "hho/hho_stiffness_matrix.hpp"
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
   projector_bq(const BQData& bqd) : m_bqd(bqd) {}

   matrix_type cell_mm;
   matrix_type face_mm;
   matrix_type grad_mm;

   template<typename Function>
   vector_type
   projectOnCell(const mesh_type& msh, const cell_type& cl, const Function& func, int degree = -1)
   {
      if (degree < 0) degree = m_bqd.cell_degree();

      const auto cell_basis_size = m_bqd.cell_basis.range(0, degree).size();

      cell_mm               = mass_matrix(msh, cl, m_bqd, 0, degree);
      const vector_type rhs = compute_rhs(msh, cl, func, m_bqd, degree);

      assert(cell_mm.rows() == cell_mm.cols() && cell_mm.rows() == cell_basis_size);
      assert(rhs.rows() == cell_basis_size);

      return cell_mm.llt().solve(rhs);
   }

   template<typename Function>
   vector_type
   projectOnFace(const mesh_type& msh, const face_type& fc, const Function& func, int degree = -1)
   {
      if (degree < 0) degree = m_bqd.face_degree();

      const auto face_basis_size = m_bqd.face_basis.range(0, degree).size();

      face_mm               = mass_matrix(msh, fc, m_bqd, 0, degree);
      const vector_type rhs = compute_rhs(msh, fc, func, m_bqd, degree);

      assert(face_mm.rows() == face_mm.cols() && face_mm.rows() == face_basis_size);
      assert(rhs.rows() == face_basis_size);

      return face_mm.llt().solve(rhs);
   }

   template<typename Function>
   vector_type
   projectOnStiffnessSpace(const mesh_type& msh,
                           const cell_type& cl,
                           const Function&  f,
                           int              degree = -1)
   {
      typedef typename BQData::cell_basis_type::gradient_value_type gvt;

      if (degree < 0) degree = m_bqd.cell_degree() + 1;
      const auto grad_basis_size = m_bqd.cell_basis.range(1, degree).size();

      grad_mm = stiffness_matrix(msh, cl, m_bqd, 1, degree);
      assert(grad_mm.rows() == grad_mm.cols() && grad_mm.rows() == grad_basis_size);

      vector_type rhs = vector_type::Zero(grad_basis_size);

      const auto grad_quadpoints = m_bqd.cell_quadrature.integrate(msh, cl);
      assert(2 * degree <= m_bqd.cell_quadrature.order());

      for (auto& qp : grad_quadpoints) {
         const auto gphi = m_bqd.cell_basis.eval_gradients(msh, cl, qp.point(), 1, degree);
         assert(gphi.size() == grad_basis_size);

         const gvt f_qp = qp.weight() * f(qp.point());
         for (size_t i = 0; i < grad_basis_size; i++) {
            rhs(i) += mm_prod(f_qp, gphi[i]);
         }
      }

      return grad_mm.llt().solve(rhs);
   }

   template<typename Function>
   vector_type
   projectOnSymStiffnessSpace(const mesh_type& msh,
                              const cell_type& cl,
                              const Function&  f,
                              int              degree = -1)
   {
      if (degree < 0) degree = m_bqd.cell_degree() + 1;
      const auto grad_basis_size = m_bqd.cell_basis.range(1, degree).size();

      grad_mm = sym_stiffness_matrix(msh, cl, m_bqd, 1, degree);
      assert(grad_mm.rows() == grad_mm.cols() && grad_mm.rows() == grad_basis_size);

      vector_type rhs = vector_type::Zero(grad_basis_size);

      const auto grad_quadpoints = m_bqd.cell_quadrature.integrate(msh, cl);
      assert(2 * degree <= m_bqd.cell_quadrature.order());

      for (auto& qp : grad_quadpoints) {
         const auto gphi = m_bqd.cell_basis.eval_sgradients(msh, cl, qp.point(), 1, degree);
         assert(gphi.size() == grad_basis_size);

         const auto f_qp = mm_prod(qp.weight(), f(qp.point()));
         for (size_t i = 0; i < grad_basis_size; i++) {
            rhs(i) += mm_prod(f_qp, gphi[i]);
         }
      }

      return grad_mm.llt().solve(rhs);
   }

   template<typename Function>
   vector_type
   projectGradOnCell(const mesh_type& msh, const cell_type& cl, const Function& f, int degree = -1)
   {
      typedef typename BQData::grad_basis_type::function_value_type gvt;

      if (degree < 0) degree = m_bqd.grad_degree();
      const auto grad_basis_size = m_bqd.grad_basis.range(0, degree).size();

      grad_mm = grad_mass_matrix(msh, cl, m_bqd, 0, degree);
      assert(grad_mm.rows() == grad_mm.cols() && grad_mm.rows() == grad_basis_size);

      vector_type rhs = vector_type::Zero(grad_basis_size);

      const auto grad_quadpoints = m_bqd.grad_quadrature.integrate(msh, cl);
      assert(2 * degree <= m_bqd.grad_quadrature.order());

      for (auto& qp : grad_quadpoints) {
         const auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point(), 0, degree);
         assert(gphi.size() == grad_basis_size);

         const gvt f_qp = qp.weight() * f(qp.point());
         for (size_t i = 0; i < grad_basis_size; i++) {
            rhs(i) += mm_prod(f_qp, gphi[i]);
         }
      }

      return grad_mm.llt().solve(rhs);
   }
};

} // hho
} // disk
