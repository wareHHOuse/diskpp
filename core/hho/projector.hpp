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

   const static size_t dimension = mesh_type::dimension;

   std::vector<static_matrix<scalar_type, dimension, dimension>>
   compute_SS()
   {
      assert(dimension == 2 || dimension == 3);

      std::vector<static_matrix<scalar_type, dimension, dimension>> ret;

      if (dimension == 2) {
         ret.reserve(1);
         static_matrix<scalar_type, dimension, dimension> mat =
           static_matrix<scalar_type, dimension, dimension>::Zero();

         mat << 0., 1., -1., 0.;

         ret.push_back(mat);
      } else if (dimension == 3) {
         ret.reserve(3);
         static_matrix<scalar_type, dimension, dimension> mat =
           static_matrix<scalar_type, dimension, dimension>::Zero();

         mat(0, 2) = 1;
         mat(2, 0) = -1;

         ret.push_back(mat);

         mat = static_matrix<scalar_type, dimension, dimension>::Zero();

         mat(0, 1) = 1;
         mat(1, 0) = -1;

         ret.push_back(mat);

         mat = static_matrix<scalar_type, dimension, dimension>::Zero();

         mat(1, 2) = 1;
         mat(2, 1) = -1;

         ret.push_back(mat);
      } else
         throw std::logic_error("dimension have to be =2 or =3");

      return ret;
   }

   size_t
   num_lag()
   {
      if (dimension == 2)
         return 1;
      else if (dimension == 3)
         return 3;

      throw std::logic_error("dimension have to be =2 or =3");
   }

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

   // project gradient on  gradient space
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

      // RHS
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

   // project gradient (! not the symetric gradient) on symetric gradient space
   template<typename Function>
   vector_type
   projectOnSymStiffnessSpace(const mesh_type& msh,
                              const cell_type& cl,
                              const Function&  f,
                              int              degree = -1)
   {
      typedef typename BQData::cell_basis_type::gradient_value_type gvt;

      if (degree < 0) degree = m_bqd.cell_degree() + 1;
      const auto grad_basis_size = m_bqd.cell_basis.range(1, degree).size();

      grad_mm = sym_stiffness_matrix(msh, cl, m_bqd, 1, degree);
      assert(grad_mm.rows() == grad_mm.cols() && grad_mm.rows() == grad_basis_size);

      /* LHS: take basis functions derivatives from degree 1 to K+1 */
      const auto  MG_rowcol_range                      = m_bqd.cell_basis.range(1, degree);
      const auto  MG_size                              = MG_rowcol_range.size() + num_lag();
      matrix_type MG                                   = matrix_type::Zero(MG_size, MG_size);
      MG.block(0, 0, grad_basis_size, grad_basis_size) = grad_mm;

      // LHS: impose zero_average condition for the skew-symetric part
      const auto cell_quadpoints = m_bqd.cell_quadrature.integrate(msh, cl);
      assert(2 * (degree - 1) <= m_bqd.cell_quadrature.order());

      const auto SS = compute_SS();

      for (auto& qp : cell_quadpoints) {
         const auto c_dphi = m_bqd.cell_basis.eval_gradients(msh, cl, qp.point(), 1, degree);
         assert(c_dphi.size() == MG_rowcol_range.size());

         for (size_t j = 0; j < num_lag(); j++) {
            for (size_t i = 0; i < MG_rowcol_range.size(); i++) {
               const scalar_type qp_ss_cdphi = qp.weight() * mm_prod(SS[j], c_dphi[i]);
               MG(i, MG_rowcol_range.size() + j) += qp_ss_cdphi;
               MG(MG_rowcol_range.size() + j, i) += qp_ss_cdphi;
            }
         }
      }

      // RHS
      vector_type rhs = vector_type::Zero(grad_basis_size + num_lag());

      const auto grad_quadpoints = m_bqd.cell_quadrature.integrate(msh, cl);
      assert(2 * degree <= m_bqd.cell_quadrature.order());

      for (auto& qp : grad_quadpoints) {
         const auto gphi = m_bqd.cell_basis.eval_sgradients(msh, cl, qp.point(), 1, degree);
         assert(gphi.size() == grad_basis_size);

         const gvt f_qp = qp.weight() * f(qp.point());
         for (size_t i = 0; i < grad_basis_size; i++) {
            rhs(i) += mm_prod(f_qp, gphi[i]);
         }
      }

      Eigen::PartialPivLU<Eigen::MatrixXd> LU;

      LU.compute(MG);

      return LU.solve(rhs).head(grad_basis_size);
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
