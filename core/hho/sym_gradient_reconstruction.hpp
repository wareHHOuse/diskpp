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

#include "bases/bases_traits.hpp"
#include "common/eigen.hpp"
#include "hho/hho_bq.hpp"
#include "hho/hho_mass_matrix.hpp"
#include "hho/hho_stiffness_matrix.hpp"
#include "hho/hho_utils.hpp"
#include "quadratures/quadratures.hpp"

namespace disk {

namespace hho {

// details on implementation of the symetric gradient reconstruction

// class for the gradient reconstruction

template<typename BQData>
struct sym_gradient_reconstruction_bq
{
   typedef typename BQData::mesh_type       mesh_type;
   typedef typename mesh_type::scalar_type  scalar_type;
   typedef typename mesh_type::cell         cell_type;
   typedef typename BQData::cell_basis_type cell_basis_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   const static size_t dimension = mesh_type::dimension;

   const BQData& m_bqd;

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
   matrix_type oper;
   matrix_type data;

   sym_gradient_reconstruction_bq() = delete;

   sym_gradient_reconstruction_bq(const BQData& bqd) : m_bqd(bqd) {}

   void
   compute(const mesh_type& msh, const cell_type& cl)
   {
      const size_t cell_degree     = m_bqd.cell_degree();
      const size_t num_cell_dofs   = howmany_dofs(m_bqd.cell_basis);
      const size_t cell_basis_size = m_bqd.cell_basis.computed_size();

      const size_t face_degree   = m_bqd.face_degree();
      const size_t num_face_dofs = howmany_dofs(m_bqd.face_basis);

      const matrix_type stiff_mat = sym_stiffness_matrix(msh, cl, m_bqd, 0, cell_degree + 1);
      assert(stiff_mat.rows() == stiff_mat.cols() && stiff_mat.rows() == cell_basis_size);

      /* LHS: take basis functions derivatives from degree 1 to K+1 */
      const auto  MG_rowcol_range = m_bqd.cell_basis.range(1, cell_degree + 1);
      const auto  MG_size         = MG_rowcol_range.size() + num_lag();
      matrix_type MG              = matrix_type::Zero(MG_size, MG_size);
      MG.block(0, 0, MG_rowcol_range.size(), MG_rowcol_range.size()) =
        take(stiff_mat, MG_rowcol_range, MG_rowcol_range);

      // LHS: impose zero_average condition for the skew-symetric part
      const auto cell_quadpoints = m_bqd.cell_quadrature.integrate(msh, cl);
      assert(2 * cell_degree <= m_bqd.cell_quadrature.order());

      const auto SS = compute_SS();

      for (auto& qp : cell_quadpoints) {
         const auto c_dphi =
           m_bqd.cell_basis.eval_gradients(msh, cl, qp.point(), 1, cell_degree + 1);
         assert(c_dphi.size() == MG_rowcol_range.size());

         for (size_t j = 0; j < num_lag(); j++) {
            for (size_t i = 0; i < MG_rowcol_range.size(); i++) {
               const scalar_type qp_ss_cdphi = qp.weight() * mm_prod(SS[j], c_dphi[i]);
               MG(i, MG_rowcol_range.size() + j) += qp_ss_cdphi;
               MG(MG_rowcol_range.size() + j, i) += qp_ss_cdphi;
            }
         }
      }

      /* RHS, volumetric part. */
      const auto BG_row_range = m_bqd.cell_basis.range(1, cell_degree + 1);
      const auto BG_col_range = m_bqd.cell_basis.range(0, cell_degree);

      const auto fcs       = faces(msh, cl);
      const auto num_faces = fcs.size();

      const dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

      assert(dsr.total_size() == (num_cell_dofs + num_faces * num_face_dofs));

      matrix_type BG = matrix_type::Zero(BG_row_range.size() + num_lag(), dsr.total_size());

      BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
        take(stiff_mat, BG_row_range, BG_col_range);

      for (size_t face_i = 0; face_i < num_faces; face_i++) {
         const auto current_face_range = dsr.face_range(face_i);
         const auto fc                 = fcs[face_i];
         const auto n                  = normal(msh, cl, fc);

         const auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);
         assert((msh.dimension == 1 && m_bqd.face_quadrature.order() == 0) ||
                (cell_degree + face_degree) <= m_bqd.face_quadrature.order());

         for (auto& qp : face_quadpoints) {
            const auto c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);
            assert(c_phi.size() == num_cell_dofs);

            const auto c_sdphi =
              m_bqd.cell_basis.eval_sgradients(msh, cl, qp.point(), 1, cell_degree + 1);
            assert(c_sdphi.size() == BG_row_range.size());

            std::vector<static_vector<scalar_type, dimension>> qp_c_sdphi_n;
            qp_c_sdphi_n.reserve(BG_row_range.size());

            for (size_t i = 0; i < BG_row_range.size(); i++) {
               qp_c_sdphi_n.push_back(qp.weight() * mm_prod(c_sdphi[i], n));
            }
            assert(qp_c_sdphi_n.size() == BG_row_range.size());

            matrix_type T = matrix_type::Zero(BG_row_range.size(), num_cell_dofs);

            for (size_t j = 0; j < num_cell_dofs; j++) {
               for (size_t i = 0; i < BG_row_range.size(); i++) {
                  T(i, j) = mm_prod(qp_c_sdphi_n[i], c_phi[j]);
               }
            }

            BG.block(0, 0, BG_row_range.size(), num_cell_dofs) -= T;

            const auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point(), 0, face_degree);

            assert(f_phi.size() == num_face_dofs);

            matrix_type F = matrix_type::Zero(BG_row_range.size(), num_face_dofs);

            for (size_t j = 0; j < num_face_dofs; j++) {
               for (size_t i = 0; i < BG_row_range.size(); i++) {

                  F(i, j) = mm_prod(qp_c_sdphi_n[i], f_phi[j]);
               }
            }

            BG.block(0, current_face_range.min(), BG_row_range.size(), num_face_dofs) += F;
         }
      }

      assert(MG.rows() == MG.cols());
      assert(MG.cols() == BG.rows());

      Eigen::PartialPivLU<Eigen::MatrixXd> LU;

      LU.compute(MG);

      oper = (LU.solve(BG)).block(0, 0, MG_rowcol_range.size(), dsr.total_size());       // GT
      data = (BG.block(0, 0, BG_row_range.size(), dsr.total_size())).transpose() * oper; // A

      assert(oper.rows() == MG_rowcol_range.size() && oper.cols() == dsr.total_size());
      assert(data.rows() == dsr.total_size() && data.cols() == dsr.total_size());
   }
};

template<typename BQData>
struct sym_gradient_reconstruction_full_bq
{
   typedef typename BQData::mesh_type                    mesh_type;
   typedef typename mesh_type::scalar_type               scalar_type;
   typedef typename mesh_type::cell                      cell_type;
   typedef typename BQData::cell_basis_type              cell_basis_type;
   typedef typename cell_basis_type::function_value_type fvt;

   typedef dynamic_matrix<scalar_type> matrix_type;

   const BQData& m_bqd;

 public:
   matrix_type oper;
   matrix_type data;

   sym_gradient_reconstruction_full_bq() = delete;

   sym_gradient_reconstruction_full_bq(const BQData& bqd) : m_bqd(bqd) {}

   void
   compute(const mesh_type& msh, const cell_type& cl)
   {
      static_assert(have_grad_space<BQData>::value, "need to have gradient space");
      const size_t cell_degree     = m_bqd.cell_degree();
      const size_t face_degree     = m_bqd.face_degree();
      const size_t grad_degree     = m_bqd.grad_degree();
      const size_t cell_basis_size = (m_bqd.cell_basis.range(0, cell_degree)).size();
      const size_t face_basis_size = (m_bqd.face_basis.range(0, face_degree)).size();
      const size_t grad_basis_size = (m_bqd.grad_basis.range(0, grad_degree)).size();

      const auto   fcs       = faces(msh, cl);
      const size_t num_faces = fcs.size();

      const dofspace_ranges dsr(cell_basis_size, face_basis_size, num_faces);

      assert(dsr.total_size() == (cell_basis_size + num_faces * face_basis_size));

      // Compute mass matrix (gphi_i,ghi_j)
      const matrix_type MG = grad_mass_matrix(msh, cl, m_bqd, 0, grad_degree);
      assert(MG.rows() == MG.cols() && MG.rows() == grad_basis_size);

      // Compute lhs
      matrix_type BG = matrix_type::Zero(grad_basis_size, dsr.total_size());

      // Compute (ghpi_i,grad(cphi_j))
      const auto grad_cell_quadpoints = m_bqd.grad_cell_quadrature.integrate(msh, cl);
      assert((cell_degree + grad_degree) <= m_bqd.grad_cell_quadrature.order());

      for (auto& qp : grad_cell_quadpoints) {
         const auto gphi  = m_bqd.grad_basis.eval_functions(msh, cl, qp.point(), 0, grad_degree);
         const auto cdphi = m_bqd.cell_basis.eval_sgradients(msh, cl, qp.point(), 0, cell_degree);
         assert(grad_basis_size == gphi.size());
         assert(cell_basis_size == cdphi.size());

         for (size_t j = 0; j < cell_basis_size; j++) {
            for (size_t i = 0; i < grad_basis_size; i++) {
               BG(i, j) += qp.weight() * mm_prod(gphi[i], cdphi[j]);
            }
         }
      } // end qp

      for (size_t face_i = 0; face_i < num_faces; face_i++) {
         const auto current_face_range = dsr.face_range(face_i);
         const auto fc                 = fcs[face_i];
         const auto n                  = normal(msh, cl, fc);

         const auto grad_face_quadpoints = m_bqd.grad_face_max_quadrature.integrate(msh, fc);
         assert((msh.dimension == 1 && m_bqd.grad_face_max_quadrature.order() == 0) ||
                (std::max(cell_degree, face_degree) + grad_degree) <=
                  m_bqd.grad_face_max_quadrature.order());

         for (auto& qp : grad_face_quadpoints) {
            const auto c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);
            const auto gphi  = m_bqd.grad_basis.eval_functions(msh, cl, qp.point(), 0, grad_degree);

            assert(grad_basis_size == gphi.size());
            assert(cell_basis_size == c_phi.size());
            // tau.n

            std::vector<fvt> qp_gphi_n;

            qp_gphi_n.reserve(gphi.size());

            for (size_t i = 0; i < gphi.size(); i++) {
               qp_gphi_n.push_back(qp.weight() * mm_prod(gphi[i], n));
            }

            assert(qp_gphi_n.size() == grad_basis_size);

            matrix_type T = matrix_type::Zero(grad_basis_size, cell_basis_size);

            for (size_t j = 0; j < cell_basis_size; j++) {
               for (size_t i = 0; i < grad_basis_size; i++) {
                  T(i, j) = mm_prod(qp_gphi_n[i], c_phi[j]);
               }
            }

            BG.block(0, 0, grad_basis_size, cell_basis_size) -= T;

            const auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point(), 0, face_degree);
            assert(f_phi.size() == face_basis_size);

            matrix_type F = matrix_type::Zero(grad_basis_size, face_basis_size);

            for (size_t j = 0; j < face_basis_size; j++) {
               for (size_t i = 0; i < grad_basis_size; i++) {
                  F(i, j) = mm_prod(qp_gphi_n[i], f_phi[j]);
               }
            }

            BG.block(0, current_face_range.min(), grad_basis_size, face_basis_size) += F;
         }
      }

      oper = MG.llt().solve(BG);    // GT
      data = BG.transpose() * oper; // A

      assert(oper.rows() == grad_basis_size);
      assert(oper.cols() == dsr.total_size());

      assert(data.rows() == dsr.total_size());
      assert(data.cols() == dsr.total_size());
   }
};

} // hho
} // disk
