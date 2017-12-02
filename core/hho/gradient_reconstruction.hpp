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
#include "hho/hho_stiffness_matrix.hpp"

namespace disk {

namespace hho {

// details on implementation of the gradient reconstruction
namespace priv {

template<bool b, typename BQData>
struct scalar_basis_and_gradient
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   static void impl(const mesh_type& msh,
                    const cell_type& cl,
                    const BQData&    bqd,
                    matrix_type&     oper,
                    matrix_type&     data)
   {}
};

template<typename BQData>
struct scalar_basis_and_gradient<true, BQData>
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   typedef material_tensor<scalar_type, mesh_type::dimension, mesh_type::dimension>
     material_tensor_type;

   template<typename TensorField>
   static void impl(const mesh_type&   msh,
                    const cell_type&   cl,
                    const TensorField& mtens,
                    const BQData&      bqd,
                    matrix_type&       oper,
                    matrix_type&       data)
   {
      const auto cell_degree     = bqd.cell_degree();
      const auto cell_basis_size = bqd.cell_basis.computed_size();
      const auto num_cell_dofs   = howmany_dofs(bqd.cell_basis);
      const auto face_basis_size = howmany_dofs(bqd.face_basis);

      matrix_type stiff_mat = matrix_type::Zero(cell_basis_size, cell_basis_size);

      const auto cell_quadpoints = bqd.cell_quadrature.integrate(msh, cl);
      for (auto& qp : cell_quadpoints) {
         const matrix_type dphi =
           bqd.cell_basis.eval_gradients(msh, cl, qp.point(), 0, cell_degree + 1);
         assert(dphi.rows() == cell_basis_size);

         stiff_mat += qp.weight() * dphi * (dphi * mtens(qp.point())).transpose();
      }

      /* LHS: take basis functions derivatives from degree 1 to K+1 */
      const auto        MG_rowcol_range = bqd.cell_basis.range(1, cell_degree + 1);
      const matrix_type MG              = take(stiff_mat, MG_rowcol_range, MG_rowcol_range);

      /* RHS, volumetric part. */
      const auto BG_row_range = bqd.cell_basis.range(1, cell_degree + 1);
      const auto BG_col_range = bqd.cell_basis.range(0, cell_degree);

      const auto fcs       = faces(msh, cl);
      const auto num_faces = fcs.size();

      const dofspace_ranges dsr(num_cell_dofs, face_basis_size, num_faces);

      matrix_type BG = matrix_type::Zero(BG_row_range.size(), dsr.total_size());

      BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
        take(stiff_mat, BG_row_range, BG_col_range);

      for (size_t face_i = 0; face_i < num_faces; face_i++) {
         const auto current_face_range = dsr.face_range(face_i);
         const auto fc                 = fcs[face_i];
         const auto n                  = normal(msh, cl, fc);
         const auto face_quadpoints    = bqd.face_quadrature.integrate(msh, fc);

         for (auto& qp : face_quadpoints) {
            const matrix_type c_phi =
              bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);
            const matrix_type c_dphi =
              bqd.cell_basis.eval_gradients(msh, cl, qp.point(), 1, cell_degree + 1);
            assert(c_phi.rows() == num_cell_dofs);

            const matrix_type c_dphi_n = (c_dphi * mtens(qp.point()).transpose()) * n;
            const matrix_type T        = qp.weight() * c_dphi_n * c_phi.transpose();

            BG.block(0, 0, BG.rows(), BG_col_range.size()) -= T;

            const matrix_type f_phi = bqd.face_basis.eval_functions(msh, fc, qp.point());
            assert(f_phi.size() == face_basis_size);
            const matrix_type F = qp.weight() * c_dphi_n * f_phi.transpose();

            BG.block(0, current_face_range.min(), BG.rows(), current_face_range.size()) += F;
         }
      }

      oper = MG.llt().solve(BG);    // GT
      data = BG.transpose() * oper; // A
   }

   static void impl(const mesh_type& msh,
                    const cell_type& cl,
                    const BQData&    bqd,
                    matrix_type&     oper,
                    matrix_type&     data)
   {
      material_tensor_type id_tens;
      const auto           tf = [](const typename mesh_type::point_type& pt) -> auto
      {
         return material_tensor_type::Identity();
      };

      impl(msh, cl, tf, bqd, oper, data);
   }
};

template<bool b, typename BQData>
struct scalar_basis_and_gradient_full
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   static void impl(const mesh_type& msh,
                    const cell_type& cl,
                    const BQData&    bqd,
                    matrix_type&     oper,
                    matrix_type&     data)
   {
      const size_t cell_degree     = bqd.cell_degree();
      const size_t face_degree     = bqd.face_degree();
      const size_t grad_degree     = bqd.grad_degree();
      const size_t cell_basis_size = (bqd.cell_basis.range(0, cell_degree)).size();
      const size_t face_basis_size = (bqd.face_basis.range(0, face_degree)).size();
      const size_t grad_basis_size = (bqd.grad_basis.range(0, grad_degree)).size();

      const auto   fcs       = faces(msh, cl);
      const size_t num_faces = fcs.size();

      const dofspace_ranges dsr(cell_basis_size, face_basis_size, num_faces);

      assert(dsr.total_size() == (cell_basis_size + num_faces * face_basis_size));

      // Compute mass matrix (gphi_i,ghi_j)
      matrix_type MG = matrix_type::Zero(grad_basis_size, grad_basis_size);

      const auto grad_quadpoints = bqd.grad_quadrature.integrate(msh, cl);
      for (auto& qp : grad_quadpoints) {
         const auto gphi = bqd.grad_basis.eval_functions(msh, cl, qp.point());
         assert(grad_basis_size == gphi.size());

         for (size_t j = 0; j < grad_basis_size; j++) {
            for (size_t i = j; i < grad_basis_size; i++) {
               MG(i, j) += qp.weight() * mm_prod(gphi[i], gphi[j]);
            }
         }
      }

      // lower part MG
      for (size_t i = 0; i < grad_basis_size; i++) {
         for (size_t j = i; j < grad_basis_size; j++) {
            MG(i, j) = MG(j, i);
         }
      }

      // Compute lhs
      matrix_type BG = matrix_type::Zero(grad_basis_size, dsr.total_size());

      // Compute (ghpi_i,grad(cphi_j))
      const auto grad_cell_quadpoints = bqd.grad_cell_quadrature.integrate(msh, cl);
      for (auto& qp : grad_cell_quadpoints) {
         const auto gphi  = bqd.grad_basis.eval_functions(msh, cl, qp.point());
         const auto cdphi = bqd.cell_basis.eval_gradients(msh, cl, qp.point());
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

         const auto grad_face_quadpoints = bqd.grad_face_max_quadrature.integrate(msh, fc);

         for (auto& qp : grad_face_quadpoints) {
            const auto c_phi = bqd.cell_basis.eval_functions(msh, cl, qp.point());
            const auto gphi  = bqd.grad_basis.eval_functions(msh, cl, qp.point());

            assert(grad_basis_size == gphi.size());
            assert(cell_basis_size == c_phi.size());
            // tau.n

            std::vector<static_vector<scalar_type, mesh_type::dimension>> gphi_n;

            gphi_n.reserve(gphi.size());

            for (size_t i = 0; i < gphi.size(); i++) {
               gphi_n.push_back(mm_prod(gphi[i], n));
            }

            assert(gphi_n.size() == grad_basis_size);

            matrix_type T = matrix_type::Zero(grad_basis_size, cell_basis_size);

            for (size_t j = 0; j < cell_basis_size; j++) {
               for (size_t i = 0; i < grad_basis_size; i++) {
                  T(i, j) = qp.weight() * mm_prod(gphi_n[i], c_phi[j]);
               }
            }

            assert(T.rows() == grad_basis_size);
            assert(T.cols() == cell_basis_size);

            BG.block(0, 0, grad_basis_size, cell_basis_size) -= T;

            const auto f_phi = bqd.face_basis.eval_functions(msh, fc, qp.point());
            assert(f_phi.size() == face_basis_size);
            matrix_type F = matrix_type::Zero(BG.rows(), current_face_range.size());

            for (size_t j = 0; j < current_face_range.size(); j++) {
               for (size_t i = 0; i < grad_basis_size; i++) {
                  F(i, j) = qp.weight() * mm_prod(gphi_n[i], f_phi[j]);
               }
            }

            assert(F.rows() == grad_basis_size);
            assert(F.cols() == current_face_range.size());

            BG.block(0, current_face_range.min(), grad_basis_size, current_face_range.size()) += F;
         }
      }

      oper = MG.llt().solve(BG);    // GT
      data = BG.transpose() * oper; // A

      assert(oper.rows() == grad_basis_size);
      assert(oper.cols() == dsr.total_size());
   }
};

template<typename BQData>
struct scalar_basis_and_gradient_full<true, BQData>
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   static void impl(const mesh_type& msh,
                    const cell_type& cl,
                    const BQData&    bqd,
                    matrix_type&     oper,
                    matrix_type&     data)
   {
      const size_t cell_degree     = bqd.cell_degree();
      const size_t face_degree     = bqd.face_degree();
      const size_t grad_degree     = bqd.grad_degree();
      const size_t cell_basis_size = (bqd.cell_basis.range(0, cell_degree)).size();
      const size_t face_basis_size = (bqd.face_basis.range(0, face_degree)).size();
      const size_t grad_basis_size = (bqd.grad_basis.range(0, grad_degree)).size();

      const auto   fcs       = faces(msh, cl);
      const size_t num_faces = fcs.size();

      const dofspace_ranges dsr(cell_basis_size, face_basis_size, num_faces);

      assert(dsr.total_size() == (cell_basis_size + num_faces * face_basis_size));

      // Compute mass matrix (gphi_i,ghi_j)
      matrix_type MG = matrix_type::Zero(grad_basis_size, grad_basis_size);

      const auto grad_quadpoints = bqd.grad_quadrature.integrate(msh, cl);
      for (auto& qp : grad_quadpoints) {
         const auto gphi = bqd.grad_basis.eval_functions(msh, cl, qp.point());
         assert(grad_basis_size == gphi.size());

         for (size_t j = 0; j < grad_basis_size; j++) {
            for (size_t i = j; i < grad_basis_size; i++) {
               MG(i, j) += qp.weight() * mm_prod(gphi[i], gphi[j]);
            }
         }
      }

      // lower part MG
      for (size_t i = 0; i < grad_basis_size; i++) {
         for (size_t j = i; j < grad_basis_size; j++) {
            MG(i, j) = MG(j, i);
         }
      }

      // Compute lhs
      matrix_type BG = matrix_type::Zero(grad_basis_size, dsr.total_size());

      // Compute (ghpi_i,grad(cphi_j))
      const auto grad_cell_quadpoints = bqd.grad_cell_quadrature.integrate(msh, cl);
      for (auto& qp : grad_cell_quadpoints) {
         const auto gphi = bqd.grad_basis.eval_functions(msh, cl, qp.point());
         const auto dphi = bqd.cell_basis.eval_gradients(msh, cl, qp.point());
         assert(grad_basis_size == gphi.size());
         assert(cell_basis_size == dphi.rows());
         assert(dphi.cols() == mesh_type::dimension);

         for (size_t j = 0; j < cell_basis_size; j++) {
            // on converti dphi_j
            static_vector<scalar_type, mesh_type::dimension> dphi_j =
              static_vector<scalar_type, mesh_type::dimension>::Zero();
            for (size_t k = 0; k < mesh_type::dimension; k++) {
               dphi_j(k) = dphi(j, k);
            }

            for (size_t i = 0; i < grad_basis_size; i++) {
               BG(i, j) += qp.weight() * mm_prod(gphi[i], dphi_j);
            }
         }
      } // end qp

      for (size_t face_i = 0; face_i < num_faces; face_i++) {
         const auto current_face_range = dsr.face_range(face_i);
         const auto fc                 = fcs[face_i];
         const auto n                  = normal(msh, cl, fc);

         const auto grad_face_quadpoints = bqd.grad_face_max_quadrature.integrate(msh, fc);

         for (auto& qp : grad_face_quadpoints) {
            const matrix_type c_phi = bqd.cell_basis.eval_functions(msh, cl, qp.point());
            const auto        gphi  = bqd.grad_basis.eval_functions(msh, cl, qp.point());

            assert(grad_basis_size == gphi.size());
            assert(cell_basis_size == c_phi.size());
            // tau.n
            matrix_type gphi_n = matrix_type::Zero(grad_basis_size, 1);
            for (size_t i = 0; i < grad_basis_size; i++)
               gphi_n(i, 0) = mm_prod(gphi[i], n);

            const matrix_type T = qp.weight() * gphi_n * c_phi.transpose();

            assert(T.rows() == grad_basis_size);
            assert(T.cols() == cell_basis_size);

            BG.block(0, 0, grad_basis_size, cell_basis_size) -= T;

            const matrix_type f_phi = bqd.face_basis.eval_functions(msh, fc, qp.point());
            assert(f_phi.size() == face_basis_size);
            const matrix_type F = qp.weight() * gphi_n * f_phi.transpose();

            assert(F.rows() == grad_basis_size);
            assert(F.cols() == current_face_range.size());

            BG.block(0, current_face_range.min(), grad_basis_size, current_face_range.size()) += F;
         }
      }

      oper = MG.llt().solve(BG);    // GT
      data = BG.transpose() * oper; // A

      assert(oper.rows() == grad_basis_size);
      assert(oper.cols() == dsr.total_size());
   }
};

template<bool b, typename BQData>
struct algo_selector_compute_gradient_full
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef typename BQData::cell_basis_type cell_basis_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   template<typename TensorField>
   static void impl(const mesh_type&   msh,
                    const cell_type&   cl,
                    const TensorField& mtens,
                    const BQData&      bqd,
                    matrix_type&       oper,
                    matrix_type&       data)
   {
      scalar_basis_and_gradient<is_scalar_basis<cell_basis_type>::value, BQData>::impl(
        msh, cl, mtens, bqd, oper, data);
   }

   static void impl(const mesh_type& msh,
                    const cell_type& cl,
                    const BQData&    bqd,
                    matrix_type&     oper,
                    matrix_type&     data)
   {
      scalar_basis_and_gradient<is_scalar_basis<cell_basis_type>::value, BQData>::impl(
        msh, cl, bqd, oper, data);
   }
};

template<typename BQData>
struct algo_selector_compute_gradient_full<true, BQData>
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef typename BQData::cell_basis_type cell_basis_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   static void impl(const mesh_type& msh,
                    const cell_type& cl,
                    const BQData&    bqd,
                    matrix_type&     oper,
                    matrix_type&     data)
   {
      scalar_basis_and_gradient_full<is_scalar_basis<cell_basis_type>::value, BQData>::impl(
        msh, cl, bqd, oper, data);
   }
};

} // priv

// class for the gradient reconstruction

template<typename BQData>
struct gradient_reconstruction_bq
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   const BQData& m_bqd;

 public:
   matrix_type oper;
   matrix_type data;

   gradient_reconstruction_bq() == delete;

   gradient_reconstruction_bq(const BQData& bqd)
     : m_bqd(bqd)
   {}

   template<typename TensorField>
   void compute(const mesh_type& msh, const cell_type& cl, const TensorField& mtens)
   {
      priv::algo_selector_compute_gradient_full<have_grad_space<BQData>::value, BQData>::impl(
        msh, cl, mtens, m_bqd, oper, data);
   }

   void compute(const mesh_type& msh, const cell_type& cl)
   {
      priv::algo_selector_compute_gradient_full<have_grad_space<BQData>::value, BQData>::impl(
        msh, cl, m_bqd, oper, data);
   }
};

} // hho
} // disk
