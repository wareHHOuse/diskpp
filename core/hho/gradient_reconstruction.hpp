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

// details on implementation of the gradient reconstruction
namespace priv {

template<bool b, typename BQData>
struct use_vector_container_and_gradient // false for scalar basis functions
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   typedef material_tensor<scalar_type, mesh_type::dimension, mesh_type::dimension>
     material_tensor_type;

   template<typename TensorField>
   static void
   impl(const mesh_type&   msh,
        const cell_type&   cl,
        const TensorField& mtens,
        const BQData&      bqd,
        matrix_type&       oper,
        matrix_type&       data)
   {
      const auto cell_degree     = bqd.cell_degree();
      const auto cell_basis_size = bqd.cell_basis.computed_size();
      const auto num_cell_dofs   = howmany_dofs(bqd.cell_basis);

      const auto face_degree     = bqd.face_degree();
      const auto face_basis_size = howmany_dofs(bqd.face_basis);

      const matrix_type stiff_mat = stiffness_matrix(msh, cl, mtens, bqd, 0, cell_degree + 1);
      assert(stiff_mat.rows() == stiff_mat.cols() && stiff_mat.rows() == cell_basis_size);

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

         const auto face_quadpoints = bqd.face_quadrature.integrate(msh, fc);
         assert((msh.dimension == 1 && bqd.face_quadrature.order() == 0) ||
                (cell_degree + face_degree) <= bqd.face_quadrature.order());

         for (auto& qp : face_quadpoints) {
            const matrix_type c_phi =
              bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);
            assert(c_phi.rows() == num_cell_dofs);

            const matrix_type c_dphi =
              bqd.cell_basis.eval_gradients(msh, cl, qp.point(), 1, cell_degree + 1);
            assert(c_dphi.rows() == BG_row_range.size());

            const matrix_type qp_c_dphi_n =
              qp.weight() * (c_dphi * mtens(qp.point()).transpose()) * n;
            const matrix_type T = qp_c_dphi_n * c_phi.transpose();

            BG.block(0, 0, BG.rows(), BG_col_range.size()) -= T;

            const matrix_type f_phi =
              bqd.face_basis.eval_functions(msh, fc, qp.point(), 0, face_degree);
            assert(f_phi.size() == face_basis_size);

            const matrix_type F = qp_c_dphi_n * f_phi.transpose();

            BG.block(0, current_face_range.min(), BG.rows(), current_face_range.size()) += F;
         }
      }

      oper = MG.llt().solve(BG);    // GT
      data = BG.transpose() * oper; // A

      assert(oper.rows() == MG_rowcol_range.size() && oper.cols() == dsr.total_size());
      assert(data.rows() == dsr.total_size() && data.cols() == dsr.total_size());
   }

   static void
   impl(const mesh_type& msh,
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

template<typename BQData>
struct use_vector_container_and_gradient<true, BQData>
{
   typedef typename BQData::mesh_type                            mesh_type;
   typedef typename mesh_type::scalar_type                       scalar_type;
   typedef typename mesh_type::cell                              cell_type;
   typedef typename BQData::cell_basis_type::function_value_type fvt;

   const static size_t                 dimension = mesh_type::dimension;
   typedef dynamic_matrix<scalar_type> matrix_type;

   static void
   impl(const mesh_type& msh,
        const cell_type& cl,
        const BQData&    bqd,
        matrix_type&     oper,
        matrix_type&     data)
   {
      const size_t cell_degree     = bqd.cell_degree();
      const size_t num_cell_dofs   = howmany_dofs(bqd.cell_basis);
      const size_t cell_basis_size = bqd.cell_basis.computed_size();

      const size_t face_degree   = bqd.face_degree();
      const size_t num_face_dofs = howmany_dofs(bqd.face_basis);

      const matrix_type stiff_mat = stiffness_matrix(msh, cl, bqd, 0, cell_degree + 1);
      assert(stiff_mat.rows() == stiff_mat.cols() && stiff_mat.rows() == cell_basis_size);

      /* LHS: take basis functions derivatives from degree 1 to K+1 */
      const auto        MG_rowcol_range = bqd.cell_basis.range(1, cell_degree + 1);
      const matrix_type MG              = take(stiff_mat, MG_rowcol_range, MG_rowcol_range);

      /* RHS, volumetric part. */
      const auto BG_row_range = bqd.cell_basis.range(1, cell_degree + 1);
      const auto BG_col_range = bqd.cell_basis.range(0, cell_degree);

      const auto fcs       = faces(msh, cl);
      const auto num_faces = fcs.size();

      const dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

      assert(dsr.total_size() == (num_cell_dofs + num_faces * num_face_dofs));

      matrix_type BG = matrix_type::Zero(BG_row_range.size(), dsr.total_size());

      BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
        take(stiff_mat, BG_row_range, BG_col_range);

      for (size_t face_i = 0; face_i < num_faces; face_i++) {
         const auto current_face_range = dsr.face_range(face_i);
         const auto fc                 = fcs[face_i];
         const auto n                  = normal(msh, cl, fc);

         const auto face_quadpoints = bqd.face_quadrature.integrate(msh, fc);
         assert((msh.dimension == 1 && bqd.face_quadrature.order() == 0) ||
                (cell_degree + face_degree) <= bqd.face_quadrature.order());

         for (auto& qp : face_quadpoints) {
            const auto c_phi = bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);
            assert(c_phi.size() == num_cell_dofs);

            const auto c_dphi =
              bqd.cell_basis.eval_gradients(msh, cl, qp.point(), 1, cell_degree + 1);
            assert(c_dphi.size() == BG_row_range.size());

            std::vector<fvt> qp_c_dphi_n;
            qp_c_dphi_n.reserve(BG_row_range.size());

            for (size_t i = 0; i < BG_row_range.size(); i++) {
               qp_c_dphi_n.push_back(qp.weight() * mm_prod(c_dphi[i], n));
            }
            assert(qp_c_dphi_n.size() == BG.rows());

            matrix_type T = matrix_type::Zero(BG.rows(), num_cell_dofs);

            for (size_t j = 0; j < num_cell_dofs; j++) {
               for (size_t i = 0; i < BG.rows(); i++) {
                  T(i, j) = mm_prod(qp_c_dphi_n[i], c_phi[j]);
               }
            }

            BG.block(0, 0, BG.rows(), num_cell_dofs) -= T;

            const auto f_phi = bqd.face_basis.eval_functions(msh, fc, qp.point(), 0, face_degree);

            assert(f_phi.size() == num_face_dofs);

            matrix_type F = matrix_type::Zero(BG.rows(), num_face_dofs);

            for (size_t j = 0; j < num_face_dofs; j++) {
               for (size_t i = 0; i < BG.rows(); i++) {

                  F(i, j) = mm_prod(qp_c_dphi_n[i], f_phi[j]);
               }
            }

            BG.block(0, current_face_range.min(), BG.rows(), num_face_dofs) += F;
         }
      }

      assert(MG.rows() == MG.cols());
      assert(MG.cols() == BG.rows());

      oper = MG.ldlt().solve(BG);   // GT
      data = BG.transpose() * oper; // A

      assert(oper.rows() == MG_rowcol_range.size() && oper.cols() == dsr.total_size());
      assert(data.rows() == dsr.total_size() && data.cols() == dsr.total_size());
   }
};

template<bool b, typename BQData>
struct use_vector_container_and_gradient_full // false for scalar basis functions
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef dynamic_matrix<scalar_type> matrix_type;
   typedef dynamic_vector<scalar_type> vector_type;

   static void
   impl(const mesh_type& msh,
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
      const matrix_type MG = grad_mass_matrix(msh, cl, bqd, 0, grad_degree);
      assert(MG.rows() == MG.cols() && MG.rows() == grad_basis_size);

      // Compute lhs
      matrix_type BG = matrix_type::Zero(grad_basis_size, dsr.total_size());

      // Compute (ghpi_i,grad(cphi_j))
      const auto grad_cell_quadpoints = bqd.grad_cell_quadrature.integrate(msh, cl);
      assert((cell_degree + grad_degree) <= bqd.grad_cell_quadrature.order());

      for (auto& qp : grad_cell_quadpoints) {
         const auto gphi = bqd.grad_basis.eval_functions(msh, cl, qp.point(), 0, grad_degree);
         const auto dphi = bqd.cell_basis.eval_gradients(msh, cl, qp.point(), 0, cell_degree);
         assert(grad_basis_size == gphi.size());
         assert(cell_basis_size == dphi.rows());
         assert(dphi.cols() == mesh_type::dimension);

         // we covert dphi since we can have mismatch type
         const auto dphi_j = convert_to_vector(dphi);

         for (size_t j = 0; j < cell_basis_size; j++) {
            const auto qp_dphi_j = mm_prod(qp.weight(), dphi_j[j]);
            for (size_t i = 0; i < grad_basis_size; i++) {
               BG(i, j) += mm_prod(gphi[i], qp_dphi_j);
            }
         }
      } // end qp

      for (size_t face_i = 0; face_i < num_faces; face_i++) {
         const auto current_face_range = dsr.face_range(face_i);
         const auto fc                 = fcs[face_i];
         const auto n                  = normal(msh, cl, fc);

         const auto grad_face_quadpoints = bqd.grad_face_max_quadrature.integrate(msh, fc);
         assert((msh.dimension == 1 && bqd.face_quadrature.order() == 0) ||
                ((std::max(cell_degree, face_degree) + grad_degree) <=
                 bqd.grad_face_max_quadrature.order()));

         for (auto& qp : grad_face_quadpoints) {
            const matrix_type c_phi =
              bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);
            const auto gphi = bqd.grad_basis.eval_functions(msh, cl, qp.point(), 0, grad_degree);

            assert(grad_basis_size == gphi.size());
            assert(cell_basis_size == c_phi.size());
            // tau.n
            vector_type qp_gphi_n = vector_type::Zero(grad_basis_size);
            for (size_t i = 0; i < grad_basis_size; i++)
               qp_gphi_n(i) = qp.weight() * mm_prod(gphi[i], n);

            const matrix_type T = qp_gphi_n * c_phi.transpose();

            BG.block(0, 0, grad_basis_size, cell_basis_size) -= T;

            const matrix_type f_phi =
              bqd.face_basis.eval_functions(msh, fc, qp.point(), 0, face_degree);
            assert(f_phi.size() == face_basis_size);
            const matrix_type F = qp_gphi_n * f_phi.transpose();

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

template<typename BQData>
struct use_vector_container_and_gradient_full<true, BQData>
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   typedef typename BQData::cell_basis_type::function_value_type fvt;

   static void
   impl(const mesh_type& msh,
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
      const matrix_type MG = grad_mass_matrix(msh, cl, bqd, 0, grad_degree);
      assert(MG.rows() == MG.cols() && MG.rows() == grad_basis_size);

      // Compute lhs
      matrix_type BG = matrix_type::Zero(grad_basis_size, dsr.total_size());

      // Compute (ghpi_i,grad(cphi_j))
      const auto grad_cell_quadpoints = bqd.grad_cell_quadrature.integrate(msh, cl);
      assert((cell_degree + grad_degree) <= bqd.grad_cell_quadrature.order());

      for (auto& qp : grad_cell_quadpoints) {
         const auto gphi  = bqd.grad_basis.eval_functions(msh, cl, qp.point(), 0, grad_degree);
         const auto cdphi = bqd.cell_basis.eval_gradients(msh, cl, qp.point(), 0, cell_degree);
         assert(grad_basis_size == gphi.size());
         assert(cell_basis_size == cdphi.size());

         for (size_t j = 0; j < cell_basis_size; j++) {
            const auto qp_cdphi_j = mm_prod(qp.weight(), cdphi[j]);
            for (size_t i = 0; i < grad_basis_size; i++) {
               BG(i, j) += mm_prod(gphi[i], qp_cdphi_j);
            }
         }
      } // end qp

      for (size_t face_i = 0; face_i < num_faces; face_i++) {
         const auto current_face_range = dsr.face_range(face_i);
         const auto fc                 = fcs[face_i];
         const auto n                  = normal(msh, cl, fc);

         const auto grad_face_quadpoints = bqd.grad_face_max_quadrature.integrate(msh, fc);
         assert((msh.dimension == 1 && bqd.grad_face_max_quadrature.order() == 0) ||
                (std::max(cell_degree, face_degree) + grad_degree) <=
                  bqd.grad_face_max_quadrature.order());

         for (auto& qp : grad_face_quadpoints) {
            const auto c_phi = bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);
            const auto gphi  = bqd.grad_basis.eval_functions(msh, cl, qp.point(), 0, grad_degree);

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

            const auto f_phi = bqd.face_basis.eval_functions(msh, fc, qp.point(), 0, face_degree);
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

} // priv

// class for the gradient reconstruction

template<typename BQData>
struct gradient_reconstruction_bq
{
   typedef typename BQData::mesh_type       mesh_type;
   typedef typename mesh_type::scalar_type  scalar_type;
   typedef typename mesh_type::cell         cell_type;
   typedef typename BQData::cell_basis_type cell_basis_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   const BQData& m_bqd;

 public:
   matrix_type oper;
   matrix_type data;

   gradient_reconstruction_bq() = delete;

   gradient_reconstruction_bq(const BQData& bqd) : m_bqd(bqd) {}

   template<typename TensorField>
   void
   compute(const mesh_type& msh, const cell_type& cl, const TensorField& mtens)
   {
      priv::use_vector_container_and_gradient<use_vector_container<cell_basis_type>::value,
                                              BQData>::impl(msh, cl, mtens, m_bqd, oper, data);
   }

   void
   compute(const mesh_type& msh, const cell_type& cl)
   {
      priv::use_vector_container_and_gradient<use_vector_container<cell_basis_type>::value,
                                              BQData>::impl(msh, cl, m_bqd, oper, data);
   }
};

template<typename BQData>
struct gradient_reconstruction_full_bq
{
   typedef typename BQData::mesh_type       mesh_type;
   typedef typename mesh_type::scalar_type  scalar_type;
   typedef typename mesh_type::cell         cell_type;
   typedef typename BQData::cell_basis_type cell_basis_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   const BQData& m_bqd;

 public:
   matrix_type oper;
   matrix_type data;

   gradient_reconstruction_full_bq() = delete;

   gradient_reconstruction_full_bq(const BQData& bqd) : m_bqd(bqd) {}

   void
   compute(const mesh_type& msh, const cell_type& cl)
   {
      static_assert(have_grad_space<BQData>::value, "need to have gradient space");
      priv::use_vector_container_and_gradient_full<use_vector_container<cell_basis_type>::value,
                                                   BQData>::impl(msh, cl, m_bqd, oper, data);
   }
};

} // hho
} // disk
