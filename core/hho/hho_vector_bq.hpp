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

#include "bases/bases_ranges.hpp"
#include "bases/bases_utils.hpp"
#include "common/eigen.hpp"
#include "mechanics/BoundaryConditions.hpp"
#include "timecounter.h"

//#define USE_BLAS
#define FILL_COLMAJOR

namespace disk {

namespace hho {

template<typename BQData>
class gradient_reconstruction_vector_bq
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef typename BQData::cell_basis_type cell_basis_type;
   typedef typename BQData::cell_quad_type  cell_quad_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   const BQData& m_bqd;

   cell_basis_type cell_basis;
   cell_quad_type  cell_quadrature;

 public:
   matrix_type oper;
   matrix_type data;

   gradient_reconstruction_vector_bq(const BQData& bqd)
     : m_bqd(bqd)
   {
      cell_basis      = cell_basis_type(bqd.cell_degree() + 1);
      cell_quadrature = cell_quad_type(2 * (bqd.cell_degree() + 1));
   }

   void compute(const mesh_type& msh, const cell_type& cl)
   {
      const size_t cell_degree     = m_bqd.cell_degree();
      const size_t face_degree     = m_bqd.face_degree();
      const size_t cell_basis_size = (cell_basis.range(0, cell_degree + 1)).size();
      const size_t face_basis_size = m_bqd.face_basis.range(0, face_degree).size();

      matrix_type stiff_mat = matrix_type::Zero(cell_basis_size, cell_basis_size);

      const auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
      for (auto& qp : cell_quadpoints) {
         const auto dphi = cell_basis.eval_gradients(msh, cl, qp.point());
         assert(cell_basis_size == dphi.size());

         for (size_t i = 0; i < cell_basis_size; i++) {
            for (size_t j = i; j < cell_basis_size; j++) {
               stiff_mat(i, j) += qp.weight() * mm_prod(dphi[i], dphi[j]);
            }
         }
      }

      // lower part
      for (size_t i = 1; i < cell_basis_size; i++)
         for (size_t j = 0; j < i; j++)
            stiff_mat(i, j) = stiff_mat(j, i);

      /* LHS: take basis functions derivatives from degree 1 to K+1 */
      const auto        MG_rowcol_range = cell_basis.range(1, cell_degree + 1);
      const matrix_type MG              = take(stiff_mat, MG_rowcol_range, MG_rowcol_range);

      /* RHS, volumetric part. */
      const auto BG_row_range = cell_basis.range(1, cell_degree + 1);
      const auto BG_col_range = cell_basis.range(0, cell_degree);

      const auto   fcs       = faces(msh, cl);
      const size_t num_faces = fcs.size();

      const size_t num_cell_dofs = BG_col_range.size();

      const dofspace_ranges dsr(num_cell_dofs, face_basis_size, num_faces);

      assert(dsr.total_size() == (num_cell_dofs + num_faces * face_basis_size));

      matrix_type BG = matrix_type::Zero(BG_row_range.size(), dsr.total_size());

      BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
        take(stiff_mat, BG_row_range, BG_col_range);

      for (size_t face_i = 0; face_i < num_faces; face_i++) {
         const auto current_face_range = dsr.face_range(face_i);
         const auto fc                 = fcs[face_i];
         const auto n                  = normal(msh, cl, fc);

         const auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);

         for (auto& qp : face_quadpoints) {
            auto       c_phi  = cell_basis.eval_functions(msh, cl, qp.point()); // 0, m_degree);
            const auto c_dphi = cell_basis.eval_gradients(msh, cl, qp.point()); // 1, m_degree+1);

            assert(c_phi.size() == cell_basis_size);

            decltype(c_phi) c_dphi_n;

            c_dphi_n.reserve(BG_row_range.to() - BG_row_range.from());

            for (size_t i = BG_row_range.from(); i < BG_row_range.to(); i++) {
               c_dphi_n.push_back(mm_prod(c_dphi[i], n));
            }

            assert(c_dphi_n.size() == BG.rows());

            matrix_type T = matrix_type::Zero(BG.rows(), BG_col_range.size());

            for (size_t i = 0; i < BG.rows(); i++) {
               for (size_t j = 0; j < BG_col_range.size(); j++) {
                  T(i, j) = qp.weight() * mm_prod(c_dphi_n[i], c_phi[j]);
               }
            }

            BG.block(0, 0, BG.rows(), BG_col_range.size()) -= T;

            const auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());

            assert(f_phi.size() == face_basis_size);

            matrix_type F = matrix_type::Zero(BG.rows(), current_face_range.size());

            for (size_t i = 0; i < BG.rows(); i++) {
               for (size_t j = 0; j < current_face_range.size(); j++) {
                  F(i, j) = qp.weight() * mm_prod(c_dphi_n[i], f_phi[j]);
               }
            }

            BG.block(0, current_face_range.min(), BG.rows(), current_face_range.size()) += F;
         }
      }

      assert(MG.rows() == MG.cols());
      assert(MG.cols() == BG.rows());

      oper = MG.ldlt().solve(BG);   // GT
      data = BG.transpose() * oper; // A
   }
};

template<typename BQData>
class gradient_reconstruction_vector_full_bq
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef dynamic_matrix<scalar_type> matrix_type;
   typedef dynamic_vector<scalar_type> vector_type;

   const BQData& m_bqd;

   template<int DIM>
   std::vector<static_vector<scalar_type, DIM>> compute_gphi_n(
     const std::vector<static_matrix<scalar_type, DIM, DIM>> gphi,
     const static_vector<scalar_type, DIM>                   n) const
   {
      const size_t grad_basis_size = gphi.size();
      const size_t DIM2            = DIM * DIM;

      std::vector<static_vector<scalar_type, DIM>> gphi_n;

      static_vector<scalar_type, DIM> vzero = static_vector<scalar_type, DIM>::Zero();

      gphi_n.reserve(gphi.size());

      for (size_t i = 0; i < grad_basis_size; i += DIM2) {
         for (size_t k = 0; k < DIM; k++) { // depend de l'ordre des bases
            scalar_type val = gphi[i](0, 0) * n(k);
            for (size_t l = 0; l < DIM; l++) { // depend de l'ordre des bases
               auto gn(vzero);
               gn(l) = val;
               gphi_n.push_back(gn);
            }
         }
      }

      //             for(size_t i = 0; i < gphi.size(); i++)
      //             {
      //                gphi_n.push_back(mm_prod(gphi[i], n));
      //             }

      return gphi_n;
   }

   template<int DIM>
   void compute_MG(matrix_type&                                            MG,
                   const std::vector<static_matrix<scalar_type, DIM, DIM>> gphi,
                   const scalar_type                                       weight) const
   {
      const size_t grad_basis_size = gphi.size();
      const size_t grad_degree     = m_bqd.grad_degree();
      const size_t poly_space =
        grad_basis_size; // DIM * DIM * binomial(grad_degree-1 + DIM, grad_degree-1);
      const size_t DIM2 = DIM * DIM;

      //             // poly classique
      //
      size_t col = 0;
      for (size_t j = 0; j < poly_space; j += DIM2) {
         for (size_t k = 0; k < DIM; k++) {    // depend de l'ordre des bases
            for (size_t l = 0; l < DIM; l++) { // depend de l'ordre des bases
               for (size_t i = col; i < poly_space; i += DIM2) {
                  MG(i, col) += weight * gphi[i](l, k) * gphi[col](l, k);
               }
               col++;
            }
         }
      }

      // RT
      //             for(std::size_t j = 0; j < grad_basis_size; j ++) {
      //                for(std::size_t i = j; i < grad_basis_size; i ++) {
      //                   MG(i,j) = mm_prod(gphi[i], gphi[j]);
      //                }
      //             }
   }

   template<int DIM>
   void compute_BG(matrix_type&                                            BG,
                   const std::vector<static_matrix<scalar_type, DIM, DIM>> gphi,
                   const std::vector<static_matrix<scalar_type, DIM, DIM>> cdphi,
                   const scalar_type&                                      weight) const
   {
      const size_t grad_basis_size = gphi.size();
      const size_t cell_basis_size = cdphi.size();
      const size_t grad_degree     = m_bqd.grad_degree();
      const size_t poly_space =
        grad_basis_size; // DIM * DIM * binomial(grad_degree-1 + DIM, grad_degree-1);
      const size_t DIM2 = DIM * DIM;

      //             // poly classique
      //
      size_t row = 0;
      for (size_t i = 0; i < poly_space; i += DIM2) {
         for (size_t k = 0; k < DIM; k++) {    // depend de l'ordre des bases
            for (size_t l = 0; l < DIM; l++) { // depend de l'ordre des bases
               for (size_t j = l; j < cell_basis_size; j += DIM) {
                  BG(row, j) += weight * gphi[row](l, k) * cdphi[j](l, k);
               }
               row++;
            }
         }
      }

      // RT
      //             for(std::size_t i = 0; i < grad_basis_size; i ++) {
      //                for(std::size_t j = 0; j < cell_basis_size; j ++) {
      //                   BG(i,j) = mm_prod(gphi[i], cdphi[j]);
      //                }
      //             }
   }

 public:
   matrix_type oper;
   matrix_type data;

   gradient_reconstruction_vector_full_bq(const BQData& bqd)
     : m_bqd(bqd)
   {}

   void compute(const mesh_type& msh, const cell_type& cl, const bool compute_data = true)
   {
      const size_t cell_degree     = m_bqd.cell_degree();
      const size_t face_degree     = m_bqd.face_degree();
      const size_t cell_basis_size = (m_bqd.cell_basis.range(0, cell_degree)).size();
      const size_t face_basis_size = (m_bqd.face_basis.range(0, face_degree)).size();
      const size_t grad_basis_size = m_bqd.grad_basis.size();

      const auto   fcs       = faces(msh, cl);
      const size_t num_faces = fcs.size();

      timecounter tc;
      double      t_base(0.0);
      double      t_cons(0.0);
      double      t_inv(0.0);

      const dofspace_ranges dsr(cell_basis_size, face_basis_size, num_faces);

      assert(dsr.total_size() == (cell_basis_size + num_faces * face_basis_size));

      matrix_type BG = matrix_type::Zero(grad_basis_size, dsr.total_size());

      matrix_type MG = matrix_type::Zero(grad_basis_size, grad_basis_size);

      const auto grad_quadpoints = m_bqd.grad_quadrature.integrate(msh, cl);
      for (auto& qp : grad_quadpoints) {
         tc.tic();
         const auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());
         assert(grad_basis_size == gphi.size());

         tc.toc();
         t_base += tc.to_double();

         tc.tic();

         compute_MG(MG, gphi, qp.weight());
      }

      const auto grad_cell_quadpoints = m_bqd.grad_cell_quadrature.integrate(msh, cl);
      for (auto& qp : grad_cell_quadpoints) {
         const auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());
         assert(grad_basis_size == gphi.size());

         const auto dphi = m_bqd.cell_basis.eval_gradients(msh, cl, qp.point());
         assert(cell_basis_size == dphi.size());

         compute_BG(BG, gphi, dphi, qp.weight());

         tc.toc();
         t_cons += tc.to_double();
      }

      tc.tic();
      // lower part MG
      for (size_t i = 0; i < grad_basis_size; i++)
         for (size_t j = i; j < grad_basis_size; j++)
            MG(i, j) = MG(j, i);

      tc.toc();
      t_cons += tc.to_double();

      for (size_t face_i = 0; face_i < num_faces; face_i++) {
         const auto current_face_range   = dsr.face_range(face_i);
         const auto fc                   = fcs[face_i];
         const auto n                    = normal(msh, cl, fc);
         const auto grad_face_quadpoints = m_bqd.grad_face_quadrature.integrate(msh, fc);

         for (auto& qp : grad_face_quadpoints) {
            tc.tic();
            const auto c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point());
            const auto gphi  = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());
            tc.toc();
            t_base += tc.to_double();

            // tau.n
            tc.tic();
            const auto gphi_n = compute_gphi_n(gphi, n);

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
            tc.toc();
            t_cons += tc.to_double();

            tc.tic();
            const auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());
            tc.toc();
            t_base += tc.to_double();

            tc.tic();
            matrix_type F = matrix_type::Zero(BG.rows(), current_face_range.size());

            for (size_t j = 0; j < current_face_range.size(); j++) {
               for (size_t i = 0; i < grad_basis_size; i++) {
                  F(i, j) = qp.weight() * mm_prod(gphi_n[i], f_phi[j]);
               }
            }

            assert(F.rows() == grad_basis_size);
            assert(F.cols() == current_face_range.size());

            BG.block(0, current_face_range.min(), grad_basis_size, current_face_range.size()) += F;
            tc.toc();
            t_cons += tc.to_double();
         }
      }

      tc.tic();
      oper = MG.llt().solve(BG); // GT
      tc.toc();
      t_inv += tc.to_double();
      if (compute_data) data = BG.transpose() * oper; // A

      assert(oper.rows() == grad_basis_size);
      assert(oper.cols() == dsr.total_size());
   }
};

template<typename BQData>
class projector_vector_bq
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef typename mesh_type::face        face_type;

   typedef typename BQData::cell_basis_type cell_basis_type;
   typedef typename BQData::cell_quad_type  cell_quad_type;

   typedef dynamic_matrix<scalar_type> matrix_type;
   typedef dynamic_vector<scalar_type> vector_type;

   const BQData& m_bqd;

 public:
   projector_vector_bq(const BQData& bqd)
     : m_bqd(bqd)
   {}

   matrix_type cell_mm;
   matrix_type face_mm;
   matrix_type whole_mm;
   matrix_type grad_mm;
   matrix_type pot_mm;

   template<typename Function>
   vector_type compute_whole(const mesh_type& msh, const cell_type& cl, const Function& f)
   {
      const size_t cell_degree     = m_bqd.cell_degree();
      const size_t face_degree     = m_bqd.face_degree();
      const size_t cell_basis_size = (m_bqd.cell_basis.range(0, cell_degree)).size();
      const size_t face_basis_size = m_bqd.face_basis.range(0, face_degree).size();
      const auto   fcs             = faces(msh, cl);

      this->compute_cell(msh, cl, f);

      vector_type ret = vector_type::Zero(cell_basis_size + fcs.size() * face_basis_size);
      whole_mm        = matrix_type::Zero(cell_basis_size + fcs.size() * face_basis_size,
                                   cell_basis_size + fcs.size() * face_basis_size);

      ret.block(0, 0, cell_basis_size, 1)                    = compute_cell(msh, cl, f);
      whole_mm.block(0, 0, cell_basis_size, cell_basis_size) = cell_mm;

      size_t face_offset = cell_basis_size;
      for (auto& fc : fcs) {
         matrix_type mm  = matrix_type::Zero(face_basis_size, face_basis_size);
         vector_type rhs = vector_type::Zero(face_basis_size);

         auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);
         for (auto& qp : face_quadpoints) {
            auto phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());
            assert(phi.size() == face_basis_size);

            for (size_t i = 0; i < face_basis_size; i++)
               for (size_t j = 0; j < face_basis_size; j++)
                  mm(i, j) += qp.weight() * mm_prod(phi[i], phi[j]);

            for (size_t i = 0; i < face_basis_size; i++)
               rhs(i) += qp.weight() * mm_prod(f(qp.point()), phi[i]);
         }

         ret.block(face_offset, 0, face_basis_size, 1) = mm.llt().solve(rhs);
         whole_mm.block(face_offset, face_offset, face_basis_size, face_basis_size) = mm;
         face_offset += face_basis_size;
      }

      return ret;
   }

   template<typename Function>
   vector_type compute_cell_grad(const mesh_type& msh, const cell_type& cl, const Function& f)
   {
      const size_t grad_basis_size = m_bqd.grad_basis.size();

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
      return mm.llt().solve(rhs);
   }

   template<typename Function>
   vector_type compute_pot(const mesh_type& msh, const cell_type& cl, const Function& f)
   {
      const size_t DIM             = msh.dimension;
      const size_t cell_degree     = m_bqd.cell_degree();
      const size_t cell_basis_size = DIM * binomial(cell_degree + 1 + DIM, cell_degree + 1);

      cell_basis_type cell_basis      = cell_basis_type(cell_degree + 1);
      cell_quad_type  cell_quadrature = cell_quad_type(2 * (cell_degree + 1));

      assert(cell_basis.size() == cell_basis_size);

      matrix_type mm  = matrix_type::Zero(cell_basis_size, cell_basis_size);
      vector_type rhs = vector_type::Zero(cell_basis_size);

      const auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
      for (auto& qp : cell_quadpoints) {

         const auto dphi = cell_basis.eval_gradients(msh, cl, qp.point());
         assert(cell_basis_size == dphi.size());

         for (size_t j = 0; j < cell_basis_size; j++) {
            for (size_t i = j; i < cell_basis_size; i++) {
               mm(i, j) += qp.weight() * mm_prod(dphi[i], dphi[j]);
            }
         }

         for (size_t i = 0; i < cell_basis_size; i++) {
            rhs(i) += qp.weight() * mm_prod(f(qp.point()), dphi[i]);
         }
      }

      // lower part
      for (size_t i = 0; i < cell_basis_size; i++)
         for (size_t j = i; j < cell_basis_size; j++)
            mm(i, j) = mm(j, i);

      /* LHS: take basis functions derivatives from degree 1 to K+1 */
      const auto MG_rowcol_range = cell_basis.range(1, cell_degree + 1);
      pot_mm                     = take(mm, MG_rowcol_range, MG_rowcol_range);

      // std::cout << "mm " << mm << '\n';
      // std::cout << "r " << rhs<< '\n';
      return pot_mm.ldlt().solve(rhs.tail(pot_mm.cols()));
   }
};

template<typename BQData>
class displacement_reconstruction_elas_bq
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef typename BQData::cell_basis_type cell_basis_type;
   typedef typename BQData::cell_quad_type  cell_quad_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   const BQData& m_bqd;

   cell_basis_type cell_basis;
   cell_quad_type  cell_quadrature;

 public:
   matrix_type oper;
   matrix_type data;

   displacement_reconstruction_elas_bq(const BQData& bqd)
     : m_bqd(bqd)
   {
      cell_basis      = cell_basis_type(bqd.cell_degree() + 1);
      cell_quadrature = cell_quad_type(2 * (bqd.cell_degree() + 1));
   }

   void compute(const mesh_type& msh, const cell_type& cl)
   {
      const size_t cell_degree     = m_bqd.cell_degree();
      const size_t face_degree     = m_bqd.face_degree();
      const size_t cell_basis_size = (cell_basis.range(0, cell_degree + 1)).size();
      const size_t face_basis_size = (m_bqd.face_basis.range(0, face_degree)).size();

      matrix_type stiff_mat = matrix_type::Zero(cell_basis_size, cell_basis_size);

      const auto cell_quadpoints = cell_quadrature.integrate(msh, cl);

      for (auto& qp : cell_quadpoints) {
         const auto dphi = cell_basis.eval_gradients(msh, cl, qp.point());
         assert(cell_basis_size == dphi.size());

         for (size_t j = 0; j < cell_basis_size; j++) {
            for (size_t i = j; i < cell_basis_size; i++) {
               stiff_mat(i, j) += qp.weight() * mm_prod(dphi[i], dphi[j]);
            }
         }
      }

      // lower part
      for (size_t i = 0; i < cell_basis_size; i++)
         for (size_t j = i; j < cell_basis_size; j++)
            stiff_mat(i, j) = stiff_mat(j, i);

      /* LHS: take basis functions derivatives from degree 1 to K+1 */
      const auto        MG_rowcol_range = cell_basis.range(1, cell_degree + 1);
      const matrix_type MG              = take(stiff_mat, MG_rowcol_range, MG_rowcol_range);

      /* RHS, volumetric part. */
      const auto BG_row_range = cell_basis.range(1, cell_degree + 1);
      const auto BG_col_range = cell_basis.range(0, cell_degree);

      const auto   fcs       = faces(msh, cl);
      const size_t num_faces = fcs.size();

      const size_t num_cell_dofs = BG_col_range.size();

      const dofspace_ranges dsr(num_cell_dofs, face_basis_size, num_faces);

      assert(dsr.total_size() == (num_cell_dofs + num_faces * face_basis_size));

      matrix_type BG = matrix_type::Zero(BG_row_range.size(), dsr.total_size());

      BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
        take(stiff_mat, BG_row_range, BG_col_range);

      for (size_t face_i = 0; face_i < num_faces; face_i++) {
         const auto current_face_range = dsr.face_range(face_i);
         const auto fc                 = fcs[face_i];
         const auto n                  = normal(msh, cl, fc);

         const auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);

         for (auto& qp : face_quadpoints) {
            auto       c_phi  = cell_basis.eval_functions(msh, cl, qp.point()); // 0, m_degree);
            const auto c_dphi = cell_basis.eval_gradients(msh, cl, qp.point()); // 1, m_degree+1);

            decltype(c_phi) c_dphi_n;

            c_dphi_n.reserve(BG_row_range.to() - BG_row_range.from());

            for (size_t i = BG_row_range.from(); i < BG_row_range.to(); i++) {
               c_dphi_n.push_back(mm_prod(c_dphi[i], n));
            }

            matrix_type T = matrix_type::Zero(BG.rows(), BG_col_range.size());

            assert(c_dphi_n.size() == BG.rows());

            for (size_t j = 0; j < BG_col_range.size(); j++) {
               for (size_t i = 0; i < BG.rows(); i++) {
                  T(i, j) = qp.weight() * mm_prod(c_dphi_n[i], c_phi[j]);
               }
            }

            BG.block(0, 0, BG.rows(), BG_col_range.size()) -= T;

            const auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());

            assert(f_phi.size() == face_basis_size);

            matrix_type F = matrix_type::Zero(BG.rows(), current_face_range.size());

            for (size_t j = 0; j < current_face_range.size(); j++) {
               for (size_t i = 0; i < BG.rows(); i++) {
                  F(i, j) = qp.weight() * mm_prod(c_dphi_n[i], f_phi[j]);
               }
            }

            BG.block(0, current_face_range.min(), BG.rows(), current_face_range.size()) += F;
         }
      }

      assert(MG.rows() == MG.cols());
      assert(MG.cols() == BG.rows());

      oper = MG.ldlt().solve(BG);   // GT
      data = BG.transpose() * oper; // A
   }
};

template<typename BQData>
class assembler_vector_bq
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef typename mesh_type::face        face_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   typedef Eigen::Triplet<scalar_type> triplet_type;

   std::vector<triplet_type> m_triplets;
   size_t                    m_num_unknowns;

   const BQData& m_bqd;

 public:
   typedef Eigen::SparseMatrix<scalar_type> sparse_matrix_type;
   typedef dynamic_vector<scalar_type>      vector_type;

   sparse_matrix_type matrix;
   vector_type        rhs;

   assembler_vector_bq() = delete;

   assembler_vector_bq(const mesh_type& msh, const BQData& bqd)
     : m_bqd(bqd)
   {
      m_num_unknowns = m_bqd.face_basis.size() * (msh.faces_size() + msh.boundary_faces_size());
      matrix         = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
      rhs            = vector_type::Zero(m_num_unknowns);
   }

   assembler_vector_bq(const mesh_type&          msh,
                       const BQData&             bqd,
                       const BoundaryConditions& boundary_conditions)
     : m_bqd(bqd)
   {
      const size_t faces_dofs = m_bqd.face_basis.size() * msh.faces_size();
      const size_t lag_dofs =
        boundary_conditions.nb_lags() * m_bqd.face_basis.size() / msh.dimension;
      m_num_unknowns = faces_dofs + lag_dofs;
      matrix         = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
      rhs            = vector_type::Zero(m_num_unknowns);
   }

   template<typename LocalContrib>
   void assemble(const mesh_type& msh, const cell_type& cl, const LocalContrib& lc)
   {
      const size_t face_basis_size = m_bqd.face_basis.size();

      const auto          fcs = faces(msh, cl);
      std::vector<size_t> l2g(fcs.size() * face_basis_size);

      assert(face_basis_size == face_basis_size);

      for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
         const auto fc  = fcs[face_i];
         const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const auto face_id = eid.second;

         const auto face_offset = face_id * face_basis_size;

         const auto pos = face_i * face_basis_size;

         for (size_t i = 0; i < face_basis_size; i++)
            l2g[pos + i] = face_offset + i;
      }

      assert(lc.first.rows() == fcs.size() * face_basis_size);
      assert(lc.first.rows() == lc.first.cols());
      assert(lc.first.rows() == lc.second.size());
      assert(lc.second.size() == l2g.size());

      // std::cout << lc.second.size() << " " << l2g.size() << std::endl;

#ifdef FILL_COLMAJOR
      for (size_t j = 0; j < lc.first.cols(); j++) {
         for (size_t i = 0; i < lc.first.rows(); i++)
            m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));

         rhs(l2g.at(j)) += lc.second(j);
      }
#else
      for (size_t i = 0; i < lc.first.rows(); i++) {
         for (size_t j = 0; j < lc.first.cols(); j++)
            m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));

         rhs(l2g.at(i)) += lc.second(i);
      }
#endif
   }

   template<typename Function>
   void impose_boundary_conditions(const mesh_type& msh, const Function& bc)
   {
      const size_t fbs    = m_bqd.face_basis.size();
      size_t       face_i = 0;
      for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++) {
         auto bfc = *itor;

         const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const size_t face_id = eid.second;

         const size_t face_offset          = face_id * fbs;
         const size_t face_offset_lagrange = (msh.faces_size() + face_i) * fbs;

         const auto fqd = m_bqd.face_quadrature.integrate(msh, bfc);

         matrix_type MFF   = matrix_type::Zero(fbs, fbs);
         vector_type rhs_f = vector_type::Zero(fbs);

         for (auto& qp : fqd) {
            const auto f_phi = m_bqd.face_basis.eval_functions(msh, bfc, qp.point());

            assert(f_phi.size() == fbs);

            for (size_t i = 0; i < fbs; i++)
               for (size_t j = i; j < fbs; j++)
                  MFF(i, j) += qp.weight() * mm_prod(f_phi[i], f_phi[j]);

            // lower part
            for (size_t i = 1; i < fbs; i++)
               for (size_t j = 0; j < i; j++)
                  MFF(i, j) = MFF(j, i);

            for (size_t i = 0; i < fbs; i++)
               rhs_f(i) += qp.weight() * mm_prod(f_phi[i], bc(qp.point()));
         }

#ifdef FILL_COLMAJOR
         for (size_t j = 0; j < fbs; j++) {
            for (size_t i = 0; i < fbs; i++) {
               m_triplets.push_back(
                 triplet_type(face_offset + i, face_offset_lagrange + j, MFF(i, j)));
               m_triplets.push_back(
                 triplet_type(face_offset_lagrange + j, face_offset + i, MFF(i, j)));
            }
            rhs(face_offset_lagrange + j) = rhs_f(j);
         }
#else
         for (size_t i = 0; i < fbs; i++) {
            for (size_t j = 0; j < fbs; j++) {
               m_triplets.push_back(
                 triplet_type(face_offset + i, face_offset_lagrange + j, MFF(i, j)));
               m_triplets.push_back(
                 triplet_type(face_offset_lagrange + j, face_offset + i, MFF(i, j)));
            }
            rhs(face_offset_lagrange + i) = rhs_f(i);
         }
#endif

         face_i++;
      }
   }

   template<typename Function, typename NeumannFunction>
   void impose_boundary_conditions_nl(const mesh_type&                msh,
                                      const Function&                 bc,
                                      const NeumannFunction&          g,
                                      const std::vector<vector_type>& sol_faces,
                                      const std::vector<vector_type>& sol_lagr,
                                      const BoundaryConditions&       boundary_conditions)
   {
      const size_t DIM             = msh.dimension;
      const size_t face_basis_size = m_bqd.face_basis.size();
      const size_t lagr_size       = face_basis_size / DIM;
      const size_t faces_offset    = msh.faces_size() * face_basis_size;

      size_t face_i = 0;
      size_t face_dir(0);

      for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++) {
         auto bfc = *itor;

         const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const auto   face_id     = eid.second;
         const size_t face_offset = face_id * face_basis_size;

         // std::cout << "face "  << face_i << " "<<
         // boundary_conditions.is_boundary_dirichlet(face_i)   << '\n';
         if (boundary_conditions.is_boundary_dirichlet(face_i)) {
            const size_t nb_lag          = boundary_conditions.nb_lag_conditions_faceI(face_dir);
            const size_t lagr_basis_size = nb_lag * lagr_size;
            const size_t lag_pos         = boundary_conditions.begin_lag_conditions_faceI(face_dir);
            const size_t face_offset_lagrange = faces_offset + lag_pos * lagr_size;

            const auto  fqd    = m_bqd.face_quadrature.integrate(msh, bfc);
            matrix_type MFF    = matrix_type::Zero(face_basis_size, face_basis_size);
            vector_type rhs_bc = vector_type::Zero(face_basis_size);
            vector_type rhs_l  = vector_type::Zero(face_basis_size);
            for (auto& qp : fqd) {
               const auto f_phi = m_bqd.face_basis.eval_functions(msh, bfc, qp.point());
               assert(f_phi.size() == face_basis_size);
               for (size_t i = 0; i < face_basis_size; i++)
                  for (size_t j = i; j < face_basis_size; j++)
                     MFF(i, j) += qp.weight() * mm_prod(f_phi[i], f_phi[j]);

               for (size_t i = 0; i < face_basis_size; i++)
                  rhs_bc(i) += qp.weight() * mm_prod(f_phi[i], bc(qp.point()));
            }

            // lower part
            for (size_t i = 1; i < face_basis_size; i++)
               for (size_t j = 0; j < i; j++)
                  MFF(i, j) = MFF(j, i);
            // impose displacement
            const vector_type rhs_fc = rhs_bc - MFF * sol_faces.at(face_id);
            const vector_type rhs_f  = -MFF * sol_faces.at(face_id);

            vector_type rhs_f2;
            matrix_type MFF2;

            bool dirichlet_standart = true;
            switch (boundary_conditions.boundary_type(face_i)) {
               case CLAMPED: {
                  dirichlet_standart = false;
                  rhs_f2             = rhs_f;
                  MFF2               = MFF;
                  rhs_l              = MFF * sol_lagr.at(face_dir);
                  break;
               }
               case DX: {
                  dirichlet_standart = false;
                  rhs_f2.resize(lagr_basis_size, 1);
                  MFF2.resize(face_basis_size, lagr_basis_size);

                  assert(rhs_f2.rows() == lagr_basis_size);
                  assert(MFF2.rows() == face_basis_size);
                  assert(MFF2.cols() == lagr_basis_size);

                  size_t ind(0);
                  for (size_t i = 0; i < face_basis_size; i += DIM) {
                     MFF2.col(ind) = MFF.col(i);
                     rhs_f2(ind)   = rhs_f(i);
                     ind++;
                  }

                  assert(ind == lagr_basis_size);

                  vector_type rhs_l2 = MFF2 * sol_lagr.at(face_dir);
                  ind                = 0;
                  for (size_t i = 0; i < face_basis_size; i += DIM) {
                     rhs_l(i) = rhs_l2(ind);
                     ind++;
                  }

                  assert(ind == lagr_basis_size);

                  break;
               }
               //
               //                   case DY:
               //                      nb_lag_conditions -= DIM - 1;
               //                      break;
               //                   case DZ:
               //                      if( DIM != 3){
               //                         std::cout << "Invalid condition for face:" << face_id <<
               //                         std::endl; throw std::invalid_argument(" ONLY DIM = 3 for
               //                         this Dirichlet Conditions");
               //                      }
               //                      else
               //                         nb_lag_conditions -= 2;
               //                      break;
               //                   case DXDY:
               //                      if( DIM != 3){
               //                         std::cout << "Invalid condition for face:" << face_id <<
               //                         std::endl; throw std::invalid_argument(" ONLY DIM = 3 for
               //                         this Dirichlet Conditions");
               //                      }
               //                      else
               //                         nb_lag_conditions -= 1;
               //                      break;
               //                   case DXDZ:
               //                      if( DIM != 3){
               //                         std::cout << "Invalid condition for face:" << face_id <<
               //                         std::endl; throw std::invalid_argument(" ONLY DIM = 3 for
               //                         this Dirichlet Conditions");
               //                      }
               //                      else
               //                         nb_lag_conditions -= 1;
               //                      break;
               //                   case DYDZ:
               //                      if( DIM != 3){
               //                         std::cout << "Invalid condition for face:" << face_id <<
               //                         std::endl; throw std::invalid_argument(" ONLY DIM = 3 for
               //                         this Dirichlet Conditions");
               //                      }
               //                      else
               //                         nb_lag_conditions -= 1;
               //                      break;
               case OTHER: break;
               default:
                  std::cout << "Unknown Dirichlet Conditions: we do nothing" << std::endl;
                  break;
            }

            if (dirichlet_standart) {
               rhs_f2 = rhs_fc;
               MFF2   = MFF;
               rhs_l  = MFF * sol_lagr.at(face_dir);
            }

            assert(MFF2.rows() == face_basis_size);
#ifdef FILL_COLMAJOR
            for (size_t j = 0; j < MFF2.cols(); j++) {
               for (size_t i = 0; i < MFF2.rows(); i++) {
                  m_triplets.push_back(
                    triplet_type(face_offset + i, face_offset_lagrange + j, MFF2(i, j)));
                  m_triplets.push_back(
                    triplet_type(face_offset_lagrange + j, face_offset + i, MFF2(i, j)));
               }
               rhs(face_offset_lagrange + j) = rhs_f2(j);
            }

            for (size_t i = 0; i < face_basis_size; i++)
               rhs(face_offset + i) -= rhs_l(i);

#else
            for (size_t i = 0; i < MFF2.rows(); i++) {
               for (size_t j = 0; j < MFF2.cols(); j++) {
                  m_triplets.push_back(
                    triplet_type(face_offset + i, face_offset_lagrange + j, MFF2(i, j)));
                  m_triplets.push_back(
                    triplet_type(face_offset_lagrange + j, face_offset + i, MFF2(i, j)));
               }
               rhs(face_offset_lagrange + i) = rhs_f2(i);
            }
#endif
            face_dir++;
         } else {
            if (boundary_conditions.boundary_type(face_i) == NEUMANN) {
               vector_type TF = vector_type::Zero(face_basis_size);

               const auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, bfc);
               for (auto& qp : face_quadpoints) {
                  auto fphi = m_bqd.face_basis.eval_functions(msh, bfc, qp.point());
                  assert(fphi.size() == face_basis_size);

                  for (size_t i = 0; i < face_basis_size; i++)
                     TF(i) += qp.weight() * mm_prod(g(qp.point()), fphi[i]);
               }

               for (size_t i = 0; i < TF.rows(); i++) {
                  rhs(face_offset + i) += TF(i);
               }
            }
         }
         face_i++;
      }
   }

   void finalize()
   {
      matrix.setFromTriplets(m_triplets.begin(), m_triplets.end());
      m_triplets.clear();
   }

   void finalize(sparse_matrix_type& mat, vector_type& vec)
   {
      mat = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
      mat.setFromTriplets(m_triplets.begin(), m_triplets.end());
      m_triplets.clear();
      vec = rhs;
   }
};

template<typename BQData>
class assembler_full_vector_bq
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef dynamic_matrix<scalar_type> matrix_type;

   typedef Eigen::Triplet<scalar_type> triplet_type;

   std::vector<triplet_type> m_triplets;
   size_t                    m_num_unknowns;

   const BQData& m_bqd;

 public:
   typedef Eigen::SparseMatrix<scalar_type> sparse_matrix_type;
   typedef dynamic_vector<scalar_type>      vector_type;

   sparse_matrix_type matrix;
   vector_type        rhs;

   assembler_full_vector_bq() = delete;

   assembler_full_vector_bq(const mesh_type& msh, const BQData& bqd)
     : m_bqd(bqd)
   {
      const size_t cells_dofs = m_bqd.cell_basis.size() * msh.cells_size();
      const size_t faces_dofs =
        m_bqd.face_basis.size() * (msh.faces_size() + msh.boundary_faces_size());
      m_num_unknowns = cells_dofs + faces_dofs;
      matrix         = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
      rhs            = vector_type::Zero(m_num_unknowns);
   }

   assembler_full_vector_bq(const mesh_type&          msh,
                            const BQData&             bqd,
                            const BoundaryConditions& boundary_conditions)
     : m_bqd(bqd)
   {
      const size_t cells_dofs = m_bqd.cell_basis.size() * msh.cells_size();
      const size_t faces_dofs = m_bqd.face_basis.size() * msh.faces_size();
      const size_t lag_dofs =
        boundary_conditions.nb_lags() * m_bqd.face_basis.size() / msh.dimension;
      m_num_unknowns = cells_dofs + faces_dofs + lag_dofs;
      matrix         = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
      rhs            = vector_type::Zero(m_num_unknowns);
   }

   template<typename LocalContrib>
   void assemble(const mesh_type& msh, const cell_type& cl, const LocalContrib& lc)
   {
      // assemble cell Parameters
      const size_t cell_basis_size = m_bqd.cell_basis.size();
      const size_t face_basis_size = m_bqd.face_basis.size();
      const auto   fcs             = faces(msh, cl);
      const size_t total_dof       = cell_basis_size + fcs.size() * face_basis_size;

      std::vector<size_t> l2g(total_dof);

      const auto cid = find_element_id(msh.cells_begin(), msh.cells_end(), cl);
      if (!cid.first) throw std::invalid_argument("This is a bug: cell not found");

      const auto cell_id     = cid.second;
      const auto cell_offset = cell_id * cell_basis_size;

      for (size_t i = 0; i < cell_basis_size; i++) {
         l2g[i] = cell_offset + i;
      }

      const size_t cells_offset = msh.cells_size() * cell_basis_size;
      for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
         const auto fc  = fcs[face_i];
         const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const auto face_id     = eid.second;
         const auto face_offset = cells_offset + face_id * face_basis_size;
         const auto pos         = cell_basis_size + face_i * face_basis_size;

         for (size_t i = 0; i < face_basis_size; i++)
            l2g[pos + i] = face_offset + i;
      }

      assert(lc.first.rows() == lc.first.cols());
      assert(lc.first.rows() == lc.second.size());
      assert(lc.second.size() == l2g.size());

      // std::cout << lc.second.size() << " " << l2g.size() << std::endl;

#ifdef FILL_COLMAJOR
      for (size_t j = 0; j < lc.first.cols(); j++) {
         for (size_t i = 0; i < lc.first.rows(); i++)
            m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));

         rhs(l2g.at(j)) += lc.second(j);
      }
#else
      for (size_t i = 0; i < lc.first.rows(); i++) {
         for (size_t j = 0; j < lc.first.cols(); j++)
            m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));

         rhs(l2g.at(i)) += lc.second(i);
      }
#endif
   }

   template<typename Function>
   void impose_boundary_conditions(const mesh_type& msh, const Function& bc)
   {
      const size_t face_basis_size = m_bqd.face_basis.size();
      const size_t cell_basis_size = m_bqd.cell_basis.size();
      const size_t cells_offset    = msh.cells_size() * cell_basis_size;

      size_t face_i = 0;
      for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++) {
         auto bfc = *itor;

         const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const size_t face_id = eid.second;

         const size_t face_offset = cells_offset + face_id * face_basis_size;
         const size_t face_offset_lagrange =
           cells_offset + (msh.faces_size() + face_i) * face_basis_size;

         const auto fqd = m_bqd.face_quadrature.integrate(msh, bfc);

         matrix_type MFF   = matrix_type::Zero(face_basis_size, face_basis_size);
         vector_type rhs_f = vector_type::Zero(face_basis_size);

         for (auto& qp : fqd) {
            const auto f_phi = m_bqd.face_basis.eval_functions(msh, bfc, qp.point());

            assert(f_phi.size() == face_basis_size);

            for (size_t i = 0; i < face_basis_size; i++) {
               for (size_t j = i; j < face_basis_size; j++) {
                  MFF(i, j) += qp.weight() * mm_prod(f_phi[i], f_phi[j]);
               }
            }

            // lower part
            for (size_t i = 1; i < face_basis_size; i++) {
               for (size_t j = 0; j < i; j++) {
                  MFF(i, j) = MFF(j, i);
               }
            }

            for (size_t i = 0; i < face_basis_size; i++) {
               rhs_f(i) += qp.weight() * mm_prod(f_phi[i], bc(qp.point()));
            }
         }

#ifdef FILL_COLMAJOR
         for (size_t j = 0; j < face_basis_size; j++) {
            for (size_t i = 0; i < face_basis_size; i++) {
               m_triplets.push_back(
                 triplet_type(face_offset + i, face_offset_lagrange + j, MFF(i, j)));
               m_triplets.push_back(
                 triplet_type(face_offset_lagrange + j, face_offset + i, MFF(i, j)));
            }
            rhs(face_offset_lagrange + j) = rhs_f(j);
         }
#else
         for (size_t i = 0; i < face_basis_size; i++) {
            for (size_t j = 0; j < face_basis_size; j++) {
               m_triplets.push_back(
                 triplet_type(face_offset + i, face_offset_lagrange + j, MFF(i, j)));
               m_triplets.push_back(
                 triplet_type(face_offset_lagrange + j, face_offset + i, MFF(i, j)));
            }
            rhs(face_offset_lagrange + i) = rhs_f(i);
         }
#endif

         face_i++;
      }
   }

   // non-linear case
   template<typename Function, typename NeumannFunction>
   void impose_boundary_conditions_nl(const mesh_type&                msh,
                                      const Function&                 bc,
                                      const NeumannFunction&          g,
                                      const std::vector<vector_type>& sol_faces,
                                      const std::vector<vector_type>& sol_lagr,
                                      const BoundaryConditions&       boundary_conditions)
   {
      const size_t DIM             = msh.dimension;
      const size_t face_basis_size = m_bqd.face_basis.size();
      const size_t lagr_size       = face_basis_size / DIM;
      size_t       face_i          = 0;
      size_t       face_dir(0);
      const size_t cells_offset = msh.cells_size() * m_bqd.cell_basis.size();

      for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++) {
         auto bfc = *itor;

         const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const auto   face_id     = eid.second;
         const size_t face_offset = cells_offset + face_id * face_basis_size;

         // std::cout << "face "  << face_i << " "<<
         // boundary_conditions.is_boundary_dirichlet(face_i)   << '\n';
         if (boundary_conditions.is_boundary_dirichlet(face_i)) {
            const size_t nb_lag          = boundary_conditions.nb_lag_conditions_faceI(face_dir);
            const size_t lagr_basis_size = nb_lag * lagr_size;
            const size_t lag_pos         = boundary_conditions.begin_lag_conditions_faceI(face_dir);
            const size_t face_offset_lagrange =
              cells_offset + msh.faces_size() * face_basis_size + lag_pos * lagr_size;

            const auto fqd = m_bqd.face_quadrature.integrate(msh, bfc);

            matrix_type MFF    = matrix_type::Zero(face_basis_size, face_basis_size);
            vector_type rhs_bc = vector_type::Zero(face_basis_size);
            vector_type rhs_l  = vector_type::Zero(face_basis_size);
            for (auto& qp : fqd) {
               const auto f_phi = m_bqd.face_basis.eval_functions(msh, bfc, qp.point());
               assert(f_phi.size() == face_basis_size);
               for (size_t i = 0; i < face_basis_size; i++)
                  for (size_t j = i; j < face_basis_size; j++)
                     MFF(i, j) += qp.weight() * mm_prod(f_phi[i], f_phi[j]);

               for (size_t i = 0; i < face_basis_size; i++)
                  rhs_bc(i) += qp.weight() * mm_prod(f_phi[i], bc(qp.point()));
            }

            // lower part
            for (size_t i = 1; i < face_basis_size; i++)
               for (size_t j = 0; j < i; j++)
                  MFF(i, j) = MFF(j, i);

            // impose displacement
            const vector_type rhs_fc = rhs_bc - MFF * sol_faces.at(face_id);
            const vector_type rhs_f  = -MFF * sol_faces.at(face_id);

            vector_type rhs_f2;
            matrix_type MFF2;

            bool dirichlet_standart = true;
            switch (boundary_conditions.boundary_type(face_i)) {
               case CLAMPED: {
                  dirichlet_standart = false;
                  rhs_f2             = rhs_f;
                  MFF2               = MFF;
                  rhs_l              = MFF * sol_lagr.at(face_dir);
                  break;
               }
               case DX: {
                  dirichlet_standart = false;
                  rhs_f2.resize(lagr_basis_size, 1);
                  MFF2.resize(face_basis_size, lagr_basis_size);

                  assert(rhs_f2.rows() == lagr_basis_size);
                  assert(MFF2.rows() == face_basis_size);
                  assert(MFF2.cols() == lagr_basis_size);

                  size_t ind(0);
                  for (size_t i = 0; i < face_basis_size; i += DIM) {
                     MFF2.col(ind) = MFF.col(i);
                     rhs_f2(ind)   = rhs_f(i);
                     ind++;
                  }

                  assert(ind == lagr_basis_size);

                  vector_type rhs_l2 = MFF2 * sol_lagr.at(face_dir);
                  ind                = 0;
                  for (size_t i = 0; i < face_basis_size; i += DIM) {
                     rhs_l(i) = rhs_l2(ind);
                     ind++;
                  }

                  assert(ind == lagr_basis_size);

                  break;
               }
               //
               //                   case DY:
               //                      nb_lag_conditions -= DIM - 1;
               //                      break;
               //                   case DZ:
               //                      if( DIM != 3){
               //                         std::cout << "Invalid condition for face:" << face_id <<
               //                         std::endl; throw std::invalid_argument(" ONLY DIM = 3 for
               //                         this Dirichlet Conditions");
               //                      }
               //                      else
               //                         nb_lag_conditions -= 2;
               //                      break;
               //                   case DXDY:
               //                      if( DIM != 3){
               //                         std::cout << "Invalid condition for face:" << face_id <<
               //                         std::endl; throw std::invalid_argument(" ONLY DIM = 3 for
               //                         this Dirichlet Conditions");
               //                      }
               //                      else
               //                         nb_lag_conditions -= 1;
               //                      break;
               //                   case DXDZ:
               //                      if( DIM != 3){
               //                         std::cout << "Invalid condition for face:" << face_id <<
               //                         std::endl; throw std::invalid_argument(" ONLY DIM = 3 for
               //                         this Dirichlet Conditions");
               //                      }
               //                      else
               //                         nb_lag_conditions -= 1;
               //                      break;
               //                   case DYDZ:
               //                      if( DIM != 3){
               //                         std::cout << "Invalid condition for face:" << face_id <<
               //                         std::endl; throw std::invalid_argument(" ONLY DIM = 3 for
               //                         this Dirichlet Conditions");
               //                      }
               //                      else
               //                         nb_lag_conditions -= 1;
               //                      break;
               case OTHER: break;
               default:
                  std::cout << "Unknown Dirichlet Conditions: we do nothing" << std::endl;
                  break;
            }

            if (dirichlet_standart) {
               rhs_f2 = rhs_fc;
               MFF2   = MFF;
               rhs_l  = MFF * sol_lagr.at(face_dir);
            }

            assert(MFF2.rows() == face_basis_size);

#ifdef FILL_COLMAJOR
            for (size_t j = 0; j < MFF2.cols(); j++) {
               for (size_t i = 0; i < MFF2.rows(); i++) {
                  m_triplets.push_back(
                    triplet_type(face_offset + i, face_offset_lagrange + j, MFF2(i, j)));
                  m_triplets.push_back(
                    triplet_type(face_offset_lagrange + j, face_offset + i, MFF2(i, j)));
               }
               rhs(face_offset_lagrange + j) = rhs_f2(j);
            }

            for (size_t i = 0; i < face_basis_size; i++)
               rhs(face_offset + i) -= rhs_l(i);

#else
            for (size_t i = 0; i < MFF2.rows(); i++) {
               for (size_t j = 0; j < MFF2.cols(); j++) {
                  m_triplets.push_back(
                    triplet_type(face_offset + i, face_offset_lagrange + j, MFF2(i, j)));
                  m_triplets.push_back(
                    triplet_type(face_offset_lagrange + j, face_offset + i, MFF2(i, j)));
               }
               rhs(face_offset_lagrange + i) = rhs_f2(i);
               rhs(face_offset + j) -= rhs_l(j);
            }
#endif
            face_dir++;
         } else {

            if (boundary_conditions.boundary_type(face_i) == NEUMANN) {
               vector_type TF = vector_type::Zero(face_basis_size);

               const auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, bfc);
               for (auto& qp : face_quadpoints) {
                  auto fphi = m_bqd.face_basis.eval_functions(msh, bfc, qp.point());
                  assert(fphi.size() == face_basis_size);

                  for (size_t i = 0; i < face_basis_size; i++)
                     TF(i) += qp.weight() * mm_prod(g(qp.point()), fphi[i]);
               }

               for (size_t i = 0; i < TF.rows(); i++) {
                  rhs(face_offset + i) += TF(i);
               }
            }
         }
         face_i++;
      }
   }

   void finalize()
   {
      matrix.setFromTriplets(m_triplets.begin(), m_triplets.end());
      m_triplets.clear();
   }

   void finalize(sparse_matrix_type& mat, vector_type& vec)
   {
      mat = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
      mat.setFromTriplets(m_triplets.begin(), m_triplets.end());
      m_triplets.clear();
      vec = rhs;
   }
};

} // namespace hho

} // namespace disk
