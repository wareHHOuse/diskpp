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
#include "contrib/sol2/sol.hpp"
#include "hho/gradient_reconstruction.hpp"
#include "hho/hho_bq.hpp"
#include "timecounter.h"

//#define USE_BLAS
#define FILL_COLMAJOR

namespace disk {

template<typename Mesh,
         typename CellBasisType,
         typename CellQuadType,
         typename FaceBasisType,
         typename FaceQuadType>
class projector
{
   typedef Mesh                            mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef typename mesh_type::face        face_type;

   typedef CellBasisType cell_basis_type;
   typedef CellQuadType  cell_quadrature_type;
   typedef FaceBasisType face_basis_type;
   typedef FaceQuadType  face_quadrature_type;

   typedef dynamic_matrix<scalar_type> matrix_type;
   typedef dynamic_vector<scalar_type> vector_type;

   cell_basis_type      cell_basis;
   cell_quadrature_type cell_quadrature;

   face_basis_type      face_basis;
   face_quadrature_type face_quadrature;

   size_t m_degree;

 public:
   projector()
     : m_degree(1)
   {
      cell_basis      = cell_basis_type(m_degree);
      cell_quadrature = cell_quadrature_type(2 * m_degree);
      face_basis      = face_basis_type(m_degree);
      face_quadrature = face_quadrature_type(2 * m_degree);
   }

   projector(size_t degree)
     : m_degree(degree)
   {
      cell_basis      = cell_basis_type(m_degree);
      cell_quadrature = cell_quadrature_type(2 * m_degree);
      face_basis      = face_basis_type(m_degree);
      face_quadrature = face_quadrature_type(2 * m_degree);
   }

   matrix_type cell_mm;

   template<typename Function>
   vector_type compute_cell(const mesh_type& msh, const cell_type& cl, const Function& f)
   {
      matrix_type mm  = matrix_type::Zero(cell_basis.size(), cell_basis.size());
      vector_type rhs = vector_type::Zero(cell_basis.size());

      auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
      for (auto& qp : cell_quadpoints) {
         auto phi = cell_basis.eval_functions(msh, cl, qp.point());

         mm += qp.weight() * phi * phi.transpose();
         rhs += qp.weight() * f(qp.point()) * phi;
      }

      cell_mm = mm;
      return mm.llt().solve(rhs);
   }

   template<typename Function>
   vector_type compute_whole(const mesh_type& msh, const cell_type& cl, const Function& f)
   {
      auto        fcs = faces(msh, cl);
      vector_type ret(cell_basis.size() + fcs.size() * face_basis.size());

      ret.block(0, 0, cell_basis.size(), 1) = compute_cell(msh, cl, f);

      size_t face_offset = cell_basis.size();
      for (auto& fc : fcs) {
         matrix_type mm  = matrix_type::Zero(face_basis.size(), face_basis.size());
         vector_type rhs = vector_type::Zero(face_basis.size());

         auto face_quadpoints = face_quadrature.integrate(msh, fc);
         for (auto& qp : face_quadpoints) {
            auto phi = face_basis.eval_functions(msh, fc, qp.point());

            mm += qp.weight() * phi * phi.transpose();
            rhs += qp.weight() * f(qp.point()) * phi;
         }

         ret.block(face_offset, 0, face_basis.size(), 1) = mm.llt().solve(rhs);
         face_offset += face_basis.size();
      }

      return ret;
   }
};

template<typename BQData>
class projector_bq
{
   typedef typename BQData::mesh_type       mesh_type;
   typedef typename mesh_type::scalar_type  scalar_type;
   typedef typename mesh_type::cell         cell_type;
   typedef typename mesh_type::face         face_type;
   typedef typename BQData::cell_basis_type cell_basis_type;
   typedef typename BQData::face_basis_type face_basis_type;
   typedef typename BQData::cell_quad_type  cell_quadrature_type;
   typedef typename BQData::face_quad_type  face_quadrature_type;
   typedef dynamic_matrix<scalar_type>      matrix_type;
   typedef dynamic_vector<scalar_type>      vector_type;

   const BQData& m_bqd;

 public:
   projector_bq(const BQData& bqd)
     : m_bqd(bqd)
   {}

   matrix_type cell_mm;
   matrix_type grad_mm;

   template<typename Function>
   vector_type compute_cell(const mesh_type& msh, const cell_type& cl, const Function& f)
   {
      auto cell_degree     = m_bqd.cell_degree();
      auto cell_basis_size = m_bqd.cell_basis.range(0, cell_degree).size();

      matrix_type mm  = matrix_type::Zero(cell_basis_size, cell_basis_size);
      vector_type rhs = vector_type::Zero(cell_basis_size);

      auto cell_quadpoints = m_bqd.cell_quadrature.integrate(msh, cl);
      for (auto& qp : cell_quadpoints) {
         auto phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);

         mm += qp.weight() * phi * phi.transpose();
         rhs += qp.weight() * f(qp.point()) * phi;
      }

      cell_mm = mm;

      return mm.llt().solve(rhs);
   }

   template<typename Function>
   vector_type compute_face(const mesh_type& msh, const face_type& fc, const Function& f)
   {
      auto face_basis_size = howmany_dofs(m_bqd.face_basis);

      matrix_type mm  = matrix_type::Zero(face_basis_size, face_basis_size);
      vector_type rhs = vector_type::Zero(face_basis_size);

      auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);
      for (auto& qp : face_quadpoints) {
         auto phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());

         mm += qp.weight() * phi * phi.transpose();
         rhs += qp.weight() * f(qp.point()) * phi;
      }

      return mm.llt().solve(rhs);
   }

   template<typename Function>
   vector_type compute_whole(const mesh_type& msh, const cell_type& cl, const Function& f)
   {
      auto cell_degree     = m_bqd.cell_degree();
      auto face_degree     = m_bqd.face_degree();
      auto cell_basis_size = m_bqd.cell_basis.range(0, cell_degree).size();
      auto face_basis_size = m_bqd.face_basis.range(0, face_degree).size();

      auto        fcs = faces(msh, cl);
      vector_type ret(cell_basis_size + fcs.size() * face_basis_size);

      ret.block(0, 0, cell_basis_size, 1) = compute_cell(msh, cl, f);

      size_t face_offset = cell_basis_size;
      for (auto& fc : fcs) {
         matrix_type mm  = matrix_type::Zero(face_basis_size, face_basis_size);
         vector_type rhs = vector_type::Zero(face_basis_size);

         auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);
         for (auto& qp : face_quadpoints) {
            auto phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());

            mm += qp.weight() * phi * phi.transpose();
            rhs += qp.weight() * f(qp.point()) * phi;
         }

         ret.block(face_offset, 0, face_basis_size, 1) = mm.llt().solve(rhs);
         face_offset += face_basis_size;
      }

      return ret;
   }

   template<typename Function>
   vector_type compute_cell_grad(const mesh_type& msh, const cell_type& cl, const Function& f)
   {
      const size_t cell_degree = m_bqd.cell_degree();

      const size_t grad_basis_size = (m_bqd.grad_basis.range(0, cell_degree + 1)).size();

      matrix_type mm  = matrix_type::Zero(grad_basis_size, grad_basis_size);
      vector_type rhs = vector_type::Zero(grad_basis_size);

      auto grad_quadpoints = m_bqd.grad_quadrature.integrate(msh, cl);
      for (auto& qp : grad_quadpoints) {
         auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());
         assert(gphi.size() == grad_basis_size);

         for (size_t i = 0; i < grad_basis_size; i++)
            for (size_t j = i; j < grad_basis_size; j++)
               mm(i, j) += qp.weight() * mm_prod(gphi[i], gphi[j]);

         // lower part
         for (size_t i = 1; i < grad_basis_size; i++)
            for (size_t j = 0; j < i; j++)
               mm(i, j) = mm(j, i);

         for (size_t i = 0; i < grad_basis_size; i++) {
            rhs(i) += qp.weight() * mm_prod(f(qp.point()), gphi[i]);
         }
      }

      grad_mm = mm;
      return mm.llt().solve(rhs);
   }
};

template<typename BQData>
class eigval_mass_matrix_bq
{
   typedef typename BQData::mesh_type       mesh_type;
   typedef typename mesh_type::scalar_type  scalar_type;
   typedef typename mesh_type::cell         cell_type;
   typedef typename mesh_type::face         face_type;
   typedef typename BQData::cell_basis_type cell_basis_type;
   typedef typename BQData::face_basis_type face_basis_type;
   typedef typename BQData::cell_quad_type  cell_quadrature_type;
   typedef typename BQData::face_quad_type  face_quadrature_type;
   typedef dynamic_matrix<scalar_type>      matrix_type;
   typedef dynamic_vector<scalar_type>      vector_type;

   const BQData& m_bqd;

 public:
   matrix_type data;

   eigval_mass_matrix_bq(const BQData& bqd)
     : m_bqd(bqd)
   {}

   void compute(const mesh_type& msh, const cell_type& cl)
   {
      auto cell_basis_size = howmany_dofs(m_bqd.cell_basis);
      auto cell_degree     = m_bqd.cell_degree();

      data = matrix_type::Zero(cell_basis_size, cell_basis_size);

      auto cell_quadpoints = m_bqd.cell_quadrature.integrate(msh, cl);
      for (auto& qp : cell_quadpoints) {
         matrix_type phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);
         data += qp.weight() * phi * phi.transpose();
      }
   }
};

template<typename Mesh,
         typename CellBasisType,
         typename CellQuadType,
         typename FaceBasisType,
         typename FaceQuadType>
class gradient_reconstruction
{
   typedef Mesh                            mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef typename mesh_type::face        face_type;

   typedef CellBasisType cell_basis_type;
   typedef CellQuadType  cell_quadrature_type;
   typedef FaceBasisType face_basis_type;
   typedef FaceQuadType  face_quadrature_type;

   typedef dynamic_matrix<scalar_type> matrix_type;
   typedef dynamic_vector<scalar_type> vector_type;

   typedef material_tensor<scalar_type, mesh_type::dimension, mesh_type::dimension>
     material_tensor_type;

   cell_basis_type      cell_basis;
   cell_quadrature_type cell_quadrature;

   face_basis_type      face_basis;
   face_quadrature_type face_quadrature;

   size_t m_degree;

 public:
   matrix_type oper;
   matrix_type data;

   gradient_reconstruction()
     : m_degree(1)
   {
      cell_basis      = cell_basis_type(m_degree + 1);
      cell_quadrature = cell_quadrature_type(2 * (m_degree + 1));
      face_basis      = face_basis_type(m_degree);
      face_quadrature = face_quadrature_type(2 * m_degree);
   }

   gradient_reconstruction(size_t degree)
     : m_degree(degree)
   {
      cell_basis      = cell_basis_type(m_degree + 1);
      cell_quadrature = cell_quadrature_type(2 * (m_degree + 1));
      face_basis      = face_basis_type(m_degree);
      face_quadrature = face_quadrature_type(2 * m_degree);
   }

   void compute(const mesh_type& msh, const cell_type& cl)
   {
      material_tensor_type id_tens;
      id_tens = material_tensor_type::Identity();
      compute(msh, cl, id_tens);
   }

   void compute(const mesh_type& msh, const cell_type& cl, const material_tensor_type& mtens)
   {
      matrix_type stiff_mat = matrix_type::Zero(cell_basis.size(), cell_basis.size());

      auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
      for (auto& qp : cell_quadpoints) {
         matrix_type dphi = cell_basis.eval_gradients(msh, cl, qp.point());
         stiff_mat += qp.weight() * dphi * (/*mtens **/ dphi.transpose());
      }

      /* LHS: take basis functions derivatives from degree 1 to K+1 */
      auto        MG_rowcol_range = cell_basis.range(1, m_degree + 1);
      matrix_type MG              = take(stiff_mat, MG_rowcol_range, MG_rowcol_range);

      /* RHS, volumetric part. */
      auto BG_row_range = cell_basis.range(1, m_degree + 1);
      auto BG_col_range = cell_basis.range(0, m_degree);

      auto fcs       = faces(msh, cl);
      auto num_faces = fcs.size();

      auto num_cell_dofs = cell_basis.range(0, m_degree).size();
      auto num_face_dofs = face_basis.size();

      dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

      matrix_type BG = matrix_type::Zero(BG_row_range.size(), dsr.total_size());

      BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
        take(stiff_mat, BG_row_range, BG_col_range);

      for (size_t face_i = 0; face_i < num_faces; face_i++) {
         auto current_face_range = dsr.face_range(face_i);
         auto fc                 = fcs[face_i];
         auto n                  = normal(msh, cl, fc);
         auto face_quadpoints    = face_quadrature.integrate(msh, fc);

         for (auto& qp : face_quadpoints) {
            matrix_type c_phi  = cell_basis.eval_functions(msh, cl, qp.point(), 0, m_degree);
            matrix_type c_dphi = cell_basis.eval_gradients(msh, cl, qp.point(), 1, m_degree + 1);

            matrix_type c_dphi_n = (c_dphi /** mtens*/) * n;
            matrix_type T        = qp.weight() * c_dphi_n * c_phi.transpose();

            BG.block(0, 0, BG.rows(), BG_col_range.size()) -= T;

            matrix_type f_phi = face_basis.eval_functions(msh, fc, qp.point());
            matrix_type F     = qp.weight() * c_dphi_n * f_phi.transpose();

            BG.block(0, current_face_range.min(), BG.rows(), current_face_range.size()) += F;
         }
      }

      oper = MG.llt().solve(BG);    // GT
      data = BG.transpose() * oper; // A
   }
};

template<typename Mesh,
         typename CellBasisType,
         typename CellQuadType,
         typename FaceBasisType,
         typename FaceQuadType,
         typename DivCellBasisType,
         typename DivCellQuadType>
class divergence_reconstruction
{
   typedef Mesh                            mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef typename mesh_type::face        face_type;

   typedef CellBasisType    cell_basis_type;
   typedef CellQuadType     cell_quadrature_type;
   typedef FaceBasisType    face_basis_type;
   typedef FaceQuadType     face_quadrature_type;
   typedef DivCellBasisType div_cell_basis_type;
   typedef DivCellQuadType  div_cell_quadrature_type;

   typedef dynamic_matrix<scalar_type> matrix_type;
   typedef dynamic_vector<scalar_type> vector_type;

   cell_basis_type      cell_basis;
   cell_quadrature_type cell_quadrature;

   face_basis_type      face_basis;
   face_quadrature_type face_quadrature;

   div_cell_basis_type      div_cell_basis;
   div_cell_quadrature_type div_cell_quadrature;

   size_t m_degree;

 public:
   matrix_type oper;
   matrix_type data;

   divergence_reconstruction()
     : m_degree(1)
   {
      cell_basis          = cell_basis_type(m_degree + 1);
      cell_quadrature     = cell_quadrature_type(2 * (m_degree + 1));
      face_basis          = face_basis_type(m_degree);
      face_quadrature     = face_quadrature_type(2 * m_degree);
      div_cell_basis      = div_cell_basis_type(m_degree);
      div_cell_quadrature = div_cell_quadrature_type(2 * m_degree);
   }

   divergence_reconstruction(size_t degree)
     : m_degree(degree)
   {
      cell_basis          = cell_basis_type(m_degree + 1);
      cell_quadrature     = cell_quadrature_type(2 * (m_degree + 1));
      face_basis          = face_basis_type(m_degree);
      face_quadrature     = face_quadrature_type(2 * m_degree);
      div_cell_basis      = div_cell_basis_type(m_degree);
      div_cell_quadrature = div_cell_quadrature_type(2 * m_degree);
   }

   void compute(const mesh_type& msh, const cell_type& cl)
   {
      auto        dcbs = div_cell_basis.size();
      matrix_type MD   = matrix_type::Zero(dcbs, dcbs);

      auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
      for (auto& qp : cell_quadpoints) {
         auto phi_d = div_cell_basis.eval_functions(msh, cl, qp.point());
         for (size_t i = 0; i < dcbs; i++)
            for (size_t j = 0; j < dcbs; j++)
               MD(i, j) += qp.weight() * mm_prod(phi_d[i], phi_d[j]);
      }

      /* RHS, volumetric part. */
      auto fcs           = faces(msh, cl);
      auto num_cell_dofs = cell_basis.range(0, m_degree).size();
      auto num_face_dofs = face_basis.size();
      auto num_faces     = fcs.size();

      dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

      matrix_type BD = matrix_type::Zero(dcbs, dsr.total_size());
      for (auto& qp : cell_quadpoints) {
         auto phi    = cell_basis.eval_functions(msh, cl, qp.point());
         auto dphi_d = div_cell_basis.eval_gradients(msh, cl, qp.point());

         for (size_t i = 0; i < dphi_d.size(); i++)
            for (size_t j = 0; j < num_cell_dofs; j++)
               BD(i, j) -= qp.weight() * mm_prod(dphi_d[i], phi[j]);
      }

      size_t face_offset = num_cell_dofs;

      for (size_t face_i = 0; face_i < num_faces; face_i++) {
         auto fc              = fcs[face_i];
         auto n               = normal(msh, cl, fc);
         auto face_quadpoints = face_quadrature.integrate(msh, fc);

         for (auto& qp : face_quadpoints) {
            auto phi_d = div_cell_basis.eval_functions(msh, cl, qp.point());
            auto phi   = face_basis.eval_functions(msh, fc, qp.point());
            for (size_t i = 0; i < phi_d.size(); i++) {
               for (size_t j = 0; j < face_basis.size(); j++) {
                  auto        p1 = mm_prod(phi[j], n);
                  scalar_type p2 = mm_prod(p1, phi_d[i]);
                  BD(i, face_offset + j) += qp.weight() * p2;
               }
            }
         }

         face_offset += face_basis.size();
      }

      oper = MD.partialPivLu().solve(BD);
      data = BD.transpose() * oper;
   }
};

} // namespace disk
