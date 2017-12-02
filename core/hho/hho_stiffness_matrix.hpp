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

namespace disk {

namespace hho {

// compute stiffness matrix
namespace priv {
template<typename BQData, typename BasisType>
struct stiffness_matrix_face_F
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::face        face_type;

   static void impl(const mesh_type& msh,
                    const face_type& fc,
                    const BQData&    bqd,
                    const size_t     degree)
   {
      static_assert(sizeof(BQData) == -1, "BQData not known in stiffness_matrix_face_F");
   }
};

// scaled_monomial_scalar_basis
template<typename BQData>
struct stiffness_matrix_face_F<
  BQData,
  scaled_monomial_scalar_basis<typename BQData::mesh_type, typename BQData::face_type>>
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::face        face_type;
   typedef dynamic_matrix<scalar_type>     matrix_type;

   static matrix_type impl(const mesh_type& msh,
                           const face_type& fc,
                           const BQData&    bqd,
                           const size_t     degree)
   {
      const auto face_basis_size = bqd.face_basis.range(0, degree).size();

      matrix_type stiffness = matrix_type::Zero(face_basis_size, face_basis_size);

      const auto face_quadpoints = bqd.face_quadrature.integrate(msh, fc);
      for (auto& qp : face_quadpoints) {
         const matrix_type dfphi = bqd.face_basis.eval_gradient(msh, fc, qp.point(), 0, degree);
         assert(dfphi.rows() == face_basis_size);

         stiffness += qp.weight() * dfphi * dfphi.transpose();
      }

      return stiffness;
   }
};

// scaled_monomial_vector_basis
template<typename BQData>
struct stiffness_matrix_face_F<
  BQData,
  scaled_monomial_vector_basis<typename BQData::mesh_type, typename BQData::face_type>>
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::face        face_type;
   typedef dynamic_matrix<scalar_type>     matrix_type;

   static matrix_type impl(const mesh_type& msh,
                           const face_type& fc,
                           const BQData&    bqd,
                           const size_t     degree)
   {
      const auto face_basis_size = bqd.face_basis.range(0, degree).size();

      matrix_type stiffness = matrix_type::Zero(face_basis_size, face_basis_size);

      const auto face_quadpoints = bqd.face_quadrature.integrate(msh, fc);
      for (auto& qp : face_quadpoints) {
         const auto dfphi = bqd.face_basis.eval_gradients(msh, fc, qp.point(), 0, degree);
         assert(dfphi.size() == face_basis_size);

         for (size_t j = 0; j < face_basis_size; j++) {
            for (size_t i = j; i < face_basis_size; i++) {
               stiffness(i, j) += qp.weight() * mm_prod(dfphi[i], dfphi[j]);
            }
         }
      }

      // lower part
      for (size_t j = 0; j < face_basis_size; j++) {
         for (size_t i = 0; i < j; i++) {
            stiffness(i, j) = stiffness(j, i);
         }
      }

      return stiffness;
   }
};

template<typename BQData, typename BasisType>
struct stiffness_matrix_cell_F
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   static void impl(const mesh_type& msh,
                    const cell_type& cl,
                    const BQData&    bqd,
                    const size_t     degree)
   {
      static_assert(sizeof(BQData) == -1, "BQData not known in stiffness_matrix_cell_F");
   }
};

// scaled_monomial_scalar_basis
template<typename BQData>
struct stiffness_matrix_cell_F<
  BQData,
  scaled_monomial_scalar_basis<typename BQData::mesh_type, typename BQData::cell_type>>
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef dynamic_matrix<scalar_type>     matrix_type;

   static matrix_type impl(const mesh_type& msh,
                           const cell_type& cl,
                           const BQData&    bqd,
                           const size_t     degree)
   {
      const auto cell_basis_size = bqd.cell_basis.range(0, degree).size();

      matrix_type stiffness = matrix_type::Zero(cell_basis_size, cell_basis_size);

      const auto cell_quadpoints = bqd.cell_quadrature.integrate(msh, cl);
      for (auto& qp : cell_quadpoints) {
         const matrix_type dcphi = bqd.cell_basis.eval_gradients(msh, cl, qp.point(), 0, degree);
         assert(dcphi.rows() == cell_basis_size);

         stiffness += qp.weight() * dcphi * dcphi.transpose();
      }

      return stiffness;
   }
};

// scaled_monomial_vector_basis
template<typename BQData>
struct stiffness_matrix_cell_F<
  BQData,
  scaled_monomial_vector_basis<typename BQData::mesh_type, typename BQData::cell_type>>
{

   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef dynamic_matrix<scalar_type>     matrix_type;

   static matrix_type impl(const mesh_type& msh,
                           const cell_type& cl,
                           const BQData&    bqd,
                           const size_t     degree)
   {
      const auto cell_basis_size = bqd.cell_basis.range(0, degree).size();

      matrix_type stiffness = matrix_type::Zero(cell_basis_size, cell_basis_size);

      const auto cell_quadpoints = bqd.cell_quadrature.integrate(msh, cl);
      for (auto& qp : cell_quadpoints) {
         const auto dcphi = bqd.cell_basis.eval_gradients(msh, cl, qp.point(), 0, degree);
         assert(dcphi.size() == cell_basis_size);

         for (size_t j = 0; j < cell_basis_size; j++) {
            for (size_t i = j; i < cell_basis_size; i++) {
               stiffness(i, j) += qp.weight() * mm_prod(dcphi[i], dcphi[j]);
            }
         }
      }

      // lower part
      for (size_t j = 0; j < cell_basis_size; j++) {
         for (size_t i = 0; i < j; i++) {
            stiffness(i, j) = stiffness(j, i);
         }
      }

      return stiffness;
   }
};

} // namespace priv

// face stiffness matrix
template<typename BQData>
dynamic_matrix<typename BQData::mesh_type::scalar_type>
stiffness_matrix(const typename BQData::mesh_type& msh,
                 const typename BQData::face_type& fc,
                 const BQData&                     bqd,
                 int                               degree = -1)
{
   if (degree < 0) degree = bqd.face_degree();
   return priv::stiffness_matrix_face_F<BQData, typename BQData::face_basis_type>::impl(
     msh, fc, bqd, degree);
}

// cell stiffness matrix
template<typename BQData>
dynamic_matrix<typename BQData::mesh_type::scalar_type>
stiffness_matrix(const typename BQData::mesh_type& msh,
                 const typename BQData::cell_type& cl,
                 const BQData&                     bqd,
                 int                               degree = -1)
{
   if (degree < 0) degree = bqd.cell_degree();
   return priv::stiffness_matrix_cell_F<BQData, typename BQData::cell_basis_type>::impl(
     msh, cl, bqd, degree);
}

} // hho namespace
} // disk namespace