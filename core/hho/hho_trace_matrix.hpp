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

// compute mass matrix
namespace priv {
template<typename BQData, typename BasisType>
struct trace_matrix_F
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef typename mesh_type::face        face_type;

   static void
   impl(const mesh_type& msh,
        const cell_type& cl,
        const face_type& fc,
        const BQData&    bqd,
        const size_t&    cell_degree,
        const size_t&    face_degree)
   {
      static_assert(sizeof(BQData) == -1, "BQData not known in trace_matrix_F");
   }
};

// scaled_monomial_scalar_basis
template<typename BQData>
struct trace_matrix_F<
  BQData,
  scaled_monomial_scalar_basis<typename BQData::mesh_type, typename BQData::face_type>>
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef typename mesh_type::face        face_type;
   typedef dynamic_matrix<scalar_type>     matrix_type;

   static matrix_type
   impl(const mesh_type& msh,
        const cell_type& cl,
        const face_type& fc,
        const BQData&    bqd,
        const size_t&    cell_degree,
        const size_t&    face_degree)
   {
      const auto cell_basis_size = bqd.cell_basis.range(0, cell_degree).size();
      const auto face_basis_size = bqd.face_basis.range(0, face_degree).size();

      matrix_type mat = matrix_type::Zero(face_basis_size, cell_basis_size);

      const auto face_quadpoints = bqd.face_trace_quadrature.integrate(msh, fc);
      assert((msh.dimension == 1 && bqd.face_trace_quadrature.order() == 0) ||
             (cell_degree + face_degree) <= bqd.face_trace_quadrature.order());

      for (auto& qp : face_quadpoints) {
         const matrix_type fphi =
           bqd.face_basis.eval_functions(msh, fc, qp.point(), 0, face_degree);
         assert(fphi.rows() == face_basis_size);

         const auto cphi = bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);
         assert(cphi.rows() == cell_basis_size);

         mat += qp.weight() * fphi * cphi.transpose();
      }

      return mat;
   }
};

// scaled_monomial_vector_basis
template<typename BQData>
struct trace_matrix_F<
  BQData,
  scaled_monomial_vector_basis<typename BQData::mesh_type, typename BQData::face_type>>
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef typename mesh_type::face        face_type;
   typedef dynamic_matrix<scalar_type>     matrix_type;

   const static size_t dimension = mesh_type::dimension;

   static matrix_type
   impl(const mesh_type& msh,
        const cell_type& cl,
        const face_type& fc,
        const BQData&    bqd,
        const size_t&    cell_degree,
        const size_t&    face_degree)
   {
      const auto cell_basis_size = bqd.cell_basis.range(0, cell_degree).size();
      const auto face_basis_size = bqd.face_basis.range(0, face_degree).size();

      matrix_type mat = matrix_type::Zero(face_basis_size, cell_basis_size);

      const auto face_quadpoints = bqd.face_trace_quadrature.integrate(msh, fc);
      assert((msh.dimension == 1 && bqd.face_trace_quadrature.order() == 0) ||
             (cell_degree + face_degree) <= bqd.face_trace_quadrature.order());

      for (auto& qp : face_quadpoints) {
         const auto fphi = bqd.face_basis.eval_functions(msh, fc, qp.point(), 0, face_degree);
         assert(fphi.size() == face_basis_size);

         const auto cphi = bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);
         assert(cphi.size() == cell_basis_size);

         for (size_t j = 0; j < cell_basis_size; j += dimension) {
            for (size_t col = 0; col < dimension; col++) {
               const auto qp_cphi_j = mm_prod(qp.weight(), cphi[j + col]);
               for (size_t i = col; i < face_basis_size; i += dimension) {
                  mat(i, j + col) += mm_prod(fphi[i], qp_cphi_j);
               }
            }
         }
      }

      return mat;
   }
};

} // namespace priv

// face mass matrix
template<typename BQData>
dynamic_matrix<typename BQData::mesh_type::scalar_type>
trace_matrix(const typename BQData::mesh_type& msh,
             const typename BQData::cell_type& cl,
             const typename BQData::face_type& fc,
             const BQData&                     bqd,
             int                               cell_degree = -1,
             int                               face_degree = -1)
{
   if (cell_degree < 0 or face_degree < 0) {
      cell_degree = bqd.cell_degree();
      face_degree = bqd.face_degree();
   }
   return priv::trace_matrix_F<BQData, typename BQData::face_basis_type>::impl(
     msh, cl, fc, bqd, cell_degree, face_degree);
}

} // hho namespace
} // disk namespace