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
#include "quadratures/quadratures.hpp"

namespace disk {

namespace hho {

namespace priv {

template<typename BQData, typename BasisType>
struct stiffness_matrix_cell_F
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   static void
   impl(const mesh_type& msh,
        const cell_type& cl,
        const BQData&    bqd,
        const size_t&    min_degree,
        const size_t&    max_degree)
   {
      static_assert(sizeof(BQData) == -1, "BQData not known in stiffness_matrix_cell_F");
   }

   template<typename TensorField>
   static void
   impl(const mesh_type&   msh,
        const cell_type&   cl,
        const TensorField& mtens,
        const BQData&      bqd,
        const size_t&      min_degree,
        const size_t&      max_degree)
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

   static matrix_type
   impl(const mesh_type& msh,
        const cell_type& cl,
        const BQData&    bqd,
        const size_t&    min_degree,
        const size_t&    max_degree)
   {
      const auto cell_basis_size = bqd.cell_basis.range(min_degree, max_degree).size();

      matrix_type stiffness = matrix_type::Zero(cell_basis_size, cell_basis_size);

      const auto cell_quadpoints = bqd.cell_quadrature.integrate(msh, cl);
      assert(2 * max_degree <= bqd.cell_quadrature.order());

      for (auto& qp : cell_quadpoints) {
         const matrix_type dcphi =
           bqd.cell_basis.eval_gradients(msh, cl, qp.point(), min_degree, max_degree);
         assert(dcphi.rows() == cell_basis_size);

         stiffness += qp.weight() * dcphi * dcphi.transpose();
      }

      return stiffness;
   }

   template<typename TensorField>
   static matrix_type
   impl(const mesh_type&   msh,
        const cell_type&   cl,
        const TensorField& mtens,
        const BQData&      bqd,
        const size_t&      min_degree,
        const size_t&      max_degree)
   {
      const auto cell_basis_size = bqd.cell_basis.range(min_degree, max_degree).size();

      matrix_type stiffness = matrix_type::Zero(cell_basis_size, cell_basis_size);

      const auto cell_quadpoints = bqd.cell_quadrature.integrate(msh, cl);
      assert(2 * max_degree <= bqd.cell_quadrature.order());

      for (auto& qp : cell_quadpoints) {
         const matrix_type dcphi =
           bqd.cell_basis.eval_gradients(msh, cl, qp.point(), min_degree, max_degree);
         assert(dcphi.rows() == cell_basis_size);

         stiffness += qp.weight() * dcphi * (dcphi * mtens(qp.point())).transpose();
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

   static matrix_type
   impl(const mesh_type& msh,
        const cell_type& cl,
        const BQData&    bqd,
        const size_t&    min_degree,
        const size_t&    max_degree)
   {
      const size_t dimension       = mesh_type::dimension;
      const auto   cell_basis_size = bqd.cell_basis.range(min_degree, max_degree).size();

      matrix_type stiffness = matrix_type::Zero(cell_basis_size, cell_basis_size);

      const auto cell_quadpoints = bqd.cell_quadrature.integrate(msh, cl);
      assert(2 * max_degree <= bqd.cell_quadrature.order());

      for (auto& qp : cell_quadpoints) {
         const auto dcphi =
           bqd.cell_basis.eval_gradients(msh, cl, qp.point(), min_degree, max_degree);
         assert(dcphi.size() == cell_basis_size);

         for (size_t j = 0; j < cell_basis_size; j++) {
            const auto qp_dcphi_j = mm_prod(qp.weight(), dcphi[j]);
            for (size_t i = j; i < cell_basis_size; i += dimension) {
               stiffness(i, j) += mm_prod(dcphi[i], qp_dcphi_j);
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

// cell stiffness matrix
template<typename BQData>
dynamic_matrix<typename BQData::mesh_type::scalar_type>
stiffness_matrix(const typename BQData::mesh_type& msh,
                 const typename BQData::cell_type& cl,
                 const BQData&                     bqd,
                 size_t                            min_degree = 0,
                 int                               max_degree = -1)
{
   if (max_degree < 0) max_degree = bqd.cell_degree();
   return priv::stiffness_matrix_cell_F<BQData, typename BQData::cell_basis_type>::impl(
     msh, cl, bqd, min_degree, max_degree);
}

template<typename BQData, typename TensorField>
dynamic_matrix<typename BQData::mesh_type::scalar_type>
stiffness_matrix(const typename BQData::mesh_type& msh,
                 const typename BQData::cell_type& cl,
                 const TensorField&                mtens,
                 const BQData&                     bqd,
                 size_t                            min_degree = 0,
                 int                               max_degree = -1)
{
   if (max_degree < 0) max_degree = bqd.cell_degree();
   return priv::stiffness_matrix_cell_F<BQData, typename BQData::cell_basis_type>::impl(
     msh, cl, mtens, bqd, min_degree, max_degree);
}

// symetric stiffness matrix
template<typename BQData>
dynamic_matrix<typename BQData::mesh_type::scalar_type>
sym_stiffness_matrix(const typename BQData::mesh_type& msh,
                     const typename BQData::cell_type& cl,
                     const BQData&                     bqd,
                     size_t                            min_degree = 0,
                     int                               max_degree = -1)
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef dynamic_matrix<scalar_type>     matrix_type;

   if (max_degree < 0) max_degree = bqd.cell_degree();

   const auto cell_basis_size = bqd.cell_basis.range(min_degree, max_degree).size();

   matrix_type stiffness = matrix_type::Zero(cell_basis_size, cell_basis_size);

   const auto cell_quadpoints = bqd.cell_quadrature.integrate(msh, cl);
   assert(2 * max_degree <= bqd.cell_quadrature.order());

   for (auto& qp : cell_quadpoints) {
      const auto dcphi =
        bqd.cell_basis.eval_sgradients(msh, cl, qp.point(), min_degree, max_degree);
      assert(dcphi.size() == cell_basis_size);

      for (size_t j = 0; j < cell_basis_size; j++) {
         const auto qp_dcphi_j = mm_prod(qp.weight(), dcphi[j]);
         for (size_t i = j; i < cell_basis_size; i++) {
            stiffness(i, j) += mm_prod(dcphi[i], qp_dcphi_j);
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

} // hho namespace
} // disk namespace