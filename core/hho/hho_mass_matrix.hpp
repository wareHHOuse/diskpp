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

#include "bases/bases.hpp"
#include "bases/bases_traits.hpp"
#include "common/eigen.hpp"
#include "quadratures/quadratures.hpp"
#include <cassert>

namespace disk {

namespace hho {

// compute mass matrix
namespace priv {
template<typename BQData, typename BasisType>
struct mass_matrix_face_F
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::face        face_type;

   static void
   impl(const mesh_type& msh,
        const face_type& fc,
        const BQData&    bqd,
        const size_t&    min_degree,
        const size_t&    max_degree)
   {
      static_assert(sizeof(BQData) == -1, "BQData not known in mass_matrix_face_F");
   }
};

// scaled_monomial_scalar_basis
template<typename BQData>
struct mass_matrix_face_F<
  BQData,
  scaled_monomial_scalar_basis<typename BQData::mesh_type, typename BQData::face_type>>
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::face        face_type;
   typedef dynamic_matrix<scalar_type>     matrix_type;

   static matrix_type
   impl(const mesh_type& msh,
        const face_type& fc,
        const BQData&    bqd,
        const size_t&    min_degree,
        const size_t&    max_degree)
   {
      const auto face_basis_size = bqd.face_basis.range(min_degree, max_degree).size();

      matrix_type mass = matrix_type::Zero(face_basis_size, face_basis_size);

      const auto face_quadpoints = bqd.face_quadrature.integrate(msh, fc);
      assert((msh.dimension == 1 && bqd.face_quadrature.order() == 0) ||
             2 * max_degree <= bqd.face_quadrature.order());
      for (auto& qp : face_quadpoints) {
         const matrix_type fphi =
           bqd.face_basis.eval_functions(msh, fc, qp.point(), min_degree, max_degree);
         assert(fphi.rows() == face_basis_size);

         mass += qp.weight() * fphi * fphi.transpose();
      }

      return mass;
   }
};

// scaled_monomial_vector_basis
template<typename BQData>
struct mass_matrix_face_F<
  BQData,
  scaled_monomial_vector_basis<typename BQData::mesh_type, typename BQData::face_type>>
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::face        face_type;
   typedef dynamic_matrix<scalar_type>     matrix_type;

   const static size_t dimension = mesh_type::dimension;

   static matrix_type
   impl(const mesh_type& msh,
        const face_type& fc,
        const BQData&    bqd,
        const size_t&    min_degree,
        const size_t&    max_degree)
   {
      const auto face_basis_size = bqd.face_basis.range(min_degree, max_degree).size();

      matrix_type mass = matrix_type::Zero(face_basis_size, face_basis_size);

      const auto face_quadpoints = bqd.face_quadrature.integrate(msh, fc);
      assert((msh.dimension == 1 && bqd.face_quadrature.order() == 0) ||
             (2 * max_degree <= bqd.face_quadrature.order()));
      for (auto& qp : face_quadpoints) {
         const auto fphi =
           bqd.face_basis.eval_functions(msh, fc, qp.point(), min_degree, max_degree);
         assert(fphi.size() == face_basis_size);

         for (size_t j = 0; j < face_basis_size; j++) {
            const auto qp_fhi_j = mm_prod(qp.weight(), fphi[j]);
            for (size_t i = j; i < face_basis_size; i += dimension) {
               mass(i, j) += mm_prod(fphi[i], qp_fhi_j);
            }
         }
      }

      // upper part
      for (size_t j = 0; j < face_basis_size; j++) {
         for (size_t i = 0; i < j; i++) {
            mass(i, j) = mass(j, i);
         }
      }

      return mass;
   }
};

template<typename BQData, typename BasisType>
struct mass_matrix_cell_F
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
      static_assert(sizeof(BQData) == -1, "BQData not known in mass_matrix_cell_F");
   }
};

// scaled_monomial_scalar_basis
template<typename BQData>
struct mass_matrix_cell_F<
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

      matrix_type mass = matrix_type::Zero(cell_basis_size, cell_basis_size);

      const auto cell_quadpoints = bqd.cell_quadrature.integrate(msh, cl);
      assert(2 * max_degree <= bqd.cell_quadrature.order());
      for (auto& qp : cell_quadpoints) {
         const matrix_type cphi =
           bqd.cell_basis.eval_functions(msh, cl, qp.point(), min_degree, max_degree);
         assert(cphi.rows() == cell_basis_size);

         mass += qp.weight() * cphi * cphi.transpose();
      }

      return mass;
   }
};

// scaled_monomial_vector_basis
template<typename BQData>
struct mass_matrix_cell_F<
  BQData,
  scaled_monomial_vector_basis<typename BQData::mesh_type, typename BQData::cell_type>>
{

   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef dynamic_matrix<scalar_type>     matrix_type;

   const static size_t dimension = mesh_type::dimension;

   static matrix_type
   impl(const mesh_type& msh,
        const cell_type& cl,
        const BQData&    bqd,
        const size_t&    min_degree,
        const size_t&    max_degree)
   {
      const auto cell_basis_size = bqd.cell_basis.range(min_degree, max_degree).size();

      matrix_type mass = matrix_type::Zero(cell_basis_size, cell_basis_size);

      const auto cell_quadpoints = bqd.cell_quadrature.integrate(msh, cl);
      assert(2 * max_degree <= bqd.cell_quadrature.order());
      for (auto& qp : cell_quadpoints) {
         const auto cphi =
           bqd.cell_basis.eval_functions(msh, cl, qp.point(), min_degree, max_degree);
         assert(cphi.size() == cell_basis_size);

         for (size_t j = 0; j < cell_basis_size; j++) {
            const auto qp_cphi_j = mm_prod(qp.weight(), cphi[j]);
            for (size_t i = j; i < cell_basis_size; i += dimension) {
               mass(i, j) += mm_prod(cphi[i], qp_cphi_j);
            }
         }
      }

      // upper part
      for (size_t j = 0; j < cell_basis_size; j++) {
         for (size_t i = 0; i < j; i++) {
            mass(i, j) = mass(j, i);
         }
      }

      return mass;
   }
};

} // namespace priv

// face mass matrix
template<typename BQData>
dynamic_matrix<typename BQData::mesh_type::scalar_type>
mass_matrix(const typename BQData::mesh_type& msh,
            const typename BQData::face_type& fc,
            const BQData&                     bqd,
            size_t                            min_degree = 0,
            size_t                            max_degree = VERY_HIGH_DEGREE)
{
   max_degree = (max_degree == VERY_HIGH_DEGREE) ? bqd.face_degree() : max_degree;

   return priv::mass_matrix_face_F<BQData, typename BQData::face_basis_type>::impl(
     msh, fc, bqd, min_degree, max_degree);
}

// cell mass matrix
template<typename BQData>
dynamic_matrix<typename BQData::mesh_type::scalar_type>
mass_matrix(const typename BQData::mesh_type& msh,
            const typename BQData::cell_type& cl,
            const BQData&                     bqd,
            size_t                            min_degree = 0,
            size_t                            max_degree = VERY_HIGH_DEGREE)
{
   max_degree = (max_degree == VERY_HIGH_DEGREE) ? bqd.cell_degree() : max_degree;

   return priv::mass_matrix_cell_F<BQData, typename BQData::cell_basis_type>::impl(
     msh, cl, bqd, min_degree, max_degree);
}

// gradient mass_matrix

namespace priv {

template<typename BQData, typename GradBasis>
struct grad_mass_matrix_cell_F
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef dynamic_matrix<scalar_type>     matrix_type;

   const static size_t dimension = mesh_type::dimension;

   static matrix_type
   impl(const mesh_type& msh,
        const cell_type& cl,
        const BQData&    bqd,
        const size_t&    min_degree,
        const size_t&    max_degree)
   {
      static_assert(use_vector_container<typename BQData::grad_basis_type>::value,
                    "cell basic function have to use vector container");
      const auto grad_basis_size = bqd.grad_basis.range(min_degree, max_degree).size();

      matrix_type mass = matrix_type::Zero(grad_basis_size, grad_basis_size);

      const auto grad_quadpoints = bqd.grad_quadrature.integrate(msh, cl);
      assert(2 * max_degree <= bqd.grad_quadrature.order());

      for (auto& qp : grad_quadpoints) {
         const auto gphi =
           bqd.grad_basis.eval_functions(msh, cl, qp.point(), min_degree, max_degree);
         assert(gphi.size() == grad_basis_size);

         for (size_t j = 0; j < grad_basis_size; j++) {
            const auto qp_gphi_j = mm_prod(qp.weight(), gphi[j]);
            for (size_t i = j; i < grad_basis_size; i++) {
               mass(i, j) += mm_prod(gphi[i], qp_gphi_j);
            }
         }
      }

      // lower part
      for (size_t j = 0; j < grad_basis_size; j++) {
         for (size_t i = 0; i < j; i++) {
            mass(i, j) = mass(j, i);
         }
      }

      return mass;
   }
};

template<typename BQData>
struct grad_mass_matrix_cell_F<
  BQData,
  scaled_monomial_matrix_basis<typename BQData::mesh_type, typename BQData::cell_type>>
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef dynamic_matrix<scalar_type>     matrix_type;

   const static size_t dimension = mesh_type::dimension;

   static matrix_type
   impl(const mesh_type& msh,
        const cell_type& cl,
        const BQData&    bqd,
        const size_t&    min_degree,
        const size_t&    max_degree)
   {
      static_assert(use_vector_container<typename BQData::grad_basis_type>::value,
                    "cell basic function have to use vector container");
      const auto grad_basis_size = bqd.grad_basis.range(min_degree, max_degree).size();

      matrix_type mass = matrix_type::Zero(grad_basis_size, grad_basis_size);

      const auto grad_quadpoints = bqd.grad_quadrature.integrate(msh, cl);
      assert(2 * max_degree <= bqd.grad_quadrature.order());

      const size_t dim2 = dimension * dimension;

      for (auto& qp : grad_quadpoints) {
         const auto gphi =
           bqd.grad_basis.eval_functions(msh, cl, qp.point(), min_degree, max_degree);
         assert(gphi.size() == grad_basis_size);

         for (size_t j = 0; j < grad_basis_size; j += dim2) {
            const auto qp_gphi_j = mm_prod(qp.weight(), gphi[j]);
            for (size_t i = j; i < grad_basis_size; i += dim2) {
               mass(i, j) += mm_prod(gphi[i], qp_gphi_j);
            }
         }
      }

      // upper part
      for (size_t j = 0; j < grad_basis_size; j++) {
         for (size_t i = 0; i < j; i++) {
            mass(i, j) = mass(j, i);
         }
      }

      // copy of each cols
      for (size_t j = 0; j < grad_basis_size; j += dim2) {
         for (size_t col = 1; col < dim2; col++) {
            mass.block(col, j + col, grad_basis_size - dim2 + 1, 1) =
              mass.block(0, j, grad_basis_size - dim2 + 1, 1);
         }
      }

      return mass;
   }
};

template<typename BQData>
struct grad_mass_matrix_cell_F<
  BQData,
  scaled_monomial_sym_matrix_basis<typename BQData::mesh_type, typename BQData::cell_type>>
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef dynamic_matrix<scalar_type>     matrix_type;

   const static size_t dimension = mesh_type::dimension;

   static matrix_type
   impl(const mesh_type& msh,
        const cell_type& cl,
        const BQData&    bqd,
        const size_t&    min_degree,
        const size_t&    max_degree)
   {
      static_assert(use_vector_container<typename BQData::grad_basis_type>::value,
                    "cell basic function have to use vector container");
      const auto grad_basis_size = bqd.grad_basis.range(min_degree, max_degree).size();

      matrix_type mass = matrix_type::Zero(grad_basis_size, grad_basis_size);

      const auto grad_quadpoints = bqd.grad_quadrature.integrate(msh, cl);
      assert(2 * max_degree <= bqd.grad_quadrature.order());

      size_t dec = 0;
      if (dimension == 3) {
         dec = 6;
      } else if (dimension == 2) {
         dec = 3;
      } else
         assert(false);

      for (auto& qp : grad_quadpoints) {
         const auto gphi =
           bqd.grad_basis.eval_functions(msh, cl, qp.point(), min_degree, max_degree);
         assert(gphi.size() == grad_basis_size);

         for (size_t j = 0; j < grad_basis_size; j++) {
            const auto qp_gphi_j = mm_prod(qp.weight(), gphi[j]);
            for (size_t i = j; i < grad_basis_size; i += dec) {
               mass(i, j) += mm_prod(gphi[i], qp_gphi_j);
            }
         }
      }

      // upper part
      for (size_t j = 0; j < grad_basis_size; j++) {
         for (size_t i = 0; i < j; i++) {
            mass(i, j) = mass(j, i);
         }
      }

      return mass;
   }
};

} // end priv

template<typename BQData>
dynamic_matrix<typename BQData::mesh_type::scalar_type>
grad_mass_matrix(const typename BQData::mesh_type& msh,
                 const typename BQData::cell_type& cl,
                 const BQData&                     bqd,
                 size_t                            min_degree = 0,
                 size_t                            max_degree = VERY_HIGH_DEGREE)
{
   static_assert(have_grad_space<BQData>::value, "You need to define a gradient space");
   max_degree = (max_degree == VERY_HIGH_DEGREE) ? bqd.cell_degree() : max_degree;

   return priv::grad_mass_matrix_cell_F<BQData, typename BQData::grad_basis_type>::impl(
     msh, cl, bqd, min_degree, max_degree);
}

} // hho namespace
} // disk namespace