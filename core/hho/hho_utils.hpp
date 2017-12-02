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

// scalar case
template<typename T>
T
eval(const dynamic_vector<T>& tab_coeff, const dynamic_vector<T>& base)
{
   assert(tab_coeff.rows() == base.rows());

   return tab_coeff.dot(base);
}

template<typename T, int DIM>
static_vector<T, DIM>
eval(const dynamic_vector<T>& tab_coeff, const Eigen::Matrix<T, Eigen::Dynamic, DIM>& base)
{
   assert(tab_coeff.rows() == base.rows());

   const auto prod = tab_coeff * base;

   static_vector<T, DIM> ret = static_vector<T, DIM>::Zero();

   for (size_t i = 0; i < DIM; i++) {
      ret(i) = prod(i);
   }

   return ret;
}

// vectorial case
template<typename T, int DIM>
static_matrix<T, DIM, DIM>
eval(const dynamic_vector<T>& tab_coeff, const std::vector<static_matrix<T, DIM, DIM>>& base)
{
   static_matrix<T, DIM, DIM> ret = static_matrix<T, DIM, DIM>::Zero();
   assert(tab_coeff.size() == base.size());

   for (size_t i = 0; i < base.size(); i++) {
      ret += tab_coeff(i) * base[i];
   }

   return ret;
}

// matricial case
template<typename T, int DIM>
static_vector<T, DIM>
eval(const dynamic_vector<T>& tab_coeff, const std::vector<static_vector<T, DIM>>& base)
{
   static_vector<T, DIM> ret = static_vector<T, DIM>::Zero();
   assert(tab_coeff.size() == base.size());

   for (size_t i = 0; i < base.size(); i++) {
      ret += tab_coeff(i) * base[i];
   }

   return ret;
}

namespace priv {
// compute_rhs
template<typename BQData, typename BasisType, typename Function>
struct compute_rhs_face_F
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::face        face_type;

   static void impl(const mesh_type& msh,
                    const face_type& fc,
                    const Function&  func,
                    const BQData&    bqd,
                    const size_t&    degree)
   {
      static_assert(sizeof(BQData) == -1, "BQData not known in compute_rhs_face_F");
   }
};

// scaled_monomial_scalar_basis
template<typename BQData, typename Function>
struct compute_rhs_face_F<
  BQData,
  scaled_monomial_scalar_basis<typename BQData::mesh_type, typename BQData::face_type>,
  Function>
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::face        face_type;
   typedef dynamic_vector<scalar_type>     vector_type;
   typedef dynamic_matrix<scalar_type>     matrix_type;

   static vector_type impl(const mesh_type& msh,
                           const face_type& fc,
                           const Function&  func,
                           const BQData&    bqd,
                           const size_t&    degree)
   {
      const auto face_basis_size = bqd.face_basis.range(0, degree).size();

      vector_type vec = vector_type::Zero(face_basis_size);

      const auto face_quadpoints = bqd.face_quadrature.integrate(msh, fc);
      for (auto& qp : face_quadpoints) {
         const matrix_type fphi = bqd.face_basis.eval_functions(msh, fc, qp.point(), 0, degree);
         assert(fphi.rows() == face_basis_size);

         vec += qp.weight() * fphi * func(qp.point());
      }

      return vec;
   }
};

// scaled_monomial_vector_basis
template<typename BQData, typename Function>
struct compute_rhs_face_F<
  BQData,
  scaled_monomial_vector_basis<typename BQData::mesh_type, typename BQData::face_type>,
  Function>
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::face        face_type;
   typedef dynamic_vector<scalar_type>     vector_type;

   static vector_type impl(const mesh_type& msh,
                           const face_type& fc,
                           const Function&  func,
                           const BQData&    bqd,
                           const size_t&    degree)
   {
      const auto face_basis_size = bqd.face_basis.range(0, degree).size();

      vector_type vec = vector_type::Zero(face_basis_size);

      const auto face_quadpoints = bqd.face_quadrature.integrate(msh, fc);
      for (auto& qp : face_quadpoints) {
         const auto fphi = bqd.face_basis.eval_functions(msh, fc, qp.point(), 0, degree);
         assert(fphi.size() == face_basis_size);

         const auto fval = func(qp.point());
         for (size_t i = 0; i < face_basis_size; i++) {
            vec(i) += qp.weight() * mm_prod(fphi[i], fval);
         }
      }

      return vec;
   }
};

template<typename BQData, typename BasisType, typename Function>
struct compute_rhs_cell_F
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   static void impl(const mesh_type& msh,
                    const cell_type& cl,
                    const Function&  func,
                    const BQData&    bqd,
                    const size_t&    degree)
   {
      static_assert(sizeof(BQData) == -1, "BQData not known in compute_rhs_cell_F");
   }
};

// scaled_monomial_scalar_basis
template<typename BQData, typename Function>
struct compute_rhs_cell_F<
  BQData,
  scaled_monomial_scalar_basis<typename BQData::mesh_type, typename BQData::cell_type>,
  Function>
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef dynamic_vector<scalar_type>     vector_type;

   static vector_type impl(const mesh_type& msh,
                           const cell_type& cl,
                           const Function&  func,
                           const BQData&    bqd,
                           const size_t&    degree)
   {
      const auto cell_basis_size = bqd.cell_basis.range(0, degree).size();

      vector_type vec = vector_type::Zero(cell_basis_size);

      const auto cell_quadpoints = bqd.cell_quadrature.integrate(msh, cl);
      for (auto& qp : cell_quadpoints) {
         const vector_type cphi = bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, degree);
         assert(cphi.rows() == cell_basis_size);

         vec += qp.weight() * cphi * func(qp.point());
      }

      return vec;
   }
};

// scaled_monomial_vector_basis
template<typename BQData, typename Function>
struct compute_rhs_cell_F<
  BQData,
  scaled_monomial_vector_basis<typename BQData::mesh_type, typename BQData::cell_type>,
  Function>
{
   typedef typename BQData::mesh_type      mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;
   typedef dynamic_vector<scalar_type>     vector_type;

   static vector_type impl(const mesh_type& msh,
                           const cell_type& cl,
                           const Function&  func,
                           const BQData&    bqd,
                           const size_t&    degree)
   {
      const auto cell_basis_size = bqd.cell_basis.range(0, degree).size();

      vector_type vec = vector_type::Zero(cell_basis_size);

      const auto cell_quadpoints = bqd.cell_quadrature.integrate(msh, cl);
      for (auto& qp : cell_quadpoints) {
         const auto cphi = bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, degree);
         assert(cphi.size() == cell_basis_size);

         const auto feval = func(qp.point());
         for (size_t i = 0; i < cell_basis_size; i++) {
            vec(i) += qp.weight() * mm_prod(cphi[i], feval);
         }
      }

      return vec;
   }
};

} // namespace priv

// face compute_rhs
template<typename BQData, typename Function>
dynamic_vector<typename BQData::mesh_type::scalar_type>
compute_rhs(const typename BQData::mesh_type& msh,
            const typename BQData::face_type& fc,
            const Function&                   func,
            const BQData&                     bqd,
            int                               degree = -1)
{
   if (degree < 0) degree = bqd.face_degree();
   return priv::compute_rhs_face_F<BQData, typename BQData::face_basis_type, Function>::impl(
     msh, fc, func, bqd, degree);
}

// cell compute_rhs
template<typename BQData, typename Function>
dynamic_vector<typename BQData::mesh_type::scalar_type>
compute_rhs(const typename BQData::mesh_type& msh,
            const typename BQData::cell_type& cl,
            const Function&                   func,
            const BQData&                     bqd,
            int                               degree = -1)
{
   if (degree < 0) degree = bqd.cell_degree();
   return priv::compute_rhs_cell_F<BQData, typename BQData::cell_basis_type, Function>::impl(
     msh, cl, func, bqd, degree);
}

} // hho namespace
} // disk namespace