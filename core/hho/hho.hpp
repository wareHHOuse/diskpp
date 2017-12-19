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
#include "hho/hho_bq.hpp"
#include "hho/hho_mass_matrix.hpp"

namespace disk {
namespace hho {

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

   eigval_mass_matrix_bq(const BQData& bqd) : m_bqd(bqd) {}

   void
   compute(const mesh_type& msh, const cell_type& cl)
   {
      auto cell_basis_size = howmany_dofs(m_bqd.cell_basis);
      auto cell_degree     = m_bqd.cell_degree();

      data = mass_matrix(msh, cl, m_bqd, cell_degree);
      assert(data.rows() == data.cols() && data.rows() == cell_basis_size);
   }
};

} // end hho
} // namespace disk
