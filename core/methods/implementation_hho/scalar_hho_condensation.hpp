/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
 * Nicolas Pignet  (C) 2018, 2019               nicolas.pignet@enpc.fr
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

#include <vector>

#include "adaptivity/adaptivity.hpp"
#include "bases/bases.hpp"
#include "boundary_conditions/boundary_conditions.hpp"
#include "common/eigen.hpp"
#include "quadratures/quadratures.hpp"
#include "utils_hho.hpp"

namespace disk
{

namespace priv
{
// static condensation
template<typename T>
auto
static_condensation_impl(const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                         const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs,
                         const size_t                                                     num_cell_dofs,
                         const size_t                                                     num_faces_dofs)
{
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    const auto num_total_dofs = num_cell_dofs + num_faces_dofs;

    assert(lhs.rows() == lhs.cols());
    assert(lhs.cols() == num_total_dofs);

    if ((rhs.size() != num_cell_dofs) && (rhs.size() != num_total_dofs))
    {
        throw std::invalid_argument("static condensation: incorrect size of the rhs");
    }

    const matrix_type K_TT = lhs.topLeftCorner(num_cell_dofs, num_cell_dofs);
    const matrix_type K_TF = lhs.topRightCorner(num_cell_dofs, num_faces_dofs);
    const matrix_type K_FT = lhs.bottomLeftCorner(num_faces_dofs, num_cell_dofs);
    const matrix_type K_FF = lhs.bottomRightCorner(num_faces_dofs, num_faces_dofs);

    assert(K_TT.cols() == num_cell_dofs);
    assert(K_TT.cols() + K_TF.cols() == lhs.cols());
    assert(K_TT.rows() + K_FT.rows() == lhs.rows());
    assert(K_TF.rows() + K_FF.rows() == lhs.rows());
    assert(K_FT.cols() + K_FF.cols() == lhs.cols());

    const vector_type cell_rhs  = rhs.head(num_cell_dofs);
    vector_type       faces_rhs = vector_type::Zero(num_faces_dofs);

    if (rhs.size() == num_total_dofs)
    {
        faces_rhs = rhs.tail(num_faces_dofs);
    }

    const auto K_TT_ldlt = K_TT.ldlt();
    if (K_TT_ldlt.info() != Eigen::Success)
    {
        throw std::invalid_argument("static condensation: K_TT is not positive definite");
    }

    const matrix_type AL = K_TT_ldlt.solve(K_TF);
    const vector_type bL = K_TT_ldlt.solve(cell_rhs);

    const matrix_type AC = K_FF - K_FT * AL;
    const vector_type bC = faces_rhs - K_FT * bL;

    return std::make_tuple(std::make_pair(AC, bC), AL, bL);
}

// static decondensation for primal scalar problem
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
static_decondensation_impl(const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                           const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs,
                           const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              solF,
                           const size_t                                                     num_cell_dofs,
                           const size_t                                                     num_faces_dofs)
{
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    const auto num_total_dofs = num_cell_dofs + num_faces_dofs;

    assert(lhs.rows() == lhs.cols());
    assert(lhs.cols() == num_total_dofs);

    if ((rhs.size() < num_cell_dofs))
    {
        throw std::invalid_argument("static condensation: incorrect size of the rhs");
    }

    const matrix_type K_TT = lhs.topLeftCorner(num_cell_dofs, num_cell_dofs);
    const matrix_type K_TF = lhs.topRightCorner(num_cell_dofs, num_faces_dofs);

    const vector_type solT = K_TT.ldlt().solve(rhs.head(num_cell_dofs) - K_TF * solF);

    vector_type ret          = vector_type::Zero(num_total_dofs);
    ret.head(num_cell_dofs)  = solT;
    ret.tail(num_faces_dofs) = solF;

    return ret;
}
} // namespace priv

// static condensation for primal scalar problem like diffusion
template<typename Mesh, typename T>
auto
make_scalar_static_condensation_withMatrix(const Mesh&                                                      msh,
                                           const typename Mesh::cell_type&                                  cl,
                                           const hho_degree_info&                                           hdi,
                                           const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                                           const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs)
{
    const auto num_cell_dofs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    const auto num_faces     = howmany_faces(msh, cl);
    const auto num_face_dofs = scalar_basis_size(hdi.face_degree(), Mesh::dimension - 1);

    return priv::static_condensation_impl(lhs, rhs, num_cell_dofs, num_faces * num_face_dofs);
}

template<typename Mesh, typename T>
auto
make_scalar_static_condensation(const Mesh&                                                      msh,
                                const typename Mesh::cell_type&                                  cl,
                                const hho_degree_info&                                           hdi,
                                const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                                const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs)
{
    return std::get<0>(make_scalar_static_condensation_withMatrix(hdi, lhs, rhs));
}

// static decondensation for primal scalar problem
template<typename Mesh, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
make_scalar_static_decondensation(const Mesh&                                                      msh,
                                  const typename Mesh::cell_type&                                  cl,
                                  const hho_degree_info                                            hdi,
                                  const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                                  const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs,
                                  const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              solF)
{
    const auto num_cell_dofs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    const auto num_faces     = howmany_faces(msh, cl);
    const auto num_face_dofs = scalar_basis_size(hdi.face_degree(), Mesh::dimension - 1);

    return static_decondensation_impl(lhs, rhs, solF, num_cell_dofs, num_faces * num_face_dofs);
}

// static decondensation for primal scalar problem
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
make_scalar_static_decondensation_withMatrix(const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& AL,
                                             const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              bL,
                                             const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              solF)
{
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    vector_type ret       = vector_type::Zero(bL.size() + solF.size());
    ret.head(bL.size())   = bL - AL * solF;
    ret.tail(solF.size()) = solF;

    return ret;
}

} // end diskpp