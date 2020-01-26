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
#include "scalar_hho_condensation.hpp"

namespace disk
{

// static condensation for primal vectorial problem like elasticity
template<typename Mesh, typename T>
auto
make_vector_static_condensation_withMatrix(const Mesh&                                                      msh,
                                           const typename Mesh::cell_type&                                  cl,
                                           const hho_degree_info&                                           hdi,
                                           const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                                           const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs)
{
    const auto num_cell_dofs = vector_basis_size(hdi.cell_degree(), Mesh::dimension, Mesh::dimension);
    const auto num_faces     = howmany_faces(msh, cl);
    const auto num_face_dofs = vector_basis_size(hdi.face_degree(), Mesh::dimension - 1, Mesh::dimension);

    return priv::static_condensation_impl(lhs, rhs, num_cell_dofs, num_faces * num_face_dofs);
}

template<typename Mesh, typename T>
auto
make_vector_static_condensation_withMatrix(const Mesh&                                                      msh,
                                           const typename Mesh::cell_type&                                  cl,
                                           const MeshDegreeInfo<Mesh>&                                      msh_infos,
                                           const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                                           const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs,
                                           bool check_size = true)
{
    const auto cell_infos     = msh_infos.cellDegreeInfo(msh, cl);
    const auto num_cell_dofs  = vector_basis_size(cell_infos.cell_degree(), Mesh::dimension, Mesh::dimension);

    if(check_size)
    {
        const auto num_faces_dofs = vector_faces_dofs(msh, cell_infos.facesDegreeInfo());
        return priv::static_condensation_impl(lhs, rhs, num_cell_dofs, num_faces_dofs);
    }

    return priv::static_condensation_impl(lhs, rhs, num_cell_dofs, rhs.size()-num_cell_dofs);
}

template<typename Mesh, typename T>
auto
make_vector_static_condensation(const Mesh&                                                      msh,
                                const typename Mesh::cell_type&                                  cl,
                                const hho_degree_info&                                           hdi,
                                const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                                const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs)
{
    return std::get<0>(make_vector_static_condensation_withMatrix(msh, cl, hdi, lhs, rhs));
}

// static decondensation for primal vector problem
template<typename Mesh, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
make_vector_static_decondensation(const Mesh&                                                      msh,
                                  const typename Mesh::cell_type&                                  cl,
                                  const hho_degree_info                                            hdi,
                                  const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                                  const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs,
                                  const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              solF)
{
    const auto num_cell_dofs = vector_basis_size(hdi.cell_degree(), Mesh::dimension, Mesh::dimension);
    const auto num_faces     = howmany_faces(msh, cl);
    const auto num_face_dofs = vector_basis_size(hdi.face_degree(), Mesh::dimension - 1, Mesh::dimension);

    return static_decondensation_impl(lhs, rhs, solF, num_cell_dofs, num_faces * num_face_dofs);
}

// static decondensation for primal vector problem
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
make_vector_static_decondensation_withMatrix(const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& AL,
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