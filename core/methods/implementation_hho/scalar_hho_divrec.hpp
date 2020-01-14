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

/**
 * @brief compute the divergence reconstruction for vectorial HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param cell_infos cell degree informations
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the divergence term
 */
template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_hho_divergence_reconstruction(const Mesh&                     msh,
                                   const typename Mesh::cell_type& cl,
                                   const CellDegreeInfo<Mesh>&     cell_infos)
{
    using T = typename Mesh::coordinate_type;
    typedef dynamic_matrix<T> matrix_type;

    const auto cbas_s = make_scalar_monomial_basis(msh, cl, cell_infos.grad_degree());

    const auto dr_lhs = make_mass_matrix(msh, cl, cbas_s);
    const auto dr_rhs = make_hho_divergence_reconstruction_rhs(msh, cl, cell_infos);

    matrix_type oper = dr_lhs.ldlt().solve(dr_rhs);
    matrix_type data = dr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

/**
 * @brief compute the right hand-side of divergence reconstruction for vectorial HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param cell_infos cell degree informations
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the divergence term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_hho_divergence_reconstruction_rhs(const Mesh&                     msh,
                                       const typename Mesh::cell_type& cl,
                                       const CellDegreeInfo<Mesh>&     cell_infos)
{
    using T = typename Mesh::coordinate_type;
    typedef dynamic_matrix<T> matrix_type;

    const auto celdeg = cell_infos.cell_degree();
    const auto recdeg = cell_infos.grad_degree();

    const auto cbas_v = make_vector_monomial_basis(msh, cl, celdeg);
    const auto cbas_s = make_scalar_monomial_basis(msh, cl, recdeg);

    const auto rbs = scalar_basis_size(recdeg, Mesh::dimension);
    const auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);

    const auto faces_infos = cell_infos.facesDegreeInfo();

    const auto num_faces_dofs = vector_faces_dofs(msh, faces_infos);

    matrix_type dr_rhs = matrix_type::Zero(rbs, cbs + num_faces_dofs);

    if (recdeg > 0)
    {
        const auto qps = integrate(msh, cl, celdeg + recdeg - 1);
        for (auto& qp : qps)
        {
            const auto s_dphi = cbas_s.eval_gradients(qp.point());
            const auto v_phi  = cbas_v.eval_functions(qp.point());

            dr_rhs.block(0, 0, rbs, cbs) -= qp.weight() * priv::outer_product(s_dphi, v_phi);
        }
    }

    const auto fcs    = faces(msh, cl);
    size_t     offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fdi = faces_infos[i];

        if (fdi.hasUnknowns())
        {
            const auto fc     = fcs[i];
            const auto facdeg = fdi.degree();
            const auto fbs    = vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension);

            const auto n      = normal(msh, cl, fc);
            const auto fbas_v = make_vector_monomial_basis(msh, fc, facdeg);

            const auto qps_f = integrate(msh, fc, facdeg + recdeg);
            for (auto& qp : qps_f)
            {
                const auto s_phi = cbas_s.eval_functions(qp.point());
                const auto f_phi = fbas_v.eval_functions(qp.point());

                const auto qp_f_phi_n = priv::inner_product(f_phi, priv::inner_product(qp.weight(), n));
                dr_rhs.block(0, offset, rbs, fbs) += priv::outer_product(s_phi, qp_f_phi_n);
            }

            offset += fbs;
        }
    }

    return dr_rhs;
}

/**
 * @brief compute the divergence reconstruction for vectorial HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param di hho degree information
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the divergence term
 */
template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_hho_divergence_reconstruction(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    const CellDegreeInfo<Mesh> cell_infos(msh, cl, di.cell_degree(), di.face_degree(), di.grad_degree());

    return make_hho_divergence_reconstruction(msh, cl, cell_infos);
}

/**
 * @brief compute the divergence reconstruction for vectorial HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param cell_infos mesh degree information
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the divergence term
 */
template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_hho_divergence_reconstruction(const Mesh&                     msh,
                                   const typename Mesh::cell_type& cl,
                                   const MeshDegreeInfo<Mesh>&     msh_infos)
{
    return make_hho_divergence_reconstruction(msh, cl, msh_infos.cellDegreeInfo(msh, cl));
}

/**
 * @brief compute the divergence reconstruction for vectorial HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param di hho degree information
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the divergence term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_hho_divergence_reconstruction_rhs(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    const CellDegreeInfo<Mesh> cell_infos(msh, cl, di.cell_degree(), di.face_degree(), di.grad_degree());

    return make_hho_divergence_reconstruction_rhs(msh, cl, cell_infos);
}

/**
 * @brief compute the right hand-side of divergence reconstruction for vectorial HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param cell_infos mesh degree information
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the divergence term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_hho_divergence_reconstruction_rhs(const Mesh&                     msh,
                                       const typename Mesh::cell_type& cl,
                                       const MeshDegreeInfo<Mesh>&     msh_infos)
{
    return make_hho_divergence_reconstruction_rhs(msh, cl, msh_infos.cellDegreeInfo(msh, cl));
}

} // end diskpp