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

#include "diskpp/adaptivity/adaptivity.hpp"
#include "diskpp/bases/bases.hpp"
#include "diskpp/boundary_conditions/boundary_conditions.hpp"
#include "diskpp/common/eigen.hpp"
#include "diskpp/quadratures/quadratures.hpp"
#include "utils_hho.hpp"

namespace disk
{

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_matrix_symmetric_gradrec(const Mesh&                     msh,
                              const typename Mesh::cell_type& cl,
                              const CellDegreeInfo<Mesh>&     cell_infos)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const size_t N = Mesh::dimension;

    const auto graddeg = cell_infos.grad_degree();
    const auto celdeg  = cell_infos.cell_degree();

    const auto gb = make_sym_matrix_monomial_basis(msh, cl, graddeg);
    const auto cb = make_vector_monomial_basis(msh, cl, celdeg);

    const auto gbs = sym_matrix_basis_size(graddeg, Mesh::dimension, Mesh::dimension);
    const auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);

    const auto faces_infos    = cell_infos.facesDegreeInfo();
    const auto num_faces_dofs = vector_faces_dofs(msh, faces_infos);

    matrix_type gr_lhs = matrix_type::Zero(gbs, gbs);
    matrix_type gr_rhs = matrix_type::Zero(gbs, cbs + num_faces_dofs);

    // this is very costly to build it
    const auto qps = integrate(msh, cl, 2 * graddeg);

    size_t dec = 0;
    if (N == 3)
        dec = 6;
    else if (N == 2)
        dec = 3;
    else
        std::logic_error("Expected 3 >= dim > 1");

    for (auto& qp : qps)
    {
        const auto gphi = gb.eval_functions(qp.point());

        for (size_t j = 0; j < gbs; j++)
        {
            const auto qp_gphi_j = priv::inner_product(qp.weight(), gphi[j]);
            for (size_t i = j; i < gbs; i += dec)
                gr_lhs(i, j) += priv::inner_product(gphi[i], qp_gphi_j);
        }
    }

    // upper part
    for (size_t j = 0; j < gbs; j++)
        for (size_t i = 0; i < j; i++)
            gr_lhs(i, j) = gr_lhs(j, i);

    // compute rhs
    if (celdeg > 0)
    {
        const auto qpc = integrate(msh, cl, graddeg + celdeg - 1);
        for (auto& qp : qpc)
        {
            const auto gphi    = gb.eval_functions(qp.point());
            const auto dphi    = cb.eval_sgradients(qp.point());
            const auto qp_dphi = priv::inner_product(qp.weight(), dphi);

            gr_rhs.block(0, 0, gbs, cbs) += priv::outer_product(gphi, qp_dphi);

        } // end qp
    }

    const auto fcs    = faces(msh, cl);
    size_t     offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fdi = faces_infos[i];
        if (fdi.hasUnknowns())
        {
            const auto fc = fcs[i];
            const auto n  = normal(msh, cl, fc);

            const auto facdeg = fdi.degree();
            const auto fb     = make_vector_monomial_basis(msh, fc, facdeg);
            const auto fbs    = vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension);

            const auto qps_f = integrate(msh, fc, graddeg + std::max(celdeg, facdeg));
            for (auto& qp : qps_f)
            {
                const auto cphi = cb.eval_functions(qp.point());
                const auto gphi = gb.eval_functions(qp.point());
                const auto fphi = fb.eval_functions(qp.point());

                const auto qp_gphi_n = priv::inner_product(gphi, priv::inner_product(qp.weight(), n));
                gr_rhs.block(0, offset, gbs, fbs) += priv::outer_product(qp_gphi_n, fphi);
                gr_rhs.block(0, 0, gbs, cbs) -= priv::outer_product(qp_gphi_n, cphi);
            }

            offset += fbs;
        }
    }

    matrix_type oper = gr_lhs.ldlt().solve(gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_matrix_symmetric_gradrec(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    const CellDegreeInfo<Mesh> cell_infos(msh, cl, di.cell_degree(), di.face_degree(), di.grad_degree());

    return make_matrix_symmetric_gradrec(msh, cl, cell_infos);
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_matrix_symmetric_gradrec(const Mesh&                     msh,
                              const typename Mesh::cell_type& cl,
                              const MeshDegreeInfo<Mesh>&     msh_infos)
{
    return make_matrix_symmetric_gradrec(msh, cl, msh_infos.cellDegreeInfo(msh, cl));
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_matrix_symmetric_gradrec_RT(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const size_t N = Mesh::dimension;

    const auto graddeg = di.grad_degree();
    const auto celdeg  = di.cell_degree();
    const auto facdeg  = di.face_degree();

    const auto gb = make_sym_matrix_monomial_basis_RT(msh, cl, graddeg);
    const auto cb = make_vector_monomial_basis(msh, cl, celdeg);

    const auto gbs = sym_matrix_basis_size_RT(graddeg, Mesh::dimension, Mesh::dimension);
    const auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
    const auto fbs = vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension);

    const auto num_faces = howmany_faces(msh, cl);

    const matrix_type gr_lhs = make_mass_matrix(msh, cl, gb);
    matrix_type       gr_rhs = matrix_type::Zero(gbs, cbs + num_faces * fbs);

    // compute rhs
    if (celdeg > 0)
    {
        const auto qpc = integrate(msh, cl, graddeg + celdeg - 1);
        for (auto& qp : qpc)
        {
            const auto gphi    = gb.eval_functions(qp.point());
            const auto dphi    = cb.eval_sgradients(qp.point());
            const auto qp_dphi = priv::inner_product(qp.weight(), dphi);

            gr_rhs.block(0, 0, gbs, cbs) += priv::outer_product(gphi, qp_dphi);

        } // end qp
    }

    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = normal(msh, cl, fc);
        const auto fb = make_vector_monomial_basis(msh, fc, facdeg);

        const auto qps_f = integrate(msh, fc, graddeg + std::max(celdeg, facdeg));
        for (auto& qp : qps_f)
        {
            const auto cphi = cb.eval_functions(qp.point());
            const auto gphi = gb.eval_functions(qp.point());
            const auto fphi = fb.eval_functions(qp.point());

            const auto qp_gphi_n = priv::inner_product(gphi, priv::inner_product(qp.weight(), n));
            gr_rhs.block(0, cbs + i * fbs, gbs, fbs) += priv::outer_product(qp_gphi_n, fphi);
            gr_rhs.block(0, 0, gbs, cbs) -= priv::outer_product(qp_gphi_n, cphi);
        }
    }

    matrix_type oper = gr_lhs.ldlt().solve(gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_marix_hho_gradrec(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& hdi)

{
    const auto hho_gradrec_vector = make_vector_hho_gradrec(msh, cl, hdi);
    const auto lhs                = priv::compute_lhs_vector(msh, cl, hdi, hho_gradrec_vector.second);
    const auto oper               = priv::compute_grad_matrix(msh, cl, hdi, hho_gradrec_vector.first);

    return std::make_pair(oper, lhs);
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_marix_hho_gradrec(const Mesh& msh, const typename Mesh::cell_type& cl, const MeshDegreeInfo<Mesh>& msh_infos)

{
    const auto cell_infos = msh_infos.cellDegreeInfo(msh, cl);

    const auto hho_gradrec_vector = make_vector_hho_gradrec(msh, cl, cell_infos);
    const auto lhs                = priv::compute_lhs_vector(msh, cl, cell_infos, hho_gradrec_vector.second);
    const auto oper               = priv::compute_grad_matrix(msh, cl, cell_infos, hho_gradrec_vector.first);

    return std::make_pair(oper, lhs);
}

} // end diskpp