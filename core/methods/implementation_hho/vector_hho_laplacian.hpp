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
size_t
nb_lag(const size_t dim)
{
    size_t lag = 1;
    if (dim == 3)
        lag = 3;
    return lag;
}
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_symmetric_laplacian(const Mesh&                     msh,
                                    const typename Mesh::cell_type& cl,
                                    const CellDegreeInfo<Mesh>&     cell_infos)
{
    using T        = typename Mesh::coordinate_type;
    const size_t N = Mesh::dimension;

    const auto recdeg = cell_infos.reconstruction_degree();
    const auto celdeg = cell_infos.cell_degree();

    const auto rb = make_vector_monomial_basis(msh, cl, recdeg);
    const auto cb = make_vector_monomial_basis(msh, cl, recdeg);

    const auto rbs = vector_basis_size(recdeg, N, N);
    const auto cbs = vector_basis_size(celdeg, N, N);

    const auto faces_infos    = cell_infos.facesDegreeInfo();
    const auto num_faces_dofs = vector_faces_dofs(msh, faces_infos);

    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, N, N>             gradient_type;
    typedef Matrix<T, Dynamic, N>       function_type;

    const size_t rbs_ho         = rbs - N;
    const size_t num_total_dofs = cbs + num_faces_dofs;
    const size_t nb_lag         = priv::nb_lag(N);

    matrix_type stiff  = matrix_type::Zero(rbs, rbs);
    matrix_type gr_lhs = matrix_type::Zero(rbs_ho + nb_lag, rbs_ho + nb_lag);
    matrix_type gr_rhs = matrix_type::Zero(rbs_ho + nb_lag, num_total_dofs);

    const auto qps = integrate(msh, cl, 2 * (recdeg - 1));
    for (auto& qp : qps)
    {
        const auto dphi    = rb.eval_sgradients(qp.point());
        const auto qp_dphi = priv::inner_product(qp.weight(), dphi);
        stiff += priv::outer_product(qp_dphi, dphi);
    }

    gr_lhs.block(0, 0, rbs_ho, rbs_ho) = stiff.block(N, N, rbs_ho, rbs_ho);
    gr_rhs.block(0, 0, rbs_ho, cbs)    = stiff.block(N, 0, rbs_ho, cbs);

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
            const auto fbs    = vector_basis_size(facdeg, N - 1, N);
            const auto fb     = make_vector_monomial_basis(msh, fc, facdeg);

            const auto qps_f = integrate(msh, fc, std::max(facdeg, celdeg) + recdeg - 1);
            for (auto& qp : qps_f)
            {
                eigen_compatible_stdvector<gradient_type> r_dphi_tmp = rb.eval_sgradients(qp.point());

                auto                                      begin_iter = std::next(r_dphi_tmp.begin(), N);
                eigen_compatible_stdvector<gradient_type> r_dphi(rbs_ho);
                std::copy(begin_iter, r_dphi_tmp.end(), r_dphi.begin());

                const function_type c_phi       = cb.eval_functions(qp.point());
                const function_type f_phi       = fb.eval_functions(qp.point());
                const function_type qp_r_dphi_n = qp.weight() * priv::inner_product(r_dphi, n);
                gr_rhs.block(0, offset, rbs_ho, fbs) += priv::outer_product(qp_r_dphi_n, f_phi);
                gr_rhs.block(0, 0, rbs_ho, cbs) -= priv::outer_product(qp_r_dphi_n, c_phi);
            }

            offset += fbs;
        }
    }

    const auto qps_2 = integrate(msh, cl, recdeg);

    matrix_type rot = matrix_type::Zero(rbs, nb_lag);
    for (auto& qp : qps_2)
    {
        const auto rphi = rb.eval_curls(qp.point());
        rot += qp.weight() * rphi;
    }
    gr_lhs.block(0, rbs_ho, rbs_ho, nb_lag) += rot.bottomLeftCorner(rbs_ho, nb_lag);
    gr_lhs.block(rbs_ho, 0, nb_lag, rbs_ho) += rot.bottomLeftCorner(rbs_ho, nb_lag).transpose();

    // use LU solver because lhs is only symmetric and positive
    matrix_type sol  = gr_lhs.lu().solve(gr_rhs);
    matrix_type oper = sol.block(0, 0, rbs_ho, num_total_dofs);
    matrix_type gr   = gr_rhs.block(0, 0, rbs_ho, num_total_dofs);
    matrix_type data = gr.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_symmetric_laplacian(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    const CellDegreeInfo<Mesh> cell_infos(msh, cl, di.cell_degree(), di.face_degree(), di.grad_degree());

    return make_vector_hho_symmetric_laplacian(msh, cl, cell_infos);
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_symmetric_laplacian(const Mesh&                     msh,
                                    const typename Mesh::cell_type& cl,
                                    const MeshDegreeInfo<Mesh>&     msh_infos)
{
    return make_vector_hho_symmetric_laplacian(msh, cl, msh_infos.cellDegreeInfo(msh, cl));
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_laplacian(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& hdi)
{
    const auto hho_scalar_laplacian = make_scalar_hho_laplacian(msh, cl, hdi);

    return make_vector_hho_laplacian(msh, cl, hdi, hho_scalar_laplacian);
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_laplacian(const Mesh&                                                      msh,
                          const typename Mesh::cell_type&                                  cl,
                          const hho_degree_info&                                           hdi,
                          const std::pair<dynamic_matrix<typename Mesh::coordinate_type>,
                                          dynamic_matrix<typename Mesh::coordinate_type>>& hho_scalar_laplacian)
{
    const auto oper = priv::compute_grad_vector(msh, cl, hdi, hho_scalar_laplacian.first);
    const auto data = priv::compute_lhs_vector(msh, cl, hdi, hho_scalar_laplacian.second);

    return std::make_pair(oper, data);
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_laplacian(const Mesh& msh, const typename Mesh::cell_type& cl, const MeshDegreeInfo<Mesh>& mesh_infos)
{
    const auto hho_scalar_laplacian = make_scalar_hho_laplacian(msh, cl, mesh_infos);

    return make_vector_hho_laplacian(msh, cl, mesh_infos.cellDegreeInfo(msh, cl), hho_scalar_laplacian);
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_laplacian(const Mesh&                                                      msh,
                          const typename Mesh::cell_type&                                  cl,
                          const MeshDegreeInfo<Mesh>&                                      mesh_infos,
                          const std::pair<dynamic_matrix<typename Mesh::coordinate_type>,
                                          dynamic_matrix<typename Mesh::coordinate_type>>& hho_scalar_laplacian)
{
    const auto oper =
      priv::compute_grad_vector(msh, cl, mesh_infos.cellDegreeInfo(msh, cl), hho_scalar_laplacian.first);
    const auto data =
      priv::compute_lhs_vector(msh, cl, mesh_infos.cellDegreeInfo(msh, cl), hho_scalar_laplacian.second);

    return std::make_pair(oper, data);
}

} // end diskpp