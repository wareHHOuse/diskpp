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

/**
* @brief Compute the term \f$ (\nabla R^{k+1}_T(\hat{u}), \nabla R^{k+1}_T(\hat{v}) )_T   \f$  where \f$ R^{k+1}_T \f$
is the high-order
* reconstruction operator which if solution of
*
* \f$ (\nabla R^{k+1}_T(\hat{u}), \nabla w )_T = (\nabla u_T, \nabla w )_T + (u_{\partial T} - u_T, \nabla w
n_T)_{\partial T}, \quad \forall w \f$
*
*  The term \f$ (\nabla R^{k+1}_T(\hat{u}), \nabla R^{k+1}_T(\hat{v}) )_T   \f$ is the HHO conterpart of \f$ (\nabla u,
\nabla v )_T   \f$
*
* @tparam Mesh type of the mesh
* @param msh mesh
* @param cl  cell
* @param cell_infos  cell degree informations
* @return std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
dynamic_matrix<typename Mesh::coordinate_type>  > the first term is \f$ \nabla R^{k+1}_T \f$ and the second term is \f$
(\nabla R^{k+1}_T(\hat{u}), \nabla R^{k+1}_T(\hat{v}) )_T   \f$
*/
template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_scalar_hho_laplacian(const Mesh& msh, const typename Mesh::cell_type& cl, const CellDegreeInfo<Mesh>& cell_infos)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    const size_t DIM = Mesh::dimension;

    const auto cb = make_scalar_monomial_basis(msh, cl, cell_infos.reconstruction_degree());

    const auto rbs = scalar_basis_size(cell_infos.reconstruction_degree(), Mesh::dimension);
    const auto cbs = scalar_basis_size(cell_infos.cell_degree(), Mesh::dimension);

    const auto fcs_di = cell_infos.facesDegreeInfo();

    const auto num_faces_dofs = scalar_faces_dofs(msh, fcs_di);
    const auto num_total_dofs = cbs + num_faces_dofs;

    const matrix_type stiff  = make_stiffness_matrix(msh, cl, cb);
    matrix_type       gr_lhs = matrix_type::Zero(rbs - 1, rbs - 1);
    matrix_type       gr_rhs = matrix_type::Zero(rbs - 1, num_total_dofs);

    gr_lhs                           = stiff.block(1, 1, rbs - 1, rbs - 1);
    gr_rhs.block(0, 0, rbs - 1, cbs) = stiff.block(1, 0, rbs - 1, cbs);

    const auto fcs    = faces(msh, cl);
    size_t     offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc  = fcs[i];
        const auto fdi = fcs_di[i];

        if (fdi.hasUnknowns())
        {
            const auto face_degree = fdi.degree();
            const auto n           = normal(msh, cl, fc);
            const auto fb          = make_scalar_monomial_basis(msh, fc, face_degree);
            const auto fbs         = fb.size();

            auto qps_f = integrate(
              msh, fc, cell_infos.reconstruction_degree() - 1 + std::max(face_degree, cell_infos.cell_degree()));
            for (auto& qp : qps_f)
            {
                vector_type             c_phi_tmp  = cb.eval_functions(qp.point());
                vector_type             c_phi      = c_phi_tmp.head(cbs);
                Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic, DIM> c_dphi     = c_dphi_tmp.block(1, 0, rbs - 1, DIM);
                vector_type             f_phi      = fb.eval_functions(qp.point());
                gr_rhs.block(0, offset, rbs - 1, fbs) += qp.weight() * (c_dphi * n) * f_phi.transpose();
                gr_rhs.block(0, 0, rbs - 1, cbs) -= qp.weight() * (c_dphi * n) * c_phi.transpose();
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
make_scalar_hho_laplacian(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& hdi)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    const size_t DIM = Mesh::dimension;

    const auto rb = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
    const auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());

    const auto rbs = rb.size();
    const auto cbs = cb.size();
    const auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);

    const auto fcs = faces(msh, cl);

    const auto num_faces_dofs = fbs * fcs.size();
    const auto num_total_dofs = cbs + num_faces_dofs;

    const matrix_type stiff  = make_stiffness_matrix(msh, cl, rb);
    matrix_type       gr_lhs = matrix_type::Zero(rbs - 1, rbs - 1);
    matrix_type       gr_rhs = matrix_type::Zero(rbs - 1, num_total_dofs);

    gr_lhs                           = stiff.block(1, 1, rbs - 1, rbs - 1);
    gr_rhs.block(0, 0, rbs - 1, cbs) = stiff.block(1, 0, rbs - 1, cbs);

    size_t     offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc  = fcs[i];

        const auto face_degree = hdi.face_degree();
        const auto n           = normal(msh, cl, fc);
        const auto fb          = make_scalar_monomial_basis(msh, fc, face_degree);
        const auto fbs         = fb.size();

        auto rd = hdi.reconstruction_degree();
        auto cd = hdi.face_degree();
        auto fd = hdi.face_degree();

        auto qps_f = integrate(msh, fc, rd + std::max(cd, fd) );
        for (auto& qp : qps_f)
        {
            Matrix<T, Dynamic, DIM> r_dphi_tmp = rb.eval_gradients(qp.point());
            Matrix<T, Dynamic, DIM> r_dphi     = r_dphi_tmp.block(1, 0, rbs - 1, DIM);
            vector_type             c_phi      = cb.eval_functions(qp.point());
            vector_type             f_phi      = fb.eval_functions(qp.point());
            gr_rhs.block(0, offset, rbs - 1, fbs) += qp.weight() * (r_dphi * n) * f_phi.transpose();
            gr_rhs.block(0, 0, rbs - 1, cbs) -= qp.weight() * (r_dphi * n) * c_phi.transpose();
        }

        offset += fbs;
    }

    matrix_type oper = gr_lhs.ldlt().solve(gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}


/**
 * @brief Compute the term \f$ (\nabla R^{k+1}_T(\hat{u}), \nabla R^{k+1}_T(\hat{v}) )_T   \f$  where \f$ R^{k+1}_T \f$
is the high-order
 * reconstruction operator which if solution of
 *
 * \f$ (\nabla R^{k+1}_T(\hat{u}), \nabla w )_T = (\nabla u_T, \nabla w )_T + (u_{\partial T} - u_T, \nabla w
n_T)_{\partial T}, \quad \forall w \f$
 *
 *  The term \f$ (\nabla R^{k+1}_T(\hat{u}), \nabla R^{k+1}_T(\hat{v}) )_T   \f$ is the HHO conterpart of \f$ (\nabla u,
\nabla v )_T   \f$
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl  cell
 * @param di  hho degree information
 * @return std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
dynamic_matrix<typename Mesh::coordinate_type>  > the first term is \f$ \nabla R^{k+1}_T \f$ and the second term is \f$
(\nabla R^{k+1}_T(\hat{u}), \nabla R^{k+1}_T(\hat{v}) )_T   \f$
 */
#if 0
template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_scalar_hho_laplacian(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    const CellDegreeInfo<Mesh> cell_infos(msh, cl, di.cell_degree(), di.face_degree(), di.grad_degree());

    return make_scalar_hho_laplacian(msh, cl, cell_infos);
}
#endif

/**
 * @brief Compute the term \f$ (\nabla R^{k+1}_T(\hat{u}), \nabla R^{k+1}_T(\hat{v}) )_T   \f$  where \f$ R^{k+1}_T \f$
is the high-order
 * reconstruction operator which if solution of
 *
 * \f$ (\nabla R^{k+1}_T(\hat{u}), \nabla w )_T = (\nabla u_T, \nabla w )_T + (u_{\partial T} - u_T, \nabla w
n_T)_{\partial T}, \quad \forall w \f$
 *
 *  The term \f$ (\nabla R^{k+1}_T(\hat{u}), \nabla R^{k+1}_T(\hat{v}) )_T   \f$ is the HHO conterpart of \f$ (\nabla u,
\nabla v )_T   \f$
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl  cell
 * @param mesh_infos  hho degree informations
 * @return std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
dynamic_matrix<typename Mesh::coordinate_type>  > the first term is \f$ \nabla R^{k+1}_T \f$ and the second term is \f$
(\nabla R^{k+1}_T(\hat{u}), \nabla R^{k+1}_T(\hat{v}) )_T   \f$
 */
template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_scalar_hho_laplacian(const Mesh& msh, const typename Mesh::cell_type& cl, const MeshDegreeInfo<Mesh>& mesh_infos)
{
    return make_scalar_hho_laplacian(msh, cl, mesh_infos.cellDegreeInfo(msh, cl));
}

/**
 * @briefApply Dirichlet boundary conditions via Nitsche trick, version for diffusion.
 * @param rec_op the reconstruction operator
 * @param lhs the element lhs contribution
 * @param rhs the element rhs contribution
 *
 * lhs and rhs get updated in order to apply the boundary conditions.
 */

template<typename T>
using DM = Matrix<T, Dynamic, Dynamic>;

template<typename T>
using DV = Matrix<T, Dynamic, 1>;

template<typename Mesh, typename Function>
void
apply_dirichlet_via_nitsche(const Mesh& msh,
                            const typename Mesh::cell_type& cl,
                            DM<typename Mesh::coordinate_type>& rec_op,
                            DM<typename Mesh::coordinate_type>& lhs,
                            DV<typename Mesh::coordinate_type>& rhs,
                            const hho_degree_info& hdi,
                            const Function& bcs_fun,
                            typename Mesh::coordinate_type penalization = 1.0)
{
    const size_t DIM = Mesh::dimension;
    typedef typename Mesh::coordinate_type              T;
    typedef Matrix<T, Dynamic, Dynamic>       matrix_type;
    typedef Matrix<T, Dynamic, 1>             vector_type;
    typedef Matrix<T, Dynamic, DIM>           grad_type;

    auto rb = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
    auto cbs = scalar_basis_size(hdi.cell_degree(), DIM);
    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    for (size_t face_i = 0; face_i < fcs.size(); face_i++)
    {
        auto fc = fcs[face_i];

        if ( !msh.is_boundary(fc) )
            continue;

        auto fb = make_scalar_monomial_basis(msh, fc, hdi.face_degree());

        auto rbs = rb.size();
        auto fbs = fb.size();

        matrix_type mass = matrix_type::Zero(fbs, fbs);
        matrix_type trace = matrix_type::Zero(fbs, rbs-1);
        vector_type vec1 = vector_type::Zero(rbs-1);
        auto hf = measure(msh, fc);
        auto blkofs = cbs+fbs*face_i;
        auto qps = integrate(msh, fc, 2*hdi.reconstruction_degree());
        for (auto& qp : qps)
        {
            auto n = normal(msh, cl, fc);
            grad_type   r_dphi = rb.eval_gradients( qp.point() );
            vector_type r_dphi_n = r_dphi.block(1,0,rbs-1,DIM)*n;
            vector_type f_phi = fb.eval_functions( qp.point() );

            trace += qp.weight() * f_phi * r_dphi_n.transpose();
            mass += qp.weight() * f_phi * f_phi.transpose();
            vec1 += qp.weight() * r_dphi_n * bcs_fun(qp.point());

            rhs.block(blkofs, 0, fbs, 1) += qp.weight() * f_phi * (penalization/hf) * bcs_fun(qp.point());
        }

        matrix_type l = trace * rec_op;

        assert(l.cols() == lhs.cols());
        lhs.block(blkofs, 0, l.rows(), l.cols()) -= l;
        lhs.block(0, blkofs, l.cols(), l.rows()) -= l.transpose();
        lhs.block(blkofs, blkofs, mass.rows(), mass.cols()) += (penalization/hf)*mass;

        rhs -= rec_op.transpose() * vec1;
    }
}

} // end diskpp