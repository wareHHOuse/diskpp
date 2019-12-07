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
template<typename Mesh, typename GradBasis>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_gradrec_impl(const Mesh&                     msh,
                             const typename Mesh::cell_type& cl,
                             const CellDegreeInfo<Mesh>&     cell_infos,
                             const GradBasis&                gb)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    const auto celdeg  = cell_infos.cell_degree();
    const auto graddeg = gb.degree();

    const auto cb = make_scalar_monomial_basis(msh, cl, celdeg);

    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    const auto gbs = gb.size();

    const auto fcs_di = cell_infos.facesDegreeInfo();

    const auto num_faces_dofs = scalar_faces_dofs(msh, fcs_di);

    const matrix_type gr_lhs = make_mass_matrix(msh, cl, gb);
    matrix_type       gr_rhs = matrix_type::Zero(gbs, cbs + num_faces_dofs);

    // (vT, div(tau))_T
    if (graddeg > 0)
    {
        const auto qps = integrate(msh, cl, celdeg + graddeg - 1);
        for (auto& qp : qps)
        {
            const auto        c_phi     = cb.eval_functions(qp.point());
            const auto        g_dphi    = gb.eval_divergences(qp.point());
            const vector_type qp_g_dphi = qp.weight() * g_dphi;

            gr_rhs.block(0, 0, gbs, cbs) -= priv::outer_product(qp_g_dphi, c_phi);
        }
    }

    // (vF, tau.n)_F
    const auto fcs    = faces(msh, cl);
    size_t     offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fdi = fcs_di[i];
        if (fdi.hasUnknowns())
        {
            const auto fc     = fcs[i];
            const auto facdeg = fdi.degree();
            const auto n      = normal(msh, cl, fc);
            const auto fbs    = scalar_basis_size(facdeg, Mesh::dimension - 1);
            const auto fb     = make_scalar_monomial_basis(msh, fc, facdeg);

            const auto qps_f = integrate(msh, fc, graddeg + facdeg);
            for (auto& qp : qps_f)
            {
                const vector_type f_phi      = fb.eval_functions(qp.point());
                const auto        g_phi      = gb.eval_functions(qp.point());
                const vector_type qp_g_phi_n = g_phi * (qp.weight() * n);

                gr_rhs.block(0, offset, gbs, fbs) += priv::outer_product(qp_g_phi_n, f_phi);
            }
            offset += fbs;
        }
    }

    matrix_type oper = gr_lhs.ldlt().solve(gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}
}

/**
 * @brief Compute the term \f$ (G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$  where \f$ G^{k}_T \f$ is the hho gradient
reconstruction which if solution of
 *
 * \f$ (G^{k}_T(\hat{u})(\hat{u}), \tau )_T = (\nabla u_T, \tau )_T + (u_{\partial T} - u_T, \tau n_T)_{\partial T},
\quad \forall w \f$
 *
 *  The term \f$ (G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$ is the HHO conterpart of \f$ (\nabla u, \nabla v )_T \f$
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl  cell
 * @param di  hho degree information
 * @return std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
dynamic_matrix<typename Mesh::coordinate_type>  > the first term is \f$ G^{k}_T \f$ and the second term is \f$
(G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$
 */
template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_gradrec(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    const auto graddeg = di.grad_degree();
    const auto gb      = make_vector_monomial_basis(msh, cl, graddeg);

    const CellDegreeInfo<Mesh> cell_infos(msh, cl, di.cell_degree(), di.face_degree(), di.grad_degree());

    return priv::make_vector_hho_gradrec_impl(msh, cl, cell_infos, gb);
}

/**
 * @brief Compute the term \f$ (G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$  where \f$ G^{k}_T \f$ is the hho gradient
reconstruction which if solution of
 *
 * \f$ (G^{k}_T(\hat{u})(\hat{u}), \tau )_T = (\nabla u_T, \tau )_T + (u_{\partial T} - u_T, \tau n_T)_{\partial T},
\quad \forall w \f$
 *
 *  The term \f$ (G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$ is the HHO conterpart of \f$ (\nabla u, \nabla v )_T \f$
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl  cell
 * @param msh_infos  mesh degree informations
 * @return std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
dynamic_matrix<typename Mesh::coordinate_type>  > the first term is \f$ G^{k}_T \f$ and the second term is \f$
(G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$
 */
template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_gradrec(const Mesh& msh, const typename Mesh::cell_type& cl, const MeshDegreeInfo<Mesh>& msh_infos)
{
    const auto cell_infos = msh_infos.cellDegreeInfo(msh, cl);
    const auto gb         = make_vector_monomial_basis(msh, cl, cell_infos.grad_degree());

    return priv::make_vector_hho_gradrec_impl(msh, cl, cell_infos, gb);
}

/**
 * @brief Compute the term \f$ (G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$  where \f$ G^{k}_T \f$ is the hho gradient
reconstruction which if solution of
 *
 * \f$ (G^{k}_T(\hat{u})(\hat{u}), \tau )_T = (\nabla u_T, \tau )_T + (u_{\partial T} - u_T, \tau n_T)_{\partial T},
\quad \forall w \f$
 *
 *  The term \f$ (G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$ is the HHO conterpart of \f$ (\nabla u, \nabla v )_T \f$
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl  cell
 * @param cell_infos  cell degree informations
 * @return std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
dynamic_matrix<typename Mesh::coordinate_type>  > the first term is \f$ G^{k}_T \f$ and the second term is \f$
(G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$
 */
template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_gradrec(const Mesh& msh, const typename Mesh::cell_type& cl, const CellDegreeInfo<Mesh>& cell_infos)
{
    const auto gb = make_vector_monomial_basis(msh, cl, cell_infos.grad_degree());

    return priv::make_vector_hho_gradrec_impl(msh, cl, cell_infos, gb);
}

/**
 * @brief Compute the term \f$ (K G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$  where \f$ G^{k}_T \f$ is the hho
gradient reconstruction which if solution of
 *
 * \f$ (G^{k}_T(\hat{u})(\hat{u}), \tau )_T = (\nabla u_T, \tau )_T + (u_{\partial T} - u_T, \tau n_T)_{\partial T},
\quad \forall w \f$
 *
 *  The term \f$ (K G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$ is the HHO conterpart of \f$ (K \nabla u, \nabla v )_T
\f$
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl  cell
 * @param di  hho degree informations
 * @param mfield is a material fiel which can depend on the coordiantes (mfield is symmetric positive definite matrix)
 * @return std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
dynamic_matrix<typename Mesh::coordinate_type>  > the first term is \f$ G^{k}_T \f$ and the second term is \f$ (K
G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$
 */
template<typename Mesh, typename TensorField>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_gradrec(const Mesh&                     msh,
                        const typename Mesh::cell_type& cl,
                        const hho_degree_info&          di,
                        const TensorField&              mfield)
{
    const auto gradrec = make_vector_hho_gradrec(msh, cl, di);

    const auto graddeg = di.grad_degree();
    const auto gb      = make_vector_monomial_basis(msh, cl, graddeg);

    const auto mass = make_mass_matrix(msh, cl, gb, mfield);
    const auto LHS  = gradrec.first.transpose() * mass * gradrec.first;

    return std::make_pair(gradrec.first, LHS);
}

/**
 * @brief Compute the term \f$ (K G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$  where \f$ G^{k}_T \f$ is the hho
gradient reconstruction which if solution of
 *
 * \f$ (G^{k}_T(\hat{u})(\hat{u}), \tau )_T = (\nabla u_T, \tau )_T + (u_{\partial T} - u_T, \tau n_T)_{\partial T},
\quad \forall w \f$
 *
 *  The term \f$ (K G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$ is the HHO conterpart of \f$ (K \nabla u, \nabla v )_T
\f$
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl  cell
 * @param msh_infos  mesh degree informations
 * @param mfield is a material fiel which can depend on the coordiantes (mfield is symmetric positive definite matrix)
 * @return std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
dynamic_matrix<typename Mesh::coordinate_type>  > the first term is \f$ G^{k}_T \f$ and the second term is \f$ (K
G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$
 */
template<typename Mesh, typename TensorField>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_gradrec(const Mesh&                     msh,
                        const typename Mesh::cell_type& cl,
                        const MeshDegreeInfo<Mesh>&     msh_infos,
                        const TensorField&              mfield)
{
    const auto cell_infos = msh_infos.cellDegreeInfo(msh, cl);
    const auto gradrec    = make_vector_hho_gradrec(msh, cl, msh_infos);

    const auto graddeg = cell_infos.grad_degree();
    const auto gb      = make_vector_monomial_basis(msh, cl, graddeg);

    const auto mass = make_mass_matrix(msh, cl, gb, mfield);
    const auto LHS  = gradrec.first.transpose() * mass * gradrec.first;

    return std::make_pair(gradrec.first, LHS);
}

/**
 * @brief Compute the term \f$ (G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$  where \f$ G^{k}_T \f$ is the hho gradient
reconstruction which if solution of
 *
 * \f$ (G^{k}_T(\hat{u})(\hat{u}), \tau )_T = (\nabla u_T, \tau )_T + (u_{\partial T} - u_T, \tau n_T)_{\partial T},
\quad \forall w \f$
 *
 *  The term \f$ (G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$ is the HHO conterpart of \f$ (K \nabla u, \nabla v )_T
\f$
 *
 *  Here the reconstructed gradient \f$ G^{k}_T \f$ is reconstructed in the Raviart-Thomas space (works only for
simplical cell)
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl  simplicial cell (triangle or tetrahedron)
 * @param di  hho degree information
 * @return std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
dynamic_matrix<typename Mesh::coordinate_type>  > the first term is \f$ G^{k}_T \f$ and the second term is \f$ (
G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$
 */

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_gradrec_RT(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    const auto graddeg = di.grad_degree();
    const auto gb      = make_vector_monomial_basis_RT(msh, cl, graddeg);

    const CellDegreeInfo<Mesh> cell_infos(msh, cl, di.cell_degree(), di.face_degree(), di.grad_degree());

    return make_vector_hho_gradrec_impl(msh, cl, cell_infos, gb);
}

/**
 * @brief Compute the term \f$ (G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$  where \f$ G^{k}_T \f$ is the hho gradient
reconstruction which if solution of
 *
 * \f$ (G^{k}_T(\hat{u})(\hat{u}), \tau )_T = (\nabla u_T, \tau )_T + (u_{\partial T} - u_T, \tau n_T)_{\partial T},
\quad \forall w \f$
 *
 *  The term \f$ (G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$ is the HHO conterpart of \f$ (K \nabla u, \nabla v )_T
\f$
 *
 *  Here the reconstructed gradient \f$ G^{k}_T \f$ is reconstructed in the Raviart-Thomas space (works only for
simplical cell)
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl  simplicial cell (triangle or tetrahedron)
 * @param msh_infos  mesh degree information
 * @return std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
dynamic_matrix<typename Mesh::coordinate_type>  > the first term is \f$ G^{k}_T \f$ and the second term is \f$ (
G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$
 */

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_gradrec_RT(const Mesh& msh, const typename Mesh::cell_type& cl, const MeshDegreeInfo<Mesh>& msh_infos)
{
    const auto cell_infos = msh_infos.cellDegreeInfo(msh, cl);
    const auto gb         = make_vector_monomial_basis(msh, cl, cell_infos.grad_degree());

    return make_vector_hho_gradrec_impl(msh, cl, cell_infos, gb);
}

/**
 * @brief Compute the term \f$ (K G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$  where \f$ G^{k}_T \f$ is the hho
gradient reconstruction which if solution of
 *
 * \f$ (G^{k}_T(\hat{u})(\hat{u}), \tau )_T = (\nabla u_T, \tau )_T + (u_{\partial T} - u_T, \tau n_T)_{\partial T},
\quad \forall w \f$
 *
 *  The term \f$ (K G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$ is the HHO conterpart of \f$ (K \nabla u, \nabla v )_T
\f$
 *
 *  Here the reconstructed gradient \f$ G^{k}_T \f$ is reconstructed in the Raviart-Thomas space (works only for
simplical cell)
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl  simplicial cell (triangle or tetrahedron)
 * @param di  hho degree information
 * @param mfield is a material fiel which can depend on the coordiantes (mfield is symmetric positive definite matrix)
 * @return std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
dynamic_matrix<typename Mesh::coordinate_type>  > the first term is \f$ G^{k}_T \f$ and the second term is \f$ (K
G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$
 */
template<typename Mesh, typename TensorField>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_gradrec_RT(const Mesh&                     msh,
                           const typename Mesh::cell_type& cl,
                           const hho_degree_info&          di,
                           const TensorField&              mfield)
{
    const auto gradrec_RT = make_vector_hho_gradrec_RT(msh, cl, di);

    const auto graddeg = di.grad_degree();
    const auto gb      = make_vector_monomial_basis_RT(msh, cl, graddeg);

    const auto mass = make_mass_matrix(msh, cl, gb, mfield);
    const auto LHS  = gradrec_RT.first.transpose() * mass * gradrec_RT.first;

    return std::make_pair(gradrec_RT.first, LHS);
}

/**
 * @brief Compute the term \f$ (K G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$  where \f$ G^{k}_T \f$ is the hho
gradient reconstruction which if solution of
 *
 * \f$ (G^{k}_T(\hat{u})(\hat{u}), \tau )_T = (\nabla u_T, \tau )_T + (u_{\partial T} - u_T, \tau n_T)_{\partial T},
\quad \forall w \f$
 *
 *  The term \f$ (K G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$ is the HHO conterpart of \f$ (K \nabla u, \nabla v )_T
\f$
 *
 *  Here the reconstructed gradient \f$ G^{k}_T \f$ is reconstructed in the Raviart-Thomas space (works only for
simplical cell)
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl  simplicial cell (triangle or tetrahedron)
 * @param msh_infos  mesh degree informations
 * @param mfield is a material fiel which can depend on the coordiantes (mfield is symmetric positive definite matrix)
 * @return std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
dynamic_matrix<typename Mesh::coordinate_type>  > the first term is \f$ G^{k}_T \f$ and the second term is \f$ (K
G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$
 */
template<typename Mesh, typename TensorField>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_gradrec_RT(const Mesh&                     msh,
                           const typename Mesh::cell_type& cl,
                           const MeshDegreeInfo<Mesh>&     msh_infos,
                           const TensorField&              mfield)
{
    const auto gradrec_RT = make_vector_hho_gradrec_RT(msh, cl, msh_infos);

    const auto graddeg = msh_infos.cellDegreeInfo(msh, cl).grad_degree();
    const auto gb      = make_vector_monomial_basis_RT(msh, cl, graddeg);

    const auto mass = make_mass_matrix(msh, cl, gb, mfield);
    const auto LHS  = gradrec_RT.first.transpose() * mass * gradrec_RT.first;

    return std::make_pair(gradrec_RT.first, LHS);
}

} // end diskpp