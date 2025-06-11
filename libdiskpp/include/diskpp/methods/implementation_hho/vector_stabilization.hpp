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
#include "scalar_stabilization.hpp"

namespace disk
{

/**
 * @brief compute the stabilization term \f$\sum_{F \in F_T} 1/h_F(u_F - u_T + \Pi^k_T R^{k+1}_T(\hat{u}_T) -
 * R^{k+1}_T(\hat{u}_T), v_F - v_T + \Pi^k_T R^{k+1}_T(\hat{v}_T) - R^{k+1}_T(\hat{v}_T))_F \f$ for vectorial HHO
 * unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param reconstruction reconstruction operator \f$ R^{k+1}_T \f$
 * @param di hho degree information
 * @param hF use diameter of face for scaling if true (or cell diameter if false)
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_vector_hho_stabilization(const Mesh&                                           msh,
                              const typename Mesh::cell_type&                       cl,
                              const dynamic_matrix<typename Mesh::coordinate_type>& reconstruction,
                              const CellDegreeInfo<Mesh>&                           cell_infos,
                              const bool hF = true)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto recdeg = cell_infos.reconstruction_degree();
    const auto celdeg = cell_infos.cell_degree();

    const auto rbs = vector_basis_size(recdeg, Mesh::dimension, Mesh::dimension);
    const auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);

    const size_t N = Mesh::dimension;

    const auto cb = make_vector_monomial_basis(msh, cl, recdeg);

    const matrix_type mass_mat = make_mass_matrix(msh, cl, cb);

    T h = 0;
    if (!hF)
    {
        h = diameter(msh, cl);
    }

    // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

    // Step 1: compute \pi_T^k p_T^k v (third term).
    const matrix_type M1    = mass_mat.block(0, 0, cbs, cbs);
    const matrix_type M2    = mass_mat.block(0, N, cbs, rbs - N);
    matrix_type       proj1 = -M1.ldlt().solve(M2 * reconstruction);

    assert(M2.cols() == reconstruction.rows());

    // Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
    proj1.block(0, 0, cbs, cbs) += matrix_type::Identity(cbs, cbs);

    const auto fcs = faces(msh, cl);

    const auto fcs_di         = cell_infos.facesDegreeInfo();
    const auto num_faces_dofs = vector_faces_dofs(msh, fcs_di);

    matrix_type data = matrix_type::Zero(cbs + num_faces_dofs, cbs + num_faces_dofs);

    size_t offset = cbs;
    // Step 3: project on faces (eqn. 21)
    for (size_t face_i = 0; face_i < fcs.size(); face_i++)
    {
        const auto fdi = fcs_di[face_i];

        if (fdi.hasUnknowns())
        {
            const auto fc = fcs[face_i];
            if(hF)
                h = diameter(msh, fc);

            const auto facdeg = fdi.degree();
            const auto fb     = make_vector_monomial_basis(msh, fc, facdeg);
            const auto fbs    = vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension);

            matrix_type face_mass_matrix  = make_mass_matrix(msh, fc, fb);
            matrix_type face_trace_matrix = matrix_type::Zero(fbs, rbs);

            const auto face_quadpoints = integrate(msh, fc, recdeg + facdeg);
            for (auto& qp : face_quadpoints)
            {
                const auto f_phi = fb.eval_functions(qp.point());
                const auto c_phi = cb.eval_functions(qp.point());
                face_trace_matrix += priv::outer_product(priv::inner_product(qp.weight(), f_phi), c_phi);
            }

            LLT<matrix_type> piKF;
            piKF.compute(face_mass_matrix);

            // Step 3a: \pi_F^k( v_F - p_T^k v )
            const matrix_type MR1 = face_trace_matrix.block(0, N, fbs, rbs - N);

            matrix_type proj2 = piKF.solve(MR1 * reconstruction);
            proj2.block(0, offset, fbs, fbs) -= matrix_type::Identity(fbs, fbs);

            // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
            const matrix_type MR2   = face_trace_matrix.block(0, 0, fbs, cbs);
            const matrix_type proj3 = piKF.solve(MR2 * proj1);
            const matrix_type BRF   = proj2 + proj3;

            data += BRF.transpose() * face_mass_matrix * BRF / h;

            offset += fbs;
        }
    }

    return data;
}

/**
 * @brief compute the stabilization term \f$\sum_{F \in F_T} 1/h_F(u_F - u_T + \Pi^k_T R^{k+1}_T(\hat{u}_T) -
 * R^{k+1}_T(\hat{u}_T), v_F - v_T + \Pi^k_T R^{k+1}_T(\hat{v}_T) - R^{k+1}_T(\hat{v}_T))_F \f$ for vectorial HHO
 * unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param reconstruction reconstruction operator \f$ R^{k+1}_T \f$
 * @param di hho degree information
 * @param hF use diameter of face for scaling if true (or cell diameter if false)
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_vector_hho_stabilization(const Mesh&                                           msh,
                              const typename Mesh::cell_type&                       cl,
                              const dynamic_matrix<typename Mesh::coordinate_type>& reconstruction,
                              const MeshDegreeInfo<Mesh>&                           msh_infos,
                              const bool                                            hF = true)
{
    return make_vector_hho_stabilization(msh, cl, reconstruction, msh_infos.cellDegreeInfo(msh, cl), hF);
}

/**
 * @brief compute the stabilization term \f$\sum_{F \in F_T} 1/h_F(u_F - u_T + \Pi^k_T R^{k+1}_T(\hat{u}_T) -
 * R^{k+1}_T(\hat{u}_T), v_F - v_T + \Pi^k_T R^{k+1}_T(\hat{v}_T) - R^{k+1}_T(\hat{v}_T))_F \f$ for vectorial HHO
 * unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param reconstruction reconstruction operator \f$ R^{k+1}_T \f$
 * @param di hho degree information
 * @param hF use diameter of face for scaling if true (or cell diameter if false)
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_vector_hho_stabilization(const Mesh&                                           msh,
                              const typename Mesh::cell_type&                       cl,
                              const dynamic_matrix<typename Mesh::coordinate_type>& reconstruction,
                              const hho_degree_info&                                di,
                              const bool                                            hF = true)
{
    const CellDegreeInfo<Mesh> cell_infos(
      msh, cl, di.cell_degree(), di.face_degree(), di.grad_degree());

    return make_vector_hho_stabilization(msh, cl, reconstruction, cell_infos, hF);
}

/**
 * @brief compute the stabilization term \f$ \sum_{F \in F_T} 1/h_F(u_F-\Pi^k_F(u_T), v_F-\Pi^k_F(v_T))_F \f$
 * for vectorial HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param di hho degree information
 * @param hF use diameter of face for scaling if true (or cell diameter if false)
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_vector_hdg_stabilization(const Mesh&                     msh,
                              const typename Mesh::cell_type& cl,
                              const hho_degree_info&          di,
                              bool                            hF = true,
                              bool                            scaling = true)
{
    const auto hdg_scalar_stab = make_scalar_hdg_stabilization(msh, cl, di, hF, scaling);

    return priv::compute_lhs_vector(msh, cl, di, hdg_scalar_stab);
}

/**
 * @brief compute the stabilization term \f$ \sum_{F \in F_T} 1/h_F(u_F-\Pi^k_F(u_T), v_F-\Pi^k_F(v_T))_F \f$
 * for vectorial HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param msh_infos mesh degree information
 * @param hF use diameter of face for scaling if true (or cell diameter if false)
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_vector_hdg_stabilization(const Mesh&                     msh,
                              const typename Mesh::cell_type& cl,
                              const MeshDegreeInfo<Mesh>&     msh_infos,
                              const bool                      hF = true,
                              bool                            scaling = true)
{
    const auto hdg_scalar_stab = make_scalar_hdg_stabilization(msh, cl, msh_infos, hF, scaling);

    return priv::compute_lhs_vector(msh, cl, msh_infos.cellDegreeInfo(msh, cl), hdg_scalar_stab);
}

/**
 * @brief compute the stabilization term \f$ \sum_{F \in F_T} 1/h_F(u_F-u_T, v_F-v_T)_F \f$
 * for vectorial HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param di hho degree information
 * @param hF use diameter of face for scaling if true (or cell diameter if false)
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_vector_dg_stabilization(const Mesh&                     msh,
                             const typename Mesh::cell_type& cl,
                             const hho_degree_info&          di,
                             const bool                      hF = true)
{
    const auto dg_scalar_stab = make_scalar_dg_stabilization(msh, cl, di, hF);

    return priv::compute_lhs_vector(msh, cl, di, dg_scalar_stab);
}

/**
 * @brief compute the stabilization term \f$ \sum_{F \in F_T} 1/h_F(u_F-u_T, v_F-v_T)_F \f$
 * for vectorial HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param msh_infos mesh degree information
 * @param hF use diameter of face for scaling if true (or cell diameter if false)
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_vector_dg_stabilization(const Mesh&                     msh,
                             const typename Mesh::cell_type& cl,
                             const MeshDegreeInfo<Mesh>&     msh_infos,
                             const bool                      hF = true)
{
    const auto dg_scalar_stab = make_scalar_dg_stabilization(msh, cl, msh_infos, hF);

    return priv::compute_lhs_vector(msh, cl, msh_infos.cellDegreeInfo(msh, cl), dg_scalar_stab);
}

// doesn't work for symmetric gradient
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_vector_hho_stabilization_optim(const Mesh&                                           msh,
                                    const typename Mesh::cell_type&                       cl,
                                    const dynamic_matrix<typename Mesh::coordinate_type>& reconstruction_scalar,
                                    const hho_degree_info&                                hdi,
                                    const bool                                            hF = true)
{
    const auto hho_scalar_stab = make_scalar_hho_stabilization(msh, cl, reconstruction_scalar, hdi, hF);

    return priv::compute_lhs_vector(msh, cl, hdi, hho_scalar_stab);
}

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_vector_hho_stabilization_optim(const Mesh&                                           msh,
                                    const typename Mesh::cell_type&                       cl,
                                    const dynamic_matrix<typename Mesh::coordinate_type>& reconstruction_scalar,
                                    const MeshDegreeInfo<Mesh>&                           msh_infos,
                                    const bool                                            hF = true)
{
    const auto hho_scalar_stab = make_scalar_hho_stabilization(msh, cl, reconstruction_scalar, msh_infos, hF);

    return priv::compute_lhs_vector(msh, cl, msh_infos.cellDegreeInfo(msh, cl), hho_scalar_stab);
}

} // end diskpp