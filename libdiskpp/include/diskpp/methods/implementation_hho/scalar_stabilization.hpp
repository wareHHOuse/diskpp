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

#include <algorithm>
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
 * @brief compute the difference operator \f$ \Delta := u_{\partial T}-u_T \f$
 * for scalar HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param cell_infos cell degree informations
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the differenc operator
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_hho_difference(const Mesh& msh, const typename Mesh::cell_type& cl, const CellDegreeInfo<Mesh>& cell_infos)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto celdeg = cell_infos.cell_degree();

    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);

    const auto faces_infos    = cell_infos.facesDegreeInfo();
    const auto num_faces_dofs = scalar_faces_dofs(msh, faces_infos);
    const auto total_dofs     = cbs + num_faces_dofs;
    const auto num_diff_dofs  = scalar_diff_dofs(msh, cell_infos);

    matrix_type data = matrix_type::Zero(num_diff_dofs, total_dofs);

    const auto cb = make_scalar_monomial_basis(msh, cl, celdeg);

    const auto fcs          = faces(msh, cl);
    size_t     offset_faces = cbs;
    size_t     offset_diff  = 0;

    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fdi = faces_infos[i];

        if (fdi.hasUnknowns())
        {
            const auto fc     = fcs[i];
            const auto facdeg = fdi.degree();
            const auto ddeg   = std::max(facdeg, celdeg);
            const auto db     = make_scalar_monomial_basis(msh, fc, ddeg);
            const auto fbs    = scalar_basis_size(facdeg, Mesh::dimension - 1);
            const auto dbs    = scalar_basis_size(ddeg, Mesh::dimension - 1);

            matrix_type mass_F = matrix_type::Zero(dbs, dbs);
            matrix_type trace  = matrix_type::Zero(dbs, cbs);

            const auto qps = integrate(msh, fc, 2 * ddeg);
            for (auto& qp : qps)
            {
                const auto c_phi = cb.eval_functions(qp.point());
                const auto d_phi = db.eval_functions(qp.point());

                mass_F += priv::outer_product(priv::inner_product(qp.weight(), d_phi), d_phi);
                trace += priv::outer_product(priv::inner_product(qp.weight(), d_phi), c_phi);
            }

            const matrix_type proj = mass_F.ldlt().solve(trace);

            data.block(offset_diff, 0, dbs, cbs)            = -proj;
            data.block(offset_diff, offset_faces, fbs, fbs) = matrix_type::Identity(fbs, fbs);

            offset_faces += fbs;
            offset_diff += dbs;
        }
    }

    return data;
}

/**
 * @brief compute the difference operator \f$ \Delta := u_{\partial T}-u_T \f$
 * for scalar HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param di hho degree information
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the difference operator
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_hho_difference(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    const CellDegreeInfo<Mesh> cell_infos(msh, cl, di.cell_degree(), di.face_degree(), di.grad_degree());

    return make_scalar_hho_difference(msh, cl, cell_infos);
}

/**
 * @brief compute the stabilization term \f$ \sum_{F \in F_T} 1/h_F(u_F-\Pi^k_F(u_T), v_F-\Pi^k_F(v_T))_F \f$
 * for scalar HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param cell_infos cell degree informations
 * @param hF use diameter of face for scaling if true (or cell diameter if false)
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_hdg_stabilization(const Mesh&                     msh,
                              const typename Mesh::cell_type& cl,
                              const CellDegreeInfo<Mesh>&     cell_infos,
                              bool hF = true,
                              bool scaling = true)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto celdeg = cell_infos.cell_degree();

    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);

    const auto faces_infos    = cell_infos.facesDegreeInfo();
    const auto num_faces_dofs = scalar_faces_dofs(msh, faces_infos);
    const auto total_dofs     = cbs + num_faces_dofs;

    matrix_type data = matrix_type::Zero(total_dofs, total_dofs);

    T h = 0;
    if(!hF){
        h = diameter(msh, cl);
    }

    const auto cb     = make_scalar_monomial_basis(msh, cl, celdeg);
    const auto fcs    = faces(msh, cl);
    size_t     offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fdi = faces_infos[i];

        if (fdi.hasUnknowns())
        {
            const auto fc     = fcs[i];
            const auto facdeg = fdi.degree();

            if (hF) {
                h = diameter(msh, fc);
            }

            const auto fb  = make_scalar_monomial_basis(msh, fc, facdeg);
            const auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);

            const matrix_type If    = matrix_type::Identity(fbs, fbs);
            matrix_type       oper  = matrix_type::Zero(fbs, total_dofs);
            matrix_type       tr    = matrix_type::Zero(fbs, total_dofs);
            matrix_type       mass  = make_mass_matrix(msh, fc, fb);
            matrix_type       trace = matrix_type::Zero(fbs, cbs);

            oper.block(0, offset, fbs, fbs) = -If;

            const auto qps = integrate(msh, fc, facdeg + celdeg);
            for (auto& qp : qps)
            {
                const auto c_phi = cb.eval_functions(qp.point());
                const auto f_phi = fb.eval_functions(qp.point());

                assert(c_phi.rows() == cbs);
                assert(f_phi.rows() == fbs);
                assert(c_phi.cols() == f_phi.cols());

                trace += priv::outer_product(priv::inner_product(qp.weight(), f_phi), c_phi);
            }

            tr.block(0, offset, fbs, fbs) = -mass;
            tr.block(0, 0, fbs, cbs)      = trace;

            oper.block(0, 0, fbs, cbs) = mass.ldlt().solve(trace);
            if (scaling) {
                data += oper.transpose() * tr * (1. / h);
            }
            else {
                data += oper.transpose() * tr;
            }

            offset += fbs;
        }
    }

    return data;
}

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_hdg_stabilization2(const Mesh&                      msh,
    const typename Mesh::cell_type&                         cl,
    const hho_degree_info&                                  hdi,
    const static_matrix<typename Mesh::coordinate_type, Mesh::dimension, Mesh::dimension>& kappaT)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto celdeg = hdi.cell_degree();
    const auto cb = make_scalar_monomial_basis(msh, cl, celdeg);
    const auto cbs = cb.size();

    const auto fcs = faces(msh, cl);
    const auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    const auto num_faces_dofs = fbs*fcs.size();
    const auto total_dofs     = cbs + num_faces_dofs;

    matrix_type data = matrix_type::Zero(total_dofs, total_dofs);

    T h = diameter(msh, cl);

    size_t offset   = cbs;

    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto facdeg = hdi.face_degree();
        const auto fb  = make_scalar_monomial_basis(msh, fc, facdeg);
        const auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);

        const matrix_type If    = matrix_type::Identity(fbs, fbs);
        matrix_type       oper  = matrix_type::Zero(fbs, total_dofs);
        matrix_type       tr    = matrix_type::Zero(fbs, total_dofs);
        matrix_type       mass  = make_mass_matrix(msh, fc, fb);
        matrix_type       trace = matrix_type::Zero(fbs, cbs);

        oper.block(0, offset, fbs, fbs) = -If;

        const auto qps = integrate(msh, fc, facdeg + celdeg);
        for (auto& qp : qps)
        {
            const auto c_phi = cb.eval_functions(qp.point());
            const auto f_phi = fb.eval_functions(qp.point());

            assert(c_phi.rows() == cbs);
            assert(f_phi.rows() == fbs);
            assert(c_phi.cols() == f_phi.cols());

            trace += priv::outer_product(priv::inner_product(qp.weight(), f_phi), c_phi);
        }

        tr.block(0, offset, fbs, fbs) = -mass;
        tr.block(0, 0, fbs, cbs)      = trace;

        auto n = normal(msh, cl, fc);
        T stabparam = n.dot(kappaT*n)/h;
        oper.block(0, 0, fbs, cbs) = mass.ldlt().solve(trace);
        data += oper.transpose() * tr * stabparam;

        offset += fbs;

    }

    return data;
}

/**
 * @brief compute the stabilization term \f$ \sum_{F \in F_T} 1/h_F(u_F-\Pi^k_F(u_T), v_F-\Pi^k_F(v_T))_F \f$
 * for scalar HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param cell_infos cell degree informations
 * @param hF use diameter of face for scaling if true (or cell diameter if false)
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_hdg_stabilization_diff(const Mesh&                     msh,
                                   const typename Mesh::cell_type& cl,
                                   const CellDegreeInfo<Mesh>&     cell_infos,
                                   bool hF = true)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto celdeg = cell_infos.cell_degree();

    const auto faces_infos = cell_infos.facesDegreeInfo();

    const auto num_diff_dofs  = scalar_diff_dofs(msh, cell_infos);
    const auto num_faces_dofs = scalar_faces_dofs(msh, faces_infos);

    matrix_type data = matrix_type::Zero(num_diff_dofs, num_diff_dofs);

    const auto fcs         = faces(msh, cl);
    size_t     offset      = 0;
    size_t     offset_diff = 0;

    T h = 0;
    if (!hF)
    {
        h = diameter(msh, cl);
    }
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fdi = faces_infos[i];

        if (fdi.hasUnknowns())
        {
            const auto fc     = fcs[i];
            const auto facdeg = fdi.degree();
            if (hF)
            {
                h = diameter(msh, fc);
            }

            const auto diff_deg = std::max(celdeg, facdeg);
            const auto fb       = make_scalar_monomial_basis(msh, fc, facdeg);
            const auto db       = make_scalar_monomial_basis(msh, fc, diff_deg);

            const auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);
            const auto dbs = scalar_basis_size(diff_deg, Mesh::dimension - 1);

            matrix_type mass  = matrix_type::Zero(fbs, fbs);
            matrix_type trace = matrix_type::Zero(fbs, dbs);

            const auto qps = integrate(msh, fc, facdeg + diff_deg);
            for (auto& qp : qps)
            {
                const auto d_phi = db.eval_functions(qp.point());
                const auto f_phi = fb.eval_functions(qp.point());

                mass += priv::outer_product(priv::inner_product(qp.weight(), f_phi), f_phi);
                trace += priv::outer_product(priv::inner_product(qp.weight(), f_phi), d_phi);
            }

            matrix_type oper = matrix_type::Zero(fbs, num_diff_dofs);
            matrix_type tr   = matrix_type::Zero(fbs, num_diff_dofs);

            tr.block(0, offset_diff, fbs, dbs) = trace;

            oper.block(0, offset_diff, fbs, dbs) = mass.ldlt().solve(trace);
            data += oper.transpose() * tr * (1. / h);

            offset += fbs;
            offset_diff += dbs;
        }
    }

    return data;
}

/**
 * @brief compute the stabilization term \f$\sum_{F \in F_T} 1/h_F(u_F - u_T + \Pi^k_T R^{k+1}_T(\hat{u}_T) -
 * R^{k+1}_T(\hat{u}_T), v_F - v_T + \Pi^k_T R^{k+1}_T(\hat{v}_T) - R^{k+1}_T(\hat{v}_T))_F \f$ for scalar HHO
 * unknowns
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
make_scalar_hdg_stabilization_diff(const Mesh&                     msh,
                                   const typename Mesh::cell_type& cl,
                                   const hho_degree_info&          di,
                                   bool                            hF = true)
{
    const CellDegreeInfo<Mesh> cell_infos(msh, cl, di.cell_degree(), di.face_degree(), di.grad_degree());

    return make_scalar_hdg_stabilization_diff(msh, cl, cell_infos, hF);
}

/**
 * @brief compute the stabilization term \f$ \sum_{F \in F_T} 1/h_F(u_F-u_T, v_F-v_T)_F \f$
 * for scalar HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param cell_infos cell degree informations
 * @param hF use diameter of face for scaling if true (or cell diameter if false)
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_dg_stabilization(const Mesh&                     msh,
                             const typename Mesh::cell_type& cl,
                             const CellDegreeInfo<Mesh>&     cell_infos,
                             bool                            hF = true)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto celdeg = cell_infos.cell_degree();

    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);

    const auto faces_infos    = cell_infos.facesDegreeInfo();
    const auto num_faces_dofs = scalar_faces_dofs(msh, faces_infos);
    const auto total_dofs     = cbs + num_faces_dofs;

    matrix_type data = matrix_type::Zero(total_dofs, total_dofs);

    const auto cb = make_scalar_monomial_basis(msh, cl, celdeg);

    const auto fcs    = faces(msh, cl);
    size_t     offset = cbs;
    T          h      = 0;
    if (!hF)
    {
        h = diameter(msh, cl);
    }
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fdi = faces_infos[i];

        if (fdi.hasUnknowns())
        {
            const auto fc     = fcs[i];
            const auto facdeg = fdi.degree();
            if (hF)
            {
                h = diameter(msh, fc);
            }
            const auto fb     = make_scalar_monomial_basis(msh, fc, facdeg);
            const auto fbs    = scalar_basis_size(facdeg, Mesh::dimension - 1);

            matrix_type mass_F = make_mass_matrix(msh, fc, fb);
            matrix_type mass_T = make_mass_matrix(msh, fc, cb);
            matrix_type trace  = matrix_type::Zero(fbs, cbs);

            const auto qps = integrate(msh, fc, facdeg + celdeg);
            for (auto& qp : qps)
            {
                const auto c_phi = cb.eval_functions(qp.point());
                const auto f_phi = fb.eval_functions(qp.point());

                trace += priv::outer_product(priv::inner_product(qp.weight(), f_phi), c_phi);
            }

            data.block(0, 0, cbs, cbs) += mass_T / h;
            data.block(offset, offset, fbs, fbs) += mass_F / h;
            data.block(0, offset, cbs, fbs) -= trace.transpose() / h;
            data.block(offset, 0, fbs, cbs) -= trace / h;

            offset += fbs;
        }
    }

    return data;
}

/**
 * @brief compute the stabilization term \f$ \sum_{F \in F_T} 1/h_F(u_F-\Pi^k_F(u_T), v_F-\Pi^k_F(v_T))_F \f$
 * for scalar HHO unknowns
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
make_scalar_hdg_stabilization(const Mesh&                     msh,
                              const typename Mesh::cell_type& cl,
                              const hho_degree_info&          di,
                              bool                            hF = true,
                              bool                            scaling = true)
{
    const CellDegreeInfo<Mesh> cell_infos(msh, cl, di.cell_degree(), di.face_degree(), di.grad_degree());

    return make_scalar_hdg_stabilization(msh, cl, cell_infos, hF, scaling);
}

/**
 * @brief compute the stabilization term \f$ \sum_{F \in F_T} 1/h_F(u_F-\Pi^k_F(u_T), v_F-\Pi^k_F(v_T))_F \f$
 * for scalar HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param msh_infos mesh degree informations
 * @param hF use diameter of face for scaling if true (or cell diameter if false)
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_hdg_stabilization(const Mesh&                     msh,
                              const typename Mesh::cell_type& cl,
                              const MeshDegreeInfo<Mesh>&     msh_infos,
                              bool                            hF = true,
                              bool                            scaling = true)
{
    return make_scalar_hdg_stabilization(msh, cl, msh_infos.cellDegreeInfo(msh, cl), hF, scaling);
}

/**
 * @brief compute the stabilization term \f$ \sum_{F \in F_T} 1/h_F(u_F-u_T, v_F-v_T)_F \f$
 * for scalar HHO unknowns
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
make_scalar_dg_stabilization(const Mesh&                     msh,
                             const typename Mesh::cell_type& cl,
                             const hho_degree_info&          di,
                             bool                            hF = true)
{
    const CellDegreeInfo<Mesh> cell_infos(msh, cl, di.cell_degree(), di.face_degree(), di.grad_degree());

    return make_scalar_dg_stabilization(msh, cl, cell_infos, hF);
}

/**
 * @brief compute the stabilization term \f$ \sum_{F \in F_T} 1/h_F(u_F-u_T, v_F-v_T)_F \f$
 * for scalar HHO unknowns
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
make_scalar_dg_stabilization(const Mesh&                     msh,
                             const typename Mesh::cell_type& cl,
                             const MeshDegreeInfo<Mesh>&     msh_infos,
                             bool                            hF = true)
{
    return make_scalar_dg_stabilization(msh, cl, msh_infos.cellDegreeInfo(msh, cl), hF);
}

/**
 * @brief compute the stabilization term \f$\sum_{F \in F_T} 1/h_F(u_F - u_T + \Pi^k_T R^{k+1}_T(\hat{u}_T) -
 * R^{k+1}_T(\hat{u}_T), v_F - v_T + \Pi^k_T R^{k+1}_T(\hat{v}_T) - R^{k+1}_T(\hat{v}_T))_F \f$ for scalar HHO
 * unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param reconstruction reconstruction operator \f$ R^{k+1}_T \f$
 * @param cell_infos cell degree information
 * @param hF use diameter of face for scaling if true (or cell diameter if false)
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_hho_stabilization(const Mesh&                     msh,
        const typename Mesh::cell_type&                       cl,
        const dynamic_matrix<typename Mesh::coordinate_type>& reconstruction,
        const CellDegreeInfo<Mesh>&                           cell_infos,
        bool                                                  hF = true)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto recdeg = cell_infos.reconstruction_degree();
    const auto celdeg = cell_infos.cell_degree();

    const auto rbs = scalar_basis_size(recdeg, Mesh::dimension);
    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);

    const auto cb = make_scalar_monomial_basis(msh, cl, recdeg);

    const matrix_type mass_mat = make_mass_matrix(msh, cl, cb);

    // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

    // Step 1: compute \pi_T^k p_T^k v (third term).
    const matrix_type M1    = mass_mat.block(0, 0, cbs, cbs);
    const matrix_type M2    = mass_mat.block(0, 1, cbs, rbs - 1);
    matrix_type       proj1 = -M1.ldlt().solve(M2 * reconstruction);

    assert(M2.cols() == reconstruction.rows());

    // Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
    proj1.block(0, 0, cbs, cbs) += matrix_type::Identity(cbs, cbs);

    const auto fcs            = faces(msh, cl);
    const auto faces_infos    = cell_infos.facesDegreeInfo();
    const auto num_faces_dofs = scalar_faces_dofs(msh, faces_infos);

    matrix_type data = matrix_type::Zero(cbs + num_faces_dofs, cbs + num_faces_dofs);

    T h = 0;
    if (!hF)
    {
        h = diameter(msh, cl);
    }

    // Step 3: project on faces (eqn. 21)
    size_t offset = cbs;
    for (size_t face_i = 0; face_i < fcs.size(); face_i++)
    {
        const auto fdi = faces_infos[face_i];

        if (fdi.hasUnknowns())
        {
            const auto fc = fcs[face_i];
            if (hF)
            {
                h = diameter(msh, fc);
            }

            const auto facdeg = fdi.degree();
            const auto fb     = make_scalar_monomial_basis(msh, fc, facdeg);
            const auto fbs    = scalar_basis_size(facdeg, Mesh::dimension - 1);

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
            const matrix_type MR1 = face_trace_matrix.block(0, 1, fbs, rbs - 1);

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


template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_hho_stabilization2(const Mesh&                      msh,
        const typename Mesh::cell_type&                         cl,
        const dynamic_matrix<typename Mesh::coordinate_type>&   reconstruction,
        const hho_degree_info&                                  hdi,
        const static_matrix<typename Mesh::coordinate_type, Mesh::dimension, Mesh::dimension>& kappaT)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto celdeg = hdi.cell_degree();
    const auto recdeg = hdi.reconstruction_degree();

    const auto cb = make_scalar_monomial_basis(msh, cl, celdeg);
    const auto rb = make_scalar_monomial_basis(msh, cl, recdeg);
    
    const auto cbs = cb.size();
    const auto rbs = rb.size();

    const matrix_type mass_mat = make_mass_matrix(msh, cl, rb);

    // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

    // Step 1: compute \pi_T^k p_T^k v (third term).
    const matrix_type M1    = mass_mat.block(0, 0, cbs, cbs);
    const matrix_type M2    = mass_mat.block(0, 1, cbs, rbs - 1);
    matrix_type       proj1 = -M1.ldlt().solve(M2 * reconstruction);

    assert(M2.cols() == reconstruction.rows());

    // Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
    proj1.block(0, 0, cbs, cbs) += matrix_type::Identity(cbs, cbs);

    const auto fcs = faces(msh, cl);
    const auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    const auto num_faces_dofs = fbs*fcs.size();

    matrix_type data = matrix_type::Zero(cbs + num_faces_dofs, cbs + num_faces_dofs);

    T h = diameter(msh, cl);

    // Step 3: project on faces (eqn. 21)
    size_t offset = cbs;
    for (size_t face_i = 0; face_i < fcs.size(); face_i++)
    {
        const auto fc = fcs[face_i];

        const auto facdeg = hdi.face_degree();
        const auto fb     = make_scalar_monomial_basis(msh, fc, facdeg);
        const auto fbs    = fb.size();

        matrix_type face_mass_matrix  = matrix_type::Zero(fbs, fbs);
        matrix_type face_trace_matrix = matrix_type::Zero(fbs, rbs);

        const auto face_quadpoints = integrate(msh, fc, recdeg + facdeg);
        for (auto& qp : face_quadpoints)
        {
            const auto f_phi = fb.eval_functions(qp.point());
            const auto r_phi = rb.eval_functions(qp.point());
            face_mass_matrix += qp.weight() * f_phi * f_phi.transpose(); 
            face_trace_matrix += qp.weight() * f_phi * r_phi.transpose();
        }

        LLT<matrix_type> piKF;
        piKF.compute(face_mass_matrix);

        // Step 3a: \pi_F^k( v_F - p_T^k v )
        const matrix_type MR1 = face_trace_matrix.block(0, 1, fbs, rbs - 1);

        matrix_type proj2 = piKF.solve(MR1 * reconstruction);
        proj2.block(0, offset, fbs, fbs) -= matrix_type::Identity(fbs, fbs);

        // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
        const matrix_type MR2   = face_trace_matrix.block(0, 0, fbs, cbs);
        const matrix_type proj3 = piKF.solve(MR2 * proj1);
        const matrix_type BRF   = proj2 + proj3;

        auto n = normal(msh, cl, fc);
        T stabparam = n.dot(kappaT*n)/h;

        data += stabparam * BRF.transpose() * face_mass_matrix * BRF;

        offset += fbs;
    }

    return data;
}

/**
 * @brief compute the stabilization term \f$\sum_{F \in F_T} 1/h_F(u_F - u_T + \Pi^k_T R^{k+1}_T(\hat{u}_T) -
 * R^{k+1}_T(\hat{u}_T), v_F - v_T + \Pi^k_T R^{k+1}_T(\hat{v}_T) - R^{k+1}_T(\hat{v}_T))_F \f$ for scalar HHO
 * unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param reconstruction reconstruction operator \f$ R^{k+1}_T \f$
 * @param msh_infos mesh degree information
 * @param hF use diameter of face for scaling if true (or cell diameter if false)
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_hho_stabilization(const Mesh&                     msh,
        const typename Mesh::cell_type&                       cl,
        const dynamic_matrix<typename Mesh::coordinate_type>& reconstruction,
        const MeshDegreeInfo<Mesh>&                           msh_infos,
        bool                                                  hF = true)
{
    return make_scalar_hho_stabilization(msh, cl, reconstruction, msh_infos.cellDegreeInfo(msh, cl), hF);
}

/**
 * @brief compute the stabilization term \f$\sum_{F \in F_T} 1/h_F(u_F - u_T + \Pi^k_T R^{k+1}_T(\hat{u}_T) -
 * R^{k+1}_T(\hat{u}_T), v_F - v_T + \Pi^k_T R^{k+1}_T(\hat{v}_T) - R^{k+1}_T(\hat{v}_T))_F \f$ for scalar HHO
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
template<typename Mesh, typename ScalarReconstructionBasis>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_hho_stabilization(const Mesh&                 msh,
    const typename Mesh::cell_type&                       cl,
    const ScalarReconstructionBasis&                      rb,
    const dynamic_matrix<typename Mesh::coordinate_type>& reconstruction,
    const hho_degree_info&                                hdi)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto recdeg = hdi.reconstruction_degree();
    const auto celdeg = hdi.cell_degree();
    const auto facdeg = hdi.face_degree();

    const auto cb = make_scalar_monomial_basis(msh, cl, celdeg);

    const auto rbs = rb.size();
    const auto cbs = cb.size();
    const auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);

    // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

    // Step 1: compute \pi_T^k p_T^k v (third term).
    const matrix_type M1    = make_mass_matrix(msh, cl, cb);

    matrix_type M2    = matrix_type::Zero(cbs, rbs - 1);
    auto qps = integrate(msh, cl, recdeg+celdeg);
    for(auto& qp : qps)
    {
        auto trial_tmp = rb.eval_functions( qp.point() );
        auto trial = trial_tmp.tail(rbs-1);
        auto test = cb.eval_functions( qp.point() );
        M2 += qp.weight() * test * trial.transpose();
    }

    matrix_type       proj1 = -M1.ldlt().solve(M2 * reconstruction);

    assert(M2.cols() == reconstruction.rows());

    // Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
    proj1.block(0, 0, cbs, cbs) += matrix_type::Identity(cbs, cbs);

    const auto fcs            = faces(msh, cl);
    const auto num_faces_dofs = fcs.size() * fbs;

    matrix_type data = matrix_type::Zero(cbs + num_faces_dofs, cbs + num_faces_dofs);

    T h = diameter(msh, cl);

    // Step 3: project on faces (eqn. 21)
    size_t offset = cbs;
    for (size_t face_i = 0; face_i < fcs.size(); face_i++)
    {
        const auto& fc = fcs[face_i];
        const auto fb     = make_scalar_monomial_basis(msh, fc, facdeg);
        const auto fbs    = scalar_basis_size(facdeg, Mesh::dimension - 1);

        matrix_type face_mass_matrix  = make_mass_matrix(msh, fc, fb);
        matrix_type face_trace_matrix = matrix_type::Zero(fbs, rbs-1);
        matrix_type face_trace_matrix2 = matrix_type::Zero(fbs, cbs);

        const auto face_quadpoints = integrate(msh, fc, recdeg + facdeg);
        for (auto& qp : face_quadpoints)
        {
            const auto f_phi = fb.eval_functions(qp.point());
            const auto c_phi = cb.eval_functions(qp.point());
            const auto r_phi_tmp = rb.eval_functions(qp.point());
            const auto r_phi = r_phi_tmp.tail(rbs-1);
            face_trace_matrix += qp.weight() * f_phi * r_phi.transpose();
            face_trace_matrix2 += qp.weight() * f_phi * c_phi.transpose();
        }

        LLT<matrix_type> piKF;
        piKF.compute(face_mass_matrix);

        // Step 3a: \pi_F^k( v_F - p_T^k v )
        const matrix_type MR1 = face_trace_matrix;

        matrix_type proj2 = piKF.solve(MR1 * reconstruction);
        proj2.block(0, offset, fbs, fbs) -= matrix_type::Identity(fbs, fbs);

        // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
        const matrix_type MR2   = face_trace_matrix2;
        const matrix_type proj3 = piKF.solve(MR2 * proj1);
        const matrix_type BRF   = proj2 + proj3;

        data += BRF.transpose() * face_mass_matrix * BRF / h;

        offset += fbs;
    }

    return data;
}

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_hho_stabilization(const Mesh&                     msh,
        const typename Mesh::cell_type&                       cl,
        const dynamic_matrix<typename Mesh::coordinate_type>& reconstruction,
        const hho_degree_info&                                hdi,
        bool x = false)
{
    auto rb = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
    return make_scalar_hho_stabilization(msh, cl, rb, reconstruction, hdi);
}

/**
 * @brief compute the stabilization term \f$\sum_{F \in F_T} 1/h_F(u_F - u_T + \Pi^k_T R^{k+1}_T(\hat{u}_T) -
 * R^{k+1}_T(\hat{u}_T), v_F - v_T + \Pi^k_T R^{k+1}_T(\hat{v}_T) - R^{k+1}_T(\hat{v}_T))_F \f$ for scalar HHO
 * unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param reconstruction reconstruction operator \f$ R^{k+1}_T \f$
 * @param cell_infos cell degree information
 * @param hF use diameter of face for scaling if true (or cell diameter if false)
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_hho_stabilization_fluxes(const Mesh&              msh,
        const typename Mesh::cell_type&                       cl,
        const dynamic_matrix<typename Mesh::coordinate_type>& reconstruction,
        const CellDegreeInfo<Mesh>&                           cell_infos,
        bool                                                  hF = true)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const matrix_type stab = make_scalar_hho_stabilization(msh, cl, reconstruction, cell_infos);

    const auto celdeg = cell_infos.cell_degree();

    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);

    const auto cb = make_scalar_monomial_basis(msh, cl, celdeg);

    const auto fcs            = faces(msh, cl);
    const auto faces_infos    = cell_infos.facesDegreeInfo();
    const auto num_faces_dofs = scalar_faces_dofs(msh, faces_infos);

    matrix_type matA = matrix_type::Zero(num_faces_dofs, cbs + num_faces_dofs);

    T h = 0;
    if (!hF)
    {
        h = diameter(msh, cl);
    }

    // Step 3: project on faces (eqn. 21)
    size_t offset = cbs;
    size_t row    = 0;
    for (size_t face_i = 0; face_i < fcs.size(); face_i++)
    {
        const auto fdi = faces_infos[face_i];

        if (fdi.hasUnknowns())
        {
            const auto fc = fcs[face_i];
            if (hF)
            {
                h = diameter(msh, fc);
            }

            const auto facdeg = fdi.degree();
            const auto fb     = make_scalar_monomial_basis(msh, fc, facdeg);
            const auto fbs    = scalar_basis_size(facdeg, Mesh::dimension - 1);

            matrix_type face_mass_matrix  = make_mass_matrix(msh, fc, fb);
            matrix_type face_trace_matrix = matrix_type::Zero(fbs, cbs);

            const auto face_quadpoints = integrate(msh, fc, celdeg + facdeg);
            for (auto& qp : face_quadpoints)
            {
                const auto f_phi = fb.eval_functions(qp.point());
                const auto c_phi = cb.eval_functions(qp.point());
                face_trace_matrix += priv::outer_product(priv::inner_product(qp.weight(), f_phi), c_phi);
            }

            matA.block(row, offset, fbs, fbs) += face_mass_matrix / h;
            matA.block(row, 0, fbs, cbs) -= face_trace_matrix / h;

            offset += fbs;
            row += fbs;
        }
    }

    return (matA.transpose()).fullPivHouseholderQr().solve(stab);
}

/**
 * @brief compute the stabilization term \f$\sum_{F \in F_T} 1/h_F(u_F - u_T + \Pi^k_T R^{k+1}_T(\hat{u}_T) -
 * R^{k+1}_T(\hat{u}_T), v_F - v_T + \Pi^k_T R^{k+1}_T(\hat{v}_T) - R^{k+1}_T(\hat{v}_T))_F \f$ for scalar HHO
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
make_scalar_hho_stabilization_fluxes(const Mesh&              msh,
        const typename Mesh::cell_type&                       cl,
        const dynamic_matrix<typename Mesh::coordinate_type>& reconstruction,
        const hho_degree_info&                                di,
        bool                                                  hF = true)
{
    const CellDegreeInfo<Mesh> cell_infos(msh, cl, di.cell_degree(), di.face_degree(), di.grad_degree());

    return make_scalar_hho_stabilization_fluxes(msh, cl, reconstruction, cell_infos, hF);
}

/**
 * @brief compute the stabilization term \f$\sum_{F \in F_T} 1/h_F(u_F - u_T + \Pi^k_T R^{k+1}_T(\hat{u}_T) -
 * R^{k+1}_T(\hat{u}_T), v_F - v_T + \Pi^k_T R^{k+1}_T(\hat{v}_T) - R^{k+1}_T(\hat{v}_T))_F \f$ for scalar HHO
 * unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param cell_infos cell degree information
 * @param hF use diameter of face for scaling if true (or cell diameter if false)
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_dg_stabilization_fluxes(const Mesh&                     msh,
                                    const typename Mesh::cell_type& cl,
                                    const CellDegreeInfo<Mesh>&     cell_infos,
                                    bool                            hF = true)
{
    const auto stab    = make_scalar_dg_stabilization_diff(msh, cl, cell_infos, hF);
    const auto adjoint = make_scalar_stabilization_fluxes(msh, cl, cell_infos, stab, hF);
    return adjoint * make_scalar_hho_difference(msh, cl, cell_infos, hF);
}

/**
 * @brief compute the stabilization term \f$\sum_{F \in F_T} 1/h_F(u_F - u_T + \Pi^k_T R^{k+1}_T(\hat{u}_T) -
 * R^{k+1}_T(\hat{u}_T), v_F - v_T + \Pi^k_T R^{k+1}_T(\hat{v}_T) - R^{k+1}_T(\hat{v}_T))_F \f$ for scalar HHO
 * unknowns
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
make_scalar_dg_stabilization_fluxes(const Mesh&                     msh,
                                    const typename Mesh::cell_type& cl,
                                    const hho_degree_info&          di,
                                    bool                            hF = true)
{
    const CellDegreeInfo<Mesh> cell_infos(msh, cl, di.cell_degree(), di.face_degree(), di.grad_degree());

    return make_scalar_dg_stabilization_fluxes(msh, cl, cell_infos, hF);
}

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_stabilization_fluxes(const Mesh&                                          msh,
                                 const typename Mesh::cell_type&                      cl,
                                 const CellDegreeInfo<Mesh>&                          cell_infos,
                                 const dynamic_matrix<typename Mesh::coordinate_type> stab,
                                 bool                                                 hF = true)
{
    const auto mass = make_scalar_dg_stabilization_diff(msh, cl, cell_infos, hF);

    return mass.ldlt().solve(stab);
}

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_dg_stabilization_diff(const Mesh&                     msh,
                                  const typename Mesh::cell_type& cl,
                                  const CellDegreeInfo<Mesh>&     cell_infos,
                                  bool                            hF = true)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto celdeg = cell_infos.cell_degree();

    const auto faces_infos = cell_infos.facesDegreeInfo();

    const auto num_diff_dofs  = scalar_diff_dofs(msh, cell_infos);
    const auto num_faces_dofs = scalar_faces_dofs(msh, faces_infos);

    matrix_type data = matrix_type::Zero(num_diff_dofs, num_diff_dofs);

    T h = 0;
    if (!hF)
    {
        h = diameter(msh, cl);
    }

    const auto fcs         = faces(msh, cl);
    size_t     offset_diff = 0;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fdi = faces_infos[i];

        if (fdi.hasUnknowns())
        {
            const auto fc     = fcs[i];
            const auto facdeg = fdi.degree();
            if (hF)
            {
                h = diameter(msh, fc);
            }

            const auto diff_deg = std::max(celdeg, facdeg);
            const auto db       = make_scalar_monomial_basis(msh, fc, diff_deg);
            const auto dbs      = scalar_basis_size(diff_deg, Mesh::dimension - 1);

            matrix_type mass = make_mass_matrix(msh, fc, db);

            data.block(offset_diff, offset_diff, dbs, dbs) = mass / h;

            offset_diff += dbs;
        }
    }

    return data;
}

} // end diskpp