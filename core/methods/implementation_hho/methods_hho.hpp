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
 * Karol Cascavita (C) 2018                     karol.cascavita@enpc.fr
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

#include <iterator>

#include "adaptivity/adaptivity.hpp"
#include "bases/bases.hpp"
#include "boundary_conditions/boundary_conditions.hpp"
#include "common/eigen.hpp"
#include "quadratures/quadratures.hpp"
#include "utils_hho.hpp"

namespace disk
{

/**
 * @brief Compute the interpolation operator \f$ \hat{I}(u) := (\Pi^{cd}_T(u), \Pi^{fd}_{\partial T}(u)) \f$
 * of a function \f$ u \in H^1(T;\mathbb{R}) \f$ where \f$ \Pi^{cd}_T \f$, respc. \f$ \Pi^{fd}_{\partial T} \f$ are the
 * L2-orthogonal projector on the cell and faces
 *
 * @tparam Mesh type of the mesh
 * @param msh  mesh
 * @param cl  cell
 * @param hdi  hho degree informations
 * @param u  function  \f$ u \in H^1(T;\mathbb{R}) \f$ to project
 * @param di additonal quadratue degree
 * @return dynamic_vector<typename Mesh::coordinate_type> coefficient of the the interpolation operator \f$ \hat{I}(u)
 * \f$
 */
template<typename Mesh>
dynamic_vector<typename Mesh::coordinate_type>
project_function(const Mesh&                      msh,
                 const typename Mesh::cell_type&  cl,
                 const hho_degree_info&           hdi,
                 const scalar_rhs_function<Mesh>& u,
                 size_t                           di = 0)
{
    typedef dynamic_vector<typename Mesh::coordinate_type> vector_type;

    const auto cbs       = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    const auto fbs       = scalar_basis_size(hdi.face_degree(), Mesh::dimension - 1);
    const auto num_faces = howmany_faces(msh, cl);

    vector_type ret = vector_type::Zero(cbs + num_faces * fbs);

    ret.block(0, 0, cbs, 1) = project_function(msh, cl, hdi.cell_degree(), u, di);

    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < num_faces; i++)
    {
        ret.segment(cbs + i * fbs, fbs) = project_function(msh, fcs[i], hdi.face_degree(), u, di);
    }

    return ret;
}

/**
 * @brief Compute the interpolation operator \f$ \hat{I}(u) := (\Pi^{cd}_T(u), \Pi^{fd}_{\partial T}(u)) \f$
 * of a function \f$ u \in H^1(T;\mathbb{R}^d) \f$ where \f$ \Pi^{cd}_T \f$, respc. \f$ \Pi^{fd}_{\partial T} \f$ are
 * the L2-orthogonal projector on the cell and faces
 *
 * @tparam Mesh type of the mesh
 * @param msh  mesh
 * @param cl  cell
 * @param hdi  hho degree informations
 * @param u  function  \f$ u \in H^1(T;\mathbb{R}^d) \f$ to project
 * @param di additonal quadratue degree
 * @return dynamic_vector<typename Mesh::coordinate_type> coefficient of the the interpolation operator \f$ \hat{I}(u)
 * \f$
 */
template<typename Mesh>
dynamic_vector<typename Mesh::coordinate_type>
project_function(const Mesh&                      msh,
                 const typename Mesh::cell_type&  cl,
                 const hho_degree_info&           hdi,
                 const vector_rhs_function<Mesh>& u,
                 size_t                           di = 0)
{
    typedef dynamic_vector<typename Mesh::coordinate_type> vector_type;

    const auto cbs       = vector_basis_size(hdi.cell_degree(), Mesh::dimension, Mesh::dimension);
    const auto fbs       = vector_basis_size(hdi.face_degree(), Mesh::dimension - 1, Mesh::dimension);
    const auto num_faces = howmany_faces(msh, cl);

    vector_type ret = vector_type::Zero(cbs + num_faces * fbs);

    ret.block(0, 0, cbs, 1) = project_function(msh, cl, hdi.cell_degree(), u, di);

    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < num_faces; i++)
    {
        ret.segment(cbs + i * fbs, fbs) = project_function(msh, fcs[i], hdi.face_degree(), u, di);
    }

    return ret;
}

/**
 * @brief Compute the interpolation operator \f$ \hat{I}(u) := (\Pi^{cd}_T(u), \Pi^{fd}_{\partial T}(u)) \f$
 * of a function \f$ u \in H^1(T;\mathbb{R}^d) \f$ where \f$ \Pi^{cd}_T \f$, respc. \f$ \Pi^{fd}_{\partial T} \f$ are
 * the L2-orthogonal projector on the cell and faces
 *
 * @tparam Mesh type of the mesh
 * @param msh  mesh
 * @param cl  cell
 * @param hdi  hho degree informations
 * @param u  function  \f$ u \in H^1(T;\mathbb{R}^d) \f$ to project
 * @param di additonal quadratue degree
 * @return dynamic_vector<typename Mesh::coordinate_type> coefficient of the the interpolation operator \f$ \hat{I}(u)
 * \f$
 */
template<typename Mesh>
dynamic_vector<typename Mesh::coordinate_type>
project_tangent(const Mesh&                      msh,
                const typename Mesh::cell_type&  cl,
                const hho_degree_info&           hdi,
                const vector_rhs_function<Mesh>& u,
                size_t                           di = 0)
{
    typedef dynamic_vector<typename Mesh::coordinate_type> vector_type;

    const auto cbs       = vector_basis_size(hdi.cell_degree(), Mesh::dimension, Mesh::dimension);
    const auto fbs       = vector_basis_size(hdi.face_degree(), Mesh::dimension - 1, Mesh::dimension-1);
    const auto num_faces = howmany_faces(msh, cl);

    vector_type ret = vector_type::Zero(cbs + num_faces * fbs);

    ret.block(0, 0, cbs, 1) = project_function(msh, cl, hdi.cell_degree(), u, di);

    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < num_faces; i++)
    {
        auto n = normal(msh, cl, fcs[i]);
        auto trace_u = [&](const typename Mesh::point_type& pt) -> auto {
            auto uval = u(pt);
            return n.cross(uval.cross(n));
        };

        auto fc = fcs[i];
        const auto fb = make_vector_monomial_tangential_basis(msh, fc, hdi.face_degree());
        ret.segment(cbs + i * fbs, fbs) = project_function(msh, fcs[i], fb, trace_u, di);
    }

    return ret;
}

template<typename Mesh>
dynamic_vector<typename Mesh::coordinate_type>
project_tangent_nedelec(const Mesh&                      msh,
                const typename Mesh::cell_type&  cl,
                const hho_degree_info&           hdi,
                const vector_rhs_function<Mesh>& u,
                size_t                           di = 0)
{
    typedef dynamic_vector<typename Mesh::coordinate_type> vector_type;

    const auto cbs       = vector_basis_size(hdi.cell_degree(), Mesh::dimension, Mesh::dimension);
    const auto fbs       = nedelec_tangential_basis_size(hdi.face_degree());
    const auto num_faces = howmany_faces(msh, cl);

    vector_type ret = vector_type::Zero(cbs + num_faces * fbs);

    ret.block(0, 0, cbs, 1) = project_function(msh, cl, hdi.cell_degree(), u, di);

    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < num_faces; i++)
    {
        auto n = normal(msh, cl, fcs[i]);
        auto trace_u = [&](const typename Mesh::point_type& pt) -> auto {
            auto uval = u(pt);
            return n.cross(uval.cross(n));
        };

        auto fc = fcs[i];
        const auto fb = make_vector_monomial_nedelec_tangential_basis(msh, fc, hdi.face_degree());
        ret.segment(cbs + i * fbs, fbs) = project_function(msh, fcs[i], fb, trace_u, di);
    }

    return ret;
}

/**
 * @brief Compute the interpolation operator \f$ \hat{I}(u) := (\Pi^{cd}_T(u), \Pi^{fd}_{\partial T}(u)) \f$
 * of a function \f$ u \in H^1(T;\mathbb{R}) \f$ where \f$ \Pi^{cd}_T \f$, respc. \f$ \Pi^{fd}_{\partial T} \f$ are the
 * L2-orthogonal projector on the cell and faces
 *
 * @tparam Mesh type of the mesh
 * @param msh  mesh
 * @param cl  cell
 * @param hdi  hho degree informations
 * @param u  function  \f$ u \in H^1(T;\mathbb{R}) \f$ to project
 * @param di additonal quadratue degree
 * @return dynamic_vector<typename Mesh::coordinate_type> coefficient of the the interpolation operator \f$ \hat{I}(u)
 * \f$
 */
template<typename Mesh>
dynamic_vector<typename Mesh::coordinate_type>
project_function(const Mesh&                      msh,
                 const typename Mesh::cell_type&  cl,
                 const MeshDegreeInfo<Mesh>&      degree_infos,
                 const scalar_rhs_function<Mesh>& u,
                 size_t                           di = 0)
{
    typedef dynamic_vector<typename Mesh::coordinate_type> vector_type;

    const auto cdi    = degree_infos.cellDegreeInfo(msh, cl);
    const auto fcs_di = cdi.facesDegreeInfo();

    const auto num_faces_dofs = scalar_faces_dofs(msh, cl, fcs_di);

    const auto cbs = scalar_basis_size(cdi.cell_degree(), Mesh::dimension);

    vector_type ret = vector_type::Zero(cbs + num_faces_dofs);

    ret.head(cbs) = project_function(msh, cl, cdi.cell_degree(), u, di);

    const auto fcs    = faces(msh, cl);
    size_t     offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const DegreeInfo fdi = fcs_di[i];

        if (fdi.hasUnknowns())
        {
            const size_t num_face_dofs         = scalar_face_dofs(msh, fdi);
            ret.segment(offset, num_face_dofs) = project_function(msh, fcs[i], fdi.degree(), u, di);
            offset += num_face_dofs;
        }
    }

    return ret;
}

/**
 * @brief Compute the interpolation operator \f$ \hat{I}(u) := (\Pi^{cd}_T(u), \Pi^{fd}_{\partial T}(u)) \f$
 * of a function \f$ u \in H^1(T;\mathbb{R}^d) \f$ where \f$ \Pi^{cd}_T \f$, respc. \f$ \Pi^{fd}_{\partial T} \f$ are
 * the L2-orthogonal projector on the cell and faces
 *
 * @tparam Mesh type of the mesh
 * @param msh  mesh
 * @param cl  cell
 * @param hdi  hho degree informations
 * @param u  function  \f$ u \in H^1(T;\mathbb{R}^d) \f$ to project
 * @param di additonal quadratue degree
 * @return dynamic_vector<typename Mesh::coordinate_type> coefficient of the the interpolation operator \f$ \hat{I}(u)
 * \f$
 */
template<typename Mesh>
dynamic_vector<typename Mesh::coordinate_type>
project_function(const Mesh&                      msh,
                 const typename Mesh::cell_type&  cl,
                 const MeshDegreeInfo<Mesh>&      degree_infos,
                 const vector_rhs_function<Mesh>& u,
                 size_t                           di = 0)
{
    typedef dynamic_vector<typename Mesh::coordinate_type> vector_type;

    const auto cdi    = degree_infos.cellDegreeInfo(msh, cl);
    const auto fcs_di = cdi.facesDegreeInfo();

    const auto num_faces_dofs = vector_faces_dofs(msh, fcs_di);

    const auto cbs = vector_basis_size(cdi.cell_degree(), Mesh::dimension, Mesh::dimension);

    vector_type ret = vector_type::Zero(cbs + num_faces_dofs);

    ret.head(cbs) = project_function(msh, cl, cdi.cell_degree(), u, di);

    const auto fcs    = faces(msh, cl);
    size_t     offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const DegreeInfo fdi = fcs_di[i];

        if (fdi.hasUnknowns())
        {
            const size_t num_face_dofs         = vector_face_dofs(msh, fdi);
            ret.segment(offset, num_face_dofs) = project_function(msh, fcs[i], fdi.degree(), u, di);
            offset += num_face_dofs;
        }
    }

    return ret;
}

template<typename Mesh>
auto
make_hho_stokes(const Mesh&                     msh,
                const typename Mesh::cell_type& cl,
                const hho_degree_info&          hdi,
                const bool&                     use_sym_grad)
{
    if (use_sym_grad)
        return make_vector_hho_symmetric_laplacian(msh, cl, hdi);
    else
        return make_vector_hho_laplacian(msh, cl, hdi);
}

template<typename Mesh>
auto
make_hlow_stokes(const Mesh&                     msh,
                 const typename Mesh::cell_type& cl,
                 const hho_degree_info&          hdi,
                 const bool&                     use_sym_grad)
{
    if (use_sym_grad)
        return make_matrix_symmetric_gradrec(msh, cl, hdi);
    else
        return make_marix_hho_gradrec(msh, cl, hdi);
}


} // end diskpp