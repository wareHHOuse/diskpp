/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
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

#include <cassert>

#include "bases/bases.hpp"
#include "common/eigen.hpp"
#include "mechanics/behaviors/maths_tensor.hpp"
#include "methods/hho"
#include "quadratures/quadratures.hpp"

#include "timecounter.h"

namespace MK
{

template<typename Mesh>
std::pair<Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
          Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>>
make_matrix_symmetric_gradrec(const Mesh&                                   msh,
                              const typename Mesh::cell_type&               cl,
                              const disk::hho_degree_info&                  di,
                              const disk::vector_boundary_conditions<Mesh>& bnd)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const size_t N = Mesh::dimension;

    const auto ci = disk::contact_info<Mesh>(msh, cl, di, bnd);

    const auto graddeg = di.grad_degree();
    const auto celdeg  = ci.cell_degree();

    const auto gb = make_sym_matrix_monomial_basis(msh, cl, graddeg);
    const auto cb = make_vector_monomial_basis(msh, cl, celdeg);

    const auto gbs = disk::sym_matrix_basis_size(graddeg, Mesh::dimension, Mesh::dimension);
    const auto cbs = ci.num_cell_dofs();

    const auto fcs       = ci.faces();
    const auto num_faces = ci.num_faces();

    matrix_type gr_lhs = matrix_type::Zero(gbs, gbs);
    matrix_type gr_rhs = matrix_type::Zero(gbs, ci.num_total_dofs());

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
            const auto qp_gphi_j = disk::priv::inner_product(qp.weight(), gphi[j]);
            for (size_t i = j; i < gbs; i += dec)
                gr_lhs(i, j) += disk::priv::inner_product(gphi[i], qp_gphi_j);
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
            const auto qp_dphi = disk::priv::inner_product(qp.weight(), dphi);

            gr_rhs.block(0, 0, gbs, cbs) += disk::priv::outer_product(gphi, qp_dphi);

        } // end qp
    }

    size_t offset = cbs;

    for (size_t i = 0; i < num_faces; i++)
    {
        const auto fc     = fcs[i];
        const auto facdeg = ci.face_degree(bnd, fc);
        const auto n      = normal(msh, cl, fc);
        const auto fb     = make_vector_monomial_basis(msh, fc, facdeg);
        const auto fbs    = fb.size();

        const auto qps_f = integrate(msh, fc, graddeg + std::max(celdeg, facdeg));
        for (auto& qp : qps_f)
        {
            const auto cphi = cb.eval_functions(qp.point());
            const auto gphi = gb.eval_functions(qp.point());
            const auto fphi = fb.eval_functions(qp.point());

            const auto qp_gphi_n = disk::priv::inner_product(gphi, disk::priv::inner_product(qp.weight(), n));
            gr_rhs.block(0, offset, gbs, fbs) += disk::priv::outer_product(qp_gphi_n, fphi);
            gr_rhs.block(0, 0, gbs, cbs) -= disk::priv::outer_product(qp_gphi_n, cphi);
        }

        offset += fbs;
    }

    matrix_type oper = gr_lhs.ldlt().solve(gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
std::pair<Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
          Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>>
make_vector_hho_symmetric_laplacian(const Mesh&                                   msh,
                                    const typename Mesh::cell_type&               cl,
                                    const disk::hho_degree_info&                  di,
                                    const disk::vector_boundary_conditions<Mesh>& bnd)
{
    using T        = typename Mesh::coordinate_type;
    const size_t N = Mesh::dimension;

    const auto ci = disk::contact_info<Mesh>(msh, cl, di, bnd);

    const auto recdeg = di.reconstruction_degree();
    const auto celdeg = ci.cell_degree();

    const auto rb = make_vector_monomial_basis(msh, cl, recdeg);
    const auto cb = make_vector_monomial_basis(msh, cl, celdeg);

    const auto rbs = disk::vector_basis_size(recdeg, N, N);
    const auto cbs = disk::vector_basis_size(celdeg, N, N);

    const auto fcs       = ci.faces();
    const auto num_faces = fcs.size();

    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, N, N>             gradient_type;
    typedef Matrix<T, Dynamic, N>       function_type;

    const size_t rbs_ho         = rbs - N;
    const size_t num_total_dofs = ci.num_total_dofs();
    const size_t nb_lag         = disk::priv::nb_lag(N);

    matrix_type stiff  = matrix_type::Zero(rbs, rbs);
    matrix_type gr_lhs = matrix_type::Zero(rbs_ho + nb_lag, rbs_ho + nb_lag);
    matrix_type gr_rhs = matrix_type::Zero(rbs_ho + nb_lag, num_total_dofs);

    const auto qps = integrate(msh, cl, 2 * (recdeg - 1));
    for (auto& qp : qps)
    {
        const auto dphi    = rb.eval_sgradients(qp.point());
        const auto qp_dphi = disk::priv::inner_product(qp.weight(), dphi);
        stiff += disk::priv::outer_product(qp_dphi, dphi);
    }

    gr_lhs.block(0, 0, rbs_ho, rbs_ho) = stiff.block(N, N, rbs_ho, rbs_ho);
    gr_rhs.block(0, 0, rbs_ho, cbs)    = stiff.block(N, 0, rbs_ho, cbs);

    size_t offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc     = fcs[i];
        const auto facdeg = ci.face_degree(bnd, fc);
        const auto fbs    = disk::vector_basis_size(facdeg, N - 1, N);
        const auto n      = normal(msh, cl, fc);
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
            const function_type qp_r_dphi_n = qp.weight() * disk::priv::inner_product(r_dphi, n);
            gr_rhs.block(0, offset, rbs_ho, fbs) += disk::priv::outer_product(qp_r_dphi_n, f_phi);
            gr_rhs.block(0, 0, rbs_ho, cbs) -= disk::priv::outer_product(qp_r_dphi_n, c_phi);
        }

        offset += fbs;
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
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_vector_hho_stabilization(const Mesh&                                                     msh,
                              const typename Mesh::cell_type&                                 cl,
                              const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>& reconstruction,
                              const disk::hho_degree_info&                                    di,
                              const disk::vector_boundary_conditions<Mesh>&                   bnd)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    const size_t                        N = Mesh::dimension;

    const auto ci = disk::contact_info<Mesh>(msh, cl, di, bnd);

    const auto recdeg = di.reconstruction_degree();
    const auto celdeg = ci.cell_degree();

    const auto rbs = disk::vector_basis_size(recdeg, Mesh::dimension, Mesh::dimension);
    const auto cbs = disk::vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);

    const auto num_total_dofs = ci.num_total_dofs();

    const auto cb = make_vector_monomial_basis(msh, cl, recdeg);

    const matrix_type mass_mat = make_mass_matrix(msh, cl, cb);

    // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

    // Step 1: compute \pi_T^k p_T^k v (third term).
    const matrix_type M1    = mass_mat.block(0, 0, cbs, cbs);
    const matrix_type M2    = mass_mat.block(0, N, cbs, rbs - N);
    matrix_type       proj1 = -M1.ldlt().solve(M2 * reconstruction);

    assert(M2.cols() == reconstruction.rows());

    // Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
    proj1.block(0, 0, cbs, cbs) += matrix_type::Identity(cbs, cbs);

    const auto fcs       = ci.faces();
    const auto num_faces = fcs.size();

    matrix_type data = matrix_type::Zero(num_total_dofs, num_total_dofs);

    // Step 3: project on faces (eqn. 21)
    size_t offset = cbs;
    for (size_t face_i = 0; face_i < num_faces; face_i++)
    {
        const auto fc     = fcs[face_i];
        const auto facdeg = ci.face_degree(bnd, fc);
        const auto fbs    = disk::vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension);
        const auto hf     = diameter(msh, fc);
        const auto fb     = make_vector_monomial_basis(msh, fc, facdeg);

        matrix_type face_mass_matrix  = make_mass_matrix(msh, fc, fb);
        matrix_type face_trace_matrix = matrix_type::Zero(fbs, rbs);

        const auto face_quadpoints = integrate(msh, fc, recdeg + facdeg);
        for (auto& qp : face_quadpoints)
        {
            const auto f_phi = fb.eval_functions(qp.point());
            const auto c_phi = cb.eval_functions(qp.point());
            face_trace_matrix += disk::priv::outer_product(disk::priv::inner_product(qp.weight(), f_phi), c_phi);
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

        data += BRF.transpose() * face_mass_matrix * BRF / hf;

        offset += fbs;
    }

    return data;
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_scalar_hdg_stabilization(const Mesh&                                   msh,
                              const typename Mesh::cell_type&               cl,
                              const disk::hho_degree_info&                  di,
                              const disk::vector_boundary_conditions<Mesh>& bnd)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto ci = disk::contact_info<Mesh>(msh, cl, di, bnd);

    const auto celdeg = ci.cell_degree();

    const auto cbs = disk::scalar_basis_size(celdeg, Mesh::dimension);

    const auto fcs        = ci.faces();
    const auto num_faces  = fcs.size();
    const auto total_dofs = ci.num_total_dofs();

    matrix_type data = matrix_type::Zero(total_dofs, total_dofs);

    auto cb = make_scalar_monomial_basis(msh, cl, celdeg);

    size_t offset = cbs;

    for (size_t i = 0; i < num_faces; i++)
    {
        const auto fc     = fcs[i];
        const auto h      = diameter(msh, fc);
        const auto facdeg = ci.face_degree(bnd, fc);
        const auto fbs    = disk::scalar_basis_size(facdeg, Mesh::dimension - 1);
        const auto fb     = make_scalar_monomial_basis(msh, fc, facdeg);

        const matrix_type If = matrix_type::Identity(fbs, fbs);

        matrix_type oper  = matrix_type::Zero(fbs, total_dofs);
        matrix_type tr    = matrix_type::Zero(fbs, total_dofs);
        matrix_type mass  = make_mass_matrix(msh, fc, fb);
        matrix_type trace = matrix_type::Zero(fbs, cbs);

        oper.block(0, offset, fbs, fbs) = -If;

        const auto qps = integrate(msh, fc, facdeg + celdeg);
        for (auto& qp : qps)
        {
            const auto c_phi = cb.eval_functions(qp.point());
            const auto f_phi = fb.eval_functions(qp.point());

            assert(c_phi.rows() == cbs);
            assert(f_phi.rows() == fbs);
            assert(c_phi.cols() == f_phi.cols());

            trace += disk::priv::outer_product(disk::priv::inner_product(qp.weight(), f_phi), c_phi);
        }

        tr.block(0, offset, fbs, fbs) = -mass;
        tr.block(0, 0, fbs, cbs)      = trace;

        offset += fbs;

        oper.block(0, 0, fbs, cbs) = mass.ldlt().solve(trace);
        data += oper.transpose() * tr * (1. / h);
    }

    return data;
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_vector_hdg_stabilization(const Mesh&                                   msh,
                              const typename Mesh::cell_type&               cl,
                              const disk::hho_degree_info&                  di,
                              const disk::vector_boundary_conditions<Mesh>& bnd)
{
    const auto hdg_scalar_stab = make_scalar_hdg_stabilization(msh, cl, di, bnd);

    return disk::priv::compute_lhs_vector(msh, cl, di, hdg_scalar_stab, bnd);
}

template<typename Mesh, typename T>
auto
make_vector_static_condensation_withMatrix(const Mesh&                                                      msh,
                                           const typename Mesh::cell_type&                                  cl,
                                           const disk::hho_degree_info&                                     di,
                                           const disk::vector_boundary_conditions<Mesh>&                    bnd,
                                           const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                                           const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs)
{
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    const auto ci = disk::contact_info<Mesh>(msh, cl, di, bnd);

    const auto fcs            = ci.faces();
    const auto num_faces      = fcs.size();
    const auto num_cell_dofs  = ci.num_cell_dofs();
    const auto num_faces_dofs = ci.num_faces_dofs();
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

    const auto K_TT_lu = K_TT.lu();
    // if (K_TT_lu.info() != Eigen::Success)
    // {
    //     throw std::invalid_argument("static condensation: K_TT is not positive definite");
    // }

    const matrix_type AL = K_TT_lu.solve(K_TF);
    const vector_type bL = K_TT_lu.solve(cell_rhs);

    const matrix_type AC = K_FF - K_FT * AL;
    const vector_type bC = faces_rhs - K_FT * bL;

    return std::make_tuple(std::make_pair(AC, bC), AL, bL);
}

namespace priv
{
template<typename T>
T
compute_g0(const disk::point<T, 2>& pt, const static_vector<T, 2>& n)
{
    // compute the distance to the plane y = 0

    if (std::abs(n(1)) < T(1E-12))
        return 10E12;

    const T dist = std::abs(pt.y() / n(1));

    if (pt.y() < T(0))
        return -dist;

    return dist;

    // compute the distance to the plane y = -0.05(x-0.5)^2
    // return pt.y() + 0.05 * (pt.x() - 0.5) * (pt.x() - 0.5);
}

template<typename T>
T
compute_g0(const disk::point<T, 3>& pt, const static_vector<T, 3>& n)
{
    // distance to the plane z = 0

    // the normal is orthogonal to the plane z = 0
    if (std::abs(n(2)) < T(1E-12))
        return 10E12;

    const T dist = std::abs(pt.z() / n(2));

    if (pt.z() < T(0))
        return -dist;

    return dist;
}

template<typename Basis, typename T, size_t DIM>
point<T, DIM>
new_pt(const Basis& base, const disk::dynamic_vector<T>& tab_coeff, const point<T, DIM>& pt)
{
    const auto beval = base.eval_functions(pt);

    const auto depl = disk::eval(tab_coeff, beval);

    return pt + depl;
}

template<typename T>
static_vector<T, 2>
compute_normal(const disk::point<T, 2>& a, const disk::point<T, 2>& b)
{
    const static_vector<T, 2> t1 = (b - a).to_vector();
    static_vector<T, 2>       nor;
    nor(0) = -t1(1);
    nor(1) = t1(0);

    return nor / nor.norm();
}

template<typename T>
static_vector<T, 3>
compute_normal(const disk::point<T, 3>& p1, const disk::point<T, 3>& p2, const disk::point<T, 3>& p3)
{
    static_vector<T, 3> t1 = (p2 - p1).to_vector();
    static_vector<T, 3> t2 = (p3 - p1).to_vector();

    t1 /= t1.norm();
    t2 /= t2.norm();

    static_vector<T, 3> nor = cross(t1, t2);

    return nor / nor.norm();
}

template<typename Mesh, typename Elem, typename ElemBasis, typename T>
T
compute_gap_fb(const Mesh&                msh,
               const Elem&                elem,
               const ElemBasis&           eb,
               const disk::dynamic_vector<T>&   tab_coeff,
               const disk::point<T, 2>&         pt,
               const static_vector<T, 2>& n)
{
    const disk::point<T, 2>         pt_def = new_pt(eb, tab_coeff, pt);
    const auto                pts    = points(msh, elem);
    const static_vector<T, 2> n_ref  = compute_normal(pts[0], pts[1]);

    const disk::point<T, 2> pta_def = new_pt(eb, tab_coeff, pts[0]);
    const disk::point<T, 2> ptb_def = new_pt(eb, tab_coeff, pts[1]);

    const T sign = std::copysign(T(1), n.dot(n_ref));

    const static_vector<T, 2> n_def = sign * compute_normal(pta_def, ptb_def);

    // std::cout << "pt: " << pt << std::endl;
    // std::cout << "pta: " << pts[0] << std::endl;
    // std::cout << "ptb: " << pts[1] << std::endl;
    // std::cout << "pt def: " << pt_def << std::endl;
    // std::cout << "pta def: " << pta_def << std::endl;
    // std::cout << "ptb def: " << ptb_def << std::endl;

    // std::cout << "normal: " << n.transpose() << std::endl;
    // std::cout << "normal ref: " << n_ref.transpose() << std::endl;
    // std::cout << "normal def: " << n_def.transpose() << std::endl;

    return compute_g0(pt_def, n_def);
}

template<typename Mesh, typename Elem, typename ElemBasis, typename T>
T
compute_gap_fb(const Mesh&                msh,
               const Elem&                elem,
               const ElemBasis&           eb,
               const disk::dynamic_vector<T>&   tab_coeff,
               const disk::point<T, 3>&         pt,
               const static_vector<T, 3>& n)
{
    const disk::point<T, 3>         pt_def = new_pt(eb, tab_coeff, pt);
    const auto                pts    = points(msh, elem);
    const static_vector<T, 3> n_ref  = compute_normal(pts[0], pts[1], pts[2]);

    const disk::point<T, 3> pt0_def = new_pt(eb, tab_coeff, pts[0]);
    const disk::point<T, 3> pt1_def = new_pt(eb, tab_coeff, pts[1]);
    const disk::point<T, 3> pt2_def = new_pt(eb, tab_coeff, pts[2]);

    const T sign = std::copysign(T(1), n.dot(n_ref));

    const static_vector<T, 3> n_def = sign * compute_normal(pt0_def, pt1_def, pt2_def);

    return compute_g0(pt_def, n_def);
}
}

template<typename MeshType>
class tresca
{
  private:
    typedef MeshType                                    mesh_type;
    typedef typename mesh_type::coordinate_type         scalar_type;
    typedef typename mesh_type::cell                    cell_type;
    typedef typename mesh_type::face                    face_type;
    typedef point<scalar_type, mesh_type::dimension>    point_type;
    typedef disk::MaterialData<scalar_type>             material_type;
    typedef typename disk::hho_degree_info              hdi_type;
    typedef ParamRun<scalar_type>                       param_type;
    typedef disk::vector_boundary_conditions<mesh_type> bnd_type;
    typedef typename disk::contact_info<mesh_type>      ci_type;

    const static int dimension = mesh_type::dimension;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    typedef static_matrix<scalar_type, dimension, dimension> matrix_static;
    typedef static_vector<scalar_type, dimension>            vector_static;

    const mesh_type&     m_msh;
    const hdi_type&      m_hdi;
    const material_type& m_material_data;
    const param_type&    m_rp;
    const bnd_type&      m_bnd;

    size_t cell_basis_size, grad_basis_size, num_total_dofs;

    bool two_dim;

    // contact contrib;

    int
    num_dofs_dim() const
    {
        if (two_dim)
            return 3;
        else
            return 6;
    }

    matrix_type
    make_hho_aT(const cell_type& cl, const matrix_type& ET, const matrix_type& ST) const
    {
        // grad_sym part
        matrix_type mass_grad = matrix_type::Zero(grad_basis_size, grad_basis_size);
        const auto  gb        = make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());

        // this is very costly to build it
        const auto qps = integrate(m_msh, cl, 2 * m_hdi.grad_degree());

        size_t dec = 0;
        if (dimension == 3)
            dec = 6;
        else if (dimension == 2)
            dec = 3;
        else
            std::logic_error("Expected 3 >= dim > 1");

        for (auto& qp : qps)
        {
            const auto gphi = gb.eval_functions(qp.point());

            for (size_t j = 0; j < grad_basis_size; j++)
            {
                const auto qp_gphi_j = disk::priv::inner_product(qp.weight(), gphi[j]);
                for (size_t i = j; i < grad_basis_size; i += dec)
                    mass_grad(i, j) += disk::priv::inner_product(gphi[i], qp_gphi_j);
            }
        }

        // upper part
        for (size_t j = 0; j < grad_basis_size; j++)
            for (size_t i = 0; i < j; i++)
                mass_grad(i, j) = mass_grad(j, i);

        matrix_type aT = 2.0 * m_material_data.getMu() * ET.transpose() * mass_grad * ET;

        // divergence part
        matrix_type mass_div = matrix_type::Zero(grad_basis_size, grad_basis_size);

        for (auto& qp : qps)
        {
            const auto gphi       = gb.eval_functions(qp.point());
            const auto gphi_trace = disk::trace(gphi);

            const auto qp_gphi_trace = disk::priv::inner_product(qp.weight(), gphi_trace);

            mass_div += disk::priv::outer_product(qp_gphi_trace, gphi_trace);
        }

        aT += m_material_data.getLambda() * ET.transpose() * mass_div * ET;

        // stabilization part
        if (m_rp.m_stab)
        {
            aT += m_rp.m_beta * ST;
        }

        return aT;
    }

    scalar_type
    make_hho_distance(const point_type& pt, const vector_static& n) const
    {
        return priv::compute_g0(pt, n);
    }

    // normal part of u : u_n = u.n
    template<typename TraceBasis>
    vector_type
    make_hho_u_n(const vector_static& n, const TraceBasis& tb, const point_type& pt) const
    {
        const auto t_phi = tb.eval_functions(pt);

        //(phi_T . n)
        return disk::priv::inner_product(t_phi, n);
    }

    // tangential part of u : u_t = u - u_n*n
    template<typename TraceBasis>
    Matrix<scalar_type, Dynamic, dimension>
    make_hho_u_t(const vector_static& n, const TraceBasis& tb, const point_type& pt) const
    {
        const auto t_phi = tb.eval_functions(pt);
        const auto u_n   = make_hho_u_n(n, tb, pt);

        // phi_T - (phi_T . n)n
        return t_phi - disk::priv::inner_product(u_n, n);
    }

    // cauchy traction : sigma_n = sigma * n
    template<typename GradBasis>
    Matrix<scalar_type, Dynamic, dimension>
    make_hho_sigma_n(const matrix_type& ET, const vector_static& n, const GradBasis& gb, const point_type& pt) const
    {
        Matrix<scalar_type, Dynamic, dimension> sigma_n =
          Matrix<scalar_type, Dynamic, dimension>::Zero(num_total_dofs, dimension);

        const auto gphi   = gb.eval_functions(pt);
        const auto gphi_n = disk::priv::inner_product(gphi, n);

        sigma_n = 2.0 * m_material_data.getMu() * ET.transpose() * gphi_n;

        const auto gphi_trace_n = disk::priv::inner_product(disk::trace(gphi), n);

        sigma_n += m_material_data.getLambda() * ET.transpose() * gphi_trace_n;

        // sigma. n
        return sigma_n;
    }

    template<typename GradBasis>
    vector_type
    make_hho_sigma_nn(const matrix_type& ET, const vector_static& n, const GradBasis& gb, const point_type& pt) const
    {
        const auto sigma_n = make_hho_sigma_n(ET, n, gb, pt);

        // sigma_n . n
        return disk::priv::inner_product(sigma_n, n);
    }

    vector_type
    make_hho_sigma_nn(const Matrix<scalar_type, Dynamic, dimension>& sigma_n, const vector_static& n) const
    {
        // sigma_n . n
        return disk::priv::inner_product(sigma_n, n);
    }

    template<typename GradBasis>
    Matrix<scalar_type, Dynamic, dimension>
    make_hho_sigma_nt(const matrix_type& ET, const vector_static& n, const GradBasis& gb, const point_type& pt) const
    {
        const auto sigma_n  = make_hho_sigma_n(ET, n, gb, pt);
        const auto sigma_nn = make_hho_sigma_nn(sigma_n, n);

        // sigma_n - sigma_nn * n
        return sigma_n - disk::priv::inner_product(sigma_nn, n);
    }

    vector_type
    make_hho_phi_n_uT(const vector_type& sigma_nn,
                      const vector_type& uT_n,
                      const scalar_type& theta,
                      const scalar_type& gamma_F) const
    {
        vector_type phi_n = theta * sigma_nn;
        phi_n.head(cell_basis_size) -= gamma_F * uT_n;

        // theta * sigma_nn - gamma uT_n
        return phi_n;
    }

    Matrix<scalar_type, Dynamic, dimension>
    make_hho_phi_t_uT(const Matrix<scalar_type, Dynamic, dimension>& sigma_nt,
                      const Matrix<scalar_type, Dynamic, dimension>& uT_t,
                      const scalar_type&                             theta,
                      const scalar_type&                             gamma_F) const
    {
        Matrix<scalar_type, Dynamic, dimension> phi_t = theta * sigma_nt;
        phi_t.block(0, 0, cell_basis_size, dimension) -= gamma_F * uT_t;

        // theta * sigma_nt - gamma uT_t
        return phi_t;
    }

    vector_type
    make_hho_phi_n_uF(const vector_type& sigma_nn,
                      const vector_type& uF_n,
                      const scalar_type& theta,
                      const scalar_type& gamma_F,
                      const size_t&      offset) const
    {
        vector_type phi_n = theta * sigma_nn;

        assert(offset + uF_n.size() <= phi_n.size());

        phi_n.segment(offset, uF_n.size()) -= gamma_F * uF_n;

        // theta * sigma_nn - gamma uF_n
        return phi_n;
    }

    Matrix<scalar_type, Dynamic, dimension>
    make_hho_phi_t_uF(const Matrix<scalar_type, Dynamic, dimension>& sigma_nt,
                      const Matrix<scalar_type, Dynamic, dimension>& uF_t,
                      const scalar_type&                             theta,
                      const scalar_type&                             gamma_F,
                      const size_t&                                  offset) const
    {
        Matrix<scalar_type, Dynamic, dimension> phi_t = theta * sigma_nt;

        phi_t.block(offset, 0, uF_t.rows(), dimension) -= gamma_F * uF_t;

        // theta * sigma_nt - gamma uF_t
        return phi_t;
    }

    // projection on the ball of radius alpha centered on 0
    vector_static
    make_proj_alpha(const vector_static& x, const scalar_type& alpha) const
    {
        const scalar_type x_norm = x.norm();

        if (x_norm <= alpha)
        {
            return x;
        }

        return alpha * x / x_norm;
    }

    // derivative of the projection on the ball of radius alpha centered on 0
    matrix_static
    make_d_proj_alpha(const vector_static& x, const scalar_type& alpha) const
    {
        const scalar_type x_norm = x.norm();

        if (x_norm <= alpha)
        {
            return matrix_static::Identity();
        }

        return alpha / x_norm * (matrix_static::Identity() - disk::Kronecker(x, x) / (x_norm * x_norm));
    }

    // compute theta/gamma *(sigma_n, sigma_n)_Fc
    matrix_type
    make_hho_nitsche(const cell_type& cl, const matrix_type& ET, const ci_type& ci) const
    {
        const auto gb = make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());

        matrix_type nitsche = matrix_type::Zero(num_total_dofs, num_total_dofs);

        const auto fcs = m_bnd.faces_with_contact(cl);
        for (auto& fc : fcs)
        {
            const auto n       = normal(m_msh, cl, fc);
            const auto qps     = integrate(m_msh, fc, 2 * m_hdi.grad_degree() + 2);
            const auto hF      = diameter(m_msh, fc);
            const auto gamma_F = m_rp.m_gamma_0 / hF;

            for (auto& qp : qps)
            {
                const auto sigma_n    = make_hho_sigma_n(ET, n, gb, qp.point());
                const auto qp_sigma_n = disk::priv::inner_product(qp.weight() / gamma_F, sigma_n);

                nitsche += disk::priv::outer_product(qp_sigma_n, sigma_n);
            }
        }

        return m_rp.m_theta * nitsche;
    }

    // compute (phi_n_theta, H(-phi_n_1(u))*phi_n_1)_FC / gamma
    matrix_type
    make_hho_heaviside_contact(const cell_type&   cl,
                               const matrix_type& ET,
                               const vector_type& stress_coeff,
                               const vector_type& uTF,
                               const ci_type&     ci) const
    {
        const auto cb = make_vector_monomial_basis(m_msh, cl, m_hdi.cell_degree());
        const auto gb = make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());

        matrix_type lhs = matrix_type::Zero(num_total_dofs, num_total_dofs);

        const auto fcs    = faces(m_msh, cl);
        size_t     offset = cell_basis_size;
        for (auto& fc : fcs)
        {
            const auto facedeg = ci.face_degree(m_bnd, fc);
            const auto fbs     = ci.num_face_dofs(m_bnd, fc);

            if (m_bnd.is_contact_face(fc))
            {
                const auto contact_type = m_bnd.contact_boundary_type(fc);
                const auto n            = normal(m_msh, cl, fc);
                const auto qp_deg       = std::max(m_hdi.cell_degree(), m_hdi.grad_degree());
                const auto qps          = integrate(m_msh, fc, 2 * qp_deg + 2);
                const auto hF           = diameter(m_msh, fc);
                const auto gamma_F      = m_rp.m_gamma_0 / hF;

                const auto fb = make_vector_monomial_basis(m_msh, fc, facedeg);

                for (auto& qp : qps)
                {
                    const vector_type sigma_nn = make_hho_sigma_nn(ET, n, gb, qp.point());

                    if (contact_type == disk::SIGNORINI_CELL)
                    {
                        const vector_type uT_n = make_hho_u_n(n, cb, qp.point());

                        const scalar_type phi_n_1_u = eval_phi_n_uT(stress_coeff, gb, cb, uTF, n, gamma_F, qp.point());

                        // Heaviside(-phi_n_1(u))
                        if (phi_n_1_u <= scalar_type(0))
                        {
                            const vector_type phi_n_theta = make_hho_phi_n_uT(sigma_nn, uT_n, m_rp.m_theta, gamma_F);
                            const vector_type phi_n_1     = make_hho_phi_n_uT(sigma_nn, uT_n, scalar_type(1), gamma_F);
                            const auto qp_phi_n_theta = disk::priv::inner_product(qp.weight() / gamma_F, phi_n_theta);

                            lhs += disk::priv::outer_product(qp_phi_n_theta, phi_n_1);
                        }
                    }
                    else
                    {
                        const vector_type uF_n = make_hho_u_n(n, fb, qp.point());

                        const scalar_type phi_n_1_u =
                          eval_phi_n_uF(fc, stress_coeff, gb, fb, uTF, offset, n, gamma_F, qp.point());

                        // Heaviside(-phi_n_1(u))
                        if (phi_n_1_u <= scalar_type(0))
                        {
                            const vector_type phi_n_theta =
                              make_hho_phi_n_uF(sigma_nn, uF_n, m_rp.m_theta, gamma_F, offset);
                            const vector_type phi_n_1 =
                              make_hho_phi_n_uF(sigma_nn, uF_n, scalar_type(1), gamma_F, offset);
                            const auto qp_phi_n_theta = disk::priv::inner_product(qp.weight() / gamma_F, phi_n_theta);

                            lhs += disk::priv::outer_product(qp_phi_n_theta, phi_n_1);
                        }
                    }
                }
            }
            offset += fbs;
        }

        return lhs;
    }

    // compute (phi_n_theta, [phi_n_1(u)]R-)_FC / gamma
    vector_type
    make_hho_negative_contact(const cell_type&   cl,
                              const matrix_type& ET,
                              const vector_type& stress_coeff,
                              const vector_type& uTF,
                              const ci_type&     ci) const
    {
        const auto cb = make_vector_monomial_basis(m_msh, cl, m_hdi.cell_degree());
        const auto gb = make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());

        vector_type rhs = vector_type::Zero(num_total_dofs);

        const auto fcs    = faces(m_msh, cl);
        size_t     offset = cell_basis_size;
        for (auto& fc : fcs)
        {
            const auto facedeg = ci.face_degree(m_bnd, fc);
            const auto fbs     = ci.num_face_dofs(m_bnd, fc);

            if (m_bnd.is_contact_face(fc))
            {
                const auto contact_type = m_bnd.contact_boundary_type(fc);
                const auto n            = normal(m_msh, cl, fc);
                const auto qp_deg       = std::max(m_hdi.cell_degree(), m_hdi.grad_degree());
                const auto qps          = integrate(m_msh, fc, 2 * qp_deg + 2);
                const auto hF           = diameter(m_msh, fc);
                const auto gamma_F      = m_rp.m_gamma_0 / hF;

                const auto fb = make_vector_monomial_basis(m_msh, fc, facedeg);

                for (auto& qp : qps)
                {
                    const vector_type sigma_nn = make_hho_sigma_nn(ET, n, gb, qp.point());

                    if (contact_type == disk::SIGNORINI_CELL)
                    {
                        const vector_type uT_n = make_hho_u_n(n, cb, qp.point());

                        const scalar_type phi_n_1_u = eval_phi_n_uT(stress_coeff, gb, cb, uTF, n, gamma_F, qp.point());

                        // [phi_n_1_u]_R-
                        if (phi_n_1_u <= scalar_type(0))
                        {
                            const vector_type phi_n_theta = make_hho_phi_n_uT(sigma_nn, uT_n, m_rp.m_theta, gamma_F);

                            rhs += (qp.weight() / gamma_F * phi_n_1_u) * phi_n_theta;
                        }
                    }
                    else
                    {
                        const vector_type uF_n = make_hho_u_n(n, fb, qp.point());
                        const scalar_type phi_n_1_u =
                          eval_phi_n_uF(fc, stress_coeff, gb, fb, uTF, offset, n, gamma_F, qp.point());

                        // [phi_n_1_u]_R-
                        if (phi_n_1_u <= scalar_type(0))
                        {
                            const vector_type phi_n_theta =
                              make_hho_phi_n_uF(sigma_nn, uF_n, m_rp.m_theta, gamma_F, offset);

                            rhs += (qp.weight() / gamma_F * phi_n_1_u) * phi_n_theta;
                        }
                    }
                }
            }
            offset += fbs;
        }
        return rhs;
    }

    // compute (phi_t_theta, [phi_t_1(u)]_(s))_FC / gamma
    vector_type
    make_hho_threshold_tresca(const cell_type&   cl,
                              const matrix_type& ET,
                              const vector_type& stress_coeff,
                              const vector_type& uTF,
                              const ci_type&     ci) const
    {
        const auto cb = make_vector_monomial_basis(m_msh, cl, m_hdi.cell_degree());
        const auto gb = make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());

        vector_type rhs = vector_type::Zero(num_total_dofs);

        const auto fcs    = faces(m_msh, cl);
        size_t     offset = cell_basis_size;
        for (auto& fc : fcs)
        {
            const auto facedeg = ci.face_degree(m_bnd, fc);
            const auto fbs     = ci.num_face_dofs(m_bnd, fc);

            if (m_bnd.is_contact_face(fc))
            {
                const auto contact_type = m_bnd.contact_boundary_type(fc);
                const auto n            = normal(m_msh, cl, fc);
                const auto qp_deg       = std::max(m_hdi.cell_degree(), m_hdi.grad_degree());
                const auto qps          = integrate(m_msh, fc, 2 * qp_deg + 2);
                const auto hF           = diameter(m_msh, fc);
                const auto gamma_F      = m_rp.m_gamma_0 / hF;

                const auto fb     = make_vector_monomial_basis(m_msh, fc, facedeg);
                const auto s_func = m_bnd.contact_boundary_func(fc);

                for (auto& qp : qps)
                {
                    const auto sigma_nt = make_hho_sigma_nt(ET, n, gb, qp.point());

                    if (contact_type == disk::SIGNORINI_CELL)
                    {
                        const auto uT_t = make_hho_u_t(n, cb, qp.point());

                        const auto phi_t_theta = make_hho_phi_t_uT(sigma_nt, uT_t, m_rp.m_theta, gamma_F);

                        const vector_static phi_t_1_u_proj =
                          eval_proj_phi_t_uT(stress_coeff, gb, cb, uTF, n, gamma_F, s_func(qp.point()), qp.point());

                        const vector_static qp_phi_t_1_u_pro = qp.weight() * phi_t_1_u_proj / gamma_F;

                        rhs += disk::priv::inner_product(phi_t_theta, qp_phi_t_1_u_pro);
                    }
                    else
                    {
                        const auto uF_t = make_hho_u_t(n, fb, qp.point());

                        const auto phi_t_theta = make_hho_phi_t_uF(sigma_nt, uF_t, m_rp.m_theta, gamma_F, offset);

                        const vector_static phi_t_1_u_proj = eval_proj_phi_t_uF(
                          stress_coeff, gb, fb, uTF, offset, n, gamma_F, s_func(qp.point()), qp.point());

                        const vector_static qp_phi_t_1_u_pro = qp.weight() * phi_t_1_u_proj / gamma_F;

                        rhs += disk::priv::inner_product(phi_t_theta, qp_phi_t_1_u_pro);
                    }
                }
            }
            offset += fbs;
        }
        return rhs;
    }

    // compute (phi_t_theta, (d_proj_alpha(u)) phi_t_1)_FC / gamma
    matrix_type
    make_hho_matrix_tresca(const cell_type&   cl,
                           const matrix_type& ET,
                           const vector_type& stress_coeff,
                           const vector_type& uTF,
                           const ci_type&     ci) const
    {
        const auto cb = make_vector_monomial_basis(m_msh, cl, m_hdi.cell_degree());
        const auto gb = make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());

        matrix_type lhs = matrix_type::Zero(num_total_dofs, num_total_dofs);

        const auto fcs    = faces(m_msh, cl);
        size_t     offset = cell_basis_size;
        for (auto& fc : fcs)
        {
            const auto facedeg = ci.face_degree(m_bnd, fc);
            const auto fbs     = ci.num_face_dofs(m_bnd, fc);

            if (m_bnd.is_contact_face(fc))
            {
                const auto contact_type = m_bnd.contact_boundary_type(fc);
                const auto n            = normal(m_msh, cl, fc);
                const auto qp_deg       = std::max(m_hdi.cell_degree(), m_hdi.grad_degree());
                const auto qps          = integrate(m_msh, fc, 2 * qp_deg + 2);
                const auto hF           = diameter(m_msh, fc);
                const auto gamma_F      = m_rp.m_gamma_0 / hF;

                const auto fb     = make_vector_monomial_basis(m_msh, fc, facedeg);
                const auto s_func = m_bnd.contact_boundary_func(fc);

                for (auto& qp : qps)
                {
                    const auto sigma_nt = make_hho_sigma_nt(ET, n, gb, qp.point());

                    if (contact_type == disk::SIGNORINI_CELL)
                    {
                        const auto uT_t = make_hho_u_t(n, cb, qp.point());

                        const auto phi_t_1     = make_hho_phi_t_uT(sigma_nt, uT_t, scalar_type(1), gamma_F);
                        const auto phi_t_theta = make_hho_phi_t_uT(sigma_nt, uT_t, m_rp.m_theta, gamma_F);

                        const auto phi_t_1_u      = eval_phi_t_uT(stress_coeff, gb, cb, uTF, n, gamma_F, qp.point());
                        const auto d_proj_phi_t_u = make_d_proj_alpha(phi_t_1_u, s_func(qp.point()));

                        const auto d_proj_u_phi_t_1 = disk::priv::inner_product(d_proj_phi_t_u, phi_t_1);

                        const auto qp_phi_t_theta = disk::priv::inner_product(qp.weight() / gamma_F, phi_t_theta);

                        lhs += disk::priv::outer_product(qp_phi_t_theta, d_proj_u_phi_t_1);
                    }
                    else
                    {
                        const auto uF_t = make_hho_u_t(n, fb, qp.point());

                        const auto phi_t_1     = make_hho_phi_t_uF(sigma_nt, uF_t, scalar_type(1), gamma_F, offset);
                        const auto phi_t_theta = make_hho_phi_t_uF(sigma_nt, uF_t, m_rp.m_theta, gamma_F, offset);

                        const auto phi_t_1_u = eval_phi_t_uF(stress_coeff, gb, fb, uTF, offset, n, gamma_F, qp.point());
                        const auto d_proj_phi_t_u = make_d_proj_alpha(phi_t_1_u, s_func(qp.point()));

                        const auto d_proj_u_phi_t_1 = disk::priv::inner_product(d_proj_phi_t_u, phi_t_1);

                        const auto qp_phi_t_theta = disk::priv::inner_product(qp.weight() / gamma_F, phi_t_theta);

                        lhs += disk::priv::outer_product(qp_phi_t_theta, d_proj_u_phi_t_1);
                    }
                }
            }
            offset += fbs;
        }

        return lhs;
    }

  public:
    matrix_type K_int;
    vector_type RTF;
    vector_type F_int;
    double      time_contact;

    tresca(const mesh_type&     msh,
           const hdi_type&      hdi,
           const material_type& material_data,
           const param_type&    rp,
           const bnd_type&      bnd) :
      m_msh(msh),
      m_hdi(hdi), m_material_data(material_data), m_rp(rp), m_bnd(bnd)
    {
        cell_basis_size = disk::vector_basis_size(m_hdi.cell_degree(), dimension, dimension);
        grad_basis_size = disk::sym_matrix_basis_size(m_hdi.grad_degree(), dimension, dimension);
        num_total_dofs  = 0;

        if (dimension == 2)
            two_dim = true;
        else if (dimension == 3)
            two_dim = false;
        else
            assert(false);
    }

    template<typename Function, typename LawCell>
    void
    compute(const cell_type&   cl,
            const Function&    load,
            const matrix_type& ET,
            const matrix_type& ST,
            const vector_type& uTF,
            const bool&        has_vector_face,
            LawCell&           law)
    {
        timecounter tc;

        const auto ci  = ci_type(m_msh, cl, m_hdi, m_bnd);
        num_total_dofs = ci.num_total_dofs();

        assert(uTF.size() == num_total_dofs);

        const auto cb = make_vector_monomial_basis(m_msh, cl, ci.cell_degree());

        const auto dim_dofs = num_dofs_dim();

        RTF   = vector_type::Zero(num_total_dofs);
        F_int = vector_type::Zero(num_total_dofs);
        K_int = matrix_type::Zero(num_total_dofs, num_total_dofs);

        matrix_type AT = matrix_type::Zero(grad_basis_size, grad_basis_size);
        vector_type aT = vector_type::Zero(grad_basis_size);

        const vector_type GsT_uTF = ET * uTF;

        assert(GsT_uTF.size() == grad_basis_size);

        auto& law_quadpoints = law.getQPs();

        auto gb = disk::make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());

        for (auto& qp : law_quadpoints)
        {
            // std::cout << "qp: " << qp.point() << std::endl;
            const auto gphi = gb.eval_functions(qp.point());

            assert(gphi.size() == grad_basis_size);

            // Compute local gradient and norm
            // std::cout << "GT_utf: " << GsT_uTF << std::endl;
            const auto GsT_iqn = disk::eval(GsT_uTF, gphi);
            // std::cout << "Em" << std::endl;
            // std::cout << qp.getTotalStrainPrev() << std::endl;
            // std::cout << "dE" << std::endl;
            // std::cout << incr_strain << std::endl;

            // Compute bahavior
            tc.tic();
            const auto tensor_behavior = qp.compute_whole(GsT_iqn, m_material_data, true);
            tc.toc();

            //  std::cout << "module " << tensor_behavior.second << std::endl;

            for (int j = 0; j < grad_basis_size; j++)
            {
                const matrix_static Agphi_j = qp.weight() * disk::tm_prod(tensor_behavior.second, gphi[j]);
                // std::cout << j << std::endl;
                // std::cout << gphi[j] << std::endl;
                // std::cout << Agphi_j << std::endl;
                for (int i = 0; i <= j; i += dim_dofs)
                {
                    // compute (Gkt v, A(u) : Gkt du)
                    if (two_dim)
                    {
                        AT(i, j) += Agphi_j(0, 0) * gphi[i](0, 0);
                        AT(i + 1, j) += 2 * Agphi_j(0, 1) * gphi[i + 1](0, 1);
                        AT(i + 2, j) += Agphi_j(1, 1) * gphi[i + 2](1, 1);
                    }
                    else
                    {
                        AT(i, j) += Agphi_j(0, 0) * gphi[i](0, 0);
                        AT(i + 1, j) += 2 * Agphi_j(0, 1) * gphi[i + 1](0, 1);
                        AT(i + 2, j) += Agphi_j(1, 1) * gphi[i + 2](1, 1);
                        AT(i + 3, j) += 2 * Agphi_j(0, 2) * gphi[i + 3](0, 2);
                        AT(i + 4, j) += 2 * Agphi_j(1, 2) * gphi[i + 4](1, 2);
                        AT(i + 5, j) += Agphi_j(2, 2) * gphi[i + 5](2, 2);
                    }
                    // AT(i, j) += disk::mm_prod(gphi[i], Agphi_j);
                }
            }

            // compute (PK1(u), G^k_T v)_T
            const auto stress_qp = disk::priv::inner_product(qp.weight(), tensor_behavior.first);
            //  std::cout << "stress" << std::endl;
            //  std::cout << tensor_behavior.first << std::endl;

            for (int i = 0; i < grad_basis_size; i += dim_dofs)
            {
                if (two_dim)
                {
                    aT(i) += stress_qp(0, 0) * gphi[i](0, 0);
                    aT(i + 1) += 2 * stress_qp(0, 1) * gphi[i + 1](0, 1);
                    aT(i + 2) += stress_qp(1, 1) * gphi[i + 2](1, 1);
                }
                else
                {
                    aT(i) += stress_qp(0, 0) * gphi[i](0, 0);
                    aT(i + 1) += 2 * stress_qp(0, 1) * gphi[i + 1](0, 1);
                    aT(i + 2) += stress_qp(1, 1) * gphi[i + 2](1, 1);
                    aT(i + 3) += 2 * stress_qp(0, 2) * gphi[i + 3](0, 2);
                    aT(i + 4) += 2 * stress_qp(1, 2) * gphi[i + 4](1, 2);
                    aT(i + 5) += stress_qp(2, 2) * gphi[i + 5](2, 2);
                }
            }
        }

        // lower part AT
        for (int j = 0; j < grad_basis_size; j++)
            for (int i = j; i < grad_basis_size; i++)
                AT(i, j) = AT(j, i);

        // add volumic term
        RTF.head(cell_basis_size) = make_rhs(m_msh, cl, cb, load, 2);

        // contact contribution
        time_contact = 0.0;
        if (has_vector_face)
        {
            tc.tic();
            const vector_type stress_coeff = law.projectStressSymOnCell(m_msh, cl, m_hdi, m_material_data);

            // compute theta/gamma *(sigma_n, sigma_n)_Fc
            if (m_rp.m_theta != 0.0)
            {
                const matrix_type K_sig = make_hho_nitsche(cl, ET, ci);
                K_int -= K_sig;
                F_int -= K_sig * uTF;
            }

            //  std::cout << "Nitche: " << std::endl;
            // std::cout << make_hho_nitsche(cl, ET) << std::endl;

            // compute (phi_n_theta, H(-phi_n_1(u))*phi_n_1)_FC / gamma
            K_int += make_hho_heaviside_contact(cl, ET, stress_coeff, uTF, ci);

            //  std::cout << "Heaviside: " << std::endl;
            // std::cout << make_hho_heaviside_contact(cl, ET, uTF) << std::endl;

            // compute (phi_n_theta, [phi_n_1(u)]R-)_FC / gamma
            F_int += make_hho_negative_contact(cl, ET, stress_coeff, uTF, ci);

            //  std::cout << "Negative: " << std::endl;
            // std::cout << make_hho_negative_contact(cl, ET, uTF).transpose() << std::endl;

            // friction contribution
            if (m_rp.m_frot)
            {
                // compute (phi_t_theta, [phi_t_1(u)]_s)_FC / gamma
                F_int += make_hho_threshold_tresca(cl, ET, stress_coeff, uTF, ci);

                // std::cout << "Threshold: " << std::endl;
                // std::cout << make_hho_threshold_tresca(cl, ET, uTF).transpose() << std::endl;

                // compute (phi_t_theta, (d_proj_alpha(u)) phi_t_1)_FC / gamma
                K_int += make_hho_matrix_tresca(cl, ET, stress_coeff, uTF, ci);
            }

            tc.toc();
            time_contact = tc.to_double();
        }

        K_int += ET.transpose() * AT * ET;
        F_int += ET.transpose() * aT;
        if (m_rp.m_stab)
        {
            K_int += m_rp.m_beta * ST;
            F_int += m_rp.m_beta * ST * uTF;
        }

        RTF -= F_int;

        //  std::cout << "K: " << K_int.norm() << std::endl;
        //  std::cout << K_int << std::endl;
        //  std::cout << "R: " << RTF.norm() << std::endl;
        //  std::cout << RTF.transpose() << std::endl;

        assert(K_int.rows() == num_total_dofs);
        assert(K_int.cols() == num_total_dofs);
        assert(RTF.rows() == num_total_dofs);
    }

    template<typename CellBasis>
    scalar_type
    eval_uT_n(const CellBasis& cb, const vector_type& uTF, const vector_static& n, const point_type& pt) const
    {
        const vector_type uT_n = make_hho_u_n(n, cb, pt);

        return uT_n.dot(uTF.head(cell_basis_size));
    }

    template<typename CellBasis>
    vector_static
    eval_uT_t(const CellBasis& cb, const vector_type& uTF, const vector_static& n, const point_type& pt) const
    {
        const auto uT_t = make_hho_u_t(n, cb, pt);

        return uT_t.transpose() * (uTF.head(cell_basis_size));
    }

    template<typename FaceBasis>
    scalar_type
    eval_uF_n(const FaceBasis& fb, const vector_type& uF, const vector_static& n, const point_type& pt) const
    {
        const vector_type uF_n = make_hho_u_n(n, fb, pt);
        assert(uF_n.size() == uF.size());
        return uF_n.dot(uF);
    }

    template<typename FaceBasis>
    vector_static
    eval_uF_t(const FaceBasis& fb, const vector_type& uF, const vector_static& n, const point_type& pt) const
    {
        const auto uF_t = make_hho_u_t(n, fb, pt);
        assert(uF_t.rows() == uF.size());
        return uF_t.transpose() * uF;
    }

    template<typename GradBasis>
    scalar_type
    eval_stress_nn(const vector_type&   stress_coeff,
                   const GradBasis&     gb,
                   const vector_static& n,
                   const point_type&    pt) const
    {
        const auto stress = eval_stress(stress_coeff, gb, pt);
        return (stress * n).dot(n);
    }

    template<typename GradBasis>
    vector_static
    eval_stress_nt(const vector_type&   stress_coeff,
                   const GradBasis&     gb,
                   const vector_static& n,
                   const point_type&    pt) const
    {
        const auto s_n  = eval_stress(stress_coeff, gb, pt) * n;
        const auto s_nn = s_n.dot(n);
        return s_n - s_nn * n;
    }

    template<typename GradBasis, typename CellBasis>
    scalar_type
    eval_phi_n_uT(const vector_type&   stress_coeff,
                  const GradBasis&     gb,
                  const CellBasis&     cb,
                  const vector_type&   uTF,
                  const vector_static& n,
                  const scalar_type&   gamma_F,
                  const point_type&    pt) const
    {
        const scalar_type sigma_nn = eval_stress_nn(stress_coeff, gb, n, pt);
        const scalar_type uT_n     = eval_uT_n(cb, uTF, n, pt);
        const scalar_type g0       = make_hho_distance(pt, n);

        return sigma_nn - gamma_F * (uT_n - g0);
    }

    template<typename GradBasis, typename CellBasis>
    scalar_type
    eval_proj_phi_n_uT(const vector_type&   stress_coeff,
                       const GradBasis&     gb,
                       const CellBasis&     cb,
                       const vector_type&   uTF,
                       const vector_static& n,
                       const scalar_type&   gamma_F,
                       const point_type&    pt) const
    {
        const scalar_type phi_n_1_u = eval_phi_n_uT(stress_coeff, gb, cb, uTF, n, gamma_F, pt);

        if (phi_n_1_u <= scalar_type(0))
            return phi_n_1_u;

        return scalar_type(0);
    }

    template<typename GradBasis, typename CellBasis>
    vector_static
    eval_phi_t_uT(const vector_type&   stress_coeff,
                  const GradBasis&     gb,
                  const CellBasis&     cb,
                  const vector_type&   uTF,
                  const vector_static& n,
                  const scalar_type&   gamma_F,
                  const point_type&    pt) const
    {
        const auto sigma_nt = eval_stress_nt(stress_coeff, gb, n, pt);
        const auto uT_t     = eval_uT_t(cb, uTF, n, pt);

        return sigma_nt - gamma_F * uT_t;
    }

    template<typename GradBasis, typename CellBasis>
    vector_static
    eval_proj_phi_t_uT(const vector_type&   stress_coeff,
                       const GradBasis&     gb,
                       const CellBasis&     cb,
                       const vector_type&   uTF,
                       const vector_static& n,
                       const scalar_type&   gamma_F,
                       const scalar_type&   s,
                       const point_type&    pt) const
    {
        const vector_static phi_t_1_u = eval_phi_t_uT(stress_coeff, gb, cb, uTF, n, gamma_F, pt);

        return make_proj_alpha(phi_t_1_u, s);
    }

    template<typename GradBasis, typename FaceBasis>
    scalar_type
    eval_phi_n_uF(const face_type&     fc,
                  const vector_type&   stress_coeff,
                  const GradBasis&     gb,
                  const FaceBasis&     fb,
                  const vector_type&   uTF,
                  const size_t         offset,
                  const vector_static& n,
                  const scalar_type&   gamma_F,
                  const point_type&    pt) const
    {
        const vector_type uF       = uTF.segment(offset, fb.size());
        const scalar_type sigma_nn = eval_stress_nn(stress_coeff, gb, n, pt);
        const scalar_type uF_n     = eval_uF_n(fb, uF, n, pt);
        const scalar_type g0       = make_hho_distance(pt, n);
        const scalar_type gap      = priv::compute_gap_fb(m_msh, fc, fb, uF, pt, n);

        // std::cout << gap << std::endl;

        return sigma_nn + gamma_F * gap;
    }

    template<typename GradBasis, typename FaceBasis>
    scalar_type
    eval_proj_phi_n_uF(const face_type&     fc,
                       const vector_type&   stress_coeff,
                       const GradBasis&     gb,
                       const FaceBasis&     fb,
                       const vector_type&   uTF,
                       const size_t         offset,
                       const vector_static& n,
                       const scalar_type&   gamma_F,
                       const point_type&    pt) const
    {
        const scalar_type phi_n_1_u = eval_phi_n_uF(fc, stress_coeff, gb, fb, uTF, offset, n, gamma_F, pt);

        if (phi_n_1_u <= scalar_type(0))
            return phi_n_1_u;

        return scalar_type(0);
    }

    template<typename GradBasis, typename FaceBasis>
    vector_static
    eval_phi_t_uF(const vector_type&   stress_coeff,
                  const GradBasis&     gb,
                  const FaceBasis&     fb,
                  const vector_type&   uTF,
                  const size_t&        offset,
                  const vector_static& n,
                  const scalar_type&   gamma_F,
                  const point_type&    pt) const
    {
        const vector_type uF       = uTF.segment(offset, fb.size());
        const auto        sigma_nt = eval_stress_nt(stress_coeff, gb, n, pt);
        const auto        uF_t     = eval_uF_t(fb, uF, n, pt);

        return sigma_nt - gamma_F * uF_t;
    }

    template<typename GradBasis, typename FaceBasis>
    vector_static
    eval_proj_phi_t_uF(const vector_type&   stress_coeff,
                       const GradBasis&     gb,
                       const FaceBasis&     fb,
                       const vector_type&   uTF,
                       const size_t&        offset,
                       const vector_static& n,
                       const scalar_type&   gamma_F,
                       const scalar_type&   s,
                       const point_type&    pt) const
    {
        const vector_static phi_t_1_u = eval_phi_t_uF(stress_coeff, gb, fb, uTF, offset, n, gamma_F, pt);

        return make_proj_alpha(phi_t_1_u, s);
    }

    template<typename GradBasis>
    matrix_static
    eval_stress(const vector_type& stress_coeff, const GradBasis& gb, const point_type pt) const
    {
        const auto gphi = gb.eval_functions(pt);

        return disk::eval(stress_coeff, gphi);
    }
};

} // end namespace MK