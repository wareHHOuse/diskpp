/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018, 2019
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
 * Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework.
 * M. Abbas, A. Ern, N. Pignet.
 * International Journal of Numerical Methods in Engineering (2019)
 * 120(3), 303-327
 * DOI: 10.1002/nme.6137
 */

#pragma once

#include <cassert>

#include "bases/bases.hpp"
#include "common/eigen.hpp"
#include "mechanics/behaviors/maths_tensor.hpp"
#include "methods/hho"
#include "quadratures/quadratures.hpp"

#include "timecounter.h"

namespace disk
{

namespace mechanics
{

namespace priv
{

template<typename T>
T
compute_g0_gv(const point<T, 3>& pt, const static_vector<T, 3>& n)
{
    // distance to the cylinder z^2+y^2 = r^2 (axe x)

    const T r = 8.77;

    // eq in a* t^2 + b *t + c =0
    const T a = n(1) * n(1) + n(2) * n(2);
    const T b = 2 * (n(1) * pt.y() + n(2) * pt.z());
    const T c = pt.y() * pt.y() + pt.z() * pt.z() - r * r;

    const T delta = b * b - 4.0 * a * c;

    if (abs(delta) <= 1E-12)
    {
        throw std::invalid_argument("wrong prjoection for GV");
    }

    const T t1 = (-b + sqrt(delta)) / (2.0 * a);
    const T t2 = (-b - sqrt(delta)) / (2.0 * a);

    const static_vector<T, 3> p_0 = pt.to_vector();
    const static_vector<T, 3> p_1 = p_0 + t1 * n;
    const static_vector<T, 3> p_2 = p_0 + t2 * n;

    const T gap_1 = (p_1 - p_0).dot(n);
    const T gap_2 = (p_2 - p_0).dot(n);

    //    std::cout << "pt : " << pt << std::endl;
    //    std::cout << "n : " << n.transpose() << std::endl;
    //    std::cout << "p1 : " << p_1.transpose() << std::endl;
    //    std::cout << "p2 : " << p_2.transpose() << std::endl;
    //    std::cout << sqrt(p_0(1)*p_0(1) + p_0(2)*p_0(2)) << " " << sqrt(p_1(1)*p_1(1) + p_1(2)*p_1(2)) << " " <<
    //    sqrt(p_2(1)*p_2(1) + p_2(2)*p_2(2)) << std::endl; std::cout << (p_1 - p_0).norm() << " " << (p_2 - p_0).norm()
    //    << std::endl; std::cout << gap_1 << " " << gap_2 << std::endl;

    if (abs(gap_1) < abs(gap_2))
    {
        return gap_1;
    }

    return gap_2;
}

template<typename T>
T
compute_g0(const point<T, 2>& pt, const static_vector<T, 2>& n)
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
compute_g0(const point<T, 3>& pt, const static_vector<T, 3>& n)
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
new_pt(const Basis& base, const dynamic_vector<T>& tab_coeff, const point<T, DIM>& pt)
{
    const auto beval = base.eval_functions(pt);

    const auto depl = disk::eval(tab_coeff, beval);

    return pt + depl;
}

template<typename T>
static_vector<T, 2>
compute_normal(const point<T, 2>& a, const point<T, 2>& b)
{
    const static_vector<T, 2> t1 = (b - a).to_vector();
    static_vector<T, 2>       nor;
    nor(0) = -t1(1);
    nor(1) = t1(0);

    return nor / nor.norm();
}

template<typename T>
static_vector<T, 3>
compute_normal(const point<T, 3>& p1, const point<T, 3>& p2, const point<T, 3>& p3)
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
               const dynamic_vector<T>&   tab_coeff,
               const point<T, 2>&         pt,
               const static_vector<T, 2>& n)
{
    const point<T, 2>         pt_def = new_pt(eb, tab_coeff, pt);
    const auto                pts    = points(msh, elem);
    const static_vector<T, 2> n_ref  = compute_normal(pts[0], pts[1]);

    const point<T, 2> pta_def = new_pt(eb, tab_coeff, pts[0]);
    const point<T, 2> ptb_def = new_pt(eb, tab_coeff, pts[1]);

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
               const dynamic_vector<T>&   tab_coeff,
               const point<T, 3>&         pt,
               const static_vector<T, 3>& n)
{
    const point<T, 3>         pt_def = new_pt(eb, tab_coeff, pt);
    const auto                pts    = points(msh, elem);
    const static_vector<T, 3> n_ref  = compute_normal(pts[0], pts[1], pts[2]);

    const point<T, 3> pt0_def = new_pt(eb, tab_coeff, pts[0]);
    const point<T, 3> pt1_def = new_pt(eb, tab_coeff, pts[1]);
    const point<T, 3> pt2_def = new_pt(eb, tab_coeff, pts[2]);

    const T sign = std::copysign(T(1), n.dot(n_ref));

    const static_vector<T, 3> n_def = sign * compute_normal(pt0_def, pt1_def, pt2_def);

    return compute_g0_gv(pt_def, n_def);
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



    scalar_type
    make_hho_distance(const point_type& pt, const vector_static& n) const
    {
        return priv::compute_g0(pt, n);
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
    double      time_contact, time_law;

    tresca(const mesh_type&     msh,
           const hdi_type&      hdi,
           const material_type& material_data,
           const param_type&    rp,
           const bnd_type&      bnd) :
      m_msh(msh),
      m_hdi(hdi), m_material_data(material_data), m_rp(rp), m_bnd(bnd), time_contact(0.0), time_law(0.0)
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

        time_contact = 0.0;
        time_law     = 0.0;

        // contact contribution
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
        // return sigma_nn - gamma_F * (uF_n - g0);
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

};
}