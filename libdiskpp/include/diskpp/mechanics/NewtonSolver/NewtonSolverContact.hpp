/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2024                     nicolas.pignet@enpc.fr
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

#include "diskpp/bases/bases.hpp"
#include "diskpp/boundary_conditions/boundary_conditions.hpp"
#include "diskpp/common/eigen.hpp"
#include "diskpp/common/timecounter.hpp"
#include "diskpp/mechanics/NewtonSolver/NewtonSolverParameters.hpp"
#include "diskpp/mechanics/behaviors/laws/materialData.hpp"
#include "diskpp/mechanics/behaviors/maths_tensor.hpp"
#include "diskpp/methods/hho"
#include "diskpp/quadratures/quadratures.hpp"

namespace disk
{

namespace mechanics
{

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

    static_vector<T, 3> nor = t1.cross(t2);

    return nor / nor.norm();
}

template<typename Mesh, typename Elem, typename ElemBasis, typename T>
T
compute_gap_fb(const Mesh&                    msh,
               const Elem&                    elem,
               const ElemBasis&               eb,
               const disk::dynamic_vector<T>& tab_coeff,
               const disk::point<T, 2>&       pt,
               const static_vector<T, 2>&     n)
{
    const disk::point<T, 2>   pt_def = new_pt(eb, tab_coeff, pt);
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
compute_gap_fb(const Mesh&                    msh,
               const Elem&                    elem,
               const ElemBasis&               eb,
               const disk::dynamic_vector<T>& tab_coeff,
               const disk::point<T, 3>&       pt,
               const static_vector<T, 3>&     n)
{
    const disk::point<T, 3>   pt_def = new_pt(eb, tab_coeff, pt);
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
class contact_contribution
{
  private:
    typedef MeshType                                    mesh_type;
    typedef typename mesh_type::coordinate_type         scalar_type;
    typedef typename mesh_type::cell                    cell_type;
    typedef typename mesh_type::face                    face_type;
    typedef point<scalar_type, mesh_type::dimension>    point_type;
    typedef disk::MaterialData<scalar_type>             material_type;
    typedef NewtonSolverParameter<scalar_type>          param_type;
    typedef disk::vector_boundary_conditions<mesh_type> bnd_type;

    const static int dimension = mesh_type::dimension;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    typedef static_matrix<scalar_type, dimension, dimension> matrix_static;
    typedef static_vector<scalar_type, dimension>            vector_static;

    typedef Matrix<scalar_type, Dynamic, dimension> matrix_d_type;

    const mesh_type&     m_msh;
    const material_type& m_material_data;
    const param_type&    m_rp;
    const bnd_type&      m_bnd;

    // contact contrib;

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

        // std::cout << "t_phi: " << t_phi.transpose() << std::endl;

        //(phi_T . n)
        return disk::priv::inner_product(t_phi, n);
    }

    // tangential part of u : u_t = u - u_n*n
    template<typename TraceBasis>
    matrix_d_type
    make_hho_u_t(const vector_static& n, const TraceBasis& tb, const point_type& pt) const
    {
        const auto t_phi = tb.eval_functions(pt);
        const auto u_n   = make_hho_u_n(n, tb, pt);

        // phi_T - (phi_T . n)n
        return t_phi - disk::priv::inner_product(u_n, n);
    }

    // cauchy traction : sigma_n = sigma * n
    template<typename GradBasis>
    matrix_d_type
    make_hho_sigma_n(const matrix_type& ET, const vector_static& n, const GradBasis& gb, const point_type& pt) const
    {
        matrix_d_type sigma_n = matrix_d_type::Zero(ET.cols(), dimension);

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
    make_hho_sigma_nn(const matrix_d_type& sigma_n, const vector_static& n) const
    {
        // sigma_n . n
        return disk::priv::inner_product(sigma_n, n);
    }

    template<typename GradBasis>
    matrix_d_type
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
        phi_n.head(uT_n.size()) -= gamma_F * uT_n;

        // theta * sigma_nn - gamma uT_n
        return phi_n;
    }

    matrix_d_type
    make_hho_phi_t_uT(const matrix_d_type& sigma_nt,
                      const matrix_d_type& uT_t,
                      const scalar_type&   theta,
                      const scalar_type&   gamma_F) const
    {
        matrix_d_type phi_t = theta * sigma_nt;
        phi_t.block(0, 0, uT_t.rows(), dimension) -= gamma_F * uT_t;

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

    matrix_d_type
    make_hho_phi_t_uF(const matrix_d_type& sigma_nt,
                      const matrix_d_type& uF_t,
                      const scalar_type&   theta,
                      const scalar_type&   gamma_F,
                      const size_t&        offset) const
    {
        matrix_d_type phi_t = theta * sigma_nt;

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
    make_hho_nitsche(const cell_type& cl, const matrix_type& ET, const CellDegreeInfo<MeshType>& cell_infos) const
    {
        const auto gb = make_sym_matrix_monomial_basis(m_msh, cl, cell_infos.grad_degree());

        matrix_type nitsche = matrix_type::Zero(ET.cols(), ET.cols());

        const auto fcs = m_bnd.faces_with_contact(cl);
        for (auto& fc : fcs)
        {
            const auto n       = normal(m_msh, cl, fc);
            const auto qps     = integrate(m_msh, fc, 2 * cell_infos.grad_degree() + 2);
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
    make_hho_heaviside_contact(const cell_type&                cl,
                               const matrix_type&              ET,
                               const vector_type&              uTF,
                               const CellDegreeInfo<MeshType>& cell_infos) const
    {
        const auto cb = make_vector_monomial_basis(m_msh, cl, cell_infos.cell_degree());
        const auto gb = make_sym_matrix_monomial_basis(m_msh, cl, cell_infos.grad_degree());

        matrix_type lhs = matrix_type::Zero(uTF.size(), uTF.size());

        const auto fcs    = faces(m_msh, cl);
        size_t     offset = cb.size();

        const auto fcs_di = cell_infos.facesDegreeInfo();
        size_t     face_i = 0;

        const vector_type ET_uTF = ET * uTF;

        for (auto& fc : fcs)
        {
            const auto fdi     = fcs_di[face_i++];
            const auto facedeg = fdi.degree();
            const auto fb      = make_vector_monomial_basis(m_msh, fc, facedeg);
            const auto fbs     = fb.size();

            if (m_bnd.is_contact_face(fc))
            {
                const auto contact_type = m_bnd.contact_boundary_type(fc);
                const auto n            = normal(m_msh, cl, fc);
                const auto qp_deg       = std::max(cell_infos.cell_degree(), cell_infos.grad_degree());
                const auto qps          = integrate(m_msh, fc, 2 * qp_deg + 2);
                const auto hF           = diameter(m_msh, fc);
                const auto gamma_F      = m_rp.m_gamma_0 / hF;

                for (auto& qp : qps)
                {
                    const vector_type sigma_nn = make_hho_sigma_nn(ET, n, gb, qp.point());

                    if (contact_type == disk::SIGNORINI_CELL)
                    {
                        const vector_type uT_n = make_hho_u_n(n, cb, qp.point());

                        const scalar_type phi_n_1_u = eval_phi_n_uT(ET_uTF, gb, cb, uTF, n, gamma_F, qp.point());

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
                          eval_phi_n_uF(fc, ET_uTF, gb, fb, uTF, offset, n, gamma_F, qp.point());

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
    make_hho_negative_contact(const cell_type&                cl,
                              const matrix_type&              ET,
                              const vector_type&              uTF,
                              const CellDegreeInfo<MeshType>& cell_infos) const
    {
        const auto cb = make_vector_monomial_basis(m_msh, cl, cell_infos.cell_degree());
        const auto gb = make_sym_matrix_monomial_basis(m_msh, cl, cell_infos.grad_degree());

        vector_type rhs = vector_type::Zero(uTF.size());

        const auto fcs    = faces(m_msh, cl);
        size_t     offset = cb.size();

        const auto fcs_di = cell_infos.facesDegreeInfo();
        size_t     face_i = 0;

        const vector_type ET_uTF = ET * uTF;

        for (auto& fc : fcs)
        {
            const auto fdi     = fcs_di[face_i++];
            const auto facedeg = fdi.degree();
            const auto fb      = make_vector_monomial_basis(m_msh, fc, facedeg);
            const auto fbs     = fb.size();

            if (m_bnd.is_contact_face(fc))
            {
                const auto contact_type = m_bnd.contact_boundary_type(fc);
                const auto n            = normal(m_msh, cl, fc);
                const auto qp_deg       = std::max(cell_infos.cell_degree(), cell_infos.grad_degree());
                const auto qps          = integrate(m_msh, fc, 2 * qp_deg + 2);
                const auto hF           = diameter(m_msh, fc);
                const auto gamma_F      = m_rp.m_gamma_0 / hF;

                for (auto& qp : qps)
                {
                    const vector_type sigma_nn = make_hho_sigma_nn(ET, n, gb, qp.point());

                    if (contact_type == disk::SIGNORINI_CELL)
                    {
                        const vector_type uT_n = make_hho_u_n(n, cb, qp.point());

                        const scalar_type phi_n_1_u = eval_phi_n_uT(ET_uTF, gb, cb, uTF, n, gamma_F, qp.point());

                        // [phi_n_1_u]_R-
                        if (phi_n_1_u <= scalar_type(0))
                        {
                            const vector_type phi_n_theta = make_hho_phi_n_uT(sigma_nn, uT_n, m_rp.m_theta, gamma_F);

                            rhs += (qp.weight() / gamma_F * phi_n_1_u) * phi_n_theta;
                        }
                    }
                    else
                    {
                        const scalar_type phi_n_1_u =
                          eval_phi_n_uF(fc, ET_uTF, gb, fb, uTF, offset, n, gamma_F, qp.point());

                        // std::cout << "qp: " << qp.point() << std::endl;
                        // std::cout << "phi_n_1_u: " << phi_n_1_u << std::endl;

                        // [phi_n_1_u]_R-
                        if (phi_n_1_u <= scalar_type(0))
                        {
                            const vector_type uF_n = make_hho_u_n(n, fb, qp.point());
                            const vector_type phi_n_theta =
                              make_hho_phi_n_uF(sigma_nn, uF_n, m_rp.m_theta, gamma_F, offset);

                            // std::cout << "sigma_nn: " << sigma_nn.transpose() << std::endl;
                            // std::cout << "uF_n: " << uF_n.transpose() << std::endl;
                            // std::cout << "phi_n_theta: " << phi_n_theta.transpose() << std::endl;

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
    make_hho_threshold_tresca(const cell_type&                cl,
                              const matrix_type&              ET,
                              const vector_type&              uTF,
                              const CellDegreeInfo<MeshType>& cell_infos) const
    {
        const auto cb = make_vector_monomial_basis(m_msh, cl, cell_infos.cell_degree());
        const auto gb = make_sym_matrix_monomial_basis(m_msh, cl, cell_infos.grad_degree());

        vector_type rhs = vector_type::Zero(uTF.size());

        const auto fcs    = faces(m_msh, cl);
        size_t     offset = cb.size();

        const auto fcs_di = cell_infos.facesDegreeInfo();
        size_t     face_i = 0;

        const vector_type ET_uTF = ET * uTF;

        for (auto& fc : fcs)
        {
            const auto fdi     = fcs_di[face_i++];
            const auto facedeg = fdi.degree();
            const auto fb      = make_vector_monomial_basis(m_msh, fc, facedeg);
            const auto fbs     = fb.size();

            if (m_bnd.is_contact_face(fc))
            {
                const auto contact_type = m_bnd.contact_boundary_type(fc);
                const auto n            = normal(m_msh, cl, fc);
                const auto qp_deg       = std::max(cell_infos.cell_degree(), cell_infos.grad_degree());
                const auto qps          = integrate(m_msh, fc, 2 * qp_deg + 2);
                const auto hF           = diameter(m_msh, fc);
                const auto gamma_F      = m_rp.m_gamma_0 / hF;

                const auto s_func = m_bnd.contact_boundary_func(fc);

                for (auto& qp : qps)
                {
                    const auto sigma_nt = make_hho_sigma_nt(ET, n, gb, qp.point());

                    if (contact_type == disk::SIGNORINI_CELL)
                    {
                        const auto uT_t = make_hho_u_t(n, cb, qp.point());

                        const auto phi_t_theta = make_hho_phi_t_uT(sigma_nt, uT_t, m_rp.m_theta, gamma_F);

                        const vector_static phi_t_1_u_proj =
                          eval_proj_phi_t_uT(ET_uTF, gb, cb, uTF, n, gamma_F, s_func(qp.point()), qp.point());

                        const vector_static qp_phi_t_1_u_pro = qp.weight() * phi_t_1_u_proj / gamma_F;

                        rhs += disk::priv::inner_product(phi_t_theta, qp_phi_t_1_u_pro);
                    }
                    else
                    {
                        const auto uF_t = make_hho_u_t(n, fb, qp.point());

                        const auto phi_t_theta = make_hho_phi_t_uF(sigma_nt, uF_t, m_rp.m_theta, gamma_F, offset);

                        // std::cout << "sigma_nt: " << sigma_nt.transpose() << std::endl;
                        //                         std::cout << "uF_t: " << uF_t.transpose() << std::endl;
                        //                         std::cout << "phi_t_theta: " << phi_t_theta.transpose() << std::endl;

                        const vector_static phi_t_1_u_proj = eval_proj_tresca_phi_t_uF(
                          ET_uTF, gb, fb, uTF, offset, n, gamma_F, s_func(qp.point()), qp.point());

                        const vector_static qp_phi_t_1_u_pro = qp.weight() * phi_t_1_u_proj / gamma_F;

                        // std::cout << "phi_t_1_u_proj: " << phi_t_1_u_proj.transpose() << std::endl;

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
    make_hho_matrix_tresca(const cell_type&                cl,
                           const matrix_type&              ET,
                           const vector_type&              uTF,
                           const CellDegreeInfo<MeshType>& cell_infos) const
    {
        const auto cb = make_vector_monomial_basis(m_msh, cl, cell_infos.cell_degree());
        const auto gb = make_sym_matrix_monomial_basis(m_msh, cl, cell_infos.grad_degree());

        matrix_type lhs = matrix_type::Zero(uTF.size(), uTF.size());

        const auto fcs    = faces(m_msh, cl);
        size_t     offset = cb.size();

        const auto fcs_di = cell_infos.facesDegreeInfo();
        size_t     face_i = 0;

        const vector_type ET_uTF = ET * uTF;

        for (auto& fc : fcs)
        {
            const auto fdi     = fcs_di[face_i++];
            const auto facedeg = fdi.degree();
            const auto fb      = make_vector_monomial_basis(m_msh, fc, facedeg);
            const auto fbs     = fb.size();

            if (m_bnd.is_contact_face(fc))
            {
                const auto contact_type = m_bnd.contact_boundary_type(fc);
                const auto n            = normal(m_msh, cl, fc);
                const auto qp_deg       = std::max(cell_infos.cell_degree(), cell_infos.grad_degree());
                const auto qps          = integrate(m_msh, fc, 2 * qp_deg + 2);
                const auto hF           = diameter(m_msh, fc);
                const auto gamma_F      = m_rp.m_gamma_0 / hF;

                const auto s_func = m_bnd.contact_boundary_func(fc);

                for (auto& qp : qps)
                {
                    const auto sigma_nt = make_hho_sigma_nt(ET, n, gb, qp.point());

                    if (contact_type == disk::SIGNORINI_CELL)
                    {
                        const auto uT_t = make_hho_u_t(n, cb, qp.point());

                        const auto phi_t_1     = make_hho_phi_t_uT(sigma_nt, uT_t, scalar_type(1), gamma_F);
                        const auto phi_t_theta = make_hho_phi_t_uT(sigma_nt, uT_t, m_rp.m_theta, gamma_F);

                        const auto phi_t_1_u      = eval_phi_t_uT(ET_uTF, gb, cb, uTF, n, gamma_F, qp.point());
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

                        const auto phi_t_1_u      = eval_phi_t_uF(ET_uTF, gb, fb, uTF, offset, n, gamma_F, qp.point());
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

    // compute (phi_t_theta, [phi_t_1(u)]_(s))_FC / gamma
    vector_type
    make_hho_threshold_coulomb(const cell_type&                cl,
                               const matrix_type&              ET,
                               const vector_type&              uTF,
                               const CellDegreeInfo<MeshType>& cell_infos) const
    {
        const auto cb = make_vector_monomial_basis(m_msh, cl, cell_infos.cell_degree());
        const auto gb = make_sym_matrix_monomial_basis(m_msh, cl, cell_infos.grad_degree());

        vector_type rhs = vector_type::Zero(uTF.size());

        const auto fcs    = faces(m_msh, cl);
        size_t     offset = cb.size();

        const auto fcs_di = cell_infos.facesDegreeInfo();
        size_t     face_i = 0;

        const vector_type ET_uTF = ET * uTF;

        for (auto& fc : fcs)
        {
            const auto fdi     = fcs_di[face_i++];
            const auto facedeg = fdi.degree();
            const auto fb      = make_vector_monomial_basis(m_msh, fc, facedeg);
            const auto fbs     = fb.size();

            if (m_bnd.is_contact_face(fc))
            {
                const auto contact_type = m_bnd.contact_boundary_type(fc);
                const auto n            = normal(m_msh, cl, fc);
                const auto qp_deg       = std::max(cell_infos.cell_degree(), cell_infos.grad_degree());
                const auto qps          = integrate(m_msh, fc, 2 * qp_deg + 2);
                const auto hF           = diameter(m_msh, fc);
                const auto gamma_F      = m_rp.m_gamma_0 / hF;

                const auto s_func = m_bnd.contact_boundary_func(fc);

                for (auto& qp : qps)
                {
                    const auto sigma_nt = make_hho_sigma_nt(ET, n, gb, qp.point());

                    if (contact_type == disk::SIGNORINI_CELL)
                    {
                        throw std::runtime_error("Not implemented");
                    }
                    else
                    {
                        const auto uF_t = make_hho_u_t(n, fb, qp.point());

                        const auto phi_t_theta = make_hho_phi_t_uF(sigma_nt, uF_t, m_rp.m_theta, gamma_F, offset);

                        const vector_static phi_t_1_u_proj = eval_proj_coulomb_phi_t_uF(
                          fc, ET_uTF, gb, fb, uTF, offset, n, gamma_F, s_func(qp.point()), qp.point());

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
    make_hho_matrix_coulomb(const cell_type&                cl,
                            const matrix_type&              ET,
                            const vector_type&              uTF,
                            const CellDegreeInfo<MeshType>& cell_infos) const
    {
        const auto cb = make_vector_monomial_basis(m_msh, cl, cell_infos.cell_degree());
        const auto gb = make_sym_matrix_monomial_basis(m_msh, cl, cell_infos.grad_degree());

        matrix_type lhs = matrix_type::Zero(uTF.size(), uTF.size());

        const auto fcs    = faces(m_msh, cl);
        size_t     offset = cb.size();

        const auto fcs_di = cell_infos.facesDegreeInfo();
        size_t     face_i = 0;

        const vector_type ET_uTF = ET * uTF;

        for (auto& fc : fcs)
        {
            const auto fdi     = fcs_di[face_i++];
            const auto facedeg = fdi.degree();
            const auto fb      = make_vector_monomial_basis(m_msh, fc, facedeg);
            const auto fbs     = fb.size();

            if (m_bnd.is_contact_face(fc))
            {
                const auto contact_type = m_bnd.contact_boundary_type(fc);
                const auto n            = normal(m_msh, cl, fc);
                const auto qp_deg       = std::max(cell_infos.cell_degree(), cell_infos.grad_degree());
                const auto qps          = integrate(m_msh, fc, 2 * qp_deg + 2);
                const auto hF           = diameter(m_msh, fc);
                const auto gamma_F      = m_rp.m_gamma_0 / hF;

                const auto s_func = m_bnd.contact_boundary_func(fc);

                for (auto& qp : qps)
                {
                    const auto sigma_nt = make_hho_sigma_nt(ET, n, gb, qp.point());

                    if (contact_type == disk::SIGNORINI_CELL)
                    {
                        assert(false);
                    }
                    else
                    {
                        const auto uF_t = make_hho_u_t(n, fb, qp.point());

                        const auto phi_t_1     = make_hho_phi_t_uF(sigma_nt, uF_t, scalar_type(1), gamma_F, offset);
                        const auto phi_t_theta = make_hho_phi_t_uF(sigma_nt, uF_t, m_rp.m_theta, gamma_F, offset);

                        const auto phi_t_1_u = eval_phi_t_uF(ET_uTF, gb, fb, uTF, offset, n, gamma_F, qp.point());
                        const scalar_type phi_n_1_u =
                          eval_phi_n_uF(fc, ET_uTF, gb, fb, uTF, offset, n, gamma_F, qp.point());

                        const scalar_type proj_phi_n_1_u = std::min(scalar_type(0), phi_n_1_u);
                        const scalar_type fric_bound     = -s_func(qp.point()) * proj_phi_n_1_u;
                        const auto        d_proj_phi_t_u = make_d_proj_alpha(phi_t_1_u, fric_bound);

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
    matrix_type K_cont;
    vector_type F_cont;
    double      time_contact;

    contact_contribution(const mesh_type&     msh,
                         const material_type& material_data,
                         const param_type&    rp,
                         const bnd_type&      bnd) :
      m_msh(msh),
      m_material_data(material_data), m_rp(rp), m_bnd(bnd)
    {
    }

    void
    compute(const cell_type&                 cl,
            const CellDegreeInfo<mesh_type>& cell_infos,
            const matrix_type&               ET,
            const vector_type&               uTF)
    {
        timecounter tc;
        tc.tic();

        // contact contribution
        time_contact = 0.0;

        const auto cb = disk::make_vector_monomial_basis(m_msh, cl, cell_infos.cell_degree());
        const auto gb = disk::make_sym_matrix_monomial_basis(m_msh, cl, cell_infos.grad_degree());

        F_cont = vector_type::Zero(uTF.size());
        K_cont = matrix_type::Zero(uTF.size(), uTF.size());

        // compute theta/gamma *(sigma_n, sigma_n)_Fc
        if (m_rp.m_theta != 0.0)
        {
            const matrix_type K_sig = make_hho_nitsche(cl, ET, cell_infos);
            K_cont -= K_sig;
            F_cont -= K_sig * uTF;
        }

        // std::cout << "Nitche: " << std::endl;
        // std::cout << make_hho_nitsche(cl, ET, cell_infos) << std::endl;

        // compute (phi_n_theta, H(-phi_n_1(u))*phi_n_1)_FC / gamma
        K_cont += make_hho_heaviside_contact(cl, ET, uTF, cell_infos);

        // std::cout << "Heaviside: " << std::endl;
        // std::cout << make_hho_heaviside_contact(cl, ET, uTF, cell_infos) << std::endl;

        // compute (phi_n_theta, [phi_n_1(u)]R-)_FC / gamma
        F_cont += make_hho_negative_contact(cl, ET, uTF, cell_infos);

        // auto Fc1 = make_hho_negative_contact(cl, ET, uTF, cell_infos);
        // std::cout << "Negative: " << Fc1.norm() << std::endl;
        // std::cout << Fc1.transpose() << std::endl;

        // friction contribution
        if (m_rp.m_frot_type != NO_FRICTION)
        {
            if (m_rp.m_frot_type == TRESCA)
            {
                // compute (phi_t_theta, [phi_t_1(u)]_s)_FC / gamma
                F_cont += make_hho_threshold_tresca(cl, ET, uTF, cell_infos);

                // auto Ff1 = make_hho_threshold_tresca(cl, ET, uTF, cell_infos);
                // std::cout << "Threshold: " << Ff1.norm() << std::endl;
                // std::cout << Ff1.transpose() << std::endl

                // compute (phi_t_theta, (d_proj_alpha(u)) phi_t_1)_FC / gamma
                K_cont += make_hho_matrix_tresca(cl, ET, uTF, cell_infos);
            }
            else if (m_rp.m_frot_type == COULOMB)
            {
                // compute (phi_t_theta, [phi_t_1(u)]_s)_FC / gamma
                F_cont += make_hho_threshold_coulomb(cl, ET, uTF, cell_infos);

                // auto Ff1 = make_hho_threshold_tresca(cl, ET, uTF, cell_infos);
                // std::cout << "Threshold: " << Ff1.norm() << std::endl;
                // std::cout << Ff1.transpose() << std::endl

                // compute (phi_t_theta, (d_proj_alpha(u)) phi_t_1)_FC / gamma
                K_cont += make_hho_matrix_coulomb(cl, ET, uTF, cell_infos);
            }
        }

        tc.toc();
        time_contact = tc.elapsed();

        //  std::cout << "K_cont: " << K_cont.norm() << std::endl;
        //  std::cout << K_cont << std::endl;
        // std::cout << "F_cont: " << F_cont.norm() << std::endl;
        // std::cout << F_cont.transpose() << std::endl;
    }

    template<typename CellBasis>
    scalar_type
    eval_uT_n(const CellBasis& cb, const vector_type& uTF, const vector_static& n, const point_type& pt) const
    {
        const vector_type uT_n = make_hho_u_n(n, cb, pt);

        return uT_n.dot(uTF.head(cb.size()));
    }

    template<typename CellBasis>
    vector_static
    eval_uT_t(const CellBasis& cb, const vector_type& uTF, const vector_static& n, const point_type& pt) const
    {
        const auto uT_t = make_hho_u_t(n, cb, pt);

        return uT_t.transpose() * (uTF.head(cb.size()));
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
    eval_stress_nn(const vector_type& ET_uTF, const GradBasis& gb, const vector_static& n, const point_type& pt) const
    {
        const auto stress = eval_stress(ET_uTF, gb, pt);
        return (stress * n).dot(n);
    }

    template<typename GradBasis>
    vector_static
    eval_stress_nt(const vector_type& ET_uTF, const GradBasis& gb, const vector_static& n, const point_type& pt) const
    {
        const matrix_static sig  = eval_stress(ET_uTF, gb, pt);
        const vector_static s_n  = sig * n;
        const auto          s_nn = s_n.dot(n);
        return s_n - s_nn * n;
    }

    template<typename GradBasis, typename CellBasis>
    scalar_type
    eval_phi_n_uT(const vector_type&   ET_uTF,
                  const GradBasis&     gb,
                  const CellBasis&     cb,
                  const vector_type&   uTF,
                  const vector_static& n,
                  const scalar_type&   gamma_F,
                  const point_type&    pt) const
    {
        const scalar_type sigma_nn = eval_stress_nn(ET_uTF, gb, n, pt);
        const scalar_type uT_n     = eval_uT_n(cb, uTF, n, pt);
        const scalar_type g0       = make_hho_distance(pt, n);

        return sigma_nn - gamma_F * (uT_n - g0);
    }

    template<typename GradBasis, typename CellBasis>
    scalar_type
    eval_proj_phi_n_uT(const vector_type&   ET_uTF,
                       const GradBasis&     gb,
                       const CellBasis&     cb,
                       const vector_type&   uTF,
                       const vector_static& n,
                       const scalar_type&   gamma_F,
                       const point_type&    pt) const
    {
        const scalar_type phi_n_1_u = eval_phi_n_uT(ET_uTF, gb, cb, uTF, n, gamma_F, pt);

        if (phi_n_1_u <= scalar_type(0))
            return phi_n_1_u;

        return scalar_type(0);
    }

    template<typename GradBasis, typename CellBasis>
    vector_static
    eval_phi_t_uT(const vector_type&   ET_uTF,
                  const GradBasis&     gb,
                  const CellBasis&     cb,
                  const vector_type&   uTF,
                  const vector_static& n,
                  const scalar_type&   gamma_F,
                  const point_type&    pt) const
    {
        const auto sigma_nt = eval_stress_nt(ET_uTF, gb, n, pt);
        const auto uT_t     = eval_uT_t(cb, uTF, n, pt);

        return sigma_nt - gamma_F * uT_t;
    }

    template<typename GradBasis, typename CellBasis>
    vector_static
    eval_proj_phi_t_uT(const vector_type&   ET_uTF,
                       const GradBasis&     gb,
                       const CellBasis&     cb,
                       const vector_type&   uTF,
                       const vector_static& n,
                       const scalar_type&   gamma_F,
                       const scalar_type&   s,
                       const point_type&    pt) const
    {
        const vector_static phi_t_1_u = eval_phi_t_uT(ET_uTF, gb, cb, uTF, n, gamma_F, pt);

        return make_proj_alpha(phi_t_1_u, s);
    }

    template<typename GradBasis, typename FaceBasis>
    scalar_type
    eval_phi_n_uF(const face_type&     fc,
                  const vector_type&   ET_uTF,
                  const GradBasis&     gb,
                  const FaceBasis&     fb,
                  const vector_type&   uTF,
                  const size_t         offset,
                  const vector_static& n,
                  const scalar_type&   gamma_F,
                  const point_type&    pt) const
    {
        const vector_type uF       = uTF.segment(offset, fb.size());
        const scalar_type sigma_nn = eval_stress_nn(ET_uTF, gb, n, pt);
        const scalar_type uF_n     = eval_uF_n(fb, uF, n, pt);
        const scalar_type g0       = make_hho_distance(pt, n);
        const scalar_type gap      = priv::compute_gap_fb(m_msh, fc, fb, uF, pt, n);

        // std::cout << gap << std::endl;

        return sigma_nn + gamma_F * gap;
    }

    template<typename GradBasis, typename FaceBasis>
    scalar_type
    eval_proj_phi_n_uF(const face_type&     fc,
                       const vector_type&   ET_uTF,
                       const GradBasis&     gb,
                       const FaceBasis&     fb,
                       const vector_type&   uTF,
                       const size_t         offset,
                       const vector_static& n,
                       const scalar_type&   gamma_F,
                       const point_type&    pt) const
    {
        const scalar_type phi_n_1_u = eval_phi_n_uF(ET_uTF, gb, fb, uTF, offset, n, gamma_F, pt);

        if (phi_n_1_u <= scalar_type(0))
            return phi_n_1_u;

        return scalar_type(0);
    }

    template<typename GradBasis, typename FaceBasis>
    vector_static
    eval_phi_t_uF(const vector_type&   ET_uTF,
                  const GradBasis&     gb,
                  const FaceBasis&     fb,
                  const vector_type&   uTF,
                  const size_t&        offset,
                  const vector_static& n,
                  const scalar_type&   gamma_F,
                  const point_type&    pt) const
    {
        const vector_type uF       = uTF.segment(offset, fb.size());
        const auto        sigma_nt = eval_stress_nt(ET_uTF, gb, n, pt);
        const auto        uF_t     = eval_uF_t(fb, uF, n, pt);

        return sigma_nt - gamma_F * uF_t;
    }

    template<typename GradBasis, typename FaceBasis>
    vector_static
    eval_proj_tresca_phi_t_uF(const vector_type&   ET_uTF,
                              const GradBasis&     gb,
                              const FaceBasis&     fb,
                              const vector_type&   uTF,
                              const size_t&        offset,
                              const vector_static& n,
                              const scalar_type&   gamma_F,
                              const scalar_type&   s,
                              const point_type&    pt) const
    {
        const vector_static phi_t_1_u = eval_phi_t_uF(ET_uTF, gb, fb, uTF, offset, n, gamma_F, pt);

        return make_proj_alpha(phi_t_1_u, s);
    }

    template<typename GradBasis, typename FaceBasis>
    vector_static
    eval_proj_coulomb_phi_t_uF(const face_type&     fc,
                               const vector_type&   ET_uTF,
                               const GradBasis&     gb,
                               const FaceBasis&     fb,
                               const vector_type&   uTF,
                               const size_t&        offset,
                               const vector_static& n,
                               const scalar_type&   gamma_F,
                               const scalar_type&   Fc,
                               const point_type&    pt) const
    {
        const vector_static phi_t_1_u      = eval_phi_t_uF(ET_uTF, gb, fb, uTF, offset, n, gamma_F, pt);
        const scalar_type   phi_n_1_u      = eval_phi_n_uF(fc, ET_uTF, gb, fb, uTF, offset, n, gamma_F, pt);
        const scalar_type   proj_phi_n_1_u = std::min(scalar_type(0), phi_n_1_u);
        const scalar_type   fric_bound     = -Fc * proj_phi_n_1_u;

        return make_proj_alpha(phi_t_1_u, fric_bound);
    }

    template<typename GradBasis>
    matrix_static
    eval_stress(const vector_type& ET_uTF, const GradBasis& gb, const point_type pt) const
    {
        const auto gphi = gb.eval_functions(pt);

        const matrix_static Eu = disk::eval(ET_uTF, gphi);

        return 2 * m_material_data.getMu() * Eu + m_material_data.getLambda() * Eu.trace() * matrix_static::Identity();
    }
};

} // end namespace MK

} // end namespace disk
