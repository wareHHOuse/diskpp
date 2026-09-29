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

#include "diskpp/bases/bases.hpp"
#include "diskpp/boundary_conditions/boundary_conditions.hpp"
#include "diskpp/common/eigen.hpp"
#include "diskpp/common/timecounter.hpp"
#include "diskpp/mechanics/NewtonSolver/ContactLinearization.hpp"
#include "diskpp/mechanics/NewtonSolver/NonLinearParameters.hpp"
#include "diskpp/mechanics/behaviors/laws/materialData.hpp"
#include "diskpp/mechanics/behaviors/maths_tensor.hpp"
#include "diskpp/methods/hho"
#include "diskpp/quadratures/quadratures.hpp"

#include <cassert>

namespace disk {

namespace mechanics {

template < typename MeshType >
class contact_contribution {
  private:
    typedef MeshType mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell cell_type;
    typedef typename mesh_type::face face_type;
    typedef point< scalar_type, mesh_type::dimension > point_type;
    typedef MaterialData< scalar_type > material_type;
    typedef NonLinearParameters< scalar_type > param_type;
    typedef vector_boundary_conditions< mesh_type > bnd_type;

    const static int dimension = mesh_type::dimension;

    typedef dynamic_matrix< scalar_type > matrix_type;
    typedef dynamic_vector< scalar_type > vector_type;

    typedef static_matrix< scalar_type, dimension, dimension > matrix_static;
    typedef static_vector< scalar_type, dimension > vector_static;

    typedef Matrix< scalar_type, Dynamic, dimension > matrix_d_type;

    const mesh_type &m_msh;
    const material_type &m_material_data;
    const param_type &m_rp;
    const bnd_type &m_bnd;

    scalar_type m_time;
    // dv/du of the time integrator; 1 when the friction law is displacement-based.
    scalar_type m_cN = scalar_type( 1 );
    ContactKinematics m_cont_kine;

    // contact contrib;
    // normal part of u : u_n = u.n
    template < typename TraceBasis >
    vector_type make_hho_u_n( const vector_static &n, const TraceBasis &tb,
                              const point_type &pt ) const {
        const auto t_phi = tb.eval_functions( pt );

        // std::cout << "t_phi: " << t_phi.transpose() << std::endl;

        //(phi_T . n)
        return disk::priv::inner_product( t_phi, n );
    }

    // tangential part of u : u_t = u - u_n*n
    template < typename TraceBasis >
    matrix_d_type make_hho_u_t( const vector_static &n, const TraceBasis &tb,
                                const point_type &pt ) const {
        const auto t_phi = tb.eval_functions( pt );
        const auto u_n = make_hho_u_n( n, tb, pt );

        // phi_T - (phi_T . n)n
        return t_phi - disk::priv::inner_product( u_n, n );
    }

    // cauchy traction : sigma_n = sigma * n
    template < typename GradBasis >
    matrix_d_type make_hho_sigma_n( const matrix_type &ET, const vector_static &n,
                                    const GradBasis &gb, const point_type &pt ) const {
        matrix_d_type sigma_n = matrix_d_type::Zero( ET.cols(), dimension );

        const auto gphi = gb.eval_functions( pt );
        const auto gphi_n = disk::priv::inner_product( gphi, n );

        sigma_n = 2.0 * m_material_data.getMu() * ET.transpose() * gphi_n;

        const auto gphi_trace_n = disk::priv::inner_product( disk::trace( gphi ), n );

        sigma_n += m_material_data.getLambda() * ET.transpose() * gphi_trace_n;

        // sigma. n
        return sigma_n;
    }

    template < typename GradBasis >
    vector_type make_hho_sigma_nn( const matrix_type &ET, const vector_static &n,
                                   const GradBasis &gb, const point_type &pt ) const {
        const auto sigma_n = make_hho_sigma_n( ET, n, gb, pt );

        // sigma_n . n
        return disk::priv::inner_product( sigma_n, n );
    }

    vector_type make_hho_sigma_nn( const matrix_d_type &sigma_n, const vector_static &n ) const {
        // sigma_n . n
        return disk::priv::inner_product( sigma_n, n );
    }

    template < typename GradBasis >
    matrix_d_type make_hho_sigma_nt( const matrix_type &ET, const vector_static &n,
                                     const GradBasis &gb, const point_type &pt ) const {
        const auto sigma_n = make_hho_sigma_n( ET, n, gb, pt );
        const auto sigma_nn = make_hho_sigma_nn( sigma_n, n );

        // sigma_n - sigma_nn * n
        return sigma_n - disk::priv::inner_product( sigma_nn, n );
    }

    // Cell-version contribution

    vector_type
    make_hho_phi_n_uT( const vector_type &sigma_nn,
                       const vector_type &uT_n,
                       scalar_type theta,
                       scalar_type gamma_n_F ) const {
        vector_type phi_n = theta * sigma_nn;
        phi_n.head( uT_n.size() ) -= gamma_n_F * uT_n;

        // theta * sigma_nn - gamma uT_n
        return phi_n;
    }

    matrix_d_type
    make_hho_phi_t_uT( const matrix_d_type &sigma_nt,
                       const matrix_d_type &uT_t,
                       scalar_type theta,
                       scalar_type gamma_t_F ) const {
        matrix_d_type phi_t = theta * sigma_nt;
        phi_t.block( 0, 0, uT_t.rows(), dimension ) -= gamma_t_F * uT_t;

        // theta * sigma_nt - gamma uT_t
        return phi_t;
    }

    vector_type
    make_hho_phi_n_uF( const vector_type &sigma_nn,
                       const vector_type &uF_n,
                       scalar_type theta,
                       scalar_type gamma_n_F,
                       size_t offset ) const {
        vector_type phi_n = theta * sigma_nn;

        assert( offset + uF_n.size() <= phi_n.size() );

        phi_n.segment( offset, uF_n.size() ) -= gamma_n_F * uF_n;

        // theta * sigma_nn - gamma uF_n
        return phi_n;
    }

    matrix_d_type
    make_hho_phi_t_uF( const matrix_d_type &sigma_nt,
                       const matrix_d_type &uF_t,
                       scalar_type theta,
                       scalar_type gamma_t_F,
                       size_t offset ) const {
        matrix_d_type phi_t = theta * sigma_nt;

        phi_t.block( offset, 0, uF_t.rows(), dimension ) -= gamma_t_F * uF_t;

        // theta * sigma_nt - gamma uF_t
        return phi_t;
    }

    vector_type
    make_hho_dphi_n_uT( const vector_type &sigma_nn_derivative,
                        const vector_type &gap_derivative,
                        const scalar_type gamma_n_F ) const {
        vector_type derivative = sigma_nn_derivative;

        assert( gap_derivative.size() <= derivative.size() );

        derivative.head( gap_derivative.size() ) += gamma_n_F * gap_derivative;

        /*
         * D Phi_n =
         *
         * D sigma_nn + gamma_n_F D gap.
         */
        return derivative;
    }

    vector_type
    make_hho_dphi_n_uF( const vector_type &sigma_nn_derivative,
                        const vector_type &gap_derivative,
                        const scalar_type gamma_n_F,
                        const size_t offset ) const {
        vector_type derivative = sigma_nn_derivative;

        assert( offset + gap_derivative.size() <= static_cast< size_t >( derivative.size() ) );

        derivative.segment( offset, gap_derivative.size() ) += gamma_n_F * gap_derivative;

        return derivative;
    }

    // projection on [-inf;0]
    scalar_type
    make_proj_Rmin( const scalar_type x ) const {
        return std::min( scalar_type( 0. ), x );
    }

    // projection on [-inf;0]
    scalar_type
    coulomb_bound( const scalar_type fric, const scalar_type normal_bound ) const {
        return fric * std::abs( normal_bound );
    }

    // projection on the ball of radius alpha centered on 0
    vector_static make_proj_alpha( const vector_static &x, scalar_type alpha ) const {
        const scalar_type x_norm = x.norm();

        if ( x_norm <= alpha ) {
            return x;
        }

        return alpha * x / x_norm;
    }

    // derivative of the projection on the ball of radius alpha centered on 0 versus x
    matrix_static
    make_du_proj_alpha( const vector_static &x, scalar_type alpha ) const {
        const scalar_type x_norm = x.norm();

        if ( alpha <= std::numeric_limits< scalar_type >::epsilon() )
            return matrix_static::Zero();

        if ( x_norm <= alpha ) {
            return matrix_static::Identity();
        }

        return alpha / x_norm *
               ( matrix_static::Identity() - disk::Kronecker( x, x ) / ( x_norm * x_norm ) );
    }

    // derivative of the projection on the ball of radius alpha centered on 0 versus alpha
    vector_static
    make_da_proj_alpha( const vector_static &x, scalar_type alpha ) const {
        const scalar_type x_norm = x.norm();

        if ( alpha <= std::numeric_limits< scalar_type >::epsilon() || x_norm <= alpha )
            return vector_static::Zero();

        return x / x_norm;
    }

    // Compute:
    //
    // theta * [
    // (sigma_nn, sigma_nn)_Fc / gamma_n
    // + (sigma_nt, sigma_nt)_Fc / gamma_t
    // ].
    //
    // In the frictionless case, only the normal contribution is assembled.
    matrix_type
    make_hho_nitsche( const cell_type &cl,
                      const matrix_type &ET,
                      const CellDegreeInfo< MeshType > &cell_infos ) const {
        const auto gb = make_sym_matrix_monomial_basis( m_msh, cl, cell_infos.grad_degree() );

        matrix_type nitsche = matrix_type::Zero( ET.cols(), ET.cols() );

        const auto contact_faces = m_bnd.faces_with_contact( cl );

        for ( const auto &fc : contact_faces ) {
            const vector_static reference_normal = normal( m_msh, cl, fc );

            const auto quadrature_points = integrate( m_msh, fc, 2 * cell_infos.grad_degree() + 2 );

            const scalar_type hF = diameter( m_msh, fc );

            const scalar_type gamma_n_F = m_rp.gamma_0_n() / hF;
            const scalar_type gamma_t_F = m_rp.gamma_0_t() / hF;

            for ( const auto &qp : quadrature_points ) {

                const matrix_d_type sigma_n =
                    make_hho_sigma_n( ET, reference_normal, gb, qp.point() );

                /*
                 * sigma_nn = (sigma(.) N) . N.
                 *
                 */
                const vector_type sigma_nn = make_hho_sigma_nn( sigma_n, reference_normal );

                /*
                 * Normal contribution:
                 *
                 * theta / gamma_n* integral sigma_nn(z) sigma_nn(u).
                 */
                const vector_type weighted_sigma_nn = qp.weight() / gamma_n_F * sigma_nn;

                nitsche += disk::priv::outer_product( weighted_sigma_nn, sigma_nn );

                /*
                 * In the frictionless formulation, there is no
                 * tangential contact constraint to enforce.
                 */
                if ( m_rp.m_frot_type == NO_FRICTION )
                    continue;

                /*
                 * Tangential traction operator:
                 *
                 * sigma_nt = sigma_n- sigma_nn N = P_t(N) sigma(.) N.
                 */
                const matrix_d_type sigma_nt =
                    sigma_n - disk::priv::inner_product( sigma_nn, reference_normal );

                /*
                 * Tangential contribution:
                 *
                 * theta / gamma_t*integral sigma_nt(z) . sigma_nt(u).
                 */
                const matrix_d_type weighted_sigma_nt = qp.weight() / gamma_t_F * sigma_nt;

                nitsche += disk::priv::outer_product( weighted_sigma_nt, sigma_nt );
            }
        }

        return m_rp.m_theta * nitsche;
    }

    // compute (phi_n_theta, H(-phi_n_1(u))*phi_n_1)_FC / gamma
    matrix_type
    make_hho_heaviside_contact( const cell_type &cl,
                                const matrix_type &ET,
                                const vector_type &uTF,
                                const CellDegreeInfo< MeshType > &cell_infos ) const {
        const auto cb = make_vector_monomial_basis( m_msh, cl, cell_infos.cell_degree() );

        const auto gb = make_sym_matrix_monomial_basis( m_msh, cl, cell_infos.grad_degree() );

        matrix_type lhs = matrix_type::Zero( uTF.size(), uTF.size() );

        const auto fcs = faces( m_msh, cl );

        size_t offset = cb.size();

        const auto fcs_di = cell_infos.facesDegreeInfo();

        size_t face_i = 0;

        const vector_type ET_uTF = ET * uTF;

        for ( const auto &fc : fcs ) {
            const auto fdi = fcs_di[face_i++];

            const auto face_degree = fdi.degree();

            const auto fb = make_vector_monomial_basis( m_msh, fc, face_degree );

            const auto face_basis_size = fb.size();

            if ( m_bnd.is_contact_face( fc ) ) {
                const auto contact_type = m_bnd.contact_boundary_type( fc );

                const vector_type uF = uTF.segment( offset, fb.size() );
                const auto ni = normal( m_msh, cl, fc );
                const auto nc = NormalKinematics( m_cont_kine, m_msh, fc, fb, uF, ni );

                const auto qp_degree =
                    std::max( cell_infos.cell_degree(), cell_infos.grad_degree() );

                const auto quadrature_points = integrate( m_msh, fc, 2 * qp_degree + 2 );

                const auto hF = diameter( m_msh, fc );
                const auto gamma_n_F = m_rp.gamma_0_n() / hF;
                const auto gap_function = m_bnd.contact_boundary_gap( fc );

                for ( const auto &qp : quadrature_points ) {
                    /*
                     * D sigma_nn.
                     */
                    const vector_type sigma_nn_derivative =
                        make_hho_sigma_nn( ET, nc.reference_normal, gb, qp.point() );

                    scalar_type phi_n_value = scalar_type( 0 );
                    vector_type phi_n_theta, dphi_n;

                    if ( contact_type == disk::SIGNORINI_CELL ) {
                        const vector_type uT = uTF.head( cb.size() );

                        const auto gap_data = priv::linearize_gap_fb( cb,
                                                                      uT,
                                                                      gap_function,
                                                                      qp.point(),
                                                                      nc.cont_kine,
                                                                      nc.kinematic_normal,
                                                                      m_time );

                        if ( !gap_data.valid )
                            continue;

                        phi_n_value = eval_phi_n_uT( fc,
                                                     ET_uTF,
                                                     gb,
                                                     cb,
                                                     uTF,
                                                     nc.reference_normal,
                                                     nc.kinematic_normal,
                                                     gamma_n_F,
                                                     qp.point() );

                        // always reference normal for test function
                        const vector_type uT_n_test =
                            make_hho_u_n( nc.reference_normal, cb, qp.point() );

                        phi_n_theta = make_hho_phi_n_uT(
                            sigma_nn_derivative, uT_n_test, m_rp.m_theta, gamma_n_F );

                        /*
                         * Vraie dérivée numérique de Phi_n.
                         */
                        dphi_n = make_hho_dphi_n_uT(
                            sigma_nn_derivative, gap_data.derivative, gamma_n_F );
                    } else {

                        const auto gap_data = priv::linearize_gap_fb( fb,
                                                                      uF,
                                                                      gap_function,
                                                                      qp.point(),
                                                                      nc.cont_kine,
                                                                      nc.kinematic_normal,
                                                                      m_time );

                        if ( !gap_data.valid )
                            continue;

                        phi_n_value = eval_phi_n_uF( fc,
                                                     ET_uTF,
                                                     gb,
                                                     fb,
                                                     uTF,
                                                     offset,
                                                     nc.reference_normal,
                                                     nc.kinematic_normal,
                                                     gamma_n_F,
                                                     qp.point() );

                        // always reference normal for test function
                        const vector_type uF_n_test =
                            make_hho_u_n( nc.reference_normal, fb, qp.point() );

                        phi_n_theta = make_hho_phi_n_uF(
                            sigma_nn_derivative, uF_n_test, m_rp.m_theta, gamma_n_F, offset );

                        dphi_n = make_hho_dphi_n_uF(
                            sigma_nn_derivative, gap_data.derivative, gamma_n_F, offset );
                    }

                    /*
                     * D [Phi_n]_- = H(-Phi_n) D Phi_n.
                     */
                    if ( phi_n_value <= scalar_type( 0 ) ) {
                        const auto weighted_phi_n_theta =
                            disk::priv::inner_product( qp.weight() / gamma_n_F, phi_n_theta );

                        lhs += disk::priv::outer_product( weighted_phi_n_theta, dphi_n );
                    }
                }
            }

            offset += face_basis_size;
        }

        return lhs;
    }

    // compute (phi_n_theta, [phi_n_1(u)]R-)_FC / gamma
    vector_type make_hho_negative_contact( const cell_type &cl, const matrix_type &ET,
                                           const vector_type &uTF,
                                           const CellDegreeInfo< MeshType > &cell_infos ) const {
        const auto cb = make_vector_monomial_basis( m_msh, cl, cell_infos.cell_degree() );
        const auto gb = make_sym_matrix_monomial_basis( m_msh, cl, cell_infos.grad_degree() );

        vector_type rhs = vector_type::Zero( uTF.size() );

        const auto fcs = faces( m_msh, cl );
        size_t offset = cb.size();

        const auto fcs_di = cell_infos.facesDegreeInfo();
        size_t face_i = 0;

        const vector_type ET_uTF = ET * uTF;

        for ( auto &fc : fcs ) {
            const auto fdi = fcs_di[face_i++];
            const auto facedeg = fdi.degree();
            const auto fb = make_vector_monomial_basis( m_msh, fc, facedeg );
            const auto fbs = fb.size();

            if ( m_bnd.is_contact_face( fc ) ) {
                const auto contact_type = m_bnd.contact_boundary_type( fc );
                const vector_type uF = uTF.segment( offset, fb.size() );
                const auto ni = normal( m_msh, cl, fc );
                const auto nc = NormalKinematics( m_cont_kine, m_msh, fc, fb, uF, ni );
                const auto qp_deg = std::max( cell_infos.cell_degree(), cell_infos.grad_degree() );
                const auto qps = integrate( m_msh, fc, 2 * qp_deg + 2 );
                const auto hF = diameter( m_msh, fc );
                const auto gamma_n_F = m_rp.gamma_0_n() / hF;

                for ( auto &qp : qps ) {
                    const vector_type sigma_nn =
                        make_hho_sigma_nn( ET, nc.reference_normal, gb, qp.point() );

                    if ( contact_type == disk::SIGNORINI_CELL ) {

                        const scalar_type phi_n_1_u = eval_phi_n_uT( fc,
                                                                     ET_uTF,
                                                                     gb,
                                                                     cb,
                                                                     uTF,
                                                                     nc.reference_normal,
                                                                     nc.kinematic_normal,
                                                                     gamma_n_F,
                                                                     qp.point() );

                        // [phi_n_1_u]_R-
                        if ( phi_n_1_u <= scalar_type( 0 ) ) {
                            // always reference normal for test function
                            const vector_type uT_n =
                                make_hho_u_n( nc.reference_normal, cb, qp.point() );
                            const vector_type phi_n_theta =
                                make_hho_phi_n_uT( sigma_nn, uT_n, m_rp.m_theta, gamma_n_F );

                            rhs += ( qp.weight() / gamma_n_F * phi_n_1_u ) * phi_n_theta;
                        }
                    } else {
                        const scalar_type phi_n_1_u = eval_phi_n_uF( fc,
                                                                     ET_uTF,
                                                                     gb,
                                                                     fb,
                                                                     uTF,
                                                                     offset,
                                                                     nc.reference_normal,
                                                                     nc.kinematic_normal,
                                                                     gamma_n_F,
                                                                     qp.point() );

                        // std::cout << "qp: " << qp.point() << std::endl;
                        // std::cout << "phi_n_1_u: " << phi_n_1_u << std::endl;

                        // [phi_n_1_u]_R-
                        if ( phi_n_1_u <= scalar_type( 0 ) ) {
                            // always reference normal for test function
                            const vector_type uF_n =
                                make_hho_u_n( nc.reference_normal, fb, qp.point() );
                            const vector_type phi_n_theta = make_hho_phi_n_uF(
                                sigma_nn, uF_n, m_rp.m_theta, gamma_n_F, offset );

                            // std::cout << "sigma_nn: " << sigma_nn.transpose() << std::endl;
                            // std::cout << "uF_n: " << uF_n.transpose() << std::endl;
                            // std::cout << "phi_n_theta: " << phi_n_theta.transpose() << std::endl;

                            rhs += ( qp.weight() / gamma_n_F * phi_n_1_u ) * phi_n_theta;
                        }
                    }
                }
            }
            offset += fbs;
        }
        return rhs;
    }

    // compute (phi_t_theta, [phi_t_1(u,w)]_(s))_FC / gamma
    vector_type
    make_hho_threshold_tresca( const cell_type &cl,
                               const matrix_type &ET,
                               const vector_type &uTF,
                               const vector_type &vuTF,
                               const CellDegreeInfo< MeshType > &cell_infos ) const {
        const auto cb = make_vector_monomial_basis( m_msh, cl, cell_infos.cell_degree() );
        const auto gb = make_sym_matrix_monomial_basis( m_msh, cl, cell_infos.grad_degree() );

        vector_type rhs = vector_type::Zero( uTF.size() );

        const auto fcs = faces( m_msh, cl );
        size_t offset = cb.size();

        const auto fcs_di = cell_infos.facesDegreeInfo();
        size_t face_i = 0;

        const vector_type ET_uTF = ET * uTF;

        for ( auto &fc : fcs ) {
            const auto fdi = fcs_di[face_i++];
            const auto facedeg = fdi.degree();
            const auto fb = make_vector_monomial_basis( m_msh, fc, facedeg );
            const auto fbs = fb.size();

            if ( m_bnd.is_contact_face( fc ) ) {
                const auto contact_type = m_bnd.contact_boundary_type( fc );

                const vector_type uF = uTF.segment( offset, fb.size() );
                const auto ni = normal( m_msh, cl, fc );
                const auto nc = NormalKinematics( m_cont_kine, m_msh, fc, fb, uF, ni );

                const auto qp_deg = std::max( cell_infos.cell_degree(), cell_infos.grad_degree() );
                const auto qps = integrate( m_msh, fc, 2 * qp_deg + 2 );
                const auto hF = diameter( m_msh, fc );
                const auto gamma_t_F = m_rp.gamma_0_t() / hF;

                const auto s_func = m_bnd.contact_boundary_func( fc );

                for ( auto &qp : qps ) {
                    const auto sigma_nt =
                        make_hho_sigma_nt( ET, nc.reference_normal, gb, qp.point() );

                    if ( contact_type == disk::SIGNORINI_CELL ) {
                        // always reference normal for test function
                        const auto uT_t_test = make_hho_u_t( nc.reference_normal, cb, qp.point() );

                        const auto phi_t_theta =
                            make_hho_phi_t_uT( sigma_nt, uT_t_test, m_rp.m_theta, gamma_t_F );

                        const vector_static phi_t_1_u_proj =
                            eval_proj_phi_t_uT( ET_uTF,
                                                gb,
                                                cb,
                                                vuTF,
                                                nc.reference_normal,
                                                nc.kinematic_normal,
                                                gamma_t_F,
                                                s_func( qp.point() ),
                                                qp.point() );

                        const vector_static qp_phi_t_1_u_pro =
                            qp.weight() * phi_t_1_u_proj / gamma_t_F;

                        rhs += disk::priv::inner_product( phi_t_theta, qp_phi_t_1_u_pro );
                    } else {
                        const auto uF_t_test = make_hho_u_t( nc.reference_normal, fb, qp.point() );

                        const auto phi_t_theta = make_hho_phi_t_uF(
                            sigma_nt, uF_t_test, m_rp.m_theta, gamma_t_F, offset );

                        // std::cout << "sigma_nt: " << sigma_nt.transpose() << std::endl;
                        //                         std::cout << "uF_t: " << uF_t.transpose() <<
                        //                         std::endl; std::cout << "phi_t_theta: " <<
                        //                         phi_t_theta.transpose() << std::endl;

                        const vector_static phi_t_1_u_proj =
                            eval_proj_tresca_phi_t_uF( ET_uTF,
                                                       gb,
                                                       fb,
                                                       vuTF,
                                                       offset,
                                                       nc.reference_normal,
                                                       nc.kinematic_normal,
                                                       gamma_t_F,
                                                       s_func( qp.point() ),
                                                       qp.point() );

                        const vector_static qp_phi_t_1_u_pro =
                            qp.weight() * phi_t_1_u_proj / gamma_t_F;

                        // std::cout << "phi_t_1_u_proj: " << phi_t_1_u_proj.transpose() <<
                        // std::endl;

                        rhs += disk::priv::inner_product( phi_t_theta, qp_phi_t_1_u_pro );
                    }
                }
            }
            offset += fbs;
        }
        return rhs;
    }

    // compute (phi_t_theta, (d_proj_alpha(u,w)) phi_t_1)_FC / gamma
    matrix_type
    make_hho_matrix_tresca( const cell_type &cl,
                            const matrix_type &ET,
                            const vector_type &uTF,
                            const vector_type &vuTF,
                            const CellDegreeInfo< MeshType > &cell_infos ) const {
        const auto cb = make_vector_monomial_basis( m_msh, cl, cell_infos.cell_degree() );
        const auto gb = make_sym_matrix_monomial_basis( m_msh, cl, cell_infos.grad_degree() );

        matrix_type lhs = matrix_type::Zero( uTF.size(), uTF.size() );

        const auto fcs = faces( m_msh, cl );
        size_t offset = cb.size();

        const auto fcs_di = cell_infos.facesDegreeInfo();
        size_t face_i = 0;

        const vector_type ET_uTF = ET * uTF;

        for ( auto &fc : fcs ) {
            const auto fdi = fcs_di[face_i++];
            const auto facedeg = fdi.degree();
            const auto fb = make_vector_monomial_basis( m_msh, fc, facedeg );
            const auto fbs = fb.size();

            if ( m_bnd.is_contact_face( fc ) ) {
                const auto contact_type = m_bnd.contact_boundary_type( fc );
                const vector_type uF = uTF.segment( offset, fb.size() );
                const auto ni = normal( m_msh, cl, fc );
                const auto nc = NormalKinematics( m_cont_kine, m_msh, fc, fb, uF, ni );
                const auto qp_deg = std::max( cell_infos.cell_degree(), cell_infos.grad_degree() );
                const auto qps = integrate( m_msh, fc, 2 * qp_deg + 2 );
                const auto hF = diameter( m_msh, fc );
                const auto gamma_t_F = m_rp.gamma_0_t() / hF;

                const auto s_func = m_bnd.contact_boundary_func( fc );

                for ( auto &qp : qps ) {
                    const auto sigma_nt =
                        make_hho_sigma_nt( ET, nc.reference_normal, gb, qp.point() );

                    if ( contact_type == disk::SIGNORINI_CELL ) {
                        // always reference normal for test function
                        const auto uT_t = make_hho_u_t( nc.kinematic_normal, cb, qp.point() );
                        const auto uT_t_test = make_hho_u_t( nc.reference_normal, cb, qp.point() );

                        const auto phi_t_1 =
                            make_hho_phi_t_uT( sigma_nt, m_cN * uT_t, scalar_type( 1 ), gamma_t_F );
                        const auto phi_t_theta =
                            make_hho_phi_t_uT( sigma_nt, uT_t_test, m_rp.m_theta, gamma_t_F );

                        const auto phi_t_1_u = eval_phi_t_uT( ET_uTF,
                                                              gb,
                                                              cb,
                                                              vuTF,
                                                              nc.reference_normal,
                                                              nc.kinematic_normal,
                                                              gamma_t_F,
                                                              qp.point() );
                        const auto d_proj_phi_t_u =
                            make_du_proj_alpha( phi_t_1_u, s_func( qp.point() ) );

                        const auto d_proj_u_phi_t_1 =
                            disk::priv::inner_product( d_proj_phi_t_u, phi_t_1 );

                        const auto qp_phi_t_theta =
                            disk::priv::inner_product( qp.weight() / gamma_t_F, phi_t_theta );

                        lhs += disk::priv::outer_product( qp_phi_t_theta, d_proj_u_phi_t_1 );
                    } else {
                        // always reference normal for test function
                        const auto uF_t = make_hho_u_t( nc.kinematic_normal, fb, qp.point() );
                        const auto uF_t_test = make_hho_u_t( nc.reference_normal, fb, qp.point() );

                        const auto phi_t_1 = make_hho_phi_t_uF(
                            sigma_nt, m_cN * uF_t, scalar_type( 1 ), gamma_t_F, offset );
                        const auto phi_t_theta = make_hho_phi_t_uF(
                            sigma_nt, uF_t_test, m_rp.m_theta, gamma_t_F, offset );

                        const auto phi_t_1_u = eval_phi_t_uF( ET_uTF,
                                                              gb,
                                                              fb,
                                                              vuTF,
                                                              offset,
                                                              nc.reference_normal,
                                                              nc.kinematic_normal,
                                                              gamma_t_F,
                                                              qp.point() );
                        const auto d_proj_phi_t_u =
                            make_du_proj_alpha( phi_t_1_u, s_func( qp.point() ) );

                        const auto d_proj_u_phi_t_1 =
                            disk::priv::inner_product( d_proj_phi_t_u, phi_t_1 );

                        const auto qp_phi_t_theta =
                            disk::priv::inner_product( qp.weight() / gamma_t_F, phi_t_theta );

                        lhs += disk::priv::outer_product( qp_phi_t_theta, d_proj_u_phi_t_1 );
                    }
                }
            }
            offset += fbs;
        }

        return lhs;
    }

    // compute (phi_t_theta, [phi_t_1(u,w)]_(s))_FC / gamma
    vector_type
    make_hho_threshold_coulomb( const cell_type &cl,
                                const matrix_type &ET,
                                const vector_type &uTF,
                                const vector_type &vuTF,
                                const CellDegreeInfo< MeshType > &cell_infos ) const {
        const auto cb = make_vector_monomial_basis( m_msh, cl, cell_infos.cell_degree() );
        const auto gb = make_sym_matrix_monomial_basis( m_msh, cl, cell_infos.grad_degree() );

        vector_type rhs = vector_type::Zero( uTF.size() );

        const auto fcs = faces( m_msh, cl );
        size_t offset = cb.size();

        const auto fcs_di = cell_infos.facesDegreeInfo();
        size_t face_i = 0;

        const vector_type ET_uTF = ET * uTF;

        for ( auto &fc : fcs ) {
            const auto fdi = fcs_di[face_i++];
            const auto facedeg = fdi.degree();
            const auto fb = make_vector_monomial_basis( m_msh, fc, facedeg );
            const auto fbs = fb.size();

            if ( m_bnd.is_contact_face( fc ) ) {
                const auto contact_type = m_bnd.contact_boundary_type( fc );
                const vector_type uF = uTF.segment( offset, fb.size() );
                const auto ni = normal( m_msh, cl, fc );
                const auto nc = NormalKinematics( m_cont_kine, m_msh, fc, fb, uF, ni );
                const auto qp_deg = std::max( cell_infos.cell_degree(), cell_infos.grad_degree() );
                const auto qps = integrate( m_msh, fc, 2 * qp_deg + 2 );
                const auto hF = diameter( m_msh, fc );
                const auto gamma_t_F = m_rp.gamma_0_t() / hF;
                const auto gamma_n_F = m_rp.gamma_0_n() / hF;

                const auto s_func = m_bnd.contact_boundary_func( fc );

                for ( auto &qp : qps ) {
                    const auto sigma_nt =
                        make_hho_sigma_nt( ET, nc.reference_normal, gb, qp.point() );

                    if ( contact_type == disk::SIGNORINI_CELL ) {
                        // always reference normal for test function
                        const auto uT_t_test = make_hho_u_t( nc.reference_normal, cb, qp.point() );

                        const auto phi_t_theta =
                            make_hho_phi_t_uT( sigma_nt, uT_t_test, m_rp.m_theta, gamma_t_F );

                        const vector_static phi_t_1_u_proj =
                            eval_proj_coulomb_phi_t_uT( fc,
                                                        ET_uTF,
                                                        gb,
                                                        cb,
                                                        uTF,
                                                        vuTF,
                                                        nc.reference_normal,
                                                        nc.kinematic_normal,
                                                        gamma_t_F,
                                                        gamma_n_F,
                                                        s_func( qp.point() ),
                                                        qp.point() );

                        const vector_static qp_phi_t_1_u_pro =
                            qp.weight() * phi_t_1_u_proj / gamma_t_F;

                        rhs += disk::priv::inner_product( phi_t_theta, qp_phi_t_1_u_pro );
                    } else {
                        // always reference normal for test function
                        const auto uF_t_test = make_hho_u_t( nc.reference_normal, fb, qp.point() );

                        const auto phi_t_theta = make_hho_phi_t_uF(
                            sigma_nt, uF_t_test, m_rp.m_theta, gamma_t_F, offset );

                        const vector_static phi_t_1_u_proj =
                            eval_proj_coulomb_phi_t_uF( fc,
                                                        ET_uTF,
                                                        gb,
                                                        fb,
                                                        uTF,
                                                        vuTF,
                                                        offset,
                                                        nc.reference_normal,
                                                        nc.kinematic_normal,
                                                        gamma_t_F,
                                                        gamma_n_F,
                                                        s_func( qp.point() ),
                                                        qp.point() );

                        const vector_static qp_phi_t_1_u_pro =
                            qp.weight() * phi_t_1_u_proj / gamma_t_F;

                        rhs += disk::priv::inner_product( phi_t_theta, qp_phi_t_1_u_pro );
                    }
                }
            }
            offset += fbs;
        }
        return rhs;
    }

    // compute (phi_t_theta, (d_proj_alpha(u,w)) phi_t_1)_FC / gamma
    matrix_type
    make_hho_matrix_coulomb( const cell_type &cl,
                             const matrix_type &ET,
                             const vector_type &uTF,
                             const vector_type &vuTF,
                             const CellDegreeInfo< MeshType > &cell_infos ) const {
        const auto cb = make_vector_monomial_basis( m_msh, cl, cell_infos.cell_degree() );
        const auto gb = make_sym_matrix_monomial_basis( m_msh, cl, cell_infos.grad_degree() );

        matrix_type lhs = matrix_type::Zero( uTF.size(), uTF.size() );

        const auto fcs = faces( m_msh, cl );
        size_t offset = cb.size();

        const auto fcs_di = cell_infos.facesDegreeInfo();
        size_t face_i = 0;

        const vector_type ET_uTF = ET * uTF;

        for ( auto &fc : fcs ) {
            const auto fdi = fcs_di[face_i++];
            const auto facedeg = fdi.degree();
            const auto fb = make_vector_monomial_basis( m_msh, fc, facedeg );
            const auto fbs = fb.size();

            if ( m_bnd.is_contact_face( fc ) ) {
                const auto contact_type = m_bnd.contact_boundary_type( fc );
                const vector_type uF = uTF.segment( offset, fbs );
                const auto ni = normal( m_msh, cl, fc );
                const auto nc = NormalKinematics( m_cont_kine, m_msh, fc, fb, uF, ni );
                const auto qp_deg = std::max( cell_infos.cell_degree(), cell_infos.grad_degree() );
                const auto qps = integrate( m_msh, fc, 2 * qp_deg + 2 );
                const auto hF = diameter( m_msh, fc );
                const auto gamma_t_F = m_rp.gamma_0_t() / hF;
                const auto gamma_n_F = m_rp.gamma_0_n() / hF;

                const auto s_func = m_bnd.contact_boundary_func( fc );
                const auto gap_function = m_bnd.contact_boundary_gap( fc );

                for ( auto &qp : qps ) {
                    const auto sigma_nt =
                        make_hho_sigma_nt( ET, nc.reference_normal, gb, qp.point() );

                    if ( contact_type == disk::SIGNORINI_CELL ) {
                        // always reference normal for test function
                        const auto uT_t = make_hho_u_t( nc.kinematic_normal, cb, qp.point() );
                        const auto uT_t_test = make_hho_u_t( nc.reference_normal, cb, qp.point() );

                        const auto phi_t_1 =
                            make_hho_phi_t_uT( sigma_nt, m_cN * uT_t, scalar_type( 1 ), gamma_t_F );
                        const auto phi_t_theta =
                            make_hho_phi_t_uT( sigma_nt, uT_t_test, m_rp.m_theta, gamma_t_F );

                        const auto phi_t_1_u = eval_phi_t_uT( ET_uTF,
                                                              gb,
                                                              cb,
                                                              vuTF,
                                                              nc.reference_normal,
                                                              nc.kinematic_normal,
                                                              gamma_t_F,
                                                              qp.point() );
                        const scalar_type phi_n_1_u = eval_phi_n_uT( fc,
                                                                     ET_uTF,
                                                                     gb,
                                                                     cb,
                                                                     uTF,
                                                                     nc.reference_normal,
                                                                     nc.kinematic_normal,
                                                                     gamma_n_F,
                                                                     qp.point() );

                        const scalar_type proj_phi_n_1_u = make_proj_Rmin( phi_n_1_u );
                        const scalar_type fric_bound =
                            coulomb_bound( s_func( qp.point() ), proj_phi_n_1_u );
                        const auto d_proj_phi_t_u = make_du_proj_alpha( phi_t_1_u, fric_bound );

                        const auto d_proj_u_phi_t_1 =
                            disk::priv::inner_product( d_proj_phi_t_u, phi_t_1 );

                        const auto qp_phi_t_theta =
                            disk::priv::inner_product( qp.weight() / gamma_t_F, phi_t_theta );

                        // compute (phi_t_theta, (d_proj_alpha(u,w)) phi_t_1)_FC / gamma
                        lhs += disk::priv::outer_product( qp_phi_t_theta, d_proj_u_phi_t_1 );

                        // d/du of the friction radius s(u) = -F [phi_n_1(u)]_-, slip regime
                        if ( m_rp.consistentFrictionTangent() && phi_n_1_u < scalar_type( 0 ) &&
                             phi_t_1_u.norm() > fric_bound ) {

                            const vector_static q_hat = phi_t_1_u / phi_t_1_u.norm();

                            const vector_type uT = uTF.head( cb.size() );

                            const auto gap_data = priv::linearize_gap_fb( cb,
                                                                          uT,
                                                                          gap_function,
                                                                          qp.point(),
                                                                          nc.cont_kine,
                                                                          nc.kinematic_normal,
                                                                          m_time );

                            if ( gap_data.valid ) {
                                const vector_type sigma_nn_derivative =
                                    make_hho_sigma_nn( ET, nc.reference_normal, gb, qp.point() );

                                const vector_type dphi_n = make_hho_dphi_n_uT(
                                    sigma_nn_derivative, gap_data.derivative, gamma_n_F );

                                const vector_type phi_t_theta_qhat = phi_t_theta * q_hat;

                                lhs -= ( qp.weight() / gamma_t_F * s_func( qp.point() ) ) *
                                       disk::priv::outer_product( phi_t_theta_qhat, dphi_n );
                            }
                        }
                    } else {
                        const scalar_type phi_n_1_u = eval_phi_n_uF( fc,
                                                                     ET_uTF,
                                                                     gb,
                                                                     fb,
                                                                     uTF,
                                                                     offset,
                                                                     nc.reference_normal,
                                                                     nc.kinematic_normal,
                                                                     gamma_n_F,
                                                                     qp.point() );

                        /*
                         * Inactive normal contact:
                         *
                         * alpha = -mu [Phi_n]_- = 0.
                         *
                         * The Coulomb contribution and its generalized derivative
                         * are taken equal to zero.
                         */
                        if ( phi_n_1_u > scalar_type( 0 ) ) {
                            continue;
                        }

                        const auto uF_t = make_hho_u_t( nc.kinematic_normal, fb, qp.point() );

                        const auto uF_t_test = make_hho_u_t( nc.reference_normal, fb, qp.point() );

                        const auto phi_t_1 = make_hho_phi_t_uF(
                            sigma_nt, m_cN * uF_t, scalar_type( 1 ), gamma_t_F, offset );

                        const auto phi_t_theta = make_hho_phi_t_uF(
                            sigma_nt, uF_t_test, m_rp.m_theta, gamma_t_F, offset );

                        /*
                         * Phi_t(u,w) = sigma_nt(u) - gamma_t P_t(N) w.
                         */
                        const vector_static phi_t_1_u = eval_phi_t_uF( ET_uTF,
                                                                       gb,
                                                                       fb,
                                                                       vuTF,
                                                                       offset,
                                                                       nc.reference_normal,
                                                                       nc.kinematic_normal,
                                                                       gamma_t_F,
                                                                       qp.point() );

                        /*
                         * Coulomb radius:
                         *
                         * alpha(u) = -mu [Phi_n(u)]_-.
                         */
                        const scalar_type proj_phi_n_1_u = make_proj_Rmin( phi_n_1_u );

                        const scalar_type friction_coefficient = s_func( qp.point() );

                        const scalar_type fric_bound =
                            coulomb_bound( friction_coefficient, proj_phi_n_1_u );

                        /*
                         * First contribution:
                         *
                         * D_{Phi_t} Proj_alpha(Phi_t) D Phi_t.
                         */
                        const matrix_static du_proj_phi_t_u =
                            make_du_proj_alpha( phi_t_1_u, fric_bound );

                        const auto du_proj_u_phi_t_1 =
                            disk::priv::inner_product( du_proj_phi_t_u, phi_t_1 );

                        const auto qp_phi_t_theta =
                            disk::priv::inner_product( qp.weight() / gamma_t_F, phi_t_theta );

                        lhs += disk::priv::outer_product( qp_phi_t_theta, du_proj_u_phi_t_1 );

                        /*
                         * Second contribution:
                         *
                         * D_alpha Proj_alpha(Phi_t) D alpha,
                         *
                         * with:
                         *
                         * D alpha = -mu D Phi_n.
                         *
                         * This term is nonzero only in active contact and
                         * in the sliding regime.
                         */
                        const scalar_type phi_t_norm = phi_t_1_u.norm();

                        const bool active_contact = phi_n_1_u < 0.0;

                        const bool sliding = phi_t_norm > fric_bound;

                        if ( m_rp.consistentFrictionTangent() && active_contact && sliding &&
                             phi_t_norm > 0.0 ) {
                            /*
                             * In sliding:
                             *
                             * D_alpha Proj_alpha(Phi_t) = Phi_t / ||Phi_t||.
                             *
                             * Use the dedicated utility instead of rebuilding q_hat
                             * independently.
                             */
                            const vector_static da_proj_phi_t_u =
                                make_da_proj_alpha( phi_t_1_u, fric_bound );

                            /*
                             * Frozen-normal finite-difference derivative of the gap.
                             */
                            const auto gap_data = priv::linearize_gap_fb( fb,
                                                                          uF,
                                                                          gap_function,
                                                                          qp.point(),
                                                                          nc.cont_kine,
                                                                          nc.kinematic_normal,
                                                                          m_time );

                            if ( gap_data.valid ) {
                                /*
                                 * D sigma_nn, with the reference normal held fixed.
                                 */
                                const vector_type sigma_nn_derivative =
                                    make_hho_sigma_nn( ET, nc.reference_normal, gb, qp.point() );

                                /*
                                 * D Phi_n = D sigma_nn + gamma_n D gap.
                                 */
                                const vector_type dphi_n = make_hho_dphi_n_uF(
                                    sigma_nn_derivative, gap_data.derivative, gamma_n_F, offset );

                                /*
                                 * phi_t_theta * D_alpha Proj_alpha.
                                 */
                                const vector_type phi_t_theta_da_proj =
                                    phi_t_theta * da_proj_phi_t_u;

                                /*
                                 * Since:
                                 *
                                 * D alpha = -mu D Phi_n,
                                 *
                                 * the contribution is:
                                 *
                                 * - weight * mu / gamma_t (Phi_{t,theta} . D_alpha Proj) tensor
                                 * D Phi_n.
                                 */
                                lhs -= ( qp.weight() / gamma_t_F * friction_coefficient ) *
                                       disk::priv::outer_product( phi_t_theta_da_proj, dphi_n );
                            }
                        }
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
    double time_contact;

    contact_contribution( const mesh_type &msh,
                          const material_type &material_data,
                          const param_type &rp,
                          const bnd_type &bnd )
        : m_msh( msh ),
          m_material_data( material_data ),
          m_rp( rp ),
          m_bnd( bnd ),
          m_time( scalar_type( 0 ) ),
          m_cont_kine( rp.getContactKinematics() ) {}

    // One quadrature point of a contact face, with every quantity projected on that
    // face's own discrete normal n = normal( m_msh, cl, fc ).
    struct trace_point {
        point_type pt;           // reference position
        vector_static n;         // discrete facet normal (outward)
        vector_static u;         // displacement trace
        scalar_type gap0;        // initial gap g0( pt, n )
        scalar_type u_n;         // u . n
        scalar_type u_t;         // u . t, t = tangent (2D only)
        scalar_type sigma_nn;    // (sigma n) . n
        scalar_type sigma_nt;    // (sigma n) . t (2D only)
        scalar_type phi_n;       // Nitsche normal indicator, < 0 <=> active contact
        scalar_type phi_t;       // Nitsche tangential trial stress projected on t
        scalar_type fric_bound;  // Fc * [-phi_n]_+
        scalar_type Fc;          // Coulomb coefficient at pt
        scalar_type weight;      // quadrature weight
    };

    // 2D tangent associated with the outward normal n.
    static vector_static tangent_of( const vector_static &n ) {
        vector_static t = vector_static::Zero();
        static_assert( dimension == 2, "contact trace output is 2D only" );
        t( 0 ) = -n( 1 );
        t( 1 ) = n( 0 );
        return t;
    }

    // Per-quadrature-point trace of the contact boundary of one cell. Mirrors the
    // quadrature, gamma_F and friction bound used by nitsche_friction_energy so the
    // reported state matches the state the solver actually enforced.
    std::vector< trace_point >
    contact_boundary_trace( const cell_type &cl, const CellDegreeInfo< MeshType > &cell_infos,
                            const matrix_type &ET, const vector_type &uTF ) const {
        std::vector< trace_point > out;

        const auto cb = make_vector_monomial_basis( m_msh, cl, cell_infos.cell_degree() );
        const auto gb = make_sym_matrix_monomial_basis( m_msh, cl, cell_infos.grad_degree() );

        const auto fcs = faces( m_msh, cl );
        size_t offset = cb.size();
        const auto fcs_di = cell_infos.facesDegreeInfo();
        size_t face_i = 0;

        const vector_type ET_uTF = ET * uTF;

        for ( auto &fc : fcs ) {
            const auto fdi = fcs_di[face_i++];
            const auto fb = make_vector_monomial_basis( m_msh, fc, fdi.degree() );
            const auto fbs = fb.size();

            if ( m_bnd.is_contact_face( fc ) ) {
                const auto contact_type = m_bnd.contact_boundary_type( fc );
                const auto n = normal( m_msh, cl, fc );
                const auto t = tangent_of( n );
                const auto qp_deg = std::max( cell_infos.cell_degree(), cell_infos.grad_degree() );
                const auto qps = integrate( m_msh, fc, 2 * qp_deg + 2 );
                const auto hF = diameter( m_msh, fc );
                const auto gamma_n_F = m_rp.gamma_0_n() / hF;
                const auto gamma_t_F = m_rp.gamma_0_t() / hF;
                const auto gap_func = m_bnd.contact_boundary_gap( fc );

                const vector_type uF = uTF.segment( offset, fbs );

                for ( auto &qp : qps ) {
                    trace_point tp;
                    tp.pt = qp.point();
                    tp.n = n;
                    tp.weight = qp.weight();
                    tp.gap0 = gap_func( qp.point(), n );

                    // displacement trace: cell trace for SIGNORINI_CELL, face trace otherwise
                    if ( contact_type == disk::SIGNORINI_CELL ) {
                        const auto t_phi = cb.eval_functions( qp.point() );
                        tp.u = t_phi.transpose() * uTF.head( cb.size() );
                    } else {
                        const auto t_phi = fb.eval_functions( qp.point() );
                        tp.u = t_phi.transpose() * uF;
                    }

                    tp.u_n = tp.u.dot( n );
                    tp.u_t = tp.u.dot( t );

                    const vector_static sig_n_vec = eval_stress( ET_uTF, gb, qp.point() ) * n;
                    tp.sigma_nn = sig_n_vec.dot( n );
                    tp.sigma_nt = sig_n_vec.dot( t );

                    // Nitsche indicators, each on its own penalty
                    tp.phi_n = tp.sigma_nn + gamma_n_F * ( tp.gap0 - tp.u_n );

                    const vector_static u_tang = tp.u - tp.u_n * n;
                    const vector_static sigma_tang = sig_n_vec - tp.sigma_nn * n;
                    const vector_static phi_t_vec = sigma_tang - gamma_t_F * u_tang;
                    tp.phi_t = phi_t_vec.dot( t );

                    tp.Fc = m_bnd.contact_boundary_func( fc )( qp.point() );
                    tp.fric_bound = -tp.Fc * make_proj_Rmin( tp.phi_n );

                    out.push_back( tp );
                }
            }
            offset += fbs;
        }
        return out;
    }

    // wTF carries the tangential kinematics, uTF always the stress
    void
    compute( const cell_type &cl,
             const CellDegreeInfo< mesh_type > &cell_infos,
             const matrix_type &ET,
             const vector_type &uTF,
             const vector_type &vTF,
             bool tangent_matix,
             const scalar_type c_N = scalar_type( 1 ) ) {
        timecounter tc;
        tc.tic();

        m_cN = c_N;

        // contact contribution
        time_contact = 0.0;

        const auto cb = disk::make_vector_monomial_basis( m_msh, cl, cell_infos.cell_degree() );
        const auto gb = disk::make_sym_matrix_monomial_basis( m_msh, cl, cell_infos.grad_degree() );

        F_cont = vector_type::Zero( uTF.size() );
        K_cont = matrix_type::Zero( uTF.size(), uTF.size() );

        // compute theta/gamma *(sigma_n, sigma_n)_Fc
        if ( m_rp.m_theta != 0.0 ) {
            const matrix_type K_sig = make_hho_nitsche( cl, ET, cell_infos );
            if ( tangent_matix ) {
                K_cont -= K_sig;
            }
            F_cont -= K_sig * uTF;
        }

        // std::cout << "Nitche: " << std::endl;
        // std::cout << make_hho_nitsche(cl, ET, cell_infos) << std::endl;

        // compute (phi_n_theta, H(-phi_n_1(u))*phi_n_1)_FC / gamma
        if ( tangent_matix ) {
            K_cont += make_hho_heaviside_contact( cl, ET, uTF, cell_infos );
        }

        // std::cout << "Heaviside: " << std::endl;
        // std::cout << make_hho_heaviside_contact(cl, ET, uTF, cell_infos) << std::endl;

        // compute (phi_n_theta, [phi_n_1(u)]R-)_FC / gamma
        F_cont += make_hho_negative_contact( cl, ET, uTF, cell_infos );

        // auto Fc1 = make_hho_negative_contact(cl, ET, uTF, cell_infos);
        // std::cout << "Negative: " << Fc1.norm() << std::endl;
        // std::cout << Fc1.transpose() << std::endl;

        // friction contribution
        if ( m_rp.m_frot_type != NO_FRICTION ) {
            if ( m_rp.m_frot_type == TRESCA ) {
                if ( m_rp.isUnsteady() ) {
                    // compute (phi_t_theta, [phi_t_1(v)]_s)_FC / gamma
                    F_cont += make_hho_threshold_tresca( cl, ET, uTF, vTF, cell_infos );

                    // compute (phi_t_theta, (d_proj_alpha(v)) phi_t_1)_FC / gamma
                    if ( tangent_matix ) {
                        K_cont += make_hho_matrix_tresca( cl, ET, uTF, vTF, cell_infos );
                    }
                } else {
                    // compute (phi_t_theta, [phi_t_1(u)]_s)_FC / gamma
                    F_cont += make_hho_threshold_tresca( cl, ET, uTF, uTF, cell_infos );

                    // auto Ff1 = make_hho_threshold_tresca(cl, ET, uTF, cell_infos);
                    // std::cout << "Threshold: " << Ff1.norm() << std::endl;
                    // std::cout << Ff1.transpose() << std::endl

                    // compute (phi_t_theta, (d_proj_alpha(u)) phi_t_1)_FC / gamma
                    if ( tangent_matix ) {
                        K_cont += make_hho_matrix_tresca( cl, ET, uTF, uTF, cell_infos );
                    }
                }
            } else if ( m_rp.m_frot_type == COULOMB ) {
                if ( m_rp.isUnsteady() ) {
                    // compute (phi_t_theta, [phi_t_1(v)]_s)_FC / gamma
                    F_cont += make_hho_threshold_coulomb( cl, ET, uTF, vTF, cell_infos );

                    // compute (phi_t_theta, (d_proj_alpha(v)) phi_t_1)_FC / gamma
                    if ( tangent_matix ) {
                        K_cont += make_hho_matrix_coulomb( cl, ET, uTF, vTF, cell_infos );
                    }
                } else {
                    // compute (phi_t_theta, [phi_t_1(u)]_s)_FC / gamma
                    F_cont += make_hho_threshold_coulomb( cl, ET, uTF, uTF, cell_infos );

                    // auto Ff1 = make_hho_threshold_tresca(cl, ET, uTF, cell_infos);
                    // std::cout << "Threshold: " << Ff1.norm() << std::endl;
                    // std::cout << Ff1.transpose() << std::endl

                    // compute (phi_t_theta, (d_proj_alpha(u)) phi_t_1)_FC / gamma
                    if ( tangent_matix ) {
                        K_cont += make_hho_matrix_coulomb( cl, ET, uTF, uTF, cell_infos );
                    }
                }
            }
        }

        tc.toc();
        time_contact = tc.elapsed();

        //  std::cout << "K_cont: " << K_cont.norm() << std::endl;
        //  std::cout << K_cont << std::endl;
        // std::cout << "F_cont: " << F_cont.norm() << std::endl;
        // std::cout << F_cont.transpose() << std::endl;
    }

    // Nitsche contact/friction energies, for the discrete-energy diagnostic.
    scalar_type
    nitsche_contact_energy( const cell_type &cl,
                            const CellDegreeInfo< MeshType > &cell_infos,
                            const matrix_type &ET,
                            const vector_type &uTF ) const {
        const auto cb = make_vector_monomial_basis( m_msh, cl, cell_infos.cell_degree() );
        const auto gb = make_sym_matrix_monomial_basis( m_msh, cl, cell_infos.grad_degree() );

        scalar_type energy = scalar_type( 0 );

        const auto fcs = faces( m_msh, cl );
        size_t offset = cb.size();
        const auto fcs_di = cell_infos.facesDegreeInfo();
        size_t face_i = 0;

        const vector_type ET_uTF = ET * uTF;

        for ( auto &fc : fcs ) {
            const auto fdi = fcs_di[face_i++];
            const auto fb = make_vector_monomial_basis( m_msh, fc, fdi.degree() );
            const auto fbs = fb.size();

            if ( m_bnd.is_contact_face( fc ) ) {
                const auto n = normal( m_msh, cl, fc );
                const auto qp_deg = std::max( cell_infos.cell_degree(), cell_infos.grad_degree() );
                const auto qps = integrate( m_msh, fc, 2 * qp_deg + 2 );
                const auto hF = diameter( m_msh, fc );
                const auto gamma_n_F = m_rp.gamma_0_n() / hF;

                for ( auto &qp : qps ) {
                    const scalar_type phi_n = eval_phi_n_uF(
                        fc, ET_uTF, gb, fb, uTF, offset, n, n, gamma_n_F, qp.point() );
                    const scalar_type sig_nn = eval_stress_nn( ET_uTF, gb, n, qp.point() );
                    const scalar_type neg = make_proj_Rmin( phi_n );
                    energy += qp.weight() / ( scalar_type( 2 ) * gamma_n_F ) *
                              ( neg * neg - m_rp.m_theta * sig_nn * sig_nn );
                }
            }
            offset += fbs;
        }
        return energy;
    }

    scalar_type
    nitsche_friction_energy( const cell_type &cl,
                             const CellDegreeInfo< MeshType > &cell_infos,
                             const matrix_type &ET,
                             const vector_type &uTF ) const {
        if ( m_rp.m_frot_type == NO_FRICTION )
            return scalar_type( 0 );

        const auto cb = make_vector_monomial_basis( m_msh, cl, cell_infos.cell_degree() );
        const auto gb = make_sym_matrix_monomial_basis( m_msh, cl, cell_infos.grad_degree() );

        scalar_type energy = scalar_type( 0 );

        const auto fcs = faces( m_msh, cl );
        size_t offset = cb.size();
        const auto fcs_di = cell_infos.facesDegreeInfo();
        size_t face_i = 0;

        const vector_type ET_uTF = ET * uTF;

        for ( auto &fc : fcs ) {
            const auto fdi = fcs_di[face_i++];
            const auto fb = make_vector_monomial_basis( m_msh, fc, fdi.degree() );
            const auto fbs = fb.size();

            if ( m_bnd.is_contact_face( fc ) ) {
                const auto n = normal( m_msh, cl, fc );
                const auto qp_deg = std::max( cell_infos.cell_degree(), cell_infos.grad_degree() );
                const auto qps = integrate( m_msh, fc, 2 * qp_deg + 2 );
                const auto hF = diameter( m_msh, fc );
                const auto gamma_t_F = m_rp.gamma_0_t() / hF;
                const auto gamma_n_F = m_rp.gamma_0_n() / hF;

                for ( auto &qp : qps ) {
                    const vector_static sigma_nt = eval_stress_nt( ET_uTF, gb, n, qp.point() );

                    const scalar_type phi_n = eval_phi_n_uF(
                        fc, ET_uTF, gb, fb, uTF, offset, n, n, gamma_n_F, qp.point() );

                    vector_static proj = vector_static::Zero();
                    if ( phi_n < scalar_type( 0 ) ) {

                        //   u(pt) = sum_i uF_i phi_i(pt) = t_phi^T uF ;  u_t = u - (u.n) n .
                        const auto t_phi = fb.eval_functions( qp.point() );
                        const vector_type uF = uTF.segment( offset, fb.size() );
                        const vector_static u_full = t_phi.transpose() * uF;
                        const vector_static u_t = u_full - u_full.dot( n ) * n;

                        // Displacement-based tangential trial stress and Coulomb projection.
                        const vector_static phi_t = sigma_nt - gamma_t_F * u_t;
                        const scalar_type Fc = m_bnd.contact_boundary_func( fc )( qp.point() );
                        const scalar_type fric_bound = -Fc * make_proj_Rmin( phi_n );
                        const scalar_type ptn = phi_t.norm();
                        proj = ( ptn <= fric_bound ) ? phi_t : ( fric_bound / ptn ) * phi_t;
                    }

                    energy += qp.weight() / ( scalar_type( 2 ) * gamma_n_F ) *
                              ( proj.squaredNorm() - m_rp.m_theta * sigma_nt.squaredNorm() );
                }
            }
            offset += fbs;
        }
        return energy;
    }

    template < typename CellBasis >
    scalar_type
    eval_uT_n( const CellBasis &cb,
               const vector_type &uTF,
               const vector_static &kinematic_normal,
               const point_type &pt ) const {
        const vector_type uT_n = make_hho_u_n( kinematic_normal, cb, pt );

        return uT_n.dot( uTF.head( cb.size() ) );
    }

    template < typename CellBasis >
    vector_static
    eval_uT_t( const CellBasis &cb,
               const vector_type &uTF,
               const vector_static &kinematic_normal,
               const point_type &pt ) const {
        const auto uT_t = make_hho_u_t( kinematic_normal, cb, pt );

        return uT_t.transpose() * ( uTF.head( cb.size() ) );
    }

    template < typename FaceBasis >
    scalar_type
    eval_uF_n( const FaceBasis &fb,
               const vector_type &uF,
               const vector_static &kinematic_normal,
               const point_type &pt ) const {
        const vector_type uF_n = make_hho_u_n( kinematic_normal, fb, pt );
        assert( uF_n.size() == uF.size() );
        return uF_n.dot( uF );
    }

    template < typename FaceBasis >
    vector_static
    eval_uF_t( const FaceBasis &fb,
               const vector_type &uF,
               const vector_static &kinematic_normal,
               const point_type &pt ) const {
        const auto uF_t = make_hho_u_t( kinematic_normal, fb, pt );
        assert( uF_t.rows() == uF.size() );
        return uF_t.transpose() * uF;
    }

    template < typename GradBasis >
    scalar_type
    eval_stress_nn( const vector_type &ET_uTF,
                    const GradBasis &gb,
                    const vector_static &reference_normal,
                    const point_type &pt ) const {
        const auto stress = eval_stress( ET_uTF, gb, pt );
        return ( stress * reference_normal ).dot( reference_normal );
    }

    template < typename GradBasis >
    vector_static
    eval_stress_nt( const vector_type &ET_uTF,
                    const GradBasis &gb,
                    const vector_static &reference_normal,
                    const point_type &pt ) const {
        const matrix_static sig = eval_stress( ET_uTF, gb, pt );
        const vector_static s_n = sig * reference_normal;
        const auto s_nn = s_n.dot( reference_normal );
        return s_n - s_nn * reference_normal;
    }

    template < typename GradBasis, typename CellBasis >
    scalar_type
    eval_phi_n_uT( const face_type &fc,
                   const vector_type &ET_uTF,
                   const GradBasis &gb,
                   const CellBasis &cb,
                   const vector_type &uTF,
                   const vector_static &reference_normal,
                   const vector_static &kinematic_normal,
                   scalar_type gamma_n_F,
                   const point_type &pt ) const {
        const scalar_type sigma_nn = eval_stress_nn( ET_uTF, gb, reference_normal, pt );
        const vector_type uT = uTF.head( cb.size() );

        const auto gap_func = m_bnd.contact_boundary_gap( fc );
        const scalar_type gap =
            priv::compute_gap_fb( m_msh, fc, cb, uT, gap_func, pt, kinematic_normal, m_time );

        return sigma_nn + gamma_n_F * gap;
    }

    template < typename GradBasis, typename CellBasis >
    scalar_type
    eval_proj_phi_n_uT( const face_type &fc,
                        const vector_type &ET_uTF,
                        const GradBasis &gb,
                        const CellBasis &cb,
                        const vector_type &uTF,
                        const vector_static &reference_normal,
                        const vector_static &kinematic_normal,
                        scalar_type gamma_n_F,
                        const point_type &pt ) const {
        const scalar_type phi_n_1_u = eval_phi_n_uT(
            fc, ET_uTF, gb, cb, uTF, reference_normal, kinematic_normal, gamma_n_F, pt );

        if ( phi_n_1_u <= scalar_type( 0 ) )
            return phi_n_1_u;

        return scalar_type( 0 );
    }

    template < typename GradBasis, typename CellBasis >
    vector_static
    eval_phi_t_uT( const vector_type &ET_uTF,
                   const GradBasis &gb,
                   const CellBasis &cb,
                   const vector_type &uTF,
                   const vector_static &reference_normal,
                   const vector_static &kinematic_normal,
                   scalar_type gamma_t_F,
                   const point_type &pt ) const {
        const auto sigma_nt = eval_stress_nt( ET_uTF, gb, reference_normal, pt );
        const auto uT_t = eval_uT_t( cb, uTF, kinematic_normal, pt );

        return sigma_nt - gamma_t_F * uT_t;
    }

    template < typename GradBasis, typename CellBasis >
    vector_static
    eval_proj_phi_t_uT( const vector_type &ET_uTF,
                        const GradBasis &gb,
                        const CellBasis &cb,
                        const vector_type &uTF,
                        const vector_static &reference_normal,
                        const vector_static &kinematic_normal,
                        scalar_type gamma_t_F,
                        scalar_type s,
                        const point_type &pt ) const {
        const vector_static phi_t_1_u =
            eval_phi_t_uT( ET_uTF, gb, cb, uTF, reference_normal, kinematic_normal, gamma_t_F, pt );

        return make_proj_alpha( phi_t_1_u, s );
    }

    template < typename GradBasis, typename FaceBasis >
    scalar_type
    eval_phi_n_uF( const face_type &fc,
                   const vector_type &ET_uTF,
                   const GradBasis &gb,
                   const FaceBasis &fb,
                   const vector_type &uTF,
                   const size_t offset,
                   const vector_static &reference_normal,
                   const vector_static &kinematic_normal,
                   scalar_type gamma_n_F,
                   const point_type &pt ) const {
        const vector_type uF = uTF.segment( offset, fb.size() );
        const scalar_type sigma_nn = eval_stress_nn( ET_uTF, gb, reference_normal, pt );
        const auto gap_func = m_bnd.contact_boundary_gap( fc );
        const scalar_type gap =
            priv::compute_gap_fb( m_msh, fc, fb, uF, gap_func, pt, kinematic_normal, m_time );

        return sigma_nn + gamma_n_F * gap;
    }

    template < typename GradBasis, typename FaceBasis >
    scalar_type
    eval_proj_phi_n_uF( const face_type &fc,
                        const vector_type &ET_uTF,
                        const GradBasis &gb,
                        const FaceBasis &fb,
                        const vector_type &uTF,
                        const size_t offset,
                        const vector_static &reference_normal,
                        const vector_static &kinematic_normal,
                        scalar_type gamma_n_F,
                        const point_type &pt ) const {
        const scalar_type phi_n_1_u = eval_phi_n_uF(
            fc, ET_uTF, gb, fb, uTF, offset, reference_normal, kinematic_normal, gamma_n_F, pt );

        if ( phi_n_1_u <= scalar_type( 0 ) )
            return phi_n_1_u;

        return scalar_type( 0 );
    }

    template < typename GradBasis, typename FaceBasis >
    vector_static
    eval_phi_t_uF( const vector_type &ET_uTF,
                   const GradBasis &gb,
                   const FaceBasis &fb,
                   const vector_type &uTF,
                   size_t offset,
                   const vector_static &reference_normal,
                   const vector_static &kinematic_normal,
                   scalar_type gamma_t_F,
                   const point_type &pt ) const {
        const vector_type uF = uTF.segment( offset, fb.size() );
        const auto sigma_nt = eval_stress_nt( ET_uTF, gb, reference_normal, pt );
        const auto uF_t = eval_uF_t( fb, uF, kinematic_normal, pt );

        return sigma_nt - gamma_t_F * uF_t;
    }

    template < typename GradBasis, typename FaceBasis >
    vector_static
    eval_proj_tresca_phi_t_uF( const vector_type &ET_uTF,
                               const GradBasis &gb,
                               const FaceBasis &fb,
                               const vector_type &uTF,
                               size_t offset,
                               const vector_static &reference_normal,
                               const vector_static &kinematic_normal,
                               scalar_type gamma_t_F,
                               scalar_type s,
                               const point_type &pt ) const {
        const vector_static phi_t_1_u = eval_phi_t_uF(
            ET_uTF, gb, fb, uTF, offset, reference_normal, kinematic_normal, gamma_t_F, pt );

        return make_proj_alpha( phi_t_1_u, s );
    }

    template < typename GradBasis, typename FaceBasis >
    vector_static
    eval_proj_coulomb_phi_t_uF( const face_type &fc,
                                const vector_type &ET_uTF,
                                const GradBasis &gb,
                                const FaceBasis &fb,
                                const vector_type &uTF,
                                const vector_type &vuTF,
                                size_t offset,
                                const vector_static &reference_normal,
                                const vector_static &kinematic_normal,
                                scalar_type gamma_t_F,
                                scalar_type gamma_n_F,
                                scalar_type Fc,
                                const point_type &pt ) const {
        const vector_static phi_t_1_u = eval_phi_t_uF(
            ET_uTF, gb, fb, vuTF, offset, reference_normal, kinematic_normal, gamma_t_F, pt );
        const scalar_type phi_n_1_u = eval_phi_n_uF(
            fc, ET_uTF, gb, fb, uTF, offset, reference_normal, kinematic_normal, gamma_n_F, pt );
        const scalar_type proj_phi_n_1_u = make_proj_Rmin( phi_n_1_u );
        const scalar_type fric_bound = coulomb_bound( Fc, proj_phi_n_1_u );

        return make_proj_alpha( phi_t_1_u, fric_bound );
    }

    // cell-trace counterpart of eval_proj_coulomb_phi_t_uF, hence no offset
    template < typename GradBasis, typename CellBasis >
    vector_static
    eval_proj_coulomb_phi_t_uT( const face_type &fc,
                                const vector_type &ET_uTF,
                                const GradBasis &gb,
                                const CellBasis &cb,
                                const vector_type &uTF,
                                const vector_type &vuTF,
                                const vector_static &reference_normal,
                                const vector_static &kinematic_normal,
                                scalar_type gamma_t_F,
                                scalar_type gamma_n_F,
                                scalar_type Fc,
                                const point_type &pt ) const {
        const vector_static phi_t_1_u = eval_phi_t_uT(
            ET_uTF, gb, cb, vuTF, reference_normal, kinematic_normal, gamma_t_F, pt );
        const scalar_type phi_n_1_u = eval_phi_n_uT(
            fc, ET_uTF, gb, cb, uTF, reference_normal, kinematic_normal, gamma_n_F, pt );
        const scalar_type proj_phi_n_1_u = make_proj_Rmin( phi_n_1_u );
        const scalar_type fric_bound = coulomb_bound( Fc, proj_phi_n_1_u );

        return make_proj_alpha( phi_t_1_u, fric_bound );
    }

    template < typename GradBasis >
    matrix_static eval_stress( const vector_type &ET_uTF, const GradBasis &gb,
                               const point_type pt ) const {
        const auto gphi = gb.eval_functions( pt );

        const matrix_static Eu = disk::eval( ET_uTF, gphi );

        return 2 * m_material_data.getMu() * Eu +
               m_material_data.getLambda() * Eu.trace() * matrix_static::Identity();
    }

    void
    setTime( const scalar_type time ) {
        m_time = time;
    }
};

} // namespace mechanics

} // end namespace disk
