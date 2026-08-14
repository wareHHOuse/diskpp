/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019, 2025               nicolas.pignet@enpc.fr
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

#include "diskpp/bases/bases.hpp"
#include "diskpp/common/eigen.hpp"
#include "diskpp/common/timecounter.hpp"
#include "diskpp/mechanics/NewtonSolver/Fields.hpp"
#include "diskpp/mechanics/NewtonSolver/NonLinearParameters.hpp"
#include "diskpp/mechanics/NewtonSolver/TimeManager.hpp"
#include "diskpp/methods/hho"
#include "diskpp/quadratures/quadratures.hpp"

#include <cassert>

namespace disk {

namespace mechanics {

template < typename T >
int getNumberOfStepToSave( const NonLinearParameters< T > &rp ) {
    if ( rp.isUnsteady() ) {
        switch ( rp.getUnsteadyScheme() ) {
        case DynamicType::NEWMARK:
        case DynamicType::THETA:
        case DynamicType::BACKWARD_EULER:
        case DynamicType::CRANK_NICOLSON: {
            return 2;
            break;
        }
        case DynamicType::LEAP_FROG: {
            return 3;
            break;
        }
        default:
            throw std::invalid_argument( "Unsupported dynamic scheme" );
            break;
        }
    }

    return 2;
}

template < typename T >
void reformulation_dynamic( NonLinearParameters< T > &rp ) {
    if ( rp.isUnsteady() ) {
        switch ( rp.getUnsteadyScheme() ) {
        case DynamicType::BACKWARD_EULER: {
            rp.setUnsteadyScheme( DynamicType::THETA );
            std::map< std::string, T > dyna_para;
            dyna_para["theta"] = 1.0;
            rp.setUnsteadyParameters( dyna_para );
            break;
        }
        case DynamicType::CRANK_NICOLSON: {
            rp.setUnsteadyScheme( DynamicType::THETA );
            std::map< std::string, T > dyna_para;
            dyna_para["theta"] = 0.5;
            rp.setUnsteadyParameters( dyna_para );
            break;
        }
        default:
            break;
        }
    }
}

template < typename MeshType >
class dynamic_computation {
    typedef MeshType mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell cell_type;

    typedef dynamic_matrix< scalar_type > matrix_type;
    typedef dynamic_vector< scalar_type > vector_type;

    std::map< std::string, scalar_type > m_param;
    DynamicType m_scheme;

    std::vector< vector_type > m_acce_cells_pred, m_acce_pred;

  public:
    matrix_type K_iner;
    vector_type R_iner;

    double time_dyna;

    dynamic_computation() : m_scheme( DynamicType::STATIC ) {}

    dynamic_computation( const NonLinearParameters< scalar_type > &rp ) {
        m_param = rp.getUnsteadyParameters();
        m_scheme = rp.getUnsteadyScheme();

        // for (auto &[key, val] : m_param) {
        //     std::cout << key << " : " << val << std::endl;
        // }
    }

    auto
    getParam( std::string param ) {
        return m_param.at( param );
    }

    bool enable( void ) const { return m_scheme != DynamicType::STATIC; }

    bool isExplicit() const { return m_scheme == DynamicType::LEAP_FROG; }

    void prediction( const mesh_type &msh, const MeshDegreeInfo< mesh_type > &degree_infos,
                     const TimeStep< scalar_type > &time_step,
                     MultiTimeField< scalar_type > &fields ) {
        if ( this->enable() ) {
            m_acce_cells_pred.clear();
            m_acce_cells_pred.resize( msh.cells_size() );
            m_acce_pred.clear();
            m_acce_pred.resize( msh.cells_size() );

            switch ( m_scheme ) {
            case DynamicType::NEWMARK: {
                std::vector< vector_type > acce, acce_cells;
                acce.resize( msh.cells_size() );
                acce_cells.resize( msh.cells_size() );

                const auto depl_prev = fields.getField( -1, FieldName::DEPL );
                const auto velo_prev = fields.getField( -1, FieldName::VITE );
                const auto acce_prev = fields.getField( -1, FieldName::ACCE );

                const auto beta = m_param.at( "beta" );
                const scalar_type dt = time_step.increment_time();
                const scalar_type cd = 1.0 / ( beta * dt * dt );
                const scalar_type cv = 1.0 / ( beta * dt );
                const scalar_type ca = ( 1.0 - 2.0 * beta ) / ( 2.0 * beta );

                for ( auto &cl : msh ) {
                    const auto c_id = msh.lookup( cl );
                    const auto c_idnfos = degree_infos.cellDegreeInfo( msh, cl );
                    const auto num_cell_dofs = vector_cell_dofs( msh, c_idnfos );

                    const vector_type acce_curr = cv * velo_prev[c_id] + ca * acce_prev[c_id];
                    const vector_type acce_pred = cd * depl_prev[c_id] + acce_curr;
                    m_acce_pred[c_id] = acce_pred;
                    m_acce_cells_pred[c_id] = acce_pred.head( num_cell_dofs );
                    acce[c_id] = -acce_curr;
                    acce_cells[c_id] = -acce_curr.head( num_cell_dofs );
                }

                fields.setCurrentField( FieldName::ACCE, acce );
                fields.setCurrentField( FieldName::ACCE_CELLS, acce_cells );
                break;
            }
            case DynamicType::THETA: {
                std::vector< vector_type > acce, acce_cells;
                acce.resize( msh.cells_size() );
                acce_cells.resize( msh.cells_size() );

                const auto depl_prev = fields.getField( -1, FieldName::DEPL );
                const auto velo_prev = fields.getField( -1, FieldName::VITE );
                const auto acce_prev = fields.getField( -1, FieldName::ACCE );

                const scalar_type dt = time_step.increment_time();
                const auto theta = m_param.at( "theta" );

                const scalar_type cd = 1.0 / ( theta * theta * dt * dt );
                const scalar_type cv = 1.0 / ( theta * theta * dt );
                const scalar_type ca = ( 1.0 - theta ) / theta;

                for ( auto &cl : msh ) {
                    const auto c_idnfos = degree_infos.cellDegreeInfo( msh, cl );
                    const auto num_cell_dofs = vector_cell_dofs( msh, c_idnfos );

                    const auto c_id = msh.lookup( cl );
                    const vector_type acce_curr = cv * velo_prev[c_id] + ca * acce_prev[c_id];
                    const vector_type acce_pred = cd * depl_prev[c_id] + acce_curr;
                    m_acce_pred[c_id] = acce_pred;
                    m_acce_cells_pred[c_id] = acce_pred.head( num_cell_dofs );
                    acce[c_id] = -acce_curr;
                    acce_cells[c_id] = -acce_curr.head( num_cell_dofs );
                }

                fields.setCurrentField( FieldName::ACCE, acce );
                fields.setCurrentField( FieldName::ACCE_CELLS, acce_cells );
                break;
            }
            case DynamicType::LEAP_FROG: {

                std::vector< vector_type > vite, acce, depl;
                vite.resize( msh.cells_size() );
                acce.resize( msh.cells_size() );
                depl.resize( msh.cells_size() );

                const auto tf2 = fields.getTimeField( -2 );

                const auto depl_prev = fields.getField( -1, FieldName::DEPL_CELLS );
                const auto vite_prev = fields.getField( -1, FieldName::VITE_CELLS );

                auto depl_curr = fields.getCurrentField( FieldName::DEPL );

                const scalar_type dt = time_step.increment_time();
                const scalar_type dt2 = dt * dt;
                const scalar_type dt2s2 = dt * dt / 2.0;
                const scalar_type un_dt = 1.0 / dt;
                const scalar_type un_dt2 = 1.0 / ( dt * dt );

                if ( tf2.empty() ) {
                    const auto acce_prev = fields.getField( -1, FieldName::ACCE_CELLS );

                    for ( auto &cl : msh ) {
                        const auto c_id = msh.lookup( cl );
                        const auto uT = depl_prev.at( c_id ) + dt * vite_prev.at( c_id ) +
                                        dt2s2 * acce_prev.at( c_id );
                        const auto vT = vite_prev.at( c_id ) + dt * acce_prev.at( c_id );
                        const auto aT = acce_prev.at( c_id );

                        depl_curr[c_id].head( uT.size() ) = uT;

                        acce[c_id] = aT;
                        vite[c_id] = vT;
                        depl[c_id] = uT;
                    }
                } else {
                    const auto resi_prev = fields.getField( -1, FieldName::RESI_CELLS );
                    const auto depl_pprev = tf2.getField( FieldName::DEPL_CELLS );

                    for ( auto &cl : msh ) {
                        const auto c_id = msh.lookup( cl );
                        const matrix_type mm = this->mass_matrix( msh, cl, degree_infos );

                        const vector_type depl_pred =
                            2.0 * depl_prev.at( c_id ) - depl_pprev.at( c_id );

                        const vector_type uT =
                            dt2 * ( mm.ldlt().solve( resi_prev.at( c_id ) ) ) + depl_pred;
                        const auto vT = un_dt * ( uT - depl_prev.at( c_id ) );
                        const auto aT = un_dt * ( vT - vite_prev.at( c_id ) );

                        acce[c_id] = aT;
                        vite[c_id] = vT;
                        depl[c_id] = uT;
                        depl_curr[c_id].head( uT.size() ) = uT;
                    }
                }

                fields.setCurrentField( FieldName::VITE_CELLS, vite );
                fields.setCurrentField( FieldName::ACCE_CELLS, acce );
                fields.setCurrentField( FieldName::DEPL_CELLS, depl );
                fields.setCurrentField( FieldName::DEPL, depl_curr );
                break;
            }
            default: {
                std::runtime_error( "Scheme not implemented for prediction" );
                break;
            }
            }
        }
    }

    scalar_type postprocess( const mesh_type &msh, const TimeStep< scalar_type > &time_step,
                             MultiTimeField< scalar_type > &fields ) const {
        timecounter tc;
        tc.tic();

        if ( this->enable() ) {
            switch ( m_scheme ) {
            case DynamicType::NEWMARK: {
                std::vector< vector_type > vite, vite_cells, acce, acce_cells;
                vite.resize( msh.cells_size() );
                acce.resize( msh.cells_size() );
                vite_cells.resize( msh.cells_size() );
                acce_cells.resize( msh.cells_size() );

                const auto vite_prev = fields.getField( -1, FieldName::VITE );
                const auto acce_prev = fields.getField( -1, FieldName::ACCE );

                const auto depl_cells = fields.getCurrentField( FieldName::DEPL_CELLS );
                const auto depl = fields.getCurrentField( FieldName::DEPL );

                auto beta = m_param.at( "beta" );
                auto gamma = m_param.at( "gamma" );
                scalar_type dt = time_step.increment_time();

                const scalar_type g0 = ( 1.0 - gamma ) * dt, g1 = gamma * dt;
                const scalar_type cd = 1.0 / ( beta * dt * dt );

                for ( auto &cl : msh ) {
                    const auto c_id = msh.lookup( cl );
                    const auto num_cell_dofs = depl_cells.at( c_id ).size();

                    auto aT = cd * depl.at( c_id ) - m_acce_pred[c_id];
                    auto vT = vite_prev.at( c_id ) + g0 * acce_prev.at( c_id ) + g1 * aT;

                    vite[c_id] = vT;
                    acce[c_id] = aT;
                    vite_cells[c_id] = vT.head( num_cell_dofs );
                    acce_cells[c_id] = aT.head( num_cell_dofs );
                }

                fields.setCurrentField( FieldName::VITE, vite );
                fields.setCurrentField( FieldName::ACCE, acce );
                fields.setCurrentField( FieldName::VITE_CELLS, vite_cells );
                fields.setCurrentField( FieldName::ACCE_CELLS, acce_cells );
                break;
            }
            case DynamicType::THETA: {
                std::vector< vector_type > vite, vite_cells, acce, acce_cells;
                vite.resize( msh.cells_size() );
                acce.resize( msh.cells_size() );
                vite_cells.resize( msh.cells_size() );
                acce_cells.resize( msh.cells_size() );

                const auto vite_prev = fields.getField( -1, FieldName::VITE );
                const auto acce_prev = fields.getField( -1, FieldName::ACCE );

                const auto depl_cells = fields.getCurrentField( FieldName::DEPL_CELLS );
                const auto depl = fields.getCurrentField( FieldName::DEPL );

                scalar_type dt = time_step.increment_time();
                const auto theta = m_param.at( "theta" );

                const scalar_type cd = 1.0 / ( theta * theta * dt * dt );

                const scalar_type v0 = ( 1.0 - theta ) * dt;
                const scalar_type v1 = theta * dt;

                for ( auto &cl : msh ) {
                    const auto c_id = msh.lookup( cl );
                    const auto num_cell_dofs = depl_cells.at( c_id ).size();

                    auto aT = cd * depl.at( c_id ) - m_acce_pred[c_id];
                    auto vT = vite_prev.at( c_id ) + v0 * acce_prev.at( c_id ) + v1 * aT;

                    vite[c_id] = vT;
                    acce[c_id] = aT;
                    vite_cells[c_id] = vT.head( num_cell_dofs );
                    acce_cells[c_id] = aT.head( num_cell_dofs );
                }

                fields.setCurrentField( FieldName::VITE, vite );
                fields.setCurrentField( FieldName::ACCE, acce );
                fields.setCurrentField( FieldName::VITE_CELLS, vite_cells );
                fields.setCurrentField( FieldName::ACCE_CELLS, acce_cells );
                break;
            }
            case DynamicType::LEAP_FROG: {
                break;
            }
            default: {
                std::runtime_error( "Scheme not implemented for post-processing" );
                break;
            }
            }
        }

        tc.toc();
        return tc.elapsed();
    }

    void compute( const mesh_type &msh, const cell_type &cl,
                  const MeshDegreeInfo< mesh_type > &degree_infos, const vector_type &uTF,
                  const vector_type &aT, const TimeStep< scalar_type > &time_step ) {
        // Unsteady Computation

        time_dyna = 0.0;
        timecounter tc;

        tc.tic();
        if ( this->enable() ) {
            const auto c_id = msh.lookup( cl );

            const auto c_idnfos = degree_infos.cellDegreeInfo( msh, cl );
            const auto cell_degree = c_idnfos.cell_degree();

            const auto faces_infos = c_idnfos.facesDegreeInfo();

            const auto num_cell_dofs =
                vector_basis_size( cell_degree, mesh_type::dimension, mesh_type::dimension );
            const auto num_faces_dofs = vector_faces_dofs( msh, faces_infos );
            const auto num_total_dofs = num_cell_dofs + num_faces_dofs;

            const vector_type uT = uTF.head( num_cell_dofs );

            R_iner = vector_type::Zero( num_cell_dofs );
            K_iner = matrix_type::Zero( num_cell_dofs, num_cell_dofs );

            switch ( m_scheme ) {
            case DynamicType::NEWMARK:
            case DynamicType::THETA: {

                auto dt = time_step.increment_time();
                scalar_type c0 = 1.0;
                if ( m_scheme == DynamicType::NEWMARK ) {
                    c0 = m_param.at( "beta" );
                } else if ( m_scheme == DynamicType::THETA ) {
                    const scalar_type theta = m_param.at( "theta" );
                    c0 = theta * theta;
                }

                const matrix_type mass_mat = this->mass_matrix( msh, cl, degree_infos );

                const auto coeff = 1.0 / ( c0 * dt * dt );

                K_iner = coeff * mass_mat;

                R_iner -= mass_mat * aT;
                break;
            }
            case DynamicType::LEAP_FROG: {
                break;
            }
            default: {
                std::runtime_error( "Scheme not implemented for inertial forces" );
                break;
            }
            }
        }
        tc.toc();
        time_dyna += tc.elapsed();
    }

    matrix_type mass_matrix( const mesh_type &msh, const cell_type &cl,
                             const MeshDegreeInfo< mesh_type > &degree_infos ) const {

        const auto c_idnfos = degree_infos.cellDegreeInfo( msh, cl );
        const auto cell_degree = c_idnfos.cell_degree();

        const auto cb = make_vector_monomial_basis( msh, cl, cell_degree );

        auto rho = m_param.at( "rho" );
        return rho * make_mass_matrix( msh, cl, cb );
    };
};

} // namespace mechanics

} // namespace disk