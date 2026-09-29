/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
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

// Explicit iteration

#pragma once

#include "diskpp/mechanics/NewtonSolver/GenericIteration.hpp"

#include <set>

namespace disk {

namespace mechanics {

/**
 * @brief Newton-Raphson iteration for nonlinear solid mechanics
 *
 *  Specialized for HHO methods
 *
 *  Options :  - small and finite deformations
 *             - plasticity, hyperelasticity (various laws)
 *
 * @tparam MeshType type of the mesh
 */
template < typename MeshType >
class ExplicitIteration : public GenericIteration< MeshType > {
    typedef typename GenericIteration< MeshType >::mesh_type mesh_type;
    typedef typename GenericIteration< MeshType >::cell_type cell_type;
    typedef typename GenericIteration< MeshType >::scalar_type scalar_type;

    typedef typename GenericIteration< MeshType >::matrix_type matrix_type;
    typedef typename GenericIteration< MeshType >::vector_type vector_type;

    typedef typename GenericIteration< MeshType >::param_type param_type;
    typedef typename GenericIteration< MeshType >::bnd_type bnd_type;
    typedef typename GenericIteration< MeshType >::behavior_type behavior_type;

    typedef typename GenericIteration< MeshType >::elem_type elem_type;
    typedef typename GenericIteration< MeshType >::assembler_type assembler_type;
    typedef typename GenericIteration< MeshType >::func_type func_type;

    std::vector< vector_type > _vite_pred;

    // Boundary condition to assemble a problem with acceleration variable and not
    // displacement - only CLAMPED dirichlet are used.
    bnd_type _bnd_acce;

  public:
    ExplicitIteration( const mesh_type &msh,
                       const bnd_type &bnd,
                       const param_type &rp,
                       const MeshDegreeInfo< mesh_type > &degree_infos,
                       const std::shared_ptr< solvers::sparse_solver< scalar_type > > lin_solv,
                       const TimeStep< scalar_type > &current_step )
        : GenericIteration< MeshType >( msh, bnd, rp, degree_infos, lin_solv, current_step ),
          _bnd_acce( msh ) {
        if ( rp.getUnsteadyScheme() != DynamicType::LEAP_FROG ) {
            throw std::invalid_argument( "Time scheme not supported by Explicit solver." );
        }

        std::set< std::size_t > bc_id;

        for ( auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++ ) {
            const auto bfc = *itor;
            const auto f_id = msh.lookup( bfc );

            if ( bnd.dirichlet_boundary_type( f_id ) == DirichletType::CLAMPED ) {
                const auto b_id = bnd.dirichlet_boundary_id( f_id );
                if ( !bc_id.contains( b_id ) ) {
                    _bnd_acce.addDirichletBC(
                        DirichletType::CLAMPED, b_id, bnd.dirichlet_boundary_func( f_id ) );
                    bc_id.insert( b_id );
                }
            }
        }

        this->m_assembler = assembler_type( msh, degree_infos, _bnd_acce );
    }

    InitInfo
    initialize( const mesh_type &msh,
                const bnd_type &bnd,
                const param_type &rp,
                const MeshDegreeInfo< mesh_type > &degree_infos,
                NonLinearData< scalar_type > &data,
                behavior_type &behavior,
                const StabCoeffManager< scalar_type > &stab_manager,
                MultiTimeField< scalar_type > &fields ) override {
        InitInfo ii;

        if ( !this->m_lin_solv->is_factorized() ) {
            this->m_assembler.initialize();

            timecounter tc;
            tc.tic();

            const bool mixed_order = rp.m_cell_degree > rp.m_face_degree;
            const bool small_def = ( behavior.getDeformation() == SMALL_DEF );

            data.m_lhs_loc = std::make_shared< std::map< int, matrix_type > >();

            const auto rho = this->m_dyna.getParam( "rho" );

            for ( auto &cl : msh ) {
                const auto cell_i = msh.lookup( cl );

                const auto cell_infos = degree_infos.cellDegreeInfo( msh, cl );
                const auto num_cell_dofs = vector_cell_dofs( msh, cell_infos );

                const auto faces_infos = cell_infos.facesDegreeInfo();
                const auto num_faces_dofs = vector_faces_dofs( msh, faces_infos );

                matrix_type lhs;

                switch ( rp.getNonLinearSolver() ) {
                case NonLinearSolverType::EXPLICIT: {

                    const auto hT = diameter( msh, cl );

                    // matrice de mass + stab
                    // maybe compute an other stabilization without hF^-1 weight
                    lhs = ( rho * hT * hT ) *
                          _stab( msh, cl, rp, degree_infos, data.m_stab_precomputed );

                    lhs.topLeftCorner( num_cell_dofs, num_cell_dofs ) +=
                        this->m_dyna.mass_matrix( msh, cl, degree_infos );
                    break;
                }
                default: {
                    throw std::invalid_argument( "Explicit option is unknown." );
                    break;
                }
                }

                const vector_type rhs = vector_type::Zero( num_cell_dofs + num_faces_dofs );

                ( *data.m_lhs_loc )[cell_i] = lhs;

                const auto scnp = make_vector_static_condensation_withMatrix(
                    msh, cl, degree_infos, lhs, rhs, true );
                const auto &lc = std::get< 0 >( scnp );

                this->m_assembler.assemble( msh, cl, _bnd_acce, lc.first, lc.second );
            }
            tc.toc();
            ii.m_time_dyna += tc.elapsed();

            this->m_assembler.finalize();

            tc.tic();
            const auto status = this->m_lin_solv->factorize( this->m_assembler.LHS );
            tc.toc();
            ii.m_time_solve += tc.elapsed();

            if ( status != solvers::direct_solver_status::ok ) {
                throw std::runtime_error( "Fail to factorize the matrix." );
            }
        }

        // Compute current displacement and predicted velocity

        const auto depl_prev = fields.getField( -1, FieldName::DEPL );
        const auto vite_prev = fields.getField( -1, FieldName::VITE );
        const auto acce_prev = fields.getField( -1, FieldName::ACCE );

        const scalar_type dt = this->m_time_step.increment_time();
        const scalar_type dts2 = 0.5 * dt;

        std::vector< vector_type > depl, depl_cells;
        depl.resize( msh.cells_size() );
        depl_cells.resize( msh.cells_size() );
        _vite_pred.clear();
        _vite_pred.resize( msh.cells_size() );

        vector_type huT, hvT_pred;

        for ( auto &cl : msh ) {
            const auto c_id = msh.lookup( cl );
            const auto cell_infos = degree_infos.cellDegreeInfo( msh, cl );
            const auto num_cell_dofs = vector_cell_dofs( msh, cell_infos );

            switch ( rp.getUnsteadyScheme() ) {
            case DynamicType::LEAP_FROG: {
                hvT_pred = vite_prev[c_id] + dts2 * acce_prev[c_id];
                apply_dirichlet( msh, cl, cell_infos, _bnd_acce, hvT_pred );
                huT = depl_prev[c_id] + dt * hvT_pred;
                apply_dirichlet( msh, cl, cell_infos, bnd, huT );
            } break;

            default:
                throw std::invalid_argument( "Time scheme not implemented." );
                break;
            }

            _vite_pred[c_id] = hvT_pred;
            depl[c_id] = huT;
            depl_cells[c_id] = huT.head( num_cell_dofs );
        }

        fields.setCurrentField( FieldName::DEPL, depl );
        fields.setCurrentField( FieldName::DEPL_CELLS, depl_cells );

        return ii;
    }

    AssemblyInfo
    assemble( const mesh_type &msh,
              const bnd_type &bnd,
              const param_type &rp,
              const MeshDegreeInfo< mesh_type > &degree_infos,
              const std::unique_ptr< func_type > &lf,
              NonLinearData< scalar_type > &data,
              behavior_type &behavior,
              StabCoeffManager< scalar_type > &stab_manager,
              MultiTimeField< scalar_type > &fields ) override {
        elem_type elem;
        AssemblyInfo ai;

        // set RHS to zero
        this->m_assembler.initialize();
        this->m_F_int = 0.0;

        const bool small_def = ( behavior.getDeformation() == SMALL_DEF );

        // Like if it is an implicit scheme
        auto current_time = this->m_time_step.end_time();
        auto depl = fields.getCurrentField( FieldName::DEPL );

        const auto lhs_loc = data.m_lhs_loc;

        const auto rlf = this->_getLoad( lf, current_time );

        timecounter tc, ttot;

        ttot.tic();

        for ( auto &cl : msh ) {
            const auto cell_i = msh.lookup( cl );

            const auto huT = depl.at( cell_i );
            const auto num_tot_dofs = huT.size();

            vector_type hvT = vector_type::Zero( num_tot_dofs );

            if ( this->m_dyna.enable() ) {
                switch ( rp.getUnsteadyScheme() ) {
                case DynamicType::LEAP_FROG: {
                    hvT = _vite_pred[cell_i];
                } break;

                default:
                    throw std::invalid_argument( "Time scheme not implemented." );
                    break;
                }
            }

            // Gradient Reconstruction
            // std::cout << "Grad" << std::endl;
            tc.tic();
            matrix_type GT =
                _gradrec( msh, cl, rp, degree_infos, small_def, data.m_gradient_precomputed );
            tc.toc();
            ai.m_time_gradrec += tc.elapsed();

            // Mechanical Computation

            tc.tic();
            // std::cout << "Elem" << std::endl;
            elem.compute( msh,
                          cl,
                          bnd,
                          rp,
                          degree_infos,
                          rlf,
                          GT,
                          huT,
                          hvT,
                          this->m_time_step,
                          behavior,
                          stab_manager,
                          small_def,
                          false );

            vector_type rhs = elem.RTF;

            tc.toc();
            ai.m_time_elem += tc.elapsed();

            tc.tic();
            if ( rp.m_stab ) {
                const auto beta_s = stab_manager.getValue( msh, cl );

                vector_type suT = _stab( msh, cl, rp, degree_infos, data.m_stab_precomputed ) * huT;

                rhs -= beta_s * suT;
            }
            tc.toc();
            ai.m_time_stab += tc.elapsed();

            tc.tic();
            const auto scnp = make_vector_static_condensation_withMatrix(
                msh, cl, degree_infos, ( *lhs_loc )[cell_i], rhs, true );

            const auto &lc = std::get< 0 >( scnp );
            this->m_AL[cell_i] = std::get< 1 >( scnp );
            this->m_bL[cell_i] = std::get< 2 >( scnp );

            tc.toc();
            ai.m_time_statcond += tc.elapsed();

            tc.tic();
            this->m_assembler.assemble_rhs( msh, cl, _bnd_acce, lc.first, lc.second );
            tc.toc();
            ai.m_time_assembler += tc.elapsed();
        }

        ai.m_time_law += elem.time_law;
        ai.m_time_contact += elem.time_contact;
        ai.m_time_load += elem.time_load;
        ai.m_time_rigi += elem.time_rigi;
        ai.m_time_fint += elem.time_fint;

        tc.tic();
        this->m_assembler.impose_neumann_boundary_conditions( msh, bnd );
        this->m_assembler.finalize();
        tc.toc();
        ai.m_time_assembler += tc.elapsed();

        ttot.toc();
        ai.m_time_assembly = ttot.elapsed();
        ai.m_linear_system_size = this->m_assembler.LHS.rows();
        return ai;
    }

    SolveInfo
    solve() override {
        timecounter tc;

        // std::cout << "RHS" << this->m_assembler.RHS.transpose() << std::endl;
        // The soluton is acce_curr and not depl_curr.
        tc.tic();
        const auto status = this->m_lin_solv->solve( this->m_assembler.RHS, this->m_system_displ );
        tc.toc();

        if ( status != solvers::direct_solver_status::ok ) {
            throw std::runtime_error( "Error during linear system resolution" );
        }

        // std::cout << "SOL" << this->m_system_displ.transpose() << std::endl;

        return SolveInfo(
            this->m_assembler.LHS.rows(), this->m_assembler.LHS.nonZeros(), tc.elapsed() );
    }

    scalar_type
    postprocess( const mesh_type &msh,
                 const bnd_type &bnd,
                 const param_type &rp,
                 const MeshDegreeInfo< mesh_type > &degree_infos,
                 const std::unique_ptr< func_type > &lf,
                 NonLinearData< scalar_type > &data,
                 behavior_type &behavior,
                 StabCoeffManager< scalar_type > &stab_manager,
                 MultiTimeField< scalar_type > &fields ) override {
        timecounter tc;
        tc.tic();

        // Get acceleration
        const auto [aFh, idx] =
            this->m_assembler.expand_solution( msh, _bnd_acce, this->m_system_displ );

        const scalar_type dt = this->m_time_step.increment_time();
        const scalar_type dts2 = 0.5 * dt;

        std::vector< vector_type > acce, acce_cells, vite, vite_cells;
        acce.resize( msh.cells_size() );
        acce_cells.resize( msh.cells_size() );
        vite.resize( msh.cells_size() );
        vite_cells.resize( msh.cells_size() );

        vector_type hvT, haT;

        for ( auto &cl : msh ) {
            const auto c_id = msh.lookup( cl );
            const auto cell_infos = degree_infos.cellDegreeInfo( msh, cl );
            const auto faces_infos = cell_infos.facesDegreeInfo();
            const auto num_cell_dofs = vector_cell_dofs( msh, cell_infos );
            const auto num_faces_dofs = vector_faces_dofs( msh, faces_infos );

            vector_type adT = vector_type( num_faces_dofs );

            const auto fcs_id = faces_id( msh, cl );
            size_t face_offset = 0;
            for ( size_t face_i = 0; face_i < fcs_id.size(); face_i++ ) {
                const size_t face_id = fcs_id[face_i];
                const auto num_face_dofs = idx( face_id + 1 ) - idx( face_id );

                adT.segment( face_offset, num_face_dofs ) =
                    aFh.segment( idx( face_id ), num_face_dofs );
                face_offset += num_face_dofs;
            }

            // static decondensation
            const vector_type aT = this->m_bL[c_id] - this->m_AL[c_id] * adT;

            vector_type haT = vector_type::Zero( num_cell_dofs + num_faces_dofs );
            haT.head( num_cell_dofs ) = aT;
            haT.tail( num_faces_dofs ) = adT;

            switch ( rp.getUnsteadyScheme() ) {
            case DynamicType::LEAP_FROG: {
                hvT = _vite_pred[c_id] + dts2 * haT;
            } break;
            default:
                throw std::invalid_argument( "Time scheme not implemented." );
                break;
            }

            acce[c_id] = haT;
            acce_cells[c_id] = aT;
            vite[c_id] = hvT;
            vite_cells[c_id] = hvT.head( num_cell_dofs );
        }

        fields.setCurrentField( FieldName::VITE, vite );
        fields.setCurrentField( FieldName::VITE_CELLS, vite_cells );
        fields.setCurrentField( FieldName::ACCE, acce );
        fields.setCurrentField( FieldName::ACCE_CELLS, acce_cells );

        tc.toc();
        return tc.elapsed();
    }
};
} // namespace mechanics
} // namespace disk
