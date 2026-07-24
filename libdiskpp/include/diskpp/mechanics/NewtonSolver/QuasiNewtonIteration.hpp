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

// QuasiNewton iteration

#pragma once

#include "diskpp/mechanics/NewtonSolver/GenericIteration.hpp"

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
class QuasiNewtonIteration : public GenericIteration< MeshType > {
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

    matrix_type _mass_term( const mesh_type &msh, const cell_type &cl,
                            const MeshDegreeInfo< mesh_type > &degree_infos ) const {

        const auto cell_infos = degree_infos.cellDegreeInfo( msh, cl );
        const auto faces_infos = cell_infos.facesDegreeInfo();
        const auto num_faces_dofs = vector_faces_dofs( msh, faces_infos );

        matrix_type mm = matrix_type::Zero( num_faces_dofs, num_faces_dofs );

        const auto fcs = faces( msh, cl );

        auto to_vector = []( const matrix_type &scalar_matrix ) {
            const int scal_total_dofs = scalar_matrix.rows();
            const int vect_total_tofs = scal_total_dofs * mesh_type::dimension;

            matrix_type mm = matrix_type::Zero( vect_total_tofs, vect_total_tofs );

            for ( int i = 0; i < scal_total_dofs; i++ ) {
                const auto row = i * mesh_type::dimension;
                for ( int j = 0; j < scal_total_dofs; j++ ) {
                    const auto col = j * mesh_type::dimension;
                    for ( int k = 0; k < mesh_type::dimension; k++ ) {
                        mm( row + k, col + k ) = scalar_matrix( i, j );
                    }
                }
            }

            return mm;
        };

        int offset = 0;
        for ( size_t i = 0; i < fcs.size(); i++ ) {
            const auto fdi = faces_infos[i];

            if ( fdi.hasUnknowns() ) {
                const auto fc = fcs[i];
                const auto facdeg = fdi.degree();
                const auto hF = diameter( msh, fc );
                const auto fb = make_scalar_monomial_basis( msh, fc, facdeg );
                const auto fbs =
                    vector_basis_size( facdeg, mesh_type::dimension - 1, mesh_type::dimension );

                const matrix_type mass_F = make_mass_matrix( msh, fc, fb );

                mm.block( offset, offset, fbs, fbs ) = ( 1.0 / hF ) * to_vector( mass_F );

                offset += fbs;
            }
        }
        assert( offset == num_faces_dofs );

        return mm;
    }

  public:
    QuasiNewtonIteration( const mesh_type &msh, const bnd_type &bnd, const param_type &rp,
                          const MeshDegreeInfo< mesh_type > &degree_infos,
                          const std::shared_ptr< solvers::sparse_solver< scalar_type > > lin_solv,
                          const TimeStep< scalar_type > &current_step )
        : GenericIteration< MeshType >( msh, bnd, rp, degree_infos, lin_solv, current_step ) {
        if ( rp.getUnsteadyScheme() != DynamicType::LEAP_FROG ) {
            throw std::invalid_argument( "Sheme not supported by QuasiNewton" );
        }
    }

    InitInfo initialize( const mesh_type &msh, const bnd_type &bnd, const param_type &rp,
                         const MeshDegreeInfo< mesh_type > &degree_infos,
                         NonLinearData< scalar_type > &data, behavior_type &behavior,
                         const StabCoeffManager< scalar_type > &stab_manager,
                         MultiTimeField< scalar_type > &fields ) override {
        InitInfo ii = GenericIteration< MeshType >::initialize( msh, bnd, rp, degree_infos, data,
                                                                behavior, stab_manager, fields );

        const bool use_tangent = rp.getNonLinearSolver() == NonLinearSolverType::QNEWTON_BDIAG_JACO;
        if ( use_tangent ) {
            this->m_lin_solv->reset();
        }

        if ( !this->m_lin_solv->is_factorized() ) {
            elem_type elem;
            vector_mechanics_hho_assembler assembler( msh, degree_infos, bnd );

            timecounter tc;
            tc.tic();

            const bool mixed_order = rp.m_cell_degree > rp.m_face_degree;
            const bool small_def = ( behavior.getDeformation() == SMALL_DEF );

            data.m_lhs_loc = std::make_shared< std::map< int, matrix_type > >();
            auto depl = fields.getCurrentField( FieldName::DEPL );

            for ( auto &cl : msh ) {
                const auto cell_i = msh.lookup( cl );

                const auto cell_infos = degree_infos.cellDegreeInfo( msh, cl );
                const auto num_cell_dofs = vector_cell_dofs( msh, cell_infos );

                const auto faces_infos = cell_infos.facesDegreeInfo();
                const auto num_faces_dofs = vector_faces_dofs( msh, faces_infos );

                const auto beta_s = stab_manager.getValue( msh, cl );

                matrix_type lhs;

                switch ( rp.getNonLinearSolver() ) {
                case NonLinearSolverType::QNEWTON_BDIAG_JACO:
                case NonLinearSolverType::QNEWTON_BDIAG_ELAS: {
                    const auto huT = depl.at( cell_i );

                    matrix_type GT = _gradrec( msh, cl, rp, degree_infos, small_def,
                                               data.m_gradient_precomputed );

                    elem.compute_rigidity_matrix( msh, cl, bnd, rp, degree_infos, GT, huT,
                                                  this->m_time_step, behavior, small_def,
                                                  use_tangent );

                    const matrix_type stab =
                        beta_s * _stab( msh, cl, rp, degree_infos, data.m_stab_precomputed );

                    lhs = matrix_type::Zero( num_faces_dofs, num_faces_dofs );

                    int offset = 0;
                    const auto fcs = faces( msh, cl );
                    for ( size_t i = 0; i < fcs.size(); i++ ) {
                        const auto fdi = faces_infos[i];

                        if ( fdi.hasUnknowns() ) {
                            const auto facdeg = fdi.degree();
                            const auto fbs = vector_basis_size( facdeg, mesh_type::dimension - 1,
                                                                mesh_type::dimension );

                            lhs.block( offset, offset, fbs, fbs ) = elem.K_int.block(
                                num_cell_dofs + offset, num_cell_dofs + offset, fbs, fbs );
                            lhs.block( offset, offset, fbs, fbs ) += stab.block(
                                num_cell_dofs + offset, num_cell_dofs + offset, fbs, fbs );

                            offset += fbs;
                        }
                    }
                    assert( offset == num_faces_dofs );
                    break;
                }
                case NonLinearSolverType::QNEWTON_BDIAG_STAB: {
                    if ( mixed_order ) {
                        matrix_type stab =
                            beta_s * _stab( msh, cl, rp, degree_infos, data.m_stab_precomputed );

                        lhs = stab.bottomRightCorner( num_faces_dofs, num_faces_dofs );
                    } else {
                        lhs = beta_s * _mass_term( msh, cl, degree_infos );
                    }
                    break;
                }
                default: {
                    throw std::invalid_argument( "QuasiNewton option is unknown." );
                    break;
                }
                }

                const vector_type rhs = vector_type::Zero( num_faces_dofs );

                assembler.assemble( msh, cl, bnd, lhs, rhs );

                if ( bnd.cell_has_dirichlet_faces( cl ) ) {
                    ( *data.m_lhs_loc )[cell_i] = lhs;
                }
            }
            tc.toc();
            ii.m_time_rigi += elem.time_rigi;
            ii.m_time_dyna += tc.elapsed() - elem.time_rigi;

            assembler.finalize();

            tc.tic();
            const auto status = this->m_lin_solv->factorize( assembler.LHS );
            tc.toc();
            ii.m_time_solve += tc.elapsed();

            if ( status != solvers::direct_solver_status::ok ) {
                throw std::runtime_error( "Fail to factorize the matrix." );
            }
        }

        return ii;
    }

    AssemblyInfo assemble( const mesh_type &msh, const bnd_type &bnd, const param_type &rp,
                           const MeshDegreeInfo< mesh_type > &degree_infos,
                           const std::unique_ptr< func_type > &lf,
                           NonLinearData< scalar_type > &data, behavior_type &behavior,
                           StabCoeffManager< scalar_type > &stab_manager,
                           MultiTimeField< scalar_type > &fields ) override {
        elem_type elem;
        AssemblyInfo ai;

        // set RHS to zero
        this->m_assembler.initialize();
        this->m_F_int = 0.0;

        const bool small_def = ( behavior.getDeformation() == SMALL_DEF );

        const bool mixed_order = rp.m_cell_degree > rp.m_face_degree;

        // Like if it is an implicit scheme
        auto current_time = this->m_time_step.end_time();
        auto depl = fields.getCurrentField( FieldName::DEPL );
        auto depl_faces = fields.getCurrentField( FieldName::DEPL_FACES );

        const auto lhs_loc = data.m_lhs_loc;

        std::vector< vector_type > resi_cells;
        resi_cells.reserve( msh.cells_size() );

        const auto rlf = this->_getLoad( lf, current_time );

        timecounter tc, ttot;

        ttot.tic();

        for ( auto &cl : msh ) {
            const auto cell_i = msh.lookup( cl );

            const auto huT = depl.at( cell_i );

            const auto cell_infos = degree_infos.cellDegreeInfo( msh, cl );
            const auto num_cell_dofs = vector_cell_dofs( msh, cell_infos );

            const auto num_tot_dofs = huT.size();
            const auto num_faces_dofs = num_tot_dofs - num_cell_dofs;

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
            elem.compute( msh, cl, bnd, rp, degree_infos, rlf, GT, huT, this->m_time_step, behavior,
                          stab_manager, small_def, false );

            vector_type rhs = elem.RTF.tail( num_faces_dofs );
            this->m_F_int += elem.F_int.tail( num_faces_dofs ).squaredNorm();

            resi_cells.push_back( elem.RTF.head( num_cell_dofs ) );

            tc.toc();
            ai.m_time_elem += tc.elapsed();

            tc.tic();
            if ( rp.m_stab ) {
                const auto beta_s = stab_manager.getValue( msh, cl );

                matrix_type stab_F =
                    beta_s * ( _stab( msh, cl, rp, degree_infos, data.m_stab_precomputed )
                                   .bottomLeftCorner( num_faces_dofs, num_tot_dofs ) );

                rhs -= stab_F * huT;
            }
            tc.toc();
            ai.m_time_stab += tc.elapsed();

            tc.tic();
            this->m_assembler.assemble_nonlinear_rhs( msh, cl, bnd, ( *lhs_loc )[cell_i], rhs,
                                                      depl_faces );
            tc.toc();
            ai.m_time_assembler += tc.elapsed();
        }
        fields.setCurrentField( FieldName::RESI_CELLS, resi_cells );

        this->m_F_int = sqrt( this->m_F_int );

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

    SolveInfo solve() override {
        timecounter tc;

        // std::cout << "RHS" << this->m_assembler.RHS.transpose() << std::endl;

        tc.tic();
        const auto status = this->m_lin_solv->solve( this->m_assembler.RHS, this->m_system_displ );
        tc.toc();

        if ( status != solvers::direct_solver_status::ok ) {
            throw std::runtime_error( "Error during linear system resolution" );
        }

        // std::cout << "SOL" << this->m_system_displ.transpose() << std::endl;

        return SolveInfo( this->m_assembler.LHS.rows(), this->m_assembler.LHS.nonZeros(),
                          tc.elapsed() );
    }

    scalar_type postprocess( const mesh_type &msh, const bnd_type &bnd, const param_type &rp,
                             const MeshDegreeInfo< mesh_type > &degree_infos,
                             const std::unique_ptr< func_type > &lf,
                             NonLinearData< scalar_type > &data, behavior_type &behavior,
                             StabCoeffManager< scalar_type > &stab_manager,
                             MultiTimeField< scalar_type > &fields ) override {
        timecounter tc;
        tc.tic();

        auto depl_faces = fields.getCurrentField( FieldName::DEPL_FACES );

        auto [dudT, idx] = this->m_assembler.expand_solution_nonlinear(
            msh, bnd, this->m_system_displ, depl_faces );

        auto update_depl_cell = [&msh, &fields, &degree_infos, &idx,
                                 this]( const std::vector< vector_type > &depl_faces ) -> auto {
            auto depl = fields.getCurrentField( FieldName::DEPL );

            // Update cell
            for ( auto &cl : msh ) {
                const auto cell_i = msh.lookup( cl );

                const auto cell_infos = degree_infos.cellDegreeInfo( msh, cl );
                const auto faces_infos = cell_infos.facesDegreeInfo();
                const auto num_faces_dofs = vector_faces_dofs( msh, faces_infos );

                vector_type udT = vector_type( num_faces_dofs );

                const auto fcs_id = faces_id( msh, cl );
                size_t face_offset = 0;
                for ( size_t face_i = 0; face_i < fcs_id.size(); face_i++ ) {
                    const size_t face_id = fcs_id[face_i];
                    const auto n_face_dofs = idx( face_id + 1 ) - idx( face_id );

                    udT.segment( face_offset, n_face_dofs ) = depl_faces[face_id];
                    face_offset += n_face_dofs;
                }

                // Update element U^{i+1} = U^i + delta U^i
                depl.at( cell_i ).tail( num_faces_dofs ) = udT;

                // std::cout << "KT_F " << m_AL[cell_i].norm() << std::endl;
                // std::cout << "sol_F" << std::endl;
                // std::cout << xdT.transpose() << std::endl;
                // std::cout << "ft" << std::endl;
                // std::cout << m_bL[cell_i].transpose() << std::endl;
                // std::cout << "sol_T" << std::endl;
                // std::cout << xT.transpose() << std::endl;
                // std::cout << depl.at(cell_i).transpose() << std::endl;
            }
            fields.setCurrentField( FieldName::DEPL, depl );

            this->m_dyna.postprocess( msh, this->m_time_step, fields );
        };

        /* TODO: fix acceleration */

        if ( rp.getLineSearch() == LineSearchType::SECANT ) {

            auto _func = [&msh, &bnd, &rp, &degree_infos, &lf, &data, &behavior, &stab_manager,
                          &fields, &dudT, &idx, depl_faces, update_depl_cell,
                          this]( scalar_type rho, const bool compute = true ) -> auto {
                // update solution

                // Update  unknowns
                // Update face Uf^{i+1} = Uf^i + delta Uf^i

                auto depl_faces_new = depl_faces;

                for ( auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++ ) {
                    const auto fc = *itor;
                    const size_t face_id = msh.lookup( fc );

                    depl_faces_new.at( face_id ) +=
                        rho * dudT.segment( idx( face_id ), idx( face_id + 1 ) - idx( face_id ) );
                }

                fields.setCurrentField( FieldName::DEPL_FACES, depl_faces_new );

                update_depl_cell( depl_faces_new );

                // compute new residual
                if ( compute && std::abs( rho ) > 1e-32 )
                    const auto ai = this->assemble( msh, bnd, rp, degree_infos, lf, data, behavior,
                                                    stab_manager, fields );
                // TODO: Check sign
                return dudT.dot( this->m_assembler.RHS );
            };

            this->m_accel.secant( _func );
        } else {

            for ( auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++ ) {
                const auto fc = *itor;
                const size_t face_id = msh.lookup( fc );

                depl_faces.at( face_id ) +=
                    dudT.segment( idx( face_id ), idx( face_id + 1 ) - idx( face_id ) );
            }

            if ( rp.getLineSearch() != LineSearchType::NO_LS ) {

                vector_type udT_new;
                if ( rp.getLineSearch() == LineSearchType::RELAXATION ) {
                    udT_new = this->m_accel.relaxation( asVector( depl_faces ) );
                } else if ( rp.getLineSearch() == LineSearchType::AITKEN ) {
                    udT_new = this->m_accel.aitken( asVector( depl_faces ) );
                } else if ( rp.getLineSearch() == LineSearchType::ANDERSON ) {
                    udT_new = this->m_accel.anderson( asVector( depl_faces ) );
                } else if ( rp.getLineSearch() == LineSearchType::ANDERSON2 ) {
                    udT_new = this->m_accel.anderson( asVector( depl_faces ), 2 );
                } else if ( rp.getLineSearch() == LineSearchType::ANDERSON3 ) {
                    udT_new = this->m_accel.anderson( asVector( depl_faces ), 3 );
                } else if ( rp.getLineSearch() == LineSearchType::ANDERSON4 ) {
                    udT_new = this->m_accel.anderson( asVector( depl_faces ), 4 );
                } else if ( rp.getLineSearch() == LineSearchType::ANDERSON5 ) {
                    udT_new = this->m_accel.anderson( asVector( depl_faces ), 5 );
                } else if ( rp.getLineSearch() == LineSearchType::ANDERSON10 ) {
                    udT_new = this->m_accel.anderson( asVector( depl_faces ), 10 );
                } else {
                    throw std::invalid_argument( "LineSearch algorithm not supported." );
                }

                fromVector( udT_new, depl_faces );
            }

            fields.setCurrentField( FieldName::DEPL_FACES, depl_faces );

            update_depl_cell( depl_faces );
        }

        tc.toc();
        return tc.elapsed();
    }

    scalar_type post_convergence( const mesh_type &msh, const bnd_type &bnd, const param_type &rp,
                                  const MeshDegreeInfo< mesh_type > &degree_infos,
                                  const NonLinearData< scalar_type > &data,
                                  const StabCoeffManager< scalar_type > &stab_manager,
                                  MultiTimeField< scalar_type > &fields ) override {
        timecounter tc;
        tc.tic();

        if ( rp.m_stab ) {
            const auto depl = fields.getCurrentField( FieldName::DEPL );
            auto resi_cells = fields.getCurrentField( FieldName::RESI_CELLS );

            // Update cell

            for ( auto &cl : msh ) {
                const auto cell_i = msh.lookup( cl );

                const auto cell_infos = degree_infos.cellDegreeInfo( msh, cl );
                const auto num_cell_dofs = vector_cell_dofs( msh, cell_infos );

                const auto beta_s = stab_manager.getValue( msh, cl );
                const matrix_type stab =
                    beta_s * _stab( msh, cl, rp, degree_infos, data.m_stab_precomputed );

                const auto num_tot_dofs = stab.cols();

                resi_cells[cell_i] -=
                    stab.topLeftCorner( num_cell_dofs, num_tot_dofs ) * depl[cell_i];
            }
            fields.setCurrentField( FieldName::RESI_CELLS, resi_cells );

            // std::cout << "DEPL: " << norm( depl ) << std::endl;
            // std::cout << "DEPL_CELLS: " << norm( fields.getCurrentField( FieldName::DEPL_CELLS )
            // )
            //           << std::endl;
            // std::cout << "DEPL_FACES: " << norm( fields.getCurrentField( FieldName::DEPL_FACES )
            // )
            //           << std::endl;
            // std::cout << "RESI_CELLS: " << norm( resi_cells ) << std::endl;
        }
        tc.toc();
        return tc.elapsed();
    }
};
} // namespace mechanics
} // namespace disk
