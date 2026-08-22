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

// Generic non-linear iteration

#pragma once

#include "diskpp/adaptivity/adaptivity.hpp"
#include "diskpp/bases/bases.hpp"
#include "diskpp/boundary_conditions/boundary_conditions.hpp"
#include "diskpp/common/timecounter.hpp"
#include "diskpp/mechanics/NewtonSolver/Fields.hpp"
#include "diskpp/mechanics/NewtonSolver/LineSearch.hpp"
#include "diskpp/mechanics/NewtonSolver/NewtonSolverComput.hpp"
#include "diskpp/mechanics/NewtonSolver/NewtonSolverDynamic.hpp"
#include "diskpp/mechanics/NewtonSolver/NewtonSolverInformations.hpp"
#include "diskpp/mechanics/NewtonSolver/NonLinearData.hpp"
#include "diskpp/mechanics/NewtonSolver/NonLinearParameters.hpp"
#include "diskpp/mechanics/NewtonSolver/StabilizationManager.hpp"
#include "diskpp/mechanics/NewtonSolver/TimeManager.hpp"
#include "diskpp/mechanics/behaviors/laws/behaviorlaws.hpp"
#include "diskpp/methods/hho"
#include "diskpp/solvers/direct_solvers.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

namespace disk {

namespace mechanics {

/**
 * @brief Generic non-linear iteration for nonlinear solid mechanics
 *
 *  Specialized for HHO methods
 *
 *  Options :  - small and finite deformations
 *             - plasticity, hyperelasticity (various laws)
 *
 * @tparam MeshType type of the mesh
 */
template < typename MeshType >
class GenericIteration {
  protected:
    typedef MeshType mesh_type;
    typedef typename mesh_type::cell cell_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type point_type;

    typedef dynamic_matrix< scalar_type > matrix_type;
    typedef dynamic_vector< scalar_type > vector_type;

    typedef NonLinearParameters< scalar_type > param_type;
    typedef vector_boundary_conditions< mesh_type > bnd_type;
    typedef Behavior< mesh_type > behavior_type;

    typedef vector_mechanics_hho_assembler< mesh_type > assembler_type;
    typedef mechanical_computation< mesh_type > elem_type;
    typedef dynamic_computation< mesh_type > dyna_type;

    typedef std::function< static_vector< scalar_type, mesh_type::dimension >( const point_type &,
                                                                               scalar_type ) >
        func_type;

    using load_at_time_type =
        std::function< static_vector< scalar_type, mesh_type::dimension >( const point_type & ) >;

    vector_type m_system_displ;

    assembler_type m_assembler;

    std::vector< vector_type > m_bL;
    std::vector< matrix_type > m_AL;

    TimeStep< scalar_type > m_time_step;

    std::shared_ptr< solvers::sparse_solver< scalar_type > > m_lin_solv;

    dyna_type m_dyna;

    scalar_type m_F_int, m_resi_init;

    ConvergenceAcceleration< scalar_type > m_accel;

    bool m_verbose;

    std::unique_ptr< load_at_time_type > _getLoad( const std::unique_ptr< func_type > &lf,
                                                   scalar_type time ) {
        if ( lf ) {
            const auto load_function = *lf;
            const auto time_eval = time;
            return std::make_unique< std::function<
                static_vector< scalar_type, mesh_type::dimension >( const point_type & ) > >(
                [load_function, time_eval]( const point_type &p ) -> auto {
                    return load_function( p, time_eval );
                } );
        }

        return nullptr;
    }

  public:
    GenericIteration() : m_verbose( false ) {};

    GenericIteration( const mesh_type &msh, const bnd_type &bnd, const param_type &rp,
                      const MeshDegreeInfo< mesh_type > &degree_infos,
                      const std::shared_ptr< solvers::sparse_solver< scalar_type > > lin_solv,
                      const TimeStep< scalar_type > &current_step )
        : m_verbose( rp.m_verbose ),
          m_time_step( current_step ),
          m_dyna( rp ),
          m_lin_solv( lin_solv ) {
        m_AL.clear();
        m_AL.resize( msh.cells_size() );

        m_bL.clear();
        m_bL.resize( msh.cells_size() );

        m_assembler = assembler_type( msh, degree_infos, bnd );
    }

    virtual ~GenericIteration() = default;

    bool verbose( void ) const { return m_verbose; }

    void verbose( bool v ) { m_verbose = v; }

    virtual InitInfo initialize( const mesh_type &msh, const bnd_type &bnd, const param_type &rp,
                                 const MeshDegreeInfo< mesh_type > &degree_infos,
                                 NonLinearData< scalar_type > &data, behavior_type &behavior,
                                 const StabCoeffManager< scalar_type > &stab_manager,
                                 MultiTimeField< scalar_type > &fields ) {
        timecounter tc;
        tc.tic();
        m_dyna.prediction( msh, degree_infos, m_time_step, fields );
        tc.toc();

        InitInfo ii;
        ii.m_time_dyna = tc.elapsed();
        return ii;
    }

    virtual AssemblyInfo assemble( const mesh_type &msh, const bnd_type &bnd, const param_type &rp,
                                   const MeshDegreeInfo< mesh_type > &degree_infos,
                                   const std::unique_ptr< func_type > &lf,
                                   NonLinearData< scalar_type > &data, behavior_type &behavior,
                                   StabCoeffManager< scalar_type > &stab_manager,
                                   MultiTimeField< scalar_type > &fields ) {

        throw std::runtime_error( "GenericIteration.assemble has to be overloaded." );
        return AssemblyInfo();
    }

    virtual SolveInfo solve() {
        timecounter tc;

        // std::cout << "LHS" << m_assembler.LHS << std::endl;
        // std::cout << "RHS" << m_assembler.RHS << std::endl;

        tc.tic();
        const auto status = m_lin_solv->solve( m_assembler.LHS, m_assembler.RHS, m_system_displ );
        tc.toc();

        if ( status != solvers::direct_solver_status::ok ) {
            throw std::runtime_error( "Error during linear system resolution" );
        }

        return SolveInfo( m_assembler.LHS.rows(), m_assembler.LHS.nonZeros(), tc.elapsed() );
    }

    virtual scalar_type postprocess( const mesh_type &msh, const bnd_type &bnd,
                                     const param_type &rp,
                                     const MeshDegreeInfo< mesh_type > &degree_infos,
                                     const std::unique_ptr< func_type > &lf,
                                     NonLinearData< scalar_type > &data, behavior_type &behavior,
                                     StabCoeffManager< scalar_type > &stab_manager,
                                     MultiTimeField< scalar_type > &fields ) {
        throw std::runtime_error( "GenericIteration.postprocess has to be overloaded." );
        return 0.0;
    }

    bool convergence( const param_type &rp, const size_t iter,
                      const MultiTimeField< scalar_type > &fields ) {
        // norm of the solution
        auto norm_sol = norm( fields.getCurrentField( FieldName::DEPL_FACES ) );

        // norm of the rhs
        const scalar_type residual = m_assembler.RHS.norm();
        scalar_type max_error = 0.0;
        for ( size_t i = 0; i < m_assembler.RHS.size(); i++ ) {
            max_error = std::max( max_error, std::abs( m_assembler.RHS( i ) ) );
        }

        // norm of the increment
        const scalar_type error_incr = m_system_displ.norm();
        scalar_type relative_displ = 1.0, relative_error = 1.0;

        if ( m_F_int > 1e-12 ) {
            relative_error = residual / m_F_int;
        }
        if ( norm_sol > 1e-12 && iter > 0 ) {
            relative_displ = error_incr / norm_sol;
        }

        if ( m_verbose and rp.getNonLinearSolver() != NonLinearSolverType::EXPLICIT ) {
            std::string s_iter = "   " + std::to_string( iter ) + "               ";
            s_iter.resize( 9 );

            if ( iter == 0 ) {
                std::cout
                    << "----------------------------------------------------------------------"
                       "------------------------"
                    << std::endl;
                std::cout << "| Iteration | Norme l2 incr | Relative incr |  Residual l2  | "
                             "Relative error | Maximum error |"
                          << std::endl;
                std::cout
                    << "----------------------------------------------------------------------"
                       "------------------------"
                    << std::endl;
            }
            std::ios::fmtflags f( std::cout.flags() );
            std::cout.precision( 5 );
            std::cout.setf( std::iostream::scientific, std::iostream::floatfield );
            std::cout << "| " << s_iter << " |   " << error_incr << " |   " << relative_displ
                      << " |   " << residual << " |   " << relative_error << "  |  " << max_error
                      << "  |" << std::endl;
            std::cout << "-------------------------------------------------------------------------"
                         "---------------------"
                      << std::endl;
            std::cout.flags( f );
        }

        const scalar_type error = std::max( relative_displ, relative_error );
        // const scalar_type error = relative_error;

        if ( !std::isfinite( error ) || std::isnan( error ) || !std::isfinite( max_error ) ||
             std::isnan( max_error ) || !std::isfinite( residual ) || std::isnan( residual ) ||
             !std::isfinite( -residual ) || std::isnan( -residual ) || residual > 1e100 ) {
            throw std::runtime_error( "Norm of residual is not finite." );
        }

        if ( residual > 1e20 || relative_displ > 1e20 || relative_error > 1e20 ) {
            throw std::runtime_error( "Norm of residual is too large." );
        }

        if ( iter == 0 ) {
            this->m_resi_init = residual;
        }

        if ( residual > 1e10 * this->m_resi_init ) {
            throw std::runtime_error( "Norm of residual diverges." );
        }

        if ( error <= rp.getConvergenceCriteria() ) {
            return true;
        } else {
            return false;
        }
    }

    virtual scalar_type post_convergence( const mesh_type &msh, const bnd_type &bnd,
                                          const param_type &rp,
                                          const MeshDegreeInfo< mesh_type > &degree_infos,
                                          const NonLinearData< scalar_type > &data,
                                          const StabCoeffManager< scalar_type > &stab_manager,
                                          MultiTimeField< scalar_type > &fields ) {
        return 0.0;
    }
};
} // namespace mechanics
} // namespace disk