/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2018, 2024               nicolas.pignet@enpc.fr
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

#include "diskpp/solvers/direct_solvers.hpp"

#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>

namespace disk {

namespace mechanics {

enum StabilizationType : int { HDG = 0, HHO = 1, HHO_SYM = 4, NO = 2, DG = 3 };

enum FrictionType : int {
    NO_FRICTION = 0,
    TRESCA = 1,
    COULOMB = 2,
};

enum DynamicType : int {
    STATIC = 0,
    NEWMARK = 1,
    BACKWARD_EULER = 2,
    THETA = 3,
    CRANK_NICOLSON = 4,
    LEAP_FROG = 5,
};

enum NonLinearSolverType {
    NEWTON,
    PICARD,
    QNEWTON_BDIAG_JACO,
    QNEWTON_BDIAG_STAB,
    QNEWTON_BDIAG_ELAS,
    EXPLICIT,
};

enum LineSearchType {
    NO_LS,
    RELAXATION,
    AITKEN,
    SECANT,
    ANDERSON,
    ANDERSON2,
    ANDERSON3,
    ANDERSON4,
    ANDERSON5,
    ANDERSON10,
};

std::string StabilizationName( const StabilizationType &type ) {
    switch ( type ) {
    case StabilizationType::HDG: {
        return "HDG";
        break;
    }
    case StabilizationType::HHO: {
        return "HHO";
        break;
    }
    case StabilizationType::HHO_SYM: {
        return "HHO_SYM";
        break;
    }
    case StabilizationType::NO: {
        return "NO";
        break;
    }
    case StabilizationType::DG: {
        return "DG";
        break;
    }
    default:
        break;
    }

    throw std::invalid_argument( "Case not supported" );
}

std::string FrictionName( const FrictionType &type ) {
    switch ( type ) {
    case FrictionType::NO_FRICTION: {
        return "NO_FRICTION";
        break;
    }
    case FrictionType::TRESCA: {
        return "TRESCA";
        break;
    }
    case FrictionType::COULOMB: {
        return "COULOMB";
        break;
    }
    default:
        break;
    }

    throw std::invalid_argument( "Case not supported" );
}

std::string DynaSchemeName( const DynamicType &type ) {
    switch ( type ) {
    case DynamicType::STATIC: {
        return "STATIC";
        break;
    }
    case DynamicType::NEWMARK: {
        return "NEWMARK";
        break;
    }
    case DynamicType::BACKWARD_EULER: {
        return "BACKWARD_EULER";
        break;
    }
    case DynamicType::THETA: {
        return "THETA";
        break;
    }
    case DynamicType::CRANK_NICOLSON: {
        return "CRANK_NICOLSON";
        break;
    }
    case DynamicType::LEAP_FROG: {
        return "LEAP_FROG";
        break;
    }
    default:
        break;
    }

    throw std::invalid_argument( "Case not supported" );
}

std::string LinearSolverName( const solvers::direct_solver &type ) {
    switch ( type ) {
    case solvers::direct_solver::autosel: {
        return "AUTO";
        break;
    }
#ifdef HAVE_PARDISO
    case solvers::direct_solver::pardiso: {
        return "PARDISO_LU";
        break;
    }
#endif
#ifdef HAVE_MUMPS
    case solvers::direct_solver::mumps: {
        return "MUMPS_LU";
        break;
    }
#endif
    case solvers::direct_solver::sparselu: {
        return "SPARSE_LU";
        break;
    }
    default:
        break;
    }

    throw std::invalid_argument( "Case not supported" );
}

std::string NonLinearSolverName( const NonLinearSolverType &type ) {
    switch ( type ) {
    case NonLinearSolverType::NEWTON: {
        return "NEWTON";
        break;
    }
    case NonLinearSolverType::PICARD: {
        return "PICARD";
        break;
    }
    case NonLinearSolverType::QNEWTON_BDIAG_JACO: {
        return "QNEWTON_BDIAG_JACO";
        break;
    }
    case NonLinearSolverType::QNEWTON_BDIAG_STAB: {
        return "QNEWTON_BDIAG_STAB";
        break;
    }
    case NonLinearSolverType::QNEWTON_BDIAG_ELAS: {
        return "QNEWTON_BDIAG_ELAS";
        break;
    }
    case NonLinearSolverType::EXPLICIT: {
        return "EXPLICIT";
        break;
    }
    default:
        break;
    }

    throw std::invalid_argument( "Case not supported" );
}

std::string LineSearchName( const LineSearchType &type ) {
    switch ( type ) {
    case LineSearchType::NO_LS: {
        return "NO";
        break;
    }
    case LineSearchType::RELAXATION: {
        return "RELAXATION";
        break;
    }
    case LineSearchType::AITKEN: {
        return "AITKEN";
        break;
    }
    case LineSearchType::SECANT: {
        return "SECANT";
        break;
    }
    case LineSearchType::ANDERSON: {
        return "ANDERSON";
        break;
    }
    case LineSearchType::ANDERSON2: {
        return "ANDERSON2";
        break;
    }
    case LineSearchType::ANDERSON3: {
        return "ANDERSON3";
        break;
    }
    case LineSearchType::ANDERSON4: {
        return "ANDERSON4";
        break;
    }
    case LineSearchType::ANDERSON5: {
        return "ANDERSON5";
        break;
    }
    case LineSearchType::ANDERSON10: {
        return "ANDERSON10";
        break;
    }
    default:
        break;
    }

    throw std::invalid_argument( "Case not supported" );
}

std::string BoolName( bool value ) { return value ? "TRUE" : "FALSE"; }

template < typename T >
class NonLinearParameters {
  public:
    int m_face_degree; // face degree
    int m_cell_degree; // cell_degree
    int m_grad_degree; // grad degree

    std::vector< std::pair< T, int > > m_time_step; // number of time time_step
    bool m_has_user_end_time;                       // final time is given
    T m_user_end_time;                              // final time of the simulation
    int m_sublevel;                                 // number of sublevel if there are problems
    int m_iter_max;                                 // maximun nexton iteration
    T m_epsilon;                                    // stop criteria

    bool m_verbose; // some printing

    bool m_precomputation; // to compute the gradient before (it's memory consuption)

    StabilizationType m_stab_type; // type of stabilization
    T m_beta;                      // stabilization parameter
    bool m_stab;                   // stabilization yes or no
    bool m_adapt_stab;             // adaptative stabilization

    DynamicType m_dyna_type;                // type of dyna
    std::map< std::string, T > m_dyna_para; // list of parameters
    T m_cfl_factor;                         // CFL factor

    int m_n_time_save;          // number of saving
    std::list< T > m_time_save; // list of time where we save result;

    T m_theta;                // theta-parameter for contact
    T m_gamma_0;              // parameter for Nitsche
    T m_threshold;            // threshol for Tesca friction
    FrictionType m_frot_type; // Friction type ?

    solvers::direct_solver m_lin_solv;    // linear solver
    NonLinearSolverType m_nlin_solv;      // non-linear solver
    LineSearchType m_lsearch;             // line-search

    NonLinearParameters()
        : m_face_degree( 1 ),
          m_cell_degree( 1 ),
          m_grad_degree( 1 ),
          m_sublevel( 5 ),
          m_iter_max( 20 ),
          m_epsilon( T( 1E-6 ) ),
          m_verbose( false ),
          m_precomputation( false ),
          m_stab( true ),
          m_beta( 1 ),
          m_stab_type( StabilizationType::HHO ),
          m_n_time_save( 0 ),
          m_user_end_time( 1.0 ),
          m_has_user_end_time( false ),
          m_adapt_stab( false ),
          m_theta( 1 ),
          m_gamma_0( 1 ),
          m_threshold( 0 ),
          m_frot_type( FrictionType::NO_FRICTION ),
          m_dyna_type( DynamicType::STATIC ),
          m_lin_solv( solvers::direct_solver::autosel ),
          m_nlin_solv( NonLinearSolverType::NEWTON ),
          m_lsearch( LineSearchType::NO_LS ),
          m_cfl_factor( 0.99 ) {
        m_time_step.push_back( std::make_pair( m_user_end_time, 1 ) );
    }

    void
    error_keyword(int line, std::string keyword, std::string value)
    {
        throw std::runtime_error("Error during parsing in line: " + std::to_string(line) + "." +
                                 "  Keyword: " + keyword + " has an unexpeced value: " + value);
    }

    void infos() {
        std::cout << "Nonlinear Solver's parameters:" << std::endl;
        std::cout << " - Face degree: " << m_face_degree << std::endl;
        std::cout << " - Cell degree: " << m_cell_degree << std::endl;
        std::cout << " - Grad degree: " << m_grad_degree << std::endl;
        std::cout << " - Stabilization ?: " << BoolName( !( m_stab_type == StabilizationType::NO ) )
                  << std::endl;
        std::cout << " - AdaptativeStabilization ?: " << BoolName( m_adapt_stab ) << std::endl;
        std::cout << " - Type: " << StabilizationName( m_stab_type ) << std::endl;
        std::cout << " - Beta: " << m_beta << std::endl;
        std::cout << " - Verbose: " << BoolName( m_verbose ) << std::endl;
        std::cout << " - Sublevel: " << m_sublevel << std::endl;
        std::cout << " - IterMax: " << m_iter_max << std::endl;
        std::cout << " - Epsilon: " << m_epsilon << std::endl;
        std::cout << " - LinearSolver: " << LinearSolverName( m_lin_solv ) << std::endl;
        std::cout << " - NonLinearSolver: " << NonLinearSolverName( m_nlin_solv ) << std::endl;
        std::cout << " - LineSearch: " << LineSearchName( m_lsearch ) << std::endl;
        std::cout << " - Precomputation: " << BoolName( m_precomputation ) << std::endl;
        std::cout << " - Dynamic scheme: " << DynaSchemeName( m_dyna_type ) << std::endl;
        std::cout << " - CFL factor: " << m_cfl_factor << std::endl;
        std::cout << " - Friction ?: " << FrictionName( m_frot_type ) << std::endl;
        std::cout << " - Threshold: " << m_threshold << std::endl;
        std::cout << " - Gamma_0: " << m_gamma_0 << std::endl;
        std::cout << " - Theta: " << m_theta << std::endl;
    }

    bool readParameters( const std::string &filename ) {
        std::ifstream ifs( filename );
        std::string keyword;
        int line( 0 );

        if ( !ifs.is_open() ) {
            std::cout << "Error opening " << filename << std::endl;
            return false;
        }

        ifs >> keyword;
        line++;
        if ( keyword != "BeginParameters" ) {
            std::cout << "Expected keyword \"BeginParameters\" line: " << line << std::endl;
            return false;
        }

        ifs >> keyword;
        line++;
        while ( keyword != "EndParameters" ) {
            // std::cout << "Keyword: " << keyword << std::endl;
            if ( keyword == "FaceDegree" ) {
                ifs >> m_face_degree;
                line++;
            } else if ( keyword == "CellDegree" ) {
                ifs >> m_cell_degree;
                line++;
            } else if ( keyword == "GradDegree" ) {
                ifs >> m_grad_degree;
                line++;
            } else if ( keyword == "Sublevel" ) {
                ifs >> m_sublevel;
                line++;
            } else if ( keyword == "TimeStep" ) {
                int n_time_step( 0 );
                ifs >> n_time_step;
                line++;

                m_time_step.clear();
                m_time_step.reserve( n_time_step );
                for ( int i = 0; i < n_time_step; i++ ) {
                    T time( 0.0 );
                    int time_step( 0 );
                    ifs >> time >> time_step;
                    m_time_step.push_back( std::make_pair( time, time_step ) );
                    line++;
                }
            } else if ( keyword == "FinalTime" ) {
                ifs >> m_user_end_time;
                line++;

                m_has_user_end_time = true;
            } else if ( keyword == "TimeSave" ) {
                ifs >> m_n_time_save;
                line++;

                m_time_save.clear();
                for ( int i = 0; i < m_n_time_save; i++ ) {
                    T time( 0.0 );
                    ifs >> time;
                    m_time_save.push_back( time );
                    line++;
                }
            } else if ( keyword == "AdaptativeStabilization" ) {
                std::string logical;
                ifs >> logical;
                line++;
                m_adapt_stab = false;
                if ( logical == "true" || logical == "True" )
                    m_adapt_stab = true;
            } else if ( keyword == "StabType" ) {
                std::string type;
                ifs >> type;
                line++;
                m_stab = true;
                if ( type == "HDG" )
                    m_stab_type = StabilizationType::HDG;
                else if ( type == "HHO" )
                    m_stab_type = StabilizationType::HHO;
                else if ( type == "HHO_SYM" )
                    m_stab_type = StabilizationType::HHO_SYM;
                else if ( type == "DG" )
                    m_stab_type = StabilizationType::DG;
                else if ( type == "NO" ) {
                    m_stab = false;
                    m_stab_type = StabilizationType::NO;
                } else
                    error_keyword(line, keyword, type);
            } else if ( keyword == "Beta" ) {
                ifs >> m_beta;
                line++;
            } else if ( keyword == "Verbose" ) {
                std::string logical;
                ifs >> logical;
                line++;
                m_verbose = false;
                if ( logical == "true" || logical == "True" )
                    m_verbose = true;
            } else if ( keyword == "IterMax" ) {
                ifs >> m_iter_max;
                line++;
            } else if ( keyword == "Epsilon" ) {
                ifs >> m_epsilon;
                line++;
            } else if ( keyword == "Precomputation" ) {
                std::string logical;
                ifs >> logical;
                line++;
                m_precomputation = false;
                if ( logical == "true" || logical == "True" )
                    m_precomputation = true;
            } else if ( keyword == "Theta" ) {
                ifs >> m_theta;
                line++;
            } else if ( keyword == "Gamma0" ) {
                ifs >> m_gamma_0;
                line++;
            } else if ( keyword == "Friction" ) {
                std::string type;
                ifs >> type;
                line++;
                if ( type == "NO" )
                    m_frot_type = FrictionType::NO_FRICTION;
                else if ( type == "TRESCA" )
                    m_frot_type = FrictionType::TRESCA;
                else if ( type == "COULOMB" )
                    m_frot_type = FrictionType::COULOMB;
                else
                    error_keyword(line, keyword, type);
            } else if ( keyword == "Threshold" ) {
                ifs >> m_threshold;
            } else if ( keyword == "Dynamic" ) {
                std::string type;
                ifs >> type;
                line++;
                if ( type == "STATIC" || type == "NO" )
                    m_dyna_type = DynamicType::STATIC;
                else if ( type == "NEWMARK" )
                    m_dyna_type = DynamicType::NEWMARK;
                else if ( type == "BACKWARD_EULER" )
                    m_dyna_type = DynamicType::BACKWARD_EULER;
                else if (type == "THETA")
                    m_dyna_type = DynamicType::THETA;
                else if ( type == "CRANK_NICOLSON" )
                    m_dyna_type = DynamicType::CRANK_NICOLSON;
                else if ( type == "LEAP_FROG" )
                    m_dyna_type = DynamicType::LEAP_FROG;
                else
                    error_keyword(line, keyword, type);
            } else if ( keyword == "CFL" ) {
                ifs >> m_cfl_factor;
                line++;
            } else if ( keyword == "NLSolver" ) {
                std::string type;
                ifs >> type;
                line++;
                if ( type == "NEWTON" ) {
                    m_nlin_solv = NonLinearSolverType::NEWTON;
                } else if ( type == "PICARD" ) {
                    m_nlin_solv = NonLinearSolverType::PICARD;
                } else if ( type == "QNEWTON_BDIAG_JACO" ) {
                    m_nlin_solv = NonLinearSolverType::QNEWTON_BDIAG_JACO;
                } else if ( type == "QNEWTON_BDIAG_STAB" ) {
                    m_nlin_solv = NonLinearSolverType::QNEWTON_BDIAG_STAB;
                } else if ( type == "QNEWTON_BDIAG_ELAS" ) {
                    m_nlin_solv = NonLinearSolverType::QNEWTON_BDIAG_ELAS;
                } else if ( type == "EXPLICIT" ) {
                    m_nlin_solv = NonLinearSolverType::EXPLICIT;
                } else {
                    error_keyword(line, keyword, type);
                }
            } else if ( keyword == "LineSearch" ) {
                std::string type;
                ifs >> type;
                line++;
                if ( type == "NO" || type == "NO_LS" ) {
                    m_lsearch = LineSearchType::NO_LS;
                } else if ( type == "RELAXATION" ) {
                    m_lsearch = LineSearchType::RELAXATION;
                } else if ( type == "AITKEN" ) {
                    m_lsearch = LineSearchType::AITKEN;
                } else if ( type == "SECANT" ) {
                    m_lsearch = LineSearchType::SECANT;
                } else if ( type == "ANDERSON" ) {
                    m_lsearch = LineSearchType::ANDERSON;
                } else if ( type == "ANDERSON2" ) {
                    m_lsearch = LineSearchType::ANDERSON2;
                } else if ( type == "ANDERSON3" ) {
                    m_lsearch = LineSearchType::ANDERSON3;
                } else if ( type == "ANDERSON4" ) {
                    m_lsearch = LineSearchType::ANDERSON4;
                } else if ( type == "ANDERSON5" ) {
                    m_lsearch = LineSearchType::ANDERSON5;
                } else if ( type == "ANDERSON10" ) {
                    m_lsearch = LineSearchType::ANDERSON10;
                } else {
                    error_keyword(line, keyword, type);
                }
            } else {
                error_keyword(line, keyword, "");
            }

            ifs >> keyword;
            line++;
        }

        ifs.close();
        return true;
    }

    void setFaceDegree( const int face_degree ) { m_face_degree = face_degree; }

    int getFaceDegree() const { return m_face_degree; }

    void setCellDegree( const int cell_degree ) { m_cell_degree = cell_degree; }

    int getCellDegree() const { return m_cell_degree; }

    void setGradDegree( const int grad_degree ) { m_grad_degree = grad_degree; }

    int getGradDegree() const { return m_face_degree; }

    void setStabilizationParameter( const T stab_para ) { m_beta = stab_para; }

    T getStabilizationParameter() const { return m_beta; }

    void setVerbose( const bool verbose ) { m_verbose = verbose; }

    bool getVerbose() const { return m_verbose; }

    void setPrecomputation( const bool precomp ) { m_precomputation = precomp; }

    bool getPrecomputation() const { return m_precomputation; }

    bool isUnsteady() const { return m_dyna_type != STATIC; }

    DynamicType getUnsteadyScheme() const { return m_dyna_type; }

    void setUnsteadyScheme( const DynamicType &scheme ) { m_dyna_type = scheme; }

    auto getUnsteadyParameters() const { return m_dyna_para; }

    auto getCFLFactor() const { return m_cfl_factor; }

    void setUnsteadyParameters( const std::map< std::string, T > dyna_para ) {
        m_dyna_para = dyna_para;
    }

    void setTimeStep( const T end_time, const int n_time_step ) {
        m_time_step.clear();
        m_time_step.push_back( std::make_pair( end_time, n_time_step ) );
    }

    void setLinearSolver( const solvers::direct_solver &type ) { m_lin_solv = type; }
    solvers::direct_solver getLinearSolver() const { return m_lin_solv; }

    void setNonLinearSolver( const NonLinearSolverType &type ) { m_nlin_solv = type; }
    NonLinearSolverType getNonLinearSolver() const { return m_nlin_solv; }

    LineSearchType getLineSearch() const { return m_lsearch; }
    void setLineSearch( const LineSearchType &type ) { m_lsearch = type; }

    void setMaximumNumberNLIteration( const int &n_iter ) { m_iter_max = n_iter; }
    int getMaximumNumberNLIteration() const { return m_iter_max; }

    void setConvergenceCriteria( const T &eps ) { m_epsilon = eps; }
    T getConvergenceCriteria() const { return m_epsilon; }
};

} // namespace mechanics
} // namespace disk