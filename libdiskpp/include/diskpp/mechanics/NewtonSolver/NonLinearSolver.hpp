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

// NewtonRaphson_solver

#pragma once

#include "diskpp/adaptivity/adaptivity.hpp"
#include "diskpp/boundary_conditions/boundary_conditions.hpp"
#include "diskpp/mechanics/NewtonSolver/NewtonSolverInformations.hpp"
#include "diskpp/mechanics/NewtonSolver/NonLinearStep.hpp"
#include "diskpp/mechanics/behaviors/tensor_conversion.hpp"
#include "diskpp/mechanics/stress_tensors.hpp"
#include "diskpp/methods/hho"
#include "diskpp/output/gmshConvertMesh.hpp"
#include "diskpp/output/gmshDisk.hpp"
#include "diskpp/output/plotOverTime.hpp"
#include "diskpp/output/postMesh.hpp"

#include <iostream>
#include <sstream>
#include <vector>

#ifdef HAVE_MGIS
#include "MGIS/Behaviour/Behaviour.hxx"
#endif

#include "diskpp/common/timecounter.hpp"

namespace disk {

namespace mechanics {

/**
 * @brief Newton-Raphson solver for nonlinear solid mechanics
 *
 *  Specialized for HHO methods
 *
 *  Options :  - small and finite deformations
 *             - plasticity, hyperelasticity (various laws)
 *
 * @tparam Mesh type of the mesh
 */
template < typename Mesh >
class NonLinearSolver {
    typedef Mesh mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type point_type;

    typedef dynamic_matrix< scalar_type > matrix_type;
    typedef dynamic_vector< scalar_type > vector_type;

    typedef NonLinearParameters< scalar_type > param_type;
    typedef vector_boundary_conditions< mesh_type > bnd_type;
    typedef Behavior< mesh_type > behavior_type;
    typedef PlotPointOverTime< mesh_type > ppt_type;

    typedef std::function< static_vector< scalar_type, mesh_type::dimension >(
        const point< scalar_type, mesh_type::dimension > &, scalar_type ) >
        func_type;

    bnd_type m_bnd;
    const mesh_type &m_msh;
    param_type m_rp;
    behavior_type m_behavior;
    MeshDegreeInfo< mesh_type > m_degree_infos;
    StabCoeffManager< scalar_type > m_stab_manager;
    PostMesh< mesh_type > m_post_mesh;
    MultiTimeField< scalar_type > m_fields;
    NonLinearData< scalar_type > m_data;

    std::shared_ptr< solvers::sparse_solver< scalar_type > > m_lin_solv;

    std::vector< ppt_type > m_ppt;

    std::unique_ptr< func_type > m_load;

    bool m_verbose, m_convergence;

    void init_degree( const size_t cell_degree, const size_t face_degree,
                      const size_t grad_degree ) {
        m_degree_infos =
            MeshDegreeInfo< mesh_type >( m_msh, cell_degree, face_degree, grad_degree );

        if ( m_bnd.nb_faces_contact() > 0 ) {
            for ( auto itor = m_msh.boundary_faces_begin(); itor != m_msh.boundary_faces_end();
                  itor++ ) {
                const auto bfc = *itor;
                const auto face_id = m_msh.lookup( bfc );

                if ( m_bnd.contact_boundary_type( face_id ) == SIGNORINI_FACE ) {
                    m_degree_infos.degree( m_msh, bfc, face_degree + 1 );
                    break;
                }
            }
        }
    }

    // Initializa data structures
    void init( void ) {
        TimeField< scalar_type > tf;
        tf.createZeroField( FieldName::DEPL, m_msh, m_degree_infos );
        tf.createZeroField( FieldName::DEPL_CELLS, m_msh, m_degree_infos );
        tf.createZeroField( FieldName::DEPL_FACES, m_msh, m_degree_infos );
        if ( m_rp.isUnsteady() ) {
            tf.createZeroField( FieldName::VITE, m_msh, m_degree_infos );
            tf.createZeroField( FieldName::VITE_CELLS, m_msh, m_degree_infos );
            tf.createZeroField( FieldName::ACCE, m_msh, m_degree_infos );
            tf.createZeroField( FieldName::ACCE_CELLS, m_msh, m_degree_infos );
        }
        m_fields.setCurrentTimeField( tf );

        // compute mesh for post-processing
        m_post_mesh = PostMesh< mesh_type >( m_msh );

        if ( m_verbose ) {
            std::cout << "** Numbers of cells: " << m_msh.cells_size() << std::endl;
            std::cout << "** Numbers of faces: " << m_msh.faces_size()
                      << "  ( boundary faces: " << m_msh.boundary_faces_size() << " )" << std::endl;
            std::cout << "** Numbers of dofs after static condensation: " << this->numberOfDofs()
                      << std::endl;
            std::cout << " " << std::endl;
        }
    }

    /**
     * @brief Precompute the gradient reconstruction and stabilization operator for HHO methods
     *
     */
    void pre_computation( void ) {
        m_data.m_gradient_precomputed.clear();
        m_data.m_gradient_precomputed.reserve( m_msh.cells_size() );

        m_data.m_stab_precomputed.clear();
        m_data.m_stab_precomputed.reserve( m_msh.cells_size() );

        for ( auto &cl : m_msh ) {
            // std::cout << m_degree_infos.cellDegreeInfo(m_msh, cl) << std::endl;
            /////// Gradient Reconstruction /////////
            if ( m_behavior.getDeformation() == SMALL_DEF ) {
                const auto sgr = make_matrix_hho_symmetric_gradrec( m_msh, cl, m_degree_infos );
                m_data.m_gradient_precomputed.push_back( sgr.first );
            } else {
                const auto gr = make_matrix_hho_gradrec( m_msh, cl, m_degree_infos );
                m_data.m_gradient_precomputed.push_back( gr.first );
            }

            if ( m_rp.m_stab ) {
                switch ( m_rp.m_stab_type ) {
                case StabilizationType::HHO: {
                    const auto recons_scalar =
                        make_scalar_hho_laplacian( m_msh, cl, m_degree_infos );
                    m_data.m_stab_precomputed.push_back( make_vector_hho_stabilization_optim(
                        m_msh, cl, recons_scalar.first, m_degree_infos ) );
                    break;
                }
                case StabilizationType::HHO_SYM: {
                    const auto recons =
                        make_vector_hho_symmetric_laplacian( m_msh, cl, m_degree_infos );
                    m_data.m_stab_precomputed.push_back(
                        make_vector_hho_stabilization( m_msh, cl, recons.first, m_degree_infos ) );
                    break;
                }
                case StabilizationType::HDG: {
                    m_data.m_stab_precomputed.push_back(
                        make_vector_hdg_stabilization( m_msh, cl, m_degree_infos ) );
                    break;
                }
                case StabilizationType::DG: {
                    m_data.m_stab_precomputed.push_back(
                        make_vector_dg_stabilization( m_msh, cl, m_degree_infos ) );
                    break;
                }
                case StabilizationType::NO: {
                    break;
                }
                default:
                    throw std::invalid_argument( "Unknown stabilization" );
                }
            }
        }
    }

    // compute l2 error
    auto _eval( FieldName name, const int cell_id, const point_type &pt ) {
        const auto field = m_fields.getCurrentField( name );

        const auto cl = m_msh[cell_id];

        const auto di = m_degree_infos.cellDegreeInfo( m_msh, cl );
        const auto cb = make_vector_monomial_basis( m_msh, cl, di.cell_degree() );
        const vector_type x = field.at( cell_id );

        const auto phi = cb.eval_functions( pt );

        return eval( x, phi );
    }

  public:
    NonLinearSolver( const mesh_type &msh, const bnd_type &bnd, const param_type &rp )
        : m_msh( msh ),
          m_verbose( rp.m_verbose ),
          m_convergence( false ),
          m_rp( rp ),
          m_bnd( bnd ),
          m_stab_manager( msh, rp.m_beta ),
          m_fields( getNumberOfStepToSave( rp ) ),
          m_load( nullptr ),
          m_lin_solv(
              std::make_shared< solvers::sparse_solver< scalar_type > >( rp.getLinearSolver() ) ) {
        if ( m_verbose ) {
            std::cout << "------------------------------------------------------------------------"
                         "-------------"
                      << std::endl;
            std::cout << "|********************** Nonlinear Solver for solid mechanics "
                         "***********************|"
                      << std::endl;
            std::cout << "------------------------------------------------------------------------"
                         "-------------"
                      << std::endl;
        }
        int face_degree = rp.m_face_degree;
        if ( rp.m_face_degree < 0 ) {
            std::cout << "'face_degree' should be > 0. Reverting to 1." << std::endl;
            face_degree = 1;
        }

        m_rp.m_face_degree = face_degree;

        int cell_degree = rp.m_cell_degree;
        if ( ( face_degree - 1 > cell_degree ) or ( cell_degree > face_degree + 1 ) ) {
            std::cout << "'cell_degree' should be 'face_degree + 1' =>"
                      << "'cell_degree' => 'face_degree -1'. Reverting to 'face_degree'."
                      << std::endl;
            cell_degree = face_degree;
        }

        m_rp.m_cell_degree = cell_degree;

        int grad_degree = rp.m_grad_degree;
        if ( grad_degree < face_degree ) {
            std::cout << "'grad_degree' should be >= 'face_degree'. Reverting to 'face_degree'."
                      << std::endl;
            grad_degree = face_degree;
        }

        if ( m_verbose ) {
            m_rp.infos();
        }

        // Initialization
        if ( m_verbose ) {
            std::cout << std::endl;
            std::cout << "Initialization ..." << std::endl;
        }
        this->init_degree( cell_degree, face_degree, grad_degree );
        this->init();
    }

    /**
     * @brief return a boolean to know if the verbosity mode is activated
     *
     */
    bool verbose( void ) const { return m_verbose; }

    /**
     * @brief Set the verbosity mode
     *
     * @param v boolean to activate or desactivate the verbosity mode
     */
    void verbose( bool v ) { m_verbose = v; }

    /**
     * @brief Initialize the inital guess with a given function
     *
     * @param func given function
     */
    void initial_guess( const vector_rhs_function< mesh_type > func ) {
        size_t cell_i = 0;

        m_fields.createField( 0, FieldName::DEPL, m_msh, m_degree_infos, func );
        m_fields.createField( 0, FieldName::DEPL_CELLS, m_msh, m_degree_infos, func );
        m_fields.createField( 0, FieldName::DEPL_FACES, m_msh, m_degree_infos, func );

        /*TODO: look for contact.*/
    }

    /**
     * @brief Initialize displacement and velocity with a given function
     *
     * @param func given function
     */
    void initial_field( FieldName name, const vector_rhs_function< mesh_type > func ) {
        if ( name == FieldName::DEPL ) {
            m_fields.createField( 0, FieldName::DEPL_CELLS, m_msh, m_degree_infos, func );
            m_fields.createField( 0, FieldName::DEPL_FACES, m_msh, m_degree_infos, func );
        }
        if ( name == FieldName::VITE ) {
            m_fields.createField( 0, FieldName::VITE_CELLS, m_msh, m_degree_infos, func );
        }
        if ( name == FieldName::ACCE ) {
            m_fields.createField( 0, FieldName::ACCE_CELLS, m_msh, m_degree_infos, func );
        }
        m_fields.createField( 0, name, m_msh, m_degree_infos, func );
    }

    /**
     * @brief Add a behavior for materials
     *
     * @param deformation Type of deformation
     * @param law Type of Law
     */
    void addBehavior( const size_t deformation, const size_t law ) {
        if ( m_verbose ) {
            std::cout << std::endl;
            std::cout << "Add behavior ..." << std::endl;
        }

        const auto mater = m_behavior.getMaterialData();
        m_behavior = behavior_type( m_msh, 2 * m_rp.m_grad_degree, deformation, law );
        m_behavior.addMaterialData( mater );

        if ( m_verbose ) {
            std::cout << "** Deformations: " << m_behavior.getDeformationName() << std::endl;
            std::cout << "** Law: " << m_behavior.getLawName() << std::endl;
            std::cout << "** Number of integration points: " << m_behavior.numberOfQP()
                      << std::endl;
        }
    }

#ifdef HAVE_MGIS
    /**
     * @brief Add a behavior for materials
     *
     * @param deformation Type of deformation
     * @param law Type of Law
     */
    void addBehavior( const std::string &filename, const std::string &law,
                      const mgis::behaviour::Hypothesis h ) {
        if ( m_verbose ) {
            std::cout << std::endl;
            std::cout << "Add behavior ..." << std::endl;
        }

        const auto mater = m_behavior.getMaterialData();
        m_behavior = behavior_type( m_msh, 2 * m_rp.m_grad_degree, filename, law, h );
        m_behavior.addMaterialData( mater );

        if ( m_verbose ) {
            std::cout << "** Deformations: " << m_behavior.getDeformationName() << std::endl;
            std::cout << "** Law: " << m_behavior.getLawName() << std::endl;
            std::cout << "** Number of integration points: " << m_behavior.numberOfQP()
                      << std::endl;
        }
    }
#endif

    /**
     * @brief Add a behavior for materials (by copy)
     *
     * @param behavior Given behavior
     */
    void addBehavior( const behavior_type &behavior ) {
        m_behavior = behavior;
        if ( m_verbose ) {
            std::cout << std::endl;
            std::cout << "Add behavior ..." << std::endl;
            std::cout << "** Number of integration points: " << m_behavior.numberOfQP()
                      << std::endl;
        }
    }

    /**
     * @brief Add material properties for the behavior
     *
     * @param material_data material properties
     */
    void addMaterialData( const MaterialData< scalar_type > &material_data ) {
        m_behavior.addMaterialData( material_data );

        if ( m_verbose ) {
            std::cout << "Add material ..." << std::endl;
            m_behavior.getMaterialData().print();
        }
    }

    void addPointPlot( const point_type &pt, const std::string &filename ) {
        auto ppt = ppt_type( m_msh, pt );

        std::vector< std::string > cmps;
        cmps.push_back( "n_iter" );
        cmps.push_back( "x" );
        cmps.push_back( "y" );
        if constexpr ( mesh_type::dimension == 3 ) {
            cmps.push_back( "z" );
        }
        cmps.push_back( "ux" );
        cmps.push_back( "uy" );
        if constexpr ( mesh_type::dimension == 3 ) {
            cmps.push_back( "uz" );
        }

        if ( m_rp.isUnsteady() ) {
            cmps.push_back( "vx" );
            cmps.push_back( "vy" );
            if constexpr ( mesh_type::dimension == 3 ) {
                cmps.push_back( "vz" );
            }
            cmps.push_back( "ax" );
            cmps.push_back( "ay" );
            if constexpr ( mesh_type::dimension == 3 ) {
                cmps.push_back( "az" );
            }
        }

        // stress tensor

        cmps.push_back( "sxx" );
        cmps.push_back( "syy" );
        cmps.push_back( "szz" );
        cmps.push_back( "sxy" );

        if constexpr ( mesh_type::dimension == 3 ) {
            cmps.push_back( "sxz" );
            cmps.push_back( "syz" );
        }

        ppt.addComponents( cmps );
        ppt.setFilename( filename );

        m_ppt.push_back( ppt );
    }

    void addExternalLoad( const std::unique_ptr< func_type > &load ) { m_load = load; }
    void addExternalLoad( const func_type load ) { m_load = std::make_unique<func_type>( load ); }

    SolverInfo compute() {
        // Precomputation
        if ( m_rp.m_precomputation ) {
            timecounter t1;
            t1.tic();
            this->pre_computation();
            t1.toc();
            if ( m_verbose )
                std::cout << "Precomputation: " << t1.elapsed() << " sec" << std::endl;
        }

        if ( m_rp.isUnsteady() ) {
            reformulation_dynamic( m_rp );
            m_rp.m_dyna_para["rho"] = m_behavior.getMaterialData().getRho();
        }

        // save first state;
        for ( auto &ppt : m_ppt ) {
            std::vector< static_vector< scalar_type, mesh_type::dimension > > vals;
            auto depl = _eval( FieldName::DEPL_CELLS, ppt.getCellId(), ppt.getPoint() );
            // new coords
            static_vector< scalar_type, mesh_type::dimension > px;
            for ( int i = 0; i < mesh_type::dimension; i++ ) {
                px[i] = ppt.getPoint()[i] + depl[i];
            }
            vals.push_back( px );
            vals.push_back( depl );

            if ( m_rp.isUnsteady() ) {
                auto vite = _eval( FieldName::VITE_CELLS, ppt.getCellId(), ppt.getPoint() );
                vals.push_back( vite );
                auto acce = _eval( FieldName::ACCE_CELLS, ppt.getCellId(), ppt.getPoint() );
                vals.push_back( acce );
            }

            // stress tensor
            const auto vzero = static_vector< scalar_type, mesh_type::dimension >::Zero();
            vals.push_back( vzero );
            vals.push_back( vzero );

            ppt.addValues( 0.0, 0, vals );
        }

        SolverInfo si;
        ppt_type stat;
        stat.setFilename( "statistics.csv" );

        timecounter ttot;
        ttot.tic();

        // check CFL condition
        scalar_type h_min = minimum_diameter( m_msh );
        scalar_type vp = 1.0;
        if ( m_rp.isUnsteady() && m_rp.getUnsteadyScheme() == DynamicType::LEAP_FROG ) {
            const auto &mater = m_behavior.getMaterialData();
            vp = std::sqrt( ( mater.getLambda() + 2.0 * mater.getMu() ) / mater.getRho() );
        }

        // list of time step
        ListOfTimeStep< scalar_type > list_time_step;
        if ( m_rp.m_has_user_end_time )
            list_time_step =
                ListOfTimeStep< scalar_type >( m_rp.m_time_step, m_rp.m_user_end_time );
        else
            list_time_step = ListOfTimeStep< scalar_type >( m_rp.m_time_step );

        if ( m_verbose )
            std::cout << "** Number of time step: " << list_time_step.numberOfTimeStep()
                      << std::endl;

        // time of saving
        bool time_saving = false;
        if ( m_rp.m_n_time_save > 0 ) {
            time_saving = true;
        }

        NewtonSolverInfo newton_info;

        // update initial state;
        m_fields.update();

        // Loop on time step
        while ( !list_time_step.empty() ) {
            const auto current_step = list_time_step.getCurrentTimeStep();
            const auto current_time = current_step.end_time();
            m_fields.setCurrentTime( current_time );

            if ( m_rp.getUnsteadyScheme() == DynamicType::LEAP_FROG ) {
                const auto dt_crit = ( h_min / vp ) * m_rp.getCFLFactor();
                if ( current_step.increment_time() > dt_crit ) {
                    throw std::runtime_error(
                        "CFL is not respected. dt_crit=" + std::to_string( dt_crit ) +
                        " vs dt=" + std::to_string( current_step.increment_time() ) );
                }
            }

            if ( m_verbose ) {
                list_time_step.printCurrentTimeStep();
            }

            m_bnd.setTime( current_time );

            NonLinearStep< mesh_type > nlStep( m_rp );

            newton_info =
                nlStep.compute( m_msh, m_bnd, m_rp, m_degree_infos, m_lin_solv, m_load,
                                current_step, m_data, m_behavior, m_stab_manager, m_fields );

            // Test convergence
            m_convergence = nlStep.convergence();

            //  Newton correction
            si.updateInfo( newton_info );

            if ( m_verbose ) {
                newton_info.printInfo();
            }

            if ( !m_convergence ) {
                if ( current_step.level() + 1 > m_rp.m_sublevel ) {
                    std::cout << "***********************************************************"
                              << std::endl;
                    std::cout << "***** PROBLEM OF CONVERGENCE: We stop the calcul here *****"
                              << std::endl;
                    std::cout << "***********************************************************"
                              << std::endl;
                    break;
                } else {
                    if ( m_verbose ) {
                        std::cout << "***********************************************************"
                                  << std::endl;
                        std::cout << "*****     NO CONVERGENCE: We split the time step     ******"
                                  << std::endl;
                        std::cout << "***********************************************************"
                                  << std::endl;
                    }

                    list_time_step.splitCurrentTimeStep();
                    m_behavior.restore();
                    m_fields.restore();
                }
            } else {
                list_time_step.removeCurrentTimeStep();
                m_behavior.update();
                m_stab_manager.update();
                m_fields.update();

                if ( time_saving && ( m_rp.m_time_save.front() < current_time + 1E-5 ) ) {
                    std::cout << "** Save results" << std::endl;
                    std::string name = "result" + std::to_string( mesh_type::dimension ) + "D_t" +
                                       std::to_string( current_time ) + "_";

                    this->output_discontinuous_field( name + "depl_disc.msh",
                                                      FieldName::DEPL_CELLS );
                    this->output_continuous_field( name + "depl_cont.msh", FieldName::DEPL_CELLS );
                    this->output_CauchyStress_GP( name + "CauchyStress_GP.msh" );
                    this->output_CauchyStress_GP( name + "CauchyStress_GP_def.msh", true );
                    this->output_discontinuous_deformed( name + "deformed_disc.msh" );
                    this->output_is_plastic_GP( name + "plastic_GP.msh" );
                    this->output_stabCoeff( name + "stabCoeff.msh" );
                    this->output_equivalentPlasticStrain_GP( name +
                                                             "equivalentPlasticStrain_GP.msh" );
                    if ( m_rp.isUnsteady() ) {
                        this->output_discontinuous_field( name + "vite_disc.msh",
                                                          FieldName::VITE_CELLS );
                        this->output_continuous_field( name + "vite_cont.msh",
                                                       FieldName::VITE_CELLS );
                        this->output_discontinuous_field( name + "acce_disc.msh",
                                                          FieldName::ACCE_CELLS );
                        this->output_continuous_field( name + "acce_cont.msh",
                                                       FieldName::ACCE_CELLS );
                    }

                    m_rp.m_time_save.pop_front();
                    if ( m_rp.m_time_save.empty() )
                        time_saving = false;
                }

                // Compute observation
                for ( auto &ppt : m_ppt ) {
                    std::vector< static_vector< scalar_type, mesh_type::dimension > > vals;
                    auto depl = _eval( FieldName::DEPL_CELLS, ppt.getCellId(), ppt.getPoint() );
                    // new coords
                    static_vector< scalar_type, mesh_type::dimension > px;
                    for ( int i = 0; i < mesh_type::dimension; i++ ) {
                        px[i] = ppt.getPoint()[i] + depl[i];
                    }

                    vals.push_back( px );
                    vals.push_back( depl );

                    if ( m_rp.isUnsteady() ) {
                        auto vite = _eval( FieldName::VITE_CELLS, ppt.getCellId(), ppt.getPoint() );
                        vals.push_back( vite );
                        auto acce = _eval( FieldName::ACCE_CELLS, ppt.getCellId(), ppt.getPoint() );
                        vals.push_back( acce );
                    }

                    // stress
                    const auto cl = *std::next( m_msh.cells_begin(), ppt.getCellId() );
                    const auto di = m_degree_infos.cellDegreeInfo( m_msh, cl );
                    const auto stress =
                        m_behavior.projectStressOnCell( m_msh, cl, di.grad_degree() );

                    const auto gb = make_matrix_monomial_basis( m_msh, cl, di.grad_degree() );
                    const auto gphi = gb.eval_functions( ppt.getPoint() );
                    const auto GT_iqn = eval( stress, gphi );

                    static_vector< scalar_type, mesh_type::dimension > sdiag, sshea;

                    if constexpr ( mesh_type::dimension == 2 ) {
                        sdiag( 0 ) = GT_iqn( 0, 0 );
                        sdiag( 1 ) = GT_iqn( 1, 1 );
                        sshea( 0 ) = 0.0;
                        sshea( 1 ) = GT_iqn( 0, 1 );
                    } else {
                        sdiag( 0 ) = GT_iqn( 0, 0 );
                        sdiag( 1 ) = GT_iqn( 1, 1 );
                        sdiag( 2 ) = GT_iqn( 2, 2 );
                        sshea( 0 ) = GT_iqn( 0, 1 );
                        sshea( 1 ) = GT_iqn( 0, 2 );
                        sshea( 2 ) = GT_iqn( 1, 2 );
                    }
                    vals.push_back( sdiag );
                    vals.push_back( sshea );

                    ppt.addValues( current_time, newton_info.m_iter, vals );
                }

                // Update stats
                stat.addValues( current_time, newton_info.getValues() );
            }
        }

        for ( auto &ppt : m_ppt ) {
            ppt.write();
        }
        stat.write();

        si.m_time_step = list_time_step.numberOfTimeStep();

        ttot.toc();
        si.m_time_solver = ttot.elapsed();

        return si;
    }

    bool convergence() const { return m_convergence; }

    size_t numberOfDofs() {
        const auto dimension = mesh_type::dimension;
        size_t num_faces_dofs = 0;
        for ( auto itor = m_msh.faces_begin(); itor != m_msh.faces_end(); itor++ ) {
            const auto fc = *itor;
            const auto di = m_degree_infos.degreeInfo( m_msh, fc );

            if ( di.hasUnknowns() ) {
                num_faces_dofs += vector_basis_size( di.degree(), dimension - 1, dimension );
            }
        }
        return num_faces_dofs;
    }

    void printSolutionCell() const {
        size_t cell_i = 0;
        const auto depl_cells = m_fields.getCurrentField( FieldName::DEPL_CELLS );

        std::cout << "Solution at the cells:" << std::endl;
        for ( auto &cl : m_msh ) {
            std::cout << "cell " << cell_i << ": " << std::endl;
            std::cout << depl_cells.at( cell_i++ ).transpose() << std::endl;
        }
    }

    // compute l2 error
    template < typename AnalyticalSolution >
    long double compute_l2_error( FieldName name, const AnalyticalSolution &as ) {
        using quad_type = long double;
        quad_type err_dof = 0.;

        size_t cell_i = 0;
        const auto field = m_fields.getCurrentField( name );

        for ( auto &cl : m_msh ) {
            const auto cdi = m_degree_infos.degreeInfo( m_msh, cl );
            const vector_type comp_dof = field.at( cell_i++ );
            const vector_type true_dof = project_function( m_msh, cl, cdi.degree(), as, 2 );

            const auto cb = make_vector_monomial_basis( m_msh, cl, cdi.degree() );
            const matrix_type mass = make_mass_matrix( m_msh, cl, cb );

            const vector_type diff_dof = ( true_dof - comp_dof );
            const vector_type mass_diff = mass * diff_dof;
            const auto size = diff_dof.size();

            for ( int i = 0; i < size; i++ ) {
                const quad_type x = static_cast< quad_type >( diff_dof( i ) );
                const quad_type y = static_cast< quad_type >( mass_diff( i ) );
                err_dof = std::fma( x, y, err_dof );
            }
        }

        return std::sqrt( err_dof );
    }

    // compute l2 error
    template < typename AnalyticalSolution >
    scalar_type compute_l2_displacement_error( const AnalyticalSolution &as ) {
        return compute_l2_error( FieldName::DEPL_CELLS, as );
    }

    // compute l2 error
    template < typename AnalyticalSolution >
    long double compute_H1_error( const AnalyticalSolution &as ) {
        using quad_type = long double;
        quad_type err_dof = 0.;

        const auto depl = m_fields.getCurrentField( FieldName::DEPL );

        matrix_type grad;
        matrix_type stab;

        for ( auto &cl : m_msh ) {
            // std::cout << m_degree_infos.cellDegreeInfo(m_msh, cl) << std::endl;
            /////// Gradient Reconstruction /////////
            if ( m_behavior.getDeformation() == SMALL_DEF ) {
                grad = make_matrix_hho_symmetric_gradrec( m_msh, cl, m_degree_infos ).second;
            } else {
                grad = make_matrix_hho_gradrec( m_msh, cl, m_degree_infos ).second;
            }

            if ( m_rp.m_stab ) {
                switch ( m_rp.m_stab_type ) {
                case StabilizationType::HHO_SYM: {
                    const auto recons =
                        make_vector_hho_symmetric_laplacian( m_msh, cl, m_degree_infos );
                    stab = make_vector_hho_stabilization( m_msh, cl, recons.first, m_degree_infos );
                    break;
                }
                case StabilizationType::HHO: {
                    const auto recons_scalar =
                        make_scalar_hho_laplacian( m_msh, cl, m_degree_infos );
                    stab = make_vector_hho_stabilization_optim( m_msh, cl, recons_scalar.first,
                                                                m_degree_infos );
                    break;
                }
                case StabilizationType::HDG: {
                    stab = make_vector_hdg_stabilization( m_msh, cl, m_degree_infos );
                    break;
                }
                case StabilizationType::DG: {
                    stab = make_vector_dg_stabilization( m_msh, cl, m_degree_infos );
                    break;
                }
                case StabilizationType::NO: {
                    break;
                    stab.setZero();
                }
                default:
                    throw std::invalid_argument( "Unknown stabilization" );
                }
            }

            const auto Ah = grad + stab;

            const auto cell_i = m_msh.lookup( cl );

            const vector_type comp_dof = depl.at( cell_i );
            const vector_type true_dof = project_function( m_msh, cl, m_degree_infos, as, 2 );

            const vector_type diff_dof = ( true_dof - comp_dof );
            const vector_type Ah_diff_dof = Ah * diff_dof;

            const auto size = diff_dof.size();
            for ( int i = 0; i < size; i++ ) {
                const quad_type x = static_cast< quad_type >( diff_dof( i ) );
                const quad_type y = static_cast< quad_type >( Ah_diff_dof( i ) );
                err_dof = std::fma( x, y, err_dof );
            }
        }

        return std::sqrt( err_dof );
    }

    void output_discontinuous_field( const std::string &filename, FieldName name ) const {
        gmsh::Gmesh gmsh( mesh_type::dimension );

        std::vector< gmsh::Data > data;             // create data (not used)
        const std::vector< gmsh::SubData > subdata; // create subdata to save soution at gauss point

        const auto depl_cells = m_fields.getCurrentField( name );

        int cell_i = 0;
        int nb_nodes = 0;
        for ( auto &cl : m_msh ) {
            const auto di = m_degree_infos.cellDegreeInfo( m_msh, cl );
            const auto cb = make_vector_monomial_basis( m_msh, cl, di.cell_degree() );
            const vector_type x = depl_cells.at( cell_i++ );
            auto cell_nodes = points( m_msh, cl );
            std::vector< gmsh::Node > new_nodes;

            // loop on the nodes of the cell
            for ( auto &pt : cell_nodes ) {
                nb_nodes++;

                const auto phi = cb.eval_functions( pt );
                const auto depl = eval( x, phi );

                const std::vector< double > deplv = convertToVectorGmsh( depl );
                const std::array< double, 3 > coor = init_coor( pt );

                // Add a node
                const gmsh::Node tmp_node( coor, nb_nodes, 0 );
                new_nodes.push_back( tmp_node );
                gmsh.addNode( tmp_node );

                const gmsh::Data datatmp( nb_nodes, deplv );
                data.push_back( datatmp );
            }
            // Add new element
            add_element( gmsh, new_nodes );
        }

        // Create and init a nodedata view
        gmsh::NodeData nodedata( 3, 0.0, "depl_node_disc", data, subdata );

        // Save the view
        nodedata.saveNodeData( filename, gmsh );
    }

    void output_continuous_field( const std::string &filename, FieldName name ) const {
        const auto dimension = mesh_type::dimension;

        gmsh::Gmesh gmsh = convertMesh( m_post_mesh );
        auto storage = m_post_mesh.mesh().backend_storage();

        const static_vector< scalar_type, dimension > vzero =
            static_vector< scalar_type, dimension >::Zero();

        const auto depl_cells = m_fields.getCurrentField( FieldName::DEPL_CELLS );
        const size_t nb_nodes( gmsh.getNumberofNodes() );

        // first(number of data at this node), second(cumulated value)
        std::vector< std::pair< size_t, static_vector< scalar_type, dimension > > > value(
            nb_nodes, std::make_pair( 0, vzero ) );

        int cell_i = 0;
        for ( auto &cl : m_msh ) {
            const auto di = m_degree_infos.cellDegreeInfo( m_msh, cl );
            const auto cb = make_vector_monomial_basis( m_msh, cl, di.cell_degree() );
            const vector_type x = depl_cells.at( cell_i );
            auto cell_nodes = m_post_mesh.nodes_cell( cell_i );

            // Loop on the nodes of the cell
            for ( auto &point_id : cell_nodes ) {
                const auto pt = storage->points[point_id];

                const auto phi = cb.eval_functions( pt );
                const auto depl = eval( x, phi );

                // Add displacement at node
                value[point_id].first++;
                value[point_id].second += depl;
            }
            cell_i++;
        }

        std::vector< gmsh::Data > data;       // create data
        std::vector< gmsh::SubData > subdata; // create subdata
        data.reserve( nb_nodes );             // data has a size of nb_node

        // Compute the average value and save it
        for ( int i_node = 0; i_node < value.size(); i_node++ ) {
            const static_vector< scalar_type, dimension > depl_avr =
                value[i_node].second / double( value[i_node].first );

            const gmsh::Data tmp_data( i_node + 1, convertToVectorGmsh( depl_avr ) );
            data.push_back( tmp_data );
        }

        // Create and init a nodedata view
        gmsh::NodeData nodedata( 3, 0.0, "depl_node_cont", data, subdata );
        // Save the view
        nodedata.saveNodeData( filename, gmsh );
    }

    void output_CauchyStress_GP( const std::string &filename, bool def = false ) const {
        gmsh::Gmesh gmsh = convertMesh( m_post_mesh );

        std::vector< gmsh::Data > data;       // create data (not used)
        std::vector< gmsh::SubData > subdata; // create subdata to save soution at gauss point
        size_t nb_nodes( gmsh.getNumberofNodes() );

        const auto depl = m_fields.getCurrentField( FieldName::DEPL );

        int cell_i = 0;
        for ( auto &cl : m_msh ) {
            const auto di = m_degree_infos.cellDegreeInfo( m_msh, cl );

            const auto uTF = depl.at( cell_i );
            matrix_type gr;
            if ( m_rp.m_precomputation ) {
                gr = m_data.m_gradient_precomputed.at( cell_i );
            } else {
                if ( m_behavior.getDeformation() == SMALL_DEF ) {
                    gr = make_matrix_hho_symmetric_gradrec( m_msh, cl, m_degree_infos ).first;
                } else {
                    gr = make_matrix_hho_gradrec( m_msh, cl, m_degree_infos ).first;
                }
            }

            const vector_type GTuTF = gr * uTF;

            const auto gb = make_matrix_monomial_basis( m_msh, cl, di.grad_degree() );
            const auto gbs = make_sym_matrix_monomial_basis( m_msh, cl, di.grad_degree() );

            const auto cb = make_vector_monomial_basis( m_msh, cl, di.cell_degree() );
            const vector_type uT = uTF.head( cb.size() );

            // Loop on nodes
            const auto nb_qp = m_behavior.numberOfQP( cell_i );

            for ( int i_qp = 0; i_qp < nb_qp; i_qp++ ) {
                const auto qp = m_behavior.quadrature_point( cell_i, i_qp );
                std::vector< double > tens;

                if ( m_behavior.getDeformation() == SMALL_DEF ) {
                    auto stress = m_behavior.compute_stress3D( cell_i, i_qp );
                    tens = convertToVectorGmsh( stress );
                } else {
                    const auto gphi = gb.eval_functions( qp.point() );
                    const auto GT_iqn = eval( GTuTF, gphi );
                    const auto FT_iqn = convertGtoF( GT_iqn );
                    const auto FT_iqn_3D = convertMatrix3DwithOne( FT_iqn );

                    auto P = m_behavior.compute_stress3D( cell_i, i_qp );
                    auto stress = convertPK1toCauchy( P, FT_iqn_3D );
                    tens = convertToVectorGmsh( stress );
                }

                std::array< double, 3 > coor = init_coor( qp.point() );

                if ( def ) {
                    const auto cphi = cb.eval_functions( qp.point() );
                    const auto depl = eval( uT, cphi );

                    // Compute new coordinates
                    for ( int j = 0; j < mesh_type::dimension; j++ )
                        coor[j] += depl( j );
                }

                // Add GP
                // Create a node at gauss point
                nb_nodes++;
                const gmsh::Node new_node( coor, nb_nodes, 0 );
                const gmsh::SubData sdata( tens, new_node );
                subdata.push_back( sdata ); // add subdata
            }
            cell_i++;
        }

        // Save
        gmsh::NodeData nodedata( 9, 0.0, "CauchyStress_GP", data,
                                 subdata ); // create and init a nodedata view

        nodedata.saveNodeData( filename, gmsh ); // save the view
    }

    void output_is_plastic_GP( const std::string &filename ) const {
        gmsh::Gmesh gmsh = convertMesh( m_post_mesh );

        std::vector< gmsh::Data > data;       // create data (not used)
        std::vector< gmsh::SubData > subdata; // create subdata to save soution at gauss point
        size_t nb_nodes( gmsh.getNumberofNodes() );

        int cell_i = 0;
        for ( auto &cl : m_msh ) {
            // Loop on nodes
            const auto nb_qp = m_behavior.numberOfQP( cell_i );

            for ( int i_qp = 0; i_qp < nb_qp; i_qp++ ) {
                const auto qp = m_behavior.quadrature_point( cell_i, i_qp );

                scalar_type p = 0;
                if ( m_behavior.is_plastic( cell_i, i_qp ) )
                    p = 1;

                const std::vector< double > p_s = convertToVectorGmsh( p );

                // Add GP
                // Create a node at gauss point
                nb_nodes++;
                const gmsh::Node new_node = convertPoint( qp.point(), nb_nodes );
                const gmsh::SubData sdata( p_s, new_node );
                subdata.push_back( sdata ); // add subdata
            }
            cell_i++;
        }

        // Save
        gmsh::NodeData nodedata( 1, 0.0, "state_GP", data,
                                 subdata ); // create and init a nodedata view

        nodedata.saveNodeData( filename, gmsh ); // save the view
    }

    void output_equivalentPlasticStrain_GP( const std::string &filename ) const {
        gmsh::Gmesh gmsh = convertMesh( m_post_mesh );

        std::vector< gmsh::Data > data;       // create data (not used)
        std::vector< gmsh::SubData > subdata; // create subdata to save soution at gauss point
        size_t nb_nodes( gmsh.getNumberofNodes() );

        int cell_i = 0;
        for ( auto &cl : m_msh ) {
            // Loop on nodes
            const auto nb_qp = m_behavior.numberOfQP( cell_i );

            for ( int i_qp = 0; i_qp < nb_qp; i_qp++ ) {
                const auto qp = m_behavior.quadrature_point( cell_i, i_qp );

                scalar_type p = m_behavior.equivalentPlasticStrain( cell_i, i_qp );

                const std::vector< double > p_s = convertToVectorGmsh( p );

                // Add GP
                // Create a node at gauss point
                nb_nodes++;
                const gmsh::Node new_node = convertPoint( qp.point(), nb_nodes );
                const gmsh::SubData sdata( p_s, new_node );
                subdata.push_back( sdata ); // add subdata
            }
            cell_i++;
        }

        // Save
        gmsh::NodeData nodedata( 1, 0.0, "equivalentPlasticStrain_GP", data,
                                 subdata ); // create and init a nodedata view

        nodedata.saveNodeData( filename, gmsh ); // save the view
    }

    void output_discontinuous_deformed( const std::string &filename ) const {
        gmsh::Gmesh gmsh( mesh_type::dimension );
        auto storage = m_msh.backend_storage();

        const auto depl_cells = m_fields.getCurrentField( FieldName::DEPL_CELLS );
        int cell_i = 0;
        size_t nb_nodes = 0;
        for ( auto &cl : m_msh ) {
            const auto di = m_degree_infos.cellDegreeInfo( m_msh, cl );

            auto cb = make_vector_monomial_basis( m_msh, cl, di.cell_degree() );
            const vector_type x = depl_cells.at( cell_i++ );
            const auto cell_nodes = points( m_msh, cl );
            std::vector< gmsh::Node > new_nodes;

            // Loop on nodes of the cell
            for ( auto &pt : cell_nodes ) {
                nb_nodes++;

                const auto phi = cb.eval_functions( pt );
                const auto depl = eval( x, phi );

                std::array< double, 3 > coor = init_coor( pt );
                // Compute new coordinates
                for ( int j = 0; j < mesh_type::dimension; j++ )
                    coor[j] += depl( j );

                // Save node
                const gmsh::Node tmp_node( coor, nb_nodes, 0 );
                new_nodes.push_back( tmp_node );
                gmsh.addNode( tmp_node );
            }
            // Add new element
            add_element( gmsh, new_nodes );
        }
        // Save mesh
        gmsh.writeGmesh( filename, 2 );
    }

    void output_stabCoeff( const std::string &filename ) const {
        gmsh::Gmesh gmsh = convertMesh( m_post_mesh );

        std::vector< gmsh::Data > data;       // create data (not used)
        std::vector< gmsh::SubData > subdata; // create subdata to save soution at gauss point
        size_t nb_nodes( gmsh.getNumberofNodes() );

        for ( auto &cl : m_msh ) {

            std::array< double, 3 > coor = init_coor( barycenter( m_msh, cl ) );
            double beta = m_stab_manager.getValue( m_msh, cl );
            std::vector< double > tens( 1, beta );

            // Add GP
            // Create a node at gauss point
            nb_nodes++;
            const gmsh::Node new_node( coor, nb_nodes, 0 );
            const gmsh::SubData sdata( tens, new_node );
            subdata.push_back( sdata ); // add subdata
        }

        // Save
        gmsh::NodeData nodedata( 1, 0.0, "StabCoeff", data,
                                 subdata ); // create and init a nodedata view

        nodedata.saveNodeData( filename, gmsh ); // save the view
    }
};
} // namespace mechanics

} // namespace disk
