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
#include "diskpp/output/ensight/ensight_exporter.hpp"
#include "diskpp/output/gmshConvertMesh.hpp"
#include "diskpp/output/gmshDisk.hpp"
#include "diskpp/output/plotOverTime.hpp"

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

    // ensight output
    output::ensight::EnsightExporter< mesh_type > m_output;
    std::int64_t m_gauss_point_part_id = -1;
    bool m_ensight_initialized = false;

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

                // order k+1 for contact faces.
                if ( m_bnd.contact_boundary_type( face_id ) == SIGNORINI_FACE ) {
                    m_degree_infos.degree( m_msh, bfc, face_degree + 1 );
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

    auto
    _eval_stress( const int cell_id, const point_type &pt ) {

        // stress
        const auto cl = m_msh[cell_id];
        const auto di = m_degree_infos.cellDegreeInfo( m_msh, cl );
        const auto stress = m_behavior.projectStressOnCell( m_msh, cl, di.grad_degree() );

        const auto gb = make_matrix_monomial_basis( m_msh, cl, di.grad_degree() );
        const auto gphi = gb.eval_functions( pt );
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

        return std::make_pair( sdiag, sshea );
    }

    void
    _setInitialState() {
        const bool small_def = m_behavior.getDeformation() == SMALL_DEF;
        const auto depl = m_fields.getField( 0, FieldName::DEPL );

        for ( auto &cl : m_msh ) {
            const auto c_id = m_msh.lookup( cl );
            const auto di = m_degree_infos.cellDegreeInfo( m_msh, cl );

            const auto uTF = depl.at( c_id );
            matrix_type gr;
            if ( m_rp.m_precomputation ) {
                gr = m_data.m_gradient_precomputed.at( c_id );
            } else {
                if ( small_def ) {
                    gr = make_matrix_hho_symmetric_gradrec( m_msh, cl, m_degree_infos ).first;
                } else {
                    gr = make_matrix_hho_gradrec( m_msh, cl, m_degree_infos ).first;
                }
            }

            const vector_type GTuTF = gr * uTF;

            const auto gb = make_matrix_monomial_basis( m_msh, cl, di.grad_degree() );
            const auto gbs = make_sym_matrix_monomial_basis( m_msh, cl, di.grad_degree() );

            // Loop on nodes
            const auto nb_qp = m_behavior.numberOfQP( c_id );
            eigen_compatible_stdvector<
                static_matrix< scalar_type, mesh_type::dimension, mesh_type::dimension > >
                gphi;

            for ( int i_qp = 0; i_qp < nb_qp; i_qp++ ) {
                const auto qp = m_behavior.quadrature_point( c_id, i_qp );

                if ( small_def ) {
                    gphi = gbs.eval_functions( qp.point() );
                } else {
                    gphi = gb.eval_functions( qp.point() );
                }

                const auto GkT_iqn = eval( GTuTF, gphi );

                m_behavior.setInitialElasticStrain( c_id, i_qp, GkT_iqn );
            }
        }
    }

    void
    initialize_ensight_output() {
        if ( m_ensight_initialized ) {
            return;
        }

        /*
         * Compute the total number of Gauss points.
         *
         * The number of quadrature points may depend on the cell.
         */
        std::size_t total_number_of_qp = 0;

        for ( const auto &cl : m_msh ) {
            const std::size_t cell_id = m_msh.lookup( cl );

            total_number_of_qp += m_behavior.numberOfQP( cell_id );
        }

        if ( total_number_of_qp == 0 ) {
            throw std::runtime_error( "No quadrature points to export" );
        }

        /*
         * reserve() allocates memory without creating entries.
         * Each Gauss point is then appended with push_back().
         */
        std::vector< std::array< double, 3 > > gauss_points;

        gauss_points.reserve( total_number_of_qp );

        for ( const auto &cl : m_msh ) {
            const std::size_t cell_id = m_msh.lookup( cl );

            const std::size_t number_of_qp = m_behavior.numberOfQP( cell_id );

            for ( std::size_t qp_id = 0; qp_id < number_of_qp; ++qp_id ) {
                const auto qp = m_behavior.quadrature_point( cell_id, qp_id );

                std::array< double, 3 > coordinates { 0.0, 0.0, 0.0 };

                coordinates[0] = qp.point().x();

                if constexpr ( mesh_type::dimension >= 2 ) {
                    coordinates[1] = qp.point().y();
                }

                if constexpr ( mesh_type::dimension == 3 ) {
                    coordinates[2] = qp.point().z();
                }

                gauss_points.push_back( coordinates );
            }
        }

        /*
         * Internal consistency check.
         */
        if ( gauss_points.size() != total_number_of_qp ) {
            throw std::runtime_error( "Unexpected number of Gauss points: got " +
                                      std::to_string( gauss_points.size() ) + ", expected " +
                                      std::to_string( total_number_of_qp ) );
        }

        m_gauss_point_part_id =
            m_output.add_point_cloud( "Gauss points", std::move( gauss_points ) );

        m_output.write_mesh();

        m_ensight_initialized = true;
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
              std::make_shared< solvers::sparse_solver< scalar_type > >( rp.getLinearSolver() ) ),
          m_output( msh ) {
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

        // save parameters
        m_rp.writeParameters( m_output.output_directory() / "parameters.dat" );

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
            _setInitialState();
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
            const auto [sdiag, sshea] = _eval_stress( ppt.getCellId(), ppt.getPoint() );
            vals.push_back( sdiag );
            vals.push_back( sshea );

            ppt.addValues( 0.0, 0, vals );
        }

        SolverInfo si;
        ppt_type stat;
        stat.setFilename( "statistics.csv" );

        ppt_type energy_ppt;
        energy_ppt.setFilename( "energy.csv" );

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

        if ( m_rp.isUnsteady() ) {
            list_time_step.checkConstantTimeStep();
        }

        if ( m_verbose )
            std::cout << "** Number of time step: " << list_time_step.numberOfTimeStep()
                      << std::endl;

        // time of saving
        bool time_saving = false, time_saving_all = false;
        if ( m_rp.m_n_time_save > 0 ) {
            time_saving = true;
        } else if ( m_rp.m_n_time_save < 0 ) {
            time_saving_all = true;
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

            newton_info = nlStep.compute( m_msh,
                                          m_bnd,
                                          m_rp,
                                          m_degree_infos,
                                          m_lin_solv,
                                          m_load,
                                          current_step,
                                          m_data,
                                          m_behavior,
                                          m_stab_manager,
                                          m_fields );

            // Test convergence
            m_convergence = nlStep.convergence();

            //  Newton correction
            si.updateInfo( newton_info );

            if ( m_verbose and m_rp.getNonLinearSolver() != NonLinearSolverType::EXPLICIT ) {
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

                if ( time_saving_all ||
                     ( time_saving && ( m_rp.m_time_save.front() < current_time + 1E-5 ) ) ) {
                    initialize_ensight_output();
                    std::cout << "** Save results" << std::endl;
                    const auto gmsh_directory = m_output.output_directory() / "gmsh";

                    std::filesystem::create_directories( gmsh_directory );

                    std::ostringstream basename;
                    basename << "result" << mesh_type::dimension << "D_t"
                             << std::to_string( current_time ) << '_';

                    const auto filepath = gmsh_directory / basename.str();

                    m_output.begin_step( current_time );

                    this->output_discontinuous_field( filepath.string() + "depl_disc.msh",
                                                      FieldName::DEPL_CELLS );

                    this->output_continuous_field( "displacement", FieldName::DEPL_CELLS );
                    this->output_CauchyStress_GP();
                    this->output_GaussDisplacement_GP();
                    this->output_is_plastic_GP();
                    this->output_stabCoeff();
                    this->output_equivalentPlasticStrain_GP();

                    this->output_normal_stress_boundary_nodes( false );
                    if ( m_bnd.nb_faces_contact() > 0 ) {
                        this->output_normal_stress_boundary_nodes( true );
                    }
                    if ( m_rp.isUnsteady() ) {
                        this->output_continuous_field( "velocity", FieldName::VITE_CELLS );
                        this->output_continuous_field( "acceleration", FieldName::ACCE_CELLS );
                    }

                    this->output_discontinuous_deformed( filepath.string() + "deformed_disc.msh" );

                    if ( time_saving ) {
                        m_rp.m_time_save.pop_front();
                        if ( m_rp.m_time_save.empty() )
                            time_saving = false;
                    }

                    m_output.end_step();
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
                    const auto [sdiag, sshea] = _eval_stress( ppt.getCellId(), ppt.getPoint() );
                    vals.push_back( sdiag );
                    vals.push_back( sshea );

                    ppt.addValues( current_time, newton_info.m_iter, vals );
                }

                // Update stats
                stat.addValues( current_time, newton_info.getValues() );

                // Discrete HHO energy
                if ( m_rp.isUnsteady() ) {
                    const auto E = compute_discrete_energy();
                    std::map< std::string, double > emap;
                    emap["kinetic"] = E.kinetic;
                    emap["elastic"] = E.elastic;
                    emap["stab"] = E.stab;
                    emap["contact"] = E.contact;
                    emap["friction"] = E.friction;
                    emap["total"] = E.total();
                    energy_ppt.addValues( current_time, emap );
                }
            }
        }

        // CSV files
        const auto csv_directory = m_output.output_directory() / "csv_files";
        for ( auto &ppt : m_ppt ) {
            ppt.write( csv_directory );
        }
        stat.write( csv_directory );
        energy_ppt.write( csv_directory );
        if ( m_bnd.nb_faces_contact() > 0 ) {
            const auto contact_file = csv_directory / "contact.csv";
            this->output_contact_boundary( contact_file.string() );
        }

        si.m_time_step = list_time_step.numberOfTimeStep();

        ttot.toc();
        si.m_time_solver = ttot.elapsed();

        m_output.validate_and_write_case();

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

    // discrete HHO mechanical energy: 1/2 sum_T [ v^T M_T v + int_T sigma(u):eps(u)
    // + beta_s (S_T u_TF, u_TF) ], plus the Nitsche contact and friction terms

    struct DiscreteEnergy {
        scalar_type kinetic = 0;
        scalar_type elastic = 0;
        scalar_type stab = 0;
        scalar_type contact = 0;  // Nitsche, normal Signorini part
        scalar_type friction = 0; // Nitsche, tangential part
        scalar_type
        total() const {
            return kinetic + elastic + stab + contact + friction;
        }
    };

    DiscreteEnergy
    compute_discrete_energy() const {
        DiscreteEnergy E;

        const auto depl = m_fields.getCurrentField( FieldName::DEPL );

        const bool have_vite = m_rp.isUnsteady();
        std::vector< vector_type > vite;
        if ( have_vite )
            vite = m_fields.getCurrentField( FieldName::VITE_CELLS );

        const auto &mat = m_behavior.getMaterialData();
        const scalar_type mu = mat.getMu();
        const scalar_type lambda = mat.getLambda();
        const scalar_type rho = mat.getRho();

        for ( auto &cl : m_msh ) {
            const auto cell_i = m_msh.lookup( cl );
            const auto di = m_degree_infos.cellDegreeInfo( m_msh, cl );

            const auto cb = make_vector_monomial_basis( m_msh, cl, di.cell_degree() );
            const vector_type uTF = depl.at( cell_i );

            // kinetic
            if ( have_vite ) {
                const matrix_type MT = rho * make_mass_matrix( m_msh, cl, cb );
                const vector_type vT = vite.at( cell_i );
                E.kinetic += scalar_type( 0.5 ) * vT.dot( MT * vT );
            }

            // elastic strain energy (small strain)
            if ( m_behavior.getDeformation() == SMALL_DEF ) {
                matrix_type gr;
                if ( m_rp.m_precomputation )
                    gr = m_data.m_gradient_precomputed.at( cell_i );
                else
                    gr = make_matrix_hho_symmetric_gradrec( m_msh, cl, m_degree_infos ).first;

                const vector_type GTuTF = gr * uTF;
                const auto gbs = make_sym_matrix_monomial_basis( m_msh, cl, di.grad_degree() );

                const auto qps = integrate( m_msh, cl, 2 * di.grad_degree() + 1 );
                for ( auto &qp : qps ) {
                    const auto gphi = gbs.eval_functions( qp.point() );
                    const auto eps = eval( GTuTF, gphi );
                    const scalar_type tr = eps.trace();
                    E.elastic += qp.weight() *
                                 ( mu * eps.squaredNorm() + scalar_type( 0.5 ) * lambda * tr * tr );
                }
            }

            // stabilization energy
            if ( m_rp.m_stab ) {
                const matrix_type ST =
                    _stab( m_msh, cl, m_rp, m_degree_infos, m_data.m_stab_precomputed );
                const scalar_type beta_s = m_stab_manager.getValue( m_msh, cl );
                E.stab += scalar_type( 0.5 ) * beta_s * uTF.dot( ST * uTF );
            }

            // Nitsche contact and friction energy
            if ( m_behavior.getDeformation() == SMALL_DEF && m_bnd.cell_has_contact_faces( cl ) ) {
                matrix_type gr_c;
                if ( m_rp.m_precomputation )
                    gr_c = m_data.m_gradient_precomputed.at( cell_i );
                else
                    gr_c = make_matrix_hho_symmetric_gradrec( m_msh, cl, m_degree_infos ).first;

                auto cc = contact_contribution< mesh_type >( m_msh, mat, m_rp, m_bnd );
                E.contact += cc.nitsche_contact_energy( cl, di, gr_c, uTF );
                E.friction += cc.nitsche_friction_energy( cl, di, gr_c, uTF );
            }
        }
        return E;
    }

    // Trace of the whole contact boundary at the current time, one row per contact
    // quadrature point, ordered along the boundary. Every tangential quantity is
    // projected on the tangent of the face's own discrete normal, never on a fixed axis.
    void
    output_contact_boundary( const std::string &filename ) const {
        if constexpr ( mesh_type::dimension != 2 ) {
            std::cout << "output_contact_boundary: 2D only, skipped" << std::endl;
            return;
        } else {
            typedef typename contact_contribution< mesh_type >::trace_point trace_point;

            const auto depl = m_fields.getCurrentField( FieldName::DEPL );
            const auto &mat = m_behavior.getMaterialData();

            std::vector< trace_point > rows;

            int cell_i = 0;
            for ( auto &cl : m_msh ) {
                if ( m_bnd.cell_has_contact_faces( cl ) ) {
                    const auto di = m_degree_infos.cellDegreeInfo( m_msh, cl );
                    const vector_type uTF = depl.at( cell_i );

                    matrix_type gr;
                    if ( m_rp.m_precomputation )
                        gr = m_data.m_gradient_precomputed.at( cell_i );
                    else
                        gr = make_matrix_hho_symmetric_gradrec( m_msh, cl, m_degree_infos ).first;

                    auto cc = contact_contribution< mesh_type >( m_msh, mat, m_rp, m_bnd );
                    const auto tr = cc.contact_boundary_trace( cl, di, gr, uTF );
                    rows.insert( rows.end(), tr.begin(), tr.end() );
                }
                cell_i++;
            }

            // order along the boundary
            std::sort( rows.begin(), rows.end(), []( const trace_point &a, const trace_point &b ) {
                if ( a.pt.x() != b.pt.x() )
                    return a.pt.x() < b.pt.x();
                return a.pt.y() < b.pt.y();
            } );

            const scalar_type nan = std::numeric_limits< scalar_type >::quiet_NaN();

            std::ofstream ofs( filename );
            ofs << std::setprecision( 12 );
            ofs << "x,y,x_def,y_def,nx,ny,tx,ty,w,gap0,gap,u_n,u_t,ux,uy,"
                << "sigma_nn,sigma_nt,abs_sigma_nt,Fc,coulomb_limit,stress_ratio,"
                << "phi_n,phi_t,fric_bound,nitsche_ratio,state\n";

            for ( const auto &r : rows ) {
                const scalar_type tx = -r.n( 1 );
                const scalar_type ty = r.n( 0 );
                const scalar_type gap = r.gap0 - r.u_n;

                const scalar_type coulomb_limit = r.Fc * std::max( scalar_type( 0 ), -r.sigma_nn );
                const scalar_type stress_ratio =
                    coulomb_limit > scalar_type( 0 ) ? std::abs( r.sigma_nt ) / coulomb_limit : nan;
                const scalar_type nitsche_ratio =
                    r.fric_bound > scalar_type( 0 ) ? std::abs( r.phi_t ) / r.fric_bound : nan;

                std::string state;
                if ( r.phi_n >= scalar_type( 0 ) )
                    state = "SEP";
                else if ( std::abs( r.phi_t ) > r.fric_bound )
                    state = "SLIP";
                else
                    state = "STICK";

                ofs << r.pt.x() << "," << r.pt.y() << "," << r.pt.x() + r.u( 0 ) << ","
                    << r.pt.y() + r.u( 1 ) << "," << r.n( 0 ) << "," << r.n( 1 ) << "," << tx << ","
                    << ty << "," << r.weight << "," << r.gap0 << "," << gap << "," << r.u_n << ","
                    << r.u_t << "," << r.u( 0 ) << "," << r.u( 1 ) << "," << r.sigma_nn << ","
                    << r.sigma_nt << "," << std::abs( r.sigma_nt ) << "," << r.Fc << ","
                    << coulomb_limit << "," << stress_ratio << "," << r.phi_n << "," << r.phi_t
                    << "," << r.fric_bound << "," << nitsche_ratio << "," << state << "\n";
            }
            ofs.close();

            if ( m_verbose )
                std::cout << "** contact boundary trace: " << rows.size() << " points -> "
                          << filename << std::endl;
        }
    }

    void
    output_discontinuous_field( const std::string &filename, FieldName name ) const {
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
        gmsh::NodeData nodedata( 3, 0.0, "discontinuous_nodes", data, subdata );

        // Save the view
        nodedata.saveNodeData( filename, gmsh );
    }

    void
    output_continuous_field( const std::string &quantity, const FieldName name ) {
        constexpr std::size_t dimension = mesh_type::dimension;

        using nodal_vector_type = static_vector< scalar_type, dimension >;

        const nodal_vector_type zero = nodal_vector_type::Zero();

        /*
         * Avoid copying all cell fields if getCurrentField()
         * returns a persistent container.
         */
        const auto &field_cells = m_fields.getCurrentField( name );

        const std::size_t nb_nodes = m_msh.points_size();

        /*
         * Build a point table indexed by the DiSk++ point identifier.
         *
         * This follows the current DiSk++ assumption that point identifiers
         * correspond to the zero-based order returned by points_begin().
         */
        std::vector< typename mesh_type::point_type > point_of_id;

        point_of_id.reserve( nb_nodes );

        for ( auto it = m_msh.points_begin(); it != m_msh.points_end(); ++it ) {
            point_of_id.push_back( *it );
        }

        /*
         * Number of cell contributions at every node.
         */
        std::vector< std::size_t > contribution_count( nb_nodes, 0 );

        /*
         * Accumulated nodal value.
         */
        std::vector< nodal_vector_type > nodal_values( nb_nodes, zero );

        for ( const auto &cl : m_msh ) {
            const std::size_t cell_id = m_msh.lookup( cl );

            const auto degree_info = m_degree_infos.cellDegreeInfo( m_msh, cl );

            const auto cell_basis =
                make_vector_monomial_basis( m_msh, cl, degree_info.cell_degree() );

            const auto &cell_field = field_cells.at( cell_id );

            /*
             * Use point identifiers directly. This avoids looking up a point
             * identifier from its floating-point coordinates.
             */
            for ( const auto point_identifier : cl.point_ids() ) {
                const std::size_t point_id = static_cast< std::size_t >( point_identifier );

                if ( point_id >= nb_nodes ) {
                    throw std::runtime_error( "Invalid point identifier while "
                                              "reconstructing nodal field '" +
                                              quantity + "'" );
                }

                const auto &point = point_of_id.at( point_id );

                const auto basis_values = cell_basis.eval_functions( point );

                const nodal_vector_type field_value = eval( cell_field, basis_values );

                nodal_values[point_id] += field_value;

                ++contribution_count[point_id];
            }
        }

        /*
         * Average the cell reconstructions at shared nodes.
         */
        for ( std::size_t point_id = 0; point_id < nb_nodes; ++point_id ) {
            const std::size_t count = contribution_count[point_id];

            if ( count == 0 ) {
                throw std::runtime_error( "Node " + std::to_string( point_id ) +
                                          " has no cell contribution for field '" + quantity +
                                          "'" );
            }

            nodal_values[point_id] /= static_cast< scalar_type >( count );
        }

        /*
         * The EnSight writer accepts Eigen-compatible static vectors.
         *
         * In 2D, the writer automatically generates:
         *
         * ux, uy, 0
         *
         * In 3D:
         *
         * ux, uy, uz
         */
        m_output.write_vector( quantity, nodal_values );
    }

    void
    output_GaussDisplacement_GP() {
        constexpr std::size_t dimension = mesh_type::dimension;
        using displacement_type = static_vector< scalar_type, dimension >;

        const auto &displacement = m_fields.getCurrentField( FieldName::DEPL_CELLS );

        std::vector< displacement_type > values;
        values.reserve( m_output.point_cloud( m_gauss_point_part_id ).number_of_points() );

        for ( const auto &cl : m_msh ) {
            const std::size_t cell_id = m_msh.lookup( cl );
            const auto degree_info = m_degree_infos.cellDegreeInfo( m_msh, cl );

            const auto cell_basis =
                make_vector_monomial_basis( m_msh, cl, degree_info.cell_degree() );

            const auto &cell_unknowns = displacement.at( cell_id );

            const std::size_t number_of_qp = m_behavior.numberOfQP( cell_id );

            for ( std::size_t qp_id = 0; qp_id < number_of_qp; ++qp_id ) {
                const auto qp = m_behavior.quadrature_point( cell_id, qp_id );

                const auto basis_values = cell_basis.eval_functions( qp.point() );

                values.push_back( eval( cell_unknowns, basis_values ) );
            }
        }

        m_output.write_point_cloud_vector( "displacement_GP", m_gauss_point_part_id, values );
    }

    void
    output_CauchyStress_GP() {
        const auto &displacement = m_fields.getCurrentField( FieldName::DEPL );

        std::vector< static_matrix< scalar_type, 3, 3 > > stresses;
        stresses.reserve( m_output.point_cloud( m_gauss_point_part_id ).number_of_points() );

        for ( const auto &cl : m_msh ) {
            const std::size_t cell_id = m_msh.lookup( cl );
            const auto degree_info = m_degree_infos.cellDegreeInfo( m_msh, cl );

            const auto &local_displacement = displacement.at( cell_id );

            matrix_type gradient_operator;
            if ( m_rp.m_precomputation )
                gradient_operator = m_data.m_gradient_precomputed.at( cell_id );
            else if ( m_behavior.getDeformation() == SMALL_DEF )
                gradient_operator =
                    make_matrix_hho_symmetric_gradrec( m_msh, cl, m_degree_infos ).first;
            else
                gradient_operator = make_matrix_hho_gradrec( m_msh, cl, m_degree_infos ).first;

            const vector_type reconstructed_gradient = gradient_operator * local_displacement;

            const auto gradient_basis =
                make_matrix_monomial_basis( m_msh, cl, degree_info.grad_degree() );

            const std::size_t number_of_qp = m_behavior.numberOfQP( cell_id );

            for ( std::size_t qp_id = 0; qp_id < number_of_qp; ++qp_id ) {
                static_matrix< scalar_type, 3, 3 > cauchy_stress;
                cauchy_stress.setZero();

                if ( m_behavior.getDeformation() == SMALL_DEF ) {
                    cauchy_stress = m_behavior.compute_stress3D( cell_id, qp_id );
                } else {
                    const auto qp = m_behavior.quadrature_point( cell_id, qp_id );

                    const auto basis_values = gradient_basis.eval_functions( qp.point() );

                    const auto displacement_gradient = eval( reconstructed_gradient, basis_values );

                    const auto deformation_gradient = convertGtoF( displacement_gradient );

                    const auto deformation_gradient_3d =
                        convertMatrix3DwithOne( deformation_gradient );

                    const auto pk1_stress = m_behavior.compute_stress3D( cell_id, qp_id );

                    cauchy_stress = convertPK1toCauchy( pk1_stress, deformation_gradient_3d );
                }

                stresses.push_back( cauchy_stress );
            }
        }

        m_output.write_point_cloud_symmetric_tensor(
            "cauchyStress_GP", m_gauss_point_part_id, stresses );
    }

    void
    output_is_plastic_GP() {
        std::vector< double > values;
        values.reserve( m_output.point_cloud( m_gauss_point_part_id ).number_of_points() );

        for ( const auto &cl : m_msh ) {
            const std::size_t cell_id = m_msh.lookup( cl );
            const std::size_t number_of_qp = m_behavior.numberOfQP( cell_id );

            for ( std::size_t qp_id = 0; qp_id < number_of_qp; ++qp_id )
                values.push_back( m_behavior.is_plastic( cell_id, qp_id ) ? 1.0 : 0.0 );
        }

        m_output.write_point_cloud_scalar( "state_GP", m_gauss_point_part_id, values );
    }

    void
    output_equivalentPlasticStrain_GP() {
        std::vector< double > values;
        values.reserve( m_output.point_cloud( m_gauss_point_part_id ).number_of_points() );

        for ( const auto &cl : m_msh ) {
            const std::size_t cell_id = m_msh.lookup( cl );
            const std::size_t number_of_qp = m_behavior.numberOfQP( cell_id );

            for ( std::size_t qp_id = 0; qp_id < number_of_qp; ++qp_id )
                values.push_back(
                    static_cast< double >( m_behavior.equivalentPlasticStrain( cell_id, qp_id ) ) );
        }

        m_output.write_point_cloud_scalar(
            "equivalentPlasticStrain_GP", m_gauss_point_part_id, values );
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

    void
    output_stabCoeff() {
        std::vector< double > values;
        values.reserve( m_msh.cells_size() );

        for ( const auto &cl : m_msh )
            values.push_back( static_cast< double >( m_stab_manager.getValue( m_msh, cl ) ) );

        m_output.write_element_scalar( "stab_coeff", values );
    }

    void
    output_normal_stress_boundary_nodes( const bool contact_only = false ) {
        constexpr std::size_t dimension = mesh_type::dimension;

        static_assert( dimension == 2 || dimension == 3,
                       "Normal stress output is implemented only in 2D and 3D" );

        const std::size_t nb_nodes = m_msh.points_size();

        /*
         * Build a table giving direct access to a point from its
         * DiSk++ point identifier.
         *
         * This assumes that point identifiers correspond to the zero-based
         * order returned by points_begin() and points_end().
         */
        std::vector< typename mesh_type::point_type > point_of_id;

        point_of_id.reserve( nb_nodes );

        for ( auto point_it = m_msh.points_begin(); point_it != m_msh.points_end(); ++point_it ) {
            point_of_id.push_back( *point_it );
        }

        /*
         * Mark the boundary faces on which the normal stress
         * must be evaluated.
         *
         * If contact_only is true, only contact faces are selected.
         * Otherwise, all boundary faces are selected.
         */
        std::vector< bool > is_selected_face( m_msh.faces_size(), false );

        for ( auto face_it = m_msh.boundary_faces_begin(); face_it != m_msh.boundary_faces_end();
              ++face_it ) {
            const auto fc = *face_it;

            const std::size_t face_id = m_msh.lookup( fc );

            if ( contact_only ) {
                is_selected_face.at( face_id ) = m_bnd.is_contact_face( face_id );
            } else {
                is_selected_face.at( face_id ) = true;
            }
        }

        /*
         * Accumulated normal stress and number of contributions
         * at every mesh node.
         *
         * Interior nodes and nodes outside the selected boundary
         * remain equal to zero.
         */
        std::vector< scalar_type > sigma_nn_sum( nb_nodes, scalar_type { 0 } );

        std::vector< std::size_t > sigma_nn_count( nb_nodes, std::size_t { 0 } );

        const bool small_deformation = m_behavior.getDeformation() == SMALL_DEF;

        const auto &displacement = m_fields.getCurrentField( FieldName::DEPL );

        /*
         * Loop over the cells.
         *
         * Processing boundary faces from their adjacent cell gives direct
         * access to the cell polynomial used to evaluate the stress.
         */
        for ( const auto &cl : m_msh ) {
            const std::size_t cell_id = m_msh.lookup( cl );

            const auto degree_info = m_degree_infos.cellDegreeInfo( m_msh, cl );

            /*
             * Polynomial projection of the stress inside the cell.
             */
            const auto projected_stress =
                m_behavior.projectStressOnCell( m_msh, cl, degree_info.grad_degree() );

            const auto matrix_basis =
                make_matrix_monomial_basis( m_msh, cl, degree_info.grad_degree() );

            /*
             * Reconstructed displacement gradient.
             *
             * It is needed only for finite-deformation computations.
             */
            vector_type reconstructed_gradient;

            if ( !small_deformation ) {
                matrix_type gradient_operator;

                if ( m_rp.m_precomputation ) {
                    gradient_operator = m_data.m_gradient_precomputed.at( cell_id );
                } else {
                    gradient_operator = make_matrix_hho_gradrec( m_msh, cl, m_degree_infos ).first;
                }

                const auto &local_displacement = displacement.at( cell_id );

                reconstructed_gradient = gradient_operator * local_displacement;
            }

            /*
             * Loop over the faces of the current cell.
             */
            for ( const auto &fc : faces( m_msh, cl ) ) {
                const std::size_t face_id = m_msh.lookup( fc );

                /*
                 * Skip internal faces and boundary faces that do not
                 * satisfy the requested selection.
                 */
                if ( !is_selected_face.at( face_id ) ) {
                    continue;
                }

                /*
                 * Compute the unit normal oriented with respect
                 * to the current cell.
                 */
                auto normal_vector = normal( m_msh, cl, fc );

                const scalar_type normal_norm = normal_vector.norm();

                if ( normal_norm <= std::numeric_limits< scalar_type >::epsilon() ) {
                    continue;
                }

                normal_vector /= normal_norm;

                /*
                 * Evaluate the normal stress at each original vertex
                 * of the boundary face.
                 *
                 * No PostMesh node or artificial face barycentre is used.
                 */
                for ( const auto point_identifier : fc.point_ids() ) {
                    const std::size_t node_id = static_cast< std::size_t >( point_identifier );

                    if ( node_id >= nb_nodes ) {
                        throw std::out_of_range( "Invalid mesh node identifier in "
                                                 "output_normal_stress_boundary_nodes" );
                    }

                    const auto &point = point_of_id.at( node_id );

                    /*
                     * Evaluate the projected stress tensor at the node.
                     */
                    const auto stress_basis_values = matrix_basis.eval_functions( point );

                    const auto stress_tensor = eval( projected_stress, stress_basis_values );

                    scalar_type sigma_nn = scalar_type { 0 };

                    if ( small_deformation ) {
                        /*
                         * In small deformation, stress_tensor is directly
                         * interpreted as the Cauchy stress tensor.
                         */
                        sigma_nn = normal_vector.dot( stress_tensor * normal_vector );
                    } else {
                        /*
                         * Evaluate the reconstructed displacement gradient
                         * at the same node.
                         */
                        const auto displacement_gradient =
                            eval( reconstructed_gradient, stress_basis_values );

                        const auto deformation_gradient = convertGtoF( displacement_gradient );

                        if constexpr ( dimension == 3 ) {
                            /*
                             * Three-dimensional conversion:
                             *
                             * sigma = (1 / det(F)) P F^T
                             */
                            const auto cauchy_stress =
                                convertPK1toCauchy( stress_tensor, deformation_gradient );

                            sigma_nn = normal_vector.dot( cauchy_stress * normal_vector );
                        } else {
                            /*
                             * In two dimensions, embed the stress and
                             * deformation-gradient tensors in three dimensions
                             * before applying the PK1-to-Cauchy conversion.
                             */
                            static_matrix< scalar_type, 3, 3 > stress_tensor_3d;

                            stress_tensor_3d.setZero();

                            stress_tensor_3d( 0, 0 ) = stress_tensor( 0, 0 );

                            stress_tensor_3d( 0, 1 ) = stress_tensor( 0, 1 );

                            stress_tensor_3d( 1, 0 ) = stress_tensor( 1, 0 );

                            stress_tensor_3d( 1, 1 ) = stress_tensor( 1, 1 );

                            const auto deformation_gradient_3d =
                                convertMatrix3DwithOne( deformation_gradient );

                            const auto cauchy_stress_3d =
                                convertPK1toCauchy( stress_tensor_3d, deformation_gradient_3d );

                            static_vector< scalar_type, 3 > normal_vector_3d;

                            normal_vector_3d.setZero();

                            normal_vector_3d( 0 ) = normal_vector( 0 );

                            normal_vector_3d( 1 ) = normal_vector( 1 );

                            sigma_nn = normal_vector_3d.dot( cauchy_stress_3d * normal_vector_3d );
                        }
                    }

                    sigma_nn_sum.at( node_id ) += sigma_nn;

                    ++sigma_nn_count.at( node_id );
                }
            }
        }

        /*
         * Compute the arithmetic average of the contributions
         * at each boundary node.
         *
         * Interior nodes and unselected boundary nodes remain zero.
         */
        std::vector< double > nodal_values( nb_nodes, 0.0 );

        for ( std::size_t node_id = 0; node_id < nb_nodes; ++node_id ) {
            const std::size_t count = sigma_nn_count[node_id];

            if ( count == 0 ) {
                continue;
            }

            nodal_values[node_id] = static_cast< double >( sigma_nn_sum[node_id] /
                                                           static_cast< scalar_type >( count ) );
        }

        /*
         * The current step and physical time are provided by
         * EnsightExporter::begin_step(time).
         */
        if ( contact_only ) {
            m_output.write_scalar( "contact_stress", nodal_values );
        } else {
            m_output.write_scalar( "normal_stress", nodal_values );
        }
    }
};
} // namespace mechanics

} // namespace disk
