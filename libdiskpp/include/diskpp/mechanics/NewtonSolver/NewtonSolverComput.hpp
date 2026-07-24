/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019, 2024               nicolas.pignet@enpc.fr
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
#include "diskpp/mechanics/NewtonSolver/NewtonSolverContact.hpp"
#include "diskpp/mechanics/NewtonSolver/StabilizationManager.hpp"
#include "diskpp/mechanics/NewtonSolver/TimeManager.hpp"
#include "diskpp/mechanics/behaviors/laws/behaviorlaws.hpp"
#include "diskpp/mechanics/deformation_tensors.hpp"
#include "diskpp/methods/hho"
#include "diskpp/quadratures/quadratures.hpp"

#include <cassert>

namespace disk {

namespace mechanics {

template < typename mesh_type, typename T >
dynamic_matrix< T > _gradrec( const mesh_type &msh, const typename mesh_type::cell &cl,
                              const NonLinearParameters< T > &rp,
                              const MeshDegreeInfo< mesh_type > &degree_infos, bool small_def,
                              const std::vector< dynamic_matrix< T > > &gradient_precomputed ) {
    if ( rp.m_precomputation ) {
        const auto cell_i = msh.lookup( cl );
        return gradient_precomputed[cell_i];
    }

    if ( small_def ) {
        const auto gradrec_sym_full = make_matrix_hho_symmetric_gradrec( msh, cl, degree_infos );
        return gradrec_sym_full.first;
    } else {
        const auto gradrec_full = make_matrix_hho_gradrec( msh, cl, degree_infos );
        return gradrec_full.first;
    }

    return dynamic_matrix< T >();
}

template < typename mesh_type, typename T >
dynamic_matrix< T > _stab( const mesh_type &msh, const typename mesh_type::cell &cl,
                           const NonLinearParameters< T > &rp,
                           const MeshDegreeInfo< mesh_type > &degree_infos,
                           const std::vector< dynamic_matrix< T > > &stab_precomputed ) {
    if ( rp.m_precomputation ) {
        const auto cell_i = msh.lookup( cl );
        return stab_precomputed.at( cell_i );
    } else {
        switch ( rp.m_stab_type ) {
        case StabilizationType::HHO_SYM: {
            const auto recons = make_vector_hho_symmetric_laplacian( msh, cl, degree_infos );
            return make_vector_hho_stabilization( msh, cl, recons.first, degree_infos );
            break;
        }
        case StabilizationType::HHO: {
            const auto recons_scalar = make_scalar_hho_laplacian( msh, cl, degree_infos );
            return make_vector_hho_stabilization_optim( msh, cl, recons_scalar.first,
                                                        degree_infos );
            break;
        }
        case StabilizationType::HDG: {
            return make_vector_hdg_stabilization( msh, cl, degree_infos );
            break;
        }
        case StabilizationType::DG: {
            return make_vector_dg_stabilization( msh, cl, degree_infos );
            break;
        }
        case StabilizationType::NO: {
            break;
        }
        default:
            throw std::invalid_argument( "Unknown stabilization" );
        }
    }

    return dynamic_matrix< T >();
}

template < typename MeshType >
class mechanical_computation {
    typedef MeshType mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell cell_type;

    typedef NonLinearParameters< scalar_type > param_type;
    typedef Behavior< mesh_type > behavior_type;
    typedef vector_boundary_conditions< mesh_type > bnd_type;

    typedef static_matrix< scalar_type, mesh_type::dimension, mesh_type::dimension >
        static_matrix_type;
    typedef static_tensor< scalar_type, mesh_type::dimension > static_tensor_type;

    const static int dimension = mesh_type::dimension;

    typedef dynamic_matrix< scalar_type > matrix_type;
    typedef dynamic_vector< scalar_type > vector_type;

    typedef std::function< static_vector< scalar_type, mesh_type::dimension >(
        const point< scalar_type, mesh_type::dimension > & ) >
        func_type;

    bool two_dim;

    /**
     * @brief Compute A : gphi for finite deformation
     *
     * @param A fourth-order tangent modulus
     * @param gphi set of basis function for gradient reconstruction
     * @return eigen_compatible_stdvector<static_matrix_type> A : gphi
     */
    eigen_compatible_stdvector< static_matrix_type >
    compute_A_gphi( const static_tensor_type &A,
                    const eigen_compatible_stdvector< static_matrix_type > &gphi ) const {
        const auto grad_basis_size = gphi.size();
        const auto DIM2 = dimension * dimension;

        disk::eigen_compatible_stdvector< static_matrix_type > Aphi;
        Aphi.reserve( grad_basis_size );

        // poly classique
        for ( size_t i = 0; i < grad_basis_size; i += DIM2 ) {
            size_t row = i;
            for ( size_t k = 0; k < dimension; k++ ) {     // depend de l'ordre des bases
                for ( size_t l = 0; l < dimension; l++ ) { // depend de l'ordre des bases
                    Aphi.push_back( disk::tm_prod( A, gphi[row], l, k ) );
                    row++;
                }
            }
        }

        return Aphi;
    }

    /**
     * @brief Compute the stress tensor and the tangent modulus at the given qudature point
     *  - in small deformation Cauchy and Cep
     *  - in finite deformation PK1 and A
     * @tparam QPType Type of law at the quadrature point
     * @tparam LawData Type for material data
     * @param qp A behavior quadrature point where we compute the constitutive law
     * @param RkT_iqn Reconstruction operator evaluated at the quadrature point (symmetric gradient)
     * in small deformation and gradient in finite deformation)
     * @param small_def small deformation yes or no
     * @return std::pair<static_matrix_type, static_tensor_type> stress tensor and the tangent
     * modulus
     */

    std::pair< static_matrix_type, static_tensor_type >
    _compute_behavior( behavior_type &behavior, size_t cell_id, size_t qp_id,
                       const static_matrix_type &RkT_iqn, const bool small_def,
                       const bool use_tangent_modulus, const bool compute_tangent ) const {
        if ( compute_tangent ) {
            if ( small_def ) {
                return behavior.compute_whole( cell_id, qp_id, RkT_iqn, use_tangent_modulus );
            } else {
                const auto FT_iqn = convertGtoF( RkT_iqn );
                return behavior.compute_whole( cell_id, qp_id, FT_iqn, use_tangent_modulus );
            }
        } else {
            if ( small_def ) {
                return std::make_pair( behavior.compute_stress( cell_id, qp_id, RkT_iqn ),
                                       static_tensor_type() );
            } else {
                const auto FT_iqn = convertGtoF( RkT_iqn );
                return std::make_pair( behavior.compute_stress( cell_id, qp_id, FT_iqn ),
                                       static_tensor_type() );
            }
        }
    }

    void _compute_rigidity( const static_tensor_type &Cep,
                            const eigen_compatible_stdvector< static_matrix_type > &gphi,
                            bool small_def, const scalar_type weight, const size_t grad_dim_dofs,
                            const size_t grad_basis_size, matrix_type &AT ) {
        // std::cout << "module : " << Cep << std::endl;
        timecounter tc;
        tc.tic();

        if ( small_def ) {
            // upper part
            for ( size_t j = 0; j < grad_basis_size; j++ ) {
                const static_matrix_type Cgphi_j = weight * tm_prod( Cep, gphi[j] );
                // std::cout << j << std::endl;
                // std::cout << gphi[j] << std::endl;
                // std::cout << Cgphi_j << std::endl;
                for ( size_t i = 0; i <= j; i += grad_dim_dofs ) {
                    // compute (Ekt v, C(u) : Ekt du)
                    if ( two_dim ) {
                        AT( i, j ) += Cgphi_j( 0, 0 ) * gphi[i]( 0, 0 );
                        AT( i + 1, j ) += 2 * Cgphi_j( 0, 1 ) * gphi[i + 1]( 0, 1 );
                        AT( i + 2, j ) += Cgphi_j( 1, 1 ) * gphi[i + 2]( 1, 1 );
                    } else {
                        AT( i, j ) += Cgphi_j( 0, 0 ) * gphi[i]( 0, 0 );
                        AT( i + 1, j ) += 2 * Cgphi_j( 0, 1 ) * gphi[i + 1]( 0, 1 );
                        AT( i + 2, j ) += Cgphi_j( 1, 1 ) * gphi[i + 2]( 1, 1 );
                        AT( i + 3, j ) += 2 * Cgphi_j( 0, 2 ) * gphi[i + 3]( 0, 2 );
                        AT( i + 4, j ) += 2 * Cgphi_j( 1, 2 ) * gphi[i + 4]( 1, 2 );
                        AT( i + 5, j ) += Cgphi_j( 2, 2 ) * gphi[i + 5]( 2, 2 );
                    }
                    // AT(i, j) += disk::mm_prod(gphi[i], Agphi_j);
                }
            }
        } else {
            // lower part
            const auto qp_A_gphi = compute_A_gphi( weight * Cep, gphi );

            for ( size_t j = 0; j < grad_basis_size; j += grad_dim_dofs ) {
                size_t col = j;
                for ( size_t k = 0; k < dimension; k++ ) {     // depend de l'ordre des bases
                    for ( size_t l = 0; l < dimension; l++ ) { // depend de l'ordre des bases
                        for ( size_t i = col; i < grad_basis_size; i++ ) {
                            AT( i, col ) += qp_A_gphi[i]( l, k ) * gphi[col]( l, k );
                        }
                        col++;
                    }
                }
            }
        }
        tc.toc();
        time_rigi += tc.elapsed();
    }

    void _symmetrized_rigidity_matrix( const size_t grad_basis_size, matrix_type &AT,
                                       const bool small_def ) {
        timecounter tc;
        tc.tic();

        if ( small_def ) {
            // lower part AT
            for ( size_t j = 0; j < grad_basis_size; j++ )
                for ( size_t i = j; i < grad_basis_size; i++ )
                    AT( i, j ) = AT( j, i );
        } else {
            // upper part AT
            for ( size_t i = 0; i < grad_basis_size; i++ )
                for ( size_t j = i; j < grad_basis_size; j++ )
                    AT( i, j ) = AT( j, i );
        }
        tc.toc();
        time_rigi += tc.elapsed();
    }

    void compute_external_forces( const mesh_type &msh, const cell_type &cl,
                                  const std::unique_ptr< func_type > &load,
                                  const size_t cell_degree ) {
        timecounter tc;
        tc.tic();

        if ( load ) {
            // compute (f,v)_T

            const auto cb = make_vector_monomial_basis( msh, cl, cell_degree );
            RTF.head( cb.size() ) += make_rhs( msh, cl, cb, *load, 1 );
        }
        tc.toc();
        time_load += tc.elapsed();
    }

    void compute_internal_forces( const static_matrix_type &stress,
                                  const eigen_compatible_stdvector< static_matrix_type > &gphi,
                                  bool small_def, const scalar_type weight,
                                  const size_t grad_dim_dofs, const size_t grad_basis_size,
                                  vector_type &aT ) {
        timecounter tc;
        tc.tic();

        //   std::cout << "stress" << std::endl;
        //   std::cout << stress << std::endl;
        const static_matrix_type stress_qp = weight * stress;

        if ( small_def ) {
            // compute (sigma(u), E^k_T v)_T
            for ( size_t i = 0; i < grad_basis_size; i += grad_dim_dofs ) {
                if ( two_dim ) {
                    aT( i ) += stress_qp( 0, 0 ) * gphi[i]( 0, 0 );
                    aT( i + 1 ) += 2 * stress_qp( 0, 1 ) * gphi[i + 1]( 0, 1 );
                    aT( i + 2 ) += stress_qp( 1, 1 ) * gphi[i + 2]( 1, 1 );
                } else {
                    aT( i ) += stress_qp( 0, 0 ) * gphi[i]( 0, 0 );
                    aT( i + 1 ) += 2 * stress_qp( 0, 1 ) * gphi[i + 1]( 0, 1 );
                    aT( i + 2 ) += stress_qp( 1, 1 ) * gphi[i + 2]( 1, 1 );
                    aT( i + 3 ) += 2 * stress_qp( 0, 2 ) * gphi[i + 3]( 0, 2 );
                    aT( i + 4 ) += 2 * stress_qp( 1, 2 ) * gphi[i + 4]( 1, 2 );
                    aT( i + 5 ) += stress_qp( 2, 2 ) * gphi[i + 5]( 2, 2 );
                }
            }
        } else {
            // compute (PK1(u), G^k_T v)_T
            for ( size_t i = 0; i < grad_basis_size; i += grad_dim_dofs ) {
                size_t row = i;
                for ( size_t k = 0; k < dimension; k++ ) {     // depend de l'ordre des bases
                    for ( size_t l = 0; l < dimension; l++ ) { // depend de l'ordre des bases
                        // compute (PK1(u), G^k_T v)_T
                        aT( row ) += stress_qp( l, k ) * gphi[row]( l, k );
                        row++;
                    }
                }
            }
        }
        tc.toc();
        time_fint += tc.elapsed();
    }

    void compute_contact_terms( const mesh_type &msh, const cell_type &cl, const bnd_type &bnd,
                                const param_type &rp, const CellDegreeInfo< mesh_type > &cell_infos,
                                const matrix_type &RkT, const vector_type &uTF,
                                const TimeStep< scalar_type > &time_step,
                                behavior_type &behavior ) {
        if ( bnd.cell_has_contact_faces( cl ) ) {
            const auto &material_data = behavior.getMaterialData();
            auto cc = contact_contribution( msh, material_data, rp, bnd );
            cc.compute( cl, cell_infos, RkT, uTF );

            time_contact += cc.time_contact;

            K_int += cc.K_cont;
            F_int += cc.F_cont;
            RTF -= cc.F_cont;
        }
    }

    size_t num_grad_dim_dofs( bool small_def ) const {
        if ( small_def ) {
            if ( two_dim )
                return 3;
            else
                return 6;
        } else {
            if ( two_dim )
                return 4;
            else
                return 9;
        }
    }

  public:
    matrix_type K_int;
    vector_type RTF;
    vector_type F_int;
    double time_law, time_load, time_fint, time_rigi;
    double time_contact;

    mechanical_computation( void )
        : time_law( 0. ), time_load( 0. ), time_contact( 0. ), time_fint( 0. ), time_rigi( 0. ) {
        if ( dimension == 2 )
            two_dim = true;
        else if ( dimension == 3 )
            two_dim = false;
        else
            assert( false );
    }

    void compute( const mesh_type &msh, const cell_type &cl, const bnd_type &bnd,
                  const param_type &rp, const MeshDegreeInfo< mesh_type > &degree_infos,
                  const std::unique_ptr< func_type > &load, const matrix_type &RkT,
                  const vector_type &uTF, const TimeStep< scalar_type > &time_step,
                  behavior_type &behavior, StabCoeffManager< scalar_type > &stab_manager,
                  const bool small_def, const bool tangent_matix = true,
                  const bool use_tangent_modulus = true ) {
        timecounter tc;

        const auto cell_infos = degree_infos.cellDegreeInfo( msh, cl );
        const auto faces_infos = cell_infos.facesDegreeInfo();

        const auto cell_degree = cell_infos.cell_degree();
        const auto grad_degree = cell_infos.grad_degree();

        const auto cell_basis_size = vector_basis_size( cell_degree, dimension, dimension );

        size_t gb_size = 0;
        if ( small_def ) {
            gb_size = sym_matrix_basis_size( grad_degree, dimension, dimension );
        } else {
            gb_size = matrix_basis_size( grad_degree, dimension, dimension );
        }

        const auto grad_basis_size = gb_size;

        const auto num_faces_dofs = vector_faces_dofs( msh, faces_infos );
        const auto num_total_dofs = cell_basis_size + num_faces_dofs;
        const auto grad_dim_dofs = num_grad_dim_dofs( small_def );

        matrix_type AT = matrix_type::Zero( grad_basis_size, grad_basis_size );
        vector_type aT = vector_type::Zero( grad_basis_size );

        RTF = vector_type::Zero( num_total_dofs );
        F_int = vector_type::Zero( num_total_dofs );

        assert( RkT.cols() == uTF.rows() );
        assert( RkT.rows() == grad_basis_size );

        // std::cout << "sol" << std::endl;
        // std::cout << uTF.transpose() << std::endl;

        const vector_type RkT_uTF = RkT * uTF;

        // std::cout << "RkT: " << RkT.norm() << std::endl;
        // std::cout << "RkT_Utf: " << RkT_uTF.transpose() << std::endl;

        const auto gb = make_matrix_monomial_basis( msh, cl, grad_degree );
        const auto gbs = make_sym_matrix_monomial_basis( msh, cl, grad_degree );

        eigen_compatible_stdvector< static_matrix_type > gphi;

        const auto cell_id = msh.lookup( cl );
        const auto nb_qp = behavior.numberOfQP( cell_id );

        scalar_type beta_comp = 0.0, total_weight = 0.0;
        for ( int i_qp = 0; i_qp < nb_qp; i_qp++ ) {
            // Compute gradient basis function
            const auto qp = behavior.quadrature_point( cell_id, i_qp );
            // std::cout << "qp: " << qp.point() << std::endl;

            if ( small_def ) {
                gphi = gbs.eval_functions( qp.point() );
            } else {
                gphi = gb.eval_functions( qp.point() );
            }

            assert( gphi.size() == grad_basis_size );

            // Compute local gradient and norm
            // RkT_iqn = Grad_sym for small def else RkT_iqn = Grad
            const auto RkT_iqn = eval( RkT_uTF, gphi );

            // std::cout << "RkT_iqn" << std::endl;
            // std::cout << RkT_iqn << std::endl;

            // Compute behavior
            // if small_def stress = Cauchy else stress = PK1
            tc.tic();
            const auto [stress, Cep] = _compute_behavior(
                behavior, cell_id, i_qp, RkT_iqn, small_def, use_tangent_modulus, tangent_matix );
            tc.toc();
            time_law += tc.elapsed();
            // std::cout << "stress: " << stress.norm() << std::endl;
            // std::cout << stress << std::endl;
            // std::cout << "Cep: " << Cep.norm() << std::endl;
            // std::cout << Cep << std::endl;

            if ( tangent_matix ) {
                // Compute rigidity
                this->_compute_rigidity( Cep, gphi, small_def, qp.weight(), grad_dim_dofs,
                                         grad_basis_size, AT );
            }

            // Compute internal force
            this->compute_internal_forces( stress, gphi, small_def, qp.weight(), grad_dim_dofs,
                                           grad_basis_size, aT );

            // compute new possible value for stabilization
            if ( rp.m_adapt_stab ) {
                scalar_type sigma_dev_norm, eps_dev_norm;
                if ( small_def ) {
                    sigma_dev_norm = deviator( stress ).norm();
                    eps_dev_norm = deviator( RkT_iqn ).norm();
                } else {
                    const auto F = convertGtoF( RkT_iqn );
                    const auto EGL = convertFtoGreenLagrange( F );
                    eps_dev_norm = deviator( EGL ).norm();

                    const auto PK2 = convertPK1toPK2( stress, F );
                    sigma_dev_norm = deviator( PK2 ).norm();

                    // const auto Cauchy = convertPK1toCauchy(stress, F);
                    // const auto eps    = convertGtoLinearizedStrain(RkT_iqn);

                    // std::cout << sigma_dev_norm / eps_dev_norm << " vs " <<
                    // deviator(Cauchy).norm() / deviator(eps).norm()
                    //           << std::endl;
                }

                total_weight += qp.weight();
                if ( eps_dev_norm < 1E-12 )
                    beta_comp += qp.weight() * rp.m_beta;
                else
                    beta_comp += qp.weight() * sigma_dev_norm / eps_dev_norm;
            }
        }

        // save new stabilization coeff beta_s in [beta/1000, 1000*beta]
        if ( rp.m_adapt_stab ) {
            beta_comp /= total_weight;
            const auto beta_s =
                std::min( 1000. * rp.m_beta, std::max( rp.m_beta / 10000., beta_comp ) );
            // std::cout << beta_comp << " vs " << beta_s << std::endl;
            stab_manager.setValueNext( msh, cl, beta_s );
        }

        // compute external forces
        this->compute_external_forces( msh, cl, load, cell_degree );
        // std::cout << "R_ext: " << RTF.norm() << std::endl;
        // std::cout << RTF.transpose() << std::endl;

        if ( tangent_matix ) {

            // Symmetrize rigidity matrix
            this->_symmetrized_rigidity_matrix( grad_basis_size, AT, small_def );

            // std::cout << "AT: " << AT.norm() << std::endl;
            // std::cout << AT << std::endl;
            // std::cout << "aT: " << aT.norm() << std::endl;

            tc.tic();
            K_int = RkT.transpose() * AT * RkT;
            tc.toc();
            time_rigi += tc.elapsed();
        }
        tc.tic();
        F_int = RkT.transpose() * aT;
        RTF -= F_int;
        tc.toc();
        time_fint += tc.elapsed();

        // Compute contact terms

        this->compute_contact_terms( msh, cl, bnd, rp, cell_infos, RkT, uTF, time_step, behavior );

        // std::cout << "K: " << K_int.norm() << std::endl;
        // // std::cout << K_int << std::endl;
        // std::cout << "F_int: " << F_int.norm() << std::endl;
        // std::cout << F_int.transpose() << std::endl;
        // std::cout << "RTF: " << RTF.norm() << std::endl;
        // std::cout << RTF.transpose() << std::endl;

        // throw std::runtime_error("");

        assert( K_int.rows() == num_total_dofs );
        assert( K_int.cols() == num_total_dofs );
        assert( RTF.rows() == num_total_dofs );
    }

    void compute_rigidity_matrix( const mesh_type &msh, const cell_type &cl, const bnd_type &bnd,
                                  const param_type &rp,
                                  const MeshDegreeInfo< mesh_type > &degree_infos,
                                  const matrix_type &RkT, const vector_type &uTF,
                                  const TimeStep< scalar_type > &time_step, behavior_type &behavior,
                                  const bool small_def, bool use_tangente = true ) {
        timecounter tc;

        const auto cell_infos = degree_infos.cellDegreeInfo( msh, cl );
        const auto faces_infos = cell_infos.facesDegreeInfo();

        const auto cell_degree = cell_infos.cell_degree();
        const auto grad_degree = cell_infos.grad_degree();

        const auto cell_basis_size = vector_basis_size( cell_degree, dimension, dimension );

        size_t gb_size = 0;
        if ( small_def ) {
            gb_size = sym_matrix_basis_size( grad_degree, dimension, dimension );
        } else {
            gb_size = matrix_basis_size( grad_degree, dimension, dimension );
        }

        const auto grad_basis_size = gb_size;

        const auto num_faces_dofs = vector_faces_dofs( msh, faces_infos );
        const auto num_total_dofs = cell_basis_size + num_faces_dofs;
        const auto grad_dim_dofs = num_grad_dim_dofs( small_def );

        matrix_type AT = matrix_type::Zero( grad_basis_size, grad_basis_size );

        RTF = vector_type::Zero( num_total_dofs );

        assert( RkT.cols() == uTF.rows() );
        assert( RkT.rows() == grad_basis_size );

        // std::cout << "sol" << std::endl;
        // std::cout << uTF.transpose() << std::endl;

        const vector_type RkT_uTF = RkT * uTF;

        // std::cout << "RkT: " << RkT.norm() << std::endl;
        // std::cout << "RkT_Utf: " << RkT_uTF.transpose() << std::endl;

        const auto gb = make_matrix_monomial_basis( msh, cl, grad_degree );
        const auto gbs = make_sym_matrix_monomial_basis( msh, cl, grad_degree );

        eigen_compatible_stdvector< static_matrix_type > gphi;

        const auto cell_id = msh.lookup( cl );
        const auto nb_qp = behavior.numberOfQP( cell_id );

        scalar_type beta_comp = 0.0, total_weight = 0.0;
        for ( int i_qp = 0; i_qp < nb_qp; i_qp++ ) {
            // Compute gradient basis function
            const auto qp = behavior.quadrature_point( cell_id, i_qp );
            // std::cout << "qp: " << qp.point() << std::endl;

            if ( small_def ) {
                gphi = gbs.eval_functions( qp.point() );
            } else {
                gphi = gb.eval_functions( qp.point() );
            }

            assert( gphi.size() == grad_basis_size );

            // Compute local gradient and norm
            // RkT_iqn = Grad_sym for small def else RkT_iqn = Grad
            const auto RkT_iqn = eval( RkT_uTF, gphi );

            // std::cout << "RkT_iqn" << std::endl;
            // std::cout << RkT_iqn << std::endl;

            // Compute behavior
            // if small_def stress = Cauchy else stress = PK1
            tc.tic();
            const auto [stress, Cep] = this->_compute_behavior( behavior, cell_id, i_qp, RkT_iqn,
                                                                small_def, use_tangente, true );
            tc.toc();
            time_law += tc.elapsed();
            // std::cout << "stress: " << stress.norm() << std::endl;
            // std::cout << stress << std::endl;
            // std::cout << "Cep: " << Cep.norm() << std::endl;
            // std::cout << Cep << std::endl;

            // Compute rigidity
            this->_compute_rigidity( Cep, gphi, small_def, qp.weight(), grad_dim_dofs,
                                     grad_basis_size, AT );
        }

        // Symmetrize rigidity matrix
        this->_symmetrized_rigidity_matrix( grad_basis_size, AT, small_def );

        // std::cout << "AT: " << AT.norm() << std::endl;
        // std::cout << AT << std::endl;
        // std::cout << "aT: " << aT.norm() << std::endl;

        tc.tic();
        K_int = RkT.transpose() * AT * RkT;
        tc.toc();
        time_rigi += tc.elapsed();
        // Compute contact terms
        // this->compute_contact_terms( msh, cl, bnd, rp, cell_infos, RkT, uTF, time_step,
        // behavior );

        // std::cout << "K: " << K_int.norm() << std::endl;
        // // std::cout << K_int << std::endl;
        // std::cout << "F_int: " << F_int.norm() << std::endl;
        // std::cout << F_int.transpose() << std::endl;
        // std::cout << "RTF: " << RTF.norm() << std::endl;
        // std::cout << RTF.transpose() << std::endl;

        // throw std::runtime_error("");

        assert( K_int.rows() == num_total_dofs );
        assert( K_int.cols() == num_total_dofs );
    }
};
} // namespace mechanics

} // namespace disk