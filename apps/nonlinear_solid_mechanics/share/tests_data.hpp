/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019, 2024                nicolas.pignet@enpc.fr
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

#include "diskpp/mechanics/NewtonSolver/NonLinearSolver.hpp"

enum STUDY {
    COOK_ELAS,
    COOK_HPP,
    COOK_LARGE,
    COOK_DYNA,
    SPHERE_LARGE,
    TAYLOR_ROD,
    SQUARE_DYNA,
    SQUARE_MATER,
    WAVE_ELAS,
    IMPACT_2D,
    FREE_VIBR_2D,
};

/* Bibliographie */
/*
 * [1] Di Pietro, D. and Ern. A.; A hybrid high-order locking free method for linear elasticity
 * on general meshes; Comput. Methods Appl. Mech. Engrg. 203, pp1-21, (2015).
 *
 * [2] M. Abbas, A. Ern, N. Pignet. Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework; International Journal of Numerical Methods in Engineering
 * (2019) 120(3), 303-327.
 */

/*
 * COOK_ELAS: [1] Section 6.3
 *
 */

template < typename T >
auto getMaterialData( const STUDY &study ) {
    disk::mechanics::MaterialData< T > material_data;

    const T GPa = 1e9;
    const T MPa = 1e6;

    switch ( study ) {
    case STUDY::WAVE_ELAS: {

        const T E = 2.5;
        const T nu = 0.25;

        material_data.setMu( E, nu );
        material_data.setLambda( E, nu );
        material_data.setRho( 1.0 );

        break;
    }
    case STUDY::COOK_ELAS: {
        // Cook Parameters HPP (mm, MPa, kN)

        material_data.setMu( 0.375 );
        material_data.setLambda( 7.5 * 10e6 );

        break;
    }
    case STUDY::COOK_HPP: {
        // Cook Parameters HPP (mm, GPa, kN)

        const T E = 70;
        const T nu = 0.4999;

        material_data.setMu( E, nu );
        material_data.setLambda( E, nu );

        material_data.setK( 0.0 );
        material_data.setH( 0.135 );

        material_data.setSigma_y0( 0.243 );

        material_data.addMfrontParameter( "YoungModulus", material_data.getE() );
        material_data.addMfrontParameter( "PoissonRatio", material_data.getNu() );
        material_data.addMfrontParameter( "HardeningSlope", material_data.getH() );
        material_data.addMfrontParameter( "YieldStrength", material_data.getSigma_y0() );
        break;
    }
    case STUDY::COOK_LARGE: {
        // (mm, GPa, kN)

        const T E = 206.9;
        const T nu = 0.29;

        material_data.setMu( E, nu );
        material_data.setLambda( E, nu );

        material_data.addCurvePoint( 0.0, 0.45 );
        material_data.addCurvePoint( 0.003065, 0.463511 );
        material_data.addCurvePoint( 0.0061270000000000005, 0.476372 );
        material_data.addCurvePoint( 0.009187, 0.488615 );
        material_data.addCurvePoint( 0.012243, 0.500271 );
        material_data.addCurvePoint( 0.015297000000000002, 0.511369 );
        material_data.addCurvePoint( 0.018348, 0.521937 );
        material_data.addCurvePoint( 0.021396000000000002, 0.532001 );
        material_data.addCurvePoint( 0.024443, 0.541585 );
        material_data.addCurvePoint( 0.027487, 0.550714 );
        material_data.addCurvePoint( 0.030528999999999997, 0.55941 );
        material_data.addCurvePoint( 0.033569, 0.567695 );
        material_data.addCurvePoint( 0.036607, 0.575588 );
        material_data.addCurvePoint( 0.039643, 0.58311 );
        material_data.addCurvePoint( 0.042677999999999994, 0.590279 );
        material_data.addCurvePoint( 0.045711, 0.597111 );
        material_data.addCurvePoint( 0.048741999999999994, 0.603625 );
        material_data.addCurvePoint( 0.051772, 0.609835 );
        material_data.addCurvePoint( 0.054801, 0.615757 );
        material_data.addCurvePoint( 0.064887, 0.633592 );
        material_data.addCurvePoint( 0.074961, 0.648851 );
        material_data.addCurvePoint( 0.085024, 0.661934 );
        material_data.addCurvePoint( 0.095079, 0.673181 );
        material_data.addCurvePoint( 0.105126, 0.682878 );
        material_data.addCurvePoint( 0.115166, 0.691265 );
        material_data.addCurvePoint( 0.12520099999999998, 0.698548 );
        material_data.addCurvePoint( 0.135232, 0.704897 );
        material_data.addCurvePoint( 0.145259, 0.710459 );
        material_data.addCurvePoint( 0.155282, 0.715356 );
        material_data.addCurvePoint( 0.16530299999999998, 0.719691 );
        material_data.addCurvePoint( 0.17532199999999998, 0.723553 );
        material_data.addCurvePoint( 0.18533899999999998, 0.727014 );
        material_data.addCurvePoint( 0.195354, 0.730137 );
        material_data.addCurvePoint( 0.205368, 0.732975 );
        material_data.addCurvePoint( 0.21538, 0.735573 );
        material_data.addCurvePoint( 0.22539199999999998, 0.737967 );
        material_data.addCurvePoint( 0.235403, 0.740189 );
        material_data.addCurvePoint( 0.245413, 0.742267 );
        material_data.addCurvePoint( 0.25542200000000004, 0.744222 );
        material_data.addCurvePoint( 0.26543100000000003, 0.746074 );
        material_data.addCurvePoint( 0.27543900000000004, 0.747838 );
        material_data.addCurvePoint( 0.28544800000000004, 0.74953 );
        material_data.addCurvePoint( 0.295456, 0.751158 );
        material_data.addCurvePoint( 0.30546300000000004, 0.752735 );
        material_data.addCurvePoint( 0.401529, 0.766376 );
        material_data.addCurvePoint( 0.501593, 0.779544 );
        material_data.addCurvePoint( 0.6016549999999999, 0.79251 );
        material_data.addCurvePoint( 0.701718, 0.805438 );
        material_data.addCurvePoint( 0.8017799999999999, 0.81836 );
        material_data.addCurvePoint( 0.901843, 0.83128 );
        material_data.addCurvePoint( 1.001905, 0.8442 );

        material_data.checkRpCurve();

        break;
    }
    case STUDY::COOK_DYNA: {
        // Cook Parameters  (mm, GPa, kN, kg, ms)
        // https://www.dynasupport.com/howtos/general/consistent-units

        const T E = 200;
        const T nu = 0.3;

        material_data.setMu( E, nu );
        material_data.setLambda( E, nu );
        material_data.setRho( 7.800e-6 );

        material_data.setK( 0.0 );
        material_data.setH( 0.13 );

        material_data.setSigma_y0( 0.45 );

        material_data.addMfrontParameter( "YoungModulus", material_data.getE() );
        material_data.addMfrontParameter( "PoissonRatio", material_data.getNu() );
        material_data.addMfrontParameter( "HardeningSlope", material_data.getH() );
        material_data.addMfrontParameter( "YieldStrength", material_data.getSigma_y0() );

        break;
    }
    case STUDY::SQUARE_DYNA: {
        // Parameters  (m, Pa, N, kg, s)
        // https://www.dynasupport.com/howtos/general/consistent-units

        material_data.setMu( 1 );
        material_data.setLambda( 1 );
        material_data.setRho( 1 );

        material_data.setK( 0.0 );
        material_data.setH( 0.25 );

        material_data.setSigma_y0( 0.20e9 );

        material_data.addMfrontParameter( "YoungModulus", material_data.getE() );
        material_data.addMfrontParameter( "PoissonRatio", material_data.getNu() );
        material_data.addMfrontParameter( "HardeningSlope", material_data.getH() );
        material_data.addMfrontParameter( "YieldStrength", material_data.getSigma_y0() );

        break;
    }
    case STUDY::SQUARE_MATER: {
        // Parameters  (m, Pa, N, kg, s)
        // https://www.dynasupport.com/howtos/general/consistent-units

        // Steel
        // const T E = 200.0e9;
        // const T nu = 0.3;
        // const T rho = 7800;
        // const T H = 0.13e9;
        // const T Sy0 = 0.45e9;

        // Gold
        const T E = 80.0e9;
        const T nu = 0.42;
        const T rho = 18900;
        const T H = 0.2e9;
        const T Sy0 = 0.02e9;

        material_data.setMu( E, nu );
        material_data.setLambda( E, nu );
        material_data.setRho( rho );
        material_data.setK( 0.0 );
        material_data.setH( H );
        material_data.setSigma_y0( Sy0 );

        material_data.addMfrontParameter( "YoungModulus", material_data.getE() );
        material_data.addMfrontParameter( "PoissonRatio", material_data.getNu() );
        material_data.addMfrontParameter( "HardeningSlope", material_data.getH() );
        material_data.addMfrontParameter( "YieldStrength", material_data.getSigma_y0() );

        break;
    }
    case STUDY::SPHERE_LARGE: {
        // Sphere Parameters (mm, GPa, kN)

        const T E = 28.95;
        const T nu = 0.3;

        material_data.setMu( E, nu );
        material_data.setLambda( E, nu );
        material_data.setK( 0 );
        material_data.setH( 0.0 );
        material_data.setSigma_y0( 6 );
        break;
    }
    case STUDY::TAYLOR_ROD: {
        // (mm, GPa, kN, kg, ms)
        // https://www.dynasupport.com/howtos/general/consistent-units

        const T E = 120;
        const T nu = 0.35;

        material_data.setMu( E, nu );
        material_data.setLambda( E, nu );
        material_data.setRho( 8.930e-6 );

        material_data.setK( 0.0 );
        material_data.setH( 0.1 );

        material_data.setSigma_y0( 0.4 );

        material_data.addMfrontParameter( "YoungModulus", material_data.getE() );
        material_data.addMfrontParameter( "PoissonRatio", material_data.getNu() );
        material_data.addMfrontParameter( "HardeningSlope", material_data.getH() );
        material_data.addMfrontParameter( "YieldStrength", material_data.getSigma_y0() );
        break;
    }
    case STUDY::IMPACT_2D:
    case STUDY::FREE_VIBR_2D: {
        // Parameters  (m, Pa, N, kg, s)
        // https://www.dynasupport.com/howtos/general/consistent-units

        const T E = 1.0;
        const T nu = 0.0;

        material_data.setMu( E, nu );
        material_data.setLambda( E, nu );
        material_data.setRho( 1.0 );

        break;
    }
    default: {
        throw std::invalid_argument( "Unexpected study" );
        break;
    }
    }

    return material_data;
}

template < typename T >
void addAdditionalParameters( const STUDY &study, disk::mechanics::NonLinearParameters< T > &rp ) {

    switch ( study ) {
    case STUDY::COOK_ELAS:
    case STUDY::COOK_HPP:
    case STUDY::COOK_LARGE:
    case STUDY::SPHERE_LARGE: {
        break;
    }
    case STUDY::COOK_DYNA:
    case STUDY::WAVE_ELAS:
    case STUDY::SQUARE_DYNA:
    case STUDY::SQUARE_MATER:
    case STUDY::TAYLOR_ROD:
    case STUDY::IMPACT_2D:
    case STUDY::FREE_VIBR_2D: {
        std::map< std::string, T > dyna_para;
        dyna_para["beta"] = 0.25;
        dyna_para["gamma"] = 0.5;
        dyna_para["theta"] = 1.0;

        rp.setUnsteadyParameters( dyna_para );
        rp.setLinearSolver( disk::solvers::direct_solver::pardiso );

        break;
    }
    default: {
        throw std::invalid_argument( "Unexpected study" );
        break;
    }
    }
}

template < template < typename, size_t, typename > class Mesh, typename T, typename Storage >
auto getBoundaryConditions( const Mesh< T, 2, Storage > &msh,
                            const disk::mechanics::MaterialData< T > &material_data,
                            const STUDY &study ) {
    typedef Mesh< T, 2, Storage > mesh_type;
    typedef disk::static_vector< T, 2 > result_type;

    disk::vector_boundary_conditions< mesh_type > bnd( msh );

    auto zero = [material_data]( const disk::point< T, 2 > &p, const T &time ) -> result_type {
        return result_type { 0.0, 0.0 };
    };

    /* Boundary conditions */
    switch ( study ) {
    case STUDY::WAVE_ELAS: {

        auto func_space = [material_data]( const disk::point< T, 2 > &p ) -> result_type {
            T ux = -sin( M_PI * p.x() ) * cos( M_PI * p.y() );
            T uy = cos( M_PI * p.x() ) * sin( M_PI * p.y() );

            return result_type { ux, uy };
        };

        auto displacement = [material_data, func_space]( const disk::point< T, 2 > &p,
                                                         const T &time ) -> result_type {
            return time * time * func_space( p );
        };

        bnd.addDirichletEverywhere( displacement );

        break;
    }
    case STUDY::COOK_ELAS: {
        auto trac = [material_data]( const disk::point< T, 2 > &p, const T &time ) -> result_type {
            T L = 16.;
            T F = 1.;
            return time * result_type { 0.0, F / L };
        };

        /* Encast */
        bnd.addDirichletBC( disk::CLAMPED, 1, zero );
        /* Load */
        bnd.addNeumannBC( disk::NEUMANN, 2, trac );
        break;
    }
    case STUDY::COOK_HPP: {
        auto trac = [material_data]( const disk::point< T, 2 > &p, const T &time ) -> result_type {
            T L = 16.;
            T F = 1.8;
            return time * result_type { 0.0, F / L };
        };

        /* Encast */
        bnd.addDirichletBC( disk::CLAMPED, 1, zero );
        /* Load */
        bnd.addNeumannBC( disk::NEUMANN, 2, trac );
        break;
    }
    case STUDY::COOK_LARGE: {
        auto trac = [material_data]( const disk::point< T, 2 > &p, const T &time ) -> result_type {
            T L = 16.;
            T F = 5.0;
            return time * result_type { 0.0, F / L };
        };

        /* Encast */
        bnd.addDirichletBC( disk::CLAMPED, 1, zero );
        /* Load */
        bnd.addNeumannBC( disk::NEUMANN, 2, trac );
        break;
    }
    case STUDY::COOK_DYNA: {
        auto trac = [material_data]( const disk::point< T, 2 > &p, const T &time ) -> result_type {
            T L = 16.;
            T F = 3.6;
            T tref = 0.25;

            return std::min( 1.0, time / tref ) * result_type { 0.0, F / L };
        };

        /* Encast */
        bnd.addDirichletBC( disk::CLAMPED, 1, zero );
        /* Load */
        bnd.addNeumannBC( disk::NEUMANN, 2, trac );
        break;
    }
    case STUDY::SQUARE_DYNA: {
        auto trac = [material_data]( const disk::point< T, 2 > &p, const T &time ) -> result_type {
            T L = 1.;
            T F = 0.2;
            T tref = 1.0;
            const auto force = result_type { 0.0, F / L };
            if ( time <= tref ) {
                return ( time / tref ) * force;
            }

            return force;
        };

        /* BOTTOM */
        bnd.addDirichletBC( disk::CLAMPED, 1, zero );
        /* TOP */
        bnd.addNeumannBC( disk::NEUMANN, 4, trac );
        break;
    }
    case STUDY::SQUARE_MATER: {
        auto trac = [material_data]( const disk::point< T, 2 > &p, const T &time ) -> result_type {
            T L = 1.;
            T F = 5.15e8;
            T tref = 5e-3;
            const auto force = result_type { 0.0, F / L };
            if ( time <= tref ) {
                return ( time / tref ) * force;
            }

            return force;
        };

        /* BOTTOM */
        bnd.addDirichletBC( disk::CLAMPED, 1, zero );
        /* TOP */
        bnd.addNeumannBC( disk::NEUMANN, 4, trac );
        break;
    }
    case STUDY::IMPACT_2D: {

        auto s = []( const disk::point< T, 2 > &p ) -> T { return 0.0; };

        /* Encast */
        bnd.addDirichletBC( disk::CLAMPED, 3, zero );
        /* Syme */
        bnd.addDirichletBC( disk::DX, 1, zero );
        /* Contact */
        auto gap = []( const disk::point< T, 2 > &pt, const disk::static_vector< T, 2 > &n ) -> T {
            // compute the distance to the plane y = 0

            if ( std::abs( n( 1 ) ) < T( 1e-12 ) )
                return T( 1e13 );
            const auto dist = std::abs( pt.y() / n( 1 ) );
            return pt.y() < T( 0 ) ? -dist : dist;
        };

        bnd.addContactBC( disk::SIGNORINI_FACE, 0, s, gap );
        break;
    }
    case STUDY::FREE_VIBR_2D: {

        auto s = []( const disk::point< T, 2 > &p ) -> T { return 0.0; };

        /* Encast */
        bnd.addDirichletBC( disk::CLAMPED, 3, zero );
        /* Syme */
        bnd.addDirichletBC( disk::DX, 1, zero );
        break;
    }
    default: {
        throw std::invalid_argument( "Unexpected study" );
        break;
    }
    }

    return bnd;
}

template < template < typename, size_t, typename > class Mesh, typename T, typename Storage >
auto getBoundaryConditions( const Mesh< T, 3, Storage > &msh,
                            const disk::mechanics::MaterialData< T > &material_data,
                            const STUDY &study ) {
    typedef Mesh< T, 3, Storage > mesh_type;
    typedef disk::static_vector< T, 3 > result_type;

    disk::vector_boundary_conditions< mesh_type > bnd( msh );

    auto zero = [material_data]( const disk::point< T, 3 > &p, const T &time ) -> result_type {
        return result_type { 0.0, 0., 0. };
    };

    /* Boundary conditions */
    switch ( study ) {
    case STUDY::SPHERE_LARGE: {

        auto deplr = [material_data]( const disk::point< T, 3 > &p, const T &time ) -> result_type {
            result_type er = result_type::Zero();

            er( 0 ) = p.x();
            er( 1 ) = p.y();
            er( 2 ) = p.z();

            er /= er.norm();

            return time * 0.157 * er;
        };

        bnd.addDirichletBC( disk::DX, 12, zero );
        bnd.addDirichletBC( disk::DY, 24, zero );
        bnd.addDirichletBC( disk::DZ, 19, zero );
        bnd.addDirichletBC( disk::DIRICHLET, 27, deplr );
        break;
    }
    case STUDY::TAYLOR_ROD: {
        /*RIGHT*/
        bnd.addDirichletBC( disk::DX, 11, zero );
        /*LEFT*/
        bnd.addDirichletBC( disk::DY, 10, zero );
        /*BOTTOM*/
        bnd.addDirichletBC( disk::DZ, 9, zero );
        break;
    }
    default: {
        throw std::invalid_argument( "Unexpected study" );
        break;
    }
    }

    return bnd;
}

template < template < typename, size_t, typename > class Mesh, typename T, typename Storage >
void addExternalLoad( const Mesh< T, 2, Storage > &msh,
                      const disk::mechanics::MaterialData< T > &material_data, const STUDY &study,
                      disk::mechanics::NonLinearSolver< Mesh< T, 2, Storage > > &nl ) {
    typedef Mesh< T, 2, Storage > mesh_type;
    typedef disk::static_vector< T, 2 > result_type;

    auto zero = [material_data]( const disk::point< T, 2 > &p, const T &time ) -> result_type {
        return result_type { 0.0, 0 };
    };

    /* External Load */
    switch ( study ) {
    case STUDY::COOK_ELAS:
    case STUDY::COOK_HPP:
    case STUDY::COOK_LARGE:
    case STUDY::COOK_DYNA:
    case STUDY::SQUARE_DYNA:
    case STUDY::SQUARE_MATER:
    case STUDY::IMPACT_2D:
    case STUDY::FREE_VIBR_2D: {
        break;
    }
    case STUDY::WAVE_ELAS: {

        auto func_space = [material_data]( const disk::point< T, 2 > &p ) -> result_type {
            T ux = -sin( M_PI * p.x() ) * cos( M_PI * p.y() );
            T uy = cos( M_PI * p.x() ) * sin( M_PI * p.y() );

            return result_type { ux, uy };
        };

        auto load = [material_data, func_space]( const disk::point< T, 2 > &p,
                                                 const T &time ) -> result_type {
            const T mu = material_data.getMu();
            const T rho = material_data.getRho();
            const T pi2t2 = M_PI * M_PI * time * time;
            return 2.0 * ( mu * pi2t2 + rho ) * func_space( p );
        };

        nl.addExternalLoad( load );

        break;
    }
    default: {
        throw std::invalid_argument( "Unexpected study" );
        break;
    }
    }
}

template < template < typename, size_t, typename > class Mesh, typename T, typename Storage >
void addExternalLoad( const Mesh< T, 3, Storage > &msh,
                      const disk::mechanics::MaterialData< T > &material_data, const STUDY &study,
                      disk::mechanics::NonLinearSolver< Mesh< T, 3, Storage > > &nl ) {
    typedef Mesh< T, 3, Storage > mesh_type;
    typedef disk::static_vector< T, 3 > result_type;

    auto zero = [material_data]( const disk::point< T, 3 > &p, const T &time ) -> result_type {
        return result_type { 0.0, 0., 0. };
    };

    /* External Load */
    switch ( study ) {
    case STUDY::SPHERE_LARGE:
    case STUDY::TAYLOR_ROD: {
        break;
    }
    default: {
        throw std::invalid_argument( "Unexpected study" );
        break;
    }
    }
}

template < template < typename, size_t, typename > class Mesh, typename T, typename Storage >
void addNonLinearOptions( const Mesh< T, 2, Storage > &msh,
                          const disk::mechanics::MaterialData< T > &material_data,
                          const STUDY &study,
                          disk::mechanics::NonLinearSolver< Mesh< T, 2, Storage > > &nl ) {
    typedef Mesh< T, 2, Storage > mesh_type;
    typedef disk::static_vector< T, 2 > result_type;

    auto zero = [material_data]( const disk::point< T, 2 > &p ) -> result_type {
        return result_type { 0.0, 0. };
    };

    /* Non-linear parameters */
    switch ( study ) {
    case STUDY::COOK_ELAS: {
        nl.addBehavior( disk::mechanics::DeformationMeasure::SMALL_DEF,
                        disk::mechanics::LawType::ELASTIC );

        nl.addPointPlot( { 47.999, 52 }, "pointA.csv" );

        break;
    }
    case STUDY::COOK_HPP: {
        nl.addBehavior( disk::mechanics::DeformationMeasure::SMALL_DEF,
                        disk::mechanics::LawType::LINEAR_HARDENING );

        nl.addPointPlot( { 47.999, 59.999 }, "pointA.csv" );

        break;
    }
    case STUDY::COOK_LARGE: {
        nl.addBehavior( disk::mechanics::DeformationMeasure::LOGARITHMIC_DEF,
                        disk::mechanics::LawType::NONLINEAR_HARDENING );

        nl.addPointPlot( { 47.999, 59.999 }, "pointA.csv" );
        break;
    }
    case STUDY::COOK_DYNA: {
#ifdef HAVE_MGIS
        /* To compile: mfront --obuild --interface=generic LogarithmicStrainPlasticity.mfront */
        // To use a law developped with Mfront
        const auto hypo = mgis::behaviour::Hypothesis::PLANESTRAIN;
        const std::string filename = "src/libBehaviour.so";
        nl.addBehavior( filename, "LogarithmicStrainPlasticity", hypo );
#else
        nl.addBehavior( disk::mechanics::DeformationMeasure::LOGARITHMIC_DEF,
                        disk::mechanics::LawType::LINEAR_HARDENING );
#endif

        nl.addPointPlot( { 47.999, 59.999 }, "pointA.csv" );
        break;
    }
    case STUDY::SQUARE_DYNA:
    case STUDY::SQUARE_MATER: {
#ifdef HAVE_MGIS
        /* To compile: mfront --obuild --interface=generic LogarithmicStrainPlasticity.mfront */
        // To use a law developped with Mfront
        const auto hypo = mgis::behaviour::Hypothesis::PLANESTRAIN;
        const std::string filename = "src/libBehaviour.so";
        nl.addBehavior( filename, "LogarithmicStrainPlasticity", hypo );
#else
        nl.addBehavior( disk::mechanics::DeformationMeasure::LOGARITHMIC_DEF,
                        disk::mechanics::LawType::LINEAR_HARDENING );
#endif

        nl.addPointPlot( { 0.999, 0.999 }, "pointA.csv" );
        break;
    }
    case STUDY::WAVE_ELAS: {
        nl.addBehavior( disk::mechanics::DeformationMeasure::SMALL_DEF,
                        disk::mechanics::LawType::ELASTIC );

        break;
    }
    case STUDY::IMPACT_2D: {
        nl.addBehavior( disk::mechanics::DeformationMeasure::SMALL_DEF,
                        disk::mechanics::LawType::ELASTIC );

        auto u0 = []( const disk::point< T, 2 > &p ) -> result_type {
            T y = p.y();
            T x = p.x();

            return result_type { 0.0, 0.5 * ( 1.0 - y ) };
        };

        nl.initial_guess( u0 );
        nl.addPointPlot( { 0.0025, 0.0 }, "pointA.csv" );
        nl.addPointPlot( { 0.0025, 1.0 }, "pointB.csv" );
        break;
    }
    case STUDY::FREE_VIBR_2D: {
        nl.addBehavior( disk::mechanics::DeformationMeasure::SMALL_DEF,
                        disk::mechanics::LawType::ELASTIC );

        const T A = 0.5;
        auto u0 = [A]( const disk::point< T, 2 > &p ) -> result_type {
            T y = p.y();

            return result_type { 0.0, A * std::cos( 0.5 * M_PI * y ) };
        };

        auto v0 = []( const disk::point< T, 2 > &p ) -> result_type {
            return result_type { 0.0, 0.0 };
        };

        auto a0 = [A]( const disk::point< T, 2 > &p ) -> result_type {
            T y = p.y();

            return result_type { 0.0, -0.25 * A * M_PI * M_PI * std::cos( 0.5 * M_PI * y ) };
        };

        nl.initial_field( disk::mechanics::FieldName::DEPL, u0 );
        nl.initial_field( disk::mechanics::FieldName::VITE, v0 );
        nl.initial_field( disk::mechanics::FieldName::ACCE, a0 );

        nl.addPointPlot( { 0.0025, 0.0 }, "pointA.csv" );
        nl.addPointPlot( { 0.0025, 1.0 }, "pointB.csv" );

        break;
    }
    default: {
        throw std::invalid_argument( "Unexpected study" );
        break;
    }
    }

    // Add after behavior
    nl.addMaterialData( material_data );
}

template < template < typename, size_t, typename > class Mesh, typename T, typename Storage >
void addNonLinearOptions( const Mesh< T, 3, Storage > &msh,
                          const disk::mechanics::MaterialData< T > &material_data,
                          const STUDY &study,
                          disk::mechanics::NonLinearSolver< Mesh< T, 3, Storage > > &nl ) {
    typedef Mesh< T, 3, Storage > mesh_type;
    typedef disk::static_vector< T, 3 > result_type;

    auto zero = [material_data]( const disk::point< T, 3 > &p ) -> result_type {
        return result_type { 0.0, 0., 0. };
    };

    /* Non-linear parameters */
    switch ( study ) {
    case STUDY::SPHERE_LARGE: {
#ifdef HAVE_MGIS
        // To use a law developped with Mfront
        const auto hypo = mgis::behaviour::Hypothesis::TRIDIMENSIONAL;
        const std::string filename = "src/libBehaviour.so";
        nl.addBehavior( filename, "LogarithmicStrainPlasticity", hypo );
#else
        // To use a native law from DiSk++
        nl.addBehavior( disk::mechanics::DeformationMeasure::LOGARITHMIC_DEF,
                        disk::mechanics::LawType::LINEAR_HARDENING );
#endif
        break;
    }
    case STUDY::TAYLOR_ROD: {
#ifdef HAVE_MGIS
        /* To compile: mfront --obuild --interface=generic LogarithmicStrainPlasticity.mfront */
        // To use a law developped with Mfront
        const auto hypo = mgis::behaviour::Hypothesis::TRIDIMENSIONAL;
        const std::string filename = "src/libBehaviour.so";
        nl.addBehavior( filename, "LogarithmicStrainPlasticity", hypo );
#else
        // To use a native law from DiSk++
        nl.addBehavior( disk::mechanics::DeformationMeasure::LOGARITHMIC_DEF,
                        disk::mechanics::LawType::LINEAR_HARDENING );
#endif
        nl.initial_field( disk::mechanics::FieldName::VITE_CELLS,
                          []( const disk::point< T, 3 > &p ) -> auto {
                              return result_type { 0.0, 0.0, -227.0 };
                          } );
        nl.addPointPlot( { -0.00001, 3.19999, 0. }, "pointA.csv" );

        break;
    }
    default: {
        throw std::invalid_argument( "Unexpected study" );
        break;
    }
    }

    // Add after behavior
    nl.addMaterialData( material_data );
}
