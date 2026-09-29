/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2026
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
#include "diskpp/boundary_conditions/boundary_conditions.hpp"
#include "diskpp/common/eigen.hpp"
#include "diskpp/common/timecounter.hpp"
#include "diskpp/mechanics/NewtonSolver/NonLinearParameters.hpp"
#include "diskpp/mechanics/behaviors/maths_tensor.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace disk {
namespace mechanics {
namespace priv {

template < typename Basis, typename T, std::size_t DIM >
point< T, DIM >
new_pt( const Basis &basis,
        const disk::dynamic_vector< T > &coefficients,
        const point< T, DIM > &point ) {
    return point + disk::eval( coefficients, basis.eval_functions( point ) );
}

template < typename T >
static_vector< T, 2 >
compute_normal( const point< T, 2 > &a, const point< T, 2 > &b ) {
    const static_vector< T, 2 > tangent = ( b - a ).to_vector();
    const T norm = tangent.norm();
    if ( norm <= T( 100 ) * std::numeric_limits< T >::epsilon() )
        throw std::runtime_error( "Cannot compute the normal of a degenerate edge" );
    static_vector< T, 2 > normal;
    normal << -tangent( 1 ), tangent( 0 );
    return normal / norm;
}

template < typename T >
static_vector< T, 3 >
compute_normal( const point< T, 3 > &a, const point< T, 3 > &b, const point< T, 3 > &c ) {
    const static_vector< T, 3 > area = ( b - a ).to_vector().cross( ( c - a ).to_vector() );
    const T norm = area.norm();
    if ( norm <= T( 100 ) * std::numeric_limits< T >::epsilon() )
        throw std::runtime_error( "Cannot compute the normal of a degenerate face" );
    return area / norm;
}

} // namespace priv

template < typename T, int DIM >
struct NormalKinematics {
    ContactKinematics cont_kine;
    static_vector< T, DIM > reference_normal;
    static_vector< T, DIM > current_normal;
    static_vector< T, DIM > kinematic_normal;

    template < typename Mesh, typename Elem, typename Basis >
    NormalKinematics( const ContactKinematics kine,
                      const Mesh &msh,
                      const Elem &elem,
                      const Basis &basis,
                      const disk::dynamic_vector< T > &coefficients,
                      const static_vector< T, DIM > &reference )
        : cont_kine( kine ),
          reference_normal( reference ),
          current_normal( reference ),
          kinematic_normal( reference ) {
        const auto element_points = points( msh, elem );
        if constexpr ( DIM == 2 ) {
            if ( element_points.size() < 2 )
                throw std::runtime_error( "A 2D contact edge must contain at least two points" );
            const auto x0 = priv::new_pt( basis, coefficients, element_points[0] );
            const auto x1 = priv::new_pt( basis, coefficients, element_points[1] );
            const auto initial = priv::compute_normal( element_points[0], element_points[1] );
            const T orientation = std::copysign( T( 1 ), reference.dot( initial ) );
            current_normal = orientation * priv::compute_normal( x0, x1 );
        } else if constexpr ( DIM == 3 ) {
            if ( element_points.size() < 3 )
                throw std::runtime_error( "A 3D contact face must contain at least three points" );
            const auto x0 = priv::new_pt( basis, coefficients, element_points[0] );
            const auto x1 = priv::new_pt( basis, coefficients, element_points[1] );
            const auto x2 = priv::new_pt( basis, coefficients, element_points[2] );
            const auto initial =
                priv::compute_normal( element_points[0], element_points[1], element_points[2] );
            const T orientation = std::copysign( T( 1 ), reference.dot( initial ) );
            current_normal = orientation * priv::compute_normal( x0, x1, x2 );
        } else {
            static_assert( DIM == 2 || DIM == 3,
                           "NormalKinematics supports only dimensions 2 and 3" );
        }
        kinematic_normal = kine == ContactKinematics::REFERENCE ? reference_normal : current_normal;
    }
};

/* CTAD guide for the six-argument call used by NewtonSolverContact.hpp. */
template < typename Mesh, typename Elem, typename Basis, typename T, int DIM >
NormalKinematics( ContactKinematics,
                  const Mesh &,
                  const Elem &,
                  const Basis &,
                  const disk::dynamic_vector< T > &,
                  const static_vector< T, DIM > & ) -> NormalKinematics< T, DIM >;

namespace priv {

template < typename T >
struct GapLinearization {
    T gap = T( 0 );
    dynamic_vector< T > derivative;
    bool valid = false;
};

template < typename Basis, typename T, std::size_t DIM, typename FunctionGap >
T
compute_gap_fb( const Basis &basis,
                const disk::dynamic_vector< T > &coefficients,
                const FunctionGap &gap_function,
                const point< T, DIM > &point,
                const static_vector< T, static_cast< int >( DIM ) > &kinematic_normal,
                const T time ) {
    return gap_function( new_pt( basis, coefficients, point ), kinematic_normal, time );
}

/* Compatibility with old eight-argument call sites. */
template < typename Mesh,
           typename Elem,
           typename Basis,
           typename T,
           std::size_t DIM,
           typename FunctionGap >
T
compute_gap_fb( const Mesh &,
                const Elem &,
                const Basis &basis,
                const disk::dynamic_vector< T > &coefficients,
                const FunctionGap &gap_function,
                const point< T, DIM > &point,
                const static_vector< T, static_cast< int >( DIM ) > &kinematic_normal,
                const T time ) {
    return compute_gap_fb( basis, coefficients, gap_function, point, kinematic_normal, time );
}

template < typename T >
bool
is_valid_gap( const T gap ) {
    return std::isfinite( gap );
}

template < typename Basis, typename T, std::size_t DIM, typename FunctionGap >
GapLinearization< T >
linearize_gap_fb( const Basis &basis,
                  const disk::dynamic_vector< T > &coefficients,
                  const FunctionGap &gap_function,
                  const point< T, DIM > &point,
                  const ContactKinematics cont_kine,
                  const static_vector< T, static_cast< int >( DIM ) > &kinematic_normal,
                  const T time ) {
    GapLinearization< T > result;
    const Eigen::Index size = coefficients.size();
    result.derivative = dynamic_vector< T >::Zero( size );
    result.gap = compute_gap_fb( basis, coefficients, gap_function, point, kinematic_normal, time );
    result.valid = is_valid_gap( result.gap );
    if ( !result.valid )
        return result;

    if ( cont_kine == ContactKinematics::REFERENCE ) {
        // -u.n
        const auto t_phi = basis.eval_functions( point );
        result.derivative = -disk::priv::inner_product( t_phi, kinematic_normal );
    } else {
        const T relative_step = std::cbrt( std::numeric_limits< T >::epsilon() );
        for ( Eigen::Index i = 0; i < size; ++i ) {
            const T h = relative_step * std::max( T( 1 ), std::abs( coefficients( i ) ) );
            auto plus = coefficients;
            auto minus = coefficients;
            plus( i ) += h;
            minus( i ) -= h;
            const T gp = compute_gap_fb( basis, plus, gap_function, point, kinematic_normal, time );
            const T gm =
                compute_gap_fb( basis, minus, gap_function, point, kinematic_normal, time );
            const bool vp = is_valid_gap( gp );
            const bool vm = is_valid_gap( gm );
            if ( vp && vm )
                result.derivative( i ) = ( gp - gm ) / ( T( 2 ) * h );
            else if ( vp )
                result.derivative( i ) = ( gp - result.gap ) / h;
            else if ( vm )
                result.derivative( i ) = ( result.gap - gm ) / h;
        }
    }
    return result;
}

} // namespace priv
} // namespace mechanics
} // namespace disk