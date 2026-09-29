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

#include "ensight_types.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace disk::output::ensight {
namespace detail {
template < typename V >
auto
vector_component_impl( const V &v, std::size_t i, int )
    -> decltype( static_cast< double >( v( i ) ) ) {
    return static_cast< double >( v( i ) );
}
template < typename V >
auto
vector_component_impl( const V &v, std::size_t i, long )
    -> decltype( static_cast< double >( v[i] ) ) {
    return static_cast< double >( v[i] );
}
template < typename V >
double
vector_component( const V &v, std::size_t i ) {
    return vector_component_impl( v, i, 0 );
}
template < typename V >
auto
vector_size_impl( const V &v, int ) -> decltype( static_cast< std::size_t >( v.size() ) ) {
    return static_cast< std::size_t >( v.size() );
}
template < typename V >
std::size_t
vector_size( const V &v ) {
    return vector_size_impl( v, 0 );
}

template < typename M >
auto
matrix_component_impl( const M &m, std::size_t i, std::size_t j, int )
    -> decltype( static_cast< double >( m( i, j ) ) ) {
    return static_cast< double >( m( i, j ) );
}
template < typename M >
auto
matrix_component_impl( const M &m, std::size_t i, std::size_t j, long )
    -> decltype( static_cast< double >( m[i][j] ) ) {
    return static_cast< double >( m[i][j] );
}
template < typename M >
double
matrix_component( const M &m, std::size_t i, std::size_t j ) {
    return matrix_component_impl( m, i, j, 0 );
}
template < typename M >
auto
matrix_rows_impl( const M &m, int ) -> decltype( static_cast< std::size_t >( m.rows() ) ) {
    return static_cast< std::size_t >( m.rows() );
}
template < typename M >
auto
matrix_rows_impl( const M &m, long ) -> decltype( static_cast< std::size_t >( m.size() ) ) {
    return static_cast< std::size_t >( m.size() );
}
template < typename M >
std::size_t
matrix_rows( const M &m ) {
    return matrix_rows_impl( m, 0 );
}
template < typename M >
auto
matrix_cols_impl( const M &m, int ) -> decltype( static_cast< std::size_t >( m.cols() ) ) {
    return static_cast< std::size_t >( m.cols() );
}
template < typename M >
auto
matrix_cols_impl( const M &m, long ) -> decltype( static_cast< std::size_t >( m[0].size() ) ) {
    return m.size() == 0 ? 0 : static_cast< std::size_t >( m[0].size() );
}
template < typename M >
std::size_t
matrix_cols( const M &m ) {
    return matrix_cols_impl( m, 0 );
}
} // namespace detail

template < typename Getter >
void
write_nodal_field_streaming( const std::filesystem::path &path,
                             const std::string &description,
                             FieldType type,
                             std::size_t count,
                             Getter &&getter,
                             std::int64_t part_id ) {
    std::ofstream os( path );
    if ( !os )
        throw std::runtime_error( "Cannot open EnSight nodal field: " + path.string() );
    os << std::scientific << std::setprecision( 16 ) << description << '\n'
       << "part\n"
       << part_id << '\n'
       << "coordinates\n";
    const auto write_component = [&]( std::size_t c ) {
        for ( std::size_t i = 0; i < count; ++i )
            os << getter( i, c ) << '\n';
    };
    switch ( type ) {
    case FieldType::scalar:
        write_component( 0 );
        break;
    case FieldType::vector:
        write_component( 0 );
        write_component( 1 );
        write_component( 2 );
        break;
    case FieldType::symmetric_tensor: {
        for ( std::size_t c = 0; c < 6; ++c )
            write_component( c );
        break;
    }
    case FieldType::tensor: {
        for ( const auto c : { 0u, 4u, 8u, 1u, 5u, 6u, 3u, 7u, 2u } )
            write_component( c );
        break;
    }
    }
}

} // namespace disk::output::ensight
