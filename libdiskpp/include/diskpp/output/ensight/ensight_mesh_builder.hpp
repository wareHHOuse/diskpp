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

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace disk::output::ensight {
namespace detail {
template < typename Identifier >
std::size_t
point_identifier_value( const Identifier &id ) {
    return static_cast< std::size_t >( id );
}
} // namespace detail

template < typename Mesh >
Geometry
build_geometry( const Mesh &msh, const std::string &part_name = "HHO mesh" ) {
    constexpr std::size_t DIM = Mesh::dimension;
    static_assert( DIM >= 1 && DIM <= 3, "EnSight export supports dimensions 1, 2 and 3" );

    Geometry out;
    out.part_name = part_name;
    out.points.reserve( msh.points_size() );
    out.element_reference.resize( msh.cells_size() );

    std::map< std::size_t, std::int64_t > ensight_id;
    std::size_t point_index = 0;
    for ( auto it = msh.points_begin(); it != msh.points_end(); ++it, ++point_index ) {
        const auto &p = *it;
        std::array< double, 3 > xyz { 0.0, 0.0, 0.0 };
        xyz[0] = static_cast< double >( p.x() );
        if constexpr ( DIM >= 2 )
            xyz[1] = static_cast< double >( p.y() );
        if constexpr ( DIM == 3 )
            xyz[2] = static_cast< double >( p.z() );
        out.points.push_back( xyz );
        ensight_id.emplace( point_index, static_cast< std::int64_t >( point_index + 1 ) );
    }

    const auto node_id = [&]( const auto &pid ) -> std::int64_t {
        const auto found = ensight_id.find( detail::point_identifier_value( pid ) );
        if ( found == ensight_id.end() )
            throw std::runtime_error( "Unknown point identifier" );
        return found->second;
    };

    for ( const auto &cl : msh ) {
        const std::size_t cell_id = msh.lookup( cl );
        if constexpr ( DIM == 1 ) {
            const auto ids = cl.point_ids();
            if ( ids.size() != 2 )
                throw std::runtime_error( "A 1D cell must have two points" );
            const auto index = out.bars.size();
            out.bars.push_back( { node_id( ids[0] ), node_id( ids[1] ) } );
            out.element_reference.at( cell_id ) = { ElementBlock::bar2, index };
        } else if constexpr ( DIM == 2 ) {
            const auto ids = cl.point_ids();
            if ( ids.size() < 3 )
                throw std::runtime_error( "A 2D cell must have at least three points" );
            if ( ids.size() == 3 ) {
                const auto index = out.triangles.size();
                out.triangles.push_back(
                    { node_id( ids[0] ), node_id( ids[1] ), node_id( ids[2] ) } );
                out.element_reference.at( cell_id ) = { ElementBlock::tria3, index };
            } else if ( ids.size() == 4 ) {
                const auto index = out.quadrilaterals.size();
                out.quadrilaterals.push_back( { node_id( ids[0] ),
                                                node_id( ids[1] ),
                                                node_id( ids[2] ),
                                                node_id( ids[3] ) } );
                out.element_reference.at( cell_id ) = { ElementBlock::quad4, index };
            } else {
                const auto index = out.polygons.size();
                std::vector< std::int64_t > polygon;
                polygon.reserve( ids.size() );
                for ( const auto pid : ids )
                    polygon.push_back( node_id( pid ) );
                out.polygons.push_back( std::move( polygon ) );
                out.element_reference.at( cell_id ) = { ElementBlock::nsided, index };
            }
        } else {
            const auto index = out.polyhedra.size();
            NfacedCell cell;
            for ( const auto &fc : faces( msh, cl ) ) {
                std::vector< std::int64_t > face;
                for ( const auto pid : fc.point_ids() )
                    face.push_back( node_id( pid ) );
                if ( face.size() < 3 )
                    throw std::runtime_error( "A polyhedron face must have at least three points" );
                cell.faces.push_back( std::move( face ) );
            }
            out.polyhedra.push_back( std::move( cell ) );
            out.element_reference.at( cell_id ) = { ElementBlock::nfaced, index };
        }
    }
    return out;
}

} // namespace disk::output::ensight
