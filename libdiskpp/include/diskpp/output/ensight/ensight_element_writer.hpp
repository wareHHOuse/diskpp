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

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>

namespace disk::output::ensight {
inline void
write_scalar_per_element( const std::filesystem::path &path,
                          const std::string &description,
                          const Geometry &geometry,
                          const std::vector< double > &values ) {
    if ( values.size() != geometry.number_of_elements() )
        throw std::invalid_argument( "Element field size mismatch" );
    std::ofstream os( path );
    if ( !os )
        throw std::runtime_error( "Cannot open EnSight element field: " + path.string() );
    os << std::scientific << std::setprecision( 16 ) << description << '\n'
       << "part\n"
       << geometry.part_id << '\n';

    const auto write_block = [&]( const char *keyword, ElementBlock block, std::size_t count ) {
        if ( count == 0 )
            return;
        std::vector< double > ordered( count, 0.0 );
        for ( std::size_t cell = 0; cell < geometry.element_reference.size(); ++cell ) {
            const auto &ref = geometry.element_reference[cell];
            if ( ref.block == block )
                ordered.at( ref.index_in_block ) = values.at( cell );
        }
        os << keyword << '\n';
        for ( const auto v : ordered )
            os << v << '\n';
    };
    write_block( "bar2", ElementBlock::bar2, geometry.bars.size() );
    write_block( "tria3", ElementBlock::tria3, geometry.triangles.size() );
    write_block( "quad4", ElementBlock::quad4, geometry.quadrilaterals.size() );
    write_block( "nsided", ElementBlock::nsided, geometry.polygons.size() );
    write_block( "nfaced", ElementBlock::nfaced, geometry.polyhedra.size() );
}
} // namespace disk::output::ensight
