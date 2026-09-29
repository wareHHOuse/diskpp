/*
 * /\ Matteo Cicuttin (C) 2016, 2017, 2018
 * /__\ matteo.cicuttin@enpc.fr
 * /_\/_\ École Nationale des Ponts et Chaussées - CERMICS
 * /\ /\
 * /__\ /__\ DISK++, a template library for DIscontinuous SKeletal
 * /_\/_\/_\/_\ methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet (C) 2026
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
#include <map>
#include <stdexcept>
#include <string>

namespace disk::output::ensight {

inline void
write_case_file( const std::filesystem::path &path,
                 const std::filesystem::path &geometry,
                 const std::map< std::string, FieldDescription > &fields,
                 const std::map< std::size_t, double > &times,
                 const std::size_t width = 6 ) {
    if ( geometry.empty() ) {
        throw std::invalid_argument( "An EnSight model geometry filename is required" );
    }
    std::ofstream os( path );
    if ( !os ) {
        throw std::runtime_error( "Cannot open EnSight case file: " + path.string() );
    }
    os << "FORMAT\n"
       << "type: ensight gold\n\n";

    /*
     * Each case file describes one standalone EnSight model:
     *
     * simulation.case -> mesh.geo
     * gauss_points.case -> gauss_points.geo
     */
    os << "GEOMETRY\n"
       << "model: " << geometry.generic_string() << "\n\n";
    if ( !fields.empty() ) {
        os << "VARIABLE\n";
        for ( const auto &[field_key, field] : fields ) {
            static_cast< void >( field_key );
            os << case_keyword( field.type, field.location ) << ": 1 " << field.name << ' '
               << field.safe_name << '/' << field.safe_name << '.' << std::string( width, '*' )
               << '\n';
        }
        os << '\n';
    }
    if ( !times.empty() ) {
        os << "TIME\n"
           << "time set: 1\n"
           << "number of steps: " << times.size() << '\n'
           << "filename start number: 0\n"
           << "filename increment: 1\n"
           << "time values:\n"
           << std::scientific << std::setprecision( 16 );
        for ( std::size_t step = 0; step < times.size(); ++step ) {
            const auto iterator = times.find( step );
            if ( iterator == times.end() ) {
                throw std::runtime_error( "Non-contiguous EnSight time steps" );
            }
            os << iterator->second << '\n';
        }
    }
    if ( !os ) {
        throw std::runtime_error( "Error while writing EnSight case file: " + path.string() );
    }
}

} // namespace disk::output::ensight