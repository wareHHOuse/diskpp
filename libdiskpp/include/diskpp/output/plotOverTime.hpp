/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2025                     nicolas.pignet@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */

#pragma once

#include "diskpp/common/eigen.hpp"
#include "diskpp/geometry/geometry.hpp"
#include "diskpp/mesh/point.hpp"

#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace disk {

template < typename Mesh >
class PlotPointOverTime {
  private:
    typedef Mesh mesh_type;
    typedef typename mesh_type::point_type point_type;

    typedef std::pair< double, std::vector< double > > vale_type;

    int _c_id;
    point_type _point;

    std::vector< vale_type > _vale;
    std::vector< std::string > _comp;

    std::string _filename;

  public:
    PlotPointOverTime() : _c_id( -1 ) {}

    PlotPointOverTime( const Mesh &mesh, const point_type &point ) : _point( point ) {
        // Find cell id;

        bool find = false;

        for ( auto &cl : mesh ) {
            if ( is_inside( mesh, cl, _point ) ) {
                _c_id = mesh.lookup( cl );
                find = true;
                break;
            }
        }

        if ( !find ) {
            throw std::runtime_error( "Point not find in the mesh" );
        }
    }

    int getCellId() const { return _c_id; }
    auto getPoint() const { return _point; }

    void setFilename( const std::string &filename ) { _filename = filename; }

    void addComponents( const std::vector< std::string > &comp ) { _comp = comp; }
    void addValues( const double &time, const std::vector< double > &vals ) {
        _vale.push_back( std::make_pair( time, vals ) );
    }

    template < typename T, int N >
    void addValues( const double &time, const static_vector< T, N > &vals ) {
        std::vector< double > val( N );

        for ( int i = 0; i < N; i++ ) {
            val[i] = vals[i];
        }

        _vale.push_back( std::make_pair( time, val ) );
    }

    template < typename T, int N >
    void addValues( const double &time, const int &n_iter,
                    const std::vector< static_vector< T, N > > &vals ) {
        std::vector< double > val;
        val.push_back( n_iter );

        for ( auto &vec : vals ) {
            for ( int i = 0; i < N; i++ ) {
                val.push_back( vec[i] );
            }
        }

        _vale.push_back( std::make_pair( time, val ) );
    }

    template < typename T, int N >
    void addValues( const double &time, const std::vector< static_vector< T, N > > &vals ) {
        std::vector< double > val;

        for ( auto &vec : vals ) {
            for ( int i = 0; i < N; i++ ) {
                val.push_back( vec[i] );
            }
        }

        _vale.push_back( std::make_pair( time, val ) );
    }

    void addValues( const double &time, const std::map< std::string, double > &vals ) {
        std::vector< double > val;

        bool write_cmp = _comp.empty();

        for ( auto &[key, value] : vals ) {
            if ( write_cmp ) {
                _comp.push_back( key );
            }
            val.push_back( value );
        }

        _vale.push_back( std::make_pair( time, val ) );
    }

    void write() const {

        std::ofstream fio( _filename, std::ofstream::trunc );

        if ( fio.is_open() ) {

            fio << "Time ";
            for ( auto &cmp : _comp ) {
                fio << "; " << cmp;
            }
            fio << std::endl;

            for ( auto &[time, vals] : _vale ) {
                fio << time;
                for ( auto &val : vals ) {
                    fio << "; " << val;
                }
                fio << std::endl;
            }

        } else {
            throw std::runtime_error( "Error when opening the file." );
        }

        fio.close();
    }
};

} // namespace disk