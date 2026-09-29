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

#include "ensight_case_writer.hpp"
#include "ensight_element_writer.hpp"
#include "ensight_geometry_writer.hpp"
#include "ensight_mesh_builder.hpp"
#include "ensight_nodal_writer.hpp"
#include "ensight_point_cloud_writer.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace disk::output::ensight {
inline std::string
safe_filename( std::string name ) {
    for ( auto &c : name ) {
        auto u = static_cast< unsigned char >( c );
        if ( !( std::isalnum( u ) || c == '_' || c == '-' ) )
            c = '_';
    }
    return name.empty() ? "field" : name;
}
inline std::filesystem::path
create_output_directory( const std::filesystem::path &parent = ".",
                         const std::string &prefix = "results" ) {
    const auto now = std::chrono::system_clock::now();
    const auto sec = std::chrono::time_point_cast< std::chrono::seconds >( now );
    const auto ms = std::chrono::duration_cast< std::chrono::milliseconds >( now - sec ).count();
    const auto raw = std::chrono::system_clock::to_time_t( now );
    std::tm tm {};
#if defined( _WIN32 )
    localtime_s( &tm, &raw );
#else
    localtime_r( &raw, &tm );
#endif
    std::ostringstream base;
    base << prefix << '_' << std::put_time( &tm, "%Y%m%d_%H%M%S" ) << '_' << std::setw( 3 )
         << std::setfill( '0' ) << ms;
    auto dir = parent / base.str();
    std::size_t n = 0;
    while ( !std::filesystem::create_directories( dir ) )
        dir = parent / ( base.str() + "_" + std::to_string( ++n ) );
    return dir;
}

template < typename Mesh >
class EnsightExporter {
  public:
    EnsightExporter( const Mesh &msh,
                     std::filesystem::path out = create_output_directory(),
                     std::string case_name = "simulation.case",
                     std::string geo_name = "mesh.geo" )
        : m_output_directory( std::move( out ) ),
          m_case_filename( std::move( case_name ) ),
          m_geometry_filename( std::move( geo_name ) ),
          m_geometry( build_geometry( msh ) ) {
        std::filesystem::create_directories( m_output_directory );
    }

    const std::filesystem::path &
    output_directory() const noexcept {
        return m_output_directory;
    }

    std::int64_t
    add_point_cloud( const std::string &name, std::vector< std::array< double, 3 > > points ) {
        if ( m_mesh_written )
            throw std::logic_error( "Point clouds must be added before write_mesh()" );
        const auto id = static_cast< std::int64_t >( 2 + m_point_clouds.size() );
        m_point_clouds.push_back( { name, id, std::move( points ) } );
        return id;
    }

    const PointCloud &
    point_cloud( std::int64_t id ) const {
        auto it = std::find_if( m_point_clouds.begin(), m_point_clouds.end(), [&]( const auto &c ) {
            return c.part_id == id;
        } );
        if ( it == m_point_clouds.end() )
            throw std::out_of_range( "Unknown point cloud" );
        return *it;
    }

    void
    write_mesh() {
        /*
         * Write only the HHO mesh into the main EnSight geometry.
         */
        write_geometry( m_output_directory / m_geometry_filename, m_geometry );

        /*
         * Write the Gauss-point cloud as EnSight point cloud geometry.
         */
        if ( !m_point_clouds.empty() ) {
            if ( m_point_clouds.size() != 1 ) {
                throw std::logic_error( "Only one EnSight point cloud point cloud "
                                        "is currently supported" );
            }

            write_point_cloud_geometry( m_output_directory / m_point_cloud_geometry_filename,
                                        m_point_clouds.front() );
        }

        m_mesh_written = true;

        update_case_file();
    }

    std::size_t
    begin_step( double time ) {
        require_mesh();
        if ( m_current_step )
            throw std::logic_error( "A step is already active" );
        const auto step = m_times.size();
        m_times.emplace( step, time );
        m_current_step = step;
        return step;
    }

    void
    end_step() {
        require_active_step();
        m_current_step.reset();
        update_case_file();
    }

    std::size_t
    current_step() const {
        require_active_step();
        return *m_current_step;
    }

    double
    current_time() const {
        return m_times.at( current_step() );
    }

    template < typename C >
    void
    write_scalar( const std::string &n, const C &v ) {
        auto s = active();
        write_scalar( n, s.first, s.second, v );
    }

    template < typename C >
    void
    write_scalar( const std::string &n, std::size_t s, double t, const C &v ) {
        if ( v.size() != m_geometry.number_of_points() )
            throw std::invalid_argument( "Node field size mismatch" );
        write_node_field(
            n, FieldType::scalar, s, t, m_geometry.part_id, v.size(), [&]( auto i, auto ) {
                return static_cast< double >( v[i] );
            } );
    }

    template < typename C >
    void
    write_vector( const std::string &n, const C &v ) {
        auto s = active();
        write_vector( n, s.first, s.second, v );
    }

    template < typename C >
    void
    write_vector( const std::string &n, std::size_t s, double t, const C &v ) {
        if ( v.size() != m_geometry.number_of_points() || v.empty() )
            throw std::invalid_argument( "Vector field size mismatch" );
        const auto dim = detail::vector_size( v[0] );
        if ( dim < 1 || dim > 3 )
            throw std::invalid_argument( "Vector dimension must be 1..3" );
        write_node_field(
            n, FieldType::vector, s, t, m_geometry.part_id, v.size(), [&]( auto i, auto c ) {
                return c < dim ? detail::vector_component( v[i], c ) : 0.0;
            } );
    }

    template < typename C >
    void
    write_symmetric_tensor( const std::string &n, const C &v ) {
        auto s = active();
        write_symmetric_tensor( n, s.first, s.second, v );
    }

    template < typename C >
    void
    write_symmetric_tensor( const std::string &n, std::size_t s, double t, const C &v ) {
        write_matrix_node_field( n,
                                 FieldType::symmetric_tensor,
                                 s,
                                 t,
                                 m_geometry.part_id,
                                 m_geometry.number_of_points(),
                                 v,
                                 true );
    }

    template < typename C >
    void
    write_tensor( const std::string &n, const C &v ) {
        auto s = active();
        write_tensor( n, s.first, s.second, v );
    }

    template < typename C >
    void
    write_tensor( const std::string &n, std::size_t s, double t, const C &v ) {
        write_matrix_node_field( n,
                                 FieldType::tensor,
                                 s,
                                 t,
                                 m_geometry.part_id,
                                 m_geometry.number_of_points(),
                                 v,
                                 false );
    }

    template < typename C >
    void
    write_element_scalar( const std::string &n, const C &v ) {
        auto s = active();
        write_element_scalar( n, s.first, s.second, v );
    }

    template < typename C >
    void
    write_element_scalar( const std::string &n, std::size_t s, double t, const C &v ) {
        if ( v.size() != m_geometry.number_of_elements() )
            throw std::invalid_argument( "Element field size mismatch" );
        register_time( s, t );
        auto &f =
            register_field( n, FieldType::scalar, FieldLocation::element, m_geometry.part_id, s );
        std::vector< double > d;
        d.reserve( v.size() );
        for ( const auto &x : v )
            d.push_back( static_cast< double >( x ) );
        const auto path = field_path( f, s );
        write_scalar_per_element( path, n, m_geometry, d );
        f.written_steps.insert( s );
        if ( !m_current_step )
            update_case_file();
    }

    template < typename C >
    void
    write_point_cloud_scalar( const std::string &n, std::int64_t part, const C &v ) {
        auto s = active();
        write_point_cloud_scalar( n, s.first, s.second, part, v );
    }

    template < typename ScalarContainer >
    void
    write_point_cloud_scalar( const std::string &name,
                              const std::size_t step,
                              const double time,
                              const std::int64_t point_cloud_id,
                              const ScalarContainer &values ) {
        const auto &cloud = point_cloud( point_cloud_id );

        if ( values.size() != cloud.number_of_points() ) {
            throw std::invalid_argument( "Point-cloud scalar field size mismatch" );
        }

        write_point_cloud_field( name,
                                 FieldType::scalar,
                                 step,
                                 time,
                                 point_cloud_id,
                                 values.size(),
                                 [&]( const std::size_t point_id, const std::size_t ) {
                                     return static_cast< double >( values[point_id] );
                                 } );
    }

    template < typename C >
    void
    write_point_cloud_vector( const std::string &n, std::int64_t part, const C &v ) {
        auto s = active();
        write_point_cloud_vector( n, s.first, s.second, part, v );
    }

    template < typename VectorContainer >
    void
    write_point_cloud_vector( const std::string &name,
                              const std::size_t step,
                              const double time,
                              const std::int64_t point_cloud_id,
                              const VectorContainer &values ) {
        const auto &cloud = point_cloud( point_cloud_id );

        if ( values.size() != cloud.number_of_points() ) {
            throw std::invalid_argument( "Point-cloud vector field size mismatch" );
        }

        if ( values.empty() ) {
            throw std::invalid_argument( "Cannot write an empty point-cloud vector field" );
        }

        const std::size_t dimension = detail::vector_size( values.front() );

        if ( dimension == 0 || dimension > 3 ) {
            throw std::invalid_argument( "Point-cloud vector dimension must "
                                         "be between one and three" );
        }

        write_point_cloud_field( name,
                                 FieldType::vector,
                                 step,
                                 time,
                                 point_cloud_id,
                                 values.size(),
                                 [&]( const std::size_t point_id, const std::size_t component ) {
                                     if ( component >= dimension ) {
                                         return 0.0;
                                     }

                                     return detail::vector_component( values[point_id], component );
                                 } );
    }

    template < typename C >
    void
    write_point_cloud_symmetric_tensor( const std::string &n, std::int64_t part, const C &v ) {
        auto s = active();
        write_point_cloud_symmetric_tensor( n, s.first, s.second, part, v );
    }

    template < typename TensorContainer >
    void
    write_point_cloud_symmetric_tensor( const std::string &name,
                                        const std::size_t step,
                                        const double time,
                                        const std::int64_t point_cloud_id,
                                        const TensorContainer &values ) {
        const auto &cloud = point_cloud( point_cloud_id );

        if ( values.size() != cloud.number_of_points() ) {
            throw std::invalid_argument( "Point-cloud symmetric tensor field size mismatch" );
        }

        if ( values.empty() ) {
            throw std::invalid_argument( "Cannot write an empty point-cloud tensor field" );
        }

        const std::size_t rows = detail::matrix_rows( values.front() );

        const std::size_t columns = detail::matrix_cols( values.front() );

        if ( rows == 0 || columns == 0 || rows != columns || rows > 3 ) {
            throw std::invalid_argument( "Point-cloud symmetric tensor must "
                                         "have dimension one, two or three" );
        }

        static constexpr std::array< std::array< std::size_t, 2 >, 6 > indices {
            { { { 0, 0 } }, { { 1, 1 } }, { { 2, 2 } }, { { 0, 1 } }, { { 1, 2 } }, { { 0, 2 } } }
        };

        write_point_cloud_field( name,
                                 FieldType::symmetric_tensor,
                                 step,
                                 time,
                                 point_cloud_id,
                                 values.size(),
                                 [&]( const std::size_t point_id, const std::size_t component ) {
                                     const std::size_t row = indices[component][0];

                                     const std::size_t column = indices[component][1];

                                     if ( row >= rows || column >= columns ) {
                                         return 0.0;
                                     }

                                     return detail::matrix_component(
                                         values[point_id], row, column );
                                 } );
    }

    template < typename C >
    void
    write_point_cloud_tensor( const std::string &n, std::int64_t part, const C &v ) {
        auto s = active();
        write_point_cloud_tensor( n, s.first, s.second, part, v );
    }

    template < typename C >
    void
    write_point_cloud_tensor( const std::string &name,
                              const std::size_t step,
                              const double time,
                              const std::int64_t point_cloud_id,
                              const C &values ) {
        const auto &cloud = point_cloud( point_cloud_id );

        if ( values.size() != cloud.number_of_points() ) {
            throw std::invalid_argument( "Point-cloud tensor field size mismatch for '" + name +
                                         "': got " + std::to_string( values.size() ) +
                                         " values, expected " +
                                         std::to_string( cloud.number_of_points() ) );
        }

        if ( values.empty() ) {
            throw std::invalid_argument( "Cannot write an empty point-cloud tensor field" );
        }

        const std::size_t rows = detail::matrix_rows( values.front() );

        const std::size_t columns = detail::matrix_cols( values.front() );

        /*
         * EnSight tensors are embedded in a three-dimensional tensor.
         *
         * Valid input dimensions are therefore:
         *
         * 1 x 1
         * 2 x 2
         * 3 x 3
         */
        if ( rows == 0 || columns == 0 || rows != columns || rows > 3 ) {
            throw std::invalid_argument( "Point-cloud tensor field '" + name +
                                         "' must contain square matrices "
                                         "of dimension one, two or three" );
        }

        /*
         * write_point_cloud_field_streaming() uses row-major component
         * identifiers:
         *
         * 0 = xx, 1 = xy, 2 = xz
         * 3 = yx, 4 = yy, 5 = yz
         * 6 = zx, 7 = zy, 8 = zz
         *
         * It is responsible for converting this ordering to the
         * EnSight asymmetric tensor ordering.
         */
        write_point_cloud_field( name,
                                 FieldType::tensor,
                                 step,
                                 time,
                                 point_cloud_id,
                                 values.size(),
                                 [&]( const std::size_t point_id, const std::size_t component ) {
                                     const std::size_t row = component / 3;

                                     const std::size_t column = component % 3;

                                     /*
                                      * Embed one-dimensional or two-dimensional tensors
                                      * into a three-dimensional tensor by padding the
                                      * missing components with zeros.
                                      */
                                     if ( row >= rows || column >= columns ) {
                                         return 0.0;
                                     }

                                     return detail::matrix_component(
                                         values[point_id], row, column );
                                 } );
    }

    void
    validate_and_write_case() {
        require_mesh();
        if ( m_current_step )
            throw std::logic_error( "Call end_step() before validation" );
        for ( const auto &e : m_fields ) {
            if ( e.second.written_steps.size() != m_times.size() )
                throw std::runtime_error( "Field '" + e.second.name + "' is incomplete" );
        }
        update_case_file();
    }

  private:
    void
    require_mesh() const {
        if ( !m_mesh_written )
            throw std::logic_error( "write_mesh() must be called first" );
    }

    void
    require_active_step() const {
        if ( !m_current_step )
            throw std::logic_error( "No active step" );
    }

    std::pair< std::size_t, double >
    active() const {
        return { current_step(), current_time() };
    }

    void
    register_time( std::size_t s, double t ) {
        require_mesh();
        if ( m_current_step && *m_current_step != s )
            throw std::logic_error( "Explicit step differs from active step" );
        auto it = m_times.find( s );
        if ( it == m_times.end() ) {
            if ( s != m_times.size() )
                throw std::runtime_error( "Steps must be contiguous" );
            m_times.emplace( s, t );
        } else {
            const auto scale = std::max( { 1.0, std::abs( it->second ), std::abs( t ) } );
            if ( std::abs( it->second - t ) > 1e-14 * scale )
                throw std::runtime_error( "Time mismatch" );
        }
    }

    FieldDescription &
    register_field( const std::string &n,
                    FieldType type,
                    FieldLocation loc,
                    std::int64_t part,
                    std::size_t step ) {
        const auto safe = safe_filename( n );
        auto r = m_fields.emplace( safe, FieldDescription { n, safe, type, loc, part, {} } );
        auto &f = r.first->second;
        if ( !r.second &&
             ( f.type != type || f.location != loc || f.part_id != part || f.name != n ) )
            throw std::runtime_error( "Field definition changed: " + n );
        if ( f.written_steps.count( step ) )
            throw std::runtime_error( "Field already written: " + n );
        return f;
    }

    std::filesystem::path
    field_path( const FieldDescription &f, std::size_t step ) const {
        const auto dir = m_output_directory / f.safe_name;
        std::filesystem::create_directories( dir );
        std::ostringstream name;
        name << f.safe_name << '.' << std::setw( static_cast< int >( m_filename_width ) )
             << std::setfill( '0' ) << step;
        return dir / name.str();
    }

    template < typename G >
    void
    write_node_field( const std::string &n,
                      FieldType type,
                      std::size_t s,
                      double t,
                      std::int64_t part,
                      std::size_t count,
                      G &&getter ) {
        register_time( s, t );
        auto &f = register_field( n, type, FieldLocation::node, part, s );
        write_nodal_field_streaming(
            field_path( f, s ), n, type, count, std::forward< G >( getter ), part );
        f.written_steps.insert( s );
        if ( !m_current_step )
            update_case_file();
    }

    template < typename C >
    void
    write_matrix_node_field( const std::string &n,
                             FieldType type,
                             std::size_t s,
                             double t,
                             std::int64_t part,
                             std::size_t count,
                             const C &v,
                             bool symmetric ) {
        if ( v.size() != count || v.empty() )
            throw std::invalid_argument( "Matrix field size mismatch" );
        const auto rows = detail::matrix_rows( v[0] );
        const auto cols = detail::matrix_cols( v[0] );
        if ( rows < 1 || cols < 1 || rows > 3 || cols > 3 || ( symmetric && rows != cols ) )
            throw std::invalid_argument( "Invalid matrix dimensions" );

        if ( symmetric ) {
            write_node_field( n, type, s, t, part, count, [&]( auto i, auto c ) {
                static constexpr std::array< std::array< std::size_t, 2 >, 6 > indices {
                    { { { 0, 0 } },
                      { { 1, 1 } },
                      { { 2, 2 } },
                      { { 0, 1 } },
                      { { 1, 2 } },
                      { { 0, 2 } } }
                };
                const auto r = indices[c][0];
                const auto c2 = indices[c][1];
                return ( r < rows && c2 < cols ) ? detail::matrix_component( v[i], r, c2 ) : 0.0;
            } );
        } else {
            write_node_field( n, type, s, t, part, count, [&]( auto i, auto c ) {
                const auto r = c / 3;
                const auto c2 = c % 3;
                return ( r < rows && c2 < cols ) ? detail::matrix_component( v[i], r, c2 ) : 0.0;
            } );
        }
    }

    template < typename Getter >
    void
    write_point_cloud_field( const std::string &name,
                             const FieldType type,
                             const std::size_t step,
                             const double time,
                             const std::int64_t point_cloud_id,
                             const std::size_t value_count,
                             Getter &&getter ) {
        require_mesh();

        register_time( step, time );

        const auto &cloud = point_cloud( point_cloud_id );

        if ( value_count != cloud.number_of_points() ) {
            throw std::invalid_argument( "Point cloud field '" + name + "' has " +
                                         std::to_string( value_count ) +
                                         " values, but the point cloud contains " +
                                         std::to_string( cloud.number_of_points() ) + " points" );
        }

        /*
         * The point_cloud_id is retained in FieldDescription to verify
         * that the field remains attached to the same point cloud.
         *
         * It is not written as an EnSight part identifier because a
         * point cloud geometry is not a standard model part.
         */
        auto &field = register_field( name, type, FieldLocation::node, point_cloud_id, step );

        write_point_cloud_field_streaming( field_path( field, step ),
                                           name,
                                           type,
                                           cloud.number_of_points(),
                                           std::forward< Getter >( getter ) );

        field.written_steps.insert( step );

        if ( !m_current_step ) {
            update_case_file();
        }
    }

    void
    update_case_file() const {
        std::map< std::string, FieldDescription > model_fields;

        std::map< std::string, FieldDescription > point_cloud_fields;

        /*
         * Separate fields according to their logical geometry identifier.
         */
        for ( const auto &[field_name, field] : m_fields ) {
            if ( field.part_id == m_geometry.part_id ) {
                model_fields.emplace( field_name, field );
            } else {
                point_cloud_fields.emplace( field_name, field );
            }
        }

        /*
         * Main HHO mesh case.
         */
        write_case_file( m_output_directory / m_case_filename,
                         m_geometry_filename,
                         model_fields,
                         m_times,
                         m_filename_width );

        /*
         * Standalone Gauss-point case.
         */
        if ( !m_point_clouds.empty() ) {

            const std::map< std::size_t, double > empty_times;

            write_case_file( m_output_directory / m_point_cloud_case_filename,
                             m_point_cloud_geometry_filename,
                             point_cloud_fields,
                             point_cloud_fields.empty() ? empty_times : m_times,
                             m_filename_width );
        }
    }

    std::filesystem::path m_output_directory, m_case_filename, m_geometry_filename;
    std::filesystem::path m_point_cloud_geometry_filename = "gauss_points.geo";
    std::filesystem::path m_point_cloud_case_filename = "gauss_points.case";
    Geometry m_geometry;
    std::vector< PointCloud > m_point_clouds;
    std::map< std::string, FieldDescription > m_fields;
    std::map< std::size_t, double > m_times;
    std::optional< std::size_t > m_current_step;
    std::size_t m_filename_width = 6;
    bool m_mesh_written = false;
};
} // namespace disk::output::ensight
