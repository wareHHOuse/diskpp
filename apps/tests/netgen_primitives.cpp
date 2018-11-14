/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
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

#include <algorithm>

#include "core/loaders/loader.hpp"
#include "geometry/geometry.hpp"
#include "common.hpp"

#include "mesh/rational.hpp"
#include "mesh/mesh_hierarchy.hpp"

#include "output/silo.hpp"

point<rational<int>,2>
barycenter(const std::vector<point<rational<int>,2>>& mp,
           size_t a, size_t b, size_t c)
{
    return (mp[a] + mp[b] + mp[c])/3;
}

rational<int>
area(const std::vector<point<rational<int>,2>>& mp,
     size_t a, size_t b, size_t c)
{
    return abs(det(mp[b] - mp[a], mp[c] - mp[a])/2);
}

template<typename Mesh>
void x(const Mesh& msh)
{
    using R = rational<int>;
    using P = point<R,2>;

    std::vector<P> mp;
    mp.resize(14);
    mp[0]  = P( R(0,10), R(0,10) );
    mp[1]  = P( R(1,10), R(0,10) );
    mp[2]  = P( R(3,10), R(0,10) );
    mp[3]  = P( R(6,10), R(0,10) );
    mp[4]  = P( R(6,10), R(2,10) );
    mp[5]  = P( R(4,10), R(3,10) );
    mp[6]  = P( R(2,10), R(5,10) );
    mp[7]  = P( R(6,10), R(5,10) );
    mp[8]  = P( R(4,10), R(9,20) );
    mp[9]  = P( R(0,10), R(10,10) );
    mp[10] = P( R(1,10), R(10,10) );
    mp[11] = P( R(2,10), R(10,10) );
    mp[12] = P( R(4,10), R(10,10) );
    mp[13] = P( R(6,10), R(10,10) );

    for (auto& cl : msh)
    {
        auto ptids = cl.point_ids();

        auto rbar = barycenter(mp, ptids[0], ptids[1], ptids[2]);
        auto rmeas = area(mp, ptids[0], ptids[1], ptids[2]);

        auto bar = barycenter(msh, cl);
        auto meas = measure(msh, cl);

        double ULP_max = 1;
        bool success = true;

        if ( !almost_equal(bar.x(), double(rbar.x()), ULP_max) )
        {
            std::cout << "  pt.x() not accurate: ";
            std::cout << std::setprecision(16) << bar.x() << " ";
            std::cout << rbar.x() << std::endl;
            success = false;
        }

        if ( !almost_equal(bar.y(), double(rbar.y()), ULP_max) )
        {
            std::cout << "  pt.y() not accurate: ";
            std::cout << std::setprecision(16) << bar.y() << " ";
            std::cout << rbar.y() << std::endl;
            success = false;
        }
    }
}

template<typename Mesh>
void
test(const Mesh& msh)
{
    using R          = typename Mesh::coordinate_type;
    using point_type = typename Mesh::point_type;

    R tot_meas = R(0);
    point_type tot_bar( R(0), R(0) );
    for (auto& cl : msh)
    {
        auto bar = barycenter(msh, cl);
        auto meas = measure(msh, cl);
        tot_bar = tot_bar + bar*meas;
        tot_meas += meas;
    }
    std::cout << std::setprecision(16) << tot_bar << " " << tot_meas << std::endl;
}

template<typename Mesh>
void
add_triangle(Mesh& msh, size_t p0, size_t p1, size_t p2)
{
    auto storage = msh.backend_storage();

    using triangle = typename Mesh::surface_type;
    using edge = typename Mesh::edge_type;
    point_identifier<2> pi0(p0);
    point_identifier<2> pi1(p1);
    point_identifier<2> pi2(p2);
    storage->surfaces.push_back( triangle({pi0, pi1, pi2}) );
    storage->edges.push_back( edge({pi0, pi1}) );
    storage->edges.push_back( edge({pi1, pi2}) );
    storage->edges.push_back( edge({pi0, pi2}) );
}

template<typename Mesh>
void
create_geometry(Mesh& msh)
{
    auto storage = msh.backend_storage();

    add_triangle(msh,0,1,3);
    add_triangle(msh,0,2,3);
    add_triangle(msh,1,3,4);
    add_triangle(msh,2,3,7);
    add_triangle(msh,3,4,5);
    add_triangle(msh,3,5,6);
    add_triangle(msh,3,6,10);
    add_triangle(msh,3,7,8);
    add_triangle(msh,3,8,9);
    add_triangle(msh,3,9,10);

    disk::priv::sort_uniq(storage->edges);
    disk::priv::sort_uniq(storage->surfaces);

    std::vector<typename Mesh::edge_type> boundary_edges;

    auto add_bedge = [&](size_t p0, size_t p1) -> void {
        using edge = typename Mesh::edge_type;
        point_identifier<2> pi0(p0);
        point_identifier<2> pi1(p1);
        boundary_edges.push_back( edge({pi0, pi1}) );
    };

    storage->boundary_info.resize( storage->edges.size() );
    for (auto& be : boundary_edges)
    {
        auto position = find_element_id(storage->edges.begin(),
                                        storage->edges.end(), be);
        if (position.first == false)
        {
            std::cout << "Bad bug at " << __FILE__ << "("
            << __LINE__ << ")" << std::endl;
            return;
        }
        disk::bnd_info bi{1, true};
        storage->boundary_info.at(position.second) = bi;
    }
}

template<typename Mesh>
void
add_points(Mesh& msh, const std::vector<typename Mesh::point_type>& pts)
{
    auto storage = msh.backend_storage();
    storage->points.resize(pts.size());
    std::copy(pts.begin(), pts.end(), storage->points.begin());
}

template<typename RMesh, typename FMesh>
void
compute_errors(disk::silo_database& silo, const std::string& prefix,
               const RMesh& rmsh, const FMesh& fmsh)
{
    using T = typename FMesh::coordinate_type;

    std::string meshname = prefix + "_mesh";

    silo.add_mesh(rmsh, meshname.c_str());

    assert( rmsh.cells_size() == fmsh.cells_size() );

    std::vector<T> meas_err, bar_x_err, bar_y_err;
    for (size_t i = 0; i < rmsh.cells_size(); i++)
    {
        auto rcl = *std::next(rmsh.cells_begin(), i);
        auto fcl = *std::next(fmsh.cells_begin(), i);

        auto rmeas = measure(rmsh, rcl);
        auto fmeas = measure(fmsh, fcl);

        auto rbar = barycenter(rmsh, rcl);
        auto fbar = barycenter(fmsh, fcl);

        meas_err.push_back( std::abs(fmeas - T(rmeas)) );
        bar_x_err.push_back( std::abs(fbar.x() - T(rbar.x())) );
        bar_y_err.push_back( std::abs(fbar.y() - T(rbar.y())) );
    }

    disk::silo_zonal_variable<T> me(prefix + "_meas_err", meas_err);
    disk::silo_zonal_variable<T> bxe(prefix + "_bar_x_err", bar_x_err);
    disk::silo_zonal_variable<T> bye(prefix + "_bar_y_err", bar_y_err);

    silo.add_variable(meshname.c_str(), me);
    silo.add_variable(meshname.c_str(), bxe);
    silo.add_variable(meshname.c_str(), bye);
}

int main(void)
{
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    
    using T             = double;
    using R             = rational<int64_t>;
    using fmesh_type    = disk::simplicial_mesh<T, 2>;
    using rmesh_type    = disk::simplicial_mesh<R, 2>;
    using fpoint_type   = typename fmesh_type::point_type;
    using rpoint_type   = typename rmesh_type::point_type;

    int64_t bnum = 74;
    int64_t bden = 100;

    std::vector<rpoint_type>    rmesh_points;
    rmesh_points.push_back( rpoint_type(R(0), R(0)) );
    rmesh_points.push_back( rpoint_type(R(1), R(0)) );
    rmesh_points.push_back( rpoint_type(R(0), R(2,5)) );
    rmesh_points.push_back( rpoint_type(R(bnum,bden), R(bnum,bden)) );
    rmesh_points.push_back( rpoint_type(R(1), R(1,2)) );
    rmesh_points.push_back( rpoint_type(R(1), R(3,4)) );
    rmesh_points.push_back( rpoint_type(R(1), R(7,8)) );
    rmesh_points.push_back( rpoint_type( R(0), R(1) ) );
    rmesh_points.push_back( rpoint_type( R(3,10), R(1) ) );
    rmesh_points.push_back( rpoint_type( R(1,2), R(1) ) );
    rmesh_points.push_back( rpoint_type(R(1), R(1)) );

    auto transform_point = [](const rpoint_type& rpt) -> fpoint_type {
        auto x = typename fpoint_type::value_type(rpt.x());
        auto y = typename fpoint_type::value_type(rpt.y());
        return fpoint_type( x, y );
    };

    std::vector<fpoint_type>    fmesh_points;
    fmesh_points.resize( rmesh_points.size() );
    std::transform( rmesh_points.begin(), rmesh_points.end(), fmesh_points.begin(), transform_point );

    fmesh_type fmsh;
    rmesh_type rmsh;

    add_points(fmsh, fmesh_points);
    add_points(rmsh, rmesh_points);

    create_geometry(fmsh);
    create_geometry(rmsh);

    test(rmsh);
    test(fmsh);

    disk::silo_database silo;
    silo.create("netgen_primitives.silo");

    compute_errors(silo, "level_0", rmsh, fmsh);

    size_t mesh_levels = 6; /* Don't push too far with rational<> */
    disk::mesh_hierarchy<R> rhier;
    rhier.build_hierarchy(rmsh, mesh_levels);

    disk::mesh_hierarchy<T> fhier;
    fhier.build_hierarchy(fmsh, mesh_levels);

    for(size_t i = 0; i < mesh_levels; i++)
    {
        std::stringstream ss;
        ss << "level_" << i+1;
        auto rm = *std::next(rhier.meshes_begin(), i);
        auto fm = *std::next(fhier.meshes_begin(), i);
        compute_errors(silo, ss.str(), rm, fm);
    }

    silo.close();

    return 0;
}
