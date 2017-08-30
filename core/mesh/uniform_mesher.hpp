/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016,2017 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#include <iostream>
#include <iomanip>
#include <regex>

#include "core/mesh/mesh.hpp"
#include "core/geometry/element_polytopes.hpp"
#include "core/geometry/geometry_poly.hpp"
#include "contrib/timecounter.h"

namespace disk {

namespace priv {

struct hierarchy_info
{
    size_t parent_id;
};

struct mesh_refinement_info
{
    point_id_type           midpoint_id;
    bool                    is_broken;

    mesh_refinement_info() : midpoint_id(0), is_broken(false) {}
};

struct mesher_triangular_storage_class
{
    typedef triangle<0, priv::hierarchy_info>       surface_type;
    typedef edge<1, priv::mesh_refinement_info>     edge_type;
    typedef node<2, void>                           node_type;
};

std::ostream&
operator<<(std::ostream& os, const edge<1, priv::mesh_refinement_info>& e)
{
    auto ptids = e.point_identifiers();
    assert(ptids.size() == 2);
    os << ptids[0] << " " << ptids[1] << " - ";
    os << "Midpoint id: " << e.user_data.midpoint_id;
    os << ", broken: " << std::boolalpha << e.user_data.is_broken;
    auto b = e.boundary_id();
    if (b.second)
        os << " [B " << b.first << "]";
    return os;
}
    
template<typename T>
using mesher_triangular_storage = mesh_storage<T, 2, mesher_triangular_storage_class>;

template<typename T>
using mesher_triangular_mesh = mesh_v2::mesh<2, mesher_triangular_storage<T>>;

} // namespace priv

template<typename Mesh>
class uniform_mesher;

template<typename T>
class uniform_mesher< mesh_v2::triangular_mesh<T> >
{
    typedef mesh_v2::triangular_mesh<T>         user_mesh_type;
    typedef priv::mesher_triangular_mesh<T>     internal_mesh_type;

    internal_mesh_type                          internal_mesh;

    void make_first_mesh()
    {
        auto be = internal_mesh.backend_storage();
        be->points.push_back({0.0, 0.0});
        be->points.push_back({1.0, 0.0});
        be->points.push_back({1.0, 1.0});
        be->points.push_back({0.0, 1.0});
        be->points.push_back({0.5, 0.5});

        typedef typename internal_mesh_type::node_type   nodet;
        be->nodes.push_back( nodet({0}) );
        be->nodes.push_back( nodet({1}) );
        be->nodes.push_back( nodet({2}) );
        be->nodes.push_back( nodet({3}) );
        be->nodes.push_back( nodet({4}) );

        typedef typename internal_mesh_type::edge_type   edget;
        be->edges.push_back( edget({0,1}, 100) );
        be->edges.push_back( edget({0,3}, 400) );
        be->edges.push_back( edget({0,4}) );
        be->edges.push_back( edget({1,2}, 200) );
        be->edges.push_back( edget({1,4}) );
        be->edges.push_back( edget({2,3}, 300) );
        be->edges.push_back( edget({2,4}) );
        be->edges.push_back( edget({3,4}) );

        typedef typename internal_mesh_type::surface_type surft;
        be->surfaces.push_back( surft({0,1,4}, 10) );
        be->surfaces.push_back( surft({0,4,3}, 20) );
        be->surfaces.push_back( surft({1,2,4}, 30) );
        be->surfaces.push_back( surft({2,3,4}, 40) );
    }

    void refine(void)
    {
        typedef typename internal_mesh_type::edge_type        edget;
        typedef typename internal_mesh_type::surface_type     surft;

        auto be = internal_mesh.backend_storage();

        for (auto& e : be->edges)
        {
            auto ptids = e.point_identifiers();
            assert(ptids.size() == 2);
            assert(ptids[0] < ptids[1]);
            assert(ptids[1] < be->points.size());
            auto bar = (be->points[ptids[0]] + be->points[ptids[1]]) * 0.5;
            e.user_data.is_broken = true;
            e.user_data.midpoint_id = be->points.size();
            be->points.push_back(bar);
        }

        size_t triangles_to_process = be->surfaces.size();
        auto num_initial_edges = be->edges.size();
        for (size_t i = 0; i < triangles_to_process; i++)
        {
            assert( i < be->surfaces.size() );
            auto t = be->surfaces[i];
            auto ptids = t.point_identifiers();
            assert(ptids.size() == 3);

            auto e0_a = ptids[0]; auto e0_b = ptids[1];
            edget e0({e0_a, e0_b}, node_ordering::SORTED);

            auto e1_a = ptids[1]; auto e1_b = ptids[2];
            edget e1({e1_a, e1_b}, node_ordering::SORTED);

            auto e2_a = ptids[2]; auto e2_b = ptids[0];
            edget e2({e2_a, e2_b}, node_ordering::SORTED);

            /* find the actual edges of the triangle */
            auto edges_end = std::next(be->edges.begin(), num_initial_edges);
            auto t_e0_itor = std::lower_bound(be->edges.begin(), edges_end, e0);
            assert(t_e0_itor != be->edges.end());
            auto t_e0 = *t_e0_itor;
            assert(t_e0.user_data.is_broken);
            auto e0_m = t_e0.user_data.midpoint_id;

            auto t_e1 = *std::lower_bound(be->edges.begin(), edges_end, e1);
            /* ADD ASSERT */
            assert(t_e1.user_data.is_broken);
            auto e1_m = t_e1.user_data.midpoint_id;

            auto t_e2 = *std::lower_bound(be->edges.begin(), edges_end, e2);
            /* ADD ASSERT */
            assert(t_e2.user_data.is_broken);
            auto e2_m = t_e2.user_data.midpoint_id;

            /* compute the edges of the new four triangles */
            /* first triangle */
            be->edges.push_back( edget({e0_a, e0_m}, t_e0.boundary_id(), node_ordering::SORTED) );
            be->edges.push_back( edget({e0_m, e2_m}, node_ordering::SORTED) );
            be->edges.push_back( edget({e2_m, e0_a}, t_e2.boundary_id(), node_ordering::SORTED) );
            auto t1 = surft( {e0_a, e0_m, e2_m}, t.domain_id() );
            t1.user_data.parent_id = i;
            be->surfaces.push_back( t1 );

            /* second triangle */
            be->edges.push_back( edget({e1_a, e1_m}, t_e1.boundary_id(), node_ordering::SORTED) );
            be->edges.push_back( edget({e1_m, e0_m}, node_ordering::SORTED) );
            be->edges.push_back( edget({e0_m, e1_a}, t_e0.boundary_id(), node_ordering::SORTED) );
            auto t2 = surft( {e1_a, e1_m, e0_m}, t.domain_id() );
            t2.user_data.parent_id = i;
            be->surfaces.push_back( t2 );

            /* third triangle */
            be->edges.push_back( edget({e2_a, e2_m}, t_e2.boundary_id(), node_ordering::SORTED) );
            be->edges.push_back( edget({e2_m, e1_m}, node_ordering::SORTED) );
            be->edges.push_back( edget({e1_m, e2_a}, t_e1.boundary_id(), node_ordering::SORTED) );
            auto t3 = surft( {e2_a, e2_m, e1_m}, t.domain_id() );
            t3.user_data.parent_id = i;
            be->surfaces.push_back( t3 );

            /* fourth triangle */
            be->edges.push_back( edget({e0_m, e1_m}, node_ordering::SORTED) );
            be->edges.push_back( edget({e1_m, e2_m}, node_ordering::SORTED) );
            be->edges.push_back( edget({e2_m, e0_m}, node_ordering::SORTED) );
            auto t4 = surft( {e0_m, e1_m, e2_m}, t.domain_id() );
            t4.user_data.parent_id = i;
            be->surfaces.push_back( t4 );
        }

        be->surfaces.erase( be->surfaces.begin(), be->surfaces.begin()+triangles_to_process );
        std::sort(be->surfaces.begin(), be->surfaces.end());

        auto edge_is_broken = [](const edget& e) -> bool { return e.user_data.is_broken; };
        std::remove_if(be->edges.begin(), be->edges.end(), edge_is_broken);
        std::sort(be->edges.begin(), be->edges.end());
        auto last = std::unique(be->edges.begin(), be->edges.end());
        be->edges.erase(last, be->edges.end());

        be->nodes.clear();
    }

    user_mesh_type convert_to_user() const
    {
        user_mesh_type user_mesh;

        auto u_be = user_mesh.backend_storage();
        auto i_be = internal_mesh.backend_storage();

        u_be->points.resize( i_be->points.size() );
        u_be->nodes.resize( i_be->nodes.size() );
        u_be->edges.resize( i_be->edges.size() );
        u_be->surfaces.resize( i_be->surfaces.size() );

        std::copy( i_be->points.begin(), i_be->points.end(), u_be->points.begin() );
        std::copy( i_be->nodes.begin(), i_be->nodes.end(), u_be->nodes.begin() );

        typedef typename internal_mesh_type::edge_type  iet;
        typedef typename user_mesh_type::edge_type      uet;

        auto edge_cast = [](const iet& e) -> uet {
            return uet(e);
        };

        std::transform(i_be->edges.begin(), i_be->edges.end(), u_be->edges.begin(), edge_cast);

        typedef typename internal_mesh_type::surface_type   ist;
        typedef typename user_mesh_type::surface_type       ust;

        auto surface_cast = [](const ist& e) -> ust {
            return ust(e);
        };

        std::transform(i_be->surfaces.begin(), i_be->surfaces.end(), u_be->surfaces.begin(), surface_cast);

        return user_mesh;
    }

public:
    uniform_mesher()
    {
        make_first_mesh();
        dump_to_matlab(internal_mesh, "first.m");
        refine();
        dump_to_matlab(internal_mesh, "second.m");
        refine();
        dump_to_matlab(internal_mesh, "third.m");

        auto m = convert_to_user();

        dump_to_matlab(m, "user.m");


        //timecounter tc;
        //tc.tic();
        //for (size_t i = 0; i < 5; i++)
        //    refine();
        //tc.toc();
        //std::cout << internal_mesh.cells_size() << " " << tc << std::endl;

        #if 0
        for (auto& cl : internal_mesh)
        {
            std::cout << cl << " -> " << measure(internal_mesh, cl) << std::endl;
        }
        
        for (auto itor = internal_mesh.faces_begin(); itor != internal_mesh.faces_end(); itor++)
        {
            std::cout << *itor << std::endl;
            //std::cout << measure(internal_mesh, *itor) << std::endl;
        }
        #endif
        
    }
};

} // namespace disk