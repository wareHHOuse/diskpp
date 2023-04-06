/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */
/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018, 2019
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
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

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/geometry/geometry.hpp"

#include "meshgen_fvca5hex.hpp"

namespace disk {

/**
 * simple_mesher meshes the domain [0,1]^d with elements of the specified mesh.
 * To use a different bounding box, mesh the first time and then apply a
 * transform using the .transform() method.
 */

template<typename Mesh>
class simple_mesher
{
    static_assert(sizeof(Mesh) == -1, "simple_mesher does not have a specialization for the specified mesh type");
};

template<typename T>
class simple_mesher<triangular_mesh<T>>
{
    typedef triangular_mesh<T>                          mesh_type;
    static const size_t DIM = 2;
    typedef typename triangular_mesh<T>::storage_type   storage_type;
    typedef point<T,DIM>                                point_type;

    typedef typename triangular_mesh<T>::node_type      node_type;
    typedef typename triangular_mesh<T>::edge_type      edge_type;
    typedef typename triangular_mesh<T>::surface_type   surface_type;

    std::shared_ptr<storage_type>   storage;

public:
    simple_mesher(mesh_type& msh)
        : storage( msh.backend_storage() )
    {
        /* Init the first level of the mesh */
        storage->points.push_back( point_type(0.0, 0.0) );
        auto pi0 = point_identifier<2>(0);
        storage->nodes.push_back( node_type( {pi0} ) );

        storage->points.push_back( point_type(1.0, 0.0) );
        auto pi1 = point_identifier<2>(1);
        storage->nodes.push_back( node_type( {pi1} ) );

        storage->points.push_back( point_type(1.0, 1.0) );
        auto pi2 = point_identifier<2>(2);
        storage->nodes.push_back( node_type( {pi2} ) );

        storage->points.push_back( point_type(0.0, 1.0) );
        auto pi3 = point_identifier<2>(3);
        storage->nodes.push_back( node_type( {pi3} ) );

        storage->points.push_back( point_type(0.5, 0.5) );
        auto pi4 = point_identifier<2>(4);
        storage->nodes.push_back( node_type( {pi4} ) );

        storage->edges.push_back( edge_type({pi0, pi1}) ); //0 b0
        storage->edges.push_back( edge_type({pi0, pi3}) ); //1 b3
        storage->edges.push_back( edge_type({pi0, pi4}) ); //2
        storage->edges.push_back( edge_type({pi1, pi2}) ); //3 b1
        storage->edges.push_back( edge_type({pi1, pi4}) ); //4
        storage->edges.push_back( edge_type({pi2, pi3}) ); //5 b2
        storage->edges.push_back( edge_type({pi2, pi4}) ); //6
        storage->edges.push_back( edge_type({pi3, pi4}) ); //7

        storage->boundary_info.resize(8); /* Total 8 edges in the mesh */
        storage->boundary_info[0] = boundary_descriptor(0, true);
        storage->boundary_info[1] = boundary_descriptor(3, true);
        storage->boundary_info[3] = boundary_descriptor(1, true);
        storage->boundary_info[5] = boundary_descriptor(2, true);

        storage->surfaces.push_back( surface_type({pi0, pi1, pi4}) );
        storage->surfaces.push_back( surface_type({pi0, pi3, pi4}) );
        storage->surfaces.push_back( surface_type({pi1, pi2, pi4}) );
        storage->surfaces.push_back( surface_type({pi2, pi3, pi4}) );
    }

    void refine(void)
    {
        size_t node_offset = storage->nodes.size();

        typedef std::pair<edge_type, boundary_descriptor> ne_pair;

        std::vector<ne_pair> new_edges;
        new_edges.reserve( storage->edges.size()*2 );

        std::vector<surface_type> new_surfaces;
        new_surfaces.reserve( storage->surfaces.size()*4 );

        size_t edge_offset = 0;
        for (auto& e : storage->edges)
        {
            auto ptids = e.point_ids();

            assert(ptids.size() == 2);
            assert(ptids[0] < storage->points.size());
            assert(ptids[1] < storage->points.size());

            auto p0 = storage->points[ ptids[0] ];
            auto p1 = storage->points[ ptids[1] ];
            auto pm = (p0 + p1)/2.;

            storage->points.push_back(pm);
            auto pmi = point_identifier<2>( storage->nodes.size() );
            storage->nodes.push_back( node_type({pmi}) );

            assert( ptids[0] < pmi );
            assert( ptids[1] < pmi );

            /* Those are going to be boundary edges */
            assert(edge_offset < storage->boundary_info.size());
            auto ep1 = std::make_pair(edge_type({ptids[0], pmi}), storage->boundary_info[edge_offset]);
            auto ep2 = std::make_pair(edge_type({ptids[1], pmi}), storage->boundary_info[edge_offset]);
            new_edges.push_back( ep1 );
            new_edges.push_back( ep2 );

            edge_offset++;
        }

        for (auto& s : storage->surfaces)
        {
            auto ptids = s.point_ids();
            assert(ptids.size() == 3);
            assert(ptids[0] < storage->points.size());
            assert(ptids[1] < storage->points.size());
            assert(ptids[2] < storage->points.size());

            auto eofs = [&](const edge_type& e) -> auto {
                //std::cout << "Searching for edge " << e << std::endl;
                auto be = begin(storage->edges);
                auto ee = end(storage->edges);
                auto ei = std::lower_bound(be, ee, e);
                if (ei == ee or e != *ei)
                    throw std::logic_error("Edge not found. This is a bug.");
                return point_identifier<2>(std::distance(be, ei) + node_offset);
            };

            auto p_e0 = eofs( edge_type({ptids[0], ptids[1]}) );
            auto p_e1 = eofs( edge_type({ptids[1], ptids[2]}) );
            auto p_e2 = eofs( edge_type({ptids[0], ptids[2]}) );

            /* Those are going to be internal edges */
            new_edges.push_back( std::make_pair(edge_type({p_e0, p_e1}), boundary_descriptor()) );
            new_edges.push_back( std::make_pair(edge_type({p_e2, p_e1}), boundary_descriptor()) );
            new_edges.push_back( std::make_pair(edge_type({p_e0, p_e2}), boundary_descriptor()) );

            new_surfaces.push_back( surface_type( {ptids[0], p_e0, p_e2} ) );
            new_surfaces.push_back( surface_type( {ptids[1], p_e0, p_e1} ) );
            new_surfaces.push_back( surface_type( {ptids[2], p_e2, p_e1} ) );
            new_surfaces.push_back( surface_type( {p_e0, p_e2, p_e1} ) );
        }

        auto comp = [](const ne_pair& nep1, const ne_pair& nep2) -> bool {
            return nep1.first < nep2.first;
        };
        std::sort(new_edges.begin(), new_edges.end(), comp);

        storage->edges.clear(); storage->edges.reserve( new_edges.size() );
        storage->boundary_info.clear(); storage->boundary_info.reserve( new_edges.size() );
        for (auto& ne : new_edges)
        {
            storage->edges.push_back(ne.first);
            storage->boundary_info.push_back(ne.second);
        }

        std::sort(new_surfaces.begin(), new_surfaces.end());
        std::swap(storage->surfaces, new_surfaces);
    }
};

template<typename Mesh>
auto make_simple_mesher(Mesh& msh)
{
    return simple_mesher<Mesh>(msh);
}


template<typename Mesh>
class fvca5_hex_mesher
{
    static_assert(sizeof(Mesh) == -1, "fvca5_hex_mesher does not have a specialization for the specified mesh type");
};

template<typename T>
class fvca5_hex_mesher<generic_mesh<T,2>>
{
    typedef generic_mesh<T,2>                   mesh_type;
    static const size_t DIM = 2;

    typedef typename mesh_type::storage_type    storage_type;
    typedef typename mesh_type::node_type       node_type;
    typedef typename mesh_type::edge_type       edge_type;
    typedef typename mesh_type::surface_type    surface_type;
    typedef typename mesh_type::point_type      point_type;

    std::shared_ptr<storage_type>               storage;

public:
    fvca5_hex_mesher(mesh_type& msh)
        : storage( msh.backend_storage() )
    {
    }

    void make_level(size_t level)
    {
        fvca5::hexmesh hm;
        fvca5::make_hexmesh(hm, level);

        /* Points */
        storage->points.resize( hm.points.size() );
        auto point_transform = [&](const fvca5::point& pt) -> point_type {
            return point_type( T(pt.x)/hm.coord_max, T(pt.y)/hm.coord_max );
        };
        std::transform(hm.points.begin(), hm.points.end(), storage->points.begin(), point_transform);

        /* Nodes */
        std::vector<node_type> nodes(hm.points.size());
        for (size_t i = 0; i < nodes.size(); i++)
            nodes[i] = node_type(point_identifier<2>(i));
        std::swap(storage->nodes, nodes);

        /* Edges */
        storage->edges.resize(hm.edges.size());
        auto edge_transform = [&](const fvca5::edge& edge) -> edge_type {
            auto ea_ofs = fvca5::offset(hm, edge.a, false);
            auto eb_ofs = fvca5::offset(hm, edge.b, false);
            assert(ea_ofs != eb_ofs);
            auto node_a = typename node_type::id_type(ea_ofs);
            auto node_b = typename node_type::id_type(eb_ofs);
            return edge_type(node_a, node_b);
        };
        std::transform(hm.edges.begin(), hm.edges.end(), storage->edges.begin(), edge_transform);
        std::sort(storage->edges.begin(), storage->edges.end());

        /* Boundary edges */
        storage->boundary_info.resize( storage->edges.size() );
        for (const auto& [boundary, edges] : hm.boundary_edges)
        {
            for (const auto& edge : edges)
            {
                auto e = edge_transform(edge);
                auto position = find_element_id(storage->edges.begin(), storage->edges.end(), e);
                if (position.first == false)
                    throw std::logic_error("edge not found, this is a bug.");

                boundary_descriptor bi(0, true);
                storage->boundary_info.at(position.second) = bi;
            }
        }

        /* Surfaces */
        using sd_pair = std::pair<surface_type, size_t>;
        std::vector< sd_pair > surfaces;
        surfaces.reserve( hm.polygons.size() );
        for (const auto& polygon : hm.polygons)
        {
            std::vector<size_t> surface_nodes;
            std::vector<typename edge_type::id_type> surface_edges;
            for (size_t i = 0; i < polygon.points.size(); i++)
            {
                auto pa = polygon.points[i];
                auto pb = polygon.points[(i+1)%polygon.points.size()];
                auto edge = fvca5::edge(pa, pb);
                auto e = edge_transform(edge);
                auto position = find_element_id(storage->edges.begin(), storage->edges.end(), e);
                if (position.first == false)
                    throw std::logic_error("edge not found, this is a bug.");
                surface_edges.push_back(position.second);

                auto ea_ofs = fvca5::offset(hm, pa, false);
                surface_nodes.push_back(ea_ofs);
            }
            
            auto domain_id = 0;
            auto surface = surface_type(surface_edges);
            surface.set_point_ids(surface_nodes.begin(), surface_nodes.end()); /* XXX: crap */
            surfaces.push_back( std::make_pair(surface, domain_id) );
        }

        auto comp = [](const sd_pair& a, const sd_pair& b) -> bool {
            return a.first < b.first;
        };

        std::sort(surfaces.begin(), surfaces.end(), comp);

        auto get_surf = [](const sd_pair& sdp) -> auto { return sdp.first; };
        storage->surfaces.resize( surfaces.size() );
        std::transform(surfaces.begin(), surfaces.end(), storage->surfaces.begin(), get_surf);

        auto get_id = [](const sd_pair& sdp) -> auto { return subdomain_descriptor(sdp.second); };
        storage->subdomain_info.resize( surfaces.size() );
        std::transform(surfaces.begin(), surfaces.end(), storage->subdomain_info.begin(), get_id);

        //mark_internal_faces(msh);
    }
};

template<typename Mesh>
auto make_fvca5_hex_mesher(Mesh& msh)
{
    return fvca5_hex_mesher<Mesh>(msh);
}

template<typename T>
void make_single_element_mesh(generic_mesh<T,1>& msh, const T& a, const T& b)
{
    using mesh_type = generic_mesh<T,1>;
    using point_type = typename mesh_type::point_type;
    using node_type = typename mesh_type::node_type;
    using edge_type = typename mesh_type::edge_type;
    using nodeid_type = typename node_type::id_type;

    auto storage = msh.backend_storage();

    auto num_edges = 1;
    auto num_nodes = 2;

    storage->points.resize(num_nodes);
    storage->nodes.resize(num_nodes);
    storage->edges.resize(num_edges);

    storage->points.at(0)   = point_type(a);
    storage->points.at(1) = point_type(b);

    auto n0 = typename node_type::id_type(0);
    auto n1 = typename node_type::id_type(1);

    storage->nodes.at(0) = node_type(n0);
    storage->nodes.at(1) = node_type(n1);

    auto e = edge_type(n0, n1);
    storage->edges.at(0) = e;

    storage->boundary_info.resize(num_nodes);
    boundary_descriptor bi(0, true);
    storage->boundary_info.at(0) = bi;
    storage->boundary_info.at(1) = bi;    
}

template<typename T>
void make_single_element_mesh(simplicial_mesh<T,2>& msh, const point<T,2>& a, const point<T,2>& b,
    const point<T,2>& c)
{
    using mesh_type = simplicial_mesh<T,2>;
    using point_type = typename mesh_type::point_type;
    using node_type = typename mesh_type::node_type;
    using edge_type = typename mesh_type::edge_type;
    using nodeid_type = typename node_type::id_type;

    auto storage = msh.backend_storage();

    storage->points.push_back( a );
    storage->points.push_back( b );
    storage->points.push_back( c );

    auto p0 = point_identifier<2>(0);
    auto p1 = point_identifier<2>(1);
    auto p2 = point_identifier<2>(2);

    storage->nodes.push_back( { p0 } );
    storage->nodes.push_back( { p1 } );
    storage->nodes.push_back( { p2 } );

    storage->edges.push_back( {p0, p1} );
    storage->edges.push_back( {p0, p2} );
    storage->edges.push_back( {p1, p2} );

    storage->surfaces.push_back( {p0, p1, p2} );
    storage->subdomain_info.resize( 1 );
}

template<typename T>
void make_single_element_mesh(cartesian_mesh<T,2>& msh, const point<T,2>& base, const T& hx, const T& hy)
{
    using mesh_type = cartesian_mesh<T,2>;
    using point_type = typename mesh_type::point_type;
    using node_type = typename mesh_type::node_type;
    using edge_type = typename mesh_type::edge_type;
    using nodeid_type = typename node_type::id_type;

    auto storage = msh.backend_storage();

    storage->points.push_back( point_type(    base.x(),    base.y() ) );
    storage->points.push_back( point_type( base.x()+hx,    base.y() ) );
    storage->points.push_back( point_type( base.x()+hx, base.y()+hy ) );
    storage->points.push_back( point_type(    base.x(), base.y()+hy ) );

    auto p0 = point_identifier<2>(0);
    auto p1 = point_identifier<2>(1);
    auto p2 = point_identifier<2>(2);
    auto p3 = point_identifier<2>(3);

    storage->nodes.push_back( { p0 } );
    storage->nodes.push_back( { p1 } );
    storage->nodes.push_back( { p2 } );
    storage->nodes.push_back( { p3 } );

    storage->edges.push_back( {p0, p1} );
    storage->edges.push_back( {p0, p3} );
    storage->edges.push_back( {p1, p2} );
    storage->edges.push_back( {p2, p3} );

    storage->surfaces.push_back( {p0, p1, p2, p3} );
    storage->subdomain_info.resize( 1 );
}

} //namespace disk