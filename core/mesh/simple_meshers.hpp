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

#include "geometry/geometry.hpp"

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
        storage->boundary_info[0] = bnd_info(0, true);
        storage->boundary_info[1] = bnd_info(3, true);
        storage->boundary_info[3] = bnd_info(1, true);
        storage->boundary_info[5] = bnd_info(2, true);

        storage->surfaces.push_back( surface_type({pi0, pi1, pi4}) );
        storage->surfaces.push_back( surface_type({pi0, pi3, pi4}) );
        storage->surfaces.push_back( surface_type({pi1, pi2, pi4}) );
        storage->surfaces.push_back( surface_type({pi2, pi3, pi4}) );
    }

    void refine(void)
    {
        size_t node_offset = storage->nodes.size();

        typedef std::pair<edge_type, bnd_info> ne_pair;

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
                if (ei == ee)
                    throw std::logic_error("Edge not found. This is a bug.");
                return point_identifier<2>(std::distance(be, ei) + node_offset);
            };

            auto p_e0 = eofs( edge_type({ptids[0], ptids[1]}) );
            auto p_e1 = eofs( edge_type({ptids[1], ptids[2]}) );
            auto p_e2 = eofs( edge_type({ptids[0], ptids[2]}) );

            /* Those are going to be internal edges */
            new_edges.push_back( std::make_pair(edge_type({p_e0, p_e1}), bnd_info()) );
            new_edges.push_back( std::make_pair(edge_type({p_e2, p_e1}), bnd_info()) );
            new_edges.push_back( std::make_pair(edge_type({p_e0, p_e2}), bnd_info()) );

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

} //namespace disk
