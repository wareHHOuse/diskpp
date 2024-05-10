/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */
//* Karol Cascavita (C) 2023
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

#include <random>
#include <regex>

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

    void init_pattern_1(void)
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

        storage->subdomain_info.resize( storage->surfaces.size() );
    }

    void init_pattern_2(void)
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


        storage->edges.push_back( edge_type({pi0, pi1}) ); //0 b0
        storage->edges.push_back( edge_type({pi0, pi2}) ); //1
        storage->edges.push_back( edge_type({pi0, pi3}) ); //2 b3
        storage->edges.push_back( edge_type({pi1, pi2}) ); //3 b1
        storage->edges.push_back( edge_type({pi2, pi3}) ); //4 b2

        storage->boundary_info.resize(5); /* Total 5 edges in the mesh */
        storage->boundary_info[0] = boundary_descriptor(0, true);
        storage->boundary_info[2] = boundary_descriptor(3, true);
        storage->boundary_info[3] = boundary_descriptor(1, true);
        storage->boundary_info[4] = boundary_descriptor(2, true);

        storage->surfaces.push_back( surface_type({pi0, pi1, pi2}) );
        storage->surfaces.push_back( surface_type({pi0, pi2, pi3}) );

        storage->subdomain_info.resize( storage->surfaces.size() );
    }

public:
    simple_mesher(mesh_type& msh)
        : storage( msh.backend_storage() )
    {
        init_pattern_1();
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
        storage->boundary_info.clear();
        storage->boundary_info.reserve(new_edges.size());
        for (auto& ne : new_edges)
        {
            storage->edges.push_back(ne.first);
            storage->boundary_info.push_back(ne.second);
        }

        std::sort(new_surfaces.begin(), new_surfaces.end());
        std::swap(storage->surfaces, new_surfaces);
        storage->subdomain_info.resize( storage->surfaces.size() );
    }
};

template<typename T>
class simple_mesher<tetrahedral_mesh<T>>
{
    typedef tetrahedral_mesh<T>                             mesh_type;
    static const size_t DIM = 3;
    typedef typename tetrahedral_mesh<T>::storage_type      storage_type;
    typedef point<T, DIM>                                   point_type;

    typedef typename tetrahedral_mesh<T>::node_type         node_type;
    typedef typename tetrahedral_mesh<T>::edge_type         edge_type;
    typedef typename tetrahedral_mesh<T>::surface_type      surface_type;
    typedef typename tetrahedral_mesh<T>::volume_type       volume_type;
    typedef std::pair<surface_type, boundary_descriptor>    ns_pair;

    std::shared_ptr<storage_type>   storage;
    std::vector<ns_pair>            new_surfaces;
    std::vector<volume_type>        new_volumes;

    using ptid = point_identifier<3>;

    void make_tetra(const ptid& pi0, const ptid& pi1,
        const ptid& pi2, const ptid& pi3,
        const boundary_descriptor& bi0, const boundary_descriptor& bi1,
        const boundary_descriptor& bi2, const boundary_descriptor& bi3)
    {
        storage->edges.push_back( edge_type({pi0, pi1}) );
        storage->edges.push_back( edge_type({pi0, pi2}) );
        storage->edges.push_back( edge_type({pi0, pi3}) );
        storage->edges.push_back( edge_type({pi1, pi2}) );
        storage->edges.push_back( edge_type({pi1, pi3}) );
        storage->edges.push_back( edge_type({pi2, pi3}) );

        new_surfaces.push_back({surface_type({pi0, pi1, pi2}), bi0});
        new_surfaces.push_back({surface_type({pi0, pi1, pi3}), bi1});
        new_surfaces.push_back({surface_type({pi0, pi2, pi3}), bi2});
        new_surfaces.push_back({surface_type({pi1, pi2, pi3}), bi3});

        new_volumes.push_back(volume_type({ pi0, pi1, pi2, pi3 }));
    }

    public:

    simple_mesher(mesh_type& msh)
        : storage(msh.backend_storage())
    {
        /* Init the first level of the mesh */
        storage->points.push_back(point_type(0.0, 0.0, 0.0));
        auto pi0 = ptid(0);
        storage->nodes.push_back(node_type({ pi0 }));

        storage->points.push_back(point_type(1.0, 0.0, 0.0));
        auto pi1 = ptid(1);
        storage->nodes.push_back(node_type({ pi1 }));

        storage->points.push_back(point_type(1.0, 1.0, 0.0));
        auto pi2 = ptid(2);
        storage->nodes.push_back(node_type({ pi2 }));

        storage->points.push_back(point_type(0.0, 1.0, 0.0));
        auto pi3 = ptid(3);
        storage->nodes.push_back(node_type({ pi3 }));

        storage->points.push_back(point_type(0.0, 0.0, 1.0));
        auto pi4 = ptid(4);
        storage->nodes.push_back(node_type({ pi4 }));

        storage->points.push_back(point_type(1.0, 0.0, 1.0));
        auto pi5 = ptid(5);
        storage->nodes.push_back(node_type({ pi5 }));

        storage->points.push_back(point_type(1.0, 1.0, 1.0));
        auto pi6 = ptid(6);
        storage->nodes.push_back(node_type({ pi6 }));

        storage->points.push_back(point_type(0.0, 1.0, 1.0));
        auto pi7 = ptid(7);
        storage->nodes.push_back(node_type({ pi7 }));


        auto bn = boundary_descriptor(0, false);
        auto b0 = boundary_descriptor(0, true);
        auto b1 = boundary_descriptor(1, true);
        auto b2 = boundary_descriptor(2, true);
        auto b3 = boundary_descriptor(3, true);
        auto b4 = boundary_descriptor(4, true);
        auto b5 = boundary_descriptor(5, true);

        size_t num_surfaces = 6 * 4;
        new_surfaces.reserve(num_surfaces);

        make_tetra(pi0, pi1, pi3, pi7, b0, bn, b2, bn);
        make_tetra(pi0, pi1, pi4, pi7, b1, bn, b2, bn);
        make_tetra(pi1, pi2, pi3, pi7, b0, bn, bn, b4);
        make_tetra(pi1, pi2, pi6, pi7, b5, bn, bn, b4);
        make_tetra(pi1, pi4, pi5, pi7, b1, bn, bn, b3);
        make_tetra(pi1, pi5, pi6, pi7, b5, bn, bn, b3);


        auto comp_sort = [](const ns_pair& nsp1, const ns_pair& nsp2) -> bool {
            return nsp1.first < nsp2.first;
        };
        auto comp_uniq = [](const ns_pair& nsp1, const ns_pair& nsp2) -> bool {
            return nsp1.first == nsp2.first;
        };

        std::sort(new_surfaces.begin(), new_surfaces.end(), comp_sort);
        auto last_surf = std::unique(new_surfaces.begin(), new_surfaces.end(), comp_uniq);
        new_surfaces.erase(last_surf, new_surfaces.end());

        std::sort(storage->edges.begin(), storage->edges.end());
        auto last_edge = std::unique(storage->edges.begin(), storage->edges.end());
        storage->edges.erase(last_edge, storage->edges.end());

        storage->surfaces.reserve(new_surfaces.size());
        storage->boundary_info.reserve(new_surfaces.size());
        for (auto& ns : new_surfaces)
        {
            storage->surfaces.push_back(ns.first);
            storage->boundary_info.push_back(ns.second);
        }

        new_surfaces.clear();
        std::swap(storage->volumes, new_volumes);
        new_volumes.clear();
        assert(storage->surfaces.size() == 18);
        assert(storage->volumes.size() == 6);
        assert(storage->nodes.size() == 8);
        storage->subdomain_info.resize( storage->volumes.size() );
    }


    void refine(void)
    {
        size_t node_offset = storage->nodes.size();
        size_t edge_offset = storage->edges.size();
        size_t surface_offset = storage->surfaces.size();
        size_t volume_offset = storage->volumes.size();

        auto count_bnd_fcs = [](const boundary_descriptor& bi) -> bool {
            return bi.is_boundary();
        };

        size_t boundary_offset = std::count_if(storage->boundary_info.begin(),
            storage->boundary_info.end(),
            count_bnd_fcs);


        new_surfaces.clear();
        new_surfaces.reserve(storage->surfaces.size() * 4 + storage->volumes.size() * 4 );

        for (const auto& e : storage->edges)
        {
            auto ptids = e.point_ids();

            assert(ptids.size() == 2);
            assert(ptids[0] < storage->points.size());
            assert(ptids[1] < storage->points.size());

            auto p0 = storage->points[ ptids[0] ];
            auto p1 = storage->points[ ptids[1] ];
            auto pm = (p0 + p1)/2.;

            storage->points.push_back(pm);
            auto pmi = point_identifier<3>(storage->nodes.size());
            storage->nodes.push_back(node_type({ pmi }));

            assert(ptids[0] < pmi);
            assert( ptids[1] < pmi );
        }

        for (const auto& s : storage->volumes)
        {
            auto ptids = s.point_ids();

            auto eofs = [&](const edge_type& e) -> auto {
                auto be = begin(storage->edges);
                auto ee = begin(storage->edges) + edge_offset;
                auto ei = std::lower_bound(be, ee, e);

                if (ei == ee or e != *ei)
                    throw std::logic_error("Edge not found. This is a bug.");
                return point_identifier<3>(std::distance(be, ei) + node_offset);
                };

            auto ssearch = [&](const surface_type& s) -> auto {
                auto bs = begin(storage->surfaces);
                auto es = end(storage->surfaces);
                auto si = std::lower_bound(bs, es, s);
                if (si == es or s != *si)
                    throw std::logic_error("Surface not found. This is a bug.");
                return std::distance(bs, si);
            };

            auto s0 = surface_type({ ptids[0], ptids[1], ptids[2]});
            auto s1 = surface_type({ ptids[0], ptids[1], ptids[3]});
            auto s2 = surface_type({ ptids[0], ptids[2], ptids[3]});
            auto s3 = surface_type({ ptids[1], ptids[2], ptids[3]});

            auto bn = boundary_descriptor(0, false);
            auto b0 = storage->boundary_info[ssearch(s0)];
            auto b1 = storage->boundary_info[ssearch(s1)];
            auto b2 = storage->boundary_info[ssearch(s2)];
            auto b3 = storage->boundary_info[ssearch(s3)];

            auto pi4 = eofs(edge_type({ptids[0], ptids[1]}) );
            auto pi5 = eofs(edge_type({ptids[1], ptids[2]}));
            auto pi6 = eofs(edge_type({ptids[0], ptids[2]}));
            auto pi7 = eofs(edge_type({ptids[0], ptids[3]}));
            auto pi8 = eofs(edge_type({ptids[2], ptids[3]}));
            auto pi9 = eofs(edge_type({ptids[1], ptids[3]}));


            /* This is based on "Simplicial grid refinement: on Freudenthal's
             * algorithm and the optimal number of congruence classes" by
             * Jürgen Bey. Numerische Mathematik volume 85, pages 1–29 (2000).
             * See also https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/src/mesh/meshRefine.cpp#L192
             */

            make_tetra(ptids[0], pi4, pi6, pi7, b0, b1, b2, bn);
            make_tetra(pi4, ptids[1], pi5, pi9, b0, b1, bn, b3);
            make_tetra(pi6, pi5, ptids[2], pi8, b0, bn, b2, b3);
            make_tetra(pi7, pi9, pi8, ptids[3], bn, b1, b2, b3);
            make_tetra(pi4, pi6, pi7, pi9, bn, bn, b1, bn);
            make_tetra(pi4, pi9, pi5, pi6, bn, bn, b0, bn);
            make_tetra(pi6, pi7, pi9, pi8, bn, b2, bn, bn);
            make_tetra(pi6, pi8, pi9, pi5, bn, bn, bn, b3);

        }

        auto comp_sort = [](const ns_pair& nsp1, const ns_pair& nsp2) -> bool {
            return nsp1.first < nsp2.first;
        };
        auto comp_uniq = [](const ns_pair& nsp1, const ns_pair& nsp2) -> bool {
            return nsp1.first == nsp2.first;
        };

        auto n_itor = storage->edges.begin() + node_offset;
        storage->edges.erase(storage->edges.begin(), n_itor);
        std::sort(storage->edges.begin(), storage->edges.end());
        auto elast = std::unique(storage->edges.begin(), storage->edges.end());
        storage->edges.erase(elast, storage->edges.end());

        //auto v_itor = storage->volumes.begin() + volume_offset;
        //storage->volumes.erase(storage->volumes.begin(), v_itor);
        std::sort(new_volumes.begin(), new_volumes.end());
        std::swap(storage->volumes, new_volumes);
        new_volumes.clear();

        std::sort(new_surfaces.begin(), new_surfaces.end(), comp_sort);
        auto slast = std::unique(new_surfaces.begin(), new_surfaces.end(),comp_uniq);
        new_surfaces.erase(slast, new_surfaces.end());

        storage->surfaces.clear();
        storage->surfaces.reserve(new_surfaces.size());
        storage->boundary_info.clear();
        storage->boundary_info.reserve(new_surfaces.size());
        for (auto& ns : new_surfaces)
        {
            storage->surfaces.push_back(ns.first);
            storage->boundary_info.push_back(ns.second);
        }
        new_surfaces.clear();

        // expected surfaces = (init volumes * 8 tetras * 4 faces - init bnd fcs * 4 faces) /2 + init bnd fcs * 4 faces.
        size_t expected_surfaces = 16 * volume_offset + 2 * boundary_offset;
        assert(storage->volumes.size() == volume_offset * 8);
        assert(storage->surfaces.size() == expected_surfaces);
        storage->subdomain_info.resize( storage->volumes.size() );
    }
};

template<typename T>
class simple_mesher<cartesian_mesh<T,2>>
{
    typedef cartesian_mesh<T,2>                         mesh_type;
    static const size_t DIM = 2;
    typedef typename cartesian_mesh<T,2>::storage_type    storage_type;
    typedef point<T,DIM>                                point_type;

    typedef typename mesh_type::node_type       node_type;
    typedef typename mesh_type::edge_type       edge_type;
    typedef typename mesh_type::surface_type    surface_type;

    std::shared_ptr<storage_type>   storage;

public:
    simple_mesher(mesh_type& msh)
        : storage( msh.backend_storage() )
    {
        auto rot = 0.0;
        /* Init the first level of the mesh */
        storage->points.push_back( point_type(0.0, 0.0) );
        auto pi0 = point_identifier<2>(0);
        storage->nodes.push_back( node_type( {pi0} ) );

        storage->points.push_back( point_type(0.5-rot, 0.0) );
        auto pi1 = point_identifier<2>(1);
        storage->nodes.push_back( node_type( {pi1} ) );

        storage->points.push_back( point_type(1.0, 0.0) );
        auto pi2 = point_identifier<2>(2);
        storage->nodes.push_back( node_type( {pi2} ) );

        storage->points.push_back( point_type(0.0, 0.5+rot) );
        auto pi3 = point_identifier<2>(3);
        storage->nodes.push_back( node_type( {pi3} ) );

        storage->points.push_back( point_type(0.5, 0.5) );
        auto pi4 = point_identifier<2>(4);
        storage->nodes.push_back( node_type( {pi4} ) );

        storage->points.push_back( point_type(1.0, 0.5-rot) );
        auto pi5 = point_identifier<2>(5);
        storage->nodes.push_back( node_type( {pi5} ) );

        storage->points.push_back( point_type(0.0, 1.0) );
        auto pi6 = point_identifier<2>(6);
        storage->nodes.push_back( node_type( {pi6} ) );

        storage->points.push_back( point_type(0.5+rot, 1.0) );
        auto pi7 = point_identifier<2>(7);
        storage->nodes.push_back( node_type( {pi7} ) );

        storage->points.push_back( point_type(1.0, 1.0) );
        auto pi8 = point_identifier<2>(8);
        storage->nodes.push_back( node_type( {pi8} ) );

        storage->edges.push_back( edge_type({pi0, pi1}) ); //0 b0
        storage->edges.push_back( edge_type({pi0, pi3}) ); //1 b3
        storage->edges.push_back( edge_type({pi1, pi2}) ); //2 b0
        storage->edges.push_back( edge_type({pi1, pi4}) ); //3
        storage->edges.push_back( edge_type({pi2, pi5}) ); //4 b1
        storage->edges.push_back( edge_type({pi3, pi4}) ); //5
        storage->edges.push_back( edge_type({pi3, pi6}) ); //6 b3
        storage->edges.push_back( edge_type({pi4, pi5}) ); //7
        storage->edges.push_back( edge_type({pi4, pi7}) ); //8
        storage->edges.push_back( edge_type({pi5, pi8}) ); //9 b1
        storage->edges.push_back( edge_type({pi6, pi7}) ); //10 b2
        storage->edges.push_back( edge_type({pi7, pi8}) ); //11 b2

        storage->boundary_info.resize( storage->edges.size() );
        storage->boundary_info[0] = boundary_descriptor(0, true);
        storage->boundary_info[1] = boundary_descriptor(3, true);
        storage->boundary_info[2] = boundary_descriptor(0, true);
        storage->boundary_info[4] = boundary_descriptor(1, true);
        storage->boundary_info[6] = boundary_descriptor(3, true);
        storage->boundary_info[9] = boundary_descriptor(1, true);
        storage->boundary_info[10] = boundary_descriptor(2, true);
        storage->boundary_info[11] = boundary_descriptor(2, true);

        storage->surfaces.push_back( surface_type({pi0, pi1, pi3, pi4}) );
        storage->surfaces.push_back( surface_type({pi1, pi2, pi4, pi5}) );
        storage->surfaces.push_back( surface_type({pi3, pi4, pi6, pi7}) );
        storage->surfaces.push_back( surface_type({pi4, pi5, pi7, pi8}) );

        storage->subdomain_info.resize( storage->surfaces.size() );
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
            assert(ptids.size() == 4);
            assert(ptids[0] < storage->points.size());
            assert(ptids[1] < storage->points.size());
            assert(ptids[2] < storage->points.size());
            assert(ptids[3] < storage->points.size());

            auto p0 = storage->points[ ptids[0] ];
            auto p1 = storage->points[ ptids[1] ];
            auto p2 = storage->points[ ptids[2] ];
            auto p3 = storage->points[ ptids[3] ];
            auto pm = (p0 + p1 + p2 + p3) / 4.0;

            storage->points.push_back(pm);
            point_identifier<2> pmi(storage->nodes.size());
            storage->nodes.push_back( node_type({pmi}) );

            assert(ptids[0] < pmi);
            assert(ptids[1] < pmi);
            assert(ptids[2] < pmi);
            assert(ptids[3] < pmi);

            auto eofs = [&](const edge_type& e) -> auto {
                auto be = begin(storage->edges);
                auto ee = end(storage->edges);
                auto ei = std::lower_bound(be, ee, e);
                if (ei == ee or e != *ei) {
                    std::cout << s << std::endl;
                    std::cout << e << std::endl;
                    for (auto& ee : storage->edges)
                        std::cout << "  " << ee << std::endl;
                    throw std::logic_error("Edge not found. This is a bug.");
                }
                return point_identifier<2>(std::distance(be, ei) + node_offset);
            };

            auto p_e0 = eofs( edge_type({ptids[0], ptids[1]}) );
            auto p_e1 = eofs( edge_type({ptids[0], ptids[2]}) );
            auto p_e2 = eofs( edge_type({ptids[1], ptids[3]}) );
            auto p_e3 = eofs( edge_type({ptids[2], ptids[3]}) );

            /* Those are going to be internal edges */
            
            new_edges.push_back( std::make_pair(edge_type({p_e0, pmi}), boundary_descriptor()) );
            new_edges.push_back( std::make_pair(edge_type({p_e1, pmi}), boundary_descriptor()) );
            new_edges.push_back( std::make_pair(edge_type({p_e2, pmi}), boundary_descriptor()) );
            new_edges.push_back( std::make_pair(edge_type({p_e3, pmi}), boundary_descriptor()) );

            new_surfaces.push_back( surface_type( {ptids[0], p_e0, p_e1, pmi} ) );
            new_surfaces.push_back( surface_type( {ptids[1], p_e2, p_e0, pmi} ) );
            new_surfaces.push_back( surface_type( {ptids[2], p_e1, p_e3, pmi} ) );
            new_surfaces.push_back( surface_type( {ptids[3], p_e3, p_e2, pmi} ) );
        }

        auto comp = [](const ne_pair& nep1, const ne_pair& nep2) -> bool {
            return nep1.first < nep2.first;
        };
        std::sort(new_edges.begin(), new_edges.end(), comp);

        storage->edges.clear(); storage->edges.reserve( new_edges.size() );
        storage->boundary_info.clear();
        storage->boundary_info.reserve(new_edges.size());
        for (auto& ne : new_edges)
        {
            storage->edges.push_back(ne.first);
            storage->boundary_info.push_back(ne.second);
        }

        std::sort(new_surfaces.begin(), new_surfaces.end());
        std::swap(storage->surfaces, new_surfaces);
        storage->subdomain_info.resize( storage->surfaces.size() );
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

                boundary_descriptor bi(size_t(boundary), true);
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
void make_single_element_mesh(simplicial_mesh<T,3>& msh,
    const point<T,3>& a,
    const point<T,3>& b,
    const point<T,3>& c,
    const point<T,3>& d)
{
    using mesh_type = simplicial_mesh<T,3>;
    using point_type = typename mesh_type::point_type;
    using node_type = typename mesh_type::node_type;
    using edge_type = typename mesh_type::edge_type;
    using nodeid_type = typename node_type::id_type;

    auto storage = msh.backend_storage();

    storage->points.push_back( a );
    storage->points.push_back( b );
    storage->points.push_back( c );
    storage->points.push_back( d );

    auto p0 = point_identifier<3>(0);
    auto p1 = point_identifier<3>(1);
    auto p2 = point_identifier<3>(2);
    auto p3 = point_identifier<3>(3);

    // Base face
    storage->edges.push_back( {p0, p1} );
    storage->edges.push_back( {p0, p2} );
    storage->edges.push_back( {p1, p2} );
    storage->surfaces.push_back( {p0, p1, p2} );

    // Front face
    storage->edges.push_back( {p0, p3} );
    storage->edges.push_back( {p1, p3} );
    storage->surfaces.push_back( {p0, p1, p3} );

    // Lateral face
    storage->edges.push_back( {p2, p3} );
    storage->surfaces.push_back( {p1, p2, p3} );

    // Back face
    storage->surfaces.push_back( {p0, p2, p3} );

    storage->volumes.push_back( {p0, p1, p2, p3} );

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
    storage->points.push_back( point_type(    base.x(), base.y()+hy ) );
    storage->points.push_back( point_type( base.x()+hx, base.y()+hy ) );


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


template<typename T>
void make_single_element_mesh(cartesian_mesh<T,3>& msh, const point<T,3>& base, const T& hx, const T& hy, const T& hz)
{
    using mesh_type = cartesian_mesh<T,3>;
    using point_type = typename mesh_type::point_type;
    using node_type = typename mesh_type::node_type;
    using edge_type = typename mesh_type::edge_type;
    using nodeid_type = typename node_type::id_type;

    auto storage = msh.backend_storage();

    storage->points.push_back(point_type( base.x(),    base.y(), base.z()));
    storage->points.push_back(point_type( base.x()+hx, base.y(), base.z()));
    storage->points.push_back(point_type( base.x()   , base.y()+hy, base.z()));
    storage->points.push_back(point_type( base.x()+hx, base.y()+hy, base.z()));
    storage->points.push_back(point_type( base.x()   , base.y(), base.z()+hz));
    storage->points.push_back(point_type( base.x()+hx, base.y(), base.z()+hz));
    storage->points.push_back(point_type( base.x()   , base.y()+hy, base.z()+hz));
    storage->points.push_back(point_type( base.x()+hx, base.y()+hy, base.z()+hz));

    auto p0 = point_identifier<3>(0);
    auto p1 = point_identifier<3>(1);
    auto p2 = point_identifier<3>(2);
    auto p3 = point_identifier<3>(3);
    auto p4 = point_identifier<3>(4);
    auto p5 = point_identifier<3>(5);
    auto p6 = point_identifier<3>(6);
    auto p7 = point_identifier<3>(7);

    storage->nodes.push_back( { p0 } );
    storage->nodes.push_back( { p1 } );
    storage->nodes.push_back( { p2 } );
    storage->nodes.push_back( { p3 } );
    storage->nodes.push_back( { p4 } );
    storage->nodes.push_back( { p5 } );
    storage->nodes.push_back( { p6 } );
    storage->nodes.push_back( { p7 } );

    storage->edges.push_back( {p0, p1} );
    storage->edges.push_back( {p0, p2} );
    storage->edges.push_back( {p1, p3} );
    storage->edges.push_back( {p2, p3} );
    storage->edges.push_back( {p0, p4} );
    storage->edges.push_back( {p1, p5} );
    storage->edges.push_back( {p3, p7} );
    storage->edges.push_back( {p2, p6} );
    storage->edges.push_back( {p4, p6} );
    storage->edges.push_back( {p6, p7} );
    storage->edges.push_back( {p4, p5} );
    storage->edges.push_back( {p5, p7} );
    std::sort(storage->edges.begin(), storage->edges.end());

    storage->surfaces.push_back({ p0, p1, p3, p2 });
    storage->surfaces.push_back( {p4, p5, p7, p6} );
    storage->surfaces.push_back( {p0, p2, p6, p4} );
    storage->surfaces.push_back( {p1, p3, p7, p5} );
    storage->surfaces.push_back( {p0, p1, p5, p4} );
    storage->surfaces.push_back({ p2, p3, p7, p6 });
    std::sort(storage->surfaces.begin(), storage->surfaces.end());

    storage->volumes.push_back( {p0, p1, p2, p3, p4, p5, p6, p7} );

    storage->subdomain_info.resize( 1 );
}


template<typename T>
void make_single_element_mesh(generic_mesh<T,2>& msh, const T& radius, const size_t& faces)
{
    using mesh_type = generic_mesh<T,2>;
    using point_type = typename mesh_type::point_type;
    using node_type = typename mesh_type::node_type;
    using edge_type = typename mesh_type::edge_type;
    using surface_type = typename mesh_type::surface_type;
    using nodeid_type = typename node_type::id_type;

    if (faces < 3) {
        std::cout << "Can't make a polygon with less than three faces" << std::endl;
        return;
    }

    auto storage = msh.backend_storage();

    //std::random_device rd;  // Will be used to obtain a seed for the random number engine
    //std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    //std::uniform_real_distribution<> dis(-0.1, 1.0);

    for (size_t k = 0; k < faces; k++) {
        auto x = radius*std::cos(2*k*M_PI/faces);// + dis(gen);
        auto y = radius*std::sin(2*k*M_PI/faces);// + dis(gen);
        storage->points.push_back( point_type{x,y} );

        auto p = point_identifier<2>(k);
        storage->nodes.push_back( node_type(p) );
    }

    std::vector<typename node_type::id_type> surface_nodes;
    std::vector<typename edge_type::id_type> surface_edges;
    auto num_nodes = storage->nodes.size();
    for (size_t i = 0; i < num_nodes; i++) {
        auto n0 = nodeid_type(i);
        auto n1 = nodeid_type( (i+1)%num_nodes );
        storage->edges.push_back( {n0, n1} );
        surface_nodes.push_back( typename node_type::id_type(i) );
        surface_edges.push_back( typename edge_type::id_type(i) );
    }

    storage->boundary_info.resize( storage->edges.size() );
    for (auto& bi : storage->boundary_info)
        bi = boundary_descriptor(0, true);

    auto surface = surface_type(surface_edges);
    surface.set_point_ids(surface_nodes.begin(), surface_nodes.end()); /* XXX: crap */

    storage->surfaces.push_back(surface);

    storage->subdomain_info.resize( 1 );
}

template<typename T>
bool load_single_element_csv(generic_mesh<T,2>& msh, const std::string& filename)
{
    using mesh_type = generic_mesh<T,2>;
    using point_type = typename mesh_type::point_type;
    using node_type = typename mesh_type::node_type;
    using edge_type = typename mesh_type::edge_type;
    using surface_type = typename mesh_type::surface_type;
    using nodeid_type = typename node_type::id_type;

    auto storage = msh.backend_storage();


    std::ifstream ifs(filename);
    if (not ifs.is_open()) {
        std::cout << "Can't open " << filename << std::endl;
        return false;
    }

    std::string line;
    getline(ifs, line);

    size_t ptnum = 0;
    while( not ifs.eof() ) {
        getline(ifs, line);
        if (line == "")
            break;

        std::regex re(";");
        std::sregex_token_iterator iter(line.begin(), line.end(), re, -1);
        std::vector<std::string> tokens{ iter, {} };
        if (tokens.size() < 2) {
            std::cout << "Error parsing CSV input: " << tokens.size() << std::endl;
            return false;
        }
        auto x = std::stod(tokens[0]);
        auto y = std::stod(tokens[1]);

        storage->points.push_back( point_type{x,y} );
        auto p = point_identifier<2>(ptnum++);
        storage->nodes.push_back( node_type(p) );
    }

    std::vector<typename node_type::id_type> surface_nodes;
    std::vector<typename edge_type::id_type> surface_edges;
    auto num_nodes = storage->nodes.size();
    for (size_t i = 0; i < num_nodes; i++) {
        auto n0 = nodeid_type(i);
        auto n1 = nodeid_type( (i+1)%num_nodes );
        storage->edges.push_back( {n0, n1} );
        surface_nodes.push_back( typename node_type::id_type(i) );
        surface_edges.push_back( typename edge_type::id_type(i) );
    }

    storage->boundary_info.resize( storage->edges.size() );
    for (auto& bi : storage->boundary_info)
        bi = boundary_descriptor(0, true);

    auto surface = surface_type(surface_edges);
    surface.set_point_ids(surface_nodes.begin(), surface_nodes.end()); /* XXX: crap */

    storage->surfaces.push_back(surface);

    storage->subdomain_info.resize( 1 );

    return true;
}


} //namespace disk
