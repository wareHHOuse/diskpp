/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *  
 * Matteo Cicuttin (C) 2021
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group  
 */

#pragma once

#include "gmsh.h"

#include "core/mesh/mesh_storage.hpp"

#define INVALID_OFS ((size_t) (~0))

namespace disk {

template<typename Mesh>
class gmsh_geometry_loader;

template<typename T>
class gmsh_geometry_loader<simplicial_mesh<T,3>> : public mesh_loader<simplicial_mesh<T,3>>
{
    typedef simplicial_mesh<T,3>                    mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;
    typedef typename mesh_type::surface_type        surface_type;
    typedef typename mesh_type::volume_type         volume_type;

    std::vector<point_type>                         points;
    std::vector<node_type>                          nodes;
    std::vector<edge_type>                          edges;
    std::vector<surface_type>                       surfaces;
    std::vector<boundary_descriptor>                boundary_info;

    std::vector<volume_type>                                    volumes;
    std::vector<std::pair<volume_type, subdomain_descriptor>>   tmp_volumes;

    std::vector<size_t>     node_tag2ofs;

    void gmsh_get_nodes()
    {
        std::vector<size_t>     nodeTags;
        std::vector<double>     coords;
        std::vector<double>     paraCoords;

        gmsh::model::mesh::getNodes(nodeTags, coords, paraCoords, -1, -1, true, false);

        auto maxtag_pos = std::max_element(nodeTags.begin(), nodeTags.end());
        
        size_t maxtag = 0;
        if (maxtag_pos != nodeTags.end())
            maxtag = (*maxtag_pos)+1;

        points.resize( nodeTags.size() );

        for (size_t i = 0; i < nodeTags.size(); i++)
        {
            auto x = coords[3*i+0];
            auto y = coords[3*i+1];
            auto z = coords[3*i+2];
            points[i] = point_type(x,y,z);
        }

        node_tag2ofs.resize(maxtag, INVALID_OFS);
        for (size_t i = 0; i < nodeTags.size(); i++)
        {
            node_tag2ofs.at( nodeTags[i] ) = i;
            auto point_id = disk::point_identifier<3>( i );
            auto node = node_type( { point_id } );
            nodes.push_back(node);
        }
    }

    void gmsh_get_elements()
    {
        gmsh::vectorpair entities;
        gmsh::model::getEntities(entities, 3/*dimension*/);
        size_t subdom_id = 0;
        for (auto [dim, tag] : entities)
        {
            std::vector<int> elemTypes;
            gmsh::model::mesh::getElementTypes(elemTypes, dim, tag);
            for (auto& elemType : elemTypes)
            {
                std::vector<size_t> elemTags;
                std::vector<size_t> elemNodeTags;
                gmsh::model::mesh::getElementsByType(elemType, elemTags, elemNodeTags, tag);
                auto nodesPerElem = elemNodeTags.size()/elemTags.size();
                assert( elemTags.size() * nodesPerElem == elemNodeTags.size() );
            
                for (size_t i = 0; i < elemTags.size(); i++)
                {
                    auto base = nodesPerElem * i;

                    auto node0_tag = elemNodeTags[base + 0];
                    assert(node0_tag < node_tag2ofs.size());
                    auto node0_ofs = node_tag2ofs[node0_tag];
                    assert(node0_ofs != INVALID_OFS);

                    auto node1_tag = elemNodeTags[base + 1];
                    assert(node1_tag < node_tag2ofs.size());
                    auto node1_ofs = node_tag2ofs[node1_tag];
                    assert(node1_ofs != INVALID_OFS);

                    auto node2_tag = elemNodeTags[base + 2];
                    assert(node2_tag < node_tag2ofs.size());
                    auto node2_ofs = node_tag2ofs[node2_tag];
                    assert(node2_ofs != INVALID_OFS);

                    auto node3_tag = elemNodeTags[base + 3];
                    assert(node3_tag < node_tag2ofs.size());
                    auto node3_ofs = node_tag2ofs[node3_tag];
                    assert(node3_ofs != INVALID_OFS);

                    point_identifier<3>     p0(node0_ofs);
                    point_identifier<3>     p1(node1_ofs);
                    point_identifier<3>     p2(node2_ofs);
                    point_identifier<3>     p3(node3_ofs);

                    edges.push_back( edge_type( { p0, p1 } ) );
                    edges.push_back( edge_type( { p0, p2 } ) );
                    edges.push_back( edge_type( { p0, p3 } ) );
                    edges.push_back( edge_type( { p1, p2 } ) );
                    edges.push_back( edge_type( { p1, p3 } ) );
                    edges.push_back( edge_type( { p2, p3 } ) );

                    surfaces.push_back( surface_type( { p0, p1, p2 } ) );
                    surfaces.push_back( surface_type( { p0, p1, p3 } ) );
                    surfaces.push_back( surface_type( { p0, p2, p3 } ) );
                    surfaces.push_back( surface_type( { p1, p2, p3 } ) );

                    subdomain_descriptor si(i, tag);
                    auto vs = std::make_pair(volume_type( {p0, p1, p2, p3} ), si);
                    tmp_volumes.push_back( vs );
                }
            }

            subdom_id++;
        }

        using pvs = std::pair<volume_type, subdomain_descriptor>;
        auto tmpvol_comp = [](const pvs& a, const pvs& b) {
            return a.first < b.first;
        };

        priv::sort_uniq(edges);
        priv::sort_uniq(surfaces);
        boundary_info.resize(surfaces.size());
        std::sort(tmp_volumes.begin(), tmp_volumes.end(), tmpvol_comp);
    }

    void detect_boundary_faces(void)
    {
        gmsh::vectorpair entities;

        size_t b_id = 0;
        gmsh::model::getEntities(entities, 2/*dimension*/);
        for (auto [dim, tag] : entities)
        {
            std::vector<int> elemTypes;
            gmsh::model::mesh::getElementTypes(elemTypes, dim, tag);
            assert(elemTypes.size() == 1);
            for (auto& elemType : elemTypes)
            {
                std::vector<size_t> nTags;
                gmsh::model::mesh::getElementFaceNodes(elemType, 3, nTags, tag, true);
                
                std::vector<surface_type> fkeys;
                for (size_t i = 0; i < nTags.size(); i+=3)
                {
                    auto node0_tag = nTags[i + 0];
                    assert(node0_tag < node_tag2ofs.size());
                    auto node0_ofs = node_tag2ofs[node0_tag];
                    assert(node0_ofs != INVALID_OFS);

                    auto node1_tag = nTags[i + 1];
                    assert(node1_tag < node_tag2ofs.size());
                    auto node1_ofs = node_tag2ofs[node1_tag];
                    assert(node1_ofs != INVALID_OFS);

                    auto node2_tag = nTags[i + 2];
                    assert(node2_tag < node_tag2ofs.size());
                    auto node2_ofs = node_tag2ofs[node2_tag];
                    assert(node2_ofs != INVALID_OFS);

                    point_identifier<3>     p0(node0_ofs);
                    point_identifier<3>     p1(node1_ofs);
                    point_identifier<3>     p2(node2_ofs);

                    boundary_descriptor bi(b_id, tag, true);
                    surface_type surf( { p0, p1, p2 } );

                    auto itor = std::lower_bound(surfaces.begin(), surfaces.end(), surf);
                    if ( (itor == surfaces.end()) or not (*itor == surf) )
                        throw std::logic_error("Face not found");

                    auto ofs = std::distance(surfaces.begin(), itor);
                    boundary_info.at(ofs) = bi;
                }
            }

            b_id++;
        }
    }

public:
    static const char constexpr *expected_extension = "geo";
    gmsh_geometry_loader() = default;

    bool read_mesh(const std::string& s)
    {
        gmsh::initialize(0, nullptr);
        
        if (this->verbose())
            gmsh::option::setNumber("General.Terminal", 1);
        
        gmsh::open( s ); //HANDLE ERRORS!
        gmsh::model::mesh::generate( 3 /*dimension*/ );
        gmsh::model::mesh::setOrder( 1 /*element order*/ );
        gmsh_get_nodes();
        gmsh_get_elements();
        detect_boundary_faces();

        gmsh::finalize();
        return true;
    }

    bool populate_mesh(mesh_type& msh)
    {
        auto storage = msh.backend_storage();
        storage->points = std::move(points);
        storage->nodes = std::move(nodes);
        storage->edges = std::move(edges);
        storage->surfaces = std::move(surfaces);
        storage->boundary_info = std::move(boundary_info);

        std::vector<volume_type> vols;
        vols.reserve( tmp_volumes.size() );

        std::vector<subdomain_descriptor> subdoms;
        subdoms.reserve( tmp_volumes.size() );

        for (auto& [vol, sd] : tmp_volumes)
        {
            vols.push_back(vol);
            subdoms.push_back(sd);
        }

        tmp_volumes.clear();

        storage->volumes = std::move(vols);
        storage->subdomain_info = std::move(subdoms);

        mark_internal_faces(msh);

        return true;
    }
};


} //namespace disk