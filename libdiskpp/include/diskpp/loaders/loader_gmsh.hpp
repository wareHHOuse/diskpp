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

#include "diskpp/mesh/mesh_storage.hpp"

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
        else
            gmsh::option::setNumber("General.Terminal", 0);
        
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

template<typename T>
class gmsh_geometry_loader<simplicial_mesh<T,2>> : public mesh_loader<simplicial_mesh<T,2>>
{
    typedef simplicial_mesh<T,2>                    mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;
    typedef typename mesh_type::surface_type        surface_type;

    std::vector<point_type>                         points;
    std::vector<node_type>                          nodes;
    std::vector<edge_type>                          edges;
    std::vector<boundary_descriptor>                boundary_info;

    std::vector<surface_type>                       surfaces;
    std::vector<std::pair<surface_type, subdomain_descriptor>>   tmp_surfaces;

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
            points[i] = point_type(x,y);
        }

        node_tag2ofs.resize(maxtag, INVALID_OFS);
        for (size_t i = 0; i < nodeTags.size(); i++)
        {
            node_tag2ofs.at( nodeTags[i] ) = i;
            auto point_id = disk::point_identifier<2>( i );
            auto node = node_type( { point_id } );
            nodes.push_back(node);
        }
    }

    void gmsh_get_elements()
    {
        gmsh::vectorpair entities;
        gmsh::model::getEntities(entities, 2/*dimension*/);
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

                    point_identifier<2>     p0(node0_ofs);
                    point_identifier<2>     p1(node1_ofs);
                    point_identifier<2>     p2(node2_ofs);

                    edges.push_back( edge_type( { p0, p1 } ) );
                    edges.push_back( edge_type( { p0, p2 } ) );
                    edges.push_back( edge_type( { p1, p2 } ) );

                    subdomain_descriptor si(i, tag);
                    auto vs = std::make_pair(surface_type( {p0, p1, p2} ), si);
                    tmp_surfaces.push_back( vs );
                }
            }

            subdom_id++;
        }

        using pvs = std::pair<surface_type, subdomain_descriptor>;
        auto tmpsurf_comp = [](const pvs& a, const pvs& b) {
            return a.first < b.first;
        };

        priv::sort_uniq(edges);
        boundary_info.resize(edges.size());
        std::sort(tmp_surfaces.begin(), tmp_surfaces.end(), tmpsurf_comp);
    }

    void detect_boundary_edges(void)
    {
        gmsh::vectorpair entities;

        size_t b_id = 0;
        gmsh::model::getEntities(entities, 1/*dimension*/);
        for (auto [dim, tag] : entities)
        {
            std::vector<int> elemTypes;
            gmsh::model::mesh::getElementTypes(elemTypes, dim, tag);
            assert(elemTypes.size() == 1);
            for (auto& elemType : elemTypes)
            {
                std::vector<size_t> nTags;
                gmsh::model::mesh::getElementEdgeNodes(elemType, nTags, tag, true);
                
                std::vector<edge_type> fkeys;
                for (size_t i = 0; i < nTags.size(); i+=2)
                {
                    auto node0_tag = nTags[i + 0];
                    assert(node0_tag < node_tag2ofs.size());
                    auto node0_ofs = node_tag2ofs[node0_tag];
                    assert(node0_ofs != INVALID_OFS);

                    auto node1_tag = nTags[i + 1];
                    assert(node1_tag < node_tag2ofs.size());
                    auto node1_ofs = node_tag2ofs[node1_tag];
                    assert(node1_ofs != INVALID_OFS);

                    point_identifier<2>     p0(node0_ofs);
                    point_identifier<2>     p1(node1_ofs);

                    boundary_descriptor bi(b_id, tag, true);
                    edge_type edg( { p0, p1 } );

                    auto itor = std::lower_bound(edges.begin(), edges.end(), edg);
                    if ( (itor == edges.end()) or not (*itor == edg) )
                        throw std::logic_error("Edge not found");

                    auto ofs = std::distance(edges.begin(), itor);
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
        else
            gmsh::option::setNumber("General.Terminal", 0);
        
        gmsh::open( s ); //HANDLE ERRORS!
        gmsh::model::mesh::generate( 2 /*dimension*/ );
        gmsh::model::mesh::setOrder( 1 /*element order*/ );
        gmsh_get_nodes();
        gmsh_get_elements();
        detect_boundary_edges();

        gmsh::finalize();
        return true;
    }

    bool read_mesh()
    {
        //if ( ! gmsh::isInitialized() ) {
        //    std::cout << "GMSH is not initialized, can't read mesh" << std::endl;
        //    return false;
        //}
        
        if (this->verbose())
            gmsh::option::setNumber("General.Terminal", 1);
        else
            gmsh::option::setNumber("General.Terminal", 0);
        
        gmsh::model::mesh::generate( 2 /*dimension*/ );
        gmsh::model::mesh::setOrder( 1 /*element order*/ );
        gmsh_get_nodes();
        gmsh_get_elements();
        detect_boundary_edges();

        return true;
    }

    bool populate_mesh(mesh_type& msh)
    {
        auto storage = msh.backend_storage();
        storage->points = std::move(points);
        storage->nodes = std::move(nodes);
        storage->edges = std::move(edges);
        storage->boundary_info = std::move(boundary_info);

        std::vector<surface_type> surfs;
        surfs.reserve( tmp_surfaces.size() );

        std::vector<subdomain_descriptor> subdoms;
        subdoms.reserve( tmp_surfaces.size() );

        for (auto& [surf, sd] : tmp_surfaces)
        {
            surfs.push_back(surf);
            subdoms.push_back(sd);
        }

        tmp_surfaces.clear();

        storage->surfaces = std::move(surfs);
        storage->subdomain_info = std::move(subdoms);

        mark_internal_faces(msh);

        return true;
    }
};

template<disk::mesh_2D Mesh>
void submesh_via_gmsh(const Mesh& msh, disk::simplicial_mesh<typename Mesh::coordinate_type, 2>& dmsh, double scale)
{
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 0);
    gmsh::option::setNumber("Mesh.Algorithm", 1);
    gmsh::option::setNumber("Mesh.Algorithm3D", 1);

    gmsh::model::add("test");

    static const int point_ofs = 1000;

    for (int i = 0; i < msh.points_size(); i++) {
        auto pt = msh.point_at(i);
        gmsh::model::occ::addPoint(pt.x(), pt.y(), 0.0, 0.0, i+point_ofs);
    }

    std::vector<int> all_linetags;
    all_linetags.reserve(msh.faces_size());
    for (const auto& fc : faces(msh)) {
        int tag = (int) offset(msh, fc) +1;
        auto ptids = fc.point_ids();
        assert(ptids.size() == 2);
        all_linetags.push_back( gmsh::model::occ::addLine(ptids[0]+point_ofs, ptids[1]+point_ofs, tag) );
    }

    for (const auto& cl : msh) {
        std::vector<int> local_linetags;
        auto fcs = faces(msh, cl);
        for (const auto& fc : fcs)
            local_linetags.push_back( all_linetags[ offset(msh, fc) ] );

        int loop = gmsh::model::occ::addCurveLoop(local_linetags);
        gmsh::model::occ::addPlaneSurface({loop}, (int)offset(msh,cl)+1);
    
        auto h = diameter(msh, cl);
        auto ptids = cl.point_ids();
        gmsh::vectorpair vp;
        for (auto& ptid : ptids)
            vp.push_back({0, ptid+point_ofs});
        gmsh::model::occ::mesh::setSize(vp, h*scale);
    }

    gmsh::model::occ::synchronize();

    disk::gmsh_geometry_loader<disk::simplicial_mesh<double,2>> loader;
    loader.read_mesh();
    loader.populate_mesh(dmsh);

    gmsh::finalize();
}

namespace priv {

template<typename T>
struct temp_face
{
    std::vector<T>      nodes;

    temp_face() = default;
    temp_face(const temp_face&) = default;
    temp_face(temp_face&& other) = default;

    temp_face& operator=(const temp_face&) = default;

    temp_face(std::initializer_list<T> l)
    {
        nodes = std::vector<T>(l);
    }

    bool operator<(const temp_face& other) const
    {
        auto n_mine = nodes;
        auto n_other = other.nodes;
        std::sort(n_mine.begin(), n_mine.end());
        std::sort(n_other.begin(), n_other.end());
        return std::lexicographical_compare(n_mine.begin(), n_mine.end(),
                                            n_other.begin(), n_other.end());
    }

    bool operator==(const temp_face& other) const
    {
        if (nodes.size() != other.nodes.size())
            return false;

        auto n_mine = nodes;
        auto n_other = other.nodes;
        std::sort(n_mine.begin(), n_mine.end());
        std::sort(n_other.begin(), n_other.end());
        assert(n_mine.size() == n_other.size());
        for (size_t i = 0; i < n_mine.size(); i++)
            if ( n_mine[i] != n_other[i] )
                return false;

        return true;
    }
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const temp_face<T>& f)
{
    for (auto& n : f.nodes)
        os << n << " ";
    return os;
}

template<typename T>
struct temp_volume
{
    std::vector<T>              nodes;
    std::vector<temp_face<T>>   temp_faces;
    subdomain_descriptor        di;

    temp_volume(std::initializer_list<T> l, const subdomain_descriptor& p_di)
    {
        nodes = std::vector<T>(l);
        di = p_di;
    }
};

} // namespace priv

#define GMSH_SHAPE_PRISM 6

template<typename T>
class gmsh_geometry_loader<generic_mesh<T,3>> : public mesh_loader<generic_mesh<T,3>>
{
    typedef generic_mesh<T,3>                       mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;
    typedef typename mesh_type::surface_type        surface_type;
    typedef typename mesh_type::volume_type         volume_type;

    std::vector<point_type>                                     points;
    std::vector<node_type>                                      nodes;
    std::vector<edge_type>                                      edges;
    std::vector<surface_type>                                   surfaces;
    std::vector<std::pair<volume_type, subdomain_descriptor>>   tmp_mesh_vols;
    std::vector<boundary_descriptor>                            boundary_info;

    std::vector<std::pair<surface_type, size_t>>                    surf_map_table;
    std::vector<priv::temp_face<typename node_type::id_type>>       tmp_surfaces;
    std::vector<priv::temp_volume<typename node_type::id_type>>     tmp_volumes;

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

    void gmsh_get_prisms(size_t subdom_id, size_t tag, int elemType)
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

            auto node4_tag = elemNodeTags[base + 4];
            assert(node4_tag < node_tag2ofs.size());
            auto node4_ofs = node_tag2ofs[node4_tag];
            assert(node4_ofs != INVALID_OFS);

            auto node5_tag = elemNodeTags[base + 5];
            assert(node5_tag < node_tag2ofs.size());
            auto node5_ofs = node_tag2ofs[node5_tag];
            assert(node5_ofs != INVALID_OFS);

            typename node_type::id_type n0(node0_ofs);
            typename node_type::id_type n1(node1_ofs);
            typename node_type::id_type n2(node2_ofs);
            typename node_type::id_type n3(node3_ofs);
            typename node_type::id_type n4(node4_ofs);
            typename node_type::id_type n5(node5_ofs);

            edges.push_back( edge_type( { n0, n1 } ) );
            edges.push_back( edge_type( { n0, n2 } ) );
            edges.push_back( edge_type( { n0, n3 } ) );
            edges.push_back( edge_type( { n1, n2 } ) );
            edges.push_back( edge_type( { n1, n4 } ) );
            edges.push_back( edge_type( { n2, n5 } ) );
            edges.push_back( edge_type( { n3, n4 } ) );
            edges.push_back( edge_type( { n3, n5 } ) );
            edges.push_back( edge_type( { n4, n5 } ) );

            tmp_surfaces.push_back( priv::temp_face({n0, n1, n2}) );
            tmp_surfaces.push_back( priv::temp_face({n0, n1, n4, n3}) );
            tmp_surfaces.push_back( priv::temp_face({n0, n2, n5, n3}) );
            tmp_surfaces.push_back( priv::temp_face({n1, n2, n5, n4}) );
            tmp_surfaces.push_back( priv::temp_face({n3, n4, n5}) );

            subdomain_descriptor di(subdom_id, tag);
            auto tv = priv::temp_volume({n0, n1, n2, n3, n4, n5}, di);
            tv.temp_faces.push_back( priv::temp_face({n0, n1, n2}) );
            tv.temp_faces.push_back( priv::temp_face({n0, n1, n4, n3}) );
            tv.temp_faces.push_back( priv::temp_face({n0, n2, n5, n3}) );
            tv.temp_faces.push_back( priv::temp_face({n1, n2, n5, n4}) );
            tv.temp_faces.push_back( priv::temp_face({n3, n4, n5}) );
            tmp_volumes.push_back(tv);
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
                switch (elemType)
                {
                    case GMSH_SHAPE_PRISM:
                        if (this->verbose())
                            std::cout << "Entity " << tag << ": prisms" << std::endl;
                        gmsh_get_prisms(subdom_id, tag, elemType);
                        break;

                    default:
                        std::cout << "Unsupported element type " << elemType << std::endl;
                        break;
                }
            }

            subdom_id++;
        }

        make_surfaces();
        make_volumes();
    }

    void make_surfaces(void)
    {
        priv::sort_uniq(edges);
        priv::sort_uniq(tmp_surfaces);

        for (auto& ts : tmp_surfaces)
        {
            std::vector<disk::point_identifier<3>>      point_ids;
            std::vector<typename edge_type::id_type>    edge_ids;

            for (size_t i = 0; i < ts.nodes.size(); i++)
            {
                auto n0 = ts.nodes[i];
                auto n1 = ts.nodes[(i+1) % ts.nodes.size()];
                auto e = edge_type(n0, n1);
                auto itor = std::lower_bound(edges.begin(), edges.end(), e);
                if ( (itor == edges.end()) or not (*itor == e) )
                    throw("edge not found");
                auto d = std::distance(edges.begin(), itor);
                edge_ids.push_back( typename edge_type::id_type(d) );
                point_ids.push_back( disk::point_identifier<3>(n0) );
            }

            surface_type surf(edge_ids);
            surf.set_point_ids( point_ids.begin(), point_ids.end() );
            surfaces.push_back(surf);
        }

        tmp_surfaces.clear();
        std::sort(surfaces.begin(), surfaces.end());
        boundary_info.resize(surfaces.size());

        surf_map_table.resize( surfaces.size() );
        for (size_t i = 0; i < surfaces.size(); i++)
            surf_map_table[i] = std::make_pair(surfaces[i], i);

        auto comp = [](const std::pair<surface_type, size_t>& a, const std::pair<surface_type, size_t>& b) -> bool {
            auto pa = a.first.point_ids();
            auto pb = b.first.point_ids();
            std::sort(pa.begin(), pa.end());
            std::sort(pb.begin(), pb.end());
            return std::lexicographical_compare(pa.begin(), pa.end(), pb.begin(), pb.end());
        };

        std::sort(surf_map_table.begin(), surf_map_table.end(), comp);
    }

    void make_volumes(void)
    {
        using tf_t = priv::temp_face<typename node_type::id_type>;
        auto comp = [](const std::pair<surface_type, size_t>& a, const tf_t& tf)
        {
            auto pa = a.first.point_ids();
            auto pb = convert_to<disk::point_identifier<3>>(tf.nodes);
            std::sort(pa.begin(), pa.end());
            std::sort(pb.begin(), pb.end());
            return std::lexicographical_compare(pa.begin(), pa.end(), pb.begin(), pb.end());
        };

        auto eq = [](const std::pair<surface_type, size_t>& a, const tf_t& tf)
        {
            auto pa = a.first.point_ids();
            auto pb = convert_to<disk::point_identifier<3>>(tf.nodes);
            std::sort(pa.begin(), pa.end());
            std::sort(pb.begin(), pb.end());
            return std::equal(pa.begin(), pa.end(), pb.begin());
        };

        for (auto& tv : tmp_volumes)
        {
            std::vector<disk::point_identifier<3>>          point_ids;
            std::vector<typename surface_type::id_type>     surface_ids;

            for (auto& tf : tv.temp_faces)
            {
                auto itor = std::lower_bound(surf_map_table.begin(), surf_map_table.end(), tf, comp);
                if ( (itor == surf_map_table.end()) or not eq(*itor, tf) )
                    throw("face not found");
                auto d = (*itor).second;
                surface_ids.push_back( typename surface_type::id_type(d) );
            }

            point_ids = convert_to<disk::point_identifier<3>>(tv.nodes);

            volume_type vol(surface_ids);
            vol.set_point_ids( point_ids.begin(), point_ids.end() );
            tmp_mesh_vols.push_back( std::make_pair(vol, tv.di) );
        }

        tmp_volumes.clear();

        using pvs = std::pair<volume_type, subdomain_descriptor>;
        auto tmpvol_comp = [](const pvs& a, const pvs& b) {
            return a.first < b.first;
        };

        std::sort(tmp_mesh_vols.begin(), tmp_mesh_vols.end(), tmpvol_comp);
    }

    void detect_boundary_faces(void)
    {
        gmsh::vectorpair entities;

        auto comp = [](const std::pair<surface_type, size_t>& a, const std::vector<size_t>& fn)
        {
            auto pa = a.first.point_ids();
            auto pb = fn;
            std::sort(pa.begin(), pa.end());
            std::sort(pb.begin(), pb.end());
            return std::lexicographical_compare(pa.begin(), pa.end(), pb.begin(), pb.end());
        };

        auto eq = [](const std::pair<surface_type, size_t>& a, const std::vector<size_t>& fn)
        {
            auto pa = a.first.point_ids();
            auto pb = fn;
            std::sort(pa.begin(), pa.end());
            std::sort(pb.begin(), pb.end());
            return std::equal(pa.begin(), pa.end(), pb.begin());
        };


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

                    std::vector<size_t> fn{{node0_ofs, node1_ofs, node2_ofs}};

                    boundary_descriptor bi(b_id, tag, true);

                    auto itor = std::lower_bound(surf_map_table.begin(), surf_map_table.end(), fn, comp);
                    if ( (itor == surf_map_table.end()) or not eq(*itor, fn) )
                    {
                        std::cout << node0_ofs << " " << node1_ofs << " " << node2_ofs << std::endl;
                        throw std::logic_error("Face not found");
                    }

                    auto ofs = (*itor).second;
                    boundary_info.at(ofs) = bi;
                }

                nTags.clear();
                gmsh::model::mesh::getElementFaceNodes(elemType, 4, nTags, tag, true);

                for (size_t i = 0; i < nTags.size(); i+=4)
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

                    auto node3_tag = nTags[i + 3];
                    assert(node3_tag < node_tag2ofs.size());
                    auto node3_ofs = node_tag2ofs[node3_tag];
                    assert(node3_ofs != INVALID_OFS);

                    std::vector<size_t> fn{{node0_ofs, node1_ofs, node2_ofs, node3_ofs}};

                    boundary_descriptor bi(b_id, tag, true);

                    auto itor = std::lower_bound(surf_map_table.begin(), surf_map_table.end(), fn, comp);
                    if ( (itor == surf_map_table.end()) or not eq(*itor, fn) )
                        throw std::logic_error("Face not found");

                    auto ofs = (*itor).second;
                    boundary_info.at(ofs) = bi;
                }


            }

            b_id++;
        }
        surf_map_table.clear();
    }

public:
    static const char constexpr *expected_extension = "geo";
    gmsh_geometry_loader() = default;

    bool read_mesh(const std::string& s)
    {
        gmsh::initialize(0, nullptr);
        
        if (this->verbose())
            gmsh::option::setNumber("General.Terminal", 1);
        else
            gmsh::option::setNumber("General.Terminal", 0);
        
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
        vols.reserve( tmp_mesh_vols.size() );

        std::vector<subdomain_descriptor> subdoms;
        subdoms.reserve( tmp_mesh_vols.size() );

        for (auto& [vol, sd] : tmp_mesh_vols)
        {
            vols.push_back(vol);
            subdoms.push_back(sd);
        }

        tmp_mesh_vols.clear();

        storage->volumes = std::move(vols);
        storage->subdomain_info = std::move(subdoms);

        mark_internal_faces(msh);

        return true;
    }
};

/* I wrote this before discovering
 *
 *   Mesh.PartitionCreateTopology = 1;
 *
 * It essentially computes the boundaries between partitions
 * and adds them to the mesh stucture.
 * It can still be useful with other mesh formats...
 */
template<typename Mesh>
void make_interpartition_boundaries(Mesh& msh)
{
    auto cvf = disk::connectivity_via_faces(msh);

    using ipair = std::pair<size_t, size_t>;
    std::map<ipair, std::set<size_t>> newfaces;

    auto mkipair = [](size_t tag_a, size_t tag_b) {
        if (tag_a < tag_b)
            return std::pair(tag_a, tag_b);
        return std::pair(tag_b, tag_a);
    };

    for (auto& cl : msh) {
        auto di = msh.domain_info(cl);
        auto fcs = faces(msh, cl);
        for (auto& fc : fcs) {
            auto neigh = cvf.neighbour_via(msh, cl, fc);
            if (not neigh)
                continue;

            auto ncl = neigh.value();
            auto ndi = msh.domain_info(ncl);
            if ( di.tag() == ndi.tag() )
                continue;

            newfaces[mkipair(di.tag(), ndi.tag())].insert(offset(msh, fc));
        }
    }

    size_t maxbnd = 0;
    auto storage = msh.backend_storage();
    for (auto& bi : storage->boundary_info)
        maxbnd = std::max(maxbnd, bi.tag());
    
    for (auto& [ip, fcnums] : newfaces) {
        maxbnd++;
        for (auto& fcnum : fcnums) {
            if (not storage->boundary_info[fcnum].is_boundary())
                storage->boundary_info[fcnum] = {maxbnd, maxbnd, true, true};
        }
    }
}

template<disk::mesh_2D Mesh>
void dump_subdomain_boundaries(const Mesh& msh)
{
    using pt = typename Mesh::point_type;
    std::map<size_t, std::vector<std::pair<pt, pt>>> epts; 
    for (auto& cl : msh)
    {
        auto di = msh.domain_info(cl);

        auto fcs = faces(msh, cl);
        for (auto& fc : fcs)
        {
            auto bi = msh.boundary_info(fc);
            if (not bi.is_boundary())
                continue;
            //if (not bi.is_internal())
            //    continue;

            auto pts = points(msh, fc);
            epts[di.tag()].push_back({pts[0], pts[1]});
        }
    }

    std::cout << epts.size() << std::endl;

    for (auto& [tag, pts] : epts) {
        std::cout << pts.size() << std::endl;
        std::stringstream ss;
        ss << "interface_" << tag << ".txt";
        std::ofstream ofs(ss.str());
        for (auto& [p0,p1] : pts) {
            auto dp = p1 - p0;
            ofs << p0.x() << " " << p0.y() << " " << dp.x() << " " << dp.y() << std::endl;
        }
    }
}

} //namespace disk
