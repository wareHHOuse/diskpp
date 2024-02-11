/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023, 2024
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */
/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2020, 2021
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group
 */
/*
 *       /\        Matteo Cicuttin (C) 2016-2019
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

#include <tuple>
#include <fstream>
#include <vector>

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/geometry/geometry.hpp"
#include "diskpp/loaders/mesh_loader.hpp"

namespace disk
{

std::vector<size_t>
read_fvca6_line(std::ifstream& ifs)
{
    std::vector<size_t> ret;

    size_t num_entries;
    ifs >> num_entries;
    for(size_t j = 0; j < num_entries; j++)
    {
        size_t temp;
        ifs >> temp;
        ret.push_back(temp-1); //convert indices to zero-based
    }

    return ret;
}


template<typename T, size_t N>
class fvca6_mesh_loader
{
    static_assert(N == 3, "fvca6 is a 3D-only mesh format");
};

template<typename T>
class fvca6_mesh_loader<T,3> : public mesh_loader<generic_mesh<T, 3>>
{
    typedef generic_mesh<T,3>                       mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;
    typedef typename mesh_type::surface_type        surface_type;
    typedef typename mesh_type::volume_type         volume_type;

    std::vector<point_type>                         m_points;
    std::vector<node_type>                          m_nodes;
    std::vector<std::pair<size_t, edge_type>>       m_edges;

    std::vector<std::pair<size_t, std::vector<size_t>>>     vol_to_faces;
    std::vector<std::vector<size_t>>                        vol_to_vts;
    std::vector<std::pair<size_t, std::vector<size_t>>>     faces_to_edges;
    std::vector<std::vector<size_t>>                        faces_to_vts;

    bool fvca6_read(const std::string& filename)
    {
        std::ifstream   ifs(filename);
        std::string     keyword;
        size_t          lines_to_read;

        if (!ifs.is_open())
        {
            std::cout << "Error opening " << filename << std::endl;
            return false;
        }

        /* Apparently the first 16 lines of the file are comments or
         * information repeated elsewhere: throw them away */

        for (size_t i = 0; i < 16; i++)
            std::getline(ifs, keyword);

        if ( !expect(ifs, "Vertices") )
            return false;

        ifs >> lines_to_read;

        if (this->verbose())
            std::cout << "About to read " << lines_to_read << " points" << std::endl;

        m_points.reserve(lines_to_read);
        m_nodes.reserve(lines_to_read);
        for (size_t i = 0; i < lines_to_read; i++)
        {
            T x, y, z;
            ifs >> x >> y >> z;
            m_points.push_back( point_type({x, y, z}) );
            m_nodes.push_back( node_type( disk::point_identifier<3>(i) ) );
        }

        /* Volume to face data */
        if ( !expect(ifs, "Volumes->faces") )
            return false;

        ifs >> lines_to_read;

        if (this->verbose())
            std::cout << "About to read " << lines_to_read << " entries" << std::endl;

        for (size_t i = 0; i < lines_to_read; i++)
        {
            auto vol_faces = read_fvca6_line(ifs);
            vol_to_faces.push_back( std::make_pair(i, std::move(vol_faces)) );
        }

        /* Volume to vertices data */
        if ( !expect(ifs, "Volumes->Verticess") )
            return false;

        ifs >> lines_to_read;

        if (this->verbose())
            std::cout << "About to read " << lines_to_read << " entries" << std::endl;

        for (size_t i = 0; i < lines_to_read; i++)
        {
            auto vol_vts = read_fvca6_line(ifs);
            vol_to_vts.push_back( std::move(vol_vts) );
        }

        /* Faces to edges data */
        if ( !expect(ifs, "Faces->Edgess") )
            return false;

        ifs >> lines_to_read;

        if (this->verbose())
            std::cout << "About to read " << lines_to_read << " entries" << std::endl;

        for (size_t i = 0; i < lines_to_read; i++)
        {
            auto faces_edges = read_fvca6_line(ifs);
            faces_to_edges.push_back( std::make_pair(i, std::move(faces_edges)) );
        }

        /* Faces to vertices data */
        if ( !expect(ifs, "Faces->Vertices") )
            return false;

        ifs >> lines_to_read;

        if (this->verbose())
            std::cout << "About to read " << lines_to_read << " entries" << std::endl;

        for (size_t i = 0; i < lines_to_read; i++)
        {
            auto faces_vts = read_fvca6_line(ifs);
            faces_to_vts.push_back( std::move(faces_vts) );
        }

        /* Faces to cv data */
        if ( !expect(ifs, "Faces->Control") )
            return false;

        if ( !expect(ifs, "volumes") )
            return false;

        ifs >> lines_to_read;

        if (this->verbose())
            std::cout << "About to read " << lines_to_read << " entries" << std::endl;

        for (size_t i = 0; i < lines_to_read; i++)
        {
            int v1;
            int v2;
            ifs >> v1 >> v2;
        }

        /* Edges data */
        if ( !expect(ifs, "Edges") )
            return false;

        ifs >> lines_to_read;

        if (this->verbose())
            std::cout << "About to read " << lines_to_read << " entries" << std::endl;

        for (size_t i = 0; i < lines_to_read; i++)
        {
            size_t v1, v2;
            ifs >> v1 >> v2;

            if (v1 > v2)
                std::swap(v1, v2);

            auto e = edge_type({typename node_type::id_type(v1),
                                typename node_type::id_type(v2)});

            m_edges.push_back( std::make_pair(i, e) );
        }

        return true;
    }

public:
    static const char constexpr *expected_extension = "msh";

    fvca6_mesh_loader() = default;

    bool read_mesh(const std::string& s)
    {
        if (this->verbose())
            std::cout << " *** READING FVCA6 3D MESH ***" << std::endl;
        return fvca6_read(s);
    }

    bool populate_mesh(mesh_type& msh)
    {
        auto storage = msh.backend_storage();

        std::vector<size_t> conv_table;
        /* Sort the edges in lexicographical order, remember their original
         * position to convert the pointers in the faces */
        auto comp_edges = [](const std::pair<size_t, edge_type>& e1,
                             const std::pair<size_t, edge_type>& e2) {
            return e1.second < e2.second;
        };
        std::sort(m_edges.begin(), m_edges.end(), comp_edges);

        std::vector<edge_type> edges;
        edges.reserve( m_edges.size() );
        conv_table.resize( m_edges.size() );
        for (size_t i = 0; i < m_edges.size(); i++)
        {
            conv_table[m_edges[i].first] = i;   /* Make ptr conversion table */
            edges.push_back(m_edges[i].second);
        }

        /* Convert the edge pointers in the face data */
        for (auto& fe : faces_to_edges)
        {
            for (auto& ptr : fe.second)
                ptr = conv_table[ptr];
        }

        /* Sort in lexicographical order and remember original position */
        auto comp_vecs = [](const std::pair<size_t, std::vector<size_t>>& e1,
                            const std::pair<size_t, std::vector<size_t>>& e2) {
            return e1.second < e2.second;
        };
        std::sort(faces_to_edges.begin(), faces_to_edges.end(), comp_vecs);

        std::vector<surface_type> faces;
        faces.reserve( faces_to_edges.size() );
        conv_table.resize( faces_to_edges.size() );

        for (size_t i = 0; i < faces_to_edges.size(); i++)
        {
            auto fe = faces_to_edges[i];
            surface_type s( convert_to<typename edge_type::id_type>(fe.second) );
            s.set_point_ids( convert_to<disk::point_identifier<3>>(faces_to_vts.at(fe.first)) );
            faces.push_back(s);
            conv_table[fe.first] = i;
        }
        /* Now the faces are in their place and have correct ptrs */

        /* Convert the face pointers in the volume data */
        for (auto& vf : vol_to_faces )
        {
            for (auto& ptr : vf.second)
                ptr = conv_table[ptr];
        }

        //for (auto& f : faces)
        //    std::cout << f << std::endl;

        /* Sort volume data */
        std::sort(vol_to_faces.begin(), vol_to_faces.end(), comp_vecs);

        std::vector<volume_type> volumes;
        volumes.reserve( vol_to_faces.size() );

        for (size_t i = 0; i < vol_to_faces.size(); i++)
        {
            auto vf = vol_to_faces[i];
            volume_type v( convert_to<typename surface_type::id_type>(vf.second) );
            v.set_point_ids( convert_to<disk::point_identifier<3>>(vol_to_vts.at(vf.first)) );
            volumes.push_back(v);
        }

        storage->points     = std::move(m_points);
        storage->nodes      = std::move(m_nodes);
        storage->edges      = std::move(edges);
        storage->surfaces   = std::move(faces);
        storage->volumes    = std::move(volumes);

        std::vector<size_t> bf(storage->surfaces.size());
        for (auto& vol : storage->volumes)
        {
            for (auto itor = vol.subelement_id_begin(); itor != vol.subelement_id_end(); itor++)
                bf.at(*itor)++;
        }

        /* FVCA6 meshes do not provide boundary information. As usually they represent
         * the unit cube, we try to renumber according to the normals of the boundary
         * surfaces. If it fails, everything is marked as boundary 0.
         * mark_internal_faces() is not used here, as we wouldn't know which number to
         * give to the internal boundaries.
         */
        storage->boundary_info.resize(storage->surfaces.size());
        auto renumbered = renumber_hypercube_boundaries(msh);
        if (not renumbered) {
            std::cout << "Warning: unable to renumber FVCA6 mesh boundaries, ";
            std::cout << "defaulting everyting to 0.\nProbably one of the boundaries ";
            std::cout << "is not parallel to the xy, xz or yz planes.";
            std::cout << std::endl;
            boundary_descriptor bi(0, true);
            for (size_t i = 0; i < storage->surfaces.size(); i++)
                if (bf[i] == 1)
                    storage->boundary_info[i] = bi;
        }

        return true;
    }
};

}