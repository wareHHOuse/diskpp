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

 /*
  * Copyright (C) 2013-2016, Matteo Cicuttin - matteo.cicuttin@uniud.it
  * Department of Electrical Engineering, University of Udine
  * All rights reserved.
  *
  * Redistribution and use in source and binary forms, with or without
  * modification, are permitted provided that the following conditions are met:
  *     * Redistributions of source code must retain the above copyright
  *       notice, this list of conditions and the following disclaimer.
  *     * Redistributions in binary form must reproduce the above copyright
  *       notice, this list of conditions and the following disclaimer in the
  *       documentation and/or other materials provided with the distribution.
  *     * Neither the name of the University of Udine nor the
  *       names of its contributors may be used to endorse or promote products
  *       derived from this software without specific prior written permission.
  *
  * THIS SOFTWARE IS PROVIDED BY THE AUTHOR(s) ``AS IS'' AND ANY
  * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  * DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(s) BE LIABLE FOR ANY
  * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  */

#pragma once

#include <array>
#include <tuple>
#include <cassert>
#include <fstream>
#include <regex>
#include <set>
#include <thread>
#include <vector>

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/geometry/geometry.hpp"
#include "loader_utils.hpp"

#include "diskpp/common/mapped_file.h"
#include "strtot.hpp"

namespace disk {

template<typename mesh_type>
class mesh_loader
{
    bool    m_verbose;

public:
    mesh_loader()
        : m_verbose(false)
    {}

    virtual bool    read_mesh(const std::string&) { return false; }
    virtual bool    populate_mesh(mesh_type&)    = 0;

    bool    verbose(void) const     { return m_verbose; }
    void    verbose(bool v)         { m_verbose = v; }

    virtual ~mesh_loader() {}
};

template<typename T, size_t N>
class uniform_mesh_loader
{
    static_assert(N == 1, "at the moment only 1D uniform meshes are available");
};

template<typename T>
class uniform_mesh_loader<T,1> : public mesh_loader<generic_mesh<T, 1>>
{
    typedef generic_mesh<T,1>                       mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;

    T       m_h, m_x_min, m_x_max;
    size_t  m_number_of_elements;

public:
    uniform_mesh_loader()
        : m_x_min(T(0)), m_x_max(T(1)), m_number_of_elements(8)
    {
        m_h = fabs(m_x_max - m_x_min) / m_number_of_elements;
    }

    uniform_mesh_loader(T x_min, T x_max, size_t N)
        : m_x_min(x_min), m_x_max(x_max), m_number_of_elements(N)
    {
        m_h = fabs(m_x_max - m_x_min) / m_number_of_elements;
    }

    bool populate_mesh(mesh_type& msh)
    {
        if (this->verbose())
            std::cout << " *** POPULATING UNIFORM 1D MESH ***" << std::endl;

        auto storage = msh.backend_storage();

        auto num_edges = m_number_of_elements;
        auto num_nodes = m_number_of_elements + 1;

        storage->points.resize(num_nodes);
        storage->nodes.resize(num_nodes);
        storage->edges.resize(num_edges);

        for (size_t i = 0; i < num_edges; i++)
        {
            storage->points.at(i)   = point_type({m_x_min + (m_h * i)});
            storage->points.at(i+1) = point_type({m_x_min + (m_h * (i + 1))});

            auto n0 = typename node_type::id_type(i);
            auto n1 = typename node_type::id_type(i+1);
            auto e = edge_type(n0, n1);

            /*
            std::vector<point_identifier<1>> pts(2);
            pts[0] = point_identifier<1>(i);
            pts[1] = point_identifier<1>(i+1);
            e.set_point_ids(pts.begin(), pts.end());
            */
            storage->edges.at(i) = e;
        }

        for (size_t i = 0; i < num_nodes; i++)
            storage->nodes.at(i) = node_type(point_identifier<1>(i));

        storage->boundary_info.resize(num_nodes);
        boundary_descriptor bi(0, true);
        storage->boundary_info.at(0) = bi;
        storage->boundary_info.at(num_nodes - 1) = bi;

        return true;
    }
};

#include "loader_fvca5.hpp"


bool
expect(std::ifstream& ifs, const std::string& str)
{
    std::string keyword;
    ifs >> keyword;
    if ( keyword != str )
    {
        std::cout << "Expected keyword \"" << str << "\"" << std::endl;
        std::cout << "Found \"" << keyword << "\"" << std::endl;
        return false;
    }

    return true;
}

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
            size_t v1;
            int  v2;
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

        boundary_descriptor bi(0, true);
        storage->boundary_info.resize(storage->surfaces.size());
        for (size_t i = 0; i < storage->surfaces.size(); i++)
            if (bf[i] == 1)
                storage->boundary_info[i] = bi;

        return false;
    }
};










template<typename T, size_t N>
class netgen_mesh_loader
{
    static_assert(N == 2 || N == 3, "Netgen supports only 2D or 3D");
};

namespace priv {

template<typename T>
point<T,2>
read_2d_point_line(const char *str, char **endptr, T scalefactor)
{
    T t1, t2;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);

    return point<T,2>{t1*scalefactor, t2*scalefactor};
}

template<typename T>
point<T,3>
read_3d_point_line(const char *str, char **endptr, T scalefactor)
{
    T t1, t2, t3;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);

    return point<T,3>{t1*scalefactor, t2*scalefactor, t3*scalefactor};
}

template<typename T>
std::tuple<T, T, T, T, T>
read_tetrahedron_line(const char *str, char **endptr)
{
    T t1, t2, t3, t4, t5;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    t4 = strtot<T>(*endptr, endptr);
    t5 = strtot<T>(*endptr, endptr);

    return std::make_tuple(t1, t2-1, t3-1, t4-1, t5-1);
}

template<typename T>
std::tuple<T, T, T, T>
read_triangle_line(const char *str, char **endptr)
{
    T t1, t2, t3, t4;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    t4 = strtot<T>(*endptr, endptr);

    return std::make_tuple(t1, t2-1, t3-1, t4-1);
}

template<typename T>
std::tuple<T, T, T>
read_edge_line(const char *str, char **endptr)
{
    T t1, t2, t3;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);

    return std::make_tuple(t1, t2-1, t3-1);
}

} //namespace priv

template<typename T>
class netgen_mesh_loader<T,2> : public mesh_loader<simplicial_mesh<T,2>>
{
    typedef simplicial_mesh<T,2>                    mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;
    typedef typename mesh_type::surface_type        surface_type;

    std::vector<point_type>                         points;
    std::vector<node_type>                          nodes;
    std::vector<edge_type>                          edges;
    std::vector<std::pair<edge_type, size_t>>       boundary_edges;
    std::vector<surface_type>                       surfaces;

    bool netgen_read(const std::string& filename)
    {
        /* Open file */
        if (filename.size() == 0)
        {
            std::cout << "Invalid mesh file name" << std::endl;
            return false;
        }

        size_t lines, linecount;

        mapped_file mf(filename);

        //std::cout << green << " * * * Reading NETGEN format mesh * * * ";
        //std::cout << nocolor << std::endl;

        /************************ Read points ************************/
        linecount = 0;

        const char *data = mf.mem();
        char *endptr;

        lines = strtot<size_t>(data, &endptr);

        points.reserve(lines);
        nodes.reserve(lines);

        while (linecount < lines)
        {
            if ( this->verbose() && ((linecount%100000) == 0) )
            {
                std::cout << "Reading points: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto point = priv::read_2d_point_line<T>(endptr, &endptr, 1.0);

            points.push_back( point );

            auto point_id = disk::point_identifier<2>( linecount );
            auto node = node_type( { point_id } );

            nodes.push_back(node);

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading points: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        /************************ Read triangles ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        edges.reserve(lines*3);
        surfaces.reserve(lines);

        while (linecount < lines)
        {
            if ( this->verbose() && ((linecount%100000) == 0) )
            {
                std::cout << "Reading triangles: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_triangle_line<size_t>(endptr, &endptr);

            disk::point_identifier<2>     p0(std::get<1>(t));
            disk::point_identifier<2>     p1(std::get<2>(t));
            disk::point_identifier<2>     p2(std::get<3>(t));
            //domain_id_type      d(std::get<0>(t));

            edges.push_back( edge_type( { p0, p1 } ) );
            edges.push_back( edge_type( { p1, p2 } ) );
            edges.push_back( edge_type( { p0, p2 } ) );

            surfaces.push_back( surface_type( { p0, p1, p2 } ) );

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading triangles: " << linecount;
            std::cout << "/" << lines << std::endl;
        }
        /************************ Read boundary edges ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        boundary_edges.reserve(lines);

        while (linecount < lines)
        {
            if ( this->verbose() && ((linecount%50000) == 0) )
            {
                std::cout << "Reading edges: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_edge_line<size_t>(endptr, &endptr);

            disk::point_identifier<2>     p0(std::get<1>(t));
            disk::point_identifier<2>     p1(std::get<2>(t));

            edge_type   edge( { p0, p1 } );

            boundary_edges.push_back( std::make_pair(edge, std::get<0>(t)) );

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading edges: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        return true;
    }

public:
    static const char constexpr *expected_extension = "mesh2d";
    netgen_mesh_loader() = default;

    bool read_mesh(const std::string& s)
    {
        if (this->verbose())
            std::cout << " *** READING NETGEN 2D MESH ***" << std::endl;

        return netgen_read(s);
    }

    bool populate_mesh(mesh_type& msh)
    {
        auto storage = msh.backend_storage();

        if (this->verbose())
        {
            std::cout << "Sorting data...";
            std::cout.flush();
        }

        storage->points = std::move(points);
        storage->nodes = std::move(nodes);

        /* sort edges, make unique and move them in geometry */
        THREAD(edge_thread,
            priv::sort_uniq(edges);
            storage->edges = std::move(edges);
        );

        /* sort triangles, make unique and move them in geometry */
        THREAD(tri_thread,
            priv::sort_uniq(surfaces);
            storage->surfaces = std::move(surfaces);
        );

        /* wait for the threads */
        WAIT_THREAD(edge_thread);
        WAIT_THREAD(tri_thread);

        storage->boundary_info.resize(storage->edges.size());
        for (auto& be : boundary_edges)
        {
            auto position = find_element_id(storage->edges.begin(),
                                            storage->edges.end(), be.first);
            if (position.first == false)
            {
                std::cout << "Bad bug at " << __FILE__ << "("
                          << __LINE__ << ")" << std::endl;
                return false;
            }
            boundary_descriptor bi(be.second, true);
            storage->boundary_info.at(position.second) = bi;
        }

        if (this->verbose())
        {
            std::cout << "done." << std::endl;

            std::cout << "Nodes: " << storage->nodes.size() << std::endl;
            std::cout << "Edges: " << storage->edges.size() << std::endl;
            std::cout << "Faces: " << storage->surfaces.size() << std::endl;
        }

        boundary_edges.clear();

        storage->subdomain_info.resize(storage->surfaces.size());

        return true;
    }
};

template<typename T>
class netgen_mesh_loader<T,3> : public mesh_loader<simplicial_mesh<T,3>>
{
    typedef simplicial_mesh<T,3>                    mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;
    typedef typename mesh_type::surface_type        surface_type;
    typedef typename mesh_type::volume_type         volume_type;

    std::vector<point_type>                                     points;
    std::vector<node_type>                                      nodes;
    std::vector<edge_type>                                      edges;
    std::vector<surface_type>                                   surfaces;
    std::vector<boundary_descriptor>                            boundary_info;
    std::vector<std::pair<volume_type, subdomain_descriptor>>   tmp_volumes;


    bool netgen_read(const std::string& filename)
    {
        /* Open file */
        if (filename.size() == 0)
        {
            std::cout << "Invalid mesh file name" << std::endl;
            return false;
        }

        size_t lines, linecount;

        mapped_file mf(filename);

        /************************ Read points ************************/
        linecount = 0;

        const char *data = mf.mem();
        char *endptr;

        lines = strtot<size_t>(data, &endptr);

        points.reserve(lines);
        nodes.reserve(lines);

        while (linecount < lines)
        {
            if ( this->verbose() && ((linecount%100000) == 0) )
            {
                std::cout << "Reading points: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto point = priv::read_3d_point_line<T>(endptr, &endptr, 1.0);

            points.push_back( point );

            auto point_id = disk::point_identifier<3>( linecount );
            auto node = node_type( { point_id } );

            nodes.push_back(node);
            /* Do something with that point */

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading points: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        /************************ Read tetrahedra ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        edges.reserve(lines*6);
        surfaces.reserve(lines*4);
        tmp_volumes.reserve(lines);

        while (linecount < lines)
        {
            if ( this->verbose() && ((linecount%100000) == 0) )
            {
                std::cout << "Reading tetrahedra: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_tetrahedron_line<size_t>(endptr, &endptr);

            auto subdomain_num = std::get<0>(t);
            disk::point_identifier<3>     p0(std::get<1>(t));
            disk::point_identifier<3>     p1(std::get<2>(t));
            disk::point_identifier<3>     p2(std::get<3>(t));
            disk::point_identifier<3>     p3(std::get<4>(t));

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

            subdomain_descriptor si(subdomain_num);
            auto vs = std::make_pair(volume_type( {p0, p1, p2, p3} ), si);
            tmp_volumes.push_back( vs );

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading tetrahedra: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        priv::sort_uniq(edges);
        priv::sort_uniq(surfaces);
        boundary_info.resize(surfaces.size());

        using pvs = std::pair<volume_type, subdomain_descriptor>;
        auto tmpvol_comp = [](const pvs& a, const pvs& b) {
            return a.first < b.first;
        };

        std::sort(tmp_volumes.begin(), tmp_volumes.end(), tmpvol_comp);

        /************************ Read boundary surfaces ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);
        while (linecount < lines)
        {
            if ( this->verbose() && ((linecount%50000) == 0) )
            {
                std::cout << "Reading triangle: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_triangle_line<size_t>(endptr, &endptr);

            auto bnd_id = std::get<0>(t);
            disk::point_identifier<3>     p0(std::get<1>(t));
            disk::point_identifier<3>     p1(std::get<2>(t));
            disk::point_identifier<3>     p2(std::get<3>(t));

            boundary_descriptor bi(bnd_id, true);
            surface_type surf( { p0, p1, p2 } );

            auto itor = std::lower_bound(surfaces.begin(), surfaces.end(), surf);
            if ( (itor == surfaces.end()) or not (*itor == surf) )
                throw std::logic_error("Face not found");

            auto ofs = std::distance(surfaces.begin(), itor);
            boundary_info.at(ofs) = bi;

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading triangle: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        return true;
    }

public:
    static const char constexpr *expected_extension = "mesh";
    netgen_mesh_loader() = default;

    bool read_mesh(const std::string& s)
    {
        if (this->verbose())
            std::cout << " *** READING NETGEN 3D MESH ***" << std::endl;
        return netgen_read(s);
    }

    bool populate_mesh(mesh_type& msh)
    {
        auto storage = msh.backend_storage();

        if (this->verbose())
        {
            std::cout << "Sorting data...";
            std::cout.flush();
        }

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

        if (this->verbose())
        {
            std::cout << "done." << std::endl;

            std::cout << "Nodes: " << storage->nodes.size() << std::endl;
            std::cout << "Edges: " << storage->edges.size() << std::endl;
            std::cout << "Faces: " << storage->surfaces.size() << std::endl;
            std::cout << "Volumes: " << storage->volumes.size() << std::endl;
        }

        return true;
    }

};


namespace priv {

template<typename T>
std::tuple<T, T, T, T, T, T, T, T>
read_hexahedron_line(const char *str, char **endptr)
{
    T t1, t2, t3, t4, t5, t6, t7, t8;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    t4 = strtot<T>(*endptr, endptr);
    t5 = strtot<T>(*endptr, endptr);
    t6 = strtot<T>(*endptr, endptr);
    t7 = strtot<T>(*endptr, endptr);
    t8 = strtot<T>(*endptr, endptr);

    return std::make_tuple(t1, t2, t3, t4, t5, t6, t7, t8);
}

template<typename T>
std::tuple<T, T, T, T>
read_hex_face_line(const char *str, char **endptr)
{
    T t1, t2, t3, t4;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    t4 = strtot<T>(*endptr, endptr);

    return std::make_tuple(t1, t2, t3, t4);
}

template<typename T>
std::tuple<T, T, T, T>
read_quad_line(const char *str, char **endptr)
{
    T t1, t2, t3, t4;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    t4 = strtot<T>(*endptr, endptr);

    return std::make_tuple(t1, t2, t3, t4);
}

template<typename T>
std::tuple<T, T>
read_quad_face_line(const char *str, char **endptr)
{
    T t1, t2;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);

    return std::make_tuple(t1, t2);
}

} // namespace priv


template<typename T, size_t DIM>
class cartesian_mesh_loader;

template<typename T>
class cartesian_mesh_loader<T,3> : public mesh_loader<cartesian_mesh<T,3>>
{
    typedef cartesian_mesh<T,3>                     mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;
    typedef typename mesh_type::surface_type        surface_type;
    typedef typename mesh_type::volume_type         volume_type;

    std::vector<point_type>                         points;
    std::vector<node_type>                          nodes;
    std::vector<edge_type>                          edges;
    std::vector<surface_type>                       surfaces, boundary_surfaces;
    std::vector<volume_type>                        volumes;


    bool hex_read(const std::string& filename)
    {
        /* Open file */
        if (filename.size() == 0)
        {
            std::cout << "Invalid mesh file name" << std::endl;
            return false;
        }

        size_t lines, linecount;

        mapped_file mf(filename);

        //std::cout << green << " * * * Reading NETGEN format mesh * * * ";
        //std::cout << nocolor << std::endl;

        /************************ Read points ************************/
        linecount = 0;

        const char *data = mf.mem();
        char *endptr;

        lines = strtot<size_t>(data, &endptr);

        points.reserve(lines);
        nodes.reserve(lines);

        while (linecount < lines)
        {
            if ( this->verbose() && ((linecount%100000) == 0) )
            {
                std::cout << "Reading points: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto point = priv::read_3d_point_line<T>(endptr, &endptr, 1.0);

            points.push_back( point );

            auto point_id = disk::point_identifier<3>( linecount );
            auto node = node_type( { point_id } );

            nodes.push_back(node);
            /* Do something with that point */

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading points: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        /************************ Read hexahedra ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        edges.reserve(lines*12);
        surfaces.reserve(lines*6);
        volumes.reserve(lines);

        while (linecount < lines)
        {
            if ( this->verbose() && ((linecount%100000) == 0) )
            {
                std::cout << "Reading hexahedra: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_hexahedron_line<size_t>(endptr, &endptr);

            disk::point_identifier<3>     p0(std::get<0>(t));
            disk::point_identifier<3>     p1(std::get<1>(t));
            disk::point_identifier<3>     p2(std::get<2>(t));
            disk::point_identifier<3>     p3(std::get<3>(t));
            disk::point_identifier<3>     p4(std::get<4>(t));
            disk::point_identifier<3>     p5(std::get<5>(t));
            disk::point_identifier<3>     p6(std::get<6>(t));
            disk::point_identifier<3>     p7(std::get<7>(t));

            edges.push_back( edge_type( { p0, p1 } ) );
            edges.push_back( edge_type( { p0, p2 } ) );
            edges.push_back( edge_type( { p0, p4 } ) );
            edges.push_back( edge_type( { p1, p3 } ) );
            edges.push_back( edge_type( { p1, p5 } ) );
            edges.push_back( edge_type( { p2, p3 } ) );
            edges.push_back( edge_type( { p2, p6 } ) );
            edges.push_back( edge_type( { p3, p7 } ) );
            edges.push_back( edge_type( { p4, p5 } ) );
            edges.push_back( edge_type( { p4, p6 } ) );
            edges.push_back( edge_type( { p5, p7 } ) );
            edges.push_back( edge_type( { p6, p7 } ) );

            surfaces.push_back( surface_type( { p0, p2, p6, p4 } ) );
            surfaces.push_back( surface_type( { p1, p3, p7, p5 } ) );
            surfaces.push_back( surface_type( { p0, p1, p3, p2 } ) );
            surfaces.push_back( surface_type( { p4, p5, p7, p6 } ) );
            surfaces.push_back( surface_type( { p0, p4, p5, p1 } ) );
            surfaces.push_back( surface_type( { p2, p6, p7, p3 } ) );

            volumes.push_back( volume_type( { p0, p1, p2, p3, p4, p5, p6, p7 } ) );

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading hexahedra: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        /************************ Read boundary surfaces ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        boundary_surfaces.reserve(lines);

        while (linecount < lines)
        {
            if ( this->verbose() && ((linecount%50000) == 0) )
            {
                std::cout << "Reading hex face: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_hex_face_line<size_t>(endptr, &endptr);

            disk::point_identifier<3>     p0(std::get<0>(t));
            disk::point_identifier<3>     p1(std::get<1>(t));
            disk::point_identifier<3>     p2(std::get<2>(t));
            disk::point_identifier<3>     p3(std::get<3>(t));

            surface_type   quad( { p0, p1, p2, p3 } );

            boundary_surfaces.push_back( quad );

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading hex face: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        return true;
    }

public:
    static const char constexpr *expected_extension = "hex";
    cartesian_mesh_loader() = default;

    bool read_mesh(const std::string& s)
    {
        if (this->verbose())
            std::cout << " *** READING CARTESIAN 3D MESH ***" << std::endl;

        return hex_read(s);
    }

    bool populate_mesh(mesh_type& msh)
    {
        auto storage = msh.backend_storage();

        if (this->verbose())
        {
            std::cout << "Sorting data...";
            std::cout.flush();
        }

        storage->points = std::move(points);
        storage->nodes = std::move(nodes);

        /* sort edges, make unique and move them in geometry */
        THREAD(edge_thread,
            priv::sort_uniq(edges);
            storage->edges = std::move(edges);
        );

        /* sort triangles, make unique and move them in geometry */
        THREAD(quad_thread,
            priv::sort_uniq(surfaces);
            storage->surfaces = std::move(surfaces);
        );

        /* sort tetrahedra, make unique and move them in geometry */
        THREAD(hex_thread,
            std::sort(volumes.begin(), volumes.end());
            storage->volumes = std::move(volumes);
        );

        /* wait for the threads */
        WAIT_THREAD(edge_thread);
        WAIT_THREAD(quad_thread);
        WAIT_THREAD(hex_thread);

        storage->boundary_info.resize(storage->surfaces.size());
        for (auto& bs : boundary_surfaces)
        {
            auto position = find_element_id(storage->surfaces.begin(),
                                            storage->surfaces.end(), bs);
            if (position.first == false)
            {
                std::cout << "Bad bug at " << __FILE__ << "("
                          << __LINE__ << ")" << std::endl;
                return false;
            }

            boundary_descriptor bi(0, true);
            storage->boundary_info.at(position.second) = bi;
        }

        if (this->verbose())
        {
            std::cout << "done." << std::endl;

            std::cout << "Nodes: " << storage->nodes.size() << std::endl;
            std::cout << "Edges: " << storage->edges.size() << std::endl;
            std::cout << "Faces: " << storage->surfaces.size() << std::endl;
            std::cout << "Volumes: " << storage->volumes.size() << std::endl;
        }

        boundary_surfaces.clear();

        return true;
    }

};

template<typename T>
class cartesian_mesh_loader<T,2> : public mesh_loader<cartesian_mesh<T,2>>
{
    typedef cartesian_mesh<T,2>                     mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;
    typedef typename mesh_type::surface_type        surface_type;

    std::vector<point_type>                         points;
    std::vector<node_type>                          nodes;
    std::vector<edge_type>                          edges, boundary_edges;
    std::vector<surface_type>                       surfaces;


    bool quad_read(const std::string& filename)
    {
        /* Open file */
        if (filename.size() == 0)
        {
            std::cout << "Can't open '" << filename << "'" << std::endl;
            return false;
        }

        size_t lines, linecount;

        mapped_file mf(filename);
        if (!mf.is_open())
            return false;

        /************************ Read points ************************/
        linecount = 0;

        const char *data = mf.mem();
        char *endptr;

        lines = strtot<size_t>(data, &endptr);

        points.reserve(lines);
        nodes.reserve(lines);

        while (linecount < lines)
        {
            if (  this->verbose() && (linecount%100000) == 0 )
            {
                std::cout << "Reading points: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto point = priv::read_2d_point_line<T>(endptr, &endptr, 1.0);

            points.push_back( point );

            auto point_id = disk::point_identifier<2>( linecount );
            auto node = node_type( { point_id } );

            nodes.push_back(node);
            /* Do something with that point */

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading points: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        /************************ Read hexahedra ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        edges.reserve(lines*4);
        surfaces.reserve(lines);

        while (linecount < lines)
        {
            if (  this->verbose() && (linecount%100000) == 0 )
            {
                std::cout << "Reading quads: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_quad_line<size_t>(endptr, &endptr);

            disk::point_identifier<2>     p0(std::get<0>(t));
            disk::point_identifier<2>     p1(std::get<1>(t));
            disk::point_identifier<2>     p2(std::get<2>(t));
            disk::point_identifier<2>     p3(std::get<3>(t));

            edges.push_back( edge_type( { p0, p1 } ) );
            edges.push_back( edge_type( { p0, p2 } ) );
            edges.push_back( edge_type( { p1, p3 } ) );
            edges.push_back( edge_type( { p2, p3 } ) );

            surfaces.push_back( surface_type( { p0, p1, p2, p3 } ) );

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading quads: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        /************************ Read boundary surfaces ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        boundary_edges.reserve(lines);

        while (linecount < lines)
        {
            if (  this->verbose() && (linecount%50000) == 0 )
            {
                std::cout << "Reading faces: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_quad_face_line<size_t>(endptr, &endptr);

            disk::point_identifier<2>     p0(std::get<0>(t));
            disk::point_identifier<2>     p1(std::get<1>(t));

            edge_type   bnd( { p0, p1 } );

            boundary_edges.push_back( bnd );

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading faces: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        return true;
    }

public:
    static const char constexpr *expected_extension = "quad";

    cartesian_mesh_loader() = default;

    bool read_mesh(const std::string& s)
    {
       if (this->verbose())
         std::cout << " *** READING CARTESIAN 2D MESH ***" << std::endl;

        return quad_read(s);
    }

    bool populate_mesh(mesh_type& msh)
    {
        auto storage = msh.backend_storage();

        if (this->verbose())
        {
           std::cout << "Sorting data...";
           std::cout.flush();
        }

        storage->points = std::move(points);
        storage->nodes = std::move(nodes);

        /* sort edges, make unique and move them in geometry */
        THREAD(edge_thread,
            priv::sort_uniq(edges);
            storage->edges = std::move(edges);
        );

        /* sort triangles, make unique and move them in geometry */
        THREAD(quad_thread,
            priv::sort_uniq(surfaces);
            storage->surfaces = std::move(surfaces);
        );

        /* wait for the threads */
        WAIT_THREAD(edge_thread);
        WAIT_THREAD(quad_thread);

        storage->boundary_info.resize(storage->edges.size());
        for (auto& bs : boundary_edges)
        {
            auto position = find_element_id(storage->edges.begin(),
                                            storage->edges.end(), bs);
            if (position.first == false)
            {
                std::cout << "Bad bug at " << __FILE__ << "("
                          << __LINE__ << ")" << std::endl;
                return false;
            }

            boundary_descriptor bi(0, true);
            storage->boundary_info.at(position.second) = bi;
        }

        if (this->verbose())
        {
        std::cout << "done." << std::endl;

        std::cout << "Nodes: " << storage->nodes.size() << std::endl;
        std::cout << "Edges: " << storage->edges.size() << std::endl;
        std::cout << "Faces: " << storage->surfaces.size() << std::endl;
        }

        boundary_edges.clear();

        return true;
    }

};

template<typename T, size_t N>
class medit_mesh_loader
{
   static_assert(N == 2 || N == 3, "Medit supports only 2D and 3D for the moment");
};

template<typename T>
class medit_mesh_loader<T, 2> : public mesh_loader<generic_mesh<T, 2>>
{
   typedef generic_mesh<T, 2>               mesh_type;
   typedef typename mesh_type::point_type   point_type;
   typedef typename mesh_type::node_type    node_type;
   typedef typename mesh_type::edge_type    edge_type;
   typedef typename mesh_type::surface_type surface_type;

   struct medit2d_poly
   {
      std::vector<size_t>             nodes;
      size_t                          id;
      std::set<std::array<size_t, 2>> attached_edges;

      bool
      operator<(const medit2d_poly& other)
      {
         return nodes < other.nodes;
      }
   };

   std::vector<point_type>                                     m_points;
   std::vector<node_type>                                      m_nodes;
   std::vector<std::array<ident_raw_t, 2>>                    m_edges;
   std::vector<medit2d_poly>                                   m_polys;
   std::vector<std::pair<std::array<ident_raw_t, 2>, size_t>> m_boundary_edges;

   bool
   medit_read_vertices(std::ifstream& ifs)
   {
      size_t elements_to_read;
      T      x, y, z;
      T      id;

      ifs >> elements_to_read;

      if (this->verbose())
         std::cout << "Attempting to read " << elements_to_read << " points" << std::endl;

      m_points.reserve(elements_to_read);
      m_nodes.reserve(elements_to_read);

      for (size_t i = 0; i < elements_to_read; i++) {
         ifs >> x >> y >> z >> id;
         m_points.push_back(point_type{x, y});
         m_nodes.push_back(node_type(disk::point_identifier<2>(i)));
      }

      return true;
   }

   bool
   medit_read_polygons(std::ifstream& ifs, size_t polynum)
   {
      size_t elements_to_read;

      ifs >> elements_to_read;
      if (this->verbose())
         std::cout << "Reading " << elements_to_read << " " << polynum << "-angles" << std::endl;

      for (size_t i = 0; i < elements_to_read; i++) {
         medit2d_poly              p;
         std::vector<ident_raw_t> nodes(polynum + 1, 0);

         for (size_t j = 0; j < polynum; j++) {
            ident_raw_t val;
            ifs >> val;
            p.nodes.push_back(val - 1);
            nodes[j] = val - 1;
         }
         nodes[polynum] = nodes[0];

         ifs >> p.id;

         m_polys.push_back(p);

         size_t p_id = m_polys.size();

         // We have too create edges
         for (size_t j = 0; j < polynum; j++) {
            std::array<ident_raw_t, 2> b_edge = {nodes[j], nodes[j + 1]};
            assert(b_edge[0] != b_edge[1]);
            if (b_edge[0] > b_edge[1]) std::swap(b_edge[0], b_edge[1]);

            m_edges.push_back(b_edge);
            m_polys.at(p_id - 1).attached_edges.insert({b_edge[0], b_edge[1]});
         }
      }

      return true;
   }

   bool
   medit_read_boundary_edges(std::ifstream& ifs)
   {
      size_t elements_to_read;

      ifs >> elements_to_read;
      if (this->verbose())
         std::cout << "Reading " << elements_to_read << " boundary edges" << std::endl;

      m_boundary_edges.reserve(elements_to_read);

      for (size_t i = 0; i < elements_to_read; i++) {
         std::array<ident_raw_t, 2> b_edge;
         ifs >> b_edge[0];
         b_edge[0] -= 1;
         ifs >> b_edge[1];
         b_edge[1] -= 1;

         assert(b_edge[0] != b_edge[1]);

         if (b_edge[0] > b_edge[1]) std::swap(b_edge[0], b_edge[1]);

         std::array<ident_raw_t, 2> bnd = {b_edge[0], b_edge[1]};

         size_t b_id;
         ifs >> b_id;

         m_boundary_edges.push_back(std::make_pair(bnd, b_id));
      }

      return true;
   }

   bool
   face_unique(const std::array<ident_raw_t, 2>& f1, const std::array<ident_raw_t, 2>& f2)
   {
      if (f1[0] == f2[0]) {
         if (f1[1] == f2[1]) {
            return false;
         } else
            return true;
      } else
         return true;
   }

   bool
   medit_read(const std::string& filename)
   {
      std::ifstream ifs(filename);
      std::string   keyword;

      if (!ifs.is_open()) {
         std::cout << "Error opening " << filename << std::endl;
         return false;
      }

      ifs >> keyword;
      if (keyword != "MeshVersionFormatted") {
         std::cout << "Expected keyword \"MeshVersionFormatted\"" << std::endl;
         return false;
      }

      size_t format;
      ifs >> format;

      if (format != 2) {
         std::cout << "Expected format 2 (here: " << format << ")" << std::endl;
         return false;
      }

      ifs >> keyword;
      if (keyword != "Dimension") {
         std::cout << "Expected keyword \"Dimension\"" << std::endl;
         return false;
      }

      size_t dim;
      ifs >> dim;

      if (dim != 3) {
         std::cout << "Expected dimension >=2 (here: " << dim << ")" << std::endl;
         return false;
      }

      ifs >> keyword;
      while (keyword != "End") {
         if (keyword == "Vertices") {
            medit_read_vertices(ifs);
         } else if (keyword == "Triangles") {
            medit_read_polygons(ifs, 3);
         } else if (keyword == "Quadrilaterals") {
            medit_read_polygons(ifs, 4);
         } else if (keyword == "Pentagons") {
            medit_read_polygons(ifs, 5);
         } else if (keyword == "Hexagons") {
            medit_read_polygons(ifs, 6);
         } else if (keyword == "Edges") {
            m_boundary_edges.clear();
            medit_read_boundary_edges(ifs);
         } else {
            std::cout << "Error parsing Medit file" << std::endl;
            return false;
         }

         ifs >> keyword;
      }

      ifs.close();
      return true;
   }

 public:
   static const char constexpr *expected_extension = "medit2d";
   medit_mesh_loader() = default;

   bool
   read_mesh(const std::string& s)
   {
      if (this->verbose()) std::cout << " *** READING MEDIT 2D MESH ***" << std::endl;

      return medit_read(s);
   }

   bool
   populate_mesh(mesh_type& msh)
   {
      if (this->verbose()) std::cout << " *** POPULATING MEDIT MESH ***" << std::endl;
      auto storage = msh.backend_storage();

      /* Points */
      storage->points = std::move(m_points);
      storage->nodes  = std::move(m_nodes);

      /* Edges */
      /* Make the vector containing the edges */
      std::vector<edge_type> edges;
      edges.reserve(m_edges.size());
      size_t nb_edge(0);
      for (size_t i = 0; i < m_edges.size(); i++) {
         bool unique = true;
         for (size_t j = 0; j < i; j++) {
            if (!face_unique(m_edges[i], m_edges[j])) {
               unique = false;
               break;
            }
         }

         if (unique) {
            assert(m_edges[i][0] < m_edges[i][1]);
            auto node1 = typename node_type::id_type(m_edges[i][0]);
            auto node2 = typename node_type::id_type(m_edges[i][1]);

            auto e = edge_type(node1, node2);

            /* Next line not necessary anymore, see generic_element<DIM, DIM-1> */
            //e.set_point_ids(m_edges[i].begin(), m_edges[i].begin() + 2); /* XXX: crap */
            edges.push_back(e);
            nb_edge++;
         }
      }
      /* Sort them */
      edges.resize(nb_edge);
      std::sort(edges.begin(), edges.end());

      /* Detect which ones are boundary edges */
      storage->boundary_info.resize(edges.size());
      for (size_t i = 0; i < m_boundary_edges.size(); i++) {
         assert(m_boundary_edges[i].first[0] < m_boundary_edges[i].first[1]);
         auto node1 = typename node_type::id_type(m_boundary_edges[i].first[0]);
         auto node2 = typename node_type::id_type(m_boundary_edges[i].first[1]);

         auto e = edge_type(node1, node2);

         auto position = find_element_id(edges.begin(), edges.end(), e);

         if (position.first == false) {
            std::cout << "Bad bug at " << __FILE__ << "(" << __LINE__ << ")" << std::endl;
            return false;
         } else {
            boundary_descriptor bi(m_boundary_edges[i].second, true);
            storage->boundary_info.at(position.second) = bi;
         }
      }

      storage->edges = std::move(edges);

      /* Surfaces */
      std::vector<surface_type> surfaces;
      surfaces.reserve(m_polys.size());

      for (auto& p : m_polys) {
         std::vector<typename edge_type::id_type> surface_edges;
         for (auto& e : p.attached_edges) {
            assert(e[0] < e[1]);
            auto n1 = typename node_type::id_type(e[0]);
            auto n2 = typename node_type::id_type(e[1]);

            edge_type edge(n1, n2);
            auto      edge_id = find_element_id(storage->edges.begin(), storage->edges.end(), edge);
            if (!edge_id.first) {
               std::cout << "Bad bug at " << __FILE__ << "(" << __LINE__ << ")" << std::endl;
               return false;
            }

            surface_edges.push_back(edge_id.second);
         }
         auto surface = surface_type(surface_edges);
         surface.set_point_ids(p.nodes.begin(), p.nodes.end()); /* XXX: crap */
         surfaces.push_back(surface);
      }

      std::sort(surfaces.begin(), surfaces.end());
      storage->surfaces = std::move(surfaces);

      return true;
   }
};

template<typename T>
class medit_mesh_loader<T, 3> : public mesh_loader<generic_mesh<T, 3>>
{
   typedef generic_mesh<T, 3>               mesh_type;
   typedef typename mesh_type::point_type   point_type;
   typedef typename mesh_type::node_type    node_type;
   typedef typename mesh_type::edge_type    edge_type;
   typedef typename mesh_type::surface_type surface_type;
   typedef typename mesh_type::volume_type  volume_type;

   std::vector<point_type>                   m_points;
   std::vector<node_type>                    m_nodes;
   std::vector<std::pair<size_t, edge_type>> m_edges;
   std::vector<std::array<size_t, 2>>        m_boundary_edges;

   std::vector<std::pair<size_t, std::vector<size_t>>> vol_to_faces;
   std::vector<std::vector<size_t>>                    vol_to_vts;
   std::vector<std::pair<size_t, std::vector<size_t>>> faces_to_edges;
   std::vector<std::vector<size_t>>                    faces_to_vts;

   bool
   expect(std::ifstream& ifs, const std::string& str)
   {
      std::string keyword;
      ifs >> keyword;
      if (keyword != str) {
         std::cout << "Expected keyword \"" << str << "\"" << std::endl;
         std::cout << "Found \"" << keyword << "\"" << std::endl;
         return false;
      }

      return true;
   }

   std::vector<size_t>
   read_medit3d_line(std::ifstream& ifs)
   {
      std::vector<size_t> ret;

      size_t num_entries;
      ifs >> num_entries;
      for (size_t j = 0; j < num_entries; j++) {
         size_t temp;
         ifs >> temp;
         ret.push_back(temp);
      }

      return ret;
   }

   bool
   medit3d_read(const std::string& filename)
   {
      std::ifstream ifs(filename);
      std::string   keyword;
      size_t        lines_to_read;

      if (!ifs.is_open()) {
         std::cout << "Error opening " << filename << std::endl;
         return false;
      }

      ifs >> keyword;
      if (keyword != "MeshVersionFormatted") {
         std::cout << "Expected keyword \"MeshVersionFormatted\"" << std::endl;
         return false;
      }

      size_t format;
      ifs >> format;

      if (format != 2) {
         std::cout << "Expected format 2 (here: " << format << ")" << std::endl;
         return false;
      }

      ifs >> keyword;
      if (keyword != "Dimension") {
         std::cout << "Expected keyword \"Dimension\"" << std::endl;
         return false;
      }

      size_t dim;
      ifs >> dim;

      if (dim != 3) {
         std::cout << "Expected dimension == 3 (here: " << dim << ")" << std::endl;
         return false;
      }

      if (!expect(ifs, "Vertices")) return false;

      ifs >> lines_to_read;

      if (this->verbose()) std::cout << "About to read " << lines_to_read << " points" << std::endl;

      m_points.reserve(lines_to_read);
      m_nodes.reserve(lines_to_read);
      for (size_t i = 0; i < lines_to_read; i++) {
         T x, y, z;
         ifs >> x >> y >> z;
         m_points.push_back(point_type({x, y, z}));
         m_nodes.push_back(node_type(disk::point_identifier<3>(i)));
      }

      /* Volume to face data */
      if (!expect(ifs, "Volumes->Faces")) return false;

      ifs >> lines_to_read;

      if (this->verbose())
         std::cout << "About to read " << lines_to_read << " entries" << std::endl;

      for (size_t i = 0; i < lines_to_read; i++) {
         auto vol_faces = read_medit3d_line(ifs);
         vol_to_faces.push_back(std::make_pair(i, std::move(vol_faces)));
      }

      /* Volume to vertices data */
      if (!expect(ifs, "Volumes->Vertices")) return false;

      ifs >> lines_to_read;

      if (this->verbose())
         std::cout << "About to read " << lines_to_read << " entries" << std::endl;

      for (size_t i = 0; i < lines_to_read; i++) {
         auto vol_vts = read_medit3d_line(ifs);
         vol_to_vts.push_back(std::move(vol_vts));
      }

      /* Faces to edges data */
      if (!expect(ifs, "Faces->Edges")) return false;

      ifs >> lines_to_read;

      if (this->verbose())
         std::cout << "About to read " << lines_to_read << " entries" << std::endl;

      for (size_t i = 0; i < lines_to_read; i++) {
         auto faces_edges = read_medit3d_line(ifs);
         faces_to_edges.push_back(std::make_pair(i, std::move(faces_edges)));
      }

      /* Faces to vertices data */
      if (!expect(ifs, "Faces->Vertices")) return false;

      ifs >> lines_to_read;

      if (this->verbose())
         std::cout << "About to read " << lines_to_read << " entries" << std::endl;

      for (size_t i = 0; i < lines_to_read; i++) {
         auto faces_vts = read_medit3d_line(ifs);
         faces_to_vts.push_back(std::move(faces_vts));
      }

      /* Faces to cv data */
      if (!expect(ifs, "Faces->Control")) return false;

      if (!expect(ifs, "volumes")) return false;

      ifs >> lines_to_read;

      if (this->verbose())
         std::cout << "About to read " << lines_to_read << " entries" << std::endl;

      for (size_t i = 0; i < lines_to_read; i++) {
         size_t num, id;
         ifs >> num >> id;
         m_boundary_edges.push_back({num, id});
      }

      /* Edges data */
      if (!expect(ifs, "Edges")) return false;

      ifs >> lines_to_read;

      if (this->verbose())
         std::cout << "About to read " << lines_to_read << " entries" << std::endl;

      for (size_t i = 0; i < lines_to_read; i++) {
         size_t v1, v2;
         ifs >> v1 >> v2;

         if (v1 > v2) std::swap(v1, v2);

         auto e = edge_type({typename node_type::id_type(v1), typename node_type::id_type(v2)});

         m_edges.push_back(std::make_pair(i, e));
      }

      return true;
   }

 public:
   static const char constexpr *expected_extension = "medit3d";
   medit_mesh_loader() = default;

   bool
   read_mesh(const std::string& s)
   {
      if (this->verbose()) std::cout << " *** READING MEDIT 3D MESH ***" << std::endl;
      return medit3d_read(s);
   }

   bool
   populate_mesh(mesh_type& msh)
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
      edges.reserve(m_edges.size());
      conv_table.resize(m_edges.size());
      for (size_t i = 0; i < m_edges.size(); i++) {
         conv_table[m_edges[i].first] = i; /* Make ptr conversion table */
         edges.push_back(m_edges[i].second);
      }

      /* Convert the edge pointers in the face data */
      for (auto& fe : faces_to_edges) {
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
      faces.reserve(faces_to_edges.size());
      conv_table.resize(faces_to_edges.size());

      for (size_t i = 0; i < faces_to_edges.size(); i++) {
         auto         fe = faces_to_edges[i];
         surface_type s(convert_to<typename edge_type::id_type>(fe.second));
         s.set_point_ids(convert_to<disk::point_identifier<3>>(faces_to_vts.at(fe.first)));
         faces.push_back(s);
         conv_table[fe.first] = i;
      }
      /* Now the faces are in their place and have correct ptrs */

      /* Convert the face pointers in the volume data */
      for (auto& vf : vol_to_faces) {
         for (auto& ptr : vf.second)
            ptr = conv_table[ptr];
      }

      /* Sort volume data */
      std::sort(vol_to_faces.begin(), vol_to_faces.end(), comp_vecs);

      std::vector<volume_type> volumes;
      volumes.reserve(vol_to_faces.size());

      for (size_t i = 0; i < vol_to_faces.size(); i++) {
         auto        vf = vol_to_faces[i];
         volume_type v(convert_to<typename surface_type::id_type>(vf.second));
         v.set_point_ids(convert_to<disk::point_identifier<3>>(vol_to_vts.at(vf.first)));
         volumes.push_back(v);
      }

      storage->points   = std::move(m_points);
      storage->nodes    = std::move(m_nodes);
      storage->edges    = std::move(edges);
      storage->surfaces = std::move(faces);
      storage->volumes  = std::move(volumes);

      storage->boundary_info.resize(storage->surfaces.size());
      for (size_t i = 0; i < m_boundary_edges.size(); i++) {
         boundary_descriptor bi(m_boundary_edges[i][1], true);
         storage->boundary_info[conv_table[m_boundary_edges[i][0]]] = bi;
      }

      return false;
   }
};

template<typename T, size_t N>
class poly_mesh_loader
{
    static_assert(N == 2 || N == 3, "Medit supports only 2D and 3D for the moment");
};

namespace poly{


std::pair<size_t, std::vector<size_t>>
read_cell_line(std::ifstream& ifs)
{
    std::vector<size_t> elem;
    size_t              cell_id, nb_elem;

    ifs >> cell_id;
    ifs >> nb_elem;
    // std::cout << cell_id << ", " << nb_elem << std::endl;
    elem.reserve(nb_elem);
    for (size_t j = 0; j < nb_elem; j++)
    {
        size_t tmp;
        ifs >> tmp;
        elem.push_back(tmp);
    }

    return std::make_pair(cell_id, elem);
}

bool
read_cell_block(std::ifstream& ifs, const size_t nb_cells, std::vector<std::pair<size_t, std::vector<size_t>>>& cells)
{
    cells.clear();
    cells.reserve(nb_cells);

    for (size_t i = 0; i < nb_cells; i++)
        cells.push_back(read_cell_line(ifs));

    return true;
}

template<typename T>
bool
read_nodes_block(std::ifstream& ifs, const size_t nb_nodes, std::vector<std::pair<size_t, std::array<T, 3>>>& nodes)
{
    nodes.clear();
    nodes.reserve(nb_nodes);

    size_t node_id;
    std::array<T, 3> coor;

    for (size_t i = 0; i < nb_nodes; i++)
    {
        ifs >> node_id;
        for (size_t j = 0; j < 3; j++)
        {
            ifs >> coor[j];
        }
        nodes.push_back(std::make_pair(node_id, coor));
    }

    return true;
}

std::pair<size_t, std::vector<size_t>>
read_grp_line(std::ifstream& ifs)
{
    std::vector<size_t> elem;
    std::string         name;
    size_t              cell_id, nb_elem;

    ifs >> name;
    ifs >> cell_id;
    ifs >> nb_elem;
    // std::cout << cell_id << ", " << nb_elem << std::endl;
    elem.reserve(nb_elem);
    for (size_t j = 0; j < nb_elem; j++)
    {
        size_t tmp;
        ifs >> tmp;
        elem.push_back(tmp);
    }

    return std::make_pair(cell_id, elem);
}

bool
read_grp_block(std::ifstream& ifs, const size_t nb_grp, std::vector<std::pair<size_t, std::vector<size_t>>>& grps)
{
    grps.clear();
    grps.reserve(nb_grp);

    for (size_t i = 0; i < nb_grp; i++)
        grps.push_back(read_grp_line(ifs));

    return true;
}
}

template<typename T>
class poly_mesh_loader<T, 2> : public mesh_loader<generic_mesh<T, 2>>
{
    typedef generic_mesh<T, 2>               mesh_type;
    typedef typename mesh_type::point_type   point_type;
    typedef typename mesh_type::node_type    node_type;
    typedef typename mesh_type::edge_type    edge_type;
    typedef typename mesh_type::surface_type surface_type;

    std::vector<point_type>                   m_points;
    std::vector<node_type>                    m_nodes;
    std::vector<std::pair<size_t, edge_type>> m_edges;

    std::vector<std::pair<size_t, std::array<T, 3>>> vts;
    std::vector<std::pair<size_t, std::vector<size_t>>> vts_to_grp;
    std::vector<std::pair<size_t, std::vector<size_t>>> edges_to_vts;
    std::vector<std::pair<size_t, std::vector<size_t>>> edges_to_grp;
    std::vector<std::pair<size_t, std::vector<size_t>>> faces_to_edges;
    std::vector<std::pair<size_t, std::vector<size_t>>> faces_to_vts;
    std::vector<std::pair<size_t, std::vector<size_t>>> faces_to_grp;

    bool
    poly2d_read(const std::string& filename)
    {
        std::ifstream ifs(filename);
        std::string   keyword;
        size_t        dim, version, nb_elems;


        if (!ifs.is_open())
        {
            std::cout << "Error opening " << filename << std::endl;
            return false;
        }

        ifs >> keyword;
        if (keyword != "**BeginMesh")
        {
            std::cout << "Expected keyword \"**BeginMesh\"" << std::endl;
            return false;
        }

        while (keyword != "**EndMesh")
        {
            ifs >> keyword;
            if (keyword == "*Dimension")
            {
                ifs >> dim;
                if (this->verbose())
                    std::cout << "2D-Mesh" << std::endl;
                if (dim != 2)
                {
                    std::cout << "Expected dimension == 2 (here: " << dim << ")" << std::endl;
                    return false;
                }
            }
            else if (keyword == "*Version")
            {
                ifs >> version;
                if (this->verbose())
                    std::cout << "Mesh Version: " << version << std::endl;
            }
            else if (keyword == "*Nodes")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " nodes" << std::endl;

                poly::read_nodes_block(ifs, nb_elems, vts);
            }
            else if (keyword == "*Nodes->Groups")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " nodes->groups" << std::endl;

                poly::read_grp_block(ifs, nb_elems, vts_to_grp);
            }
            else if (keyword == "*Edges->Nodes")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " edges->nodes" << std::endl;

                poly::read_cell_block(ifs, nb_elems, edges_to_vts);
            }
            else if (keyword == "*Edges->Groups")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " edges->groups" << std::endl;

                poly::read_grp_block(ifs, nb_elems, edges_to_grp);
            }
            else if (keyword == "*Faces->Nodes")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " faces->nodes" << std::endl;

                poly::read_cell_block(ifs, nb_elems, faces_to_vts);
            }
            else if (keyword == "*Faces->Edges")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " faces->edges" << std::endl;

                poly::read_cell_block(ifs, nb_elems, faces_to_edges);
            }
            else if (keyword == "*Faces->Groups")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " faces->groups" << std::endl;

                poly::read_grp_block(ifs, nb_elems, faces_to_grp);
            }
            else if (keyword == "**EndMesh")
            {
                return true;
            }
            else
            {
                std::cout << "Unexpected keyword: " << keyword << std::endl;
                return false;
            }
        }

        return true;
    }

  public:
    static const char constexpr* expected_extension = "poly2d";
    poly_mesh_loader() = default;

    bool
    read_mesh(const std::string& s)
    {
        if (this->verbose())
            std::cout << " *** READING POLY2D MESH ***" << std::endl;
        return poly2d_read(s);
    }

    bool
    populate_mesh(mesh_type& msh)
    {
        if (this->verbose())
            std::cout << " *** POPULATE POLY2D MESH ***" << std::endl;

        auto storage = msh.backend_storage();

        // Create nodes and vertices
        m_points.reserve(vts.size());
        m_nodes.reserve(vts.size());
        for (auto& [id, coor] : vts)
        {
            m_points.push_back(point_type({coor[0], coor[1]}));
            m_nodes.push_back(node_type(disk::point_identifier<2>(id)));
        }

        // create edges
        m_edges.reserve(edges_to_vts.size());
        for (auto& [id, vts] : edges_to_vts)
        {
            size_t v1 = vts[0], v2 = vts[1];

            if (v1 > v2)
                std::swap(v1, v2);

            auto e = edge_type({typename node_type::id_type(v1), typename node_type::id_type(v2)});

            m_edges.push_back(std::make_pair(id, e));
        }

        std::vector<size_t> conv_table;
        /* Sort the edges in lexicographical order, remember their original
         * position to convert the pointers in the faces */
        auto comp_edges = [](const std::pair<size_t, edge_type>& e1, const std::pair<size_t, edge_type>& e2)
        {
            return e1.second < e2.second;
        };
        std::sort(m_edges.begin(), m_edges.end(), comp_edges);

        std::vector<edge_type> edges;
        edges.reserve(m_edges.size());
        conv_table.resize(m_edges.size());
        for (size_t i = 0; i < m_edges.size(); i++)
        {
            conv_table.at(m_edges.at(i).first) = i; /* Make ptr conversion table */
            edges.push_back(m_edges.at(i).second);
        }

        /* Convert the edge pointers in the face data */
        for (auto& fe : faces_to_edges)
        {
            for (auto& ptr : fe.second)
                ptr = conv_table[ptr];
        }

        /* Detect which ones are boundary edges */
        storage->boundary_info.resize(edges.size());
        for (auto& [id, edges] : edges_to_grp)
        {
            for (auto& edge_id : edges)
            {
                boundary_descriptor bi(id, true);
                storage->boundary_info.at(conv_table[edge_id]) = bi;
            }
        }

        /* Sort in lexicographical order and remember original position */
        auto comp_vecs = [](const std::pair<size_t, std::vector<size_t>>& e1,
                            const std::pair<size_t, std::vector<size_t>>& e2)
                            { return e1.second < e2.second; };
        std::sort(faces_to_edges.begin(), faces_to_edges.end(), comp_vecs);

        std::vector<surface_type> faces;
        faces.reserve(faces_to_edges.size());
        conv_table.resize(faces_to_edges.size());

        for (size_t i = 0; i < faces_to_edges.size(); i++)
        {
            auto         fe = faces_to_edges.at(i);
            surface_type s(convert_to<typename edge_type::id_type>(fe.second));
            s.set_point_ids(convert_to<disk::point_identifier<2>>(faces_to_vts.at(fe.first).second));
            faces.push_back(s);
            conv_table[fe.first] = i;
        }
        /* Now the faces are in their place and have correct ptrs */


        storage->points   = std::move(m_points);
        storage->nodes    = std::move(m_nodes);
        storage->edges    = std::move(edges);
        storage->surfaces = std::move(faces);


        return true;
    }
};

template<typename T>
class poly_mesh_loader<T, 3> : public mesh_loader<generic_mesh<T, 3>>
{
    typedef generic_mesh<T, 3>               mesh_type;
    typedef typename mesh_type::point_type   point_type;
    typedef typename mesh_type::node_type    node_type;
    typedef typename mesh_type::edge_type    edge_type;
    typedef typename mesh_type::surface_type surface_type;
    typedef typename mesh_type::volume_type  volume_type;

    std::vector<point_type>                   m_points;
    std::vector<node_type>                    m_nodes;
    std::vector<std::pair<size_t, edge_type>> m_edges;

    std::vector<std::pair<size_t, std::array<T, 3>>>    vts;
    std::vector<std::pair<size_t, std::vector<size_t>>> vts_to_grp;
    std::vector<std::pair<size_t, std::vector<size_t>>> edges_to_vts;
    std::vector<std::pair<size_t, std::vector<size_t>>> edges_to_grp;
    std::vector<std::pair<size_t, std::vector<size_t>>> faces_to_edges;
    std::vector<std::pair<size_t, std::vector<size_t>>> faces_to_vts;
    std::vector<std::pair<size_t, std::vector<size_t>>> faces_to_grp;
    std::vector<std::pair<size_t, std::vector<size_t>>> vols_to_edges;
    std::vector<std::pair<size_t, std::vector<size_t>>> vols_to_faces;
    std::vector<std::pair<size_t, std::vector<size_t>>> vols_to_vts;
    std::vector<std::pair<size_t, std::vector<size_t>>> vols_to_grp;

    bool
    poly3d_read(const std::string& filename)
    {
        std::ifstream ifs(filename);
        std::string   keyword;
        size_t        dim, version, nb_elems;

        if (!ifs.is_open())
        {
            std::cout << "Error opening " << filename << std::endl;
            return false;
        }

        ifs >> keyword;
        if (keyword != "**BeginMesh")
        {
            std::cout << "Expected keyword \"**BeginMesh\"" << std::endl;
            return false;
        }

        while (keyword != "**EndMesh")
        {
            ifs >> keyword;
            if (keyword == "*Dimension")
            {
                ifs >> dim;
                if (this->verbose())
                    std::cout << "3D-Mesh" << std::endl;
                if (dim != 3)
                {
                    std::cout << "Expected dimension == 3 (here: " << dim << ")" << std::endl;
                    return false;
                }
            }
            else if (keyword == "*Version")
            {
                ifs >> version;
                if (this->verbose())
                    std::cout << "Mesh Version: " << version << std::endl;
            }
            else if (keyword == "*Nodes")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " nodes" << std::endl;

                poly::read_nodes_block(ifs, nb_elems, vts);
            }
            else if (keyword == "*Nodes->Groups")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " nodes->groups" << std::endl;

                poly::read_grp_block(ifs, nb_elems, vts_to_grp);
            }
            else if (keyword == "*Edges->Nodes")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " edges->nodes" << std::endl;

                poly::read_cell_block(ifs, nb_elems, edges_to_vts);
            }
            else if (keyword == "*Edges->Groups")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " edges->groups" << std::endl;

                poly::read_grp_block(ifs, nb_elems, edges_to_grp);
            }
            else if (keyword == "*Faces->Nodes")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " faces->nodes" << std::endl;

                poly::read_cell_block(ifs, nb_elems, faces_to_vts);
            }
            else if (keyword == "*Faces->Edges")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " faces->edges" << std::endl;

                poly::read_cell_block(ifs, nb_elems, faces_to_edges);
            }
            else if (keyword == "*Faces->Groups")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " faces->groups" << std::endl;

                poly::read_grp_block(ifs, nb_elems, faces_to_grp);
            }
            else if (keyword == "*Volumes->Nodes")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " volumes->nodes" << std::endl;

                poly::read_cell_block(ifs, nb_elems, vols_to_vts);
            }
            else if (keyword == "*Volumes->Edges")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " volumes->edges" << std::endl;

                poly::read_cell_block(ifs, nb_elems, vols_to_edges);
            }
            else if (keyword == "*Volumes->Faces")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " volumes->faces" << std::endl;

                poly::read_cell_block(ifs, nb_elems, vols_to_faces);
            }
            else if (keyword == "*Volumes->Groups")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " volumes->groups" << std::endl;

                poly::read_grp_block(ifs, nb_elems, vols_to_grp);
            }
            else if (keyword == "**EndMesh")
            {
                return true;
            }
            else
            {
                std::cout << "Unexpected keyword: " << keyword << std::endl;
                return false;
            }
        }

        return true;
    }

  public:
    static const char constexpr* expected_extension = "poly3d";
    poly_mesh_loader()                              = default;

    bool
    read_mesh(const std::string& s)
    {
        if (this->verbose())
            std::cout << " *** READING POLY3D MESH ***" << std::endl;
        return poly3d_read(s);
    }

    bool
    populate_mesh(mesh_type& msh)
    {
        if (this->verbose())
            std::cout << " *** POPULATE POLY3D MESH ***" << std::endl;

        auto storage = msh.backend_storage();

        // Create nodes and vertices
        m_points.reserve(vts.size());
        m_nodes.reserve(vts.size());
        for (auto& [id, coor] : vts)
        {
            m_points.push_back(point_type({coor[0], coor[1], coor[2]}));
            m_nodes.push_back(node_type(disk::point_identifier<3>(id)));
        }

        // create edges
        m_edges.reserve(edges_to_vts.size());
        for (auto& [id, vts] : edges_to_vts)
        {
            size_t v1 = vts[0], v2 = vts[1];

            if (v1 > v2)
                std::swap(v1, v2);

            auto e = edge_type({typename node_type::id_type(v1), typename node_type::id_type(v2)});

            m_edges.push_back(std::make_pair(id, e));
        }

        std::vector<size_t> conv_table;
        /* Sort the edges in lexicographical order, remember their original
         * position to convert the pointers in the faces */
        auto comp_edges = [](const std::pair<size_t, edge_type>& e1, const std::pair<size_t, edge_type>& e2)
        {
            return e1.second < e2.second;
        };
        std::sort(m_edges.begin(), m_edges.end(), comp_edges);

        std::vector<edge_type> edges;
        edges.reserve(m_edges.size());
        conv_table.resize(m_edges.size());
        for (size_t i = 0; i < m_edges.size(); i++)
        {
            conv_table.at(m_edges.at(i).first) = i; /* Make ptr conversion table */
            edges.push_back(m_edges.at(i).second);
        }

        /* Convert the edge pointers in the face data */
        for (auto& fe : faces_to_edges)
        {
            for (auto& ptr : fe.second)
                ptr = conv_table[ptr];
        }

        /* Sort in lexicographical order and remember original position */
        auto comp_vecs = [](const std::pair<size_t, std::vector<size_t>>& e1,
                            const std::pair<size_t, std::vector<size_t>>& e2)
                            { return e1.second < e2.second; };
        std::sort(faces_to_edges.begin(), faces_to_edges.end(), comp_vecs);

        std::vector<surface_type> faces;
        faces.reserve(faces_to_edges.size());
        conv_table.resize(faces_to_edges.size());

        for (size_t i = 0; i < faces_to_edges.size(); i++)
        {
            auto         fe = faces_to_edges.at(i);
            surface_type s(convert_to<typename edge_type::id_type>(fe.second));
            s.set_point_ids(convert_to<disk::point_identifier<3>>(faces_to_vts.at(fe.first).second));
            faces.push_back(s);
            conv_table[fe.first] = i;
        }
        /* Now the faces are in their place and have correct ptrs */

        /* Detect which ones are boundary edges */
        storage->boundary_info.resize(faces.size());
        for (auto& [id, faces] : faces_to_grp)
        {
            for (auto& face_id : faces)
            {
                boundary_descriptor bi(id, true);
                storage->boundary_info.at(conv_table.at(face_id)) = bi;
            }
        }

        /* Convert the face pointers in the volume data */
        for (auto& vf : vols_to_faces)
        {
            for (auto& ptr : vf.second)
                ptr = conv_table[ptr];
        }

        /* Sort volume data */
        std::sort(vols_to_faces.begin(), vols_to_faces.end(), comp_vecs);

        std::vector<volume_type> volumes;
        volumes.reserve(vols_to_faces.size());

        for (size_t i = 0; i < vols_to_faces.size(); i++)
        {
            auto        vf = vols_to_faces[i];
            volume_type v(convert_to<typename surface_type::id_type>(vf.second));
            v.set_point_ids(convert_to<disk::point_identifier<3>>(vols_to_vts.at(vf.first).second));
            volumes.push_back(v);
        }

        storage->points   = std::move(m_points);
        storage->nodes    = std::move(m_nodes);
        storage->edges    = std::move(edges);
        storage->surfaces = std::move(faces);
        storage->volumes  = std::move(volumes);

        return true;
    }
};

#if 0
/* Helper to load uniform 1D meshes. */
template<typename T>
[[deprecated("DiSk++ deprecation: The load_mesh_*() functions should be preferred")]]
disk::generic_mesh<T,1>
load_uniform_1d_mesh(T min, T max, size_t cells)
{
    typedef disk::generic_mesh<T, 1>  mesh_type;

    mesh_type msh;
    disk::uniform_mesh_loader<T, 1> loader(min, max, cells);
    loader.populate_mesh(msh);

    return msh;
}


/* Helper to load 2D meshes in FVCA5 format */
template<typename T>
[[deprecated("DiSk++ deprecation: The load_mesh_*() functions should be preferred")]]
disk::generic_mesh<T,2>
load_fvca5_2d_mesh(const char *filename, bool verbose = false)
{
    typedef disk::generic_mesh<T, 2>  mesh_type;

    mesh_type msh;
    disk::fvca5_mesh_loader<T, 2> loader;
    loader.verbose(verbose);
    loader.read_mesh(filename);
    loader.populate_mesh(msh);

    return msh;
}

/* Helper to load 3D meshes in FVCA6 format */
template<typename T>
[[deprecated("DiSk++ deprecation: The load_mesh_*() functions should be preferred")]]
disk::generic_mesh<T, 3>
load_fvca6_3d_mesh(const char* filename)
{
    typedef disk::generic_mesh<T, 3> mesh_type;

    mesh_type                     msh;
    disk::fvca6_mesh_loader<T, 3> loader;
    loader.read_mesh(filename);
    loader.populate_mesh(msh);

    return msh;
}

/* Helper to load 2D meshes in Netgen format */
template<typename T>
[[deprecated("DiSk++ deprecation: The load_mesh_*() functions should be preferred")]]
disk::simplicial_mesh<T, 2>
load_netgen_2d_mesh(const char *filename)
{
    typedef disk::simplicial_mesh<T, 2>  mesh_type;

    mesh_type msh;
    disk::netgen_mesh_loader<T, 2> loader;
    loader.read_mesh(filename);
    loader.populate_mesh(msh);

    return msh;
}

/* Helper to load 2D meshes in DiSk++ format */
template<typename T>
[[deprecated("DiSk++ deprecation: The load_mesh_*() functions should be preferred")]]
disk::cartesian_mesh<T, 2>
load_cartesian_2d_mesh(const char *filename)
{
    typedef disk::cartesian_mesh<T, 2>  mesh_type;

    mesh_type msh;
    disk::cartesian_mesh_loader<T, 2> loader;
    loader.read_mesh(filename);
    loader.populate_mesh(msh);

    return msh;
}

/* Helper to load 3D meshes in Netgen format */
template<typename T>
[[deprecated("DiSk++ deprecation: The load_mesh_*() functions should be preferred")]]
disk::simplicial_mesh<T, 3>
load_netgen_3d_mesh(const char *filename)
{
    typedef disk::simplicial_mesh<T, 3>  mesh_type;

    mesh_type msh;
    disk::netgen_mesh_loader<T, 3> loader;
    loader.read_mesh(filename);
    loader.populate_mesh(msh);

    return msh;
}

/* Helper to load 3D meshes in DiSk++ format */
template<typename T>
[[deprecated("DiSk++ deprecation: The load_mesh_*() functions should be preferred")]]
disk::cartesian_mesh<T, 3>
load_cartesian_3d_mesh(const char *filename)
{
    typedef disk::cartesian_mesh<T, 3>  mesh_type;

    mesh_type msh;
    disk::cartesian_mesh_loader<T, 3> loader;
    loader.read_mesh(filename);
    loader.populate_mesh(msh);

    return msh;
}


/* Helper to load 2D meshes in Medit format */
template<typename T>
[[deprecated("DiSk++ deprecation: The load_mesh_*() functions should be preferred")]]
disk::generic_mesh<T, 2>
load_medit_2d_mesh(const char* filename)
{
   typedef disk::generic_mesh<T, 2> mesh_type;

   mesh_type                     msh;
   disk::medit_mesh_loader<T, 2> loader;
   loader.read_mesh(filename);
   loader.populate_mesh(msh);

   return msh;
}

/* Helper to load 3D meshes in Medit format */
template<typename T>
[[deprecated("DiSk++ deprecation: The load_mesh_*() functions should be preferred")]]
disk::generic_mesh<T, 3>
load_medit_3d_mesh(const char* filename)
{
   typedef disk::generic_mesh<T, 3> mesh_type;

   mesh_type                     msh;
   disk::medit_mesh_loader<T, 3> loader;
   loader.read_mesh(filename);
   loader.populate_mesh(msh);

   return msh;
}
#endif
/**************************************************************************/
// New mesh loader helpers
namespace priv {

bool
check_filename_extension(const char *filename, const char *extension)
{
    std::stringstream ss;
    ss << ".*\\." << extension << "$";

    return std::regex_match(filename, std::regex(ss.str()));
}

template<typename LT, typename MT>
bool
load_mesh(const char *filename, LT& loader, MT& msh)
{
    if ( !check_filename_extension(filename, LT::expected_extension) )
    {
        std::cout << "Warning: unexpected filename extension for ";
        std::cout << "the required mesh type" << std::endl;
    }

    bool success = loader.read_mesh(filename);
    if (!success)
        return false;

    loader.populate_mesh(msh);
    return true;
}

} //namespace priv

template<typename T>
bool
load_mesh_netgen(const char *filename, disk::simplicial_mesh<T, 2>& msh)
{
    disk::netgen_mesh_loader<T, 2> loader;
    return priv::load_mesh(filename, loader, msh);
}

template<typename T>
bool
load_mesh_netgen(const char *filename, disk::simplicial_mesh<T, 3>& msh)
{
    disk::netgen_mesh_loader<T, 3> loader;
    return priv::load_mesh(filename, loader, msh);
}

template<typename T>
bool
load_mesh_diskpp_cartesian(const char *filename, disk::cartesian_mesh<T, 2>& msh)
{
    disk::cartesian_mesh_loader<T, 2> loader;
    return priv::load_mesh(filename, loader, msh);
}

template<typename T>
bool
load_mesh_diskpp_cartesian(const char *filename, disk::cartesian_mesh<T, 3>& msh)
{
    disk::cartesian_mesh_loader<T, 3> loader;
    return priv::load_mesh(filename, loader, msh);
}

template<typename T>
bool
load_mesh_fvca5_2d(const char* filename, disk::generic_mesh<T, 2>& msh)
{
    disk::fvca5_mesh_loader<T, 2> loader;
    return priv::load_mesh(filename, loader, msh);
}

template<typename T>
bool
load_mesh_fvca6_3d(const char *filename, disk::generic_mesh<T, 3>& msh)
{
    disk::fvca6_mesh_loader<T, 3> loader;
    return priv::load_mesh(filename, loader, msh);
}

template<typename T>
bool
load_mesh_poly2d(const char* filename, disk::generic_mesh<T, 2>& msh, bool verbose = false)
{
    disk::poly_mesh_loader<T, 2> loader;
    loader.verbose(verbose);
    return priv::load_mesh(filename, loader, msh);
}

template<typename T>
bool
load_mesh_poly3d(const char* filename, disk::generic_mesh<T, 3>& msh, bool verbose = false)
{
    disk::poly_mesh_loader<T, 3> loader;
    loader.verbose(verbose);
    return priv::load_mesh(filename, loader, msh);
}


template<typename T>
bool
load_mesh_medit(const char *filename, disk::generic_mesh<T, 2>& msh)
{
    disk::medit_mesh_loader<T, 2> loader;
    return priv::load_mesh(filename, loader, msh);
}

template<typename T>
bool
load_mesh_medit(const char *filename, disk::generic_mesh<T, 3>& msh)
{
    disk::medit_mesh_loader<T, 3> loader;
    return priv::load_mesh(filename, loader, msh);
}

} // namespace disk



#ifdef HAVE_GMSH
#include "loader_gmsh.hpp"
#endif


namespace disk {

/**
 * @brief Deduce mesh type from filename extension, create appropriate mesh object
 *        and dispatch a function on it
 *
 * If you have the function `process_mesh(msh, p1, p2, ..., pn)` you can
 *
 *  ```
 *  disk::dispatch_all_meshes(mesh_filename,
 *        [](auto ...args) { process_mesh(args...); },
 *        p1, p2, ..., pn);
 *  ```
 *
 * @tparam Function
 * @tparam Args
 * @param mesh_filename Name of the file containing the mesh.
 * @param func Function to dispatch. The mesh object is passed as first parameter.
 * @param args Arguments to the function to dispatch.
 * @return int
 */

template<typename Function, typename... Args>
int
dispatch_all_meshes(const char *mesh_filename, Function func, Args && ...args)
{
    using T = double;
    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        disk::generic_mesh<T, 2> msh;
        if ( load_mesh_fvca5_2d(mesh_filename, msh) ) {
            func(msh, args...);
            return 0;
        }
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        disk::simplicial_mesh<T, 2> msh;
        if ( load_mesh_netgen(mesh_filename, msh) ) {
            func(msh, args...);
            return 0;
        }
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$") ))
    {
        disk::cartesian_mesh<T, 2> msh;
        if ( load_mesh_diskpp_cartesian(mesh_filename, msh) ) {
            func(msh, args...);
            return 0;
        }
    }


    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        disk::simplicial_mesh<T, 3> msh;
        if ( load_mesh_netgen(mesh_filename, msh) ) {
            func(msh, args...);
            return 0;
        }
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        disk::cartesian_mesh<T, 3> msh;
        if ( load_mesh_diskpp_cartesian(mesh_filename, msh) ) {
            func(msh, args...);
            return 0;
        }
    }

    /* FVCA6 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
        disk::generic_mesh<T, 3> msh;
        if ( load_mesh_fvca6_3d(mesh_filename, msh) ) {
            func(msh, args...);
            return 0;
        }
    }

#ifdef HAVE_GMSH
    /* GMSH 2D simplicials */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo2s$") ))
    {
        std::cout << "Guessed mesh format: GMSH 2D simplicials" << std::endl;
        disk::simplicial_mesh<T,2> msh;
        disk::gmsh_geometry_loader< disk::simplicial_mesh<T,2> > loader;
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);
        func(msh, args...);
        return 0;
    }

    /* GMSH 3D simplicials */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3s$") ))
    {
        std::cout << "Guessed mesh format: GMSH 3D simplicials" << std::endl;
        disk::simplicial_mesh<T,3> msh;
        disk::gmsh_geometry_loader< disk::simplicial_mesh<T,3> > loader;
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);
        func(msh, args...);
        return 0;
    }

    /* GMSH 3D generic */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3g$") ))
    {
        std::cout << "Guessed mesh format: GMSH 3D general elements" << std::endl;
        disk::generic_mesh<T,3> msh;
        disk::gmsh_geometry_loader< disk::generic_mesh<T,3> > loader;
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);
        func(msh, args...);
        return 0;
    }
#endif

    return 1;
}

} // namespace disk
