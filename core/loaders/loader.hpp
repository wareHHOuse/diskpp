/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
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

#include <vector>
#include <array>
#include <fstream>
#include <cassert>
#include <thread>

#include "geometry/geometry_generic.hpp"
#include "geometry/geometry_simplicial.hpp"

#include "mapped_file.h"
#include "strtot.hpp"

namespace disk {

template<typename mesh_type>
class mesh_loader
{
public:
    mesh_loader() {}

    virtual bool    read_mesh(const std::string&) { return false; }
    virtual bool    populate_mesh(mesh_type&)    = 0;

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
            auto e = edge_type{{n0, n1}};

            std::vector<point_identifier<1>> pts(2);
            pts[0] = point_identifier<1>(i);
            pts[1] = point_identifier<1>(i+1);
            e.set_point_ids(pts.begin(), pts.end());

            storage->edges.at(i) = e;
        }

        for (size_t i = 0; i < num_nodes; i++)
            storage->nodes.at(i) = node_type(point_identifier<1>(i));

        storage->boundary_nodes.resize(num_nodes);
        storage->boundary_nodes.at(0) = true;
        storage->boundary_nodes.at(num_nodes - 1) = true;

        return true;
    }
};

template<typename T, size_t N>
class fvca5_mesh_loader
{
    static_assert(N == 2, "fvca5 is a 2D-only mesh format");
};

template<typename T>
class fvca5_mesh_loader<T,2> : public mesh_loader<generic_mesh<T, 2>>
{
    //const char* keywords[] = {"triangles", "quadrangles", "pentagons", "hexagons",
    //                     "ennagons", "ettagons", "nonecagons" };

    typedef generic_mesh<T,2>                       mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;
    typedef typename mesh_type::surface_type        surface_type;

    template<size_t N> using polygon = std::array<ident_impl_t, N>;

    std::vector<point_type>                         m_points;
    std::vector<polygon<3>>                         m_triangles;
    std::vector<polygon<4>>                         m_quadrangles;
    std::vector<polygon<5>>                         m_pentagons;
    std::vector<polygon<6>>                         m_hexagons;
    std::vector<polygon<7>>                         m_ennagons;

    std::vector<std::array<ident_impl_t, 2>>        m_boundary_edges;
    std::vector<std::array<ident_impl_t, 4>>        m_edges;

    //std::vector<std::set<typename edge_type::id_type>>      m_edge_node_connectivity;

    bool fvca5_read_points(std::ifstream& ifs)
    {
        size_t      elements_to_read;
        T           x, y;

        ifs >> elements_to_read;
        std::cout << "Attempting to read " << elements_to_read << " points" << std::endl;
        m_points.reserve(elements_to_read);

        for (size_t i = 0; i < elements_to_read; i++)
        {
            ifs >> x >> y;
            m_points.push_back(point_type{x,y});
        }

        return true;
    }

    template<typename U, size_t N>
    bool fvca5_read_tuples(std::ifstream& ifs,
                           std::vector<std::array<U, N>>& tuples)
    {
        size_t      elements_to_read;

        ifs >> elements_to_read;
        std::cout << "Attempting to read " << elements_to_read << " " << N << "-ples " << std::endl;

        tuples.reserve(elements_to_read);

        for (size_t i = 0; i < elements_to_read; i++)
        {
            std::array<U, N> tuple;
            for (size_t j = 0; j < N; j++)
            {
                U val;
                ifs >> val;
                tuple[j] = val-1;
            }

            tuples.push_back(tuple);
        }

        return true;
    }

    bool fvca5_read(const std::string& filename)
    {
        std::ifstream   ifs(filename);
        std::string     keyword;

        if (!ifs.is_open())
        {
            std::cout << "Error opening " << filename << std::endl;
            return false;
        }

        ifs >> keyword;
        if ( keyword != "vertices" )
        {
            std::cout << "Expected keyword \"vertices\"" << std::endl;
            return false;
        }

        fvca5_read_points(ifs);

        ifs >> keyword;
        if ( keyword == "triangles" )
        {
            m_triangles.clear();
            fvca5_read_tuples(ifs, m_triangles);
            ifs >> keyword;
        }

        if ( keyword == "quadrangles" )
        {
            m_quadrangles.clear();
            fvca5_read_tuples(ifs, m_quadrangles);
            ifs >> keyword;
        }

        if ( keyword == "pentagons" )
        {
            m_pentagons.clear();
            fvca5_read_tuples(ifs, m_pentagons);
            ifs >> keyword;
        }

        if ( keyword == "hexagons" )
        {
            m_hexagons.clear();
            fvca5_read_tuples(ifs, m_hexagons);
            ifs >> keyword;
        }

        if ( keyword == "ennagons" )
        {
            m_ennagons.clear();
            fvca5_read_tuples(ifs, m_ennagons);
            ifs >> keyword;
        }

        if ( keyword == "edges" )
        {
            std::getline(ifs, keyword); //drop the rest of the line
            m_boundary_edges.clear();
            fvca5_read_tuples(ifs, m_boundary_edges);
            ifs >> keyword;
        }
        else
        {
            std::cout << "Error parsing FVCA5 file" << std::endl;
            return false;
        }

        if ( keyword == "all" )
        {
            std::getline(ifs, keyword); //drop the rest of the line
            m_edges.clear();
            fvca5_read_tuples(ifs, m_edges);
        }
        else
        {
            std::cout << "Error parsing FVCA5 file" << std::endl;
            return false;
        }

        ifs.close();
        return true;
    }

    template<size_t N>
    void put_polygons(mesh_type& msh,
                      const std::vector<polygon<N>>& polys,
                      std::vector<surface_type>& surfedg)
    {
        auto storage = msh.backend_storage();

        for (auto& p : polys)
        {
            std::vector<typename edge_type::id_type> surface_edges(N);
            for (size_t i = 0; i < N; i++)
            {
                auto pt1 = typename node_type::id_type(p[i]);
                auto pt2 = typename node_type::id_type(p[(i+1) % N]);

                if (pt2 < pt1)
                    std::swap(pt1, pt2);

                edge_type edge{{pt1, pt2}};
                auto edge_id = find_element_id(storage->edges.begin(),
                                               storage->edges.end(), edge);
                if (!edge_id.first)
                    throw std::invalid_argument("Edge not found (hanging nodes?)");

                surface_edges[i] = edge_id.second;
            }
            auto surface = surface_type(surface_edges);
            surface.set_point_ids(p.begin(), p.end()); /* XXX: crap */
            surfedg.push_back( surface );
        }

    }

public:
    fvca5_mesh_loader() = default;

    bool read_mesh(const std::string& filename)
    {
        std::cout << " *** READING FVCA5 MESH ***" << std::endl;
        return fvca5_read(filename);
    }

    bool populate_mesh(mesh_type& msh)
    {
        std::cout << " *** POPULATING FVCA5 MESH ***" << std::endl;
        auto storage = msh.backend_storage();

        /* Points */
        size_t nodes_size = m_points.size();
        storage->points = std::move(m_points);

        /* Nodes */
        std::vector<node_type> nodes(nodes_size);
        for (size_t i = 0; i < nodes_size; i++)
            nodes[i] = node_type(point_identifier<2>(i));

        storage->nodes = std::move(nodes);

        /* Edges */
        /* Make the vector containing the edges */
        std::vector<edge_type> edges;
        edges.reserve(m_edges.size());
        for (size_t i = 0; i < m_edges.size(); i++)
        {
            auto node1 = typename node_type::id_type(m_edges[i][0]);
            auto node2 = typename node_type::id_type(m_edges[i][1]);

            if (node2 < node1)
                std::swap(node1, node2);

            auto e = edge_type{{node1, node2}};

            e.set_point_ids(m_edges[i].begin(), m_edges[i].begin()+2); /* XXX: crap */
            edges.push_back(e);
        }
        /* Sort them */
        std::sort(edges.begin(), edges.end());

        /* Detect which ones are boundary edges */
        storage->boundary_edges.resize(m_edges.size());
        for (size_t i = 0; i < m_boundary_edges.size(); i++)
        {
            auto node1 = typename node_type::id_type(m_boundary_edges[i][0]);
            auto node2 = typename node_type::id_type(m_boundary_edges[i][1]);

            if (node2 < node1)
                std::swap(node1, node2);

            auto e = edge_type{{node1, node2}};

            auto position = find_element_id(edges.begin(), edges.end(), e);

            if (position.first == false)
            {
                std::cout << "Bad bug at " << __FILE__ << "("
                          << __LINE__ << ")" << std::endl;
                return false;
            }

            storage->boundary_edges.at(position.second) = true;
        }

        storage->edges = std::move(edges);

        /* Surfaces */
        std::vector<surface_type> surfaces;
        surfaces.reserve( m_triangles.size() +
                          m_quadrangles.size() +
                          m_pentagons.size() +
                          m_hexagons.size() +
                          m_ennagons.size() );

        put_polygons(msh, m_triangles, surfaces);
        m_triangles.clear();

        put_polygons(msh, m_quadrangles, surfaces);
        m_quadrangles.clear();

        put_polygons(msh, m_pentagons, surfaces);
        m_pentagons.clear();

        put_polygons(msh, m_hexagons, surfaces);
        m_hexagons.clear();

        put_polygons(msh, m_ennagons, surfaces);
        m_ennagons.clear();

        std::sort(surfaces.begin(), surfaces.end());
        storage->surfaces = std::move(surfaces);

        /* Print stats */
        storage->statistics();

        return true;
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
    T t1, t2, t3;

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

template<typename T>
void
sort_uniq(std::vector<T>& v)
{
    std::sort(v.begin(), v.end());
    auto uniq_iter = std::unique(v.begin(), v.end());
    v.erase(uniq_iter, v.end());
}

} //namespace priv

#define THREADED
#ifdef THREADED
    #define THREAD(name, body) std::thread name([&]{body})
    #define WAIT_THREAD(name) name.join()
#else
    #define THREAD(name, body) {body}
    #define WAIT_THREAD(name)
#endif

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
    std::vector<edge_type>                          edges, boundary_edges;
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
            if ( (linecount%100000) == 0 )
            {
                std::cout << "Reading points: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto point = priv::read_2d_point_line<T>(endptr, &endptr, 1.0);

            points.push_back( point );

            auto point_id = point_identifier<2>( linecount );
            auto node = node_type( { point_id } );

            nodes.push_back(node);

            linecount++;
        }

        std::cout << "Reading points: " << linecount;
        std::cout << "/" << lines << std::endl;

        /************************ Read triangles ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        edges.reserve(lines*3);
        surfaces.reserve(lines);

        while (linecount < lines)
        {
            if ( (linecount%100000) == 0 )
            {
                std::cout << "Reading triangles: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_triangle_line<size_t>(endptr, &endptr);

            point_identifier<2>     p0(std::get<1>(t));
            point_identifier<2>     p1(std::get<2>(t));
            point_identifier<2>     p2(std::get<3>(t));
            //domain_id_type      d(std::get<0>(t));

            edges.push_back( edge_type( { p0, p1 } ) );
            edges.push_back( edge_type( { p1, p2 } ) );
            edges.push_back( edge_type( { p0, p2 } ) );

            surfaces.push_back( surface_type( { p0, p1, p2 } ) );

            linecount++;
        }

        std::cout << "Reading triangles: " << linecount;
        std::cout << "/" << lines << std::endl;

        /************************ Read boundary edges ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        boundary_edges.reserve(lines);

        while (linecount < lines)
        {
            if ( (linecount%50000) == 0 )
            {
                std::cout << "Reading edges: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_edge_line<size_t>(endptr, &endptr);

            point_identifier<2>     p0(std::get<1>(t));
            point_identifier<2>     p1(std::get<2>(t));

            edge_type   edge( { p0, p1 } );

            boundary_edges.push_back( edge );

            linecount++;
        }

        std::cout << "Reading edges: " << linecount;
        std::cout << "/" << lines << std::endl;

        return true;
    }

public:
    netgen_mesh_loader() = default;

    bool read_mesh(const std::string& s)
    {
        std::cout << " *** READING NETGEN 2D MESH ***" << std::endl;
        return netgen_read(s);
    }

    bool populate_mesh(mesh_type& msh)
    {
        auto storage = msh.backend_storage();

        std::cout << "Sorting data...";
        std::cout.flush();

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

        storage->boundary_edges.resize(storage->edges.size());
        for (auto& be : boundary_edges)
        {
            auto position = find_element_id(storage->edges.begin(),
                                            storage->edges.end(), be);
            if (position.first == false)
            {
                std::cout << "Bad bug at " << __FILE__ << "("
                          << __LINE__ << ")" << std::endl;
                return false;
            }

            storage->boundary_edges[position.second] = true;
        }

        std::cout << "done." << std::endl;

        std::cout << "Nodes: " << storage->nodes.size() << std::endl;
        std::cout << "Edges: " << storage->edges.size() << std::endl;
        std::cout << "Faces: " << storage->surfaces.size() << std::endl;

        boundary_edges.clear();

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

    std::vector<point_type>                         points;
    std::vector<node_type>                          nodes;
    std::vector<edge_type>                          edges;
    std::vector<surface_type>                       surfaces, boundary_surfaces;
    std::vector<volume_type>                        volumes;


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
            if ( (linecount%100000) == 0 )
            {
                std::cout << "Reading points: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto point = priv::read_3d_point_line<T>(endptr, &endptr, 1.0);

            /*auto point = point_type( { std::get<0>(t),
                                       std::get<1>(t),
                                       std::get<2>(t) } );
*/
            points.push_back( point );

            auto point_id = point_identifier<3>( linecount );
            auto node = node_type( { point_id } );

            nodes.push_back(node);
            /* Do something with that point */

            linecount++;
        }

        std::cout << "Reading points: " << linecount;
        std::cout << "/" << lines << std::endl;

        /************************ Read tetrahedra ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        edges.reserve(lines*6);
        surfaces.reserve(lines*4);
        volumes.reserve(lines);

        while (linecount < lines)
        {
            if ( (linecount%100000) == 0 )
            {
                std::cout << "Reading tetrahedra: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_tetrahedron_line<size_t>(endptr, &endptr);

            point_identifier<3>     p0(std::get<1>(t));
            point_identifier<3>     p1(std::get<2>(t));
            point_identifier<3>     p2(std::get<3>(t));
            point_identifier<3>     p3(std::get<4>(t));
            //domain_id_type      d(std::get<0>(t));

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

            //auto tuple = std::make_tuple(volume_type(p0, p1, p2, p3), d);
            //temp_tet.push_back( tuple );

            volumes.push_back( volume_type( { p0, p1, p2, p3 } ) );

            linecount++;
        }

        std::cout << "Reading tetrahedra: " << linecount;
        std::cout << "/" << lines << std::endl;

        /************************ Read boundary surfaces ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        boundary_surfaces.reserve(lines);

        while (linecount < lines)
        {
            if ( (linecount%50000) == 0 )
            {
                std::cout << "Reading triangle: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_triangle_line<size_t>(endptr, &endptr);

            point_identifier<3>     p0(std::get<1>(t));
            point_identifier<3>     p1(std::get<2>(t));
            point_identifier<3>     p2(std::get<3>(t));

            surface_type   tri( { p0, p1, p2 } );

            boundary_surfaces.push_back( tri );

            linecount++;
        }

        std::cout << "Reading triangle: " << linecount;
        std::cout << "/" << lines << std::endl;

        return true;
    }

public:
    netgen_mesh_loader() = default;

    bool read_mesh(const std::string& s)
    {
        std::cout << " *** READING NETGEN 3D MESH ***" << std::endl;
        return netgen_read(s);
    }

    bool populate_mesh(mesh_type& msh)
    {
        auto storage = msh.backend_storage();

        std::cout << "Sorting data...";
        std::cout.flush();

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

        /* sort tetrahedra, make unique and move them in geometry */
        THREAD(tet_thread,
            std::sort(volumes.begin(), volumes.end());
            storage->volumes = std::move(volumes);
        );

        /* wait for the threads */
        WAIT_THREAD(edge_thread);
        WAIT_THREAD(tri_thread);
        WAIT_THREAD(tet_thread);

        storage->boundary_surfaces.resize(storage->surfaces.size());
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

            storage->boundary_surfaces[position.second] = true;
        }

        std::cout << "done." << std::endl;

        std::cout << "Nodes: " << storage->nodes.size() << std::endl;
        std::cout << "Edges: " << storage->edges.size() << std::endl;
        std::cout << "Faces: " << storage->surfaces.size() << std::endl;
        std::cout << "Volumes: " << storage->nodes.size() << std::endl;

        boundary_surfaces.clear();

        return true;
    }

};






} // namespace disk
