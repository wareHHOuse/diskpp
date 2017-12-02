/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016-2017 - matteo.cicuttin@enpc.fr
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

#ifndef _LOADER_HPP_WAS_INCLUDED_
#error "You must NOT include this file directly. Include loader.hpp"
#endif

#ifndef _LOADER_FVCA_HPP_
#define _LOADER_FVCA_HPP_


#include <vector>
#include <array>
#include <fstream>
#include <cassert>
#include <thread>
#include <set>

#include "geometry/geometry.hpp"
#include "loader_utils.hpp"
#include "loader_template.hpp"

#include "mapped_file.h"
#include "strtot.hpp"

namespace disk {

   template<typename T, size_t N>
   class fvca5_mesh_loader
   {
      static_assert(N == 2, "fvca5 is a 2D-only mesh format");
   };

   template<typename T>
   class fvca5_mesh_loader<T,2> : public mesh_loader<generic_mesh<T, 2>>
   {
      typedef generic_mesh<T,2>                       mesh_type;
      typedef typename mesh_type::point_type          point_type;
      typedef typename mesh_type::node_type           node_type;
      typedef typename mesh_type::edge_type           edge_type;
      typedef typename mesh_type::surface_type        surface_type;

      struct fvca5_poly
      {
         std::vector<size_t>                 nodes;
         std::set<std::array<size_t, 2>>     attached_edges;

         bool operator<(const fvca5_poly& other) {
            return nodes < other.nodes;
         }
      };

      std::vector<point_type>                         m_points;
      std::vector<fvca5_poly>                         m_polys;
      std::vector<std::array<ident_impl_t, 2>>        m_boundary_edges;
      std::vector<std::array<ident_impl_t, 4>>        m_edges;

      bool fvca5_read_points(std::ifstream& ifs)
      {
         size_t      elements_to_read;
         T           x, y;

         ifs >> elements_to_read;

         if (this->verbose())
            std::cout << "Attempting to read " << elements_to_read << " points" << std::endl;

         m_points.reserve(elements_to_read);

         for (size_t i = 0; i < elements_to_read; i++)
         {
            ifs >> x >> y;
            m_points.push_back(point_type{x,y});
         }

         return true;
      }

      bool fvca5_read_polygons(std::ifstream& ifs, size_t polynum)
      {
         size_t      elements_to_read;

         ifs >> elements_to_read;
         if (this->verbose())
            std::cout << "Reading " << elements_to_read << " " << polynum << "-angles" << std::endl;

         for (size_t i = 0; i < elements_to_read; i++)
         {
            fvca5_poly p;

            for (size_t j = 0; j < polynum; j++)
            {
               ident_impl_t val;
               ifs >> val;
               p.nodes.push_back(val-1);
            }

            m_polys.push_back(p);
         }

         return true;
      }

      bool fvca5_read_boundary_edges(std::ifstream& ifs)
      {
         size_t      elements_to_read;

         ifs >> elements_to_read;
         if (this->verbose())
            std::cout << "Reading " << elements_to_read << " boundary edges" << std::endl;

         m_boundary_edges.reserve(elements_to_read);

         for (size_t i = 0; i < elements_to_read; i++)
         {
            std::array<ident_impl_t, 2> b_edge;
            ifs >> b_edge[0]; b_edge[0] -= 1;
            ifs >> b_edge[1]; b_edge[1] -= 1;

            assert(b_edge[0] != b_edge[1]);

            if (b_edge[0] > b_edge[1])
               std::swap(b_edge[0], b_edge[1]);

            m_boundary_edges.push_back({b_edge[0], b_edge[1]});
         }

         return true;
      }

      bool fvca5_read_edges(std::ifstream& ifs)
      {
         size_t      elements_to_read;

         ifs >> elements_to_read;
         if (this->verbose())
            std::cout << "Reading " << elements_to_read << " edges" << std::endl;

         m_edges.reserve(elements_to_read);

         for (size_t i = 0; i < elements_to_read; i++)
         {
            std::array<ident_impl_t, 4> edge;
            ifs >> edge[0]; edge[0] -= 1;
            ifs >> edge[1]; edge[1] -= 1;
            ifs >> edge[2];
            ifs >> edge[3];

            assert(edge[0] != edge[1]);

            if (edge[0] > edge[1])
               std::swap(edge[0], edge[1]);

            if (edge[2] != 0)
               m_polys.at(edge[2]-1).attached_edges.insert({edge[0], edge[1]});

            if (edge[3] != 0)
               m_polys.at(edge[3]-1).attached_edges.insert({edge[0], edge[1]});

            m_edges.push_back(edge);
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
            fvca5_read_polygons(ifs, 3);
            ifs >> keyword;
         }

         if ( keyword == "quadrangles" )
         {
            fvca5_read_polygons(ifs, 4);
            ifs >> keyword;
         }

         if ( keyword == "pentagons" )
         {
            fvca5_read_polygons(ifs, 5);
            ifs >> keyword;
         }

         if ( keyword == "hexagons" )
         {
            fvca5_read_polygons(ifs, 6);
            ifs >> keyword;
         }

         if ( keyword == "ennagons" )
         {
            fvca5_read_polygons(ifs, 7);
            ifs >> keyword;
         }

         if ( keyword == "ettagons" )
         {
            fvca5_read_polygons(ifs, 8);
            ifs >> keyword;
         }

         if ( keyword == "edges" )
         {
            std::getline(ifs, keyword); //drop the rest of the line
            m_boundary_edges.clear();
            fvca5_read_boundary_edges(ifs);
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
            fvca5_read_edges(ifs);
         }
         else
         {
            std::cout << "Error parsing FVCA5 file" << std::endl;
            return false;
         }

         ifs.close();
         return true;
      }

   public:
      fvca5_mesh_loader() = default;

      bool read_mesh(const std::string& filename)
      {
         if (this->verbose())
            std::cout << " *** READING FVCA5 MESH ***" << std::endl;
         return fvca5_read(filename);
      }

      bool populate_mesh(mesh_type& msh)
      {
         if (this->verbose())
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
            assert(m_edges[i][0] < m_edges[i][1]);
            auto node1 = typename node_type::id_type(m_edges[i][0]);
            auto node2 = typename node_type::id_type(m_edges[i][1]);

            auto e = edge_type{{node1, node2}};

            e.set_point_ids(m_edges[i].begin(), m_edges[i].begin()+2); /* XXX: crap */
            edges.push_back(e);
         }
         /* Sort them */
         std::sort(edges.begin(), edges.end());

         /* Detect which ones are boundary edges */
         storage->boundary_info.resize(m_edges.size());
         for (size_t i = 0; i < m_boundary_edges.size(); i++)
         {
            assert(m_boundary_edges[i][0] < m_boundary_edges[i][1]);
            auto node1 = typename node_type::id_type(m_boundary_edges[i][0]);
            auto node2 = typename node_type::id_type(m_boundary_edges[i][1]);

            auto e = edge_type{{node1, node2}};

            auto position = find_element_id(edges.begin(), edges.end(), e);

            if (position.first == false)
            {
               std::cout << "Bad bug at " << __FILE__ << "("
               << __LINE__ << ")" << std::endl;
               return false;
            }

            bnd_info bi{0, true};
            storage->boundary_info.at(position.second) = bi;
         }

         storage->edges = std::move(edges);

         /* Surfaces */
         std::vector<surface_type> surfaces;
         surfaces.reserve( m_polys.size() );

         for (auto& p : m_polys)
         {
            std::vector<typename edge_type::id_type> surface_edges;
            for (auto& e : p.attached_edges)
            {
               assert(e[0] < e[1]);
               auto n1 = typename node_type::id_type(e[0]);
               auto n2 = typename node_type::id_type(e[1]);

               edge_type edge{{n1, n2}};
               auto edge_id = find_element_id(storage->edges.begin(),
                                              storage->edges.end(), edge);
               if (!edge_id.first)
               {
                  std::cout << "Bad bug at " << __FILE__ << "("
                  << __LINE__ << ")" << std::endl;
                  return false;
               }

               surface_edges.push_back(edge_id.second);
            }
            auto surface = surface_type(surface_edges);
            surface.set_point_ids(p.nodes.begin(), p.nodes.end()); /* XXX: crap */
            surfaces.push_back( surface );
         }

         std::sort(surfaces.begin(), surfaces.end());
         storage->surfaces = std::move(surfaces);

         return true;
      }
   };


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
            m_nodes.push_back( node_type( point_identifier<3>(i) ) );
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
                              const std::pair<size_t, edge_type>& e2) { return e1.second < e2.second;};
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
                             const std::pair<size_t, std::vector<size_t>>& e2) { return e1.second < e2.second;};

         std::sort(faces_to_edges.begin(), faces_to_edges.end(), comp_vecs);

         std::vector<surface_type> faces;
         faces.reserve( faces_to_edges.size() );
         conv_table.resize( faces_to_edges.size() );

         for (size_t i = 0; i < faces_to_edges.size(); i++)
         {
            auto fe = faces_to_edges[i];
            surface_type s( convert_to<typename edge_type::id_type>(fe.second) );
            s.set_point_ids( convert_to<point_identifier<3>>(faces_to_vts.at(fe.first)) );
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
            v.set_point_ids( convert_to<point_identifier<3>>(vol_to_vts.at(vf.first)) );
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

         bnd_info bi{0, true};
         storage->boundary_info.resize(storage->surfaces.size());
         for (size_t i = 0; i < storage->surfaces.size(); i++)
            if (bf[i] == 1)
               storage->boundary_info[i] = bi;

            return false;
      }
   };

} // namespace disk


#endif
