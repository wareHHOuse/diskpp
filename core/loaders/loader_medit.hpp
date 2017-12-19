/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
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

#ifndef _LOADER_HPP_WAS_INCLUDED_
#error "You must NOT include this file directly. Include loader.hpp"
#endif

#ifndef _LOADER_MEDIT_HPP_
#define _LOADER_MEDIT_HPP_

#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <set>
#include <thread>
#include <vector>

#include "geometry/geometry.hpp"
#include "loader_template.hpp"
#include "loader_utils.hpp"

#include "mapped_file.h"
#include "strtot.hpp"

namespace disk {

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
   std::vector<std::array<ident_impl_t, 2>>                    m_edges;
   std::vector<medit2d_poly>                                   m_polys;
   std::vector<std::pair<std::array<ident_impl_t, 2>, size_t>> m_boundary_edges;

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
         m_nodes.push_back(node_type(point_identifier<2>(i)));
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
         std::vector<ident_impl_t> nodes(polynum + 1, 0);

         for (size_t j = 0; j < polynum; j++) {
            ident_impl_t val;
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
            std::array<ident_impl_t, 2> b_edge = {nodes[j], nodes[j + 1]};
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
         std::array<ident_impl_t, 2> b_edge;
         ifs >> b_edge[0];
         b_edge[0] -= 1;
         ifs >> b_edge[1];
         b_edge[1] -= 1;

         assert(b_edge[0] != b_edge[1]);

         if (b_edge[0] > b_edge[1]) std::swap(b_edge[0], b_edge[1]);

         std::array<ident_impl_t, 2> bnd = {b_edge[0], b_edge[1]};

         size_t b_id;
         ifs >> b_id;

         m_boundary_edges.push_back(std::make_pair(bnd, b_id));
      }

      return true;
   }

   bool
   face_unique(const std::array<ident_impl_t, 2>& f1, const std::array<ident_impl_t, 2>& f2)
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

            auto e = edge_type{{node1, node2}};

            e.set_point_ids(m_edges[i].begin(), m_edges[i].begin() + 2); /* XXX: crap */
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

         auto e = edge_type{{node1, node2}};

         auto position = find_element_id(edges.begin(), edges.end(), e);

         if (position.first == false) {
            std::cout << "Bad bug at " << __FILE__ << "(" << __LINE__ << ")" << std::endl;
            return false;
         }

         bnd_info bi{m_boundary_edges[i].second, true};
         storage->boundary_info.at(position.second) = bi;
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

            edge_type edge{{n1, n2}};
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

   std::vector<std::array<size_t, 2>> tab_edges;
   std::vector<std::array<size_t, 4>> tab_faces;

   bool
   edge_exist(const std::array<size_t, 2>& f, size_t& ind)
   {
      for (size_t i = 0; i < tab_edges.size(); i++) {
         const std::array<size_t, 2> fi = tab_edges[i];

         if (f[0] == fi[0]) {
            if (f[1] == fi[1]) {
               ind = i;
               return true;
            } else
               return false;
         } else
            return false;
      }

      return false;
   }

   bool
   tria_exist(const std::array<size_t, 4>& t, size_t& ind)
   {
      for (size_t i = 0; i < tab_faces.size(); i++) {
         const std::array<size_t, 4> fi = tab_faces[i];

         if (t[0] == fi[0]) {
            if (t[1] == fi[1]) {
               if (t[2] == fi[2]) {
                  ind = i;
                  return true;
               } else
                  return false;
            } else
               return false;
         } else
            return false;
      }

      return false;
   }

   bool
   quad_exist(const std::array<size_t, 4>& q, size_t& ind)
   {
      for (size_t i = 0; i < tab_faces.size(); i++) {
         const std::array<size_t, 4> fi = tab_faces[i];

         if (q[0] == fi[0]) {
            if (q[1] == fi[1]) {
               if (q[2] == fi[2]) {
                  if (q[3] == fi[3]) {
                     ind = i;
                     return true;
                  } else
                     return false;
               } else
                  return false;
            } else
               return false;
         } else
            return false;
      }
   }

   size_t
   define_tria(const std::array<size_t, 3>& nodes)
   {
      std::array<size_t, 4> face_nodes = {0, 0, 0, 0};
      std::vector<size_t>   face_vts(3, 0);
      size_t                vol_faces;

      size_t ind(0);

      // Face  =
      face_nodes[0] = nodes[0];
      face_nodes[1] = nodes[1];
      face_nodes[2] = nodes[2];
      face_nodes[3] = nodes[2];
      std::sort(face_nodes.begin(), face_nodes.begin() + 3);

      if (tria_exist(face_nodes, ind)) {
         vol_faces = ind;
      } else {
         // We create a face
         face_vts[0] = nodes[0];
         face_vts[1] = nodes[1];
         face_vts[2] = nodes[2];
         vol_faces   = faces_to_vts.size();
         faces_to_vts.push_back(std::move(face_vts));

         tab_faces.push_back(face_nodes);

         std::array<size_t, 2> edge_nodes = {0, 0};
         std::vector<size_t>   face_edges(3, 0);
         size_t                ind2(0);

         // edge 1
         edge_nodes[0] = nodes[0];
         edge_nodes[1] = nodes[1];
         std::sort(edge_nodes.begin(), edge_nodes.end());

         if (edge_exist(edge_nodes, ind2)) {
            face_edges[0] = ind2;
         } else {
            auto e = edge_type({typename node_type::id_type(edge_nodes[0]),
                                typename node_type::id_type(edge_nodes[1])});

            tab_edges.push_back(edge_nodes);
            face_edges[0] = m_edges.size();
            m_edges.push_back(std::make_pair(m_edges.size(), e));
         }

         // edge 2
         edge_nodes[0] = nodes[1];
         edge_nodes[1] = nodes[2];
         std::sort(edge_nodes.begin(), edge_nodes.end());

         if (edge_exist(edge_nodes, ind2)) {
            face_edges[1] = ind2;
         } else {
            auto e = edge_type({typename node_type::id_type(edge_nodes[0]),
                                typename node_type::id_type(edge_nodes[1])});

            tab_edges.push_back(edge_nodes);
            face_edges[1] = m_edges.size();
            m_edges.push_back(std::make_pair(m_edges.size(), e));
         }

         // edge 3
         edge_nodes[0] = nodes[2];
         edge_nodes[1] = nodes[0];
         std::sort(edge_nodes.begin(), edge_nodes.end());

         if (edge_exist(edge_nodes, ind2)) {
            face_edges[2] = ind2;
         } else {
            auto e = edge_type({typename node_type::id_type(edge_nodes[0]),
                                typename node_type::id_type(edge_nodes[1])});

            tab_edges.push_back(edge_nodes);
            face_edges[2] = m_edges.size();
            m_edges.push_back(std::make_pair(m_edges.size(), e));
         }

         faces_to_edges.push_back(std::make_pair(vol_faces, std::move(face_edges)));
      }

      return vol_faces;
   }

   size_t
   define_quad(const std::array<size_t, 4>& nodes)
   {
      std::array<size_t, 4> face_nodes = {0, 0, 0, 0};
      std::vector<size_t>   face_vts(4, 0);
      size_t                vol_faces;

      size_t ind(0);

      // Face  =
      face_nodes[0] = nodes[0];
      face_nodes[1] = nodes[1];
      face_nodes[2] = nodes[2];
      face_nodes[3] = nodes[3];
      std::sort(face_nodes.begin(), face_nodes.end());

      if (quad_exist(face_nodes, ind)) {
         vol_faces = ind;
      } else {
         // We create a face
         face_vts[0] = nodes[0];
         face_vts[1] = nodes[1];
         face_vts[2] = nodes[2];
         face_vts[3] = nodes[3];
         vol_faces   = faces_to_vts.size();
         faces_to_vts.push_back(std::move(face_vts));
         tab_faces.push_back(face_nodes);

         std::array<size_t, 2> edge_nodes = {0, 0};
         std::vector<size_t>   face_edges(4, 0);
         size_t                ind2(0);

         // edge 1
         edge_nodes[0] = nodes[0];
         edge_nodes[1] = nodes[1];
         std::sort(edge_nodes.begin(), edge_nodes.end());

         if (edge_exist(edge_nodes, ind2)) {
            face_edges[0] = ind2;
         } else {
            auto e = edge_type({typename node_type::id_type(edge_nodes[0]),
                                typename node_type::id_type(edge_nodes[1])});

            tab_edges.push_back(edge_nodes);
            face_edges[0] = m_edges.size();
            m_edges.push_back(std::make_pair(m_edges.size(), e));
         }

         // edge 2
         edge_nodes[0] = nodes[1];
         edge_nodes[1] = nodes[2];
         std::sort(edge_nodes.begin(), edge_nodes.end());

         if (edge_exist(edge_nodes, ind2)) {
            face_edges[1] = ind2;
         } else {
            auto e = edge_type({typename node_type::id_type(edge_nodes[0]),
                                typename node_type::id_type(edge_nodes[1])});

            tab_edges.push_back(edge_nodes);
            face_edges[1] = m_edges.size();
            m_edges.push_back(std::make_pair(m_edges.size(), e));
         }

         // edge 3
         edge_nodes[0] = nodes[2];
         edge_nodes[1] = nodes[3];
         std::sort(edge_nodes.begin(), edge_nodes.end());

         if (edge_exist(edge_nodes, ind2)) {
            face_edges[2] = ind2;
         } else {
            auto e = edge_type({typename node_type::id_type(edge_nodes[0]),
                                typename node_type::id_type(edge_nodes[1])});

            tab_edges.push_back(edge_nodes);
            face_edges[2] = m_edges.size();
            m_edges.push_back(std::make_pair(m_edges.size(), e));
         }

         // edge 4
         edge_nodes[0] = nodes[3];
         edge_nodes[1] = nodes[0];
         std::sort(edge_nodes.begin(), edge_nodes.end());

         if (edge_exist(edge_nodes, ind2)) {
            face_edges[3] = ind2;
         } else {
            auto e = edge_type({typename node_type::id_type(edge_nodes[0]),
                                typename node_type::id_type(edge_nodes[1])});

            tab_edges.push_back(edge_nodes);
            face_edges[3] = m_edges.size();
            m_edges.push_back(std::make_pair(m_edges.size(), e));
         }

         faces_to_edges.push_back(std::make_pair(vol_faces, std::move(face_edges)));
      }

      return vol_faces;
   }

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
         m_points.push_back(point_type{x, y, z});
         m_nodes.push_back(node_type(point_identifier<3>(i)));
      }

      return true;
   }

   bool
   medit_read_tetra(std::ifstream& ifs)
   {
      size_t       elements_to_read;
      const size_t nb_nodes = 4;

      ifs >> elements_to_read;
      if (this->verbose())
         std::cout << "Reading " << elements_to_read << " Tetrahedra" << std::endl;

      for (size_t i = 0; i < elements_to_read; i++) {

         // Read nodes
         std::vector<size_t> nodes(nb_nodes, 0);

         for (size_t j = 0; j < nb_nodes; j++) {
            size_t val;
            ifs >> val;
            nodes.push_back(val - 1);
         }

         // Read ref of the tetra
         size_t p;
         ifs >> p;

         std::array<size_t, 3> face_nodes = {0, 0, 0};
         std::vector<size_t>   vol_faces(4, 0);

         // Face1
         face_nodes[0] = nodes[0];
         face_nodes[1] = nodes[1];
         face_nodes[2] = nodes[2];
         vol_faces[0]  = define_tria(face_nodes);

         // Face2
         face_nodes[0] = nodes[0];
         face_nodes[1] = nodes[3];
         face_nodes[2] = nodes[2];
         vol_faces[1]  = define_tria(face_nodes);

         // Face3
         face_nodes[0] = nodes[0];
         face_nodes[1] = nodes[3];
         face_nodes[2] = nodes[1];
         vol_faces[2]  = define_tria(face_nodes);

         // Face4
         face_nodes[0] = nodes[3];
         face_nodes[1] = nodes[1];
         face_nodes[2] = nodes[2];
         vol_faces[3]  = define_tria(face_nodes);

         // save nodes
         vol_to_faces.push_back(std::make_pair(vol_to_faces.size(), std::move(vol_faces)));
         vol_to_vts.push_back(std::move(nodes));
      }
      return true;
   }

   bool
   medit_read_hexa(std::ifstream& ifs)
   {
      size_t       elements_to_read;
      const size_t nb_nodes = 8;

      ifs >> elements_to_read;
      if (this->verbose()) std::cout << "Reading " << elements_to_read << " Hexahedra" << std::endl;

      for (size_t i = 0; i < elements_to_read; i++) {

         // Read nodes
         std::vector<size_t> nodes(nb_nodes, 0);

         for (size_t j = 0; j < nb_nodes; j++) {
            size_t val;
            ifs >> val;
            nodes.push_back(val - 1);
         }

         // Read ref of the tetra
         size_t p;
         ifs >> p;

         std::array<size_t, 4> face_nodes = {0, 0, 0, 0};
         std::vector<size_t>   vol_faces(6, 0);

         // Face1
         face_nodes[0] = nodes[0];
         face_nodes[1] = nodes[1];
         face_nodes[2] = nodes[2];
         face_nodes[3] = nodes[3];
         vol_faces[0]  = define_quad(face_nodes);

         // Face2
         face_nodes[0] = nodes[0];
         face_nodes[1] = nodes[4];
         face_nodes[2] = nodes[5];
         face_nodes[3] = nodes[1];
         vol_faces[1]  = define_quad(face_nodes);

         // Face3
         face_nodes[0] = nodes[1];
         face_nodes[1] = nodes[2];
         face_nodes[2] = nodes[6];
         face_nodes[3] = nodes[5];
         vol_faces[2]  = define_quad(face_nodes);

         // Face4
         face_nodes[0] = nodes[0];
         face_nodes[1] = nodes[4];
         face_nodes[2] = nodes[7];
         face_nodes[3] = nodes[3];
         vol_faces[3]  = define_quad(face_nodes);

         // Face5
         face_nodes[0] = nodes[3];
         face_nodes[1] = nodes[7];
         face_nodes[2] = nodes[6];
         face_nodes[3] = nodes[2];
         vol_faces[4]  = define_quad(face_nodes);

         // Face6
         face_nodes[0] = nodes[4];
         face_nodes[1] = nodes[5];
         face_nodes[2] = nodes[6];
         face_nodes[3] = nodes[7];
         vol_faces[5]  = define_quad(face_nodes);

         // save nodes
         vol_to_faces.push_back(std::make_pair(vol_to_faces.size(), std::move(vol_faces)));
         vol_to_vts.push_back(std::move(nodes));
      }
      return true;
   }

   bool
   medit_read_edges(std::ifstream& ifs)
   {
      size_t elements_to_read;

      ifs >> elements_to_read;

      std::string keyword;

      // We don't use this edges
      size_t v1, v2, ref;

      for (size_t i = 0; i < elements_to_read; i++)
         ifs >> v1 >> v2 >> ref;

      return true;
   }

   bool
   medit_read_boundary_tria(std::ifstream& ifs)
   {
      size_t elements_to_read;

      ifs >> elements_to_read;

      if (this->verbose())
         std::cout << "Reading " << elements_to_read << " boundary triangles" << std::endl;

      if (m_boundary_edges.empty())
         m_boundary_edges.reserve(elements_to_read);
      else
         m_boundary_edges.resize(m_boundary_edges.size() + elements_to_read);

      for (size_t i = 0; i < elements_to_read; i++) {
         std::array<size_t, 3> face_nodes = {0, 0, 0};

         for (size_t j = 0; j < 3; j++) {
            size_t val;
            ifs >> val;
            face_nodes[j] = val - 1;
         }

         // Face1
         size_t face_id = define_tria(face_nodes);

         size_t b_id;
         ifs >> b_id;

         std::array<size_t, 2> bnd = {face_id, b_id};

         m_boundary_edges.push_back(bnd);
      }

      return true;
   }

   bool
   medit_read_boundary_quad(std::ifstream& ifs)
   {
      size_t elements_to_read;

      ifs >> elements_to_read;
      if (this->verbose())
         std::cout << "Reading " << elements_to_read << " boundary quads" << std::endl;

      if (m_boundary_edges.empty())
         m_boundary_edges.reserve(elements_to_read);
      else
         m_boundary_edges.resize(m_boundary_edges.size() + elements_to_read);

      for (size_t i = 0; i < elements_to_read; i++) {
         std::array<size_t, 4> face_nodes = {0, 0, 0};

         for (size_t j = 0; j < 4; j++) {
            size_t val;
            ifs >> val;
            face_nodes[j] = val - 1;
         }

         // Face1
         size_t face_id = define_quad(face_nodes);

         size_t b_id;
         ifs >> b_id;

         std::array<size_t, 2> bnd = {face_id, b_id};

         m_boundary_edges.push_back(bnd);
      }

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
         std::cout << "Expected dimension = 3 (here: " << dim << ")" << std::endl;
         return false;
      }

      ifs >> keyword;
      while (keyword != "End") {
         std::cout << keyword << std::endl;
         if (keyword == "Vertices") {
            medit_read_vertices(ifs);
         } else if (keyword == "Triangles") {
            medit_read_boundary_tria(ifs);
         } else if (keyword == "Quadrilaterals") {
            medit_read_boundary_quad(ifs);
         } else if (keyword == "Edges") {
            medit_read_edges(ifs);
         } else if (keyword == "Tetrahedra") {
            medit_read_tetra(ifs);
         } else if (keyword == "Hexahedra") {
            medit_read_hexa(ifs);
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
   medit_mesh_loader() = default;

   bool
   read_mesh(const std::string& s)
   {
      // if (this->verbose())
      std::cout << " *** READING MEDIT 3D MESH ***" << std::endl;

      return medit_read(s);
   }

   bool
   populate_mesh(mesh_type& msh)
   {
      // if (this->verbose())
      std::cout << " *** POPULATING MEDIT MESH ***" << std::endl;
      auto storage = msh.backend_storage();

      std::vector<size_t> conv_table;
      /* Sort the edges in lexicographical order, remember their original
       * position to convert the pointers in the faces */
      auto comp_edges = [](const std::pair<size_t, edge_type>& e1,
                           const std::pair<size_t, edge_type>& e2) {
         return e1.second < e2.second;
      };
      std::sort(m_edges.begin(), m_edges.end(), comp_edges);

      std::cout << "Guessed mesh format: Medit format" << std::endl;
      std::vector<edge_type> edges;
      edges.reserve(m_edges.size());
      conv_table.resize(m_edges.size());
      for (size_t i = 0; i < m_edges.size(); i++) {
         conv_table[m_edges[i].first] = i; /* Make ptr conversion table */
         edges.push_back(m_edges[i].second);
      }

      std::cout << "Guessed mesh format: Medit format" << std::endl;
      /* Convert the edge pointers in the face data */
      for (auto& fe : faces_to_edges) {
         for (auto& ptr : fe.second)
            ptr = conv_table[ptr];
      }
      std::cout << "Guessed mesh format: Medit format" << std::endl;
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
         s.set_point_ids(convert_to<point_identifier<3>>(faces_to_vts.at(fe.first)));
         faces.push_back(s);
         conv_table[fe.first] = i;
      }
      /* Now the faces are in their place and have correct ptrs */
      std::cout << "Guessed mesh format: Medit format" << std::endl;

      /* Detect which ones are boundary edges */
      storage->boundary_info.resize(faces.size());
      for (size_t i = 0; i < m_boundary_edges.size(); i++) {
         std::array<size_t, 2> bnd      = m_boundary_edges[i];
         size_t                position = conv_table[i];

         bnd_info bi{bnd[1], true};
         storage->boundary_info.at(position) = bi;
      }

      std::cout << "Guessed mesh format: Medit format" << std::endl;

      /* Convert the face pointers in the volume data */
      for (auto& vf : vol_to_faces) {
         for (auto& ptr : vf.second)
            ptr = conv_table[ptr];
      }

      // for (auto& f : faces)
      //    std::cout << f << std::endl;

      /* Sort volume data */
      std::sort(vol_to_faces.begin(), vol_to_faces.end(), comp_vecs);

      std::vector<volume_type> volumes;
      volumes.reserve(vol_to_faces.size());

      for (size_t i = 0; i < vol_to_faces.size(); i++) {
         auto        vf = vol_to_faces[i];
         volume_type v(convert_to<typename surface_type::id_type>(vf.second));
         v.set_point_ids(convert_to<point_identifier<3>>(vol_to_vts.at(vf.first)));
         volumes.push_back(v);
      }
      std::cout << "Guessed mesh format: Medit format" << std::endl;

      storage->points   = std::move(m_points);
      storage->nodes    = std::move(m_nodes);
      storage->edges    = std::move(edges);
      storage->surfaces = std::move(faces);
      storage->volumes  = std::move(volumes);

      return true;
   }
};

} // namespace disk

#endif
