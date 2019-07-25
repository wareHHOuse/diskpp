/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2018                     nicolas.pignet@enpc.fr
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

#include <iostream>
#include <sstream>

#include <list>
#include <vector>

#include "mesh/mesh.hpp"
#include "mesh/point.hpp"

namespace disk {

// Class to create mesh for post-processing
template<typename Mesh>
struct PostMesh
{
   static_assert(sizeof(Mesh) == -1, "post_mesh: not suitable for the requested kind of mesh");
};

template<typename T>
class PostMesh<generic_mesh<T, 1>>
{
 public:
   typedef generic_mesh<T, 1>               mesh_type;
   typedef typename mesh_type::point_type   point_type;
   typedef typename mesh_type::node_type    node_type;
   typedef typename mesh_type::edge_type    edge_type;
   typedef std::vector<point_identifier<1>> list_type;

 private:
   // post-mesh
   mesh_type              post_mesh;
   std::vector<list_type> list_cell_nodes;

 public:
   PostMesh() { list_cell_nodes.clear(); }
   PostMesh(const generic_mesh<T, 1>& msh) : post_mesh(msh)
   {
      list_cell_nodes.reserve(msh.cells_size());
      for (auto& cl : msh) {
         list_cell_nodes.push_back(cell_nodes(msh, cl));
      }
   };

   mesh_type
   mesh() const
   {
      return post_mesh;
   }

   list_type
   nodes_cell(const size_t cell_id) const
   {
      return list_cell_nodes.at(cell_id);
   }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class PostMesh<Mesh<T, 2, Storage>>
{
 public:
   typedef simplicial_mesh<T, 2>            mesh_type;
   typedef typename mesh_type::point_type   point_type;
   typedef typename mesh_type::node_type    node_type;
   typedef typename mesh_type::edge_type    edge_type;
   typedef typename mesh_type::surface_type surface_type;
   typedef std::vector<disk::point_identifier<2>> list_type;

 private:
   // post-mesh
   mesh_type              post_mesh;
   std::vector<list_type> list_cell_nodes;

 public:
   PostMesh() { list_cell_nodes.clear(); }

   PostMesh(const Mesh<T, 2, Storage>& msh)
   {
      auto storage_in  = msh.backend_storage();
      auto storage_out = post_mesh.backend_storage();

      size_t num_points = storage_in->nodes.size();

      // we copy old points and nodes
      storage_out->points.reserve(num_points);
      storage_out->nodes.reserve(num_points);

      for (size_t i = 0; i < num_points; i++) {
         storage_out->points.push_back(storage_in->points[i]);
         const auto point_id = disk::point_identifier<2>(i);
         storage_out->nodes.push_back(node_type({point_id}));
      }

      // we copy old edges
      size_t num_edges = storage_in->edges.size();
      storage_out->edges.reserve(num_edges);

      for (auto& edge : storage_in->edges) {
         const auto pt = edge.point_ids();
         storage_out->edges.push_back(edge_type({pt[0], pt[1]}));
      }

      list_cell_nodes.reserve(storage_in->surfaces.size());
      assert(storage_in->surfaces.size() == msh.cells_size());
      // split all surfaces in triangles
      for (auto& cl : msh) {
         auto pts_cell = cell_nodes(msh, cl);

         list_type pts;

         for (auto pt : pts_cell)
            pts.push_back(pt);

         // If it is a triangle we save it
         if (pts.size() == 3) {
            storage_out->surfaces.push_back(surface_type({pts[0], pts[1], pts[2]}));
         } else {
            // we split element  in triangles with barycenter as node
            const auto bar = barycenter(msh, cl);
            storage_out->points.push_back(bar);
            const auto bar_id = disk::point_identifier<2>(num_points++);
            pts.push_back(bar_id);
            storage_out->nodes.push_back(node_type({bar_id}));

            auto fcs = faces(msh, cl);
            for (auto& fc : fcs) {
               const auto pts_fc = face_nodes(msh, fc);
               storage_out->surfaces.push_back(surface_type({pts_fc[0], pts_fc[1], bar_id}));
            }
         }
         list_cell_nodes.push_back(pts);
      }
   }

   mesh_type
   mesh() const
   {
      return post_mesh;
   }

   list_type
   nodes_cell(const size_t cell_id) const
   {
      return list_cell_nodes.at(cell_id);
   }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class PostMesh<Mesh<T, 3, Storage>>
{
 public:
   typedef simplicial_mesh<T, 3>    mesh_type;
   typedef typename mesh_type::cell cell_type;
   typedef typename mesh_type::face face_type;

   typedef typename mesh_type::point_type   point_type;
   typedef typename mesh_type::node_type    node_type;
   typedef typename mesh_type::edge_type    edge_type;
   typedef typename mesh_type::surface_type surface_type;
   typedef typename mesh_type::volume_type  volume_type;
   typedef disk::point_identifier<3>              point_ident;
   typedef std::vector<point_ident>         list_type;

 private:
   // post-mesh
   mesh_type              post_mesh;
   std::vector<list_type> list_cell_nodes;
   std::vector<list_type> list_face_nodes;
   size_t                 num_points;

   typedef Mesh<T, 3, Storage>        meshin_type;
   typedef typename meshin_type::cell cellin_type;
   typedef typename meshin_type::face facein_type;

   struct triangle
   {
      triangle() {}
      triangle(point_ident ap0, point_ident ap1, point_ident ap2) : p0(ap0), p1(ap1), p2(ap2) {}

      point_ident p0, p1, p2;
   };

   // a ordonner
   void
   split_in_triangles(const meshin_type& msh)
   {
      auto storage_out = post_mesh.backend_storage();

      list_face_nodes.reserve(msh.faces_size());

      for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++) {
         auto fc       = *itor;
         auto pts_face = face_nodes(msh, fc);

         list_type pts;

         for (auto pt : pts_face)
            pts.push_back(pt);

         if (pts.size() > 3) {
            const auto bar = barycenter(msh, fc);
            storage_out->points.push_back(bar);
            const auto bar_id = point_ident(num_points++);
            pts.push_back(bar_id);
            storage_out->nodes.push_back(node_type({bar_id}));
         }

         list_face_nodes.push_back(pts);
      }
   }

   std::vector<triangle>
   split_in_triangles(const meshin_type& msh, const facein_type& fc)
   {
      std::vector<triangle> ret;
      ret.clear();

      list_type pts = list_face_nodes[msh.lookup(fc)];

      if (pts.size() == 3) {
         triangle t(pts[0], pts[1], pts[2]);
         ret.push_back(t);
      } else {
         const auto bar = pts.back();
         for (size_t i = 0; i < (pts.size() - 1); i++) {
            triangle t(pts[i], pts[(i + 1) % (pts.size() - 1)], bar);
            ret.push_back(t);
         }
      }

      return ret;
   }

 public:
   PostMesh() : num_points(0)
   {
      list_cell_nodes.clear();
      list_face_nodes.clear();
   }

   PostMesh(const Mesh<T, 3, Storage>& msh)
   {
      auto storage_in  = msh.backend_storage();
      auto storage_out = post_mesh.backend_storage();

      num_points = storage_in->nodes.size();

      // we copy old points and nodes
      storage_out->points.reserve(num_points);
      storage_out->nodes.reserve(num_points);

      for (size_t i = 0; i < num_points; i++) {
         storage_out->points.push_back(storage_in->points[i]);
         const auto point_id = point_ident(i);
         storage_out->nodes.push_back(node_type({point_id}));
      }

      // we copy old edges
      size_t num_edges = storage_in->edges.size();
      storage_out->edges.reserve(num_edges);

      if (num_edges > 0) {
         for (auto& edge : storage_in->edges) {
            const auto pt = edge.point_ids();
            assert(pt.size() == 2);
            storage_out->edges.push_back(edge_type({pt[0], pt[1]}));
         }
      }

      // split all surfaces in triangles
      split_in_triangles(msh);

      list_cell_nodes.reserve(storage_in->volumes.size());
      assert(storage_in->volumes.size() == msh.cells_size());

      // split all volumes in tetra
      for (auto& cl : msh) {
         auto pts_cell = cell_nodes(msh, cl);

         list_type pts;

         for (auto pt : pts_cell)
            pts.push_back(pt);

         // If it is a tetra we save it
         if (pts.size() == 4) {
            storage_out->volumes.push_back(volume_type({pts[0], pts[1], pts[2], pts[3]}));
         } else {
            // we split element  in tetra with barycenter as node
            const auto bar = barycenter(msh, cl);
            storage_out->points.push_back(bar);
            const auto bar_id = point_ident(num_points++);
            pts.push_back(bar_id);
            storage_out->nodes.push_back(node_type({bar_id}));

            auto fcs = faces(msh, cl);
            for (auto& fc : fcs) {
               std::vector<triangle> triangles = split_in_triangles(msh, fc);

               if (triangles.size() > 1) {
                  pts.push_back(triangles[1].p2);
               }

               for (auto& tri : triangles) {
                  storage_out->volumes.push_back(volume_type({tri.p0, tri.p1, tri.p2, bar_id}));
               }
            }
         }

         list_cell_nodes.push_back(pts);
      }
   }

   mesh_type
   mesh() const
   {
      return post_mesh;
   }

   list_type
   nodes_cell(const size_t cell_id) const
   {
      return list_cell_nodes.at(cell_id);
   }
};

} // end disk