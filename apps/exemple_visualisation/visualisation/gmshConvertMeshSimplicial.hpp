/*
 *       /\
 *      /__\       Nicolas Pignet (C) 2017
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    Gmsh tools
 *  /_\/_\/_\/_\
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */


#ifndef _GmshConvertMesh_HPP
     #error "You must NOT include this file. Include gmshConvertMesh.hpp"
#endif

#ifndef _GmshConvertMeshSimplicial_HPP
#define _GmshConvertMeshSimplicial_HPP

#include<vector>
#include<string>
#include<utility>
#include <array>
#include <cassert>

#include "loaders/loader.hpp"
#include "geometry/geometry.hpp"
#include "mesh/point.hpp"
#include "gmshMesh.h"
#include "gmshElement.h"
#include "gmshDisk.h"

namespace visu{

template<typename T>
Gmesh convertMesh(const disk::simplicial_mesh<T,1> mesh)
{
      auto storage = mesh.backend_storage();

      std::vector<Node> nodes;
      nodes.reserve(storage->points.size());
      std::vector<Vertice> vertices;
      vertices.reserve(storage->points.size());
      std::vector<Edge> edges;
      edges.reserve(storage->edges.size());

      // conversion point in node
      size_t nb_node(0);
      for(auto point :  storage->points){
         nb_node +=1;
         nodes.push_back(convertPoint(point, nb_node));

         Vertice vert(nb_node, nb_node, 0, 0);
         vertices.push_back(vert);
      }
      assert(storage->points.size() == nodes.size());


      // conversion edge
      size_t nb_edge(0);
      for(auto edge :  storage->edges){
         nb_edge +=1;
         edges.push_back(convertEdge(mesh, edge, nb_edge));
      }
      assert(storage->edges.size() == edges.size());

      Gmesh msh(1, nodes, vertices, edges);

      return msh;
};


template<typename T>
Gmesh convertMesh(const disk::simplicial_mesh<T,2> mesh)
{
   auto storage = mesh.backend_storage();

   std::vector<Node> nodes;
   nodes.reserve(storage->points.size());
   std::vector<Vertice> vertices;
   vertices.reserve(storage->points.size());
   std::vector<Edge> edges;
   edges.reserve(storage->edges.size());
   std::vector<Triangle> triangles;
   triangles.reserve(storage->surfaces.size());
   std::vector<Quadrangle> quadrangles;


   // conversion point in node
   size_t nb_node(0);
   for(auto point :  storage->points){
      nb_node +=1;
      nodes.push_back(convertPoint(point, nb_node));

      Vertice vert(nb_node, nb_node, 0, 0);
      vertices.push_back(vert);
   }
   assert(storage->points.size() == nodes.size());


   // conversion edge
   size_t nb_edge(0);
   for(auto edge :  storage->edges){
      nb_edge +=1;
      edges.push_back(convertEdge(mesh, edge, nb_edge));
   }
   assert(storage->edges.size() == edges.size());

   // conversion triangle
   size_t nb_triangles(0);
   for(auto surface :  storage->surfaces){
      nb_triangles +=1;
      triangles.push_back(convertTriangle(mesh, surface, nb_triangles));
   }
   assert(storage->surfaces.size() == triangles.size());

   Gmesh msh(2, nodes, vertices, edges, triangles, quadrangles);

   return msh;

};

template<typename T>
Gmesh convertMesh(const disk::simplicial_mesh<T,3>& mesh)
{
   auto storage = mesh.backend_storage();

   std::vector<Node> nodes;
   nodes.reserve(storage->points.size());
   std::vector<Vertice> vertices;
   vertices.reserve(storage->points.size());
   std::vector<Edge> edges;
   edges.reserve(storage->edges.size());
   std::vector<Triangle> triangles;
   triangles.reserve(storage->surfaces.size());
   std::vector<Tetrahedron> tetrahedra;
   tetrahedra.reserve(storage->volumes.size());

   // conversion point in node
   size_t nb_node(0);
   for(auto point :  storage->points){
      nb_node +=1;
      nodes.push_back(convertPoint(point, nb_node));

      Vertice vert(nb_node, nb_node, 0, 0);
      vertices.push_back(vert);
   }
   assert(storage->points.size() == nodes.size());


   // conversion edge
   size_t nb_edge(0);
   for(auto edge :  storage->edges){
      nb_edge +=1;
      edges.push_back(convertEdge(mesh, edge, nb_edge));
   }
   assert(storage->edges.size() == edges.size());


   // I dont know how ordering the nodes
   // conversion triangle
   // size_t nb_triangles(0);
   // for(auto surface :  storage->surfaces){
   //    nb_triangles +=1;
   //    triangles.push_back(convertTriangle(mesh, surface, nb_triangles));
   // }
   // assert(storage->surfaces.size() == triangles.size());


   // conversion tetra
   size_t nb_tetras(0);
   for(auto volume :  storage->volumes){
      nb_tetras +=1;
      tetrahedra.push_back(convertTetrahedron(mesh, volume, nb_tetras));
   }
   assert(storage->volumes.size() == tetrahedra.size());

   Gmesh msh(3, nodes, vertices, edges, triangles, tetrahedra);

   return msh;
};


} //visu

#endif
