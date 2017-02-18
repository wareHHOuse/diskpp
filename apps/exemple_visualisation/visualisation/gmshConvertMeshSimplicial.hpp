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
      Gmesh msh;

      auto storage = mesh.backend_storage();


      // conversion point in node
      size_t nb_node(0);
      for(auto point :  storage->points){
         nb_node +=1;
         msh.addNode(convertPoint(point, nb_node));

         Vertice vert(nb_node, nb_node, 0, 0);
         msh.addVertice(vert);
      }
      assert(storage->points.size() == msh.getNumberofNodes());


      // conversion edge
      size_t nb_edge(0);
      for(auto edge :  storage->edges){
         nb_edge +=1;
         msh.addEdge(convertEdge(mesh, edge, nb_edge));
      }
      assert(storage->edges.size() == msh.getEdges().size());


      return msh;

   };


template<typename T>
Gmesh convertMesh(const disk::simplicial_mesh<T,2> mesh)
{
   Gmesh msh;

   auto storage = mesh.backend_storage();


   // conversion point in node
   size_t nb_node(0);
   for(auto point :  storage->points){
      nb_node +=1;
      msh.addNode(convertPoint(point, nb_node));

      Vertice vert(nb_node, nb_node, 0, 0);
      msh.addVertice(vert);
   }
   assert(storage->points.size() == msh.getNumberofNodes());


   // conversion edge
   size_t nb_edge(0);
   for(auto edge :  storage->edges){
      nb_edge +=1;
      msh.addEdge(convertEdge(edge, nb_edge));
   }
   assert(storage->edges.size() == msh.getEdges().size());

   // conversion triangle
   size_t nb_triangles(0);
   for(auto surface :  storage->surfaces){
      nb_triangles +=1;
      msh.addTriangle(convertTriangle(mesh, surface, nb_triangles));
   }
   assert(storage->surfaces.size() == msh.getTriangles().size());

   assert((storage->edges.size() + storage->surfaces.size() + storage->nodes.size())
         == msh.getNumberofElements());

   return msh;

};

template<typename T>
Gmesh convertMesh(const disk::simplicial_mesh<T,3>& mesh)
{
   Gmesh msh;

   auto storage = mesh.backend_storage();


   // conversion point in node
   size_t nb_node(0);
   for(auto point :  storage->points){
      nb_node +=1;
      msh.addNode(convertPoint(point, nb_node));

      Vertice vert(nb_node, nb_node, 0, 0);
      msh.addVertice(vert);
   }
   assert(storage->points.size() == msh.getNumberofNodes());


   // conversion edge
   size_t nb_edge(0);
   for(auto edge :  storage->edges){
      nb_edge +=1;
      msh.addEdge(convertEdge(edge, nb_edge));
   }
   assert(storage->edges.size() == msh.getEdges().size());

   // conversion triangle
   size_t nb_triangles(0);
   for(auto surface :  storage->surfaces){
      nb_triangles +=1;
      msh.addTriangle(convertTriangle(surface, nb_triangles));
   }
   assert(storage->surfaces.size() == msh.getTriangles().size());


   // conversion tetra
   size_t nb_tetras(0);
   for(auto volume :  storage->volumes){
      nb_tetras +=1;
      msh.addTetrahedron(convertTetrahedron(mesh, volume, nb_tetras));
   }
   assert(storage->volumes.size() == msh.getTetrahedra().size());

   assert((storage->edges.size() + storage->surfaces.size() + storage->nodes.size()
         + storage->volumes.size()) == msh.getNumberofElements());

   return msh;

};


} //visu

#endif
