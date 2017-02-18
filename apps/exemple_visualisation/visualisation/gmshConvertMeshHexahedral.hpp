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

#ifndef _GmshConvertMeshHexa_HPP
#define _GmshConvertMeshHexa_HPP

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
Gmesh convertMesh(const disk::cartesian_mesh<T,1>& mesh)
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
Gmesh convertMesh(const disk::cartesian_mesh<T,2>& mesh)
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

   // conversion quad
   size_t nb_surface(0);
   for(auto& cl : mesh ){
      auto fcs = disk::faces(mesh, cl);
      nb_surface += 1;

      auto face_id = mesh.lookup(fcs[0]);
      size_t ind0 = (storage->edges[face_id].point_ids())[0];
      size_t ind1 = (storage->edges[face_id].point_ids())[1];
      auto p0 = storage->points[ind0];
      auto p1 = storage->points[ind1];

      auto p2 = storage->points[ind0];
      auto p3 = storage->points[ind0];

      size_t ind2(0), ind3(0);

      // we search p2 and p3
      for(size_t i = 1; i < fcs.size(); i++){
         auto nodes = storage->edges[mesh.lookup(fcs[i])].point_ids();
         if( nodes[0] == ind0){
            p3 = storage->points[nodes[1]];
            ind3 = nodes[1];
         }
         else if( nodes[1] == ind0){
            p3 = storage->points[nodes[0]];
            ind3 = nodes[0];
         }
         else if( nodes[0] == ind1){
            p2 = storage->points[nodes[1]];
            ind2 = nodes[1];
         }
         else if( nodes[1] == ind1){
            p2 = storage->points[nodes[0]];
            ind2 = nodes[0];
         }
      }

      std::vector<size_t> tmpnodes(4,0);

      if(visu::test_triangle2D_direct(p0,p1,p3)){
         tmpnodes = {ind0 + 1, ind1 + 1, ind2 + 1, ind3 + 1};
      }
      else{
         tmpnodes = {ind0 + 1, ind3 + 1, ind2 + 1, ind1 + 1};
      }
      Quadrangle quad(tmpnodes, nb_surface, 0, 0);
      msh.addQuadrangle(quad);
   }
   assert(storage->surfaces.size() == msh.getQuadrangles().size());

   return msh;

};



template<typename T>
Gmesh convertMesh(const disk::cartesian_mesh<T,3>& mesh)
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
   size_t nb_quads(0);
   for(auto surface :  storage->surfaces){
      nb_quads +=1;
      msh.addQuadrangle(convertQuadrangle(surface, nb_quads));
   }
   assert(storage->surfaces.size() == msh.getQuadrangles().size());


   // conversion tetra
   size_t nb_hexas(0);
   for(auto volume :  storage->volumes){
      nb_hexas +=1;
      msh.addHexahedron(convertHexahedron(volume, nb_hexas));
   }
   assert(storage->volumes.size() == msh.getHexahedra().size());

   assert((storage->edges.size() + storage->surfaces.size() + storage->nodes.size()
         + storage->volumes.size()) == msh.getNumberofElements());

   return msh;

};


} //visu

#endif
