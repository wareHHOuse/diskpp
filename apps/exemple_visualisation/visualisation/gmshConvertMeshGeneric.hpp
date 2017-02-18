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

#ifndef _GmshConvertMeshGeneric_HPP
#define _GmshConvertMeshGeneric_HPP

#include<vector>
#include<string>
#include<utility>
#include <array>
#include <cassert>

#include "loaders/loader.hpp"
#include "geometry/geometry.hpp"
#include "core/mesh/mesh.hpp"
#include "mesh/point.hpp"
#include "gmshMesh.h"
#include "gmshElement.h"
#include "gmshDisk.h"

namespace visu{


template<typename T>
Gmesh convertMesh(const disk::generic_mesh<T,1>& mesh)
{

   Gmesh msh;

   auto storage = mesh.backend_storage();


   // conversion point in node
   size_t nb_node(0);
   for(auto& point :  storage->points){
      nb_node +=1;
      msh.addNode(convertPoint(point, nb_node));

      Vertice vert(nb_node, nb_node, 0, 0);
      msh.addVertice(vert);
   }
   assert(storage->points.size() == msh.getNumberofNodes());


   // conversion edge
   size_t nb_edge(0);
   for(auto& edge :  storage->edges){
      nb_edge +=1;
      msh.addEdge(convertEdge(mesh, edge, nb_edge));
   }
   assert(storage->edges.size() == msh.getEdges().size());

   return msh;
};


template<typename T>
Gmesh convertMesh(const disk::generic_mesh<T,2>& mesh)
{

   Gmesh msh;

   auto storage = mesh.backend_storage();


   // conversion point in node
   size_t nb_node(0);
   for(auto& point :  storage->points){
      nb_node +=1;
      msh.addNode(convertPoint(point, nb_node));

      Vertice vert(nb_node, nb_node, 0, 0);
      msh.addVertice(vert);
   }
   assert(storage->points.size() == msh.getNumberofNodes());


   // conversion edge
   size_t nb_edge(0);
   for(auto& edge :  storage->edges){
      nb_edge +=1;
      msh.addEdge(convertEdge(edge, nb_edge));
   }
   assert(storage->edges.size() == msh.getEdges().size());

   //conversion surfaces
   size_t nb_surface(0);
   for(auto& cl : mesh ){
      auto fcs = disk::faces(mesh, cl);
      if(fcs.size() == 3){
         nb_surface += 1;
         auto cell_id = mesh.lookup(cl);
         auto surface = storage->surfaces[cell_id];
         msh.addTriangle(convertTriangle(mesh, surface, nb_surface));
      }
      if(fcs.size() == 4){
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
      else{
         auto face_id = mesh.lookup(fcs[0]);
         auto ptref = (storage->edges[face_id].point_ids())[0];
         auto pref = storage->points[ptref];

         for(size_t i = 1; i < fcs.size(); i++){
            auto nodes = storage->edges[mesh.lookup(fcs[i])].point_ids();
            if( nodes[0] != ptref && nodes[1] != ptref){
               nb_surface += 1;
               auto p0 = storage->points[nodes[0]];
               auto p1 = storage->points[nodes[1]];

               std::vector<size_t> tmpnodes(3,0);

               if(visu::test_triangle2D_direct(pref,p0,p1)){
                  tmpnodes = { ptref + 1, nodes[0] + 1, nodes[1] + 1 };
               }
               else{
                  tmpnodes = { ptref + 1, nodes[1] + 1, nodes[0] + 1 };
               }

               Triangle tri(tmpnodes, nb_surface, 0, 0);
               msh.addTriangle(tri);
            }
         }
      }
   }

   return msh;
};

template<typename T>
Gmesh convertMesh(const disk::generic_mesh<T,3>& mesh)
{
   Gmesh msh;

   auto storage = mesh.backend_storage();


   // conversion point in node
   size_t nb_node(0);
   for(auto& point :  storage->points){
      nb_node +=1;
      msh.addNode(convertPoint(point, nb_node));

      Vertice vert(nb_node, nb_node, 0, 0);
      msh.addVertice(vert);
   }
   assert(storage->points.size() == msh.getNumberofNodes());


   // conversion edge
   size_t nb_edge(0);
   for(auto& edge :  storage->edges){
      nb_edge +=1;
      msh.addEdge(convertEdge(edge, nb_edge));
   }
   assert(storage->edges.size() == msh.getEdges().size());



   return msh;

};


} //visu

#endif
