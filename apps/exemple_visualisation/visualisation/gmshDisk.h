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

#ifndef GMSHDISK_H
#define GMSHDISK_H

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

namespace visu{

template< typename T>
bool test_triangle2D_direct(const T& p0, const T& p1, const T& p2)
{
   if(((p1.x()-p0.x())*(p2.y()-p0.y())-(p1.y()-p0.y())*(p2.x()-p0.x())) >= 0.0){
      return true;
   }
   else{
      return false;
   }
}


template<typename T, size_t DIM>
Node convertPoint( const point<T,DIM>& point, const size_t index)
{
   std::vector<double> coor(3,0.0);

   for(size_t i = 0; i < DIM; i++){
      coor[i] = double(point.at(i));
   }

   Node node(coor, index, 0);

   return node;
}

template< typename T>
Edge convertEdge( const T& edge, const size_t index)
{
   auto nodes = edge.point_ids();

   assert(nodes.size() == 2);

   std::vector<size_t> tmpnodes = { nodes[0] + 1, nodes[1] + 1 };

   Edge tmpedge(tmpnodes, index, 0, 0);

   return tmpedge;
}

template< typename T, typename Mesh>
Edge convertEdge(const Mesh& mesh, const T& edge, const size_t index)
{
   auto nodes = edge.point_ids();

   assert(nodes.size() == 2);

   auto storage = mesh.backend_storage();

   auto p0 = storage->points[nodes[0]];
   auto p1 = storage->points[nodes[1]];
   std::vector<size_t> tmpnodes(2, 0);

   if(p0.x() <= p1.x()){
       tmpnodes = { nodes[0] + 1, nodes[1] + 1 };
   }
   else {
       tmpnodes = { nodes[1] + 1, nodes[0] + 1 };
   }

   Edge tmpedge(tmpnodes, index, 0, 0);

   return tmpedge;
}

template< typename T>
Triangle convertTriangle( const T& triangle, const size_t index)
{
   auto nodes = triangle.point_ids();

   assert(nodes.size() == 3);

//    std::cout << " ------------------------------------------------ " << std::endl;
//    std::cout << nodes[0] + 1 << " " << nodes[1] + 1 << " "<< nodes[2] + 1 << " " << std::endl;

   std::vector<size_t> tmpnodes = { nodes[0] + 1, nodes[1] + 1, nodes[2] + 1 };

   Triangle tri(tmpnodes, index, 0, 0);

   return tri;
}

template< typename T, typename Mesh>
Triangle convertTriangle(const Mesh& mesh, const T& triangle, const size_t index)
{
   auto nodes = triangle.point_ids();

   assert(nodes.size() == 3);

   auto storage = mesh.backend_storage();

   auto p0 = storage->points[nodes[0]];
   auto p1 = storage->points[nodes[1]];
   auto p2 = storage->points[nodes[2]];
   std::vector<size_t> tmpnodes(3,0);

   if(test_triangle2D_direct(p0,p1,p2)){
      tmpnodes = { nodes[0] + 1, nodes[1] + 1, nodes[2] + 1 };
   }
   else{
      tmpnodes = { nodes[0] + 1, nodes[2] + 1, nodes[1] + 1 };
   }

//    std::cout << " ------------------------------------------------ " << std::endl;
//    std::cout << nodes[0] + 1 << " " << nodes[1] + 1 << " "<< nodes[2] + 1 << " " << std::endl;

   Triangle tri(tmpnodes, index, 0, 0);

   return tri;
}

template< typename T>
Quadrangle convertQuadrangle( const T& quadrangle, const size_t index)
{
   auto nodes = quadrangle.point_ids();

   assert(nodes.size() == 4);

   std::vector<size_t> tmpnodes = { nodes[0] + 1, nodes[1] + 1, nodes[3] + 1, nodes[2] + 1 }; //the ordering of nodes a bit strange

   Quadrangle quad(tmpnodes, index, 0, 0);

   return quad;
}

template< typename T>
Tetrahedron convertTetrahedron( const T& tetrahedron, const size_t index)
{
   auto nodes = tetrahedron.point_ids();

   assert(nodes.size() == 4);
   /*
   std::cout << " ------------------------------------------------ " << std::endl;
   std::cout << nodes[0] + 1 << " " << nodes[1] + 1 << " "<< nodes[2] + 1 << " "<< nodes[3] + 1 << " " << std::endl;*/

   std::vector<size_t> tmpnodes = { nodes[0] + 1, nodes[1] + 1, nodes[2] + 1, nodes[3] + 1 };

   Tetrahedron tetra(tmpnodes, index, 0, 0);

   return tetra;
}

template< typename T, typename Mesh>
Tetrahedron convertTetrahedron(const Mesh& mesh, const T& tetrahedron, const size_t index)
{
   auto nodes = tetrahedron.point_ids();

   assert(nodes.size() == 4);
   /*
   std::cout << " ------------------------------------------------ " << std::endl;
   std::cout << nodes[0] + 1 << " " << nodes[1] + 1 << " "<< nodes[2] + 1 << " "<< nodes[3] + 1 << " " << std::endl;*/
   auto storage = mesh.backend_storage();
   auto p0 = storage->points[nodes[0]];
   auto p1 = storage->points[nodes[1]];
   auto p2 = storage->points[nodes[2]];
   auto p3 = storage->points[nodes[3]];

   auto v1 = p1 - p0;
   auto v2 = p2 - p0;
   auto v3 = p3 - p0;

   auto det = v1.x()*v2.y()*v3.z() + v1.z()*v2.x()*v3.y() + v1.y()*v2.z()*v3.x()
            - v1.z()*v2.y()*v3.x() - v1.y()*v2.x()*v3.z() - v1.x()*v2.z()*v3.y();

   std::vector<size_t> tmpnodes(4, 0);

   // det>0 if the tetrahedra is direct
   if(det >= 0.0){
      tmpnodes = { nodes[0] + 1, nodes[1] + 1, nodes[2] + 1, nodes[3] + 1 };
   }
   else {
      tmpnodes = { nodes[0] + 1, nodes[2] + 1, nodes[1] + 1, nodes[3] + 1 };
   }

   Tetrahedron tetra(tmpnodes, index, 0, 0);

   return tetra;
}

template< typename T>
Hexahedron convertHexahedron( const T& hexahedron, const size_t index)
{
   auto nodes = hexahedron.point_ids();

   assert(nodes.size() == 8);

   /*std::cout << " ------------------------------------------------ " << std::endl;
   std::cout << nodes[0] + 1 << " " << nodes[1] + 1 << " "<< nodes[3] + 1 << " "<< nodes[2] + 1 << " " << std::endl;
   std::cout << nodes[4] + 1 << " " << nodes[5] + 1 << " "<< nodes[7] + 1 << " "<< nodes[6] + 1 << " " << std::endl;*/

   std::vector<size_t> tmpnodes = { nodes[0] + 1, nodes[1] + 1, nodes[3] + 1, nodes[2] + 1,
                                    nodes[4] + 1, nodes[5] + 1, nodes[7] + 1, nodes[6] + 1 }; //the ordering of nodes a bit strange

   Hexahedron hexa(tmpnodes, index, 0, 0);

   return hexa;
}


template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage, size_t DIM, typename Cell>
auto cell_nodes(const Mesh<T, DIM, Storage>& mesh, const Cell& cl)
{
   typedef Mesh<T, DIM, Storage>                 mesh_type;
   static_assert(sizeof(mesh_type) == -1, "convert: not suitable for the requested kind of mesh");
}

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage, typename Cell>
auto cell_nodes(const Mesh<T, 1, Storage>& mesh, const Cell& cl)
{
   auto storage = mesh.backend_storage();
   auto cell_id = mesh.lookup(cl);
   return (storage->edges[cell_id]).point_ids();
}

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage, typename Cell>
auto cell_nodes(const Mesh<T, 2, Storage>& mesh, const Cell& cl)
{
   auto storage = mesh.backend_storage();
   auto cell_id = mesh.lookup(cl);
   return (storage->surfaces[cell_id]).point_ids();
}

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage, typename Cell>
auto cell_nodes(const Mesh<T, 3, Storage>& mesh, const Cell& cl)
{
   auto storage = mesh.backend_storage();
   auto cell_id = mesh.lookup(cl);
   return (storage->volumes[cell_id]).point_ids();
}


template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage, size_t DIM, typename Face>
auto face_nodes(const Mesh<T, DIM, Storage>& mesh, const Face& fc)
{
   typedef Mesh<T, DIM, Storage>                 mesh_type;
   static_assert(sizeof(mesh_type) == -1, "convert: not suitable for the requested kind of mesh");
}

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage, typename Face>
auto face_nodes(const Mesh<T, 1, Storage>& mesh, const Face& fc)
{
   auto storage = mesh.backend_storage();
   auto face_id = mesh.lookup(fc);
   return (storage->nodes[face_id]).point_ids();
}

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage, typename Face>
auto face_nodes(const Mesh<T, 2, Storage>& mesh, const Face& fc)
{
   auto storage = mesh.backend_storage();
   auto face_id = mesh.lookup(fc);
   return (storage->edges[face_id]).point_ids();
}

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage, typename Face>
auto face_nodes(const Mesh<T, 3, Storage>& mesh, const Face& fc)
{
   auto storage = mesh.backend_storage();
   auto face_id = mesh.lookup(fc);
   return (storage->surfaces[face_id]).point_ids();
}

} //visu

#endif
