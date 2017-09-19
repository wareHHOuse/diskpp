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

#ifndef GMSHMESH_H
#define GMSHMESH_H

#include<vector>
#include<array>
#include<string>
#include "gmshElement.hpp"

namespace visu{

class Gmesh
{
   protected:
      size_t m_dim_topology;
      std::vector<Node> m_nodes;
      std::vector<Vertice> m_vertices;
      std::vector<Edge> m_edges;
      std::vector<Triangle> m_triangles;
      std::vector<Quadrangle> m_quadrangles;
      std::vector<Tetrahedron> m_tetrahedra;
      std::vector<Hexahedron> m_hexahedra;
      std::vector<Prism> m_prisms;
      std::vector<Pyramid> m_pyramids;
      size_t m_number_of_elements;
      bool m_verbose;

// private functions
      void readGmesh_MEDITformat(const std::string name_mesh);
      void writeGmesh_MEDITformat(const std::string name_mesh) const;
      void readGmesh_MSHformat(const std::string name_mesh);
      void writeGmesh_MSHformat(const std::string name_mesh) const;

      void convertInDiscontinuousMesh1();
      void convertInDiscontinuousMesh2();
      void convertInDiscontinuousMesh3();


   public:
      Gmesh();
      Gmesh(const size_t dim);

      Gmesh(const size_t dim, const bool verbose);

      Gmesh(const size_t dim, const std::vector<Node>& nodes, const std::vector<Vertice>& vertices,
            const std::vector<Edge>& edges);

      Gmesh(const size_t dim, const std::vector<Node>& nodes, const std::vector<Vertice>& vertices,
            const std::vector<Edge>& edges, const std::vector<Triangle>& triangles,
            const std::vector<Quadrangle>& quadrangles);

      Gmesh(const size_t dim, const std::vector<Node>& nodes, const std::vector<Vertice>& vertices,
            const std::vector<Edge>& edges, const std::vector<Triangle>& triangles,
            const std::vector<Quadrangle>& quadrangles, const std::vector<Tetrahedron>& tetrahedra,
            const std::vector<Hexahedron> hexahedra, const std::vector<Prism>& prisms,
            const std::vector<Pyramid>& pyramids);

      Gmesh(const size_t dim, const std::vector<Node>& nodes, const std::vector<Vertice>& vertices,
            const std::vector<Edge>& edges, const std::vector<Triangle>& triangles,
            const std::vector<Quadrangle>& quadrangles, const std::vector<Tetrahedron>& tetrahedra,
            const std::vector<Hexahedron> hexahedra);

      Gmesh(const size_t dim, const std::vector<Node>& nodes, const std::vector<Vertice>& vertices,
            const std::vector<Edge>& edges, const std::vector<Triangle>& triangles,
            const std::vector<Tetrahedron>& tetrahedra);

      Gmesh(const size_t dim, const std::vector<Node>& nodes, const std::vector<Vertice>& vertices,
            const std::vector<Edge>& edges, const std::vector<Quadrangle>& quadrangles,
            const std::vector<Hexahedron>& hexahedra);

      void readGmesh(const std::string name_mesh);
      void writeGmesh(const std::string name_mesh, const size_t format) const;
      void getInfo() const;
      size_t getDim() const;
      size_t getNumberofNodes() const;
      size_t getNumberofElements() const;
      bool verbose() const;
      void verbose(const bool verb);
      std::vector<Node> getNodes() const;
      std::vector<Vertice> getVertices() const;
      std::vector<Edge> getEdges() const;
      std::vector<Triangle> getTriangles() const;
      std::vector<Quadrangle> getQuadrangles() const;
      std::vector<Hexahedron> getHexahedra() const;
      std::vector<Tetrahedron> getTetrahedra() const;
      std::vector<Prism> getPrisms() const;
      std::vector<Pyramid> getPyramids() const;

      Node getNode(const size_t index) const;
      Vertice getVertice(const size_t index) const;
      Edge getEdge(const size_t index) const;
      Triangle getTriangle(const size_t index) const;
      Quadrangle getQuadrangle(const size_t index) const;
      Hexahedron getHexahedron(const size_t index) const;
      Tetrahedron getTetrahedron(const size_t index) const;
      Prism getPrism(const size_t index) const;
      Pyramid getPyramid(const size_t index) const;

      void addNode(const Node& node);
      void addVertice(const Vertice& vertice);
      void addEdge(const Edge& edge);
      void addTriangle(const Triangle& triangle);
      void addQuadrangle(const Quadrangle& quad);
      void addHexahedron(const Hexahedron& hexa);
      void addTetrahedron(const Tetrahedron& tetra);
      void addPrism(const Prism& prism);
      void addPyramid(const Pyramid& pyramid);

      void convertInDiscontinuousMesh();
};

} //visu

#endif
