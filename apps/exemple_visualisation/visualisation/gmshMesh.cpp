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


#include "gmshMesh.hpp"
#include "gmshElement.hpp"
#include<iostream>
#include<fstream>
#include <assert.h>

namespace visu{

//  Class Gmesh

Gmesh::Gmesh() : m_dim_topology(0), m_nodes(), m_vertices(), m_edges(), m_triangles(), m_quadrangles(),
                 m_tetrahedra() , m_hexahedra(), m_prisms(), m_pyramids(), m_number_of_elements(0), m_verbose(true) {}

Gmesh::Gmesh(const size_t dim) : m_dim_topology(dim), m_nodes(), m_vertices(), m_edges(), m_triangles(), m_quadrangles(),
                 m_tetrahedra() , m_hexahedra(), m_prisms(), m_pyramids(), m_number_of_elements(0), m_verbose(true) {}

Gmesh::Gmesh(const size_t dim, const bool verbose) : m_dim_topology(dim), m_nodes(), m_vertices(), m_edges(), m_triangles(), m_quadrangles(),
                 m_tetrahedra() , m_hexahedra(), m_prisms(), m_pyramids(), m_number_of_elements(0), m_verbose(verbose) {}

Gmesh::Gmesh(const size_t dim, const std::vector<Node>& nodes, const std::vector<Vertice>& vertices,
      const std::vector<Edge>& edges) :
      m_dim_topology(dim), m_nodes(nodes), m_vertices(vertices), m_edges(edges), m_triangles(), m_quadrangles(),
      m_tetrahedra() , m_hexahedra(), m_prisms(), m_pyramids(),
      m_number_of_elements(vertices.size() + edges.size()), m_verbose(true) {}

Gmesh::Gmesh(const size_t dim, const std::vector<Node>& nodes, const std::vector<Vertice>& vertices,
      const std::vector<Edge>& edges, const std::vector<Triangle>& triangles,
      const std::vector<Quadrangle>& quadrangles) :
      m_dim_topology(dim), m_nodes(nodes), m_vertices(vertices), m_edges(edges),
      m_triangles(triangles), m_quadrangles(quadrangles),
      m_tetrahedra() , m_hexahedra(), m_prisms(), m_pyramids(),
      m_number_of_elements(vertices.size() + edges.size() + triangles.size() + quadrangles.size()) {}

Gmesh::Gmesh(const size_t dim, const std::vector<Node>& nodes, const std::vector<Vertice>& vertices,
      const std::vector<Edge>& edges, const std::vector<Triangle>& triangles,
      const std::vector<Quadrangle>& quadrangles, const std::vector<Tetrahedron>& tetrahedra,
      const std::vector<Hexahedron> hexahedra) :
            m_dim_topology(dim), m_nodes(nodes), m_vertices(vertices), m_edges(edges),
            m_triangles(triangles), m_quadrangles(quadrangles),
            m_tetrahedra(tetrahedra) , m_hexahedra(hexahedra), m_prisms(), m_pyramids(),
            m_number_of_elements(vertices.size() + edges.size() + triangles.size() + quadrangles.size()
            + tetrahedra.size() + hexahedra.size()), m_verbose(true) {}

Gmesh::Gmesh(const size_t dim, const std::vector<Node>& nodes, const std::vector<Vertice>& vertices,
      const std::vector<Edge>& edges, const std::vector<Triangle>& triangles,
      const std::vector<Quadrangle>& quadrangles, const std::vector<Tetrahedron>& tetrahedra,
      const std::vector<Hexahedron> hexahedra, const std::vector<Prism>& prisms,
      const std::vector<Pyramid>& pyramids) :
      m_dim_topology(dim), m_nodes(nodes), m_vertices(vertices), m_edges(edges),
      m_triangles(triangles), m_quadrangles(quadrangles),
      m_tetrahedra(tetrahedra) , m_hexahedra(hexahedra), m_prisms(prisms), m_pyramids(pyramids),
      m_number_of_elements(vertices.size() + edges.size() + triangles.size() + quadrangles.size()
      + tetrahedra.size() + hexahedra.size() + prisms.size() + pyramids.size()), m_verbose(true) {}

Gmesh::Gmesh(const size_t dim, const std::vector<Node>& nodes, const std::vector<Vertice>& vertices,
      const std::vector<Edge>& edges, const std::vector<Triangle>& triangles,
      const std::vector<Tetrahedron>& tetrahedra) :
            m_dim_topology(dim), m_nodes(nodes), m_vertices(vertices), m_edges(edges),
            m_triangles(triangles), m_quadrangles(),
            m_tetrahedra(tetrahedra) , m_hexahedra(), m_prisms(), m_pyramids(),
            m_number_of_elements(vertices.size() + edges.size() + triangles.size()
            + tetrahedra.size()), m_verbose(true) {}

Gmesh::Gmesh(const size_t dim, const std::vector<Node>& nodes, const std::vector<Vertice>& vertices,
            const std::vector<Edge>& edges, const std::vector<Quadrangle>& quadrangles,
            const std::vector<Hexahedron>& hexahedra) :
               m_dim_topology(dim), m_nodes(nodes), m_vertices(vertices), m_edges(edges),
               m_triangles(), m_quadrangles(quadrangles),
               m_tetrahedra() , m_hexahedra(hexahedra), m_prisms(), m_pyramids(),
               m_number_of_elements(vertices.size() + edges.size() + quadrangles.size()
               + hexahedra.size()), m_verbose(true) {}

bool Gmesh::verbose() const
{
   return m_verbose;
}
void Gmesh::verbose(const bool verb)
{
   m_verbose = verb;
}
size_t Gmesh::getNumberofNodes() const
{
   return m_nodes.size();
}

size_t Gmesh::getNumberofElements() const
{
   return m_number_of_elements;
}

size_t Gmesh::getDim() const
{
   return m_dim_topology;
}

std::vector<Node> Gmesh::getNodes() const
{
   return m_nodes;
}

std::vector<Vertice> Gmesh::getVertices() const
{
   return m_vertices;
}

std::vector<Edge> Gmesh::getEdges() const
{
   return m_edges;
}

std::vector<Triangle> Gmesh::getTriangles() const
{
   return m_triangles;
}

std::vector<Quadrangle> Gmesh::getQuadrangles() const
{
   return m_quadrangles;
}

std::vector<Tetrahedron> Gmesh::getTetrahedra() const
{
   return m_tetrahedra;
}

std::vector<Hexahedron> Gmesh::getHexahedra() const
{
   return m_hexahedra;
}

std::vector<Prism> Gmesh::getPrisms() const
{
   return m_prisms;
}

std::vector<Pyramid> Gmesh::getPyramids() const
{
   return m_pyramids;
}

Node Gmesh::getNode(const size_t index) const
{
   return m_nodes.at(index);
}

Vertice Gmesh::getVertice(const size_t index) const
{
   return m_vertices.at(index);
}

Edge Gmesh::getEdge(const size_t index) const
{
   return m_edges.at(index);
}

Triangle Gmesh::getTriangle(const size_t index) const
{
   return m_triangles.at(index);
}

Quadrangle Gmesh::getQuadrangle(const size_t index) const
{
   return m_quadrangles.at(index);
}

Hexahedron Gmesh::getHexahedron(const size_t index) const
{
   return m_hexahedra.at(index);
}

Tetrahedron Gmesh::getTetrahedron(const size_t index) const
{
   return m_tetrahedra.at(index);
}

Prism Gmesh::getPrism(const size_t index) const
{
   return m_prisms.at(index);
}

Pyramid Gmesh::getPyramid(const size_t index) const
{
   return m_pyramids.at(index);
}

void Gmesh::readGmesh_MEDITformat(const std::string name_mesh)
{
   std::ifstream mesh_file(name_mesh.data());
   if (!mesh_file.is_open())
   {
      std::cout << "Unable to open file " << name_mesh << std::endl;
      abort();
   }

   if(m_verbose){
      std::cout << "------------------------" << std::endl;
      std::cout << "Reading mesh " << name_mesh << std::endl;
      std::cout << "MEDIT Format" << std::endl;
      std::cout << "------------------------" << std::endl;
   }

   std::string file_line;
   size_t nb_nodes(0), nb_edges(0), nb_triangles(0), nb_tetras(0), nb_quad(0), nb_hexa(0);
   std::array<double, 3> coor;
   std::vector<size_t> edg(2), tri(3), tet(4), quad(4), hexa(8);
   size_t ref(0), loopm_nodes(100), offset_elements(0);

   while (!mesh_file.eof())
   {
      getline(mesh_file, file_line);
      if ((file_line.find("Vertices") != std::string::npos)&&(loopm_nodes))
      {
         mesh_file >> nb_nodes;
         if(m_verbose)
            std::cout << "Reading vertices (" << nb_nodes << ")" << std::endl;
         for (size_t i = 1 ; i <= nb_nodes ; i++)
         {
            mesh_file >> coor[0] >> coor[1] >> coor[2] >> ref;
            Node node(coor, i, ref);
            m_nodes.push_back(node);
         }
         assert((m_nodes.size() == nb_nodes) && "Problem 1 in mesh reading.");
         loopm_nodes = 0;
      }
      if (file_line.find("Edges") != std::string::npos)
      {
         mesh_file >> nb_edges;
         if(m_verbose)
            std::cout << "Reading edges (" << nb_edges << ")" << std::endl;
         for (size_t i = 1 ; i <= nb_edges ; i++)
         {
            mesh_file >> edg[0] >> edg[1] >> ref;
            Edge edge(edg, offset_elements + i, ref, 0);
            m_edges.push_back(edge);
         }

         offset_elements += nb_edges;
         assert((m_edges.size() == nb_edges) && "Problem 2 in mesh reading.");
      }
      if (file_line.find("Triangles") != std::string::npos)
      {
         mesh_file >> nb_triangles;
         if(m_verbose)
            std::cout << "Reading triangles (" << nb_triangles << ")" << std::endl;
         for (size_t i = 1 ; i <= nb_triangles ; i++)
         {
            mesh_file >> tri[0] >> tri[1] >> tri[2] >> ref;
            Triangle triangle(tri, offset_elements + i, ref, 0);
            m_triangles.push_back(triangle);
         }

         offset_elements += nb_triangles;
         assert((m_triangles.size() == nb_triangles) && "Problem 3 in mesh reading.");
      }
      if (file_line.find("Tetrahedra") != std::string::npos)
      {
         mesh_file >> nb_tetras;
         if(m_verbose)
            std::cout << "Reading tetrahedra (" << nb_tetras << ")" << std::endl;
         for (size_t i = 0 ; i < nb_tetras ; i++)
         {
            mesh_file >> tet[0] >> tet[1] >> tet[2] >> tet[3] >> ref;
            Tetrahedron tetrahedron(tet, offset_elements + i, ref, 0);
            m_tetrahedra.push_back(tetrahedron);
         }

         offset_elements += nb_tetras;
         assert((m_tetrahedra.size() == nb_tetras) && "Problem 4 in mesh reading.");
      }
      if (file_line.find("Quadrilaterals") != std::string::npos)
      {
         mesh_file >> nb_quad;
         if(m_verbose)
            std::cout << "Reading quadrilaterals (" << nb_quad << ")" << std::endl;
         for (size_t i = 1 ; i <= nb_quad ; i++)
         {
            mesh_file >> quad[0] >> quad[1] >> quad[2] >> quad[3]  >> ref;
            Quadrangle quadrangle(quad, offset_elements + i, ref, 0);
            m_quadrangles.push_back(quadrangle);
         }

         offset_elements += nb_quad;
         assert((m_quadrangles.size() == nb_quad) && "Problem 3 in mesh reading.");
      }
      if (file_line.find("Hexahedra") != std::string::npos)
      {
         mesh_file >> nb_hexa;
         if(m_verbose)
            std::cout << "Reading hexahedra (" << nb_hexa << ")" << std::endl;
         for (size_t i = 0 ; i < nb_hexa ; i++)
         {
            mesh_file >> hexa[0] >> hexa[1] >> hexa[2] >> hexa[3] >> hexa[4]
            >> hexa[5] >> hexa[6] >> hexa[7] >> hexa[8] >> ref;
            Hexahedron hexahedron(hexa, offset_elements + i, ref, 0);
            m_hexahedra.push_back(hexahedron);
         }

         offset_elements += nb_hexa;
         assert((m_hexahedra.size() == nb_hexa) && "Problem 4 in mesh reading.");
      }
   }

   assert(offset_elements == (nb_edges + nb_triangles + nb_quad + nb_hexa + nb_tetras));
   m_number_of_elements = offset_elements;

   mesh_file.close();

   if((m_tetrahedra.size() + m_hexahedra.size() + m_prisms.size() + m_pyramids.size()) > 0){
      m_dim_topology = 3;
   }
   else if((m_triangles.size() + m_quadrangles.size()) > 0){
      m_dim_topology = 2;
   }
   else if(m_edges.size() > 0){
      m_dim_topology = 1;
   }

   std::cout << "------------------------" << std::endl;
}


void Gmesh::readGmesh_MSHformat(const std::string name_mesh)
{
   std::ifstream mesh_file(name_mesh.data());
   if (!mesh_file.is_open())
   {
      std::cout << "Unable to open file " << name_mesh << std::endl;
      abort();
   }

   if(m_verbose){
      std::cout << "------------------------" << std::endl;
      std::cout << "Reading mesh " << name_mesh << std::endl;
      std::cout << "MSH Format"  << std::endl;
      std::cout << "------------------------" << std::endl;
   }


   std::string file_line;
   size_t nb_nodes(0), nb_elements(0);
   std::array<double, 3> coor;
   size_t index(0), type_elem(0), nbinteger_tag(0), loop_nodes(100), loop_elem(100);

   while (!mesh_file.eof() && loop_elem)
   {
      getline(mesh_file, file_line);
      if ((file_line.find("$Nodes") != std::string::npos)&&(loop_nodes))
      {
         mesh_file >> nb_nodes;
         if(m_verbose)
            std::cout << "Reading nodes (" << nb_nodes << ")" << std::endl;
         for (size_t i = 0 ; i < nb_nodes ; ++i)
         {
            mesh_file >> index >> coor[0] >> coor[1] >> coor[2];
            Node node(coor, index, 1);
            m_nodes.push_back(node);
         }
         assert((m_nodes.size() == nb_nodes) && "Problem 1 in mesh reading.");
         loop_nodes = 0;
      }
      if ((file_line.find("$Elements") != std::string::npos) && (loop_elem))
      {
         mesh_file >> nb_elements;
         m_number_of_elements = nb_elements;
         for (size_t i = 0 ; i < nb_elements ; ++i)
         {
            mesh_file >> index >> type_elem >> nbinteger_tag;

            std::vector<size_t> integer_tag(nbinteger_tag);
            for(size_t jtag=0; jtag < nbinteger_tag; jtag++){
               mesh_file >> integer_tag[jtag];
            }

            std::vector<size_t> elem_nodes;
            size_t ind_node = 0;
            if(type_elem == 1){
               for(size_t inode = 0; inode < 2; inode++){
                  mesh_file >> ind_node;
                  elem_nodes.push_back(ind_node);
               }
               Edge edge(elem_nodes, i+1, integer_tag[0], integer_tag[1]);
               m_edges.push_back(edge);
            }
            else if(type_elem == 2){
               for(size_t inode = 0; inode < 3; inode++){
                  mesh_file >> ind_node;
                  elem_nodes.push_back(ind_node);
               }
               Triangle tri(elem_nodes, i+1, integer_tag[0], integer_tag[1]);
               m_triangles.push_back(tri);
            }
            else if(type_elem == 3){
               for(size_t inode = 0; inode < 4; inode++){
                  mesh_file >> ind_node;
                  elem_nodes.push_back(ind_node);
               }
               Quadrangle quad(elem_nodes, i+1, integer_tag[0], integer_tag[1]);
               m_quadrangles.push_back(quad);
            }
            else if(type_elem == 4){
               for(size_t inode = 0; inode < 4; inode++){
                  mesh_file >> ind_node;
                  elem_nodes.push_back(ind_node);
               }
               Tetrahedron tetra(elem_nodes, i+1, integer_tag[0], integer_tag[1]);
               m_tetrahedra.push_back(tetra);
            }
            else if(type_elem == 5){
               for(size_t inode = 0; inode < 8; inode++){
                  mesh_file >> ind_node;
                  elem_nodes.push_back(ind_node);
               }
               Hexahedron hexa(elem_nodes, i+1, integer_tag[0], integer_tag[1]);
               m_hexahedra.push_back(hexa);
            }
            else if(type_elem == 6){
               for(size_t inode = 0; inode < 6; inode++){
                  mesh_file >> ind_node;
                  elem_nodes.push_back(ind_node);
               }
               Prism prism(elem_nodes, i+1, integer_tag[0], integer_tag[1]);
               m_prisms.push_back(prism);
            }
            else if(type_elem == 7){
               for(size_t inode = 0; inode < 5; inode++){
                  mesh_file >> ind_node;
                  elem_nodes.push_back(ind_node);
               }
               Pyramid pyra(elem_nodes, i+1, integer_tag[0], integer_tag[1]);
               m_pyramids.push_back(pyra);
            }
            else {
               throw std::invalid_argument("Unsupported Element");
            }
         }
         loop_elem = 0;
      }
   }

   mesh_file.close();

   if((m_tetrahedra.size() + m_hexahedra.size() + m_prisms.size() + m_pyramids.size()) > 0){
      m_dim_topology = 3;
   }
   else if((m_triangles.size() + m_quadrangles.size()) > 0){
      m_dim_topology = 2;
   }
   else if(m_edges.size() > 0){
      m_dim_topology = 1;
   }

   std::cout << "------------------------" << std::endl;
}


void Gmesh::readGmesh(const std::string name_mesh)
{
   size_t pos= name_mesh.find(".m");
   std::string ext = name_mesh.substr(pos);

   if (ext == ".mesh") {
      Gmesh::readGmesh_MEDITformat(name_mesh);
   }
   else if (ext == ".msh") {
      Gmesh::readGmesh_MSHformat(name_mesh);
   }
   else {
      std::cout << "Unknown file format" << std::endl;
      abort();
   }

}

void Gmesh::writeGmesh_MEDITformat(const std::string name_mesh) const
{
   std::ofstream mesh_file(name_mesh);
   if (!mesh_file.is_open())
   {
      std::cout << "Unable to open file " << name_mesh << std::endl;
      abort();
   }
   else
   {
      if(m_verbose){
         std::cout << "------------------------" << std::endl;
         std::cout << "Writing mesh " << name_mesh << std::endl;
         std::cout << "MESH file format"  << std::endl;
         std::cout << "------------------------" << std::endl;
      }
   }

   mesh_file << "MeshVersionFormatted 2" << std::endl;
   mesh_file << "Dimension 3"  << std::endl;
   mesh_file << " " << std::endl;

   mesh_file << "Vertices" << std::endl;
   mesh_file << m_nodes.size() << std::endl;

   for (size_t i = 0; i < m_nodes.size(); i++){
      Node node(m_nodes[i]);
      mesh_file         << (node.getCoordinate())[0]
      << " "  << (node.getCoordinate())[1]
      << " "  << (node.getCoordinate())[2]
      << " "  << node.getRef()
      << std::endl;
   }

   mesh_file << " " << std::endl;
   mesh_file << "Edges" << std::endl;
   mesh_file << m_edges.size() << std::endl;

   for (size_t i = 0; i < m_edges.size(); i++){
      Edge edge(m_edges[i]);
      mesh_file           << (edge.getNodes())[0]
      << " "  << (edge.getNodes())[1]
      << " "  << edge.getPhysicalEntities()
      << std::endl;
   }

   mesh_file << " " << std::endl;
   mesh_file << "Triangles" << std::endl;
   mesh_file << m_triangles.size() << std::endl;

   for (size_t i = 0; i < m_triangles.size(); i++){
      Triangle tri(m_triangles[i]);
      mesh_file          << (tri.getNodes())[0]
      << " "  << (tri.getNodes())[1]
      << " "  << (tri.getNodes())[2]
      << " "  << tri.getPhysicalEntities()
      << std::endl;
   }

   mesh_file << " " << std::endl;
   mesh_file << "Quadrilaterals" << std::endl;
   mesh_file << m_quadrangles.size() << std::endl;

   for (size_t i = 0; i < m_quadrangles.size(); i++){
      Quadrangle quad(m_quadrangles[i]);
      mesh_file          << (quad.getNodes())[0]
      << " "  << (quad.getNodes())[1]
      << " "  << (quad.getNodes())[2]
      << " "  << (quad.getNodes())[3]
      << " "  << quad.getPhysicalEntities()
      << std::endl;
   }

   mesh_file << " " << std::endl;
   mesh_file << "Tetrahedra" << std::endl;
   mesh_file << m_tetrahedra.size() << std::endl;

   for (size_t i = 0; i < m_tetrahedra.size(); i++){
      Tetrahedron tetra(m_tetrahedra[i]);
      mesh_file          << (tetra.getNodes())[0]
      << " "  << (tetra.getNodes())[1]
      << " "  << (tetra.getNodes())[2]
      << " "  << (tetra.getNodes())[3]
      << " "  << tetra.getPhysicalEntities()
      << std::endl;
   }

   mesh_file << " " << std::endl;
   mesh_file << "Hexahedra" << std::endl;
   mesh_file << m_hexahedra.size() << std::endl;

   for (size_t i = 0; i < m_quadrangles.size(); i++){
      Hexahedron hexa(m_hexahedra[i]);
      mesh_file          << (hexa.getNodes())[0]
      << " "  << (hexa.getNodes())[1]
      << " "  << (hexa.getNodes())[2]
      << " "  << (hexa.getNodes())[3]
      << " "  << (hexa.getNodes())[4]
      << " "  << (hexa.getNodes())[5]
      << " "  << (hexa.getNodes())[6]
      << " "  << (hexa.getNodes())[7]
      << " "  << hexa.getPhysicalEntities()
      << std::endl;
   }

   mesh_file << " "    << std::endl;
   mesh_file << "End"  << std::endl;

   mesh_file.close();

   std::cout << "------------------------" << std::endl;

}

void Gmesh::writeGmesh_MSHformat(const std::string name_mesh) const
{
   std::ofstream mesh_file(name_mesh);
   if (!mesh_file.is_open())
   {
      std::cout << "Unable to open file " << name_mesh << std::endl;
      abort();
   }
   else
   {
      if(m_verbose){
         std::cout << "------------------------" << std::endl;
         std::cout << "Writing mesh " << name_mesh << std::endl;
         std::cout << "MSH ASCII file format"  << std::endl;
         std::cout << "------------------------" << std::endl;
      }
   }

   size_t offset_elements(1);

   mesh_file << "$MeshFormat" << std::endl;
   mesh_file << "2.2 0 8" << std::endl;
   mesh_file << "$EndMeshFormat" << std::endl;

   mesh_file << "$Nodes" << std::endl;
   mesh_file << m_nodes.size() << std::endl;

   for (size_t i = 0; i < m_nodes.size(); i++){
      Node node(m_nodes[i]);
      mesh_file << node.getIndex()     << " "  << (node.getCoordinate())[0]
      << " "  << (node.getCoordinate())[1]
      << " "  << (node.getCoordinate())[2]
      << std::endl;
   }

   mesh_file << "$EndNodes" << std::endl;
   mesh_file << "$Elements" << std::endl;
   mesh_file << m_number_of_elements << std::endl;

   for (size_t i = 0; i < m_vertices.size(); i++){
      Vertice vert(m_vertices[i]);
      mesh_file << offset_elements << " "  << 15
      << " "  << 2 //number of tag
      << " "  << vert.getPhysicalEntities()
      << " "  << vert.getElemTag()
      << " "  << (vert.getNodes())[0]
      << std::endl;
      offset_elements += 1;
   }

   for (size_t i = 0; i < m_edges.size(); i++){
      Edge edge(m_edges[i]);
      mesh_file << offset_elements << " "  << 1
      << " "  << 2 //number of tag
      << " "  << edge.getPhysicalEntities()
      << " "  << edge.getElemTag()
      << " "  << (edge.getNodes())[0]
      << " "  << (edge.getNodes())[1]
      << std::endl;
      offset_elements += 1;
   }

   for (size_t i = 0; i < m_triangles.size(); i++){
      Triangle tri(m_triangles[i]);
      mesh_file << offset_elements  << " "  << 2
      << " "  << 2 //number of tag
      << " "  << tri.getPhysicalEntities()
      << " "  << tri.getElemTag()
      << " "  << (tri.getNodes())[0]
      << " "  << (tri.getNodes())[1]
      << " "  << (tri.getNodes())[2]
      << std::endl;
      offset_elements += 1;
   }

   for (size_t i = 0; i < m_quadrangles.size(); i++){
      Quadrangle quad(m_quadrangles[i]);
      mesh_file << offset_elements  << " "  << 3
      << " "  << 2 //number of tag
      << " "  << quad.getPhysicalEntities()
      << " "  << quad.getElemTag()
      << " "  << (quad.getNodes())[0]
      << " "  << (quad.getNodes())[1]
      << " "  << (quad.getNodes())[2]
      << " "  << (quad.getNodes())[3]
      << std::endl;
      offset_elements += 1;
   }

   for (size_t i = 0; i < m_tetrahedra.size(); i++){
      Tetrahedron tet(m_tetrahedra[i]);
      mesh_file << offset_elements  << " "  << 4
      << " "  << 2 //number of tag
      << " "  << tet.getPhysicalEntities()
      << " "  << tet.getElemTag()
      << " "  << (tet.getNodes())[0]
      << " "  << (tet.getNodes())[1]
      << " "  << (tet.getNodes())[2]
      << " "  << (tet.getNodes())[3]
      << std::endl;
      offset_elements += 1;
   }

   for (size_t i = 0; i < m_hexahedra.size(); i++){
      Hexahedron hexa(m_hexahedra[i]);
      mesh_file << offset_elements  << " "  << 5
      << " "  << 2 //number of tag
      << " "  << hexa.getPhysicalEntities()
      << " "  << hexa.getElemTag()
      << " "  << (hexa.getNodes())[0]
      << " "  << (hexa.getNodes())[1]
      << " "  << (hexa.getNodes())[2]
      << " "  << (hexa.getNodes())[3]
      << " "  << (hexa.getNodes())[4]
      << " "  << (hexa.getNodes())[5]
      << " "  << (hexa.getNodes())[6]
      << " "  << (hexa.getNodes())[7]
      << std::endl;
      offset_elements += 1;
   }

   for (size_t i = 0; i < m_prisms.size(); i++){
      Prism pri(m_prisms[i]);
      mesh_file << offset_elements  << " "  << 6
      << " "  << 2 //number of tag
      << " "  << pri.getPhysicalEntities()
      << " "  << pri.getElemTag()
      << " "  << (pri.getNodes())[0]
      << " "  << (pri.getNodes())[1]
      << " "  << (pri.getNodes())[2]
      << " "  << (pri.getNodes())[3]
      << " "  << (pri.getNodes())[4]
      << " "  << (pri.getNodes())[5]
      << std::endl;
      offset_elements += 1;
   }

   for (size_t i = 0; i < m_pyramids.size(); i++){
      Pyramid py(m_pyramids[i]);
      mesh_file << offset_elements  << " "  << 7
      << " "  << 2 //number of tag
      << " "  << py.getPhysicalEntities()
      << " "  << py.getElemTag()
      << " "  << (py.getNodes())[0]
      << " "  << (py.getNodes())[1]
      << " "  << (py.getNodes())[2]
      << " "  << (py.getNodes())[3]
      << " "  << (py.getNodes())[4]
      << std::endl;
      offset_elements += 1;
   }

   mesh_file << "$EndElements" << std::endl;

   mesh_file.close();

   std::cout << "------------------------" << std::endl;
}

void Gmesh::writeGmesh(const std::string name_mesh, const size_t format) const
{
   if (format == 1) {
      Gmesh::writeGmesh_MEDITformat(name_mesh);
   }
   else if (format == 2) {
      Gmesh::writeGmesh_MSHformat(name_mesh);
   }
   else {
      std::cout << "Unknown file format" << std::endl;
      abort();
   }
}

void Gmesh::getInfo() const
{
   std::cout << "----------------------------" << std::endl;
   std::cout << "Informations about Gmsh mesh" << std::endl;
   std::cout << "----------------------------" << std::endl;
   std::cout << "Dimension : " << m_dim_topology << std::endl;
   std::cout << " * Numbers of vetices : " << m_nodes.size() << std::endl;
   std::cout << " * Numbers of edges : " << m_edges.size() << std::endl;
   std::cout << " * Numbers of triangles : " << m_triangles.size() << std::endl;
   std::cout << " * Numbers of quadrangles : " << m_quadrangles.size() << std::endl;
   std::cout << " * Numbers of tetrahedra : " << m_tetrahedra.size() << std::endl;
   std::cout << " * Numbers of hexahedra : " << m_hexahedra.size() << std::endl;
   std::cout << " * Numbers of prisms : " << m_prisms.size() << std::endl;
   std::cout << " * Numbers of pyramids : " << m_pyramids.size() << std::endl;
   std::cout << "----------------------------" << std::endl;
}

void Gmesh::addNode(const Node& node)
{
   m_nodes.push_back(node);
}

void Gmesh::addVertice(const Vertice& vertice)
{
   m_vertices.push_back(vertice);
   m_number_of_elements += 1;
}

void Gmesh::addEdge(const Edge& edge)
{
   m_edges.push_back(edge);
   m_number_of_elements += 1;
}

void Gmesh::addTriangle(const Triangle& triangle)
{
   m_triangles.push_back(triangle);
   m_number_of_elements += 1;
}

void Gmesh::addQuadrangle(const Quadrangle& quad)
{
   m_quadrangles.push_back(quad);
   m_number_of_elements += 1;
}

void Gmesh::addHexahedron(const Hexahedron& hexa)
{
   m_hexahedra.push_back(hexa);
   m_number_of_elements += 1;
}

void Gmesh::addTetrahedron(const Tetrahedron& tetra)
{
   m_tetrahedra.push_back(tetra);
   m_number_of_elements += 1;
}

void Gmesh::addPrism(const Prism& prism)
{
   m_prisms.push_back(prism);
   m_number_of_elements += 1;
}

void Gmesh::addPyramid(const Pyramid& pyramid)
{
   m_pyramids.push_back(pyramid);
   m_number_of_elements += 1;
}

void Gmesh::convertInDiscontinuousMesh1()
{
   const size_t numberofNewNodes = 2 * m_edges.size();
   std::vector<Node> oldNodes(m_nodes);

   m_nodes.clear();
   m_nodes.reserve(numberofNewNodes);
   size_t indexCurrentNode = 1;

   for (size_t i = 0; i < m_edges.size(); i++) {
      std::vector<size_t> indexnodes = (m_edges.at(i)).getNodes();

      for (size_t j = 0; j < 2; j++) {
         Node tmpNode((oldNodes.at(indexnodes.at(j)-1)));
         tmpNode.changeIndex(indexCurrentNode);
         (m_edges.at(i)).changeNode(j, indexCurrentNode);
         m_nodes.push_back(tmpNode);
         indexCurrentNode += 1;
      }
   }

   assert(indexCurrentNode -1 == numberofNewNodes);
}

void Gmesh::convertInDiscontinuousMesh2()
{
   const size_t numberofNewNodes = 3 * m_triangles.size()+ 4 * m_quadrangles.size();
   std::vector<Node> oldNodes(m_nodes);

   m_nodes.clear();
   m_nodes.reserve(numberofNewNodes);
   size_t indexCurrentNode = 1;

   if(m_edges.size() > 0){
      m_number_of_elements -= m_edges.size();
      m_edges.clear();
   }

   for (size_t i = 0; i < m_triangles.size(); i++) {
      std::vector<size_t> indexnodes = m_triangles[i].getNodes();
      for (size_t j = 0; j < 3; j++) {
         Node tmpNode(oldNodes[indexnodes[j]-1]);
         tmpNode.changeIndex(indexCurrentNode);
         m_triangles[i].changeNode(j, indexCurrentNode);
         m_nodes.push_back(tmpNode);
         indexCurrentNode += 1;
      }
   }

   for (size_t i = 0; i < m_quadrangles.size(); i++) {
      std::vector<size_t> indexnodes = m_quadrangles[i].getNodes();
      for (size_t j = 0; j < 4; j++) {
         Node tmpNode(oldNodes[indexnodes[j]-1]);
         tmpNode.changeIndex(indexCurrentNode);
         m_quadrangles[i].changeNode(j, indexCurrentNode);
         m_nodes.push_back(tmpNode);
         indexCurrentNode += 1;
      }
   }

   assert(indexCurrentNode -1 == numberofNewNodes);
}

void Gmesh::convertInDiscontinuousMesh3()
{
   const size_t numberofNewNodes = 4 * m_tetrahedra.size() + 5 * m_pyramids.size()
   + 6 * m_prisms.size() + 8 * m_hexahedra.size();
   std::vector<Node> oldNodes(m_nodes);

   m_nodes.clear();
   m_nodes.reserve(numberofNewNodes);
   size_t indexCurrentNode = 1;

   if(m_edges.size() > 0){
      m_number_of_elements -= m_edges.size();
      m_edges.clear();
   }

   if(m_triangles.size() > 0){
      m_number_of_elements -= m_triangles.size();
      m_triangles.clear();
   }

   if(m_quadrangles.size() > 0){
      m_number_of_elements -= m_quadrangles.size();
      m_quadrangles.clear();
   }


   for (size_t i = 0; i < m_hexahedra.size(); i++) {
      std::vector<size_t> indexnodes = m_hexahedra[i].getNodes();
      for (size_t j = 0; j < 8; j++) {
         Node tmpNode(oldNodes[indexnodes[j]-1]);
         tmpNode.changeIndex(indexCurrentNode);
         m_hexahedra[i].changeNode(j, indexCurrentNode);
         m_nodes.push_back(tmpNode);
         indexCurrentNode += 1;
      }
   }

   for (size_t i = 0; i < m_tetrahedra.size(); i++) {
      std::vector<size_t> indexnodes = m_tetrahedra[i].getNodes();
      for (size_t j = 0; j < 4; j++) {
         Node tmpNode(oldNodes[indexnodes[j]-1]);
         tmpNode.changeIndex(indexCurrentNode);
         m_tetrahedra[i].changeNode(j, indexCurrentNode);
         m_nodes.push_back(tmpNode);
         indexCurrentNode += 1;
      }
   }

   for (size_t i = 0; i < m_prisms.size(); i++) {
      std::vector<size_t> indexnodes = m_prisms[i].getNodes();
      for (size_t j = 0; j < 6; j++) {
         Node tmpNode(oldNodes[indexnodes[j]-1]);
         tmpNode.changeIndex(indexCurrentNode);
         m_prisms[i].changeNode(j, indexCurrentNode);
         m_nodes.push_back(tmpNode);
         indexCurrentNode += 1;
      }
   }

   for (size_t i = 0; i < m_pyramids.size(); i++) {
      std::vector<size_t> indexnodes = m_pyramids[i].getNodes();
      for (size_t j = 0; j < 5; j++) {
         Node tmpNode(oldNodes[indexnodes[j]-1]);
         tmpNode.changeIndex(indexCurrentNode);
         m_pyramids[i].changeNode(j, indexCurrentNode);
         m_nodes.push_back(tmpNode);
         indexCurrentNode += 1;
      }
   }
   assert(indexCurrentNode -1 == numberofNewNodes);
}

void Gmesh::convertInDiscontinuousMesh()
{
   if(m_dim_topology ==1){
      convertInDiscontinuousMesh1();
   }
   else if(m_dim_topology == 2){
      convertInDiscontinuousMesh2();
   }
   else if(m_dim_topology == 3){
      convertInDiscontinuousMesh3();
   }
}

} //visu
