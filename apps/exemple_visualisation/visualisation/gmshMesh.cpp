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


#include "gmshMesh.h"
#include "gmshElement.h"
#include<iostream>
#include<fstream>
#include <assert.h>

namespace visu{

//  Class Gmesh

Gmesh::Gmesh() : _dim_topology(0), _nodes(), _vertices(), _edges(), _triangles(), _quadrangles(),
                 _tetrahedra() , _hexahedra(), _prisms(), _pyramids(), _number_of_elements(0) {}

Gmesh::Gmesh(const size_t dim, const std::vector<Node>& nodes, const std::vector<Vertice>& vertices,
      const std::vector<Edge>& edges) :
      _dim_topology(dim), _nodes(nodes), _vertices(vertices), _edges(edges), _triangles(), _quadrangles(),
      _tetrahedra() , _hexahedra(), _prisms(), _pyramids(),
      _number_of_elements(vertices.size() + edges.size()) {}

Gmesh::Gmesh(const size_t dim, const std::vector<Node>& nodes, const std::vector<Vertice>& vertices,
      const std::vector<Edge>& edges, const std::vector<Triangle>& triangles,
      const std::vector<Quadrangle>& quadrangles) :
      _dim_topology(dim), _nodes(nodes), _vertices(vertices), _edges(edges),
      _triangles(triangles), _quadrangles(quadrangles),
      _tetrahedra() , _hexahedra(), _prisms(), _pyramids(),
      _number_of_elements(vertices.size() + edges.size() + triangles.size() + quadrangles.size()) {}

Gmesh::Gmesh(const size_t dim, const std::vector<Node>& nodes, const std::vector<Vertice>& vertices,
      const std::vector<Edge>& edges, const std::vector<Triangle>& triangles,
      const std::vector<Quadrangle>& quadrangles, const std::vector<Tetrahedron>& tetrahedra,
      const std::vector<Hexahedron> hexahedra) :
            _dim_topology(dim), _nodes(nodes), _vertices(vertices), _edges(edges),
            _triangles(triangles), _quadrangles(quadrangles),
            _tetrahedra(tetrahedra) , _hexahedra(hexahedra), _prisms(), _pyramids(),
            _number_of_elements(vertices.size() + edges.size() + triangles.size() + quadrangles.size()
            + tetrahedra.size() + hexahedra.size()) {}

Gmesh::Gmesh(const size_t dim, const std::vector<Node>& nodes, const std::vector<Vertice>& vertices,
      const std::vector<Edge>& edges, const std::vector<Triangle>& triangles,
      const std::vector<Quadrangle>& quadrangles, const std::vector<Tetrahedron>& tetrahedra,
      const std::vector<Hexahedron> hexahedra, const std::vector<Prism>& prisms,
      const std::vector<Pyramid>& pyramids) :
      _dim_topology(dim), _nodes(nodes), _vertices(vertices), _edges(edges),
      _triangles(triangles), _quadrangles(quadrangles),
      _tetrahedra(tetrahedra) , _hexahedra(hexahedra), _prisms(prisms), _pyramids(pyramids),
      _number_of_elements(vertices.size() + edges.size() + triangles.size() + quadrangles.size()
      + tetrahedra.size() + hexahedra.size() + prisms.size() + pyramids.size()) {}

Gmesh::Gmesh(const size_t dim, const std::vector<Node>& nodes, const std::vector<Vertice>& vertices,
      const std::vector<Edge>& edges, const std::vector<Triangle>& triangles,
      const std::vector<Tetrahedron>& tetrahedra) :
            _dim_topology(dim), _nodes(nodes), _vertices(vertices), _edges(edges),
            _triangles(triangles), _quadrangles(),
            _tetrahedra(tetrahedra) , _hexahedra(), _prisms(), _pyramids(),
            _number_of_elements(vertices.size() + edges.size() + triangles.size()
            + tetrahedra.size()) {}

Gmesh::Gmesh(const size_t dim, const std::vector<Node>& nodes, const std::vector<Vertice>& vertices,
            const std::vector<Edge>& edges, const std::vector<Quadrangle>& quadrangles,
            const std::vector<Hexahedron>& hexahedra) :
               _dim_topology(dim), _nodes(nodes), _vertices(vertices), _edges(edges),
               _triangles(), _quadrangles(quadrangles),
               _tetrahedra() , _hexahedra(hexahedra), _prisms(), _pyramids(),
               _number_of_elements(vertices.size() + edges.size() + quadrangles.size()
               + hexahedra.size()) {}

size_t Gmesh::getNumberofNodes() const
{
  return _nodes.size();
}

size_t Gmesh::getNumberofElements() const
{
  return _number_of_elements;
}

size_t Gmesh::getDim() const
{
   return _dim_topology;
}

std::vector<Node> Gmesh::getNodes() const
{
  return _nodes;
}

std::vector<Vertice> Gmesh::getVertices() const
{
  return _vertices;
}

std::vector<Edge> Gmesh::getEdges() const
{
  return _edges;
}

std::vector<Triangle> Gmesh::getTriangles() const
{
  return _triangles;
}

std::vector<Quadrangle> Gmesh::getQuadrangles() const
{
  return _quadrangles;
}

std::vector<Tetrahedron> Gmesh::getTetrahedra() const
{
  return _tetrahedra;
}

std::vector<Hexahedron> Gmesh::getHexahedra() const
{
  return _hexahedra;
}

std::vector<Prism> Gmesh::getPrisms() const
{
  return _prisms;
}

std::vector<Pyramid> Gmesh::getPyramids() const
{
  return _pyramids;
}

Node Gmesh::getNode(const size_t index) const
{
  return _nodes.at(index);
}

Vertice Gmesh::getVertice(const size_t index) const
{
  return _vertices.at(index);
}

Edge Gmesh::getEdge(const size_t index) const
{
  return _edges.at(index);
}

Triangle Gmesh::getTriangle(const size_t index) const
{
  return _triangles.at(index);
}

Quadrangle Gmesh::getQuadrangle(const size_t index) const
{
  return _quadrangles.at(index);
}

Hexahedron Gmesh::getHexahedron(const size_t index) const
{
  return _hexahedra.at(index);
}

Tetrahedron Gmesh::getTetrahedron(const size_t index) const
{
  return _tetrahedra.at(index);
}

Prism Gmesh::getPrism(const size_t index) const
{
  return _prisms.at(index);
}

Pyramid Gmesh::getPyramid(const size_t index) const
{
  return _pyramids.at(index);
}

void Gmesh::readGmesh_MEDITformat(const std::string name_mesh)
{
  std::ifstream mesh_file(name_mesh.data());
  if (!mesh_file.is_open())
  {
    std::cout << "Unable to open file " << name_mesh << std::endl;
    abort();
  }

  std::cout << "------------------------" << std::endl;
  std::cout << "Reading mesh " << name_mesh << std::endl;
  std::cout << "MEDIT Format" << std::endl;
  std::cout << "------------------------" << std::endl;


  std::string file_line;
  size_t nb_nodes(0), nb_edges(0), nb_triangles(0), nb_tetras(0), nb_quad(0), nb_hexa(0);
  std::vector<double> vert(3);
  std::vector<size_t> edg(2), tri(3), tet(4), quad(4), hexa(8);
  size_t ref(0), loop_nodes(100), offset_elements(0);

  while (!mesh_file.eof())
  {
    getline(mesh_file, file_line);
    if ((file_line.find("Vertices") != std::string::npos)&&(loop_nodes))
    {
        mesh_file >> nb_nodes;
        std::cout << "Reading vertices (" << nb_nodes << ")" << std::endl;
        for (size_t i = 1 ; i <= nb_nodes ; i++)
        {
          mesh_file >> vert[0] >> vert[1] >> vert[2] >> ref;
          Node node(vert, i, ref);
          _nodes.push_back(node);
        }
        assert((_nodes.size() == nb_nodes) && "Problem 1 in mesh reading.");
        loop_nodes = 0;
    }
    if (file_line.find("Edges") != std::string::npos)
    {
        mesh_file >> nb_edges;
        std::cout << "Reading edges (" << nb_edges << ")" << std::endl;
        for (size_t i = 1 ; i <= nb_edges ; i++)
        {
          mesh_file >> edg[0] >> edg[1] >> ref;
          Edge edge(edg, offset_elements + i, ref, 0);
          _edges.push_back(edge);
        }

        offset_elements += nb_edges;
        assert((_edges.size() == nb_edges) && "Problem 2 in mesh reading.");
    }
    if (file_line.find("Triangles") != std::string::npos)
    {
        mesh_file >> nb_triangles;
        std::cout << "Reading triangles (" << nb_triangles << ")" << std::endl;
        for (size_t i = 1 ; i <= nb_triangles ; i++)
        {
          mesh_file >> tri[0] >> tri[1] >> tri[2] >> ref;
          Triangle triangle(tri, offset_elements + i, ref, 0);
          _triangles.push_back(triangle);
        }

        offset_elements += nb_triangles;
        assert((_triangles.size() == nb_triangles) && "Problem 3 in mesh reading.");
    }
    if (file_line.find("Tetrahedra") != std::string::npos)
    {
        mesh_file >> nb_tetras;
        std::cout << "Reading tetrahedra (" << nb_tetras << ")" << std::endl;
        for (size_t i = 0 ; i < nb_tetras ; i++)
        {
          mesh_file >> tet[0] >> tet[1] >> tet[2] >> tet[3] >> ref;
          Tetrahedron tetrahedron(tet, offset_elements + i, ref, 0);
          _tetrahedra.push_back(tetrahedron);
        }

        offset_elements += nb_tetras;
        assert((_tetrahedra.size() == nb_tetras) && "Problem 4 in mesh reading.");
    }
    if (file_line.find("Quadrilaterals") != std::string::npos)
    {
        mesh_file >> nb_quad;
        std::cout << "Reading quadrilaterals (" << nb_quad << ")" << std::endl;
        for (size_t i = 1 ; i <= nb_quad ; i++)
        {
          mesh_file >> quad[0] >> quad[1] >> quad[2] >> quad[3]  >> ref;
          Quadrangle quadrangle(quad, offset_elements + i, ref, 0);
          _quadrangles.push_back(quadrangle);
        }

        offset_elements += nb_quad;
        assert((_quadrangles.size() == nb_quad) && "Problem 3 in mesh reading.");
    }
    if (file_line.find("Hexahedra") != std::string::npos)
    {
        mesh_file >> nb_hexa;
        std::cout << "Reading hexahedra (" << nb_hexa << ")" << std::endl;
        for (size_t i = 0 ; i < nb_hexa ; i++)
        {
          mesh_file >> hexa[0] >> hexa[1] >> hexa[2] >> hexa[3] >> hexa[4]
                    >> hexa[5] >> hexa[6] >> hexa[7] >> hexa[8] >> ref;
          Hexahedron hexahedron(hexa, offset_elements + i, ref, 0);
          _hexahedra.push_back(hexahedron);
        }

        offset_elements += nb_hexa;
        assert((_hexahedra.size() == nb_hexa) && "Problem 4 in mesh reading.");
    }
  }

   assert(offset_elements == (nb_edges + nb_triangles + nb_quad + nb_hexa + nb_tetras));
  _number_of_elements = offset_elements;

  mesh_file.close();

  if((_tetrahedra.size() + _hexahedra.size() + _prisms.size() + _pyramids.size()) > 0){
     _dim_topology = 3;
  }
  else if((_triangles.size() + _quadrangles.size()) > 0){
    _dim_topology = 2;
  }
  else if(_edges.size() > 0){
   _dim_topology = 1;
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

   std::cout << "------------------------" << std::endl;
   std::cout << "Reading mesh " << name_mesh << std::endl;
   std::cout << "MSH Format"  << std::endl;
   std::cout << "------------------------" << std::endl;


   std::string file_line;
   size_t nb_nodes(0), nb_elements(0);
   std::vector<double> vert(3);
   size_t index(0), type_elem(0), nbinteger_tag(0), loop_nodes(100), loop_elem(100);

   while (!mesh_file.eof() && loop_elem)
   {
      getline(mesh_file, file_line);
      if ((file_line.find("$Nodes") != std::string::npos)&&(loop_nodes))
      {
         mesh_file >> nb_nodes;
         std::cout << "Reading nodes (" << nb_nodes << ")" << std::endl;
         for (size_t i = 0 ; i < nb_nodes ; ++i)
         {
            mesh_file >> index >> vert[0] >> vert[1] >> vert[2];
            Node node(vert, index, 1);
            _nodes.push_back(node);
         }
         assert((_nodes.size() == nb_nodes) && "Problem 1 in mesh reading.");
         loop_nodes = 0;
      }
      if ((file_line.find("$Elements") != std::string::npos) && (loop_elem))
      {
         mesh_file >> nb_elements;
         _number_of_elements = nb_elements;
         for (size_t i = 0 ; i < nb_elements ; ++i)
         {
            std::cout << i << " / " << nb_elements - 1 << std::endl;
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
               _edges.push_back(edge);
            }
            else if(type_elem == 2){
               for(size_t inode = 0; inode < 3; inode++){
                  mesh_file >> ind_node;
                  elem_nodes.push_back(ind_node);
               }
               Triangle tri(elem_nodes, i+1, integer_tag[0], integer_tag[1]);
               _triangles.push_back(tri);
            }
            else if(type_elem == 3){
               for(size_t inode = 0; inode < 4; inode++){
                  mesh_file >> ind_node;
                  elem_nodes.push_back(ind_node);
               }
               Quadrangle quad(elem_nodes, i+1, integer_tag[0], integer_tag[1]);
               _quadrangles.push_back(quad);
            }
            else if(type_elem == 4){
               for(size_t inode = 0; inode < 4; inode++){
                  mesh_file >> ind_node;
                  elem_nodes.push_back(ind_node);
               }
               Tetrahedron tetra(elem_nodes, i+1, integer_tag[0], integer_tag[1]);
               _tetrahedra.push_back(tetra);
            }
            else if(type_elem == 5){
               for(size_t inode = 0; inode < 8; inode++){
                  mesh_file >> ind_node;
                  elem_nodes.push_back(ind_node);
               }
               Hexahedron hexa(elem_nodes, i+1, integer_tag[0], integer_tag[1]);
               _hexahedra.push_back(hexa);
            }
            else if(type_elem == 6){
               for(size_t inode = 0; inode < 6; inode++){
                  mesh_file >> ind_node;
                  elem_nodes.push_back(ind_node);
               }
               Prism prism(elem_nodes, i+1, integer_tag[0], integer_tag[1]);
               _prisms.push_back(prism);
            }
            else if(type_elem == 7){
               for(size_t inode = 0; inode < 5; inode++){
                  mesh_file >> ind_node;
                  elem_nodes.push_back(ind_node);
               }
               Pyramid pyra(elem_nodes, i+1, integer_tag[0], integer_tag[1]);
               _pyramids.push_back(pyra);
            }
            else {
               std::cout << "Unsupported element" << std::endl;
            }
         }
         loop_elem = 0;
      }
   }

   mesh_file.close();

   if((_tetrahedra.size() + _hexahedra.size() + _prisms.size() + _pyramids.size()) > 0){
     _dim_topology = 3;
   }
   else if((_triangles.size() + _quadrangles.size()) > 0){
     _dim_topology = 2;
   }
   else if(_edges.size() > 0){
    _dim_topology = 1;
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
       std::cout << "------------------------" << std::endl;
       std::cout << "Writing mesh " << name_mesh << std::endl;
       std::cout << "MESH file format"  << std::endl;
       std::cout << "------------------------" << std::endl;
    }

    mesh_file << "MeshVersionFormatted 2" << std::endl;
    mesh_file << "Dimension 3"  << std::endl;
    mesh_file << " " << std::endl;

    mesh_file << "Vertices" << std::endl;
    mesh_file << _nodes.size() << std::endl;

    for (size_t i = 0; i < _nodes.size(); i++){
      Node node(_nodes[i]);
      mesh_file         << (node.getCoordinate())[0]
                << " "  << (node.getCoordinate())[1]
                << " "  << (node.getCoordinate())[2]
                << " "  << node.getRef()
                << std::endl;
    }

    mesh_file << " " << std::endl;
    mesh_file << "Edges" << std::endl;
    mesh_file << _edges.size() << std::endl;

    for (size_t i = 0; i < _edges.size(); i++){
      Edge edge(_edges[i]);
      mesh_file           << (edge.getNodes())[0]
                  << " "  << (edge.getNodes())[1]
                  << " "  << edge.getPhysicalEntities()
                  << std::endl;
    }

    mesh_file << " " << std::endl;
    mesh_file << "Triangles" << std::endl;
    mesh_file << _triangles.size() << std::endl;

    for (size_t i = 0; i < _triangles.size(); i++){
       Triangle tri(_triangles[i]);
       mesh_file          << (tri.getNodes())[0]
                  << " "  << (tri.getNodes())[1]
                  << " "  << (tri.getNodes())[2]
                  << " "  << tri.getPhysicalEntities()
                  << std::endl;
    }

    mesh_file << " " << std::endl;
    mesh_file << "Quadrilaterals" << std::endl;
    mesh_file << _quadrangles.size() << std::endl;

    for (size_t i = 0; i < _quadrangles.size(); i++){
       Quadrangle quad(_quadrangles[i]);
       mesh_file          << (quad.getNodes())[0]
                  << " "  << (quad.getNodes())[1]
                  << " "  << (quad.getNodes())[2]
                  << " "  << (quad.getNodes())[3]
                  << " "  << quad.getPhysicalEntities()
                  << std::endl;
    }

    mesh_file << " " << std::endl;
    mesh_file << "Tetrahedra" << std::endl;
    mesh_file << _tetrahedra.size() << std::endl;

    for (size_t i = 0; i < _tetrahedra.size(); i++){
       Tetrahedron tetra(_tetrahedra[i]);
       mesh_file          << (tetra.getNodes())[0]
                  << " "  << (tetra.getNodes())[1]
                  << " "  << (tetra.getNodes())[2]
                  << " "  << (tetra.getNodes())[3]
                  << " "  << tetra.getPhysicalEntities()
                  << std::endl;
    }

    mesh_file << " " << std::endl;
    mesh_file << "Hexahedra" << std::endl;
    mesh_file << _hexahedra.size() << std::endl;

    for (size_t i = 0; i < _quadrangles.size(); i++){
       Hexahedron hexa(_hexahedra[i]);
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
      std::cout << "------------------------" << std::endl;
      std::cout << "Writing mesh " << name_mesh << std::endl;
      std::cout << "MSH ASCII file format"  << std::endl;
      std::cout << "------------------------" << std::endl;
   }

   size_t offset_elements(1);

   mesh_file << "$MeshFormat" << std::endl;
   mesh_file << "2.2 0 8" << std::endl;
   mesh_file << "$EndMeshFormat" << std::endl;

   mesh_file << "$Nodes" << std::endl;
   mesh_file << _nodes.size() << std::endl;

   for (size_t i = 0; i < _nodes.size(); i++){
     Node node(_nodes[i]);
     mesh_file << node.getIndex()     << " "  << (node.getCoordinate())[0]
                                    << " "  << (node.getCoordinate())[1]
                                    << " "  << (node.getCoordinate())[2]
                                    << std::endl;
   }

   mesh_file << "$EndNodes" << std::endl;
   mesh_file << "$Elements" << std::endl;
   mesh_file << _number_of_elements << std::endl;

   for (size_t i = 0; i < _vertices.size(); i++){
     Vertice vert(_vertices[i]);
     mesh_file << offset_elements << " "  << 15
                                  << " "  << 2 //number of tag
                                  << " "  << vert.getPhysicalEntities()
                                  << " "  << vert.getElemTag()
                                  << " "  << (vert.getNodes())[0]
                                  << std::endl;
     offset_elements += 1;
   }

   for (size_t i = 0; i < _edges.size(); i++){
     Edge edge(_edges[i]);
     mesh_file << offset_elements << " "  << 1
                                  << " "  << 2 //number of tag
                                  << " "  << edge.getPhysicalEntities()
                                  << " "  << edge.getElemTag()
                                  << " "  << (edge.getNodes())[0]
                                  << " "  << (edge.getNodes())[1]
                                  << std::endl;
     offset_elements += 1;
   }

   for (size_t i = 0; i < _triangles.size(); i++){
      Triangle tri(_triangles[i]);
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

   for (size_t i = 0; i < _quadrangles.size(); i++){
      Quadrangle quad(_quadrangles[i]);
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

   for (size_t i = 0; i < _tetrahedra.size(); i++){
      Tetrahedron tet(_tetrahedra[i]);
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

   for (size_t i = 0; i < _hexahedra.size(); i++){
      Hexahedron hexa(_hexahedra[i]);
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

   for (size_t i = 0; i < _prisms.size(); i++){
      Prism pri(_prisms[i]);
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

   for (size_t i = 0; i < _pyramids.size(); i++){
      Pyramid py(_pyramids[i]);
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
   std::cout << "Dimension : " << _dim_topology << std::endl;
   std::cout << " * Numbers of vetices : " << _nodes.size() << std::endl;
   std::cout << " * Numbers of edges : " << _edges.size() << std::endl;
   std::cout << " * Numbers of triangles : " << _triangles.size() << std::endl;
   std::cout << " * Numbers of quadrangles : " << _quadrangles.size() << std::endl;
   std::cout << " * Numbers of tetrahedra : " << _tetrahedra.size() << std::endl;
   std::cout << " * Numbers of hexahedra : " << _hexahedra.size() << std::endl;
   std::cout << " * Numbers of prisms : " << _prisms.size() << std::endl;
   std::cout << " * Numbers of pyramids : " << _pyramids.size() << std::endl;
   std::cout << "----------------------------" << std::endl;
}

void Gmesh::addNode(const Node& node)
{
   _nodes.push_back(node);
}

void Gmesh::addVertice(const Vertice& vertice)
{
   _vertices.push_back(vertice);
   _number_of_elements += 1;
}

void Gmesh::addEdge(const Edge& edge)
{
   _edges.push_back(edge);
   _number_of_elements += 1;
}

void Gmesh::addTriangle(const Triangle& triangle)
{
   _triangles.push_back(triangle);
   _number_of_elements += 1;
}

void Gmesh::addQuadrangle(const Quadrangle& quad)
{
   _quadrangles.push_back(quad);
   _number_of_elements += 1;
}

void Gmesh::addHexahedron(const Hexahedron& hexa)
{
   _hexahedra.push_back(hexa);
   _number_of_elements += 1;
}

void Gmesh::addTetrahedron(const Tetrahedron& tetra)
{
   _tetrahedra.push_back(tetra);
   _number_of_elements += 1;
}

void Gmesh::addPrism(const Prism& prism)
{
   _prisms.push_back(prism);
   _number_of_elements += 1;
}

void Gmesh::addPyramid(const Pyramid& pyramid)
{
   _pyramids.push_back(pyramid);
   _number_of_elements += 1;
}

void Gmesh::convertInDiscontinuousMesh1()
{
   const size_t numberofNewNodes = 2 * _edges.size();
   std::vector<Node> oldNodes(_nodes);

   _nodes.clear();
   _nodes.reserve(numberofNewNodes);
   size_t indexCurrentNode = 1;

   for (size_t i = 0; i < _edges.size(); i++) {
      std::vector<size_t> indexnodes = (_edges.at(i)).getNodes();

      for (size_t j = 0; j < 2; j++) {
         Node tmpNode((oldNodes.at(indexnodes.at(j)-1)));
         tmpNode.changeIndex(indexCurrentNode);
         (_edges.at(i)).changeNode(j, indexCurrentNode);
         _nodes.push_back(tmpNode);
         indexCurrentNode += 1;
      }
   }

   assert(indexCurrentNode -1 == numberofNewNodes);
}

void Gmesh::convertInDiscontinuousMesh2()
{
   const size_t numberofNewNodes = 3 * _triangles.size()
                                    + 4 * _quadrangles.size();
   std::vector<Node> oldNodes(_nodes);

   _nodes.clear();
   _nodes.reserve(numberofNewNodes);
   size_t indexCurrentNode = 1;

   if(_edges.size() > 0){
      _number_of_elements -= _edges.size();
      _edges.clear();
   }

   for (size_t i = 0; i < _triangles.size(); i++) {
      std::vector<size_t> indexnodes = _triangles[i].getNodes();
      for (size_t j = 0; j < 3; j++) {
         Node tmpNode(oldNodes[indexnodes[j]-1]);
         tmpNode.changeIndex(indexCurrentNode);
         _triangles[i].changeNode(j, indexCurrentNode);
         _nodes.push_back(tmpNode);
         indexCurrentNode += 1;
      }
   }

   for (size_t i = 0; i < _quadrangles.size(); i++) {
      std::vector<size_t> indexnodes = _quadrangles[i].getNodes();
      for (size_t j = 0; j < 4; j++) {
         Node tmpNode(oldNodes[indexnodes[j]-1]);
         tmpNode.changeIndex(indexCurrentNode);
         _quadrangles[i].changeNode(j, indexCurrentNode);
         _nodes.push_back(tmpNode);
         indexCurrentNode += 1;
      }
   }

   assert(indexCurrentNode -1 == numberofNewNodes);
}

void Gmesh::convertInDiscontinuousMesh3()
{
   const size_t numberofNewNodes = 4 * _tetrahedra.size() + 5 * _pyramids.size()
                                    + 6 * _prisms.size() + 8 * _hexahedra.size();
   std::vector<Node> oldNodes(_nodes);

   _nodes.clear();
   _nodes.reserve(numberofNewNodes);
   size_t indexCurrentNode = 1;

   if(_edges.size() > 0){
      _number_of_elements -= _edges.size();
      _edges.clear();
   }

   if(_triangles.size() > 0){
      _number_of_elements -= _triangles.size();
      _triangles.clear();
   }

   if(_quadrangles.size() > 0){
      _number_of_elements -= _quadrangles.size();
      _quadrangles.clear();
   }


   for (size_t i = 0; i < _hexahedra.size(); i++) {
      std::vector<size_t> indexnodes = _hexahedra[i].getNodes();
      for (size_t j = 0; j < 8; j++) {
         Node tmpNode(oldNodes[indexnodes[j]-1]);
         tmpNode.changeIndex(indexCurrentNode);
         _hexahedra[i].changeNode(j, indexCurrentNode);
         _nodes.push_back(tmpNode);
         indexCurrentNode += 1;
      }
   }

   for (size_t i = 0; i < _tetrahedra.size(); i++) {
      std::vector<size_t> indexnodes = _tetrahedra[i].getNodes();
      for (size_t j = 0; j < 4; j++) {
         Node tmpNode(oldNodes[indexnodes[j]-1]);
         tmpNode.changeIndex(indexCurrentNode);
         _tetrahedra[i].changeNode(j, indexCurrentNode);
         _nodes.push_back(tmpNode);
         indexCurrentNode += 1;
      }
   }

   for (size_t i = 0; i < _prisms.size(); i++) {
      std::vector<size_t> indexnodes = _prisms[i].getNodes();
      for (size_t j = 0; j < 6; j++) {
         Node tmpNode(oldNodes[indexnodes[j]-1]);
         tmpNode.changeIndex(indexCurrentNode);
         _prisms[i].changeNode(j, indexCurrentNode);
         _nodes.push_back(tmpNode);
         indexCurrentNode += 1;
      }
   }

   for (size_t i = 0; i < _pyramids.size(); i++) {
      std::vector<size_t> indexnodes = _pyramids[i].getNodes();
      for (size_t j = 0; j < 5; j++) {
         Node tmpNode(oldNodes[indexnodes[j]-1]);
         tmpNode.changeIndex(indexCurrentNode);
         _pyramids[i].changeNode(j, indexCurrentNode);
         _nodes.push_back(tmpNode);
         indexCurrentNode += 1;
      }
   }
   assert(indexCurrentNode -1 == numberofNewNodes);
}

void Gmesh::convertInDiscontinuousMesh()
{
   if(_dim_topology ==1){
      convertInDiscontinuousMesh1();
   }
   else if(_dim_topology == 2){
      convertInDiscontinuousMesh2();
   }
   else if(_dim_topology == 3){
      convertInDiscontinuousMesh3();
   }
}

void Gmesh::computeDeformed(const std::vector<Node>& newNodes)
{
   assert(newNodes.size() == _nodes.size());
   for (size_t i = 0; i < newNodes.size(); i++) {
      _nodes[i].changeCoordinates(newNodes[i].getCoordinate());
   }
}

} //visu
