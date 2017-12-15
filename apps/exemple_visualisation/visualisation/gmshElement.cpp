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


#include "gmshMesh.hpp"
#include "gmshElement.hpp"
#include<iostream>
#include<fstream>
#include <assert.h>

namespace visu{

// Class Node

Node::Node() : m_coordinate(), m_index(0), m_ref(0) {}


Node::Node(const std::array<double, 3>& coor, const size_t index, const size_t ref) :
m_coordinate(coor), m_index(index), m_ref(ref) {}

void Node::print() const
{
   std::cout << "--------------- Node --------------"  << std::endl;
   std::cout << "Index: " << m_index << std::endl;
   std::cout << "Reference: " << m_ref << std::endl;
   std::cout << "Coordinates: "  << std::endl;
   std::cout << "[x, y, z] = [ " << m_coordinate[0] << " "
   << m_coordinate[1] << " "
   << m_coordinate[2] << " ];" << std::endl;
   std::cout << "--------------------------------------" <<  std::endl;
}

std::array<double, 3> Node::getCoordinate() const
{
   return m_coordinate;
}

size_t Node::getRef() const
{
   return m_ref;
}

size_t Node::getIndex() const
{
   return m_index;
}

void Node::changeIndex(const size_t index) {m_index = index;}
void Node::changeRef(const size_t ref) {m_ref = ref;}
void Node::changeCoordinates(const std::array<double, 3>& coor) { m_coordinate = coor;}

// Class Element
Element::Element() : m_nodes(), m_index(0), m_type_elem(0),
m_physical_entities(0), m_elem_tag(0) {}

Element::Element(const size_t type_elem) : m_nodes(), m_index(0), m_type_elem(type_elem),
m_physical_entities(0), m_elem_tag(0) {}

Element::Element(const std::vector<size_t>& nodes, const size_t index, const size_t type_elem, const size_t physical_entities, const size_t elem_tag)
: m_nodes(nodes), m_index(index), m_type_elem(type_elem),
m_physical_entities(physical_entities), m_elem_tag(elem_tag) {}


size_t Element::getIndex() const
{
   return m_index;
}

size_t Element::getTypeElem() const
{
   return m_type_elem;
}

size_t Element::getPhysicalEntities() const
{
   return m_physical_entities;
}

size_t Element::getElemTag() const
{
   return m_elem_tag;
}

std::vector<size_t> Element::getNodes() const
{
   return m_nodes;
}

void Element::changeNodes(const std::vector<size_t>& nodes) {m_nodes = nodes;}
void Element::changeNode(const size_t position, const size_t node) {m_nodes[position] = node;}
void Element::changeIndex(const size_t index) {m_index = index;}
void Element::changeTypeElem(const size_t type_elem) {m_type_elem = type_elem;}
void Element::changePhysicalEntities(const size_t physical_entities) {m_physical_entities = physical_entities;}
void Element::changeElemTag(const size_t elem_tag) {m_elem_tag = elem_tag;}

void Element::print() const
{
   std::cout << "--------------- Element --------------"  << std::endl;
   std::cout << "Index: " << m_index << std::endl;
   std::cout << "Type_elem: " << m_type_elem << std::endl;
   std::cout << "Physical_entities: " << m_physical_entities << std::endl;
   std::cout << "ElemTag: " << m_elem_tag << std::endl;
   std::cout << "Number of Node: " << m_nodes.size() << std::endl;

   for (size_t i = 0; i < m_nodes.size(); i++) {
      std::cout << "Node " << i << " " << m_nodes[i]  << std::endl;
   }

   std::cout << "--------------------------------------" <<  std::endl;
}

// Class Vertice

Vertice::Vertice() : Element(15) {}


Vertice::Vertice(const size_t node, const size_t index, const size_t physical_entities, const  size_t elem_tag)
: Element((std::vector<size_t> {node}), index, 15, physical_entities, elem_tag)
{
   assert((m_nodes.size() == 1) && "The size of the vertice is not equal to 1 !");
}

// Class Edge

Edge::Edge() : Element(1) {}


Edge::Edge(const std::vector<size_t>& nodes, const size_t index, const size_t physical_entities, const  size_t elem_tag)
: Element(nodes, index, 1, physical_entities, elem_tag)
{
   assert((m_nodes.size() == 2) && "The size of the edge is not equal to 2 !");
}


// Class Triangle

Triangle::Triangle() : Element(2) {}

Triangle::Triangle(const std::vector<size_t>& nodes, const size_t index, const size_t physical_entities, const  size_t elem_tag)
: Element(nodes, index, 2, physical_entities, elem_tag)
{
   assert((m_nodes.size() == 3) && "The size of the vectices (triangle) is not equal to 3 !");
}


// Class Quadrangle

Quadrangle::Quadrangle() : Element(3) {}

Quadrangle::Quadrangle(const std::vector<size_t>& nodes, const size_t index, const size_t physical_entities, const  size_t elem_tag)
: Element(nodes, index, 3, physical_entities, elem_tag)
{
   assert((m_nodes.size() == 4) && "The size of the vectices (quadrangle) is not equal to 3 !");
}


// Class Tetrahedron

Tetrahedron::Tetrahedron() : Element(4) {}


Tetrahedron::Tetrahedron(const std::vector<size_t>& nodes, const size_t index, const size_t physical_entities, const  size_t elem_tag)
: Element(nodes, index, 4, physical_entities, elem_tag)
{
   assert((m_nodes.size() == 4) && "The size of the tetrahedron is not equal to 4 !");
}




// Class Hexahedron

Hexahedron::Hexahedron() : Element(5) {}


Hexahedron::Hexahedron(const std::vector<size_t>& nodes, const size_t index, const size_t physical_entities, const  size_t elem_tag)
: Element(nodes, index, 5, physical_entities, elem_tag)
{
   assert((m_nodes.size() == 8) && "The size of the hexahedron is not equal to 8 !");
}


// Class Prism

Prism::Prism() : Element(6) {}


Prism::Prism(const std::vector<size_t>& nodes, const size_t index, const size_t physical_entities, const  size_t elem_tag)
: Element(nodes, index, 6, physical_entities, elem_tag)
{
   assert((m_nodes.size() == 6) && "The size of the prism is not equal to 6 !");
}



// Class Pyramid

Pyramid::Pyramid() : Element(7) {}


Pyramid::Pyramid(const std::vector<size_t>& nodes, const size_t index, const size_t physical_entities, const  size_t elem_tag)
: Element(nodes, index, 7, physical_entities, elem_tag)
{
   assert((m_nodes.size() == 5) && "The size of the pyramid is not equal to 5 !");
}

} //visu
