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

#ifndef GMSHELEMENT_H
#define GMSHELEMENT_H

#include<vector>
#include<array>
#include<string>

namespace visu{

class Node
{
protected:
   std::array<double, 3> m_coordinate{ {double{0.0}, double{0.0}, double{0.0}} };
   size_t m_index;
   size_t m_ref;
public:
   Node();
   Node(const std::array<double, 3>& coor, const size_t index, const size_t ref);
   void print() const;
   std::array<double, 3> getCoordinate() const;
   size_t getRef() const;
   size_t getIndex() const;
   void changeIndex(const size_t index);
   void changeRef(const size_t ref);
   void changeCoordinates(const std::array<double, 3>& coor);
};


class Element
{
protected:
   std::vector<size_t> m_nodes;
   size_t m_index;
   size_t m_type_elem;
   size_t m_physical_entities;
   size_t m_elem_tag;

public:
   Element();
   Element(const size_t type_elem);
   Element(const std::vector<size_t>& nodes, const size_t index, const size_t type_elem, const size_t physical_entities, const  size_t elem_tag);
   size_t getIndex() const;
   size_t getTypeElem() const;
   size_t getPhysicalEntities() const;
   size_t getElemTag() const;
   std::vector<size_t> getNodes() const;
   void print() const;

   void changeNodes(const std::vector<size_t>& nodes);
   void changeNode(const size_t position, const size_t node);
   void changeIndex(const size_t index);
   void changeTypeElem(const size_t type_elem);
   void changePhysicalEntities(const size_t physical_entities);
   void changeElemTag(const size_t elem_tag);
};

class Vertice : public Element
{
private:

public:
   Vertice();
   Vertice(const size_t node, const size_t index, const size_t physical_entities, const  size_t elem_tag);
};

class Edge : public Element
{
private:

public:
   Edge();
   Edge(const std::vector<size_t>& nodes, const size_t index, const size_t physical_entities, const  size_t elem_tag);
};

class Triangle : public Element
{
private:

public:
   Triangle();
   Triangle(const std::vector<size_t>& nodes, const size_t index, const size_t physical_entities, const  size_t elem_tag);
};

class Quadrangle : public Element
{
private:

public:
   Quadrangle();
   Quadrangle(const std::vector<size_t>& nodes, const size_t index, const size_t physical_entities, const  size_t elem_tag);
   size_t getNum();
};

class Tetrahedron : public Element
{
private:

public:
   Tetrahedron();
   Tetrahedron(const std::vector<size_t>& nodes, const size_t index, const size_t physical_entities, const  size_t elem_tag);
};

class Hexahedron : public Element
{
private:

public:
   Hexahedron();
   Hexahedron(const std::vector<size_t>& nodes, const size_t index, const size_t physical_entities, const  size_t elem_tag);
};

class Prism : public Element
{
private:

public:
   Prism();
   Prism(const std::vector<size_t>& nodes, const size_t index, const size_t physical_entities, const  size_t elem_tag);
};

class Pyramid : public Element
{
private:

public:
   Pyramid();
   Pyramid(const std::vector<size_t>& nodes, const size_t index, const size_t physical_entities, const  size_t elem_tag);
};

}

#endif
