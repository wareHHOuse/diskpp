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

 #include<vector>
 #include<string>
 #include<utility>
 #include<fstream>
 #include<iostream>

 #include "gmshMesh.hpp"
 #include "gmshData.hpp"
 #include "gmshElement.hpp"


namespace visu{


// Class Data

Data::Data() : m_index(0), m_data() {}

Data::Data(const size_t index, const std::vector<double>& data)
            : m_index(index), m_data(data) {}

size_t Data::getNbComposante() const
{
   return m_index;
}

std::vector<double> Data::getData() const
{
   return m_data;
}

size_t Data::getIndex() const {return m_index;}
void Data::changeData(const size_t indice, const double value) {m_data.at(indice) = value;}
void Data::changeData(const std::vector<double>& data)
{
   m_data.clear();
   m_data.reserve(data.size());

   m_data = data;
}
void Data::changeIndex(const size_t index) {m_index = index;}

//Class Subdata
SubData::SubData() : Data(), m_node() {}

SubData::SubData(const std::vector<double>& data, const Node& node) :
Data(node.getIndex(), data), m_node(node) {}

SubData::SubData(const Data& data, const Node& node) :
SubData(data.getData(), node) {}

Node SubData::getNode() const {return m_node;}

void SubData::changeCoordinates(const std::array<double, 3>& coor)
{
   m_node.changeCoordinates(coor);
}

void SubData::changeIndex(const size_t index)
{
   m_index = index;
   m_node.changeIndex(index);
}

// Class GenericData
GenericData::GenericData() : m_composante(0), m_time_value(0.0), m_title("  "), m_datas() {}

GenericData::GenericData(const size_t nbcompo, const std::vector<Data>& datas) :
m_composante(nbcompo), m_time_value(0.0), m_title("\" Unknown title \" "), m_datas(datas) {}

GenericData::GenericData(const size_t nbcompo, const double time_value, const std::string title, const std::vector<Data>& datas)
: m_composante(nbcompo), m_time_value(time_value), m_title(title), m_datas(datas) {}

std::vector<Data> GenericData::getData() const {return m_datas;}

double GenericData::getTime() const {return m_time_value;}

std::string GenericData::getTitle() const {return m_title;}

size_t GenericData::getNbComposante() const {return m_composante;}

void GenericData::changeNbComposante(const size_t nbcompo) {m_composante = nbcompo;}

void GenericData::changeTime(const double time) {m_time_value = time;}

void GenericData::changeTitle(const std::string title) {m_title = title;}

void GenericData::addData(const Data& data) {m_datas.push_back(data);}

void GenericData::addData(const size_t index, const std::vector<double>& values)
{
   Data tmp(index, values);
   m_datas.push_back(tmp);
}


// Class NodeData
NodeData::NodeData() : GenericData(), m_subdatas() {}
NodeData::NodeData(const size_t nbcompo, const double time_value, const std::string title,
   const std::vector<Data>& datas, const std::vector<SubData>& subdatas) :
   GenericData(nbcompo, time_value, title, datas), m_subdatas(subdatas) {}

void NodeData::addSubData(const Data& subdata, const Node& node)
{
   SubData tmp(subdata, node);
   m_subdatas.push_back(tmp);
}

void NodeData::addSubData(const std::vector<double>& subdata, const Node& node)
{
   SubData tmp(subdata, node);
   m_subdatas.push_back(tmp);
}

void NodeData::writeNodeData(const std::string name_mesh) const
{
   std::ofstream mesh_file;
   mesh_file.open(name_mesh, std::ofstream::out | std::ofstream::app);
   if (!mesh_file.is_open())
   {
      std::cout << "Unable to open file " << name_mesh << std::endl;
      abort();
   }
   else
   {
      std::cout << "------------------------" << std::endl;
      std::cout << "Writing Results " << name_mesh << std::endl;
      std::cout << "MSH ASCII file format"  << std::endl;
      std::cout << "------------------------" << std::endl;
   }

   mesh_file << "$NodeData" << std::endl;
   mesh_file << 1 << std::endl;  // string tag

   mesh_file << m_title  << std::endl; // name of the view
   mesh_file << 1 << std::endl; // real tag
   mesh_file << m_time_value << std::endl; // time value
   mesh_file << 3 << std::endl; // integer tag
   mesh_file << 0 << std::endl;  // time Step

   mesh_file << m_composante << std::endl;
   mesh_file << m_datas.size() + m_subdatas.size() << std::endl;

   for (size_t i = 0; i < m_datas.size(); i++) {
      mesh_file   << m_datas[i].getIndex() << " ";

      std::vector<double> tmpresu(m_datas[i].getData());
      for (size_t j = 0; j < std::min(m_composante, tmpresu.size()) ; j++) {
         mesh_file   << tmpresu[j]   << " ";
      }
      for (size_t j = std::min(m_composante, tmpresu.size()); j < m_composante; j++) {
         mesh_file   << 0.0   << " ";
      }

      mesh_file   <<   std::endl;
   }

   for (size_t i = 0; i < m_subdatas.size(); i++) {
      mesh_file   << m_subdatas[i].getIndex() << " ";

      std::vector<double> tmpresu(m_subdatas[i].getData());
      for (size_t j = 0; j < std::min(m_composante, tmpresu.size()) ; j++) {
         mesh_file   << tmpresu[j]   << " ";
      }
      for (size_t j = std::min(m_composante, tmpresu.size()); j < m_composante; j++) {
         mesh_file   << 0.0   << " ";
      }

      mesh_file   <<   std::endl;
   }

   mesh_file << "$EndNodeData" << std::endl;
   mesh_file.close();

}

void NodeData::saveNodeData(const std::string name_mesh, const Gmesh& mesh)
{
   Gmesh meshtmp(mesh);

   if(m_subdatas.size() > 0){
      size_t indnode = mesh.getNumberofNodes();
      for (size_t i = 0; i < m_subdatas.size(); i++) {
         m_subdatas[i].changeIndex(indnode + i + 1);
         meshtmp.addNode(m_subdatas[i].getNode());
         Vertice vert(indnode + i + 1, indnode + i + 1, 1 , 1);
         meshtmp.addVertice(vert);
      }
   }
   meshtmp.writeGmesh(name_mesh, 2);
   writeNodeData(name_mesh);
}

// Class ElementData


void ElementData::writeElementData(const std::string name_mesh) const
{
   std::ofstream mesh_file;
   mesh_file.open(name_mesh, std::ofstream::out | std::ofstream::app);
   if (!mesh_file.is_open())
   {
      std::cout << "Unable to open file " << name_mesh << std::endl;
      abort();
   }
   else
   {
      std::cout << "------------------------" << std::endl;
      std::cout << "Writing Results " << name_mesh << std::endl;
      std::cout << "MSH ASCII file format"  << std::endl;
      std::cout << "------------------------" << std::endl;
   }

   mesh_file << "$ElementData" << std::endl;
   mesh_file << 1 << std::endl;  // string tag

   mesh_file << m_title  << std::endl; // name of the view
   mesh_file << 1 << std::endl; // real tag
   mesh_file << m_time_value << std::endl; // time value
   mesh_file << 3 << std::endl; // integer tag
   mesh_file << 0 << std::endl;  // time Step

   mesh_file << m_composante << std::endl;
   mesh_file << m_datas.size() << std::endl;

   for (size_t i = 0; i < m_datas.size(); i++) {
      mesh_file   << m_datas[i].getIndex() << " ";

      std::vector<double> tmpresu(m_datas[i].getData());
      for (size_t j = 0; j < std::min(m_composante, tmpresu.size()) ; j++) {
         mesh_file   << tmpresu[j]   << " ";
      }
      for (size_t j = std::min(m_composante, tmpresu.size()); j < m_composante; j++) {
         mesh_file   << 0.0   << " ";
      }

      mesh_file   <<   std::endl;
   }

   mesh_file << "$EndElementData" << std::endl;
   mesh_file.close();

}

void ElementData::saveElementData(const std::string name_mesh, const Gmesh& mesh) const
{
   mesh.writeGmesh(name_mesh, 2);
   writeElementData(name_mesh);
}

//Class ElementNodeData

void ElementNodeData::writeElementNodeData(const std::string name_mesh) const
{
   std::ofstream mesh_file;
   mesh_file.open(name_mesh, std::ofstream::out | std::ofstream::app);
   if (!mesh_file.is_open())
   {
      std::cout << "Unable to open file " << name_mesh << std::endl;
      abort();
   }
   else
   {
      std::cout << "------------------------" << std::endl;
      std::cout << "Writing Results " << name_mesh << std::endl;
      std::cout << "MSH ASCII file format"  << std::endl;
      std::cout << "------------------------" << std::endl;
   }

   mesh_file << "$ElementNodeData" << std::endl;
   mesh_file << 1 << std::endl;  // string tag

   mesh_file << m_title  << std::endl; // name of the view
   mesh_file << 1 << std::endl; // real tag
   mesh_file << m_time_value << std::endl; // time value
   mesh_file << 3 << std::endl; // integer tag
   mesh_file << 0 << std::endl;  // time Step

   mesh_file << m_composante << std::endl;
   mesh_file << m_datas.size() << std::endl;

   for (size_t i = 0; i < m_datas.size(); i++) {
      mesh_file   << m_datas[i].getIndex() << " ";

      std::vector<double> tmpresu(m_datas[i].getData());
      for (size_t j = 0; j < std::min(m_composante, tmpresu.size()) ; j++) {
         mesh_file   << tmpresu[j]   << " ";
      }
      for (size_t j = std::min(m_composante, tmpresu.size()); j < m_composante; j++) {
         mesh_file   << 0.0   << " ";
      }

      mesh_file   <<   std::endl;
   }

   mesh_file << "$EndElementNodeData" << std::endl;
   mesh_file.close();

}

void ElementNodeData::saveElementNodeData(const std::string name_mesh, const Gmesh& mesh) const
{
   mesh.writeGmesh(name_mesh, 2);
   writeElementNodeData(name_mesh);
}

} //visu
