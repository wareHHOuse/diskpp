/*
 *       /\        Nicolas Pignet (C) 2024
 *      /__\       nicolas.pignet@enpc.fr
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

#pragma once

#include <array>
#include <vector>
#include <fstream>

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/geometry/geometry.hpp"
#include "diskpp/loaders/mesh_loader.hpp"

namespace disk
{

template<typename T, size_t N>
class poly_mesh_loader
{
    static_assert(N == 2 || N == 3, "Medit supports only 2D and 3D for the moment");
};

namespace poly
{

std::pair<size_t, std::vector<size_t>>
read_cell_line(std::ifstream& ifs)
{
    std::vector<size_t> elem;
    size_t              cell_id, nb_elem;

    ifs >> cell_id;
    ifs >> nb_elem;
    // std::cout << cell_id << ", " << nb_elem << std::endl;
    elem.reserve(nb_elem);
    for (size_t j = 0; j < nb_elem; j++)
    {
        size_t tmp;
        ifs >> tmp;
        elem.push_back(tmp);
    }

    return std::make_pair(cell_id, elem);
}

bool
read_cell_block(std::ifstream& ifs, const size_t nb_cells, std::vector<std::pair<size_t, std::vector<size_t>>>& cells)
{
    cells.clear();
    cells.reserve(nb_cells);

    for (size_t i = 0; i < nb_cells; i++)
        cells.push_back(read_cell_line(ifs));

    return true;
}

template<typename T>
bool
read_nodes_block(std::ifstream& ifs, const size_t nb_nodes, std::vector<std::pair<size_t, std::array<T, 3>>>& nodes)
{
    nodes.clear();
    nodes.reserve(nb_nodes);

    size_t           node_id;
    std::array<T, 3> coor;

    for (size_t i = 0; i < nb_nodes; i++)
    {
        ifs >> node_id;
        for (size_t j = 0; j < 3; j++)
        {
            ifs >> coor[j];
        }
        nodes.push_back(std::make_pair(node_id, coor));
    }

    return true;
}

std::pair<size_t, std::vector<size_t>>
read_grp_line(std::ifstream& ifs)
{
    std::vector<size_t> elem;
    std::string         name;
    size_t              cell_id, nb_elem;

    ifs >> name;
    ifs >> cell_id;
    ifs >> nb_elem;
    // std::cout << cell_id << ", " << nb_elem << std::endl;
    elem.reserve(nb_elem);
    for (size_t j = 0; j < nb_elem; j++)
    {
        size_t tmp;
        ifs >> tmp;
        elem.push_back(tmp);
    }

    return std::make_pair(cell_id, elem);
}

bool
read_grp_block(std::ifstream& ifs, const size_t nb_grp, std::vector<std::pair<size_t, std::vector<size_t>>>& grps)
{
    grps.clear();
    grps.reserve(nb_grp);

    for (size_t i = 0; i < nb_grp; i++)
        grps.push_back(read_grp_line(ifs));

    return true;
}
}

template<typename T>
class poly_mesh_loader<T, 2> : public mesh_loader<generic_mesh<T, 2>>
{
    typedef generic_mesh<T, 2>               mesh_type;
    typedef typename mesh_type::point_type   point_type;
    typedef typename mesh_type::node_type    node_type;
    typedef typename mesh_type::edge_type    edge_type;
    typedef typename mesh_type::surface_type surface_type;

    std::vector<point_type>                   m_points;
    std::vector<node_type>                    m_nodes;
    std::vector<std::pair<size_t, edge_type>> m_edges;

    std::vector<std::pair<size_t, std::array<T, 3>>>    vts;
    std::vector<std::pair<size_t, std::vector<size_t>>> vts_to_grp;
    std::vector<std::pair<size_t, std::vector<size_t>>> edges_to_vts;
    std::vector<std::pair<size_t, std::vector<size_t>>> edges_to_grp;
    std::vector<std::pair<size_t, std::vector<size_t>>> faces_to_edges;
    std::vector<std::pair<size_t, std::vector<size_t>>> faces_to_vts;
    std::vector<std::pair<size_t, std::vector<size_t>>> faces_to_grp;

    bool
    poly2d_read(const std::string& filename)
    {
        std::ifstream ifs(filename);
        std::string   keyword;
        size_t        dim, version, nb_elems;

        if (!ifs.is_open())
        {
            std::cout << "Error opening " << filename << std::endl;
            return false;
        }

        ifs >> keyword;
        if (keyword != "**BeginMesh")
        {
            std::cout << "Expected keyword \"**BeginMesh\"" << std::endl;
            return false;
        }

        while (keyword != "**EndMesh")
        {
            ifs >> keyword;
            if (keyword == "*Dimension")
            {
                ifs >> dim;
                if (this->verbose())
                    std::cout << "2D-Mesh" << std::endl;
                if (dim != 2)
                {
                    std::cout << "Expected dimension == 2 (here: " << dim << ")" << std::endl;
                    return false;
                }
            }
            else if (keyword == "*Version")
            {
                ifs >> version;
                if (this->verbose())
                    std::cout << "Mesh Version: " << version << std::endl;
            }
            else if (keyword == "*Nodes")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " nodes" << std::endl;

                poly::read_nodes_block(ifs, nb_elems, vts);
            }
            else if (keyword == "*Nodes->Groups")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " nodes->groups" << std::endl;

                poly::read_grp_block(ifs, nb_elems, vts_to_grp);
            }
            else if (keyword == "*Edges->Nodes")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " edges->nodes" << std::endl;

                poly::read_cell_block(ifs, nb_elems, edges_to_vts);
            }
            else if (keyword == "*Edges->Groups")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " edges->groups" << std::endl;

                poly::read_grp_block(ifs, nb_elems, edges_to_grp);
            }
            else if (keyword == "*Faces->Nodes")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " faces->nodes" << std::endl;

                poly::read_cell_block(ifs, nb_elems, faces_to_vts);
            }
            else if (keyword == "*Faces->Edges")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " faces->edges" << std::endl;

                poly::read_cell_block(ifs, nb_elems, faces_to_edges);
            }
            else if (keyword == "*Faces->Groups")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " faces->groups" << std::endl;

                poly::read_grp_block(ifs, nb_elems, faces_to_grp);
            }
            else if (keyword == "**EndMesh")
            {
                return true;
            }
            else
            {
                std::cout << "Unexpected keyword: " << keyword << std::endl;
                return false;
            }
        }

        return true;
    }

  public:
    static const char constexpr* expected_extension = "poly2d";
    poly_mesh_loader()                              = default;

    bool
    read_mesh(const std::string& s)
    {
        if (this->verbose())
            std::cout << " *** READING POLY2D MESH ***" << std::endl;
        return poly2d_read(s);
    }

    bool
    populate_mesh(mesh_type& msh)
    {
        if (this->verbose())
            std::cout << " *** POPULATE POLY2D MESH ***" << std::endl;

        auto storage = msh.backend_storage();

        // Create nodes and vertices
        m_points.reserve(vts.size());
        m_nodes.reserve(vts.size());
        for (auto& [id, coor] : vts)
        {
            m_points.push_back(point_type({coor[0], coor[1]}));
            m_nodes.push_back(node_type(disk::point_identifier<2>(id)));
        }

        // create edges
        m_edges.reserve(edges_to_vts.size());
        for (auto& [id, vts] : edges_to_vts)
        {
            size_t v1 = vts[0], v2 = vts[1];

            if (v1 > v2)
                std::swap(v1, v2);

            auto e = edge_type({typename node_type::id_type(v1), typename node_type::id_type(v2)});

            m_edges.push_back(std::make_pair(id, e));
        }

        std::vector<size_t> conv_table;
        /* Sort the edges in lexicographical order, remember their original
         * position to convert the pointers in the faces */
        auto comp_edges = [](const std::pair<size_t, edge_type>& e1, const std::pair<size_t, edge_type>& e2)
        { return e1.second < e2.second; };
        std::sort(m_edges.begin(), m_edges.end(), comp_edges);

        std::vector<edge_type> edges;
        edges.reserve(m_edges.size());
        conv_table.resize(m_edges.size());
        for (size_t i = 0; i < m_edges.size(); i++)
        {
            conv_table.at(m_edges.at(i).first) = i; /* Make ptr conversion table */
            edges.push_back(m_edges.at(i).second);
        }

        /* Convert the edge pointers in the face data */
        for (auto& fe : faces_to_edges)
        {
            for (auto& ptr : fe.second)
                ptr = conv_table[ptr];
        }

        /* Detect which ones are boundary edges */
        storage->boundary_info.resize(edges.size());
        for (auto& [id, edges] : edges_to_grp)
        {
            for (auto& edge_id : edges)
            {
                boundary_descriptor bi(id, true);
                storage->boundary_info.at(conv_table[edge_id]) = bi;
            }
        }

        /* Sort in lexicographical order and remember original position */
        auto comp_vecs = [](const std::pair<size_t, std::vector<size_t>>& e1,
                            const std::pair<size_t, std::vector<size_t>>& e2) { return e1.second < e2.second; };
        std::sort(faces_to_edges.begin(), faces_to_edges.end(), comp_vecs);

        std::vector<surface_type> faces;
        faces.reserve(faces_to_edges.size());
        conv_table.resize(faces_to_edges.size());

        for (size_t i = 0; i < faces_to_edges.size(); i++)
        {
            auto         fe = faces_to_edges.at(i);
            surface_type s(convert_to<typename edge_type::id_type>(fe.second));
            s.set_point_ids(convert_to<disk::point_identifier<2>>(faces_to_vts.at(fe.first).second));
            faces.push_back(s);
            conv_table[fe.first] = i;
        }
        /* Now the faces are in their place and have correct ptrs */

        storage->points   = std::move(m_points);
        storage->nodes    = std::move(m_nodes);
        storage->edges    = std::move(edges);
        storage->surfaces = std::move(faces);

        return true;
    }
};

template<typename T>
class poly_mesh_loader<T, 3> : public mesh_loader<generic_mesh<T, 3>>
{
    typedef generic_mesh<T, 3>               mesh_type;
    typedef typename mesh_type::point_type   point_type;
    typedef typename mesh_type::node_type    node_type;
    typedef typename mesh_type::edge_type    edge_type;
    typedef typename mesh_type::surface_type surface_type;
    typedef typename mesh_type::volume_type  volume_type;

    std::vector<point_type>                   m_points;
    std::vector<node_type>                    m_nodes;
    std::vector<std::pair<size_t, edge_type>> m_edges;

    std::vector<std::pair<size_t, std::array<T, 3>>>    vts;
    std::vector<std::pair<size_t, std::vector<size_t>>> vts_to_grp;
    std::vector<std::pair<size_t, std::vector<size_t>>> edges_to_vts;
    std::vector<std::pair<size_t, std::vector<size_t>>> edges_to_grp;
    std::vector<std::pair<size_t, std::vector<size_t>>> faces_to_edges;
    std::vector<std::pair<size_t, std::vector<size_t>>> faces_to_vts;
    std::vector<std::pair<size_t, std::vector<size_t>>> faces_to_grp;
    std::vector<std::pair<size_t, std::vector<size_t>>> vols_to_edges;
    std::vector<std::pair<size_t, std::vector<size_t>>> vols_to_faces;
    std::vector<std::pair<size_t, std::vector<size_t>>> vols_to_vts;
    std::vector<std::pair<size_t, std::vector<size_t>>> vols_to_grp;

    bool
    poly3d_read(const std::string& filename)
    {
        std::ifstream ifs(filename);
        std::string   keyword;
        size_t        dim, version, nb_elems;

        if (!ifs.is_open())
        {
            std::cout << "Error opening " << filename << std::endl;
            return false;
        }

        ifs >> keyword;
        if (keyword != "**BeginMesh")
        {
            std::cout << "Expected keyword \"**BeginMesh\"" << std::endl;
            return false;
        }

        while (keyword != "**EndMesh")
        {
            ifs >> keyword;
            if (keyword == "*Dimension")
            {
                ifs >> dim;
                if (this->verbose())
                    std::cout << "3D-Mesh" << std::endl;
                if (dim != 3)
                {
                    std::cout << "Expected dimension == 3 (here: " << dim << ")" << std::endl;
                    return false;
                }
            }
            else if (keyword == "*Version")
            {
                ifs >> version;
                if (this->verbose())
                    std::cout << "Mesh Version: " << version << std::endl;
            }
            else if (keyword == "*Nodes")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " nodes" << std::endl;

                poly::read_nodes_block(ifs, nb_elems, vts);
            }
            else if (keyword == "*Nodes->Groups")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " nodes->groups" << std::endl;

                poly::read_grp_block(ifs, nb_elems, vts_to_grp);
            }
            else if (keyword == "*Edges->Nodes")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " edges->nodes" << std::endl;

                poly::read_cell_block(ifs, nb_elems, edges_to_vts);
            }
            else if (keyword == "*Edges->Groups")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " edges->groups" << std::endl;

                poly::read_grp_block(ifs, nb_elems, edges_to_grp);
            }
            else if (keyword == "*Faces->Nodes")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " faces->nodes" << std::endl;

                poly::read_cell_block(ifs, nb_elems, faces_to_vts);
            }
            else if (keyword == "*Faces->Edges")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " faces->edges" << std::endl;

                poly::read_cell_block(ifs, nb_elems, faces_to_edges);
            }
            else if (keyword == "*Faces->Groups")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " faces->groups" << std::endl;

                poly::read_grp_block(ifs, nb_elems, faces_to_grp);
            }
            else if (keyword == "*Volumes->Nodes")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " volumes->nodes" << std::endl;

                poly::read_cell_block(ifs, nb_elems, vols_to_vts);
            }
            else if (keyword == "*Volumes->Edges")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " volumes->edges" << std::endl;

                poly::read_cell_block(ifs, nb_elems, vols_to_edges);
            }
            else if (keyword == "*Volumes->Faces")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " volumes->faces" << std::endl;

                poly::read_cell_block(ifs, nb_elems, vols_to_faces);
            }
            else if (keyword == "*Volumes->Groups")
            {
                ifs >> nb_elems;
                if (this->verbose())
                    std::cout << "Reading " << nb_elems << " volumes->groups" << std::endl;

                poly::read_grp_block(ifs, nb_elems, vols_to_grp);
            }
            else if (keyword == "**EndMesh")
            {
                return true;
            }
            else
            {
                std::cout << "Unexpected keyword: " << keyword << std::endl;
                return false;
            }
        }

        return true;
    }

  public:
    static const char constexpr* expected_extension = "poly3d";
    poly_mesh_loader()                              = default;

    bool
    read_mesh(const std::string& s)
    {
        if (this->verbose())
            std::cout << " *** READING POLY3D MESH ***" << std::endl;
        return poly3d_read(s);
    }

    bool
    populate_mesh(mesh_type& msh)
    {
        if (this->verbose())
            std::cout << " *** POPULATE POLY3D MESH ***" << std::endl;

        auto storage = msh.backend_storage();

        // Create nodes and vertices
        m_points.reserve(vts.size());
        m_nodes.reserve(vts.size());
        for (auto& [id, coor] : vts)
        {
            m_points.push_back(point_type({coor[0], coor[1], coor[2]}));
            m_nodes.push_back(node_type(disk::point_identifier<3>(id)));
        }

        // create edges
        m_edges.reserve(edges_to_vts.size());
        for (auto& [id, vts] : edges_to_vts)
        {
            size_t v1 = vts[0], v2 = vts[1];

            if (v1 > v2)
                std::swap(v1, v2);

            auto e = edge_type({typename node_type::id_type(v1), typename node_type::id_type(v2)});

            m_edges.push_back(std::make_pair(id, e));
        }

        std::vector<size_t> conv_table;
        /* Sort the edges in lexicographical order, remember their original
         * position to convert the pointers in the faces */
        auto comp_edges = [](const std::pair<size_t, edge_type>& e1, const std::pair<size_t, edge_type>& e2)
        { return e1.second < e2.second; };
        std::sort(m_edges.begin(), m_edges.end(), comp_edges);

        std::vector<edge_type> edges;
        edges.reserve(m_edges.size());
        conv_table.resize(m_edges.size());
        for (size_t i = 0; i < m_edges.size(); i++)
        {
            conv_table.at(m_edges.at(i).first) = i; /* Make ptr conversion table */
            edges.push_back(m_edges.at(i).second);
        }

        /* Convert the edge pointers in the face data */
        for (auto& fe : faces_to_edges)
        {
            for (auto& ptr : fe.second)
                ptr = conv_table[ptr];
        }

        /* Sort in lexicographical order and remember original position */
        auto comp_vecs = [](const std::pair<size_t, std::vector<size_t>>& e1,
                            const std::pair<size_t, std::vector<size_t>>& e2) { return e1.second < e2.second; };
        std::sort(faces_to_edges.begin(), faces_to_edges.end(), comp_vecs);

        std::vector<surface_type> faces;
        faces.reserve(faces_to_edges.size());
        conv_table.resize(faces_to_edges.size());

        for (size_t i = 0; i < faces_to_edges.size(); i++)
        {
            auto         fe = faces_to_edges.at(i);
            surface_type s(convert_to<typename edge_type::id_type>(fe.second));
            s.set_point_ids(convert_to<disk::point_identifier<3>>(faces_to_vts.at(fe.first).second));
            faces.push_back(s);
            conv_table[fe.first] = i;
        }
        /* Now the faces are in their place and have correct ptrs */

        /* Detect which ones are boundary edges */
        storage->boundary_info.resize(faces.size());
        for (auto& [id, faces] : faces_to_grp)
        {
            for (auto& face_id : faces)
            {
                boundary_descriptor bi(id, true);
                storage->boundary_info.at(conv_table.at(face_id)) = bi;
            }
        }

        /* Convert the face pointers in the volume data */
        for (auto& vf : vols_to_faces)
        {
            for (auto& ptr : vf.second)
                ptr = conv_table[ptr];
        }

        /* Sort volume data */
        std::sort(vols_to_faces.begin(), vols_to_faces.end(), comp_vecs);

        std::vector<volume_type> volumes;
        volumes.reserve(vols_to_faces.size());

        for (size_t i = 0; i < vols_to_faces.size(); i++)
        {
            auto        vf = vols_to_faces[i];
            volume_type v(convert_to<typename surface_type::id_type>(vf.second));
            v.set_point_ids(convert_to<disk::point_identifier<3>>(vols_to_vts.at(vf.first).second));
            volumes.push_back(v);
        }

        storage->points   = std::move(m_points);
        storage->nodes    = std::move(m_nodes);
        storage->edges    = std::move(edges);
        storage->surfaces = std::move(faces);
        storage->volumes  = std::move(volumes);

        return true;
    }
};

} // end namespace disk
