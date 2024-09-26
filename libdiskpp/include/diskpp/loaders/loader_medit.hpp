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

#include <vector>
#include <fstream>
#include <array>

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/geometry/geometry.hpp"
#include "diskpp/loaders/mesh_loader.hpp"

namespace disk
{

template<typename T, size_t N>
class medit_mesh_loader
{
    static_assert(N == 2 || N == 3, "Medit supports only 2D and 3D for the moment");
};

template<typename T>
class medit_mesh_loader<T, 2> : public mesh_loader<generic_mesh<T, 2>>
{
    typedef generic_mesh<T, 2>               mesh_type;
    typedef typename mesh_type::point_type   point_type;
    typedef typename mesh_type::node_type    node_type;
    typedef typename mesh_type::edge_type    edge_type;
    typedef typename mesh_type::surface_type surface_type;

    struct medit2d_poly
    {
        std::vector<size_t>             nodes;
        size_t                          id;
        std::set<std::array<size_t, 2>> attached_edges;

        bool
        operator<(const medit2d_poly& other)
        {
            return nodes < other.nodes;
        }
    };

    std::vector<point_type>                                    m_points;
    std::vector<node_type>                                     m_nodes;
    std::vector<std::array<ident_raw_t, 2>>                    m_edges;
    std::vector<medit2d_poly>                                  m_polys;
    std::vector<std::pair<std::array<ident_raw_t, 2>, size_t>> m_boundary_edges;

    bool
    medit_read_vertices(std::ifstream& ifs)
    {
        size_t elements_to_read;
        T      x, y, z;
        T      id;

        ifs >> elements_to_read;

        if (this->verbose())
            std::cout << "Attempting to read " << elements_to_read << " points" << std::endl;

        m_points.reserve(elements_to_read);
        m_nodes.reserve(elements_to_read);

        for (size_t i = 0; i < elements_to_read; i++)
        {
            ifs >> x >> y >> z >> id;
            m_points.push_back(point_type{x, y});
            m_nodes.push_back(node_type(disk::point_identifier<2>(i)));
        }

        return true;
    }

    bool
    medit_read_polygons(std::ifstream& ifs, size_t polynum)
    {
        size_t elements_to_read;

        ifs >> elements_to_read;
        if (this->verbose())
            std::cout << "Reading " << elements_to_read << " " << polynum << "-angles" << std::endl;

        for (size_t i = 0; i < elements_to_read; i++)
        {
            medit2d_poly             p;
            std::vector<ident_raw_t> nodes(polynum + 1, 0);

            for (size_t j = 0; j < polynum; j++)
            {
                ident_raw_t val;
                ifs >> val;
                p.nodes.push_back(val - 1);
                nodes[j] = val - 1;
            }
            nodes[polynum] = nodes[0];

            ifs >> p.id;

            m_polys.push_back(p);

            size_t p_id = m_polys.size();

            // We have too create edges
            for (size_t j = 0; j < polynum; j++)
            {
                std::array<ident_raw_t, 2> b_edge = {nodes[j], nodes[j + 1]};
                assert(b_edge[0] != b_edge[1]);
                if (b_edge[0] > b_edge[1])
                    std::swap(b_edge[0], b_edge[1]);

                m_edges.push_back(b_edge);
                m_polys.at(p_id - 1).attached_edges.insert({b_edge[0], b_edge[1]});
            }
        }

        return true;
    }

    bool
    medit_read_boundary_edges(std::ifstream& ifs)
    {
        size_t elements_to_read;

        ifs >> elements_to_read;
        if (this->verbose())
            std::cout << "Reading " << elements_to_read << " boundary edges" << std::endl;

        m_boundary_edges.reserve(elements_to_read);

        for (size_t i = 0; i < elements_to_read; i++)
        {
            std::array<ident_raw_t, 2> b_edge;
            ifs >> b_edge[0];
            b_edge[0] -= 1;
            ifs >> b_edge[1];
            b_edge[1] -= 1;

            assert(b_edge[0] != b_edge[1]);

            if (b_edge[0] > b_edge[1])
                std::swap(b_edge[0], b_edge[1]);

            std::array<ident_raw_t, 2> bnd = {b_edge[0], b_edge[1]};

            size_t b_id;
            ifs >> b_id;

            m_boundary_edges.push_back(std::make_pair(bnd, b_id));
        }

        return true;
    }

    bool
    face_unique(const std::array<ident_raw_t, 2>& f1, const std::array<ident_raw_t, 2>& f2)
    {
        if (f1[0] == f2[0])
        {
            if (f1[1] == f2[1])
            {
                return false;
            }
            else
                return true;
        }
        else
            return true;
    }

    bool
    medit_read(const std::string& filename)
    {
        std::ifstream ifs(filename);
        std::string   keyword;

        if (!ifs.is_open())
        {
            std::cout << "Error opening " << filename << std::endl;
            return false;
        }

        ifs >> keyword;
        if (keyword != "MeshVersionFormatted")
        {
            std::cout << "Expected keyword \"MeshVersionFormatted\"" << std::endl;
            return false;
        }

        size_t format;
        ifs >> format;

        if (format != 2)
        {
            std::cout << "Expected format 2 (here: " << format << ")" << std::endl;
            return false;
        }

        ifs >> keyword;
        if (keyword != "Dimension")
        {
            std::cout << "Expected keyword \"Dimension\"" << std::endl;
            return false;
        }

        size_t dim;
        ifs >> dim;

        if (dim != 3)
        {
            std::cout << "Expected dimension >=2 (here: " << dim << ")" << std::endl;
            return false;
        }

        ifs >> keyword;
        while (keyword != "End")
        {
            if (keyword == "Vertices")
            {
                medit_read_vertices(ifs);
            }
            else if (keyword == "Triangles")
            {
                medit_read_polygons(ifs, 3);
            }
            else if (keyword == "Quadrilaterals")
            {
                medit_read_polygons(ifs, 4);
            }
            else if (keyword == "Pentagons")
            {
                medit_read_polygons(ifs, 5);
            }
            else if (keyword == "Hexagons")
            {
                medit_read_polygons(ifs, 6);
            }
            else if (keyword == "Edges")
            {
                m_boundary_edges.clear();
                medit_read_boundary_edges(ifs);
            }
            else
            {
                std::cout << "Error parsing Medit file" << std::endl;
                return false;
            }

            ifs >> keyword;
        }

        ifs.close();
        return true;
    }

  public:
    static const char constexpr* expected_extension = "medit2d";

    medit_mesh_loader() = default;

    bool
    read_mesh(const std::string& s)
    {
        if (this->verbose())
            std::cout << " *** READING MEDIT 2D MESH ***" << std::endl;

        return medit_read(s);
    }

    bool
    populate_mesh(mesh_type& msh)
    {
        if (this->verbose())
            std::cout << " *** POPULATING MEDIT MESH ***" << std::endl;
        auto storage = msh.backend_storage();

        /* Points */
        storage->points = std::move(m_points);
        storage->nodes  = std::move(m_nodes);

        /* Edges */
        /* Make the vector containing the edges */
        std::vector<edge_type> edges;
        edges.reserve(m_edges.size());
        size_t nb_edge(0);
        for (size_t i = 0; i < m_edges.size(); i++)
        {
            bool unique = true;
            for (size_t j = 0; j < i; j++)
            {
                if (!face_unique(m_edges[i], m_edges[j]))
                {
                    unique = false;
                    break;
                }
            }

            if (unique)
            {
                assert(m_edges[i][0] < m_edges[i][1]);
                auto node1 = typename node_type::id_type(m_edges[i][0]);
                auto node2 = typename node_type::id_type(m_edges[i][1]);

                auto e = edge_type(node1, node2);

                /* Next line not necessary anymore, see generic_element<DIM, DIM-1> */
                // e.set_point_ids(m_edges[i].begin(), m_edges[i].begin() + 2); /* XXX: crap */
                edges.push_back(e);
                nb_edge++;
            }
        }
        /* Sort them */
        edges.resize(nb_edge);
        std::sort(edges.begin(), edges.end());

        /* Detect which ones are boundary edges */
        storage->boundary_info.resize(edges.size());
        for (size_t i = 0; i < m_boundary_edges.size(); i++)
        {
            assert(m_boundary_edges[i].first[0] < m_boundary_edges[i].first[1]);
            auto node1 = typename node_type::id_type(m_boundary_edges[i].first[0]);
            auto node2 = typename node_type::id_type(m_boundary_edges[i].first[1]);

            auto e = edge_type(node1, node2);

            auto position = find_element_id(edges.begin(), edges.end(), e);

            if (position.first == false)
            {
                std::cout << "Bad bug at " << __FILE__ << "(" << __LINE__ << ")" << std::endl;
                return false;
            }
            else
            {
                boundary_descriptor bi(m_boundary_edges[i].second, true);
                storage->boundary_info.at(position.second) = bi;
            }
        }

        storage->edges = std::move(edges);

        /* Surfaces */
        std::vector<surface_type> surfaces;
        surfaces.reserve(m_polys.size());

        for (auto& p : m_polys)
        {
            std::vector<typename edge_type::id_type> surface_edges;
            for (auto& e : p.attached_edges)
            {
                assert(e[0] < e[1]);
                auto n1 = typename node_type::id_type(e[0]);
                auto n2 = typename node_type::id_type(e[1]);

                edge_type edge(n1, n2);
                auto      edge_id = find_element_id(storage->edges.begin(), storage->edges.end(), edge);
                if (!edge_id.first)
                {
                    std::cout << "Bad bug at " << __FILE__ << "(" << __LINE__ << ")" << std::endl;
                    return false;
                }

                surface_edges.push_back(edge_id.second);
            }
            auto surface = surface_type(surface_edges);
            surface.set_point_ids(p.nodes.begin(), p.nodes.end()); /* XXX: crap */
            surfaces.push_back(surface);
        }

        std::sort(surfaces.begin(), surfaces.end());
        storage->surfaces = std::move(surfaces);

        return true;
    }
};

template<typename T>
class medit_mesh_loader<T, 3> : public mesh_loader<generic_mesh<T, 3>>
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
    std::vector<std::array<size_t, 2>>        m_boundary_edges;

    std::vector<std::pair<size_t, std::vector<size_t>>> vol_to_faces;
    std::vector<std::vector<size_t>>                    vol_to_vts;
    std::vector<std::pair<size_t, std::vector<size_t>>> faces_to_edges;
    std::vector<std::vector<size_t>>                    faces_to_vts;

    bool
    expect(std::ifstream& ifs, const std::string& str)
    {
        std::string keyword;
        ifs >> keyword;
        if (keyword != str)
        {
            std::cout << "Expected keyword \"" << str << "\"" << std::endl;
            std::cout << "Found \"" << keyword << "\"" << std::endl;
            return false;
        }

        return true;
    }

    std::vector<size_t>
    read_medit3d_line(std::ifstream& ifs)
    {
        std::vector<size_t> ret;

        size_t num_entries;
        ifs >> num_entries;
        for (size_t j = 0; j < num_entries; j++)
        {
            size_t temp;
            ifs >> temp;
            ret.push_back(temp);
        }

        return ret;
    }

    bool
    medit3d_read(const std::string& filename)
    {
        std::ifstream ifs(filename);
        std::string   keyword;
        size_t        lines_to_read;

        if (!ifs.is_open())
        {
            std::cout << "Error opening " << filename << std::endl;
            return false;
        }

        ifs >> keyword;
        if (keyword != "MeshVersionFormatted")
        {
            std::cout << "Expected keyword \"MeshVersionFormatted\"" << std::endl;
            return false;
        }

        size_t format;
        ifs >> format;

        if (format != 2)
        {
            std::cout << "Expected format 2 (here: " << format << ")" << std::endl;
            return false;
        }

        ifs >> keyword;
        if (keyword != "Dimension")
        {
            std::cout << "Expected keyword \"Dimension\"" << std::endl;
            return false;
        }

        size_t dim;
        ifs >> dim;

        if (dim != 3)
        {
            std::cout << "Expected dimension == 3 (here: " << dim << ")" << std::endl;
            return false;
        }

        if (!expect(ifs, "Vertices"))
            return false;

        ifs >> lines_to_read;

        if (this->verbose())
            std::cout << "About to read " << lines_to_read << " points" << std::endl;

        m_points.reserve(lines_to_read);
        m_nodes.reserve(lines_to_read);
        for (size_t i = 0; i < lines_to_read; i++)
        {
            T x, y, z;
            ifs >> x >> y >> z;
            m_points.push_back(point_type({x, y, z}));
            m_nodes.push_back(node_type(disk::point_identifier<3>(i)));
        }

        /* Volume to face data */
        if (!expect(ifs, "Volumes->Faces"))
            return false;

        ifs >> lines_to_read;

        if (this->verbose())
            std::cout << "About to read " << lines_to_read << " entries" << std::endl;

        for (size_t i = 0; i < lines_to_read; i++)
        {
            auto vol_faces = read_medit3d_line(ifs);
            vol_to_faces.push_back(std::make_pair(i, std::move(vol_faces)));
        }

        /* Volume to vertices data */
        if (!expect(ifs, "Volumes->Vertices"))
            return false;

        ifs >> lines_to_read;

        if (this->verbose())
            std::cout << "About to read " << lines_to_read << " entries" << std::endl;

        for (size_t i = 0; i < lines_to_read; i++)
        {
            auto vol_vts = read_medit3d_line(ifs);
            vol_to_vts.push_back(std::move(vol_vts));
        }

        /* Faces to edges data */
        if (!expect(ifs, "Faces->Edges"))
            return false;

        ifs >> lines_to_read;

        if (this->verbose())
            std::cout << "About to read " << lines_to_read << " entries" << std::endl;

        for (size_t i = 0; i < lines_to_read; i++)
        {
            auto faces_edges = read_medit3d_line(ifs);
            faces_to_edges.push_back(std::make_pair(i, std::move(faces_edges)));
        }

        /* Faces to vertices data */
        if (!expect(ifs, "Faces->Vertices"))
            return false;

        ifs >> lines_to_read;

        if (this->verbose())
            std::cout << "About to read " << lines_to_read << " entries" << std::endl;

        for (size_t i = 0; i < lines_to_read; i++)
        {
            auto faces_vts = read_medit3d_line(ifs);
            faces_to_vts.push_back(std::move(faces_vts));
        }

        /* Faces to cv data */
        if (!expect(ifs, "Faces->Control"))
            return false;

        if (!expect(ifs, "volumes"))
            return false;

        ifs >> lines_to_read;

        if (this->verbose())
            std::cout << "About to read " << lines_to_read << " entries" << std::endl;

        for (size_t i = 0; i < lines_to_read; i++)
        {
            size_t num, id;
            ifs >> num >> id;
            m_boundary_edges.push_back({num, id});
        }

        /* Edges data */
        if (!expect(ifs, "Edges"))
            return false;

        ifs >> lines_to_read;

        if (this->verbose())
            std::cout << "About to read " << lines_to_read << " entries" << std::endl;

        for (size_t i = 0; i < lines_to_read; i++)
        {
            size_t v1, v2;
            ifs >> v1 >> v2;

            if (v1 > v2)
                std::swap(v1, v2);

            auto e = edge_type({typename node_type::id_type(v1), typename node_type::id_type(v2)});

            m_edges.push_back(std::make_pair(i, e));
        }

        return true;
    }

  public:

    static const char constexpr* expected_extension = "medit3d";

    medit_mesh_loader() = default;

    bool
    read_mesh(const std::string& s)
    {
        if (this->verbose())
            std::cout << " *** READING MEDIT 3D MESH ***" << std::endl;
        return medit3d_read(s);
    }

    bool
    populate_mesh(mesh_type& msh)
    {
        auto storage = msh.backend_storage();

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
            conv_table[m_edges[i].first] = i; /* Make ptr conversion table */
            edges.push_back(m_edges[i].second);
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
            auto         fe = faces_to_edges[i];
            surface_type s(convert_to<typename edge_type::id_type>(fe.second));
            s.set_point_ids(convert_to<disk::point_identifier<3>>(faces_to_vts.at(fe.first)));
            faces.push_back(s);
            conv_table[fe.first] = i;
        }
        /* Now the faces are in their place and have correct ptrs */

        /* Convert the face pointers in the volume data */
        for (auto& vf : vol_to_faces)
        {
            for (auto& ptr : vf.second)
                ptr = conv_table[ptr];
        }

        /* Sort volume data */
        std::sort(vol_to_faces.begin(), vol_to_faces.end(), comp_vecs);

        std::vector<volume_type> volumes;
        volumes.reserve(vol_to_faces.size());

        for (size_t i = 0; i < vol_to_faces.size(); i++)
        {
            auto        vf = vol_to_faces[i];
            volume_type v(convert_to<typename surface_type::id_type>(vf.second));
            v.set_point_ids(convert_to<disk::point_identifier<3>>(vol_to_vts.at(vf.first)));
            volumes.push_back(v);
        }

        storage->points   = std::move(m_points);
        storage->nodes    = std::move(m_nodes);
        storage->edges    = std::move(edges);
        storage->surfaces = std::move(faces);
        storage->volumes  = std::move(volumes);

        storage->boundary_info.resize(storage->surfaces.size());
        for (size_t i = 0; i < m_boundary_edges.size(); i++)
        {
            boundary_descriptor bi(m_boundary_edges[i][1], true);
            storage->boundary_info[conv_table[m_boundary_edges[i][0]]] = bi;
        }

        return false;
    }
};

} // end namespace disk
