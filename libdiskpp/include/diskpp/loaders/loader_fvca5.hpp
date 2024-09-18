/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023, 2024
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */
/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2020, 2021
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group
 */
/*
 *       /\        Matteo Cicuttin (C) 2016-2019
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

#include <tuple>
#include <vector>
#include <array>
#include <fstream>

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/geometry/geometry.hpp"
#include "diskpp/loaders/mesh_loader.hpp"

namespace disk
{

template<typename T, size_t N>
class fvca5_mesh_loader
{
    static_assert(N == 2, "FVCA5 is a 2D-only mesh format");
};

template<typename T>
class fvca5_mesh_loader<T, 2> : public mesh_loader<generic_mesh<T, 2>>
{
    typedef generic_mesh<T, 2>               mesh_type;
    typedef typename mesh_type::point_type   point_type;
    typedef typename mesh_type::node_type    node_type;
    typedef typename mesh_type::edge_type    edge_type;
    typedef typename mesh_type::surface_type surface_type;

    struct fvca5_poly
    {
        std::vector<size_t>             nodes;
        std::set<std::array<size_t, 2>> attached_edges;
        size_t                          domain_id;

        fvca5_poly() : domain_id(0) {}

        bool
        operator<(const fvca5_poly& other)
        {
            return nodes < other.nodes;
        }
    };

    std::ifstream ifs;
    size_t        current_line;

    std::vector<point_type>                 m_points;
    std::vector<fvca5_poly>                 m_polys;
    std::vector<std::array<ident_raw_t, 2>> m_boundary_edges;
    std::vector<std::array<ident_raw_t, 4>> m_edges;

    bool
    get_next_line(std::istringstream& iss)
    {
        std::string line;
        getline(ifs, line);
        iss = std::istringstream(line);
        current_line++;
        return true;
    }

    bool
    get_next_line(std::string& line)
    {
        getline(ifs, line);
        current_line++;
        return true;
    }

    bool
    failure(std::istringstream& ifs)
    {
        return ifs.rdstate() & std::ifstream::failbit;
    }

    bool
    fvca5_parse_vertices()
    {
        std::istringstream iss;
        get_next_line(iss);

        size_t num_vertices;
        iss >> num_vertices;
        if (failure(iss))
        {
            std::cout << "Line " << current_line << ": Error while parsing ";
            std::cout << "'vertices' block (number of vertices)" << std::endl;
            return false;
        }

        if (this->verbose())
            std::cout << "Reading " << num_vertices << " vertices" << std::endl;

        for (size_t i = 0; i < num_vertices; i++)
        {
            get_next_line(iss);
            T x, y;
            iss >> x >> y;
            if (failure(iss))
            {
                std::cout << "Line " << current_line << ": Error while parsing ";
                std::cout << "'vertices' block (vertex " << i << ")" << std::endl;
                return false;
            }

            m_points.push_back(point_type{x, y});
        }

        return true;
    }

    bool
    fvca5_parse_polygons(size_t polynum)
    {
        std::istringstream iss;
        get_next_line(iss);

        size_t num_polygons;

        iss >> num_polygons;
        if (failure(iss))
        {
            std::cout << "Line " << current_line << ": Error while parsing ";
            std::cout << "'" << polynum << "-angles' block (number of polygons)";
            std::cout << std::endl;
            return false;
        }

        if (this->verbose())
            std::cout << "Reading " << num_polygons << " " << polynum << "-angles" << std::endl;

        for (size_t i = 0; i < num_polygons; i++)
        {
            get_next_line(iss);

            fvca5_poly p;

            for (size_t j = 0; j < polynum; j++)
            {
                ident_raw_t val;
                iss >> val;
                p.nodes.push_back(val - 1);
            }

            if (failure(iss))
            {
                std::cout << "Line " << current_line << ": Error while parsing '";
                std::cout << polynum << "-angles' block (polygon ";
                std::cout << i << ")" << std::endl;
                return false;
            }

            // JMLC extension
            iss >> std::ws;
            if (std::isdigit(iss.peek()))
            {
                std::cout << "d_id" << std::endl;
                iss >> p.domain_id;
            }

            m_polys.push_back(p);
        }

        return true;
    }

    bool
    fvca5_parse_boundary_edges()
    {
        std::istringstream iss;
        get_next_line(iss);

        size_t num_edges;

        iss >> num_edges;
        if (failure(iss))
        {
            std::cout << "Line " << current_line << ": Error while parsing ";
            std::cout << "'edges of the boundary' block (number of edges)";
            std::cout << std::endl;
            return false;
        }

        if (this->verbose())
            std::cout << "Reading " << num_edges << " boundary edges" << std::endl;

        m_boundary_edges.reserve(num_edges);

        for (size_t i = 0; i < num_edges; i++)
        {
            get_next_line(iss);
            ident_raw_t node_a, node_b;
            iss >> node_a >> node_b;

            if (failure(iss))
            {
                std::cout << "Line " << current_line << ": Error while parsing ";
                std::cout << "'edges of the boundary' block (edge " << i << ")";
                std::cout << std::endl;
                return false;
            }

            if (node_a < 1 or node_b < 1)
            {
                std::cout << "Line " << current_line << ": FVCA5 format ";
                std::cout << "expects 1-based indices" << std::endl;
                return false;
            }

            node_a -= 1;
            node_b -= 1;

            if (node_a == node_b)
            {
                std::cout << "Line " << current_line << ": Edge starting and ";
                std::cout << "finishing on the same node" << std::endl;
                return false;
            }

            if (node_a > node_b)
                std::swap(node_a, node_b);

            m_boundary_edges.push_back({node_a, node_b});
        }

        return true;
    }

    bool
    fvca5_parse_all_edges()
    {
        std::istringstream iss;
        get_next_line(iss);

        size_t num_edges;

        iss >> num_edges;
        if (failure(iss))
        {
            std::cout << "Line " << current_line << ": Error while parsing ";
            std::cout << "'all edges' block (number of edges)";
            std::cout << std::endl;
            return false;
        }

        if (this->verbose())
            std::cout << "Reading " << num_edges << " edges" << std::endl;

        for (size_t i = 0; i < num_edges; i++)
        {
            get_next_line(iss);
            ident_raw_t node_a, node_b, neigh_a, neigh_b;

            iss >> node_a >> node_b >> neigh_a >> neigh_b;
            if (failure(iss))
            {
                std::cout << "Line " << current_line << ": Error while parsing ";
                std::cout << "'all edges' block (edge " << i << ")";
                std::cout << std::endl;
                return false;
            }

            if (node_a < 1 or node_b < 1)
            {
                std::cout << "Line " << current_line << ": FVCA5 format ";
                std::cout << "expects 1-based indices" << std::endl;
                return false;
            }

            node_a -= 1;
            node_b -= 1;

            if (node_a == node_b)
            {
                std::cout << "Line " << current_line << ": Edge starting and ";
                std::cout << "finishing on the same node" << std::endl;
                return false;
            }

            if (node_a > node_b)
                std::swap(node_a, node_b);

            if (neigh_a > 0)
                m_polys.at(neigh_a - 1).attached_edges.insert({node_a, node_b});

            if (neigh_b > 0)
                m_polys.at(neigh_b - 1).attached_edges.insert({node_a, node_b});

            // JMLC extension
            iss >> std::ws;
            if (std::isdigit(iss.peek()))
            {
                std::cout << "b_id" << std::endl;
                // iss >> p.domain_id;
            }

            m_edges.push_back({node_a, node_b, neigh_a, neigh_b});
        }

        return true;
    }

    bool
    fvca5_parse_block()
    {
        std::string line;
        get_next_line(line);

        if (std::regex_match(line, std::regex("^\\s*$")))
            return true;

        if (std::regex_match(line, std::regex("^\\s*vertices\\s*$")))
            return fvca5_parse_vertices();

        if (std::regex_match(line, std::regex("^\\s*triangles\\s*$")))
            return fvca5_parse_polygons(3);

        if (std::regex_match(line, std::regex("^\\s*quadrangles\\s*$")))
            return fvca5_parse_polygons(4);

        if (std::regex_match(line, std::regex("^\\s*pentagons\\s*$")))
            return fvca5_parse_polygons(5);

        if (std::regex_match(line, std::regex("^\\s*hexagons\\s*$")))
            return fvca5_parse_polygons(6);

        if (std::regex_match(line, std::regex("^\\s*ennagons\\s*$")))
            return fvca5_parse_polygons(7);

        if (std::regex_match(line, std::regex("^\\s*ettagons\\s*$")))
            return fvca5_parse_polygons(8);

        if (std::regex_match(line, std::regex("^\\s*edges\\s+of\\s+the\\s+boundary\\s*$")))
            return fvca5_parse_boundary_edges();

        /* JMLC format extension */
        if (std::regex_match(line, std::regex("^\\s*edges\\s+of\\s+transmission\\s*$")))
            return fvca5_parse_boundary_edges(); /* Yes, it is fvca5_parse_boundary_edges() */

        if (std::regex_match(line, std::regex("^\\s*all\\s+edges\\s*$")))
            return fvca5_parse_all_edges();

        std::cout << "Line " << current_line << ": Unknown block '";
        std::cout << line << "'" << std::endl;
        return false;
    }

    bool
    parse(const std::string& filename)
    {
        ifs.open(filename);
        if (not ifs.is_open())
        {
            std::cout << "Can't open " << filename << std::endl;
            return false;
        }

        while (!ifs.eof())
        {
            bool success = fvca5_parse_block();

            if (not success)
            {
                m_points.clear();
                m_polys.clear();
                m_boundary_edges.clear();
                m_edges.clear();
                return false;
            }
        }

        return true;
    }

  public:
    static const char constexpr* expected_extension = "typ1";

    fvca5_mesh_loader() : current_line(0) {}

    bool
    read_mesh(const std::string& filename)
    {
        if (this->verbose())
            std::cout << " *** READING FVCA5 MESH ***" << std::endl;
        return parse(filename);
    }

    bool
    populate_mesh(mesh_type& msh)
    {
        if (this->verbose())
            std::cout << " *** POPULATING FVCA5 MESH ***" << std::endl;
        auto storage = msh.backend_storage();

        /* Points */
        size_t nodes_size = m_points.size();
        storage->points   = std::move(m_points);

        /* Nodes */
        std::vector<node_type> nodes(nodes_size);
        for (size_t i = 0; i < nodes_size; i++)
            nodes[i] = node_type(disk::point_identifier<2>(i));

        storage->nodes = std::move(nodes);

        /* Edges */
        /* Make the vector containing the edges */
        std::vector<edge_type> edges;
        edges.reserve(m_edges.size());
        for (size_t i = 0; i < m_edges.size(); i++)
        {
            assert(m_edges[i][0] < m_edges[i][1]);
            auto node1 = typename node_type::id_type(m_edges[i][0]);
            auto node2 = typename node_type::id_type(m_edges[i][1]);

            auto e = edge_type(node1, node2);

            /* Next line not necessary anymore, see generic_element<DIM, DIM-1> */
            // e.set_point_ids(m_edges[i].begin(), m_edges[i].begin()+2); /* XXX: crap */
            edges.push_back(e);
        }
        /* Sort them */
        std::sort(edges.begin(), edges.end());

        /* Detect which ones are boundary edges */
        storage->boundary_info.resize(m_edges.size());
        for (size_t i = 0; i < m_boundary_edges.size(); i++)
        {
            assert(m_boundary_edges[i][0] < m_boundary_edges[i][1]);
            auto node1 = typename node_type::id_type(m_boundary_edges[i][0]);
            auto node2 = typename node_type::id_type(m_boundary_edges[i][1]);

            auto e = edge_type(node1, node2);

            auto position = find_element_id(edges.begin(), edges.end(), e);

            if (position.first == false)
            {
                std::cout << "Bad bug at " << __FILE__ << "(" << __LINE__ << ")" << std::endl;
                return false;
            }

            boundary_descriptor bi(0, true);
            storage->boundary_info.at(position.second) = bi;
        }

        storage->edges = std::move(edges);

        /* Surfaces */
        using sd_pair = std::pair<surface_type, size_t>;
        std::vector<sd_pair> surfaces;
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
            surfaces.push_back(std::make_pair(surface, p.domain_id));
            // storage->subdomain_info.push_back( subdomain_descriptor(p.domain_id) );
        }

        auto comp = [](const sd_pair& a, const sd_pair& b) -> bool { return a.first < b.first; };

        std::sort(surfaces.begin(), surfaces.end(), comp);

        auto get_surf = [](const sd_pair& sdp) -> auto { return sdp.first; };
        storage->surfaces.resize(surfaces.size());
        std::transform(surfaces.begin(), surfaces.end(), storage->surfaces.begin(), get_surf);

        auto get_id = [](const sd_pair& sdp) -> auto { return subdomain_descriptor(sdp.second); };
        storage->subdomain_info.resize(surfaces.size());
        std::transform(surfaces.begin(), surfaces.end(), storage->subdomain_info.begin(), get_id);

        mark_internal_faces(msh);

        /* FVCA5 meshes do not provide boundary information. As usually they represent
         * the unit square, we try to renumber according to the normals of the boundary
         * surfaces. If it fails, everything is marked as boundary 0.
         */

        std::vector<boundary_descriptor> bds_save   = storage->boundary_info;
        auto                             renumbered = renumber_hypercube_boundaries(msh);
        if (not renumbered)
        {
            std::cout << "Warning: unable to renumber FVCA5 mesh boundaries, ";
            std::cout << "defaulting everyting to 0.\nProbably one of the boundaries ";
            std::cout << "is not parallel to the x or y axes.";
            std::cout << std::endl;
            std::swap(bds_save, storage->boundary_info);
        }

        if (this->verbose())
        {
            std::cout << "Nodes:    " << storage->nodes.size() << std::endl;
            std::cout << "Edges:    " << storage->edges.size() << std::endl;
            std::cout << "Faces:    " << storage->surfaces.size() << std::endl;
            std::cout << "BND info: " << storage->boundary_info.size() << std::endl;
            std::cout << "SD info:  " << storage->subdomain_info.size() << std::endl;
        }

        return true;
    }
};

}