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

#include <tuple>
#include <vector>
#include <thread>
#include <fstream>

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/geometry/geometry.hpp"
#include "diskpp/loaders/mesh_loader.hpp"
#include "diskpp/loaders/strtot.hpp"
#include "diskpp/common/mapped_file.h"

namespace disk
{

template<typename T, size_t N>
class netgen_mesh_loader
{
    static_assert(N == 2 || N == 3, "Netgen supports only 2D or 3D");
};

namespace priv
{

template<typename T>
std::tuple<T, T, T, T, T>
read_tetrahedron_line(const char* str, char** endptr)
{
    T t1, t2, t3, t4, t5;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    t4 = strtot<T>(*endptr, endptr);
    t5 = strtot<T>(*endptr, endptr);

    return std::make_tuple(t1, t2 - 1, t3 - 1, t4 - 1, t5 - 1);
}

template<typename T>
std::tuple<T, T, T, T>
read_triangle_line(const char* str, char** endptr)
{
    T t1, t2, t3, t4;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    t4 = strtot<T>(*endptr, endptr);

    return std::make_tuple(t1, t2 - 1, t3 - 1, t4 - 1);
}

template<typename T>
std::tuple<T, T, T>
read_edge_line(const char* str, char** endptr)
{
    T t1, t2, t3;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);

    return std::make_tuple(t1, t2 - 1, t3 - 1);
}

} // namespace priv

template<typename T>
class netgen_mesh_loader<T, 2> : public mesh_loader<simplicial_mesh<T, 2>>
{
    typedef simplicial_mesh<T, 2>            mesh_type;
    typedef typename mesh_type::point_type   point_type;
    typedef typename mesh_type::node_type    node_type;
    typedef typename mesh_type::edge_type    edge_type;
    typedef typename mesh_type::surface_type surface_type;

    std::vector<point_type>                   points;
    std::vector<node_type>                    nodes;
    std::vector<edge_type>                    edges;
    std::vector<std::pair<edge_type, size_t>> boundary_edges;
    std::vector<surface_type>                 surfaces;

    bool
    netgen_read(const std::string& filename)
    {
        /* Open file */
        if (filename.size() == 0)
        {
            std::cout << "Invalid mesh file name" << std::endl;
            return false;
        }

        size_t lines, linecount;

        mapped_file mf(filename);

        // std::cout << green << " * * * Reading NETGEN format mesh * * * ";
        // std::cout << nocolor << std::endl;

        /************************ Read points ************************/
        linecount = 0;

        const char* data = mf.mem();
        char*       endptr;

        lines = strtot<size_t>(data, &endptr);

        points.reserve(lines);
        nodes.reserve(lines);

        while (linecount < lines)
        {
            if (this->verbose() && ((linecount % 100000) == 0))
            {
                std::cout << "Reading points: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto point = priv::read_2d_point_line<T>(endptr, &endptr, 1.0);

            points.push_back(point);

            auto point_id = disk::point_identifier<2>(linecount);
            auto node     = node_type({point_id});

            nodes.push_back(node);

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading points: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        /************************ Read triangles ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        edges.reserve(lines * 3);
        surfaces.reserve(lines);

        while (linecount < lines)
        {
            if (this->verbose() && ((linecount % 100000) == 0))
            {
                std::cout << "Reading triangles: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_triangle_line<size_t>(endptr, &endptr);

            disk::point_identifier<2> p0(std::get<1>(t));
            disk::point_identifier<2> p1(std::get<2>(t));
            disk::point_identifier<2> p2(std::get<3>(t));
            // domain_id_type      d(std::get<0>(t));

            edges.push_back(edge_type({p0, p1}));
            edges.push_back(edge_type({p1, p2}));
            edges.push_back(edge_type({p0, p2}));

            surfaces.push_back(surface_type({p0, p1, p2}));

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading triangles: " << linecount;
            std::cout << "/" << lines << std::endl;
        }
        /************************ Read boundary edges ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        boundary_edges.reserve(lines);

        while (linecount < lines)
        {
            if (this->verbose() && ((linecount % 50000) == 0))
            {
                std::cout << "Reading edges: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_edge_line<size_t>(endptr, &endptr);

            disk::point_identifier<2> p0(std::get<1>(t));
            disk::point_identifier<2> p1(std::get<2>(t));

            edge_type edge({p0, p1});

            boundary_edges.push_back(std::make_pair(edge, std::get<0>(t)));

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading edges: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        return true;
    }

  public:
    static const char constexpr* expected_extension = "mesh2d";
    netgen_mesh_loader()                            = default;

    bool
    read_mesh(const std::string& s)
    {
        if (this->verbose())
            std::cout << " *** READING NETGEN 2D MESH ***" << std::endl;

        return netgen_read(s);
    }

    bool
    populate_mesh(mesh_type& msh)
    {
        auto storage = msh.backend_storage();

        if (this->verbose())
        {
            std::cout << "Sorting data...";
            std::cout.flush();
        }

        storage->points = std::move(points);
        storage->nodes  = std::move(nodes);

        /* sort edges, make unique and move them in geometry */
        THREAD(edge_thread, priv::sort_uniq(edges); storage->edges = std::move(edges););

        /* sort triangles, make unique and move them in geometry */
        THREAD(tri_thread, priv::sort_uniq(surfaces); storage->surfaces = std::move(surfaces););

        /* wait for the threads */
        WAIT_THREAD(edge_thread);
        WAIT_THREAD(tri_thread);

        storage->boundary_info.resize(storage->edges.size());
        for (auto& be : boundary_edges)
        {
            auto position = find_element_id(storage->edges.begin(), storage->edges.end(), be.first);
            if (position.first == false)
            {
                std::cout << "Bad bug at " << __FILE__ << "(" << __LINE__ << ")" << std::endl;
                return false;
            }
            boundary_descriptor bi(be.second, true);
            storage->boundary_info.at(position.second) = bi;
        }

        if (this->verbose())
        {
            std::cout << "done." << std::endl;

            std::cout << "Nodes: " << storage->nodes.size() << std::endl;
            std::cout << "Edges: " << storage->edges.size() << std::endl;
            std::cout << "Faces: " << storage->surfaces.size() << std::endl;
        }

        boundary_edges.clear();

        storage->subdomain_info.resize(storage->surfaces.size());

        return true;
    }
};

template<typename T>
class netgen_mesh_loader<T, 3> : public mesh_loader<simplicial_mesh<T, 3>>
{
    typedef simplicial_mesh<T, 3>            mesh_type;
    typedef typename mesh_type::point_type   point_type;
    typedef typename mesh_type::node_type    node_type;
    typedef typename mesh_type::edge_type    edge_type;
    typedef typename mesh_type::surface_type surface_type;
    typedef typename mesh_type::volume_type  volume_type;

    std::vector<point_type>                                   points;
    std::vector<node_type>                                    nodes;
    std::vector<edge_type>                                    edges;
    std::vector<surface_type>                                 surfaces;
    std::vector<boundary_descriptor>                          boundary_info;
    std::vector<std::pair<volume_type, subdomain_descriptor>> tmp_volumes;

    bool
    netgen_read(const std::string& filename)
    {
        /* Open file */
        if (filename.size() == 0)
        {
            std::cout << "Invalid mesh file name" << std::endl;
            return false;
        }

        size_t lines, linecount;

        mapped_file mf(filename);

        /************************ Read points ************************/
        linecount = 0;

        const char* data = mf.mem();
        char*       endptr;

        lines = strtot<size_t>(data, &endptr);

        points.reserve(lines);
        nodes.reserve(lines);

        while (linecount < lines)
        {
            if (this->verbose() && ((linecount % 100000) == 0))
            {
                std::cout << "Reading points: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto point = priv::read_3d_point_line<T>(endptr, &endptr, 1.0);

            points.push_back(point);

            auto point_id = disk::point_identifier<3>(linecount);
            auto node     = node_type({point_id});

            nodes.push_back(node);
            /* Do something with that point */

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading points: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        /************************ Read tetrahedra ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        edges.reserve(lines * 6);
        surfaces.reserve(lines * 4);
        tmp_volumes.reserve(lines);

        while (linecount < lines)
        {
            if (this->verbose() && ((linecount % 100000) == 0))
            {
                std::cout << "Reading tetrahedra: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_tetrahedron_line<size_t>(endptr, &endptr);

            auto                      subdomain_num = std::get<0>(t);
            disk::point_identifier<3> p0(std::get<1>(t));
            disk::point_identifier<3> p1(std::get<2>(t));
            disk::point_identifier<3> p2(std::get<3>(t));
            disk::point_identifier<3> p3(std::get<4>(t));

            edges.push_back(edge_type({p0, p1}));
            edges.push_back(edge_type({p0, p2}));
            edges.push_back(edge_type({p0, p3}));
            edges.push_back(edge_type({p1, p2}));
            edges.push_back(edge_type({p1, p3}));
            edges.push_back(edge_type({p2, p3}));

            surfaces.push_back(surface_type({p0, p1, p2}));
            surfaces.push_back(surface_type({p0, p1, p3}));
            surfaces.push_back(surface_type({p0, p2, p3}));
            surfaces.push_back(surface_type({p1, p2, p3}));

            subdomain_descriptor si(subdomain_num);
            auto                 vs = std::make_pair(volume_type({p0, p1, p2, p3}), si);
            tmp_volumes.push_back(vs);

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading tetrahedra: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        priv::sort_uniq(edges);
        priv::sort_uniq(surfaces);
        boundary_info.resize(surfaces.size());

        using pvs        = std::pair<volume_type, subdomain_descriptor>;
        auto tmpvol_comp = [](const pvs& a, const pvs& b) { return a.first < b.first; };

        std::sort(tmp_volumes.begin(), tmp_volumes.end(), tmpvol_comp);

        /************************ Read boundary surfaces ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);
        while (linecount < lines)
        {
            if (this->verbose() && ((linecount % 50000) == 0))
            {
                std::cout << "Reading triangle: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_triangle_line<size_t>(endptr, &endptr);

            auto                      bnd_id = std::get<0>(t);
            disk::point_identifier<3> p0(std::get<1>(t));
            disk::point_identifier<3> p1(std::get<2>(t));
            disk::point_identifier<3> p2(std::get<3>(t));

            boundary_descriptor bi(bnd_id, true);
            surface_type        surf({p0, p1, p2});

            auto itor = std::lower_bound(surfaces.begin(), surfaces.end(), surf);
            if ((itor == surfaces.end()) or not(*itor == surf))
                throw std::logic_error("Face not found");

            auto ofs              = std::distance(surfaces.begin(), itor);
            boundary_info.at(ofs) = bi;

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading triangle: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        return true;
    }

  public:
    static const char constexpr* expected_extension = "mesh";
    netgen_mesh_loader()                            = default;

    bool
    read_mesh(const std::string& s)
    {
        if (this->verbose())
            std::cout << " *** READING NETGEN 3D MESH ***" << std::endl;
        return netgen_read(s);
    }

    bool
    populate_mesh(mesh_type& msh)
    {
        auto storage = msh.backend_storage();

        if (this->verbose())
        {
            std::cout << "Sorting data...";
            std::cout.flush();
        }

        storage->points        = std::move(points);
        storage->nodes         = std::move(nodes);
        storage->edges         = std::move(edges);
        storage->surfaces      = std::move(surfaces);
        storage->boundary_info = std::move(boundary_info);

        std::vector<volume_type> vols;
        vols.reserve(tmp_volumes.size());

        std::vector<subdomain_descriptor> subdoms;
        subdoms.reserve(tmp_volumes.size());

        for (auto& [vol, sd] : tmp_volumes)
        {
            vols.push_back(vol);
            subdoms.push_back(sd);
        }

        tmp_volumes.clear();

        storage->volumes        = std::move(vols);
        storage->subdomain_info = std::move(subdoms);

        mark_internal_faces(msh);

        if (this->verbose())
        {
            std::cout << "done." << std::endl;

            std::cout << "Nodes: " << storage->nodes.size() << std::endl;
            std::cout << "Edges: " << storage->edges.size() << std::endl;
            std::cout << "Faces: " << storage->surfaces.size() << std::endl;
            std::cout << "Volumes: " << storage->volumes.size() << std::endl;
        }

        return true;
    }
};

} // end namespace disk
