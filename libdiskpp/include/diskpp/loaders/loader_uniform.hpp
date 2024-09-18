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

#include "diskpp/geometry/geometry.hpp"
#include "diskpp/loaders/mesh_loader.hpp"
#include "diskpp/mesh/mesh.hpp"

namespace disk
{

template<typename T, size_t N>
class uniform_mesh_loader
{
    static_assert(N == 1, "at the moment only 1D uniform meshes are available");
};

template<typename T>
class uniform_mesh_loader<T, 1> : public mesh_loader<generic_mesh<T, 1>>
{
    typedef generic_mesh<T, 1>             mesh_type;
    typedef typename mesh_type::point_type point_type;
    typedef typename mesh_type::node_type  node_type;
    typedef typename mesh_type::edge_type  edge_type;

    T      m_h, m_x_min, m_x_max;
    size_t m_number_of_elements;

  public:
    uniform_mesh_loader() : m_x_min(T(0)), m_x_max(T(1)), m_number_of_elements(8)
    {
        m_h = fabs(m_x_max - m_x_min) / m_number_of_elements;
    }

    uniform_mesh_loader(T x_min, T x_max, size_t N) : m_x_min(x_min), m_x_max(x_max), m_number_of_elements(N)
    {
        m_h = fabs(m_x_max - m_x_min) / m_number_of_elements;
    }

    bool
    populate_mesh(mesh_type& msh)
    {
        if (this->verbose())
            std::cout << " *** POPULATING UNIFORM 1D MESH ***" << std::endl;

        auto storage = msh.backend_storage();

        auto num_edges = m_number_of_elements;
        auto num_nodes = m_number_of_elements + 1;

        storage->points.resize(num_nodes);
        storage->nodes.resize(num_nodes);
        storage->edges.resize(num_edges);

        for (size_t i = 0; i < num_edges; i++)
        {
            storage->points.at(i)     = point_type({m_x_min + (m_h * i)});
            storage->points.at(i + 1) = point_type({m_x_min + (m_h * (i + 1))});

            auto n0 = typename node_type::id_type(i);
            auto n1 = typename node_type::id_type(i + 1);
            auto e  = edge_type(n0, n1);

            storage->edges.at(i) = e;
        }

        for (size_t i = 0; i < num_nodes; i++)
            storage->nodes.at(i) = node_type(point_identifier<1>(i));

        storage->boundary_info.resize(num_nodes);
        boundary_descriptor bi(0, true);
        storage->boundary_info.at(0)             = bi;
        storage->boundary_info.at(num_nodes - 1) = bi;

        return true;
    }
};

}