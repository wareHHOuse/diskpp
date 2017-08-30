/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#pragma once

#include <vector>

#include "ident.hpp"
#include "point.hpp"

namespace disk {

/* Generic template for mesh_storage.
 *
 * This template has to be specialized for the 1D, 2D and 3D cases.
 * The function of mesh_storage is to decouple the low-level details of the
 * storage of the mesh elements (points, nodes, edges...) from the actual view
 * on the mesh. mesh_storage is part of the internal representation of the
 * mesh, for this reason is in the `priv` namespace. Internal representation
 * could change, breaking the user code, so all accesses have to be proxied by
 * `mesh_base`.
 * However, in some situations, accessing the internal representation can be
 * useful if not essential. This is the case for example of the mesh loaders.
 *
 * @DIM     Space dimension
 * @T       Type that the mesh uses to represent points.
 */
template<typename T, size_t DIM, typename StorageClass>
class mesh_storage
{
    static_assert(DIM > 0 && DIM <= 3, "mesh: Allowed dimensions are 1, 2 and 3");
};

struct bnd_info
{
    size_t  boundary_id;
    bool    is_boundary;
};

/* Template specialization of mesh_storage for 3D meshes.
 *
 * @T       Type that the mesh uses to represent points.
 */
template<typename T, typename StorageClass>
struct mesh_storage<T, 3, StorageClass>
{
    typedef StorageClass                            sc_type;
    typedef typename sc_type::volume_type           volume_type;
    typedef typename sc_type::surface_type          surface_type;
    typedef typename sc_type::edge_type             edge_type;
    typedef typename sc_type::node_type             node_type;
 
    typedef T                                       coordinate_type;
    static const size_t                             dimension = 3;
    typedef point<coordinate_type,3>                point_type;

    template<typename ET> using id = identifier<ET, ident_impl_t, 0>;
    typedef id<volume_type>                         volume_id_type;
    typedef id<surface_type>                        surface_id_type;
    typedef id<edge_type>                           edge_id_type;
    typedef id<node_type>                           node_id_type;
    //typedef id<point_type>                          point_id_type;

    std::vector<volume_type>                        volumes;
    std::vector<surface_type>                       surfaces;
    std::vector<edge_type>                          edges;
    std::vector<node_type>                          nodes;
    std::vector<point_type>                         points;

    std::vector<bnd_info>                           boundary_info;

    void statistics(void) const
    {
        std::cout << "This is a storage for a 3D mesh" << std::endl;
        std::cout << "Points: " << points.size() << std::endl;
        std::cout << "Nodes: " << nodes.size() << std::endl;
        std::cout << "Edges: " << edges.size() << std::endl;
        std::cout << "Surfaces: " << surfaces.size() << std::endl;
        std::cout << "Volumes: " << volumes.size() << std::endl;
        auto bs = std::count_if(boundary_info.begin(), boundary_info.end(), 
                    [&](const bnd_info& bi){ return bi.is_boundary; }  );
        std::cout << "Boundary surfaces: " << bs << std::endl;
    }
};

/* Template specialization of mesh_storage for 2D meshes.
 *
 * @T       Type that the mesh uses to represent points.
 */
template<typename T, typename StorageClass>
struct mesh_storage<T, 2, StorageClass>
{
    typedef StorageClass                            sc_type;
    typedef typename sc_type::surface_type          surface_type;
    typedef typename sc_type::edge_type             edge_type;
    typedef typename sc_type::node_type             node_type;

    typedef T                                       coordinate_type;
    static const size_t                             dimension = 2;
    typedef point<coordinate_type,2>                point_type;

    template<typename ET> using id = identifier<ET, ident_impl_t, 0>;
    typedef id<surface_type>                        surface_id_type;
    typedef id<edge_type>                           edge_id_type;
    typedef id<node_type>                           node_id_type;
    //typedef id<point_type>                          point_id_type;

    std::vector<surface_type>                       surfaces;
    std::vector<edge_type>                          edges;
    std::vector<node_type>                          nodes;
    std::vector<point_type>                         points;

    std::vector<bnd_info>                           boundary_info;

    void statistics(void) const
    {
        std::cout << "This is a storage for a 2D mesh" << std::endl;
        std::cout << "Points: " << points.size() << std::endl;
        std::cout << "Nodes: " << nodes.size() << std::endl;
        std::cout << "Edges: " << edges.size() << std::endl;
        std::cout << "Surfaces: " << surfaces.size() << std::endl;
        //auto be = std::count(is_boundary.begin(), is_boundary.end(), true);
        //std::cout << "Boundary edges: " << be << std::endl;
    }
};

/* Template specialization of mesh_storage for 1D meshes.
 *
 * @T       Type that the mesh uses to represent points.
 */
template<typename T, typename StorageClass>
struct mesh_storage<T, 1, StorageClass>
{
    typedef StorageClass                            sc_type;
    typedef typename sc_type::edge_type             edge_type;
    typedef typename sc_type::node_type             node_type;

    typedef T                                       coordinate_type;
    static const size_t                             dimension = 1;
    typedef point<coordinate_type,1>                point_type;

    template<typename ET> using id = identifier<ET, ident_impl_t, 0>;
    typedef id<edge_type>                           edge_id_type;
    typedef id<node_type>                           node_id_type;
    //typedef id<point_type>                          point_id_type;

    std::vector<edge_type>                          edges;
    std::vector<node_type>                          nodes;
    std::vector<point_type>                         points;

    std::vector<bnd_info>                           boundary_info;

    void statistics(void) const
    {
        std::cout << "This is a storage for a 1D mesh" << std::endl;
        std::cout << "Points: " << points.size() << std::endl;
        std::cout << "Nodes: " << nodes.size() << std::endl;
        std::cout << "Edges: " << edges.size() << std::endl;
    }
};


template<typename>
struct mesh_storage_traits {};

template<typename T, typename StorageClass>
struct mesh_storage_traits<mesh_storage<T, 3, StorageClass>>
{
    typedef mesh_storage<T, 3, StorageClass>            storage_type;
    typedef typename storage_type::node_type            node_type;
    typedef typename storage_type::node_id_type         node_id_type;
    typedef typename storage_type::edge_type            edge_type;
    typedef typename storage_type::edge_id_type         edge_id_type;
    typedef typename storage_type::surface_type         surface_type;
    typedef typename storage_type::surface_id_type      surface_id_type;
    typedef typename storage_type::volume_type          volume_type;
    typedef typename storage_type::volume_id_type       volume_id_type;
    typedef typename storage_type::coordinate_type      coordinate_type;
    typedef typename storage_type::point_type           point_type;
    static const size_t                                 dimension = 3;
};

template<typename T, typename StorageClass>
struct mesh_storage_traits<mesh_storage<T, 2, StorageClass>>
{
    typedef mesh_storage<T, 2, StorageClass>            storage_type;
    typedef typename storage_type::node_type            node_type;
    typedef typename storage_type::node_id_type         node_id_type;
    typedef typename storage_type::edge_type            edge_type;
    typedef typename storage_type::edge_id_type         edge_id_type;
    typedef typename storage_type::surface_type         surface_type;
    typedef typename storage_type::surface_id_type      surface_id_type;
    typedef typename storage_type::coordinate_type      coordinate_type;
    typedef typename storage_type::point_type           point_type;
    static const size_t                                 dimension = 2;
};

template<typename T, typename StorageClass>
struct mesh_storage_traits<mesh_storage<T, 1, StorageClass>>
{
    typedef mesh_storage<T, 1, StorageClass>            storage_type;
    typedef typename storage_type::node_type            node_type;
    typedef typename storage_type::node_id_type         node_id_type;
    typedef typename storage_type::edge_type            edge_type;
    typedef typename storage_type::edge_id_type         edge_id_type;
    typedef typename storage_type::coordinate_type      coordinate_type;
    typedef typename storage_type::point_type           point_type;
    static const size_t                                 dimension = 1;
};

} // namespace disk
