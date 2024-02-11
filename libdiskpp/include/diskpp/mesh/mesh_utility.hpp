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

#include <vector>

#include "mesh/mesh.hpp"
#include "diskpp/common/eigen.hpp"
#include "solvers/solver.hpp"


namespace disk {

/* The submesher takes a 2D mesh of any kind and computes a triangulation for its
 * elements. The triangulation is obtained by splitting the original element, so
 * the quality of the result could happen to be low. */

template<typename Mesh>
class submesher;

template<template<typename, size_t, typename> class Mesh,
            typename T, typename Storage>
class submesher<Mesh<T,2,Storage>>
{
public:
    typedef Mesh<T,2,Storage>                       outer_mesh_type;
    typedef typename outer_mesh_type::cell          outer_cell_type;

    typedef simplicial_mesh<T,2>                    inner_mesh_type;
    typedef typename inner_mesh_type::cell          inner_cell_type;

    typedef typename outer_mesh_type::point_type    point_type;

private:
    struct edge
    {
        edge(){}

        edge(size_t ap0, size_t ap1, bool bnd = false)
        {
            p0                  = ap0;
            p1                  = ap1;
            is_boundary         = bnd;
            is_broken           = false;
        }

        edge(size_t ap0, size_t ap1, size_t bid, bool bnd)
        {
            p0                  = ap0;
            p1                  = ap1;
            is_boundary         = bnd;
            is_broken           = false;
            boundary_id         = bid;
        }

        friend bool operator<(const edge& a, const edge& b) {
            return ( a.p0 < b.p0 or (a.p0 == b.p0 and a.p1 < b.p1) );
        }

        friend bool operator==(const edge& a, const edge& b) {
            return ( a.p0 == b.p0 and a.p1 == b.p1 );
        }

        size_t  p0, p1, pb;
        bool    is_boundary, is_broken;
        size_t  boundary_id;
    };


    struct triangle
    {
        triangle() {}
        triangle(size_t ap0, size_t ap1, size_t ap2) : p0(ap0), p1(ap1), p2(ap2) {}

        size_t p0, p1, p2;
    };

    std::vector<point_type>     m_points;
    std::vector<edge>           m_edges;
    std::list<triangle>         m_triangles;

    void refine(size_t steps)
    {
        if (steps == 0)
            return;

        /* break all the edges */
        for (auto& e : m_edges)
        {
            assert(!e.is_broken);
            auto bar = (m_points[e.p0] + m_points[e.p1])*0.5;
            e.pb = m_points.size();
            m_points.push_back(bar);
            e.is_broken = true;
        }

        /* refine the triangles */
        size_t triangles_to_process = m_triangles.size();
        m_edges.reserve( m_edges.size() + triangles_to_process*12 );
        auto edges_end = m_edges.end();
        for (size_t i = 0; i < triangles_to_process; i++)
        {
            /* remove the original triangle */
            triangle t = m_triangles.front();
            m_triangles.pop_front();

            /* find the edges of the triangle */
            auto t_e0 = *std::lower_bound(m_edges.begin(), edges_end, edge(t.p0, t.p1));
            assert(t_e0.is_broken);
            auto t_e1 = *std::lower_bound(m_edges.begin(), edges_end, edge(t.p1, t.p2));
            assert(t_e1.is_broken);
            auto t_e2 = *std::lower_bound(m_edges.begin(), edges_end, edge(t.p0, t.p2));
            assert(t_e2.is_broken);

            assert(t_e0.p1 == t_e1.p0);
            assert(t_e1.p1 == t_e2.p1);
            assert(t_e0.p0 == t_e2.p0);

            /* compute the edges of the new four triangles */
            /* first triangle */
            m_edges.push_back( edge(t_e0.p0, t_e0.pb, t_e0.boundary_id, t_e0.is_boundary) );
            m_edges.push_back( edge(t_e0.pb, t_e2.pb, false) );
            m_edges.push_back( edge(t_e0.p0, t_e2.pb, t_e2.boundary_id, t_e2.is_boundary) );
            triangle t0(t_e0.p0, t_e0.pb, t_e2.pb);
            m_triangles.push_back(t0);

            /* second triangle */
            m_edges.push_back( edge(t_e0.p1, t_e0.pb, t_e0.boundary_id, t_e0.is_boundary) );
            m_edges.push_back( edge(t_e0.pb, t_e1.pb, false) );
            m_edges.push_back( edge(t_e1.p0, t_e1.pb, t_e1.boundary_id, t_e1.is_boundary) );
            triangle t1(t_e0.p1, t_e0.pb, t_e1.pb);
            m_triangles.push_back(t1);

            /* third triangle */
            m_edges.push_back( edge(t_e1.p1, t_e1.pb, t_e1.boundary_id, t_e1.is_boundary) );
            m_edges.push_back( edge(t_e1.pb, t_e2.pb, false) );
            m_edges.push_back( edge(t_e2.p1, t_e2.pb, t_e2.boundary_id, t_e2.is_boundary) );
            triangle t2(t_e1.p1, t_e1.pb, t_e2.pb);
            m_triangles.push_back(t2);

            /* fourth triangle */
            triangle t3(t_e0.pb, t_e1.pb, t_e2.pb);
            m_triangles.push_back(t3);
        }

        /* we don't need the broken edges anymore, discard them */
        auto edge_is_broken = [](const edge& e) -> bool { return e.is_broken; };
        std::remove_if(m_edges.begin(), m_edges.end(), edge_is_broken);

        /* sort the edges to allow fast lookups */
        priv::sort_uniq(m_edges);

        refine(steps - 1);
    }

public:

    submesher() {};

    inner_mesh_type
    generate_mesh(const outer_mesh_type& msh, const outer_cell_type& cl,
                  size_t refinement_steps = 3)
    {
        m_points.clear();
        m_edges.clear();
        m_triangles.clear();

        /* Transform the triangle in internal representation */
        auto pts = points(msh, cl);
        m_points.insert(m_points.begin(), pts.begin(), pts.end());

        /* Here the mapping between the local boundary number and the parent
         * cell boundary number is established. */
        auto fcs = faces(msh, cl);
        assert(fcs.size() == pts.size());

        auto num_faces = fcs.size();
        /* In this case we have just a triangle. Put it as is in the new mesh. */
        if ( num_faces == 3 )
        {
            edge e0(0, 1, msh.lookup(fcs[0]), true);
            m_edges.push_back(e0);

            edge e1(1, 2, msh.lookup(fcs[1]), true);
            m_edges.push_back(e1);

            edge e2(0, 2, msh.lookup(fcs[2]), true);
            m_edges.push_back(e2);

            std::sort(m_edges.begin(), m_edges.end());

            triangle t(0,1,2);
            m_triangles.push_back(t);
        }
        else /* We split the element in triangles */
        {
            auto bar = barycenter(msh, cl);
            m_points.push_back(bar);
            for (size_t i = 0; i < num_faces; i++)
            {
                /* This is the edge corresponding to the boundary of
                 * our initial element */
                edge e0(i, (i+1)%num_faces, msh.lookup(fcs[i]), true);
                m_edges.push_back(e0);

                /* The following two are the new internal edges */
                edge e1(i, num_faces);
                m_edges.push_back(e1);

                edge e2((i+1)%num_faces, num_faces);
                m_edges.push_back(e2);

                /* Define the new triangle, */
                triangle t(i, (i+1)%num_faces, num_faces);
                m_triangles.push_back(t);
            }

            std::sort(m_edges.begin(), m_edges.end());
        }

        refine( refinement_steps );

        inner_mesh_type refined_element_submesh;
        auto storage = refined_element_submesh.backend_storage();

        storage->points = std::move(m_points);
        storage->nodes.reserve( storage->points.size() );
        for (size_t i = 0; i < storage->points.size(); i++)
        {
            auto point_id = point_identifier<2>( i );
            storage->nodes.push_back( typename inner_mesh_type::node_type({point_id}));
        }

        std::vector<std::pair<typename inner_mesh_type::edge_type, size_t>> be;
        be.reserve( 3 * m_triangles.size() );

        storage->edges.reserve( 3 * m_triangles.size() );
        storage->surfaces.reserve( m_triangles.size() );
        for (auto& t : m_triangles)
        {
            auto pi0 = point_identifier<2>(t.p0);
            auto pi1 = point_identifier<2>(t.p1);
            auto pi2 = point_identifier<2>(t.p2);

            auto e0 = typename inner_mesh_type::edge_type({pi0, pi1});
            auto e1 = typename inner_mesh_type::edge_type({pi1, pi2});
            auto e2 = typename inner_mesh_type::edge_type({pi0, pi2});

            storage->edges.push_back( e0 );
            storage->edges.push_back( e1 );
            storage->edges.push_back( e2 );
            storage->surfaces.push_back( typename inner_mesh_type::surface_type({pi0, pi1, pi2}) );

            auto t_e0 = *std::lower_bound(m_edges.begin(), m_edges.end(), edge(t.p0, t.p1));
            assert(!t_e0.is_broken);
            auto t_e1 = *std::lower_bound(m_edges.begin(), m_edges.end(), edge(t.p1, t.p2));
            assert(!t_e1.is_broken);
            auto t_e2 = *std::lower_bound(m_edges.begin(), m_edges.end(), edge(t.p0, t.p2));
            assert(!t_e2.is_broken);

            if (t_e0.is_boundary) be.push_back( std::make_pair(e0, t_e0.boundary_id) );
            if (t_e1.is_boundary) be.push_back( std::make_pair(e1, t_e1.boundary_id) );
            if (t_e2.is_boundary) be.push_back( std::make_pair(e2, t_e2.boundary_id) );
        }

        THREAD(edge_thread,
                priv::sort_uniq(storage->edges);
        );

        THREAD(tri_thread,
                priv::sort_uniq(storage->surfaces);
        );

        WAIT_THREAD(edge_thread);
        WAIT_THREAD(tri_thread);

        storage->boundary_info.resize( storage->edges.size() );
        for (auto& e : be)
        {
            auto ei = find_element_id(storage->edges.begin(),
                                      storage->edges.end(), e.first);

            if (ei.first == false)
                throw std::logic_error("Edge not found. This is a bug.");

            bnd_info bi;
            bi.boundary_id = e.second;
            bi.is_boundary = true;
            storage->boundary_info.at(ei.second) = bi;
        }

        return refined_element_submesh;
    }
};


/* Minimize the energy contained in a triangular mesh. */

template<typename T>
void relax_mesh(disk::simplicial_mesh<T,2>& msh)
{
    typedef Eigen::SparseMatrix<T>  sparse_matrix_type;
    typedef Eigen::Triplet<T>       triplet_type;

    typedef typename disk::simplicial_mesh<T,2>::point_type point_type;

    std::vector<bool> dirichlet_nodes( msh.points_size() );
    for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
    {
        auto fc = *itor;
        auto ptids = fc.point_ids();
        dirichlet_nodes.at(ptids[0]) = true;
        dirichlet_nodes.at(ptids[1]) = true;
    }

    std::vector<size_t> compress_map, expand_map;
    compress_map.resize( msh.points_size() );
    size_t system_size = std::count_if(dirichlet_nodes.begin(), dirichlet_nodes.end(), [](bool d) -> bool {return !d;});
    expand_map.resize( system_size );

    auto nnum = 0;
    for (size_t i = 0; i < msh.points_size(); i++)
    {
        if ( dirichlet_nodes.at(i) )
            continue;

        expand_map.at(nnum) = i;
        compress_map.at(i) = nnum++;
    }

    sparse_matrix_type          gA(2*system_size, 2*system_size);
    dynamic_vector<T>           gb = dynamic_vector<T>::Zero(2*system_size);
    dynamic_vector<T>           gx = dynamic_vector<T>::Zero(2*system_size);
    std::vector<triplet_type>   gtA;

    Eigen::Matrix<T,3,3> G = Eigen::Matrix<T,3,3>::Zero();
    G(0,0) = -1; G(0,1) = 1;
    G(1,1) = -1; G(1,2) = 1;
    G(2,0) = -1; G(2,2) = 1;

    auto backend_storage = msh.backend_storage();

    for (size_t i = 0; i < system_size; i++)
    {
        auto pt = backend_storage->points.at( expand_map.at(i) );
        gx(i) = pt.x();
        gx(i+system_size) = pt.y();
    }

    for (auto& cl : msh)
    {
        auto fcs = faces(msh, cl);
        auto pts = points(msh, cl);
        Eigen::Matrix<T,3,3> M = Eigen::Matrix<T,3,3>::Zero();
        M(0,0) = measure(msh, fcs[0]);
        M(1,1) = measure(msh, fcs[1]);
        M(2,2) = measure(msh, fcs[2]);

        auto K = - G.transpose() * M * G;

        auto ptids = cl.point_ids();

        for (size_t i = 0; i < K.rows(); i++)
        {
            if ( dirichlet_nodes.at(ptids[i]) )
                continue;

            for (size_t j = 0; j < K.cols(); j++)
            {
                if ( dirichlet_nodes.at(ptids[j]) )
                {
                    auto rvalx = -K(i,j)*pts[i].x();
                    auto rvaly = -K(i,j)*pts[i].y();
                    gb(compress_map.at(ptids[i])) += rvalx;
                    gb(system_size + compress_map.at(ptids[i])) += rvaly;
                    continue;
                }
                gtA.push_back( triplet_type(compress_map.at(ptids[i]),
                               compress_map.at(ptids[j]),
                               K(i,j)) );

                gtA.push_back( triplet_type(system_size + compress_map.at(ptids[i]),
                               system_size + compress_map.at(ptids[j]),
                               K(i,j)) );
            }
        }
    }

    gA.setFromTriplets(gtA.begin(), gtA.end());

    disk::solvers::conjugated_gradient_params<T> cgp;
    cgp.verbose = true;
    disk::solvers::conjugated_gradient(cgp, gA, gb, gx);



    for (size_t i = 0; i < system_size; i++)
    {
        point_type pt = point_type( {gx(i), gx(system_size+i)} );
        backend_storage->points.at( expand_map.at(i) ) = pt;
    }
}


} // namespace_disk
