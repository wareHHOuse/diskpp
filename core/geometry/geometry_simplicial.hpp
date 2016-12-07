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

 #ifndef _GEOMETRY_HPP_WAS_INCLUDED_
     #error "You must NOT include this file directly. Include geometry.hpp."
 #endif

#ifndef _GEOMETRY_SIMPLICIAL_HPP_
#define _GEOMETRY_SIMPLICIAL_HPP_

#include <set>
#include <list>

#include "geometry/element_simplicial.hpp"
#include "loaders/loader_utils.hpp"

namespace disk {

struct storage_class_simplicial;

template<size_t DIM>
struct storage_class_trait<storage_class_simplicial, DIM> {
    static_assert(DIM > 0 && DIM <= 3, "element_types: CODIM must be less than DIM");
};

template<>
struct storage_class_trait<storage_class_simplicial, 1> {
    typedef simplicial_element<1,0>    edge_type;
    typedef simplicial_element<1,1>    node_type;
};

template<>
struct storage_class_trait<storage_class_simplicial, 2> {
    typedef simplicial_element<2,0>    surface_type;
    typedef simplicial_element<2,1>    edge_type;
    typedef simplicial_element<2,2>    node_type;
};

template<>
struct storage_class_trait<storage_class_simplicial, 3> {
    typedef simplicial_element<3,0>    volume_type;
    typedef simplicial_element<3,1>    surface_type;
    typedef simplicial_element<3,2>    edge_type;
    typedef simplicial_element<3,3>    node_type;
};

template<typename T, size_t DIM>
using simplicial_mesh_storage = mesh_storage<T, DIM, storage_class_simplicial>;

template<typename T, size_t DIM>
using simplicial_mesh = mesh<T, DIM, simplicial_mesh_storage<T, DIM>>;

template<typename T>
using tetrahedral_mesh = simplicial_mesh<T, 3>;

template<typename T>
using triangular_mesh = simplicial_mesh<T, 2>;


template<typename T, size_t DIM, size_t CODIM>
std::array<typename simplicial_mesh<T,DIM>::point_type, priv::howmany<DIM, CODIM>::nodes>
points(const simplicial_mesh<T,DIM>& msh, const simplicial_element<DIM, CODIM>& elem)
{
    auto ptids = elem.point_ids();

    auto ptid_to_point = [&](const point_identifier<DIM>& pi) -> auto {
        auto itor = msh.points_begin();
        std::advance(itor, pi);
        return *itor;
    };

    typedef point_identifier<DIM>       point_id_type;

    std::array<typename simplicial_mesh<T,DIM>::point_type, priv::howmany<DIM, CODIM>::nodes> pts;
    std::transform(ptids.begin(), ptids.end(), pts.begin(), ptid_to_point);

    return pts;
}

/* Return the number of elements of the specified cell */
template<typename T, size_t DIM>
size_t
number_of_faces(const simplicial_mesh<T,DIM>& msh, const typename simplicial_mesh<T,DIM>::cell& cl)
{
    static_assert(DIM == 1 or DIM == 2 or DIM == 3, "wrong dimension");
    switch(DIM)
    {
        case 1: return 2;
        case 2: return 3;
        case 3: return 4;
    }
    /* NOTREACHED */
}

template<typename T>
std::array<typename simplicial_mesh<T, 3>::face, 4>
faces(const simplicial_mesh<T, 3>&,
      const typename simplicial_mesh<T, 3>::cell& cl)
{
    typedef typename simplicial_mesh<T, 3>::face    face_type;
    std::array<face_type, 4> ret;

    auto ptids = cl.point_ids();
    assert(ptids.size() == 4);

    ret[0] = face_type( { ptids[1], ptids[2], ptids[3] } );
    ret[1] = face_type( { ptids[0], ptids[2], ptids[3] } );
    ret[2] = face_type( { ptids[0], ptids[1], ptids[3] } );
    ret[3] = face_type( { ptids[0], ptids[1], ptids[2] } );

    return ret;
}

template<typename T>
std::array<typename simplicial_mesh<T, 2>::face, 3>
faces(const simplicial_mesh<T, 2>&,
      const typename simplicial_mesh<T, 2>::cell& cl)
{
    typedef typename simplicial_mesh<T, 2>::face    face_type;
    std::array<face_type, 3> ret;

    auto ptids = cl.point_ids();
    assert(ptids.size() == 3);

    ret[0] = face_type( { ptids[0], ptids[1] } );
    ret[1] = face_type( { ptids[1], ptids[2] } );
    ret[2] = face_type( { ptids[0], ptids[2] } );

    return ret;
}

template<typename T>
T
measure(const simplicial_mesh<T, 3>& msh,
        const typename simplicial_mesh<T, 3>::cell& cl,
        bool signed_volume = false)
{
    auto pts = points(msh, cl);
    assert(pts.size() == 4);

    auto v0 = (pts[1] - pts[0]).to_vector();
    auto v1 = (pts[2] - pts[0]).to_vector();
    auto v2 = (pts[3] - pts[0]).to_vector();

    if (signed_volume)
        return v0.dot(v1.cross(v2))/T(6);

    return std::abs( v0.dot(v1.cross(v2))/T(6) );
}

template<typename T>
T
measure(const simplicial_mesh<T, 3>& msh,
        const typename simplicial_mesh<T, 3>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 3);

    auto v0 = (pts[1] - pts[0]).to_vector();
    auto v1 = (pts[2] - pts[0]).to_vector();

    return v0.cross(v1).norm()/T(2);
}

template<typename T>
T
measure(const simplicial_mesh<T, 2>& msh,
        const typename simplicial_mesh<T, 2>::cell& cl)
{
    auto pts = points(msh, cl);
    assert(pts.size() == 3);

    auto v0 = (pts[1] - pts[0]).to_vector();
    auto v1 = (pts[2] - pts[0]).to_vector();

    return v0.cross(v1).norm()/T(2);
}

template<typename T>
T
measure(const simplicial_mesh<T, 2>& msh,
        const typename simplicial_mesh<T, 2>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 2);

    return (pts[1] - pts[0]).to_vector().norm();
}

template<typename T>
static_vector<T, 3>
normal(const simplicial_mesh<T,3>& msh,
       const typename simplicial_mesh<T,3>::cell& cl,
       const typename simplicial_mesh<T,3>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 3);

    auto v0 = (pts[1] - pts[0]).to_vector();
    auto v1 = (pts[2] - pts[0]).to_vector();
    auto n = v0.cross(v1);

    auto cell_bar = barycenter(msh, cl);
    auto face_bar = barycenter(msh, fc);
    auto outward_vector = (face_bar - cell_bar).to_vector();

    if ( n.dot(outward_vector) < T(0) )
        return -n/n.norm();

    return n/n.norm();
}

template<typename Mesh>
class submesher;

template<typename T>
class submesher<simplicial_mesh<T,2>>
{
    typedef simplicial_mesh<T,2>                mesh_type;
    typedef typename mesh_type::point_type      point_type;
    typedef typename mesh_type::cell            cell_type;

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
            return ( a.p1 < b.p1 or (a.p1 == b.p1 and a.p0 < b.p0) );
        }

        friend bool operator==(const edge& a, const edge& b) {
            return ( a.p0 == b.p0 and a.p1 == b.p0 );
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
        for (size_t i = 0; i < triangles_to_process; i++)
        {
            /* remove the original triangle */
            triangle t = m_triangles.front();
            m_triangles.pop_front();

            /* find the edges of the triangle */
            auto t_e0 = *std::lower_bound(m_edges.begin(), m_edges.end(), edge(t.p0, t.p1));
            assert(t_e0.is_broken);
            auto t_e1 = *std::lower_bound(m_edges.begin(), m_edges.end(), edge(t.p1, t.p2));
            assert(t_e1.is_broken);
            auto t_e2 = *std::lower_bound(m_edges.begin(), m_edges.end(), edge(t.p0, t.p2));
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
    submesher()
    {}

    mesh_type
    generate_mesh(const mesh_type& msh, const cell_type& cl, size_t refinement_steps = 3)
    {
        std::cout << "genmesh" << std::endl;
        m_points.clear();
        m_edges.clear();
        m_triangles.clear();

        /* Transform the triangle in internal representation */
        auto pts = points(msh, cl);
        m_points.insert(m_points.begin(), pts.begin(), pts.end());

        /* Here the mapping between the local boundary number and the parent
         * cell boundary number is established. */
        auto fcs = faces(msh, cl);
        assert(fcs.size() == 3);

        edge e0(0, 1, msh.lookup(fcs[0]), true);
        m_edges.push_back(e0);

        edge e1(1, 2, msh.lookup(fcs[1]), true);
        m_edges.push_back(e1);

        edge e2(0, 2, msh.lookup(fcs[2]), true);
        m_edges.push_back(e2);

        std::sort(m_edges.begin(), m_edges.end());

        triangle t(0,1,2);
        m_triangles.push_back(t);

        refine( refinement_steps );

        mesh_type refined_element_submesh;
        auto storage = refined_element_submesh.backend_storage();

        storage->points = std::move(m_points);
        storage->nodes.reserve( storage->points.size() );
        for (size_t i = 0; i < storage->points.size(); i++)
        {
            auto point_id = point_identifier<2>( i );
            storage->nodes.push_back( typename mesh_type::node_type({point_id}));
        }

        std::vector<std::pair<typename mesh_type::edge_type, size_t>> be;
        be.reserve( 3 * m_triangles.size() );

        storage->edges.reserve( 3 * m_triangles.size() );
        storage->surfaces.reserve( m_triangles.size() );
        for (auto& t : m_triangles)
        {
            auto pi0 = point_identifier<2>(t.p0);
            auto pi1 = point_identifier<2>(t.p1);
            auto pi2 = point_identifier<2>(t.p2);

            auto e0 = typename mesh_type::edge_type({pi0, pi1});
            auto e1 = typename mesh_type::edge_type({pi1, pi2});
            auto e2 = typename mesh_type::edge_type({pi0, pi2});

            storage->edges.push_back( e0 );
            storage->edges.push_back( e1 );
            storage->edges.push_back( e2 );
            storage->surfaces.push_back( typename mesh_type::surface_type({pi0, pi1, pi2}) );

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

} // namespace disk

#endif /* _GEOMETRY_SIMPLICIAL_HPP_ */
