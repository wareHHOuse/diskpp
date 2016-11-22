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

    struct triangle
    {
        size_t p0, p1, p2;
        bool   b0, b1, b2;
    };

    std::vector<point_type>     m_points;
    std::list<triangle>         m_triangles;

    void refine(size_t steps)
    {
        if (steps == 0)
            return;

        size_t triangles_to_process = m_triangles.size();

        for (size_t i = 0; i < triangles_to_process; i++)
        {
            triangle t = m_triangles.front();
            m_triangles.pop_front();

            auto new_point_pos = m_points.size();

            auto bar0 = (m_points[t.p0] + m_points[t.p1])/2;
            auto bar1 = (m_points[t.p1] + m_points[t.p2])/2;
            auto bar2 = (m_points[t.p0] + m_points[t.p2])/2;

            m_points.push_back(bar0); auto bar0pos = new_point_pos;
            m_points.push_back(bar1); auto bar1pos = new_point_pos+1;
            m_points.push_back(bar2); auto bar2pos = new_point_pos+2;

            triangle t0, t1, t2, t3;

            t0.p0 = t.p0; t0.p1 = bar0pos; t0.p2 = bar2pos;
            t0.b0 = t.b0; t0.b1 = false; t0.b2 = t.b2;
            m_triangles.push_back(t0);

            t1.p0 = bar0pos; t1.p1 = t.p1; t1.p2 = bar1pos;
            t1.b0 = t.b0; t1.b1 = t.b1; t1.b2 = false;
            m_triangles.push_back(t1);

            t2.p0 = bar2pos; t2.p1 = bar1pos; t2.p2 = t.p2;
            t2.b0 = false; t2.b1 = t.b1; t2.b2 = t.b2;
            m_triangles.push_back(t2);

            t3.p0 = bar0pos; t3.p1 = bar1pos; t3.p2 = bar2pos;
            t3.b0 = false; t3.b1 = false; t3.b2 = false;
            m_triangles.push_back(t3);
        }

        refine(steps-1);
    }

public:
    submesher()
    {}

    mesh_type
    generate_mesh(const mesh_type& msh, const cell_type& cl)
    {
        m_points.clear();
        m_triangles.clear();

        auto pts = points(msh, cl);
        m_points.insert(m_points.begin(), pts.begin(), pts.end());

        triangle t;
        t.p0 = 0; t.p1 = 1; t.p2 = 2;
        t.b0 = true; t.b1 = true; t.b2 = true;
        m_triangles.push_back(t);

        refine( 3 );

        mesh_type refined_element_submesh;
        auto storage = refined_element_submesh.backend_storage();

        storage->points = std::move(m_points);
        storage->nodes.reserve( storage->points.size() );
        for (size_t i = 0; i < storage->points.size(); i++)
        {
            auto point_id = point_identifier<2>( i );
            storage->nodes.push_back( typename mesh_type::node_type({point_id}));
        }

        std::vector<typename mesh_type::edge_type> be;
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

            if (t.b0) be.push_back(e0);
            if (t.b1) be.push_back(e1);
            if (t.b2) be.push_back(e2);
        }

        THREAD(edge_thread,
                priv::sort_uniq(storage->edges);
        );

        THREAD(tri_thread,
                priv::sort_uniq(storage->surfaces);
        );

        WAIT_THREAD(edge_thread);
        WAIT_THREAD(tri_thread);

        storage->boundary_edges.resize( storage->edges.size() );
        for (auto& e : be)
        {
            auto ei = find_element_id(storage->edges.begin(),
                                      storage->edges.end(), e);

            if (ei.first == false)
                throw std::logic_error("Edge not found. This is a bug.");

            storage->boundary_edges.at(ei.second) = true;
        }

        storage->statistics();

        return refined_element_submesh;
    }
};

} // namespace disk

#endif /* _GEOMETRY_SIMPLICIAL_HPP_ */
