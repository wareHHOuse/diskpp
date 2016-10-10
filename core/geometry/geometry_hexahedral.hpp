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

#ifndef _GEOMETRY_HEXAHEDRAL_HPP_
#define _GEOMETRY_HEXAHEDRAL_HPP_

#include <list>

#include "geometry/element_hexahedral.hpp"

namespace disk {

struct storage_class_hexahedral;

template<size_t DIM>
struct storage_class_trait<storage_class_hexahedral, DIM> {
    static_assert(DIM == 3, "Only DIM = 3 for hexahedra");
};

template<>
struct storage_class_trait<storage_class_hexahedral, 3> {
    typedef hexahedral_element<3,0>    volume_type;
    typedef hexahedral_element<3,1>    surface_type;
    typedef hexahedral_element<3,2>    edge_type;
    typedef hexahedral_element<3,3>    node_type;
};

template<typename T, size_t DIM>
using hexahedral_mesh_storage = mesh_storage<T, DIM, storage_class_hexahedral>;

template<typename T>
using hexahedral_mesh = mesh<T, 3, hexahedral_mesh_storage<T, 3>>;


template<typename T, size_t CODIM>
std::array<typename hexahedral_mesh<T>::point_type, hex_priv::howmany<3, CODIM>::nodes>
points(const hexahedral_mesh<T>& msh, const hexahedral_element<3, CODIM>& elem)
{
    auto ptids = elem.point_ids();

    auto ptid_to_point = [&](const point_identifier<3>& pi) -> auto {
        auto itor = msh.points_begin();
        std::advance(itor, pi);
        return *itor;
    };

    typedef point_identifier<3>       point_id_type;

    std::array<typename hexahedral_mesh<T>::point_type, hex_priv::howmany<3, CODIM>::nodes> pts;
    std::transform(ptids.begin(), ptids.end(), pts.begin(), ptid_to_point);

    return pts;
}


/* Return the number of elements of the specified cell */
template<typename T, size_t CODIM>
size_t
number_of_faces(const hexahedral_mesh<T>& msh, const hexahedral_element<3, CODIM>& e)
{
    switch(CODIM)
    {
        case 0: return 6;
        case 1: return 4;
        case 2: return 2;
        case 3: return 0;
    }
    /* NOTREACHED */
}


template<typename T>
std::array<typename hexahedral_mesh<T>::face, 6>
faces(const hexahedral_mesh<T>&,
      const typename hexahedral_mesh<T>::cell& cl)
{
    typedef typename hexahedral_mesh<T>::face    face_type;
    std::array<face_type, 6> ret;

    auto ptids = cl.point_ids();
    assert(ptids.size() == 8);

    ret[0] = face_type( { ptids[0], ptids[3], ptids[5], ptids[2] } );
    ret[1] = face_type( { ptids[1], ptids[4], ptids[7], ptids[6] } );
    ret[2] = face_type( { ptids[0], ptids[1], ptids[6], ptids[3] } );
    ret[3] = face_type( { ptids[2], ptids[5], ptids[7], ptids[4] } );
    ret[4] = face_type( { ptids[0], ptids[2], ptids[4], ptids[1] } );
    ret[5] = face_type( { ptids[3], ptids[6], ptids[7], ptids[5] } );

    return ret;
}

/* WARNING: ONLY CARTESIAN MESHES FOR NOW!! */
template<typename T>
T
measure(const hexahedral_mesh<T>& msh,
        const typename hexahedral_mesh<T>::cell& cl)
{
    auto pts = points(msh, cl);
    assert(pts.size() == 6);

    auto v0 = (pts[1] - pts[0]).to_vector().norm();
    auto v1 = (pts[2] - pts[0]).to_vector().norm();
    auto v2 = (pts[3] - pts[0]).to_vector().norm();

    return v0 * v1 * v2;
}


template<typename T>
T
measure(const hexahedral_mesh<T>& msh,
        const typename hexahedral_mesh<T>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 4);

    auto v0 = (pts[1] - pts[0]).to_vector();
    auto v1 = (pts[2] - pts[1]).to_vector();
    auto a1 = v0.cross(v1).norm();

    auto v2 = (pts[3] - pts[2]).to_vector();
    auto v3 = (pts[0] - pts[3]).to_vector();
    auto a2 = v2.cross(v3).norm();

    return (a1 + a2)*T(0.5);
}

template<typename T>
static_vector<T, 3>
normal(const hexahedral_mesh<T>& msh,
       const typename hexahedral_mesh<T>::cell& cl,
       const typename hexahedral_mesh<T>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 8);

    auto v0 = (pts[1] - pts[0]).to_vector();
    auto v1 = (pts[2] - pts[1]).to_vector();
    auto n = v0.cross(v1);

    auto cell_bar = barycenter(msh, cl);
    auto face_bar = barycenter(msh, fc);
    auto outward_vector = (face_bar - cell_bar).to_vector();

    if ( n.dot(outward_vector) < T(0) ) /* should be always positive */
        return -n/n.norm();

    return n/n.norm();
}

} // namespace disk

#endif /* _GEOMETRY_HEXAHEDRAL_HPP_ */
