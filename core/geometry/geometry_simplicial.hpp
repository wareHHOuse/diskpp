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

template<size_t DIM>
struct simplicial_storage_class {
    static_assert(DIM > 0 && DIM <= 3, "element_types: CODIM must be less than DIM");
};

template<>
struct simplicial_storage_class<1> {
    typedef simplicial_element<1,0>    edge_type;
    typedef simplicial_element<1,1>    node_type;
};

template<>
struct simplicial_storage_class<2> {
    typedef simplicial_element<2,0>    surface_type;
    typedef simplicial_element<2,1>    edge_type;
    typedef simplicial_element<2,2>    node_type;
};

template<>
struct simplicial_storage_class<3> {
    typedef simplicial_element<3,0>    volume_type;
    typedef simplicial_element<3,1>    surface_type;
    typedef simplicial_element<3,2>    edge_type;
    typedef simplicial_element<3,3>    node_type;
};

template<typename T, size_t DIM>
using simplicial_mesh_storage = mesh_storage<T, DIM, simplicial_storage_class<DIM>>;

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

    auto points_begin = msh.points_begin();

    auto ptid_to_point = [&](const point_identifier<DIM>& pi) -> typename simplicial_mesh<T,DIM>::point_type {
        return *std::next(points_begin, pi);
    };

    typedef point_identifier<DIM>       point_id_type;

    std::array<typename simplicial_mesh<T,DIM>::point_type, priv::howmany<DIM, CODIM>::nodes> pts;
    std::transform(ptids.begin(), ptids.end(), pts.begin(), ptid_to_point);

    return pts;
}

/* Return the number of elements of the specified cell */
template<typename T, size_t DIM>
constexpr size_t
howmany_faces(const simplicial_mesh<T,DIM>& msh, const typename simplicial_mesh<T,DIM>::cell& cl)
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

    assert (ptids[0] < ptids[1] && ptids[1] < ptids[2]);

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
    using namespace std;
    auto pts = points(msh, cl);
    assert(pts.size() == 3);

    auto d0 = (pts[1] - pts[0]);
    auto d1 = (pts[2] - pts[0]);
    
    return abs(d0.x()*d1.y() - d1.x()*d0.y())/T(2);
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

} // namespace disk

#endif /* _GEOMETRY_SIMPLICIAL_HPP_ */
