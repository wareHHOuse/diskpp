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

#ifndef _GEOMETRY_CARTESIAN_HPP_
#define _GEOMETRY_CARTESIAN_HPP_

#include <list>

#include "geometry/element_hexahedral.hpp"

namespace disk {

struct storage_class_cartesian;

template<size_t DIM>
struct storage_class_trait<storage_class_cartesian, DIM> {
    static_assert(DIM == 2 or DIM == 3, "Only DIM = 2 or DIM = 3 for cartesian meshes");
};

template<>
struct storage_class_trait<storage_class_cartesian, 3> {
    typedef cartesian_element<3,0>  volume_type;
    typedef cartesian_element<3,1>  surface_type;
    typedef cartesian_element<3,2>  edge_type;
    typedef cartesian_element<3,3>  node_type;
};

template<>
struct storage_class_trait<storage_class_cartesian, 2> {
    typedef cartesian_element<2,0>  surface_type;
    typedef cartesian_element<2,1>  edge_type;
    typedef cartesian_element<2,2>  node_type;
};

template<typename T, size_t DIM>
using cartesian_mesh_storage = mesh_storage<T, DIM, storage_class_cartesian>;

template<typename T, size_t DIM>
using cartesian_mesh = mesh<T, DIM, cartesian_mesh_storage<T, DIM>>;


template<typename T, size_t DIM, size_t CODIM>
std::array<typename cartesian_mesh<T,DIM>::point_type, cartesian_priv::howmany<DIM, CODIM>::nodes>
points(const cartesian_mesh<T, DIM>& msh, const cartesian_element<DIM, CODIM>& elem)
{
    auto ptids = elem.point_ids();

    auto ptid_to_point = [&](const point_identifier<DIM>& pi) -> auto {
        auto itor = msh.points_begin();
        std::advance(itor, pi);
        return *itor;
    };

    typedef point_identifier<DIM>       point_id_type;

    std::array<typename cartesian_mesh<T, DIM>::point_type, cartesian_priv::howmany<DIM, CODIM>::nodes> pts;
    std::transform(ptids.begin(), ptids.end(), pts.begin(), ptid_to_point);

    return pts;
}


/* Return the number of elements of the specified cell */
template<typename T, size_t DIM, size_t CODIM>
constexpr size_t
number_of_faces(const cartesian_mesh<T, DIM>& msh, const cartesian_element<DIM, CODIM>& e)
{
    static_assert(DIM == 2 or DIM == 3, "Wrong dimension for cartesian mesh");

    if (DIM == 2)
    {
        switch(CODIM)
        {
            case 0: return 4;
            case 1: return 2;
            case 2: return 0;
        }
    }

    if (DIM == 3)
    {
        switch(CODIM)
        {
            case 0: return 6;
            case 1: return 4;
            case 2: return 2;
            case 3: return 0;
        }
    }
    /* NOTREACHED */
}

template<typename T>
std::array<typename cartesian_mesh<T, 2>::face, 4>
faces(const cartesian_mesh<T, 2>&,
      const typename cartesian_mesh<T, 2>::cell& cl)
{
    typedef typename cartesian_mesh<T, 2>::face    face_type;
    std::array<face_type, 4> ret;

    auto ptids = cl.point_ids();
    assert(ptids.size() == 4);

    ret[0] = face_type( { ptids[0], ptids[1] } );
    ret[1] = face_type( { ptids[0], ptids[2] } );
    ret[2] = face_type( { ptids[1], ptids[3] } );
    ret[3] = face_type( { ptids[2], ptids[3] } );

    return ret;
}

template<typename T>
std::array<typename cartesian_mesh<T, 3>::face, 6>
faces(const cartesian_mesh<T, 3>&,
      const typename cartesian_mesh<T, 3>::cell& cl)
{
    typedef typename cartesian_mesh<T, 3>::face    face_type;
    std::array<face_type, 6> ret;

    auto ptids = cl.point_ids();
    assert(ptids.size() == 8);

    ret[0] = face_type( { ptids[0], ptids[2], ptids[6], ptids[4] } );
    ret[1] = face_type( { ptids[1], ptids[3], ptids[7], ptids[5] } );
    ret[2] = face_type( { ptids[0], ptids[1], ptids[3], ptids[2] } );
    ret[3] = face_type( { ptids[4], ptids[5], ptids[7], ptids[6] } );
    ret[4] = face_type( { ptids[0], ptids[4], ptids[5], ptids[1] } );
    ret[5] = face_type( { ptids[2], ptids[6], ptids[7], ptids[3] } );

    return ret;
}

template<typename T, size_t DIM>
T
measure(const cartesian_mesh<T, DIM>& msh,
        const typename cartesian_mesh<T, DIM>::cell& cl)
{
    static_assert(DIM == 2 or DIM == 3, "Wrong dimension for cartesian mesh");

    auto pts = points(msh, cl);

    if (DIM == 2)
    {
        auto v0 = (pts[1] - pts[0]).to_vector().norm();
        auto v1 = (pts[2] - pts[0]).to_vector().norm();

        return v0 * v1;
    }

    if (DIM == 3)
    {
        auto v0 = (pts[1] - pts[0]).to_vector().norm();
        auto v1 = (pts[2] - pts[0]).to_vector().norm();
        auto v2 = (pts[4] - pts[0]).to_vector().norm();

        return v0 * v1 * v2;
    }

    /* NOTREACHED */
}


template<typename T, size_t DIM>
T
measure(const cartesian_mesh<T, DIM>& msh,
        const typename cartesian_mesh<T, DIM>::face& fc)
{
    static_assert(DIM == 2 or DIM == 3, "Wrong dimension for cartesian mesh");

    auto pts = points(msh, fc);

    if (DIM == 2)
        return (pts[1] - pts[0]).to_vector().norm();

    if (DIM == 3)
    {
        auto v0 = (pts[1] - pts[0]).to_vector().norm();
        auto v1 = (pts[3] - pts[0]).to_vector().norm();

        return v0 * v1;
    }

    /* NOTREACHED */
}

} // namespace disk

#endif /* _GEOMETRY_CARTESIAN_HPP_ */
