/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
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

#ifndef _GEOMETRY_CARTESIAN_HPP_
#define _GEOMETRY_CARTESIAN_HPP_

#include <list>

#include "diskpp/geometry/element_hexahedral.hpp"

namespace disk {

template<size_t DIM>
struct cartesian_storage_class {
    static_assert(DIM == 2 or DIM == 3, "Only DIM = 2 or DIM = 3 for cartesian meshes");
};

template<>
struct cartesian_storage_class<3> {
    typedef cartesian_element<3,0>  volume_type;
    typedef cartesian_element<3,1>  surface_type;
    typedef cartesian_element<3,2>  edge_type;
    typedef cartesian_element<3,3>  node_type;
};

template<>
struct cartesian_storage_class<2> {
    typedef cartesian_element<2,0>  surface_type;
    typedef cartesian_element<2,1>  edge_type;
    typedef cartesian_element<2,2>  node_type;
};

/**
 * @brief speciliaziton for the storage class of a cartesian mesh
 * This class is more efficiant than generic_mesh for cartesain meshes since the internal data structures
 * are optimized
 *
 * @tparam T scalar type
 * @tparam DIM dimension of the mesh, i.e, 1D, 2D or 3D
 */

template<typename T, size_t DIM>
using cartesian_mesh_storage = mesh_storage<T, DIM, cartesian_storage_class<DIM>>;

/**
 * @brief speciliaziton for  cartesian mesh
 * This class is more efficiant than generic_mesh for cartesain meshes since the internal data structures
 * are optimized
 *
 * @tparam T scalar type
 * @tparam DIM dimension of the mesh, i.e, 1D, 2D or 3D
 */
template<typename T, size_t DIM>
using cartesian_mesh = mesh<T, DIM, cartesian_mesh_storage<T, DIM>>;

/**
 * @brief Return the list of point of the specified element (face or cell)
 * 
 * @tparam T scalar type
 * @tparam DIM dimension of the mesh
 * @tparam CODIM co-dimension of the element (CODIM=DIM if elem=cell and CODIM=DIM-1 if elem=face)
 * @param msh mesh
 * @param elem specified element (cell or face)
 * @return std::array<typename cartesian_mesh<T,DIM>::point_type, cartesian_priv::howmany<DIM, CODIM>::nodes> list of points of the specified element
 */
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

    std::array<typename cartesian_mesh<T, DIM>::point_type, cartesian_priv::howmany<DIM, CODIM>::nodes> pts;
    std::transform(ptids.begin(), ptids.end(), pts.begin(), ptid_to_point);

    return pts;
}

/**
 * @brief Return the number of faces of the specified cell
 *
 * @param msh mesh
 * @param cl pecified cell
 * @return template<typename T, size_t DIM> constexpr howmany_faces number of faces
 */
template<typename T, size_t DIM>
constexpr size_t
howmany_faces(const cartesian_mesh<T,DIM>& msh, const typename cartesian_mesh<T,DIM>::cell& cl)
{
    static_assert(DIM == 1 or DIM == 2 or DIM == 3, "wrong dimension");
    switch(DIM)
    {
        case 1: return 2;
        case 2: return 4;
        case 3: return 6;
    }
    /* NOTREACHED */
    assert(false);
}

/**
 * @brief Return the faces of the specified 2D-cell (quadrangle)
 *
 * @tparam T scalar_type
 * @param msh mesh
 * @param cl specified 2D-cell
 * @return std::array<typename cartesian_mesh<T, 2>::face, 4> faces of the cell
 */
template<typename T>
std::array<typename cartesian_mesh<T, 2>::face, 4>
faces(const cartesian_mesh<T, 2>& msh,
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

/**
 * @brief Return the faces of the specified 3D-cell (hexahedron)
 *
 * @tparam T scalar_type
 * @param msh mesh
 * @param cl specified 3D-cell
 * @return std::array<typename cartesian_mesh<T, 3>::face, 6> faces of the cell
 */
template<typename T>
std::array<typename cartesian_mesh<T, 3>::face, 6>
faces(const cartesian_mesh<T, 3>& msh,
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

/**
 * @brief Return the faces id of the specified 2D-cell (quadrangle)
 *
 * @tparam T scalar_type
 * @param msh mesh
 * @param cl specified 2D-cell
 * @return std::array<typename cartesian_mesh<T, 2>::face::id_type, 4> faces of the cell
 */
template<typename T>
std::array<typename cartesian_mesh<T, 2>::face::id_type, 4>
faces_id(const cartesian_mesh<T, 2>& msh, const typename cartesian_mesh<T, 2>::cell& cl)
{
    typedef typename cartesian_mesh<T, 2>::face face_type;
    std::array<typename face_type::id_type, 4>  ret;

    auto ptids = cl.point_ids();
    assert(ptids.size() == 4);

    ret[0] = msh.lookup(face_type({ptids[0], ptids[1]}));
    ret[1] = msh.lookup(face_type({ptids[0], ptids[2]}));
    ret[2] = msh.lookup(face_type({ptids[1], ptids[3]}));
    ret[3] = msh.lookup(face_type({ptids[2], ptids[3]}));

    return ret;
}

/**
 * @brief Return the faces id of the specified 3D-cell (hexahedron)
 *
 * @tparam T scalar_type
 * @param msh mesh
 * @param cl specified 3D-cell
 * @return std::array<typename cartesian_mesh<T, 3>::face::id_type, 6> faces of the cell
 */
template<typename T>
std::array<typename cartesian_mesh<T, 3>::face::id_type, 6>
faces_id(const cartesian_mesh<T, 3>& msh, const typename cartesian_mesh<T, 3>::cell& cl)
{
    typedef typename cartesian_mesh<T, 3>::face face_type;
    std::array<typename face_type::id_type, 6>   ret;

    auto ptids = cl.point_ids();
    assert(ptids.size() == 8);

    ret[0] = msh.lookup(face_type({ptids[0], ptids[2], ptids[6], ptids[4]}));
    ret[1] = msh.lookup(face_type({ptids[1], ptids[3], ptids[7], ptids[5]}));
    ret[2] = msh.lookup(face_type({ptids[0], ptids[1], ptids[3], ptids[2]}));
    ret[3] = msh.lookup(face_type({ptids[4], ptids[5], ptids[7], ptids[6]}));
    ret[4] = msh.lookup(face_type({ptids[0], ptids[4], ptids[5], ptids[1]}));
    ret[5] = msh.lookup(face_type({ptids[2], ptids[6], ptids[7], ptids[3]}));

    return ret;
}

/**
 * @brief Compute the measure (volume in 3D or area 2D) of the specified cell
 * 
 * @tparam T scalar type
 * @tparam DIM dimension of the mesh
 * @param msh mesh
 * @param cl specified cell
 * @return T measure of the cell
 */
template<typename T>
T
measure(const cartesian_mesh<T,2>& msh,
        const typename cartesian_mesh<T,2>::cell& cl)
{
    auto pts = points(msh, cl);
    assert(pts.size() == 4);
    auto v0 = (pts[1] - pts[0]).to_vector().norm();
    auto v1 = (pts[2] - pts[0]).to_vector().norm();
    return v0 * v1;
}

template<typename T>
T
measure(const cartesian_mesh<T,3>& msh,
        const typename cartesian_mesh<T,3>::cell& cl)
{
    auto pts = points(msh, cl);
    assert(pts.size() == 4);
    auto v0 = (pts[1] - pts[0]).to_vector().norm();
    auto v1 = (pts[2] - pts[0]).to_vector().norm();
    auto v2 = (pts[4] - pts[0]).to_vector().norm();
    return v0 * v1 * v2;
}

/**
 * @brief Compute the measure (area in 3D or length 2D) of the specified face
 *
 * @tparam T scalar type
 * @tparam DIM dimension of the mesh
 * @param msh mesh
 * @param fc specified face
 * @return T measure of the face
 */
template<typename T>
T
measure(const cartesian_mesh<T,2>& msh,
        const typename cartesian_mesh<T,2>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 2);
    return (pts[1] - pts[0]).to_vector().norm();
}

template<typename T>
T
measure(const cartesian_mesh<T,3>& msh,
        const typename cartesian_mesh<T,3>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 4);
    auto v0 = (pts[1] - pts[0]).to_vector().norm();
    auto v1 = (pts[3] - pts[0]).to_vector().norm();
    return v0 * v1;
}

} // namespace disk

#endif /* _GEOMETRY_CARTESIAN_HPP_ */
