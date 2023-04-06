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

#ifndef _GEOMETRY_SIMPLICIAL_HPP_
#define _GEOMETRY_SIMPLICIAL_HPP_

#include <set>
#include <list>

#include "diskpp/geometry/element_simplicial.hpp"
#include "diskpp/loaders/loader_utils.hpp"
#include "diskpp/common/simplicial_formula.hpp"
#include "diskpp/mesh/rational.hpp"

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

/**
 * @brief Return the list of point of the specified element (face or cell)
 *
 * @tparam T scalar type
 * @tparam DIM dimension of the mesh
 * @tparam CODIM co-dimension of the element (CODIM=DIM if elem=cell and CODIM=DIM-1 if elem=face)
 * @param msh (simplicial) mesh
 * @param elem specified element (cell or face)
 * @return std::array<typename cartesian_mesh<T,DIM>::point_type, cartesian_priv::howmany<DIM, CODIM>::nodes> list of
 * points of the specified element
 */
template<typename T, size_t DIM, size_t CODIM>
std::array<typename simplicial_mesh<T,DIM>::point_type, priv::howmany<DIM, CODIM>::nodes>
points(const simplicial_mesh<T,DIM>& msh, const simplicial_element<DIM, CODIM>& elem)
{
    auto ptids = elem.point_ids();

    auto points_begin = msh.points_begin();

    auto ptid_to_point = [&](const point_identifier<DIM>& pi) -> typename simplicial_mesh<T,DIM>::point_type {
        return *std::next(points_begin, pi);
    };

    std::array<typename simplicial_mesh<T,DIM>::point_type, priv::howmany<DIM, CODIM>::nodes> pts;
    std::transform(ptids.begin(), ptids.end(), pts.begin(), ptid_to_point);

    return pts;
}

/**
 * @brief Return the number of faces of the specified cell
 *
 * @param msh (simplicial) mesh
 * @param cl specified cell
 * @return template<typename T, size_t DIM> constexpr howmany_faces number of faces
 */
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

/**
 * @brief Return the faces (triangles) of the specified 3D-cell (tetrahedra)
 *
 * @tparam T scalar_type
 * @param msh (simplicial) mesh
 * @param cl specified 3D-cell (tetrahedra)
 * @return std::array<typename simplicial_mesh<T, 3>::face, 4> faces of the cell
 */
template<typename T>
std::array<typename simplicial_mesh<T, 3>::face, 4>
faces(const simplicial_mesh<T, 3>& msh,
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

/**
 * @brief Return the faces (edges) of the specified 2D-cell (triangle)
 *
 * @tparam T scalar_type
 * @param msh (simplicial) mesh
 * @param cl specified 2D-cell (triangle)
 * @return std::array<typename simplicial_mesh<T, 2>::face, 3> faces of the cell
 */
template<typename T>
std::array<typename simplicial_mesh<T, 2>::face, 3>
faces(const simplicial_mesh<T, 2>& msh,
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

/**
 * @brief Return the faces id (triangles) of the specified 3D-cell (tetrahedra)
 *
 * @tparam T scalar_type
 * @param msh (simplicial) mesh
 * @param cl specified 3D-cell (tetrahedra)
 * @return std::array<typename simplicial_mesh<T, 3>::face::id_type, 4> faces id of the cell
 */
template<typename T>
std::array<typename simplicial_mesh<T, 3>::face::id_type, 4>
faces_id(const simplicial_mesh<T, 3>& msh, const typename simplicial_mesh<T, 3>::cell& cl)
{
    typedef typename simplicial_mesh<T, 3>::face face_type;
    std::array<typename face_type::id_type, 4>   ret;

    auto ptids = cl.point_ids();
    assert(ptids.size() == 4);

    ret[0] = msh.lookup(face_type({ptids[1], ptids[2], ptids[3]}));
    ret[1] = msh.lookup(face_type({ptids[0], ptids[2], ptids[3]}));
    ret[2] = msh.lookup(face_type({ptids[0], ptids[1], ptids[3]}));
    ret[3] = msh.lookup(face_type({ptids[0], ptids[1], ptids[2]}));

    return ret;
}

/**
 * @brief Return the faces id (edges) of the specified 2D-cell (triangle)
 *
 * @tparam T scalar_type
 * @param msh (simplicial) mesh
 * @param cl specified 2D-cell (triangle)
 * @return std::array<typename simplicial_mesh<T, 2>::face::id_type, 3> faces id of the cell
 */
template<typename T>
std::array<typename simplicial_mesh<T, 2>::face::id_type, 3>
faces_id(const simplicial_mesh<T, 2>& msh, const typename simplicial_mesh<T, 2>::cell& cl)
{
    typedef typename simplicial_mesh<T, 2>::face face_type;
    std::array<typename face_type::id_type, 3>   ret;

    auto ptids = cl.point_ids();
    assert(ptids.size() == 3);

    assert(ptids[0] < ptids[1] && ptids[1] < ptids[2]);

    ret[0] = msh.lookup(face_type({ptids[0], ptids[1]}));
    ret[1] = msh.lookup(face_type({ptids[1], ptids[2]}));
    ret[2] = msh.lookup(face_type({ptids[0], ptids[2]}));

    return ret;
}

/**
 * @brief Compute the volume of the specified 3D-cell (tetrahedra)
 *
 * @tparam T scalar type
 * @param msh (simplicial) mesh
 * @param cl specified 3D-cell (tetrahedra)
 * @param signed_volume return a signed volume
 * @return T volume of the cell
 */
template<typename T>
T
measure(const simplicial_mesh<T, 3>& msh,
        const typename simplicial_mesh<T, 3>::cell& cl,
        bool signed_volume = false)
{
    auto pts = points(msh, cl);
    assert(pts.size() == 4);

    const T vol = volume_tetrahedron_kahan(pts[0], pts[1], pts[2], pts[3]);

    if (signed_volume)
    {
        const auto v0   = (pts[1] - pts[0]).to_vector();
        const auto v1   = (pts[2] - pts[0]).to_vector();
        const auto v2   = (pts[3] - pts[0]).to_vector();

        const T vol2 = v0.dot(v1.cross(v2)) / T(6);
        return vol2/std::abs(vol2) * vol;
    }

    return vol;
}

template<typename T>
rational<T>
measure(const simplicial_mesh<rational<T>, 3>& msh, const typename simplicial_mesh<rational<T>, 3>::cell& cl, bool signed_volume = false)
{
    auto pts = points(msh, cl);
    assert(pts.size() == 4);

    const auto v0 = (pts[1] - pts[0]).to_vector();
    const auto v1 = (pts[2] - pts[0]).to_vector();
    const auto v2 = (pts[3] - pts[0]).to_vector();

    const auto vol = v0.dot(v1.cross(v2)) / T(6);

    if (signed_volume)
    {
        return vol;
    }

    return abs(vol);
}

/**
 * @brief Compute the area of the specified 3D-face (triangle)
 *
 * @tparam T scalar type
 * @param msh (simplicial) mesh
 * @param fc specified 3D-face (triangle)
 * @return T area of the face
 */
template<typename T>
T
measure(const simplicial_mesh<T, 3>& msh,
        const typename simplicial_mesh<T, 3>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 3);

    return area_triangle_kahan(pts[0], pts[1], pts[2]);
}

/**
 * @brief Compute the area of the specified 2D-cell (triangle)
 *
 * @tparam T scalar type
 * @param msh (simplicial) mesh
 * @param cl specified 2D-cell (triangle)
 * @param signed_volume return a signed volume
 * @return T volume of the cell
 */
template<typename T>
T
measure(const simplicial_mesh<T, 2>& msh,
        const typename simplicial_mesh<T, 2>::cell& cl,
        bool signed_volume = false)
{
    using namespace std;
    auto pts = points(msh, cl);
    assert(pts.size() == 3);

    const T area = area_triangle_kahan(pts[0], pts[1], pts[2]);

    if (signed_volume)
    {
        const auto d0 = (pts[1] - pts[0]);
        const auto d1 = (pts[2] - pts[0]);

        const T area2 = d0.x() * d1.y() - d1.x() * d0.y();

        return area2/std::abs(area2) * area;
    }

    return area;
}

template<typename T>
rational<T>
measure(const simplicial_mesh<rational<T>, 2>&                msh,
        const typename simplicial_mesh<rational<T>, 2>::cell& cl,
        bool                                                  signed_volume = false)
{
    using namespace std;
    auto pts = points(msh, cl);
    assert(pts.size() == 3);

    const auto d0 = (pts[1] - pts[0]);
    const auto d1 = (pts[2] - pts[0]);

    const auto area = d0.x() * d1.y() - d1.x() * d0.y();

    if (signed_volume)
    {
        return area;
    }

    return abs(area);
}

/**
 * @brief Compute the length of the specified 2D-face (edge)
 *
 * @tparam T scalar type
 * @param msh (simplicial) mesh
 * @param fc specified 2D-face (edge)
 * @return T length of the face
 */
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
bool
is_inside(const simplicial_mesh<T,3>& msh,
          const typename simplicial_mesh<T,3>::cell_type& cl,
          const typename simplicial_mesh<T,3>::point_type& pt)
{
    auto pts = points(msh, cl);
    auto v0 = to_vector(pts[1] - pts[0]);
    auto v1 = to_vector(pts[2] - pts[0]);
    auto v2 = to_vector(pts[2] - pts[1]);
    auto v3 = to_vector(pts[3] - pts[0]);
    auto v4 = to_vector(pts[3] - pts[1]);
    auto v5 = to_vector(pts[3] - pts[2]);

    int count = 0;

    T t0 = v3.cross(v0).dot( to_vector(pt - pts[0]) );
    if (t0 < 0) count--; else count++;

    T t1 = v0.cross(v1).dot( to_vector(pt - pts[0]) );
    if (t1 < 0) count--; else count++;
    
    T t2 = v4.cross(v2).dot( to_vector(pt - pts[1]) );
    if (t2 < 0) count--; else count++;
    
    T t3 = v1.cross(v3).dot( to_vector(pt - pts[0]) );
    if (t3 < 0) count--; else count++;

    return (count == 4) or (count == -4);
}


} // namespace disk

#endif /* _GEOMETRY_SIMPLICIAL_HPP_ */
