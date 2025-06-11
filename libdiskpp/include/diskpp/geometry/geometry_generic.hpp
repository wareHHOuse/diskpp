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

#ifndef _GEOMETRY_GENERIC_HPP_
#define _GEOMETRY_GENERIC_HPP_

#include "diskpp/geometry/element_generic.hpp"
#include "diskpp/quadratures/bits/raw_simplices.hpp"
#include "diskpp/common/simplicial_formula.hpp"


namespace disk {

template<size_t DIM>
struct generic_storage_class;

template<size_t DIM>
struct generic_storage_class {
    static_assert(DIM > 0 && DIM <= 3, "This storage class supports DIM between 0 and 3");
};

template<>
struct generic_storage_class<1> {
    typedef generic_element<1,0>    edge_type;
    typedef generic_element<1,1>    node_type;
};

template<>
struct generic_storage_class<2> {
    typedef generic_element<2,0>    surface_type;
    typedef generic_element<2,1>    edge_type;
    typedef generic_element<2,2>    node_type;
};

template<>
struct generic_storage_class<3> {
        typedef generic_element<3,0>    volume_type;
        typedef generic_element<3,1>    surface_type;
        typedef generic_element<3,2>    edge_type;
        typedef generic_element<3,3>    node_type;
};

template<typename T, size_t DIM>
using generic_mesh_storage = mesh_storage<T, DIM, generic_storage_class<DIM>>;

template<typename T, size_t DIM>
using generic_mesh = mesh<T, DIM, generic_mesh_storage<T, DIM>>;

} // namespace disk

#include "geometry_generic_triangulations.hpp"

namespace disk {

template<typename T, size_t DIM>
size_t
howmany_faces(const generic_mesh<T,DIM>& msh,
              const typename generic_mesh<T,DIM>::cell& cl)
{
    return cl.subelement_size();
}


template<typename T, size_t DIM>
std::vector<typename generic_mesh<T, DIM>::face>
faces(const generic_mesh<T, DIM>& msh,
      const typename generic_mesh<T, DIM>::cell& cl)
{
    auto faces_begin = msh.faces_begin();
    auto id_to_face  = [&](const typename generic_mesh<T, DIM>::face::id_type& id) -> auto
    {
        return *std::next(faces_begin, id);
    };

    std::vector<typename generic_mesh<T, DIM>::face> ret;
    ret.resize(cl.subelement_size());

    std::transform(cl.subelement_id_begin(), cl.subelement_id_end(), ret.begin(), id_to_face);

    return ret;
}

template<typename T, size_t DIM>
std::vector<typename generic_mesh<T, DIM>::face::id_type>
faces_id(const generic_mesh<T, DIM>& msh, const typename generic_mesh<T, DIM>::cell& cl)
{
    return cl.faces_ids();
}

template<typename T>
std::vector<typename generic_mesh<T, 1>::face::id_type>
faces_id(const generic_mesh<T, 1>& msh, const typename generic_mesh<T, 1>::cell& cl)
{
    /* This specialization is a quick hack, the real problem should be fixed in generic_element. */
    auto face_ids = cl.faces_ids();
    return std::vector<typename generic_mesh<T, 1>::face::id_type>(face_ids.begin(), face_ids.end());
}

/* Determine if a mesh element is convex. The idea is to turn CCW and for each
 * pair of consecutive edges check the cross product. If it is positive, you're
 * turning left. If you have a right turn, then the polygon is not convex. */
template<typename T>
bool
is_convex(const disk::generic_mesh<T,2>& msh,
    const typename disk::generic_mesh<T,2>::cell_type& cl)
{
    auto pts = points(msh, cl);
    assert(pts.size() > 2);
    if (pts.size() == 3)
        return true;

    for (size_t i = 0; i < pts.size(); i++)
    {
        auto p0 = pts[i];
        auto p1 = pts[(i+1)%pts.size()];
        auto p2 = pts[(i+2)%pts.size()];
        auto v0 = p1 - p0;
        auto v1 = p2 - p1;
        auto cross = v0.x()*v1.y() - v0.y()*v1.x();
        if (cross < 0)
            return false;
    }

    return true;
}


/**
  * \brief Return the volume of the specified 3D cell
  *
  * \param msh a mesh
  * \param cl a 3D cell
  * \return Return the volume of the specified 3D cell
  *
  */

template<typename T>
T
measure(const generic_mesh<T,3>& msh, const typename generic_mesh<T,3>::cell& cl)
{
    T vol = 0.0;
    auto rss = split_in_raw_tetrahedra(msh, cl);
    for (auto& rs : rss)
        vol += measure(rs);

    return vol;
}

/**
  * \brief Return the area of the specified 3D face
  *
  * \param msh a mesh
  * \param fc a 3D face
  * \return Return the area of the specified 3D face
  *
  */

template<typename T>
T
measure(const generic_mesh<T,3>& msh, const typename generic_mesh<T,3>::face& fc)
{
    auto pts = points(msh, fc);

    T acc{};
    for (size_t i = 1; i < pts.size() - 1; i++)
    {
        acc += area_triangle_kahan(pts.at(0), pts.at(i), pts.at(i + 1));
    }

    return acc;
}

/**
 * \brief Return the area of the specified 2D cell
 *
 * \param msh Reference to the mesh
 * \param cl Reference to a cell
 * \return Return the area of the specified 2D cell
 */
template<typename T>
T
measure(const generic_mesh<T,2>& msh,
    const typename generic_mesh<T,2>::cell& cl)
{
    /* Uses the divergence theorem: this way works on
     * nonconvex elements without subtriangulating. */
    T tot_meas = 0.0;
    auto fcs = faces(msh, cl);
    for (auto& fc : fcs) {
        auto bar = barycenter(msh, fc);
        auto meas = measure(msh, fc);
        auto n = normal(msh, cl, fc);
        tot_meas += meas * (bar.x()*n[0]/2. + bar.y()*n[1]/2.);
    }
    return tot_meas;
}

template<typename T>
point<T,2>
barycenter(const generic_mesh<T,2>& msh,
    const typename generic_mesh<T,2>::cell& cl)
{
    T Cx = 0.0;
    T Cy = 0.0;
    T A = 0.0;

    auto pts = points(msh, cl);
    for (size_t i = 0; i < pts.size(); i++)
    {
        auto& p0 = pts[i];
        auto& p1 = pts[(i+1)%pts.size()];
        auto d = (p0.x()*p1.y() - p1.x()*p0.y());
        Cx += (p0.x() + p1.x())*d;
        Cy += (p0.y() + p1.y())*d;
        A += 0.5*d;
    }

    Cx /= 6*A;
    Cy /= 6*A;

    return point<T,2>(Cx, Cy);
}

template<typename Mesh, typename Element>
point<typename Mesh::coordinate_type, Mesh::dimension>
barycenter_bis(const Mesh& msh, const Element& elm)
{
    auto pts = points(msh, elm);
    auto bar = std::accumulate(std::next(pts.begin()), pts.end(), pts.front());
    return bar / typename Mesh::coordinate_type( pts.size() );
}

/**
 * \brief Return the length of the specified 2D face
 *
 * \param msh a mesh
 * \param fc a 2D face
 * \return Return the length of the specified 2D face
 */

template<typename T>
T
measure(const generic_mesh<T,2>& msh, const typename generic_mesh<T,2>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 2);
    return (pts[1] - pts[0]).to_vector().norm();
}

/**
  * \brief Return the length of the specified 1D cell
  *
  * \param msh a mesh
  * \param cl a 1D cell
  * \return Return the length of the specified 1D cell
  *
  */
template<typename T>
T
measure(const generic_mesh<T,1>& msh, const typename generic_mesh<T,1>::cell& cl)
{
    auto pts = points(msh, cl);
    assert(pts.size() == 2);
    return (pts[1] - pts[0]).to_vector().norm();
}

/**
  * \brief Return the measure, i.e. 1, of the specified 1D face
  *
  * \param msh a mesh
  * \param fc a 1D face
  * \return Return 1.0
  *
  */
template<typename T>
T
measure(const generic_mesh<T,1>& msh, const typename generic_mesh<T,1>::face& fc)
{
    return T(1);
}

/**
 * @brief
 *
 * @tparam T
 * @return T
 */
template<typename T>
T
diameter(const generic_mesh<T,1>&, const typename generic_mesh<T,1>::face&)
{
    return 1.;
}


} // namespace disk

#endif
