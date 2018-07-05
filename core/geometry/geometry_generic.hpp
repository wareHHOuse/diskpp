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

#ifndef _GEOMETRY_GENERIC_HPP_
#define _GEOMETRY_GENERIC_HPP_

#include "geometry/element_generic.hpp"
#include "quadratures/raw_simplices.hpp"

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

/* Return the number of elements of the specified cell */
template<typename T, size_t DIM>
size_t
howmany_faces(const generic_mesh<T,DIM>& msh,
              const typename generic_mesh<T,DIM>::cell& cl)
{
    return cl.subelement_size();
}

/* Return the actual faces of the specified element */
template<typename T, size_t DIM>
std::vector<typename generic_mesh<T, DIM>::face>
faces(const generic_mesh<T, DIM>& msh,
      const typename generic_mesh<T, DIM>::cell& cl)
{
    auto faces_begin = msh.faces_begin();
    auto id_to_face = [&](const typename generic_mesh<T, DIM>::face::id_type& id) -> auto {
        return *std::next(faces_begin, id);
    };

    std::vector<typename generic_mesh<T, DIM>::face> ret;
    ret.resize( cl.subelement_size() );

    std::transform(cl.subelement_id_begin(), cl.subelement_id_end(),
                   ret.begin(), id_to_face);

    return ret;
}

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

template<typename T>
T
measure(const generic_mesh<T,3>& msh, const typename generic_mesh<T,3>::face& fc)
{
    auto pts = points(msh, fc);

    T acc{};
    for (size_t i = 1; i < pts.size() - 1; i++)
    {
        auto u = (pts.at(i) - pts.at(0)).to_vector();
        auto v = (pts.at(i+1) - pts.at(0)).to_vector();
        auto n = u.cross(v);
        acc += n.norm() / T(2);
    }

    return acc;
}

/* Compute the measure of a 2-cell (= area) */
template<typename T>
T
measure(const generic_mesh<T,2>& msh, const typename generic_mesh<T,2>::cell& cl)
{
    auto pts = points(msh, cl);

    T acc{};
    for (size_t i = 1; i < pts.size() - 1; i++)
    {
        auto u = (pts.at(i) - pts.at(0)).to_vector();
        auto v = (pts.at(i+1) - pts.at(0)).to_vector();
        auto n = cross(u, v);
        acc += n.norm() / T(2);
    }

    return acc;
}

/* Compute the measure of a 2-face (= length) */
template<typename T>
T
measure(const generic_mesh<T,2>& msh, const typename generic_mesh<T,2>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 2);
    return (pts[1] - pts[0]).to_vector().norm();
}

/* Compute the measure of a 1-cell (= area) */
template<typename T>
T
measure(const generic_mesh<T,1>& msh, const typename generic_mesh<T,1>::cell& cl)
{
    auto pts = points(msh, cl);
    assert(pts.size() == 2);
    return (pts[1] - pts[0]).to_vector().norm();
}

/* Compute the measure of a 1-face (= length) */
template<typename T>
T
measure(const generic_mesh<T,1>& msh, const typename generic_mesh<T,1>::face& fc)
{
    return T(1);
}

template<typename T>
T
diameter(const generic_mesh<T,1>&, const typename generic_mesh<T,1>::face&)
{
    return 1.;
}


template<typename T>
static_vector<T, 2>
normal(const generic_mesh<T,2>& msh,
       const typename generic_mesh<T,2>::cell& cl,
       const typename generic_mesh<T,2>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 2);

    auto v = pts[1] - pts[0];
    auto n = (point<T,2>({-v.y(), v.x()})).to_vector();

    auto cell_bar = barycenter(msh, cl);
    auto face_bar = barycenter(msh, fc);
    auto outward_vector = (face_bar - cell_bar).to_vector();

    if ( n.dot(outward_vector) < T(0) )
        return -n/n.norm();

    return n/n.norm();
}

template<typename T>
T
normal(const generic_mesh<T,1>& msh,
       const typename generic_mesh<T,1>::cell& cl,
       const typename generic_mesh<T,1>::face& fc)
{
    auto fcs = faces(msh, cl);
    assert(fcs.size() == 2);

    if (fc == fcs[0])
        return -1.;

    if (fc == fcs[1])
        return 1.;

    throw std::logic_error("shouldn't have arrived here");
}


} // namespace disk

#endif
