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

#pragma once

#include "mesh/mesh.hpp"
#include "geometry/geometry_all.hpp"
#include "geometry/element_simplicial.hpp"

namespace disk {

namespace simplicial_priv {
template<size_t DIM>
struct element_types
{
    static_assert(DIM > 0 && DIM <= 3, "element_types: CODIM must be less than DIM");
};

template<>
struct element_types<3> {
        typedef simplicial_element<3,0>    volume_type;
        typedef simplicial_element<3,1>    surface_type;
        typedef simplicial_element<3,2>    edge_type;
        typedef simplicial_element<3,3>    node_type;
};

template<>
struct element_types<2> {
        typedef simplicial_element<2,0>    surface_type;
        typedef simplicial_element<2,1>    edge_type;
        typedef simplicial_element<2,2>    node_type;
};

template<>
struct element_types<1> {
        typedef simplicial_element<1,0>    edge_type;
        typedef simplicial_element<1,1>    node_type;
};


} // namespace priv

template<typename T, size_t DIM>
using simplicial_mesh_storage = mesh_storage<T, DIM, simplicial_priv::element_types<DIM>>;

template<typename T, size_t DIM>
using simplicial_mesh = mesh<T, DIM, simplicial_mesh_storage<T, DIM>>;

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

template<typename T, size_t DIM>
std::array<simplicial_element<3,1>, 4>
faces(const simplicial_mesh<T, DIM>&, const simplicial_element<3,0>& vol)
{
    std::array<simplicial_element<3,1>, 4> ret;

    auto ptids = vol.point_ids();
    assert(ptids.size() == 4);

    ret[0] = simplicial_element<3,1>( { ptids[1], ptids[2], ptids[3] } );
    ret[1] = simplicial_element<3,1>( { ptids[0], ptids[2], ptids[3] } );
    ret[2] = simplicial_element<3,1>( { ptids[0], ptids[1], ptids[3] } );
    ret[3] = simplicial_element<3,1>( { ptids[0], ptids[1], ptids[2] } );

    return ret;
}

template<typename T, size_t DIM>
T
measure(const simplicial_mesh<T, DIM>& msh, const simplicial_element<3,0>& vol, bool signed_volume = false)
{
    auto pts = points(msh, vol);
    assert(pts.size() == 4);

    auto v0 = (pts[1] - pts[0]).to_vector();
    auto v1 = (pts[2] - pts[0]).to_vector();
    auto v2 = (pts[3] - pts[0]).to_vector();

    if (signed_volume)
        return v0.dot(v1.cross(v2))/T(6);

    return std::abs( v0.dot(v1.cross(v2))/T(6) );
}

template<typename T, size_t DIM>
T
measure(const simplicial_mesh<T, DIM>& msh, const simplicial_element<3,1>& surf)
{
    auto pts = points(msh, surf);
    assert(pts.size() == 3);

    auto v0 = (pts[1] - pts[0]).to_vector();
    auto v1 = (pts[2] - pts[0]).to_vector();

    return v0.cross(v1).norm()/T(2);
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

template<typename T, size_t DIM, typename Elem>
T
diameter(const simplicial_mesh<T, DIM>& msh, const Elem& vol)
{
    auto pts = points(msh, vol);

    T diam = 0.;

    for (size_t i = 0; i < pts.size(); i++)
        for (size_t j = i; j < pts.size(); j++)
            diam = std::max((pts[i] - pts[j]).to_vector().norm(), diam);

    return diam;
}
/*
template<typename T, size_t DIM>
T
diameter(const simplicial_mesh<T, DIM>& msh, const typename simplicial_mesh<T, DIM>::cell& cl)
{
    T c_meas = measure(msh, cl);
    T af_meas = T(0);
    auto fcs = faces(msh, cl);
    for (auto& f : fcs)
        af_meas += measure(msh, f);

    return c_meas/af_meas;
}
*/

template<typename Mesh>
class submesher;

template<typename T>
class submesher<simplicial_mesh<T,2>>
{

};

template<typename T>
void
submesh(const simplicial_mesh<T,2>& msh,
        const typename simplicial_mesh<T,2>::cell& cl,
        const T& max_h)
{
    //simplicial_mesh<T,2> submesh;
    //auto storage = submesh.backend_storage();

    auto pts = points(msh, cl);

    auto b0 = (pts[0] + pts[1]) / T(2);
    auto b1 = (pts[1] + pts[2]) / T(2);
    auto b2 = (pts[0] + pts[2]) / T(2);
}




} // namespace disk
