/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2024
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */
/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
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

/*
 * This source file is part of EMT, the ElectroMagneticTool.
 *
 * Copyright (C) 2013-2015, Matteo Cicuttin - matteo.cicuttin@uniud.it
 * Department of Electrical Engineering, University of Udine
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Udine nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR(s) ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(s) BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/geometry/geometry.hpp"

namespace disk {

template<typename T>
using vec3 = Eigen::Matrix<T, 3, 1>;

namespace priv {

template<typename T>
T
dot(const vec3<T>& a, const vec3<T>& b)
{
    return a.dot(b);
}

template<typename T>
vec3<T>
cross(const vec3<T>& a, const vec3<T>& b)
{
    return a.cross(b);
}

template<typename T>
T
volume_signed(const simplicial_mesh<T, 3>& mesh,
       const typename simplicial_mesh<T, 3>::volume_type& vol)
{
    auto pts = points(mesh, vol);
    assert(pts.size() == 4);

    auto v0 = (pts[1] - pts[0]).to_vector();
    auto v1 = (pts[2] - pts[0]).to_vector();
    auto v2 = (pts[3] - pts[0]).to_vector();

    return dot(v0, cross(v1, v2))/6.0;
}

template<typename T>
T
volume_unsigned(const simplicial_mesh<T, 3>& mesh,
       const typename simplicial_mesh<T, 3>::volume_type& vol)
{
    auto pts = points(mesh, vol);
    assert(pts.size() == 4);

    auto v0 = (pts[1] - pts[0]).to_vector();
    auto v1 = (pts[2] - pts[0]).to_vector();
    auto v2 = (pts[3] - pts[0]).to_vector();

    return std::abs(dot(v0, cross(v1, v2))/6.0);
}

} //namespace priv

/***************************************************************************/
/* Get sub-elements */
template<typename T>
std::array<typename simplicial_mesh<T, 3>::edge_type, 3>
edges(const simplicial_mesh<T, 3>& msh,
      const typename simplicial_mesh<T, 3>::surface_type& surf)
{
    typedef typename simplicial_mesh<T, 3>::edge_type edge_type;
    std::array<edge_type, 3> ret;

    auto ptids = surf.point_ids();
    assert(ptids.size() == 3);

    ret[0] = edge_type({ptids[0], ptids[1]});
    ret[1] = edge_type({ptids[0], ptids[2]});
    ret[2] = edge_type({ptids[1], ptids[2]});

    return ret;
}

template<typename T>
std::array<typename simplicial_mesh<T, 3>::edge_type, 6>
edges(const simplicial_mesh<T, 3>& msh,
      const typename simplicial_mesh<T, 3>::volume_type& vol)
{ //ok
    typedef typename simplicial_mesh<T, 3>::edge_type edge_type;
    std::array<edge_type, 6> ret;

    auto ptids = vol.point_ids();
    assert(ptids.size() == 4);

    ret[0] = edge_type({ptids[0], ptids[1]});
    ret[1] = edge_type({ptids[0], ptids[2]});
    ret[2] = edge_type({ptids[0], ptids[3]});
    ret[3] = edge_type({ptids[1], ptids[2]});
    ret[4] = edge_type({ptids[1], ptids[3]});
    ret[5] = edge_type({ptids[2], ptids[3]});

    return ret;
}


template<disk::mesh_3D Mesh>
struct edge_iterproxy {
    const Mesh& mesh;
    edge_iterproxy(const Mesh& msh)
        : mesh(msh)
    {}

    auto begin() { return mesh.backend_storage()->edges.begin(); }
    auto end() { return mesh.backend_storage()->edges.end(); }
    auto begin() const { return mesh.backend_storage()->edges.begin(); }
    auto end() const { return mesh.backend_storage()->edges.end(); }
    auto size() const { return mesh.backend_storage()->edges.size(); }
};

template<disk::mesh_3D Mesh>
auto edges(const Mesh& msh) {
    return edge_iterproxy(msh);
}

template<disk::mesh_3D Mesh>
auto offset(const Mesh& msh, const typename Mesh::edge_type& edg)
{
    auto ebegin = msh.backend_storage()->edges.begin();
    auto eend = msh.backend_storage()->edges.end();
    auto itor = std::lower_bound(ebegin, eend, edg);
    if (itor == eend || *itor != edg)
        throw std::invalid_argument("offset(): edge not found");

    return std::distance(ebegin, itor);
}

template<typename T>
std::array<size_t, 3>
edge_ids(const simplicial_mesh<T, 3>& msh,
         const typename simplicial_mesh<T, 3>::surface_type& surf)
{ //ok
    auto edgs = edges(msh, surf);

    std::array<size_t, 3> ret;
    for (size_t i = 0; i < 3; i++)
        ret[i] = offset(msh, edgs[i]);

    return ret;
}

template<typename T>
std::array<size_t, 6>
edge_ids(const simplicial_mesh<T, 3>& msh,
         const typename simplicial_mesh<T, 3>::volume_type& vol)
{ //ok
    auto edgs = edges(msh, vol);

    std::array<size_t, 6> ret;
    for (size_t i = 0; i < 6; i++)
        ret[i] = offset(msh, edgs[i]);

    return ret;
}

/***************************************************************************/
/* Calculate primal edge vectors */
template<typename T>
vec3<T>
primal_edge_vectors(const simplicial_mesh<T, 3>& mesh,
                    const typename simplicial_mesh<T, 3>::edge_type& edge)
{
    auto pts = points(mesh, edge);
    assert(pts.size() == 2);

    return (pts[1] - pts[0]).to_vector();
}

template<typename T>
std::array<vec3<T>, 3>
primal_edge_vectors(const simplicial_mesh<T, 3>& mesh,
                    const typename simplicial_mesh<T, 3>::surface_type& surf)
{
    std::array<vec3<T>, 3> ret;

    auto pts = points(mesh, surf);
    assert(pts.size() == 3);

    ret[0] = (pts[1] - pts[0]).to_vector();
    ret[1] = (pts[2] - pts[0]).to_vector();
    ret[2] = (pts[2] - pts[1]).to_vector();

    return ret;
}

template<typename T>
std::array<vec3<T>, 6>
primal_edge_vectors(const simplicial_mesh<T, 3>& mesh,
                    const typename simplicial_mesh<T, 3>::volume_type& vol)
{ //ok
    std::array<vec3<T>, 6> ret;

    auto pts = points(mesh, vol);
    assert(pts.size() == 4);

    ret[0] = (pts[1] - pts[0]).to_vector();
    ret[1] = (pts[2] - pts[0]).to_vector();
    ret[2] = (pts[3] - pts[0]).to_vector();
    ret[3] = (pts[2] - pts[1]).to_vector();
    ret[4] = (pts[3] - pts[1]).to_vector();
    ret[5] = (pts[3] - pts[2]).to_vector();

    return ret;
}

/***************************************************************************/
/* Calculate primal area vectors */
template<typename T>
std::array<vec3<T>, 4>
primal_area_vectors(const simplicial_mesh<T, 3>& mesh,
                    const typename simplicial_mesh<T, 3>::volume_type& vol)
{ //ok
    typedef typename simplicial_mesh<T, 3>::edge_type edge_type;

    std::array<edge_type, 6> edgs = edges(mesh, vol);
    std::array<vec3<T>, 6> evecs;

    for (size_t i = 0; i < 6; i++)
        evecs[i] = primal_edge_vectors(mesh, edgs[i]);

    auto pts = points(mesh, vol);
    assert(pts.size() == 4);

    std::array<vec3<T>, 4> ret;

    ret[0] = priv::cross( evecs[3], evecs[4] ) / 2.0;
    ret[1] = priv::cross( evecs[1], evecs[2] ) / 2.0;
    ret[2] = priv::cross( evecs[0], evecs[2] ) / 2.0;
    ret[3] = priv::cross( evecs[0], evecs[1] ) / 2.0;

    return ret;
}

/***************************************************************************/
/* Calculate dual edge vectors */
template<typename T>
std::array<vec3<T>, 3>
dual_edge_vectors(const simplicial_mesh<T, 3>& mesh,
                  const typename simplicial_mesh<T, 3>::surface_type& t,
                  const typename simplicial_mesh<T, 3>::volume_type& parent_vol)
{
    /* the returned vectors are such that cross(pev, dev) = outward */
    std::array<vec3<T>, 3> ret;

    bool negative = priv::volume_signed(mesh, parent_vol) < 0;
    auto fcs = faces(mesh, parent_vol);

    size_t face = 42;

    if (fcs[0] == t) face = 0;
    if (fcs[1] == t) face = 1;
    if (fcs[2] == t) face = 2;
    if (fcs[3] == t) face = 3;

    if (face == 42)
        throw std::invalid_argument("The specified volume does not "
                                    "have the specified face");

    bool P, Q;
    P = !negative && ( face == 0 || face == 2 );
    Q =  negative && ( face == 1 || face == 3 );

    T sign = (P || Q) ? 1 : -1;

    auto face_bar = barycenter(mesh, t);
    auto edgs = edges(mesh, t);

    ret[0] =  ( face_bar - barycenter(mesh, edgs[0]) ).to_vector() * sign;
    ret[1] =  ( face_bar - barycenter(mesh, edgs[1]) ).to_vector() * sign;
    ret[2] = -( face_bar - barycenter(mesh, edgs[2]) ).to_vector() * sign;

    return ret;
}

#if 0
template<typename CoordT, typename IdxT>
std::array<vec3<CoordT>, 3>
dual_edge_vectors(const tetrahedral_mesh<CoordT, IdxT>& mesh,
                  const triangle<CoordT, IdxT>& t)
{
    /* the returned vectors are such that cross(pev, dev) = outward */
    std::array<vec3<CoordT>, 3> ret;

    auto vols = mesh.associated_volumes(t);
    if (vols.size() != 1)
        throw std::invalid_argument("You asked to calculate dual edges on "
                                    "an internal triangle, but without "
                                    "specifying the parent volume.");

    auto parent_vol = *vols.begin();

    return dual_edge_vectors(mesh, t, parent_vol);
}
#endif

template<typename T>
std::array<vec3<T>, 4>
dual_edge_vectors(const simplicial_mesh<T, 3>& mesh,
                  const typename simplicial_mesh<T, 3>::volume_type& vol)
{
    typedef typename simplicial_mesh<T, 3>::surface_type surface_type;
    std::array<vec3<T>, 4> ret;

    std::array<surface_type, 4> fcs = faces(mesh, vol);
    std::array<point<T,3>, 4> face_barycenters;
    for (size_t i = 0; i < 4; i++)
        face_barycenters[i] = barycenter(mesh, fcs[i]);

    auto vol_barycenter = barycenter(mesh, vol);

    T sign = (priv::volume_signed(mesh, vol) < 0) ? -1 : 1;

    ret[0] = (face_barycenters[0] - vol_barycenter).to_vector()*sign;
    ret[1] = (vol_barycenter - face_barycenters[1]).to_vector()*sign;
    ret[2] = (face_barycenters[2] - vol_barycenter).to_vector()*sign;
    ret[3] = (vol_barycenter - face_barycenters[3]).to_vector()*sign;

    return ret;
}

/***************************************************************************/
/* Calculate dual area vectors */
template<typename T>
std::array<vec3<T>, 6>
dual_area_vectors(const simplicial_mesh<T, 3>& mesh,
                  const typename simplicial_mesh<T, 3>::volume_type& vol)
{ //ok
    std::array<vec3<T>, 6> ret;

    auto edgs = edges(mesh, vol);
    std::array<point<T, 3>, 6> ebs;
    for (size_t i = 0; i < 6; i++)
        ebs[i] = barycenter(mesh, edgs[i]);

    auto fcs = faces(mesh, vol);
    std::array<point<T, 3>, 4> fbs;
    for (size_t i = 0; i < 4; i++)
        fbs[i] = barycenter(mesh, fcs[i]);

    point<T, 3> vb = barycenter(mesh, vol);

    /* Area of quadrilateral ABCD is 0.5*|AC x BD|. Orientation of the dual
     * area vector must be the same of the primal edge vector i.e. it must
     * satisfy 'dot(PEV, DAV) >= 0' */

    T sign = (priv::volume_signed(mesh, vol) < 0) ? -0.5 : 0.5;

    ret[0] = priv::cross( (vb-ebs[0]).to_vector(), (fbs[2]-fbs[3]).to_vector() )*sign;
    ret[1] = priv::cross( (vb-ebs[1]).to_vector(), (fbs[3]-fbs[1]).to_vector() )*sign;
    ret[2] = priv::cross( (vb-ebs[2]).to_vector(), (fbs[1]-fbs[2]).to_vector() )*sign;
    ret[3] = priv::cross( (vb-ebs[3]).to_vector(), (fbs[0]-fbs[3]).to_vector() )*sign;
    ret[4] = priv::cross( (vb-ebs[4]).to_vector(), (fbs[2]-fbs[0]).to_vector() )*sign;
    ret[5] = priv::cross( (vb-ebs[5]).to_vector(), (fbs[0]-fbs[1]).to_vector() )*sign;

    return ret;
}

/***************************************************************************/
/* Integrate on primal edges */
template<typename T>
Eigen::Matrix<T, 3, 1>
integrate_on_primal_edges(const simplicial_mesh<T, 3>& mesh,
                          const typename simplicial_mesh<T, 3>::surface_type& surf,
                          const vec3<T> vecf)
{
    Eigen::Matrix<T, 3, 1> ret;

    auto pev = primal_edge_vectors(mesh, surf);

    for (size_t i = 0; i < 3; i++)
        ret(i) = priv::dot(vecf, pev[i]);

    return ret;
}

template<typename T>
Eigen::Matrix<T, 6, 1>
integrate_on_primal_edges(const simplicial_mesh<T, 3>& mesh,
                          const typename simplicial_mesh<T, 3>::volume_type& vol,
                          const vec3<T> vecf)
{
    Eigen::Matrix<T, 6, 1> ret;

    auto pev = primal_edge_vectors(mesh, vol);

    for (size_t i = 0; i < 6; i++)
        ret(i) = priv::dot(vecf, pev[i]);

    return ret;
}

/***************************************************************************/
/* Integrate on primal faces */
template<typename T>
Eigen::Matrix<T, 3, 1>
integrate_on_primal_faces(const simplicial_mesh<T, 3>& mesh,
                          const typename simplicial_mesh<T, 3>::surface_type& surf,
                          const vec3<T> vecf)
{
    Eigen::Matrix<T, 3, 1> ret;

    auto pav = primal_area_vectors(mesh, surf);

    for (size_t i = 0; i < 3; i++)
        ret(i) = priv::dot(vecf, pav[i]);

    return ret;
}

template<typename T>
Eigen::Matrix<T, 4, 1>
integrate_on_primal_faces(const simplicial_mesh<T, 3>& mesh,
                          const typename simplicial_mesh<T, 3>::volume_type& vol,
                          const vec3<T> vecf)
{
    Eigen::Matrix<T, 3, 1> ret;

    auto pav = primal_area_vectors(mesh, vol);

    for (size_t i = 0; i < 4; i++)
        ret(i) = priv::dot(vecf, pav[i]);

    return ret;
}

/***************************************************************************/
/* Integrate on dual edges */
template<typename T>
Eigen::Matrix<T, 3, 1>
integrate_on_dual_edges(const simplicial_mesh<T, 3>& mesh,
                        const typename simplicial_mesh<T, 3>::surface_type& surf,
                        const vec3<T> vecf)
{
    Eigen::Matrix<T, 3, 1> ret;

    auto dev = dual_edge_vectors(mesh, surf);

    for (size_t i = 0; i < 3; i++)
        ret(i) = priv::dot(vecf, dev[i]);

    return ret;
}

template<typename T>
Eigen::Matrix<T, 4, 1>
integrate_on_dual_edges(const simplicial_mesh<T, 3>& mesh,
                        const typename simplicial_mesh<T, 3>::volume_type& vol,
                        const vec3<T> vecf)
{
    Eigen::Matrix<T, 4, 1> ret;

    auto dev = dual_edge_vectors(mesh, vol);

    for (size_t i = 0; i < 4; i++)
        ret(i) = priv::dot(vecf, dev[i]);

    return ret;
}

/***************************************************************************/
/* Integrate on dual faces */
template<typename T>
Eigen::Matrix<T, 3, 1>
integrate_on_dual_faces(const simplicial_mesh<T, 3>& mesh,
                        const typename simplicial_mesh<T, 3>::surface_type& surf,
                        const vec3<T> vecf)
{
    Eigen::Matrix<T, 3, 1> ret;

    auto dav = dual_area_vectors(mesh, surf);

    for (size_t i = 0; i < 3; i++)
        ret(i) = priv::dot(vecf, dav[i]);

    return ret;
}

template<typename T>
Eigen::Matrix<T, 6, 1>
integrate_on_dual_faces(const simplicial_mesh<T, 3>& mesh,
                        const typename simplicial_mesh<T, 3>::volume_type& vol,
                        const vec3<T> vecf)
{
    Eigen::Matrix<T, 6, 1> ret;

    auto dav = dual_area_vectors(mesh, vol);

    for (size_t i = 0; i < 6; i++)
        ret(i) = priv::dot(vecf, dav[i]);

    return ret;
}

/***************************************************************************/
/* Find normals */
namespace priv {

template<typename T>
std::array<vec3<T>, 3>
normals(const simplicial_mesh<T, 3>& mesh,
        const typename simplicial_mesh<T, 3>::surface_type& t)
{
    std::array<vec3<T>, 3> ret;

    auto pev = primal_edge_vectors(mesh, t);

    auto normal = priv::cross(pev[0], pev[2]);

    auto n0 =  priv::cross(pev[0], normal);
    auto n1 =  priv::cross(pev[1], normal);
    auto n2 = -priv::cross(pev[2], normal);

    ret[0] = n0/n0.norm();
    ret[1] = n1/n1.norm();
    ret[2] = n2/n2.norm();

    return ret;
}

template<typename T>
std::array<vec3<T>, 4>
normals(const simplicial_mesh<T, 3>& mesh,
        const typename simplicial_mesh<T, 3>::volume_type& t)
{
    std::array<vec3<T>, 4> ret;

    auto pav = primal_area_vectors(mesh, t);

    if ( priv::volume_signed(mesh, t) > 0 )
    {
        ret[0] =  pav[0]/pav[0].norm();
        ret[1] = -pav[1]/pav[1].norm();
        ret[2] =  pav[2]/pav[2].norm();
        ret[3] = -pav[3]/pav[3].norm();
    }
    else
    {
        ret[0] = -pav[0]/pav[0].norm();
        ret[1] =  pav[1]/pav[1].norm();
        ret[2] = -pav[2]/pav[2].norm();
        ret[3] =  pav[3]/pav[3].norm();
    }

    return ret;
}

} //namespace priv
/***************************************************************************/
/* Gradient matrix for tetrahedral elements */

/*
template<typename T, typename CoordT, typename IdxT>
typename arma::Mat<T>::template fixed<3,3>
tri_grad_matrix(const tetrahedral_mesh<CoordT, IdxT>&,
                const triangle<CoordT, IdxT>&)
{
    typename arma::Mat<T>::template fixed<3,3> G = { -1,  0, -1,
                                                      1, -1,  0,
                                                      0,  1,  1 };

    return G;
}
*/

template<typename T>
Eigen::Matrix<T, 6, 4>
grad_matrix(const simplicial_mesh<T, 3>&,
            const typename simplicial_mesh<T, 3>::volume_type&)
{
    Eigen::Matrix<T, 6, 4> G;
    G << -1,  1,  0,  0,
         -1,  0,  1,  0,
         -1,  0,  0,  1,
          0, -1,  1,  0,
          0, -1,  0,  1,
          0,  0, -1,  1;
    return G;
}

template<typename T>
Eigen::Matrix<T, 4, 6>
curl_matrix(const simplicial_mesh<T, 3>&,
            const typename simplicial_mesh<T, 3>::volume_type&)
{
    Eigen::Matrix<T, 4, 6> C;
    C <<  0,  0,  0, +1, -1, +1,
          0, +1, -1,  0,  0, +1,
         +1,  0, -1,  0, +1,  0,
         +1, -1,  0, +1,  0,  0;
    return C;
}

/***************************************************************************/
/* Calculate edge matrix - scalar parameter */

namespace priv {

template<typename T>
Eigen::Matrix<T, 6, 6>
edge_matrix_simplicial_element(const std::array<vec3<T>, 4>& pav,
    const std::array<vec3<T>, 4>& Ppav)
{
    Eigen::Matrix<T, 6, 6> ret;

    /* Eigen stores elements in column-major order, fill by column! */
    ret(0,0) =   dot(pav[0], Ppav[0]) + dot(pav[1], Ppav[1]);
    ret(1,0) = - dot(pav[1], Ppav[2]);
    ret(2,0) = + dot(pav[1], Ppav[3]);
    ret(3,0) = - dot(pav[0], Ppav[2]);
    ret(4,0) = + dot(pav[0], Ppav[3]);
    ret(5,0) =   0;

    ret(0,1) =   ret(1,0);
    ret(1,1) =   dot(pav[0], Ppav[0]) + dot(pav[2], Ppav[2]);
    ret(2,1) = - dot(pav[2], Ppav[3]);
    ret(3,1) = - dot(pav[0], Ppav[1]);
    ret(4,1) =   0;
    ret(5,1) = + dot(pav[0], Ppav[3]);

    ret(0,2) =   ret(2,0);
    ret(1,2) =   ret(2,1);
    ret(2,2) =   dot(pav[0], Ppav[0]) + dot(pav[3], Ppav[3]);
    ret(3,2) =   0;
    ret(4,2) = - dot(pav[0], Ppav[1]);
    ret(5,2) =   dot(pav[0], Ppav[2]);

    ret(0,3) =   ret(3,0);
    ret(1,3) =   ret(3,1);
    ret(2,3) =   0;
    ret(3,3) =   dot(pav[1], Ppav[1]) + dot(pav[2], Ppav[2]);
    ret(4,3) = - dot(pav[2], Ppav[3]);
    ret(5,3) = - dot(pav[1], Ppav[3]);

    ret(0,4) =   ret(4,0);
    ret(1,4) =   0;
    ret(2,4) =   ret(4,2);
    ret(3,4) =   ret(4,3);
    ret(4,4) =   dot(pav[1], Ppav[1]) + dot(pav[3], Ppav[3]);
    ret(5,4) = - dot(pav[1], Ppav[2]);

    ret(0,5) =   0;
    ret(1,5) =   ret(5,1);
    ret(2,5) =   ret(5,2);
    ret(3,5) =   ret(5,3);
    ret(4,5) =   ret(5,4);
    ret(5,5) =   dot(pav[2], Ppav[2]) + dot(pav[3], Ppav[3]);

    return ret;
}

template<typename T>
Eigen::Matrix<T, 4, 4>
face_matrix_simplicial_element(const std::array<vec3<T>, 6>& pev,
    const std::array<vec3<T>, 6>& Ppev)
{
    Eigen::Matrix<T, 4, 4> ret;

    /* Eigen stores elements in column-major order, fill by column! */
    ret(0,0) =   dot(pev[0], Ppev[0]) + dot(pev[1], Ppev[1]) + dot(pev[2], Ppev[2]);
    ret(1,0) = - dot(pev[3], Ppev[1]) - dot(pev[4], Ppev[2]);
    ret(2,0) = - dot(pev[3], Ppev[0]) + dot(pev[5], Ppev[2]);
    ret(3,0) = + dot(pev[4], Ppev[0]) + dot(pev[5], Ppev[1]);

    ret(0,1) =   ret(1,0);
    ret(1,1) =   dot(pev[0], Ppev[0]) + dot(pev[3], Ppev[3]) + dot(pev[4], Ppev[4]);
    ret(2,1) = - dot(pev[1], Ppev[0]) - dot(pev[5], Ppev[4]);
    ret(3,1) = + dot(pev[2], Ppev[0]) - dot(pev[3], Ppev[5]);

    ret(0,2) =   ret(2,0);
    ret(1,2) =   ret(2,1);
    ret(2,2) =   dot(pev[1], Ppev[1]) + dot(pev[3], Ppev[3]) + dot(pev[5], Ppev[5]);
    ret(3,2) = - dot(pev[2], Ppev[1]) - dot(pev[4], Ppev[3]);

    ret(0,3) =   ret(3,0);
    ret(1,3) =   ret(3,1);
    ret(2,3) =   ret(3,2);
    ret(3,3) =   dot(pev[2], Ppev[2]) + dot(pev[4], Ppev[4]) + dot(pev[5], Ppev[5]);

    return ret;
}

}

template<typename T>
Eigen::Matrix<T, 6, 6>
edge_matrix(const simplicial_mesh<T, 3>& mesh,
            const typename simplicial_mesh<T, 3>::volume_type& vol,
            const T& scalar_parameter)
{
    auto mult = scalar_parameter / ( 36 * priv::volume_unsigned(mesh, vol) );

    auto pav = primal_area_vectors(mesh, vol);
    auto Ppav = pav;
    for (auto& p : Ppav)
        p *= mult;

    return priv::edge_matrix_simplicial_element(pav, Ppav);
}


template<typename T>
Eigen::Matrix<T, 6, 6>
edge_matrix(const simplicial_mesh<T, 3>& mesh,
            const typename simplicial_mesh<T, 3>::volume_type& vol,
            const static_matrix<T, 3, 3>& tensor_parameter)
{
    auto mult = tensor_parameter / ( 36 * priv::volume_unsigned(mesh, vol) );

    auto pav = primal_area_vectors(mesh, vol);
    auto Ppav = pav;
    for (auto& p : Ppav)
        p = mult*p;

    return priv::edge_matrix_simplicial_element(pav, Ppav);
}

template<typename T>
Eigen::Matrix<T, 4, 4>
face_matrix(const simplicial_mesh<T, 3>& mesh,
            const typename simplicial_mesh<T, 3>::volume_type& vol,
            const T& scalar_parameter)
{
    auto mult = scalar_parameter / ( 36 * priv::volume_unsigned(mesh, vol) );

    auto pev = primal_edge_vectors(mesh, vol);
    auto Ppev = pev;
    for (auto& p : Ppev)
        p *= mult;

    return priv::face_matrix_simplicial_element(pev, Ppev);
}


template<typename T>
Eigen::Matrix<T, 4, 4>
face_matrix(const simplicial_mesh<T, 3>& mesh,
            const typename simplicial_mesh<T, 3>::volume_type& vol,
            const static_matrix<T, 3, 3>& tensor_parameter)
{
    auto mult = tensor_parameter / ( 36 * priv::volume_unsigned(mesh, vol) );

    auto pev = primal_edge_vectors(mesh, vol);
    auto Ppev = pev;
    for (auto& p : Ppev)
        p = mult*p;

    return priv::face_matrix_simplicial_element(pev, Ppev);
}


} // namespace disk
