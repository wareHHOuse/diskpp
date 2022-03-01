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

#ifndef _GEOMETRY_ALL_HPP_
#define _GEOMETRY_ALL_HPP_

/*
 * Here we have all the queries that don't depend on the particular mesh
 * type. For example 'barycenter()' is computed in the same way for all the
 * elements, so it goes here. It can happen that there are queries that are
 * applicable to all kinds of elements, but for only some classes of elements
 * more efficient ways to do the computation exist. Put the general one here
 * and put the specialization in the right geometry_<whatever>.hpp file.
 */

#include <algorithm>
#include <vector>

#include "common/util.h"
#include "common/eigen.hpp"

namespace disk {

/**
  * \brief Compute an estimate of the mesh discretization step 'h'
  *
  * \param msh a reference to the mesh
  * \return an estimate of the mesh discretization step 'h'
  *
  */

template<typename Mesh>
typename Mesh::coordinate_type
average_diameter(const Mesh& msh)
{
    typename Mesh::coordinate_type h{};
    for (auto& cl : msh)
    {
        h += diameter(msh, cl);
    }

    return h/msh.cells_size();
}

/**
  * \brief return the list of points of an element
  *
  * \param msh a reference to the mesh
  * \param elem a generic element (cell or face)
  * \return return a vector with the points of elem
  *
  */

template<typename Mesh, typename Element>
std::vector<typename Mesh::point_type>
points(const Mesh& msh, const Element& elem)
{
    auto ptids = elem.point_ids();

    auto points_begin = msh.points_begin();
    auto ptid_to_point = [&](const point_identifier<Mesh::dimension>& pi) -> auto {
        return *std::next(points_begin, pi);
    };

    std::vector<typename Mesh::point_type> pts(ptids.size());
    std::transform(ptids.begin(), ptids.end(), pts.begin(), ptid_to_point);

    return pts;
}

/**
  * \brief Compute the barycenter of a generic element
  *
  * \param msh a reference to the mesh
  * \param elem a generic element (cell or face)
  * \return return the barycenter of elem
  *
  */

template<typename Mesh, typename Element>
point<typename Mesh::coordinate_type, Mesh::dimension>
barycenter(const Mesh& msh, const Element& elem)
{
    auto pts = points(msh, elem);
    auto bar = std::accumulate(std::next(pts.begin()), pts.end(), pts.front());
    return bar / typename Mesh::coordinate_type( pts.size() );
}

// template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
// point<T,2>
// barycenter(const Mesh<T,2,Storage>& msh, const typename Mesh<T,2,Storage>::cell_type& cl)
// {
//     std::cout << "bar2" << std::endl;

//     using std::abs;
//     T          tot_meas{};
//     point<T,2> tot_bar{};
//     auto pts = points(msh, cl);
//     for (size_t i = 1; i < pts.size()-1; i++)
//     {
//         auto d0 = pts[i] - pts[0];
//         auto d1 = pts[i+1] - pts[0];
//         auto meas = abs(d0.x()*d1.y() - d1.x()*d0.y());
//         tot_bar = tot_bar + meas*(pts[0]+pts[i]+pts[i+1]);
//         tot_meas += meas;
//     }

//     return tot_bar/(tot_meas*T(3));
// }

/**
  * \brief Compute the diameter of a generic element, i.e, the maximum distance
  * between two different points of the element.
  *
  * \param msh a reference to the mesh
  * \param elem a generic element (cell or face)
  * \return return the diameter of elem
  *
  */

template<typename Mesh, typename Element>
typename Mesh::coordinate_type
diameter(const Mesh& msh, const Element& elem)
{
    const auto pts = points(msh, elem);

    typename Mesh::coordinate_type diam = 0.;

    for (size_t i = 0; i < pts.size(); i++)
        for (size_t j = i+1; j < pts.size(); j++)
            diam = std::max((pts[i] - pts[j]).to_vector().norm(), diam);

    return diam;
}

/**
  * \brief Compute the diameter of the bounding box af an element, i.e, the maximum distance
  * between two different points of the bounding box.
  *
  * \param msh a reference to the mesh
  * \param elem an element
  * \return compute the diameter of the bounding box af an element
  *
  */
template<typename Mesh, typename Elem>
std::array<typename Mesh::coordinate_type, Mesh::dimension>
diameter_boundingbox(const Mesh&                                                                      msh,
                     const Elem&                                                                      elem,
                     static_matrix<typename Mesh::coordinate_type, Mesh::dimension, Mesh::dimension> axes =
                       static_matrix<typename Mesh::coordinate_type, Mesh::dimension, Mesh::dimension>::Identity())
{
    using T        = typename Mesh::coordinate_type;
    const auto pts = points(msh, elem);

    std::array<T, Mesh::dimension> box_min;
    std::array<T, Mesh::dimension> box_max;
    for (size_t i = 0; i < Mesh::dimension; i++)
    {
        box_min[i] = std::numeric_limits<T>::max();
        box_max[i] = -std::numeric_limits<T>::max();
    }

    for (auto& pt : pts)
    {
        const static_vector<T, Mesh::dimension> bp = (axes.transpose()) * (pt.to_vector());

        // std::cout << bp << std::endl;

        for (size_t i = 0; i < Mesh::dimension; i++)
        {
            if (bp(i) < box_min[i])
            {
                box_min[i] = bp(i);
            }

            if (bp(i) > box_max[i])
            {
                box_max[i] = bp(i);
            }

            // std::cout << box_min[i] << " , " << box_max[i] << std::endl;
        }
    }

    std::array<T, Mesh::dimension> box_h;
    for (size_t i = 0; i < Mesh::dimension; i++)
    {
        box_h[i] = abs(box_max[i] - box_min[i]);
    }

    return box_h;
}



/**
  * \brief Allows to known if the given point is inside a 2D cell
  *
  * \param msh a reference to the mesh
  * \param cl a 3D cell
  * \param pt coordinate of a point
  *
  */

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
bool
is_inside(const Mesh<T, 2, Storage>&                      msh,
          const typename Mesh<T, 2, Storage>::cell&       cl,
          const typename Mesh<T, 2, Storage>::point_type& pt)
{
    /* Nodes MUST be ordered COUNTERCLOCKWISE and the polygon must be CONVEX */
    auto pts = points(msh, cl);

    for (size_t i = 0; i < pts.size(); i++)
    {
        auto p0 = pts[i];
        auto p1 = pts[i % pts.size()];

        auto x  = pt.x();
        auto y  = pt.y();
        auto x0 = p0.x();
        auto y0 = p0.y();
        auto x1 = p1.x();
        auto y1 = p1.y();

        if ((y - y0) * (x1 - x0) - (x - x0) * (y1 - y0) < 0.0)
            return false;
    }

    return true;
}


/**
  * \brief Allows to known if a cell has at least one face on the boundary
  *
  * \param msh a reference to the mesh
  * \param cl a  cell
  * \return return true if the cell has faces on the boundary
  *
  */

template<typename Mesh>
bool
has_faces_on_boundary(const Mesh& msh, const typename Mesh::cell& cl)
{
    auto fcs = faces(msh, cl);
    for (auto& fc : fcs)
        if ( msh.is_boundary(fc) )
            return true;

    return false;
}

/**
 * @brief Compute the outward unit normal to a 2D face
 *
 * @tparam Mesh type of mesh
 * @tparam T scalar type
 * @tparam Storage type of storage for the mesh
 * @param msh mesh
 * @param cl cell
 * @param fc face
 * @return static_vector<T, 2> Compute the outward unit normal to a 2D face
 */
template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
static_vector<T, 2>
normal(const Mesh<T,2,Storage>& msh,
       const typename Mesh<T,2,Storage>::cell& cl,
       const typename Mesh<T,2,Storage>::face& fc)
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

/**
 * @brief Compute a normal of a 2D face. No guarantee on the direction.
 *
 * @tparam Mesh type of mesh
 * @tparam T scalar type
 * @tparam Storage type of storage for the mesh
 * @param msh mesh
 * @param cl cell
 * @param fc face
 * @return static_vector<T, 2> A normal to the 2D face
 */
template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
static_vector<T, 2>
normal(const Mesh<T,2,Storage>& msh,
       const typename Mesh<T,2,Storage>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 2);

    auto v = pts[1] - pts[0];
    auto n = (point<T,2>({-v.y(), v.x()})).to_vector();

    return n/n.norm();
}

/**
 * @brief Compute the outward unit normal to a 3D face
 *
 * @tparam Mesh type of mesh
 * @tparam T scalar type
 * @tparam Storage type of storage for the mesh
 * @param msh mesh
 * @param cl cell
 * @param fc face
 * @return static_vector<T, 3> Compute the outward unit normal to a 3D face
 */
template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
static_vector<T, 3>
normal(const Mesh<T, 3, Storage>& msh,
       const typename Mesh<T, 3, Storage>::cell& cl,
       const typename Mesh<T, 3, Storage>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() >= 3);

    auto v0 = (pts[1] - pts[0]).to_vector();
    auto v1 = (pts[2] - pts[1]).to_vector();
    auto n = v0.cross(v1);

    auto cell_bar = barycenter(msh, cl);
    auto face_bar = barycenter(msh, fc);
    auto outward_vector = (face_bar - cell_bar).to_vector();

    if ( n.dot(outward_vector) < T(0) )
        return -n/n.norm();

    return n/n.norm();
}

/**
 * @brief Compute a normal to a 3D face. No guarantee on the direction.
 *
 * @tparam Mesh type of mesh
 * @tparam T scalar type
 * @tparam Storage type of storage for the mesh
 * @param msh mesh
 * @param cl cell
 * @param fc face
 * @return static_vector<T, 3> A normal to a 3D face
 */
template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
static_vector<T, 3>
normal(const Mesh<T, 3, Storage>& msh,
       const typename Mesh<T, 3, Storage>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() >= 3);

    auto v0 = (pts[1] - pts[0]).to_vector();
    auto v1 = (pts[2] - pts[1]).to_vector();
    auto n = v0.cross(v1);

    return n/n.norm();
}

/**
 * @brief Compute the outward unit normal to a 1D face (i.e a point)
 *
 * @tparam Mesh type of mesh
 * @tparam T scalar type
 * @tparam Storage type of storage for the mesh
 * @param msh mesh
 * @param cl cell
 * @param fc face
 * @return T Compute the outward unit normal to a 1D face
 */
template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
T
normal(const Mesh<T, 1, Storage>& msh,
       const typename Mesh<T, 1, Storage>::cell& cl,
       const typename Mesh<T, 1, Storage>::face& fc)
{
    auto fcs = faces(msh, cl);
    assert(fcs.size() == 2);

    if (fc == fcs[0])
        return -1.;

    if (fc == fcs[1])
        return 1.;

    throw std::logic_error("shouldn't have arrived here");
}

/**
 * \brief Compute the diameter of the bounding box af a 3D cell, i.e, the maximum distance
 * between two different points of the bounding box.
 *
 * \param msh a reference to the mesh
 * \param cl a 3D cell
 * \return compute the diameter of the bounding box af a 3D cell
 *
 */
template<typename Mesh, typename Elem>
static_matrix<typename Mesh::coordinate_type, Mesh::dimension, Mesh::dimension>
inertia_axes(const Mesh& msh, const Elem& elem)
{
    typedef static_matrix<typename Mesh::coordinate_type, Mesh::dimension, Mesh::dimension> matrix_type;

    const auto bar = barycenter(msh, elem);
    const auto qps = integrate(msh, elem, 2);

    matrix_type mass_inertia = matrix_type::Zero();
    for (auto& qp : qps)
    {
        const auto coor = (bar - qp.point()).to_vector();

        mass_inertia += qp.weight() * coor * coor.transpose();
    }

    // std::cout << "inertia matrix:" << std::endl;
    // std::cout << mass_inertia << std::endl;

    if(mass_inertia.isDiagonal())
        return matrix_type::Identity();

    Eigen::SelfAdjointEigenSolver<matrix_type> es(mass_inertia);

    // std::cout << "eigenvalues:" << std::endl;
    // std::cout << es.eigenvalues().transpose() << std::endl;
    // std::cout << "eigenvectors:" << std::endl;
    // std::cout << es.eigenvectors() << std::endl;

    return es.eigenvectors();
}

} // namespace disk

#endif /* _GEOMETRY_ALL_HPP_ */
