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
#include "quadratures/raw_simplices.hpp"

namespace disk {

/* Compute an estimate of the mesh discretization step 'h' */
template<typename Mesh>
typename Mesh::scalar_type
mesh_h(const Mesh& msh)
{
   typename Mesh::scalar_type h{};
   for (auto itor = msh.cells_begin(); itor != msh.cells_end(); itor++) {
      auto       cell         = *itor;
      const auto cell_measure = measure(msh, cell);

      const auto                 fcs = faces(msh, cell);
      typename Mesh::scalar_type face_sum{};
      for (auto& f : fcs) {
         const auto m = measure(msh, f);
         face_sum += m;
      }
      h = std::max(h, cell_measure / face_sum);
   }

   return h;
}

template<typename Mesh>
typename Mesh::scalar_type
average_diameter(const Mesh& msh)
{
   typename Mesh::scalar_type h{};
   for (auto& cl : msh) {
      h += diameter(msh, cl);
   }

   return h / msh.cells_size();
}

template<typename Mesh, typename Element>
std::vector<typename Mesh::point_type>
points(const Mesh& msh, const Element& elem)
{
   const auto ptids = elem.point_ids();

   const auto points_begin  = msh.points_begin();
   const auto ptid_to_point = [&](const point_identifier<Mesh::dimension>& pi) -> auto
   {
      return *std::next(points_begin, pi);
   };

   std::vector<typename Mesh::point_type> pts(ptids.size());
   std::transform(ptids.begin(), ptids.end(), pts.begin(), ptid_to_point);

   return pts;
}

/* Compute the barycenter of a cell */
template<typename Mesh, typename Element>
point<typename Mesh::coordinate_type, Mesh::dimension>
barycenter(const Mesh& msh, const Element& elm)
{
   const auto pts = points(msh, elm);
   const auto bar = std::accumulate(std::next(pts.begin()), pts.end(), pts.front());
   return bar / typename Mesh::coordinate_type(pts.size());
}

template<typename Mesh, typename Element>
typename Mesh::scalar_type
diameter(const Mesh& msh, const Element& elem)
{
   const auto pts = points(msh, elem);

   typename Mesh::scalar_type diam = 0.;

   for (size_t i = 0; i < pts.size(); i++)
      for (size_t j = i + 1; j < pts.size(); j++)
         diam = std::max((pts[i] - pts[j]).to_vector().norm(), diam);

   return diam;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
bool
is_inside(const Mesh<T, 2, Storage>&                      msh,
          const typename Mesh<T, 2, Storage>::cell&       cl,
          const typename Mesh<T, 2, Storage>::point_type& pt)
{
   /* Nodes MUST be ordered COUNTERCLOCKWISE and the polygon must be CONVEX */
   const auto pts = points(msh, cl);

   for (size_t i = 0; i < pts.size(); i++) {
      const auto p0 = pts[i];
      const auto p1 = pts[i % pts.size()];

      const auto x  = pt.x();
      const auto y  = pt.y();
      const auto x0 = p0.x();
      const auto y0 = p0.y();
      const auto x1 = p1.x();
      const auto y1 = p1.y();

      if ((y - y0) * (x1 - x0) - (x - x0) * (y1 - y0) < 0.0) return false;
   }

   return true;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
bool
is_inside(const Mesh<T, 3, Storage>&                      msh,
          const typename Mesh<T, 3, Storage>::cell&       cl,
          const typename Mesh<T, 3, Storage>::point_type& pt)
{
   // split in tetrahedra
   const auto rss = split_in_raw_tetrahedra(msh, cl);

   for (auto& rs : rss) {
      const auto pts = rs.points();
      assert(pts.size() == 4);

      const auto p0 = pts[0];
      const auto p1 = pts[1];
      auto       p2 = pts[2];
      auto       p3 = pts[3];

      typedef static_matrix<T, 4, 4> mat4;

      // d0
      mat4 m0;
      m0(0, 0) = p0.x();
      m0(0, 1) = p0.y();
      m0(0, 2) = p0.z();
      m0(0, 3) = 1;
      m0(1, 0) = p1.x();
      m0(1, 1) = p1.y();
      m0(1, 2) = p1.z();
      m0(1, 3) = 1;
      m0(2, 0) = p2.x();
      m0(2, 1) = p2.y();
      m0(2, 2) = p2.z();
      m0(2, 3) = 1;
      m0(3, 0) = p3.x();
      m0(3, 1) = p3.y();
      m0(3, 2) = p3.z();
      m0(3, 3) = 1;

      auto d0 = m0.determinant();
      assert(d0 != T(0));

      if (d0 < T(0)) {
         d0             = std::abs(d0);
         const auto tmp = p2;
         p2             = p3;
         p3             = tmp;
      }

      // d1
      m0(0, 0) = pt.x();
      m0(0, 1) = pt.y();
      m0(0, 2) = pt.z();
      m0(0, 3) = 1;
      m0(1, 0) = p1.x();
      m0(1, 1) = p1.y();
      m0(1, 2) = p1.z();
      m0(1, 3) = 1;
      m0(2, 0) = p2.x();
      m0(2, 1) = p2.y();
      m0(2, 2) = p2.z();
      m0(2, 3) = 1;
      m0(3, 0) = p3.x();
      m0(3, 1) = p3.y();
      m0(3, 2) = p3.z();
      m0(3, 3) = 1;

      const auto d1 = m0.determinant();

      if (d1 == T(0)) return true;
      if (d0 / d1 > T(0)) {
         // d2
         m0(0, 0) = p0.x();
         m0(0, 1) = p0.y();
         m0(0, 2) = p0.z();
         m0(0, 3) = 1;
         m0(1, 0) = pt.x();
         m0(1, 1) = pt.y();
         m0(1, 2) = pt.z();
         m0(1, 3) = 1;
         m0(2, 0) = p2.x();
         m0(2, 1) = p2.y();
         m0(2, 2) = p2.z();
         m0(2, 3) = 1;
         m0(3, 0) = p3.x();
         m0(3, 1) = p3.y();
         m0(3, 2) = p3.z();
         m0(3, 3) = 1;

         const auto d2 = m0.determinant();

         if (d2 == T(0)) return true;
         if (d0 / d2 > T(0)) {
            // d3
            m0(0, 0) = p0.x();
            m0(0, 1) = p0.y();
            m0(0, 2) = p0.z();
            m0(0, 3) = 1;
            m0(1, 0) = p1.x();
            m0(1, 1) = p1.y();
            m0(1, 2) = p1.z();
            m0(1, 3) = 1;
            m0(2, 0) = pt.x();
            m0(2, 1) = pt.y();
            m0(2, 2) = pt.z();
            m0(2, 3) = 1;
            m0(3, 0) = p3.x();
            m0(3, 1) = p3.y();
            m0(3, 2) = p3.z();
            m0(3, 3) = 1;

            const auto d3 = m0.determinant();

            if (d3 == T(0)) return true;
            if (d0 / d3 > T(0)) {
               // d4
               m0(0, 0) = p0.x();
               m0(0, 1) = p0.y();
               m0(0, 2) = p0.z();
               m0(0, 3) = 1;
               m0(1, 0) = p1.x();
               m0(1, 1) = p1.y();
               m0(1, 2) = p1.z();
               m0(1, 3) = 1;
               m0(2, 0) = p2.x();
               m0(2, 1) = p2.y();
               m0(2, 2) = p2.z();
               m0(2, 3) = 1;
               m0(3, 0) = pt.x();
               m0(3, 1) = pt.y();
               m0(3, 2) = pt.z();
               m0(3, 3) = 1;

               const auto d4 = m0.determinant();

               if (d4 == T(0)) return true;
               if (d0 / d4 > T(0)) return true;
            }
         }
      }
   }
   return false;
}

template<typename Mesh>
bool
has_faces_on_boundary(const Mesh& msh, const typename Mesh::cell& cl)
{
   const auto fcs     = faces(msh, cl);
   bool       has_bnd = false;
   for (auto& fc : fcs)
      if (msh.is_boundary(fc)) has_bnd = true;

   return has_bnd;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
static_vector<T, 2>
normal(const Mesh<T, 2, Storage>&                msh,
       const typename Mesh<T, 2, Storage>::cell& cl,
       const typename Mesh<T, 2, Storage>::face& fc)
{
   const auto pts = points(msh, fc);
   assert(pts.size() == 2);

   const auto v = pts[1] - pts[0];
   const auto n = (point<T, 2>({-v.y(), v.x()})).to_vector();

   const auto cell_bar       = barycenter(msh, cl);
   const auto face_bar       = barycenter(msh, fc);
   const auto outward_vector = (face_bar - cell_bar).to_vector();

   if (n.dot(outward_vector) < T(0)) return -n / n.norm();

   return n / n.norm();
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
static_vector<T, 3>
normal(const Mesh<T, 3, Storage>&                msh,
       const typename Mesh<T, 3, Storage>::cell& cl,
       const typename Mesh<T, 3, Storage>::face& fc)
{
   const auto pts = points(msh, fc);
   assert(pts.size() >= 3);

   const auto v0 = (pts[1] - pts[0]).to_vector();
   const auto v1 = (pts[2] - pts[1]).to_vector();
   const auto n  = v0.cross(v1);

   const auto cell_bar       = barycenter(msh, cl);
   const auto face_bar       = barycenter(msh, fc);
   const auto outward_vector = (face_bar - cell_bar).to_vector();

   if (n.dot(outward_vector) < T(0)) return -n / n.norm();

   return n / n.norm();
}

} // namespace disk

#endif /* _GEOMETRY_ALL_HPP_ */
