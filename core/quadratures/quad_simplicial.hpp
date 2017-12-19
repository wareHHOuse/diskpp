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

#ifndef _QUADRATURES_HPP_WAS_INCLUDED_
#error "You must NOT include this file. Include quadratures.hpp"
#endif

#ifndef _QUAD_SIMPLICIAL_HPP_
#define _QUAD_SIMPLICIAL_HPP_

namespace disk {

template<typename T>
point<T, 2>
map_to_reference(const simplicial_mesh<T, 3>&                msh,
                 const typename simplicial_mesh<T, 3>::face& face,
                 const point<T, 3>&                          pm)
{
   const auto pts = points(msh, face);

   const auto v0  = (pts[1] - pts[0]).to_vector();
   const auto v1  = (pts[2] - pts[0]).to_vector();
   const auto v2  = (pm - pts[0]).to_vector();
   const auto d00 = v0.dot(v0);
   const auto d01 = v0.dot(v1);
   const auto d11 = v1.dot(v1);
   const auto d20 = v2.dot(v0);
   const auto d21 = v2.dot(v1);

   const auto invden = 1. / (d00 * d11 - d01 * d01);
   const auto v      = (d11 * d20 - d01 * d21) * invden;
   const auto w      = (d00 * d21 - d01 * d20) * invden;
   const auto u      = 1 - v - w;

   const auto x = u * 0 + v * 1 + w * 0;
   const auto y = u * 0 + v * 0 + w * 1;

   return point<T, 2>{x, y};
}

template<typename T>
point<T, 3>
map_to_physical(const simplicial_mesh<T, 3>&    msh,
                const simplicial_element<3, 1>& face,
                const point<T, 2>&              pm)
{
   const auto pts = points(msh, face);
   return pts[0] + (pts[1] - pts[0]) * pm.x() + (pts[2] - pts[0]) * pm.y();
}

template<typename T>
point<T, 3>
map_to_physical(const simplicial_mesh<T, 3>&                msh,
                const typename simplicial_mesh<T, 3>::cell& cell,
                const point<T, 3>&                          p)
{
   const auto pts = points(msh, cell);

   const auto pp =
     (pts[1] - pts[0]) * p.x() + (pts[2] - pts[0]) * p.y() + (pts[3] - pts[0]) * p.z() + pts[0];

   return pp;
}

template<typename T>
class quadrature<simplicial_mesh<T, 3>, typename simplicial_mesh<T, 3>::cell>
{
   size_t                                 m_order;
   std::vector<std::pair<point<T, 3>, T>> m_quadrature_data;

 public:
   typedef simplicial_mesh<T, 3>                mesh_type;
   typedef typename simplicial_mesh<T, 3>::cell cell_type;
   typedef quadrature_point<T, 3>               quadpoint_type;
   typedef point<T, 3>                          point_type;
   typedef T                                    weight_type;

   quadrature()
     : m_order(1)
   {
      m_quadrature_data = tetrahedron_quadrature(1);
   }

   quadrature(size_t order)
     : m_order(order)
   {
      m_quadrature_data = tetrahedron_quadrature(m_order);
   }

   size_t
   order() const
   {
      return m_order;
   }

   std::vector<quadpoint_type>
   integrate(const mesh_type& msh, const cell_type& cl) const
   {
      const auto meas = measure(msh, cl);

      auto tr = [&](const std::pair<point<T, 3>, T>& qd) -> auto
      {
         const auto point  = map_to_physical(msh, cl, qd.first);
         const auto weight = qd.second * meas;
         return make_qp(point, weight);
      };

      std::vector<quadpoint_type> ret(m_quadrature_data.size());
      std::transform(m_quadrature_data.begin(), m_quadrature_data.end(), ret.begin(), tr);

      return ret;
   }
};

template<typename T>
class quadrature<simplicial_mesh<T, 3>, typename simplicial_mesh<T, 3>::face>
{
   size_t                                 m_order;
   std::vector<std::pair<point<T, 2>, T>> m_quadrature_data;

 public:
   typedef simplicial_mesh<T, 3>                mesh_type;
   typedef typename simplicial_mesh<T, 3>::face face_type;
   typedef quadrature_point<T, 3>               quadpoint_type;
   typedef point<T, 2>                          point_type;
   typedef T                                    weight_type;

   quadrature()
     : m_order(1)
   {
      m_quadrature_data = triangle_quadrature(1);
   }

   quadrature(size_t order)
     : m_order(order)
   {
      m_quadrature_data = triangle_quadrature(m_order);
   }

   size_t
   order() const
   {
      return m_order;
   }

   std::vector<quadpoint_type>
   integrate(const mesh_type& msh, const face_type& fc) const
   {
      const auto meas = measure(msh, fc);

      auto tr = [&](const std::pair<point<T, 2>, T>& qd) -> auto
      {
         const auto point  = map_to_physical(msh, fc, qd.first);
         const auto weight = qd.second * meas;
         return make_qp(point, weight);
      };

      std::vector<quadpoint_type> ret(m_quadrature_data.size());
      std::transform(m_quadrature_data.begin(), m_quadrature_data.end(), ret.begin(), tr);

      return ret;
   }
};

template<typename T>
class quadrature<simplicial_mesh<T, 2>, typename simplicial_mesh<T, 2>::cell>
{
   size_t                                 m_order;
   std::vector<std::pair<point<T, 2>, T>> m_quadrature_data;

 public:
   typedef simplicial_mesh<T, 2>          mesh_type;
   typedef typename mesh_type::cell       cell_type;
   typedef quadrature_point<T, 2>         quadpoint_type;
   typedef typename mesh_type::point_type point_type;
   typedef T                              weight_type;

   quadrature()
     : m_order(1)
   {
      m_quadrature_data = triangle_quadrature(1);
   }

   quadrature(size_t order)
     : m_order(order)
   {
      m_quadrature_data = triangle_quadrature(m_order);
   }

   size_t
   order() const
   {
      return m_order;
   }

   std::vector<quadpoint_type>
   integrate(const mesh_type& msh, const cell_type& cl) const
   {
      const auto pts = points(msh, cl);

      std::vector<quadpoint_type> ret;

      ret.resize(m_quadrature_data.size());

      const auto col1 = pts[1] - pts[0];
      const auto col2 = pts[2] - pts[0];

      /* Compute the area of the sub-triangle */
      const auto tm = (col1.x() * col2.y() - col2.x() * col1.y()) / 2.;

      auto tr = [&](const std::pair<point<T, 2>, T>& qd) -> auto
      {
         const auto point  = col1 * qd.first.x() + col2 * qd.first.y() + pts[0];
         const auto weight = qd.second * std::abs(tm);
         return make_qp(point, weight);
      };

      std::transform(m_quadrature_data.begin(), m_quadrature_data.end(), ret.begin(), tr);

      return ret;
   }
};

template<typename T>
class quadrature<simplicial_mesh<T, 2>, typename simplicial_mesh<T, 2>::face>
{
   size_t                                 m_order;
   std::vector<std::pair<point<T, 1>, T>> m_quadrature_data;

 public:
   typedef simplicial_mesh<T, 2>          mesh_type;
   typedef typename mesh_type::face       face_type;
   typedef quadrature_point<T, 2>         quadpoint_type;
   typedef typename mesh_type::point_type point_type;
   typedef T                              weight_type;

   quadrature()
     : m_order(1)
   {
      m_quadrature_data = edge_quadrature<T>(1);
   }

   quadrature(size_t order)
     : m_order(order)
   {
      m_quadrature_data = edge_quadrature<T>(m_order);
   }

   size_t
   order() const
   {
      return m_order;
   }

   std::vector<quadpoint_type>
   integrate(const mesh_type& msh, const face_type& fc) const
   {
      const auto meas = measure(msh, fc);
      const auto pts  = points(msh, fc);
      auto       tr   = [&](const std::pair<point<T, 1>, T>& qd) -> auto
      {
         const auto point  = (pts[1] - pts[0]) * qd.first.x() + pts[0];
         const auto weight = qd.second * meas;
         return make_qp(point, weight);
      };

      std::vector<quadpoint_type> ret(m_quadrature_data.size());
      std::transform(m_quadrature_data.begin(), m_quadrature_data.end(), ret.begin(), tr);

      return ret;
   }
};

} // namespace disk

#endif /* _QUAD_SIMPLICIAL_HPP_ */
