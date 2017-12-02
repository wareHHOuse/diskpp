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

#pragma once

#include "bases/bases.hpp"

//#define USE_BLAS
#define FILL_COLMAJOR

namespace disk {

namespace hho {

template<typename Mesh,
         template<typename, typename> class Basis,
         template<typename, typename> class Quadrature>
class basis_quadrature_data /* this name really sucks */
{
 public:
   typedef Mesh                     mesh_type;
   typedef typename mesh_type::cell cell_type;
   typedef typename mesh_type::face face_type;

   typedef Basis<mesh_type, cell_type>      cell_basis_type;
   typedef Basis<mesh_type, face_type>      face_basis_type;
   typedef Quadrature<mesh_type, cell_type> cell_quad_type;
   typedef Quadrature<mesh_type, face_type> face_quad_type;

   cell_basis_type cell_basis;
   face_basis_type face_basis;
   cell_quad_type  cell_quadrature;
   face_quad_type  face_quadrature;
   face_quad_type  face_trace_quadrature;

 private:
   size_t m_cell_degree, m_face_degree;

   void init(void)
   {
      cell_basis            = cell_basis_type(m_cell_degree, m_cell_degree + 1);
      face_basis            = face_basis_type(m_face_degree);
      cell_quadrature       = cell_quad_type(2 * (m_cell_degree + 1));
      face_quadrature       = face_quad_type(2 * m_face_degree);
      face_trace_quadrature = face_quad_type(m_face_degree + m_cell_degree + 1);
   }

 public:
   basis_quadrature_data()
     : m_cell_degree(1)
     , m_face_degree(1)
   {
      init();
   }

   basis_quadrature_data(size_t cell_degree, size_t face_degree)
   {
      if ((cell_degree + 1 < face_degree) or (cell_degree > face_degree + 1))
         throw std::invalid_argument("Invalid cell degree");

      m_cell_degree = cell_degree;
      m_face_degree = face_degree;

      init();
   }

   size_t cell_degree(void) const { return m_cell_degree; }
   size_t face_degree(void) const { return m_face_degree; }
};

template<typename Mesh,
         template<typename, typename> class BasisFunction,
         template<typename, typename> class BasisGradient,
         template<typename, typename> class Quadrature>
class basis_quadrature_data_full /* this name really sucks */
{
 public:
   typedef Mesh                     mesh_type;
   typedef typename mesh_type::cell cell_type;
   typedef typename mesh_type::face face_type;

   typedef BasisFunction<mesh_type, cell_type> cell_basis_type;
   typedef BasisFunction<mesh_type, face_type> face_basis_type;
   typedef BasisGradient<mesh_type, cell_type> grad_basis_type;

   typedef Quadrature<mesh_type, cell_type> cell_quad_type;
   typedef Quadrature<mesh_type, face_type> face_quad_type;

   cell_basis_type cell_basis;
   face_basis_type face_basis;
   cell_quad_type  cell_quadrature;
   face_quad_type  face_quadrature;
   face_quad_type  face_trace_quadrature;
   grad_basis_type grad_basis;
   cell_quad_type  grad_quadrature;

   cell_quad_type grad_cell_quadrature;
   face_quad_type grad_face_quadrature;
   face_quad_type grad_face_max_quadrature;

 private:
   size_t m_cell_degree, m_face_degree, m_grad_degree;

   void init(void)
   {
      cell_basis            = cell_basis_type(m_cell_degree);
      face_basis            = face_basis_type(m_face_degree);
      grad_basis            = grad_basis_type(m_grad_degree);
      cell_quadrature       = cell_quad_type(2 * m_cell_degree);
      face_quadrature       = face_quad_type(2 * m_face_degree);
      grad_quadrature       = cell_quad_type(2 * m_grad_degree);
      face_trace_quadrature = face_quad_type(m_face_degree + m_cell_degree + 1);
      grad_cell_quadrature  = cell_quad_type(m_grad_degree + m_cell_degree);
      grad_face_quadrature  = face_quad_type(m_grad_degree + m_face_degree);
      grad_face_max_quadrature =
        face_quad_type(m_grad_degree + std::max(m_face_degree, m_cell_degree));
   }

 public:
   basis_quadrature_data_full()
     : m_cell_degree(1)
     , m_face_degree(1)
     , m_grad_degree(1)
   {
      init();
   }

   basis_quadrature_data_full(const size_t face_degree,
                              const size_t cell_degree,
                              const size_t grad_degree)
   {
      if ((cell_degree + 1 < face_degree) or (cell_degree > face_degree + 1))
         throw std::invalid_argument("Invalid cell degree");

      if (grad_degree < cell_degree) throw std::invalid_argument("Invalid grad degree");

      m_cell_degree = cell_degree;
      m_face_degree = face_degree;
      m_grad_degree = grad_degree;

      init();
   }

   void info_degree() const
   {
      std::cout << "Face degree: " << m_face_degree << std::endl;
      std::cout << "Cell degree: " << m_cell_degree << std::endl;
      std::cout << "Grad degree: " << m_grad_degree << std::endl;
   }

   size_t cell_degree(void) const { return m_cell_degree; }
   size_t face_degree(void) const { return m_face_degree; }
   size_t grad_degree(void) const { return m_grad_degree; }
};

// traits
template<typename BQData>
struct have_grad_space
{
   static const bool value = false;
};

template<typename Mesh,
         template<typename, typename> class BasisFunction,
         template<typename, typename> class BasisGradient,
         template<typename, typename> class Quadrature>
struct have_grad_space<basis_quadrature_data_full<Mesh, BasisFunction, BasisGradient, Quadrature>>
{
   static const bool value = true;
};

} // namespace hho

} // namespace disk
