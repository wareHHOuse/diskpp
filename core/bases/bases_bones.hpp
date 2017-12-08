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

/* Modification needed for MsHHO:
 * A basis has a degree and a computed_degree and then a size() and a
 * computed_size(). The basis is then treated as a basis of order degree,
 * except that it is usable up to computed_degree.
 */

#ifndef _BASES_HPP_WAS_INCLUDED_
#error "You must NOT include this file directly. Include bases.hpp"
#endif

#ifndef _BASES_BONES_HPP_
#define _BASES_BONES_HPP_

#define VERY_HIGH_DEGREE 42424242 // this is shit in its purest form

#include "bases/bases_ranges.hpp"
#include "bases/monomial_generator.hpp"

namespace disk {

namespace priv {

template<size_t DIM, size_t MULT = 1>
class monomial_basis_bones
{
   monomial_generator<DIM> m_monomials;
   size_t                  m_degree, m_computed_degree;

   static_assert((not(DIM == 0)) or (MULT == 1), "Only multiplicity of one with zero-dim basis");

 protected:
   auto
   monomials_begin() const
   {
      return m_monomials.begin();
   };

   auto
   monomials_end() const
   {
      return m_monomials.end();
   };

 public:
   static const size_t multiplicity = MULT;

   monomial_basis_bones() : m_degree(1), m_computed_degree(1)
   {
      m_monomials = monomial_generator<DIM>(1);
   }

   monomial_basis_bones(size_t degree) : m_degree(degree), m_computed_degree(degree)
   {
      m_monomials = monomial_generator<DIM>(m_degree);
   }

   monomial_basis_bones(size_t degree, size_t computed_degree) :
     m_degree(degree), m_computed_degree(computed_degree)
   {
      m_monomials = monomial_generator<DIM>(m_computed_degree);
   }

   size_t
   degree_index(size_t degree) const
   {
      if (degree > m_computed_degree)
         throw std::invalid_argument("Specified degree too high for this basis");
      return MULT * m_monomials.degree_position(degree);
   }

   size_t
   size() const
   {
      return MULT * m_monomials.degree_position(m_degree + 1);
   }

   size_t
   computed_size() const
   {
      return MULT * m_monomials.size();
   }

   [[deprecated("You should use degree/computed_degree interface")]] size_t
   max_degree() const
   {
      return degree();
   }

   size_t
   degree() const
   {
      return m_degree;
   }

   size_t
   computed_degree() const
   {
      return m_computed_degree;
   }

   dof_range
   range(size_t min_degree, size_t max_degree) const
   {
      if (min_degree > m_computed_degree || max_degree > m_computed_degree)
         throw std::invalid_argument("Specified degree too high for this basis");
      auto pos_min = MULT * m_monomials.degree_position(min_degree);
      auto pos_max = MULT * m_monomials.degree_position(max_degree + 1);
      assert(pos_max > 0);
      return dof_range(pos_min, pos_max);
   }

   [[deprecated("You should use range/computed_range interface")]] dof_range
   full_range() const
   {
      return dof_range(0, size());
   }

   dof_range
   range() const
   {
      return dof_range(0, size());
   }

   dof_range
   computed_range() const
   {
      return dof_range(0, computed_size());
   }
};

template<>
class monomial_basis_bones<0, 1>
{
 public:
   static const size_t multiplicity = 1;

   monomial_basis_bones() {}
   monomial_basis_bones(size_t) {}
   monomial_basis_bones(size_t, size_t) {}
   size_t
   degree_index(size_t degree) const
   {
      return 0;
   }
   size_t
   size() const
   {
      return 1;
   }
   size_t
   computed_size() const
   {
      return 1;
   }
   dof_range range(size_t, size_t) const { return dof_range(0, 1); }
   dof_range
   full_range() const
   {
      return dof_range(0, 1);
   }
};

} // namespace priv

#define POWER_CACHE_MAX_POWER 7

template<typename T, size_t DIM, size_t MAX_POWER = POWER_CACHE_MAX_POWER>
class power_cache
{
   std::array<T, DIM*(MAX_POWER + 1)> m_powers;
   size_t                             m_max_computed_power;

 public:
   power_cache()
   {
      std::fill(m_powers.begin(), m_powers.end(), T(1));
      m_max_computed_power = 0;
   }

   power_cache(const point<T, DIM>& pt, size_t max_degree)
   {
      if (max_degree > MAX_POWER)
         throw std::invalid_argument("Power cache not sufficiently large. Recompile with the "
                                     "desired MAX_POWER or without power cache.");

      m_max_computed_power = max_degree;

      for (size_t i = 0; i < DIM; i++) {
         m_powers[i * (MAX_POWER + 1)] = T(1);
         for (size_t j = 1; j < m_max_computed_power; j++)
            m_powers[i * (MAX_POWER + 1) + j] = m_powers[i * (MAX_POWER + 1) + j - 1] * pt[i];
      }
   }

   template<typename U = T>
   typename std::enable_if<DIM == 1 || DIM == 2 || DIM == 3, U>::type
   x(size_t power) const
   {
      assert(power <= m_max_computed_power && "power too large");
      return m_powers[power];
   }

   template<typename U = T>
   typename std::enable_if<DIM == 2 || DIM == 3, U>::type
   y(size_t power) const
   {
      assert(power <= m_max_computed_power && "power too large");
      return m_powers[(MAX_POWER + 1) + power];
   }

   template<typename U = T>
   typename std::enable_if<DIM == 3, U>::type
   z(size_t power) const
   {
      assert(power <= m_max_computed_power && "power too large");
      return m_powers[2 * (MAX_POWER + 1) + power];
   }
};

/*
template<typename T, size_t M, size_t N>
class evaluated_basis
{
    std::vector<static_matrix<T, M, N>>     m_basis_data;

public:
    evaluated_basis() {}

};

template<typename T, size_t M>
class evaluated_basis<T,M,1>
{
    Eigen::Matrix<T, Eigen::Dynamic, M>     m_basis_data;

public:
    evaluated_basis() {}
};
*/

} // namespace disk

#endif /* _BASES_BONES_HPP_ */
