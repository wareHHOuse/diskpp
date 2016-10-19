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

#ifndef _BASES_HPP_WAS_INCLUDED_
    #error "You must NOT include this file directly. Include bases.hpp"
#endif

#ifndef _BASES_BONES_HPP_
#define _BASES_BONES_HPP_

#define POWER_CACHE
#define VERY_HIGH_DEGREE 42424242   // this is shit in it's purest form

#include "bases/monomial_generator.hpp"

namespace disk {

namespace priv {

template<size_t DIM, size_t MULT = 1>
class monomial_basis_bones
{
    monomial_generator<DIM>     m_monomials;
    size_t                      m_degree;

    static_assert((not (DIM == 0)) or (MULT == 1), "Only multiplicity of one with zero-dim basis");

protected:
    auto monomials_begin() const {
        return m_monomials.begin();
    };

    auto monomials_end() const {
        return m_monomials.end();
    };

public:
    static const size_t         multiplicity = MULT;

    monomial_basis_bones()
        : m_degree(1)
    {
        m_monomials = monomial_generator<DIM>(1);
    }

    monomial_basis_bones(size_t degree)
        : m_degree(degree)
    {
        m_monomials = monomial_generator<DIM>(m_degree);
    }

    size_t degree_index(size_t degree) const {
        if (degree > m_degree)
            throw std::invalid_argument("Specified degree too high for this basis");
        return MULT * m_monomials.degree_position(degree);
    }

    size_t size() const {
        return MULT * m_monomials.size();
    }

    size_t max_degree() const {
        return m_monomials.max_degree();
    }

    dof_range range(size_t min_degree, size_t max_degree) const
    {
        if (min_degree > m_degree || max_degree > m_degree)
            throw std::invalid_argument("Specified degree too high for this basis");
        auto pos_min = MULT * m_monomials.degree_position(min_degree);
        auto pos_max = MULT * m_monomials.degree_position(max_degree+1);
        assert(pos_max > 0);
        return dof_range(pos_min, pos_max);
    }

    dof_range full_range() const
    {
        return dof_range(0, size());
    }
};

template<>
class monomial_basis_bones<0, 1>
{
public:
    static const size_t         multiplicity = 1;

    monomial_basis_bones(){}
    monomial_basis_bones(size_t){}
    size_t degree_index(size_t degree) const { return 0; }
    size_t size() const { return 1; }
    dof_range range(size_t, size_t) const { return dof_range(0,1); }
    dof_range full_range() const { return dof_range(0,1); }

};

} //namespace priv
} //namespace disk

#endif /* _BASES_BONES_HPP_ */
