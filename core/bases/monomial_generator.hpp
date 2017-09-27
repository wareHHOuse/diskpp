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

#include <vector>
#include <array>
#include <type_traits>

namespace disk {

template<size_t DIM>
class monomial_generator
{
    typedef std::array<size_t, DIM>     monomial_type;
    size_t                              m_max_degree;
    std::vector<size_t>                 m_degree_indices;
    std::vector<monomial_type>          m_monomials;

    struct monomial_comparator
    {
        bool operator()(const monomial_type& a, const monomial_type& b)
        {
            auto sum_a = std::accumulate(std::next(a.begin()), a.end(), a.front());
            auto sum_b = std::accumulate(std::next(b.begin()), b.end(), b.front());

            size_t val_a = 0, val_b = 0;
            for (size_t i = 0, pow_i = 1; i < DIM; i++, pow_i *= 10)
            {
                val_a += a[i]*pow_i;
                val_b += b[i]*pow_i;
            }

            if (sum_a == sum_b)
                return val_a < val_b;

            return (sum_a < sum_b);
        }
    };

    void
    make_index_table(void)
    {
        size_t  prev_degree = 0;
        size_t  cur_deg = 1;
        m_degree_indices.resize(m_max_degree+1);
        for (size_t i = 0; i < m_monomials.size(); i++)
        {
            auto m = m_monomials[i];
            auto msum = std::accumulate(std::next(m.begin()), m.end(), m.front());
            if (msum > prev_degree)
                m_degree_indices.at(cur_deg++) = i;

            prev_degree = msum;
        }
        m_degree_indices.push_back(m_monomials.size());
    }

    template<size_t BDIM = DIM>
    typename std::enable_if<BDIM == 1>::type
    generate_monomials(size_t max_degree)
    {
        monomial_type m; /* degrees of x */

        m_degree_indices.clear();
        m_monomials.clear();

        for (size_t i = 0; i <= max_degree; i++)
        {
            m[0] = i;
            m_monomials.push_back( m );
            m_degree_indices.push_back(i);
        }
        m_degree_indices.push_back(m_monomials.size());
        assert(m_monomials.size() == max_degree+1);
    }

    template<size_t BDIM = DIM>
    typename std::enable_if<BDIM == 2>::type
    generate_monomials(size_t max_degree)
    {
        monomial_type m; /* degrees of x and y */

        for (size_t k = 0; k <= max_degree; k++)
        {
            for (size_t i = 0; i <= k; i++)
            {
                m[0] = i;      //x degree
                m[1] = k-i;    //y degree

                m_monomials.push_back( m );
            }
        }

        std::sort(m_monomials.begin(), m_monomials.end(), monomial_comparator{});
        assert( m_monomials.size() == binomial(max_degree+2,2) );
        make_index_table();
    }

    template<size_t BDIM = DIM>
    typename std::enable_if<BDIM == 3>::type
    generate_monomials(size_t max_degree)
    {
        /* The idea here is that to get *the degree K* of the 3D basis we take the
         * whole 2D basis and we multiply its components by the appropriate powers
         * of z. In particular, components of degree 0 are multiplied by z^k,
         * components of degree 1 are multiplied by z^(k-1) and so on, as
         * depicted below:
         *
         *                  1             *  z^2    = z^2
         *             x        y         *  z      = xz, yz
         *         x^2     xy     y^2     *  z^0    = x^2, xy, y^2
         *
         * In this way we got *the degree K* of the 3D basis. We must then iterate
         * on all the degrees.
         */

        monomial_type m; /* degrees of x, y and z */

        for (size_t k = 0; k <= max_degree; k++)
        {
            m[2] = max_degree - k; //z degree

            for (size_t j = 0; j <= k; j++)
            {
                for (size_t i = 0; i <= j; i++)
                {
                    m[0] = i;     //x degree
                    m[1] = j-i;   //y degree

                    m_monomials.push_back( m );
                }
            }
        }

        std::sort(m_monomials.begin(), m_monomials.end(), monomial_comparator{});
        assert( m_monomials.size() == binomial(max_degree+3,3) );
        make_index_table();
    }

public:
    monomial_generator()
        : m_max_degree(1)
    {
        generate_monomials(m_max_degree);
    }

    monomial_generator(size_t max_degree)
        : m_max_degree(max_degree)
    {
        generate_monomials(m_max_degree);
    }

    auto begin() const {
        return m_monomials.begin();
    }

    auto end() const {
        return m_monomials.end();
    }

    size_t size() const {
        return m_monomials.size();
    }

    size_t degree_position(size_t degree) const {
        return m_degree_indices.at(degree);
    }

    size_t max_degree() const {
        return m_max_degree;
    }
};

} // namespace disk
