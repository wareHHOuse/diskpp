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

#include "common/eigen.hpp"

namespace disk {

template<typename T>
T
mm_prod(const T& a, const T& b)
{
    return a*b;
}

template<typename T, int N>
T
mm_prod(const static_vector<T,N>& a, const static_vector<T,N>& b)
{
    return a.dot(b);
}

template<typename T, int N>
static_vector<T,N>
mm_prod(const static_vector<T,N>& a, const T& b)
{
    return a * b;
}

template<typename T, int N>
static_vector<T,N>
mm_prod(const T& a, const static_vector<T,N>& b)
{
    return a * b;
}

template<typename T, int N>
static_vector<T,N>
mm_prod(const static_matrix<T,N,N>& a, const static_vector<T,N>& b)
{
    return a * b;
}

template<typename T, int N>
T
mm_prod(const static_matrix<T,N,N>& a, const static_matrix<T,N,N>& b)
{
    return a.cwiseProduct(b).sum();
}

namespace priv {
struct precompute_function;
struct precompute_gradient;

template<typename Basis, typename precompute_what>
struct precomp_traits;

template<typename Basis>
struct precomp_traits<Basis, precompute_function>
{
    typedef typename Basis::function_value_type     value_type;
};

template<typename Basis>
struct precomp_traits<Basis, precompute_gradient>
{
    typedef typename Basis::gradient_value_type     value_type;
};

} // namespace priv

template<typename Mesh, typename Basis, typename what>
class precomputed_basis
{
public:
    typedef Mesh                                    mesh_type;
    typedef typename Mesh::scalar_type              scalar_type;
    typedef typename mesh_type::point_type          point_type;
    typedef Basis                                   basis_type;

    typedef typename priv::precomp_traits<Basis, what>::value_type      value_type;
    typedef std::vector<value_type>                                     vec_value_type;
    typedef std::vector<quadrature_point<scalar_type, Mesh::dimension>> vec_quadpoint_type;

private:

    struct internal_storage {
        std::shared_ptr<basis_type>     basis;
        std::vector<value_type>         data;
        vec_quadpoint_type              quad_points;

        internal_storage() {
            basis = std::make_shared<basis_type>();
            clear();
        }

        internal_storage(size_t degree) {
            basis = std::make_shared<basis_type>(degree);
            clear();
        }

        internal_storage(std::shared_ptr<basis_type> user_provided_basis) {
            basis = user_provided_basis;
        }

        void clear() {
            data.clear();
            quad_points.clear();
        }
    };

    std::shared_ptr<internal_storage>   m_storage;

    template<typename Elem, typename U = what>
    typename std::enable_if< std::is_same<U, priv::precompute_function>::value, vec_value_type >::type
    call_real_eval(const mesh_type& msh, const Elem& elem,
                   const point_type& point, const basis_type& basis)
    {
        return basis.eval_functions(msh, elem, point);
    }

    template<typename Elem, typename U = what>
    typename std::enable_if< std::is_same<U, priv::precompute_gradient>::value, vec_value_type >::type
    call_real_eval(const mesh_type& msh, const Elem& elem,
                   const point_type& point, const basis_type& basis)
    {
        return basis.eval_gradients(msh, elem, point);
    }

public:
    /* The precomputed basis uses its own basis object, but this can be
     * expensive because monomials are re-generated each time.
     */
    precomputed_basis()
    {
        m_storage = std::make_shared<internal_storage>();
    }

    /* The precomputed basis uses its own basis object, but this can be
     * expensive because monomials are re-generated each time.
     */
    precomputed_basis(size_t degree)
    {
        m_storage = std::make_shared<internal_storage>(degree);
    }

    /* In this case an user-supplied basis object is used, so monomials aren't
     * re-generated again and again.
     */
    precomputed_basis(std::shared_ptr<basis_type> user_provided_basis)
    {
        m_storage = std::make_shared<internal_storage>(user_provided_basis);
    }

    template<typename Elem>
    void
    compute(const mesh_type& msh, const Elem& elem,
            const vec_quadpoint_type& qpts)
    {
        assert(m_storage);
        assert(m_storage->basis);
        m_storage->clear();
        auto basis_size = m_storage->basis->size();
        m_storage->data.reserve( qpts.size() * basis_size );
        assert(m_storage->data.size() == 0);
        for (auto& qpt : qpts)
        {
            auto pt = qpt.point();
            auto eval_result = call_real_eval(msh, elem, pt, *m_storage->basis);
            m_storage->data.insert(m_storage->data.end(), eval_result.begin(), eval_result.end());
        }

        assert(m_storage->quad_points.size() == 0);
        m_storage->quad_points.insert(m_storage->quad_points.end(), qpts.begin(), qpts.end());
    }

    dof_range range(size_t min_degree, size_t max_degree) const {
        assert(m_storage);
        assert(m_storage->basis);
        return m_storage->basis->range(min_degree, max_degree);
    }

    dof_range full_range() const {
        assert(m_storage);
        assert(m_storage->basis);
        return m_storage->basis->full_range();
    }

    size_t degree_index(size_t degree) const {
        assert(m_storage);
        assert(m_storage->basis);
        return m_storage->basis->degree_index(degree);
    }

    size_t size() const {
        assert(m_storage);
        assert(m_storage->basis);
        return m_storage->basis->size();
    }

    size_t size_eval_points() const {
        assert(m_storage);
        return m_storage->quad_points.size();
    }

    auto data_begin() const {
        assert(m_storage);
        return m_storage->data.begin();
    }

    auto data_end() const {
        assert(m_storage);
        return m_storage->data.end();
    }

    auto eval_points_begin() const {
        assert(m_storage);
        return m_storage->quad_points.begin();
    }

    auto eval_points_end() const {
        assert(m_storage);
        return m_storage->quad_points.end();
    }
};

template<typename M, typename B>
using precomputed_basis_functions = precomputed_basis<M, B, priv::precompute_function>;

template<typename M, typename B>
using precomputed_basis_gradients = precomputed_basis<M, B, priv::precompute_gradient>;

} // namespace disk
