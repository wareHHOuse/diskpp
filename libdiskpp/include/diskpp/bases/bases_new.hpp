/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#pragma once

#include "diskpp/bases/bases_raw.hpp"

#include "diskpp/common/eigen.hpp"
#include "diskpp/mesh/mesh.hpp"

namespace disk {

namespace basis {







struct scalar_basis {};
struct vector_basis {};
struct matrix_basis {};

template<typename Basis>
struct basis_traits
{
    using category = typename Basis::category;
};





template<typename Mesh, typename Element,
    typename ScalT = typename Mesh::coordinate_type>
class scalar_monomial {
    
    scaled_monomial_basis<

public:
    using mesh_type = Mesh;
    using point_type = typename mesh_type::point_type;
    using basis_category = scalar_basis;

    scalar_type operator()(const point_type& pt) const {
        return 1.0;
    }

    scalar_type grad(const point_type& pt) const {
        return 42.0;
    }

    size_t size() const {
        return 18;
    }

    size_t degree() const {
        return 26;
    }

    size_t integration_degree() const {
        return degree;
    }
};

namespace priv {
template<typename Basis>
struct dot_evaluator
{
    const Basis& basis;
    double n;

    using point_type = typename Basis::point_type;

    dot_evaluator() = delete;
    dot_evaluator(const dot_evaluator&) = delete;

    dot_evaluator(const Basis& pb, double pn)
        : basis(pb), n(pn)
    {}

    auto operator()(const point_type& pt) const {
        return basis(pt)*n;
    }

    size_t size() const {
        return basis.size();
    }

    size_t degree() const {
        return basis.degree();
    }

    size_t integration_degree() const {
        return basis.degree();
    }
};

template<typename Basis>
struct grad_evaluator
{
    const Basis& basis;

    using mesh_type = typename Basis::mesh_type;
    using point_type = typename Basis::point_type;

    grad_evaluator() = delete;
    grad_evaluator(const grad_evaluator&) = delete;

    grad_evaluator(const Basis& pb)
        : basis(pb)
    {}

    auto operator()(const point_type& pt) const {
        return basis.grad(pt);
    }

    auto dot(const double& n) {
        return dot_evaluator(*this, n);
    }

    size_t size() const {
        return basis.size();
    }

    size_t degree() const {
        return basis.degree();
    }

    size_t integration_degree() const {
        return basis.degree()-1;
    }
};

} //namespace priv

template<typename Basis>
auto grad(const Basis& basis) {
    return priv::grad_evaluator(basis);
}

template<typename Trial, typename Test>
struct bilinear_form
{
    const Trial&   trial;
    const Test&    test;

    bilinear_form() = delete;
    bilinear_form(const bilinear_form&) = delete;

    bilinear_form(const Trial& ptrial, const Test& ptest)
        : trial(ptrial), test(ptest)
    {
        std::cout << "bf" << std::endl;
    }
};


} // namespace basis

} // namespace disk