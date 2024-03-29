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

#include "diskpp/bases/bases_traits.hpp"

namespace disk::basis {

template<typename Basis>
using piecewise_material_tensor =
    Eigen::Matrix<typename Basis::scalar_type, Basis::immersion_dimension,
        Basis::immersion_dimension>;

template<typename Basis>
using variable_material_tensor =
    std::function<piecewise_material_tensor<Basis>(typename Basis::point_type)>;

namespace priv {
template<typename Basis>
struct dot_evaluator
{
    const Basis& basis;
    using mesh_type = typename Basis::mesh_type;
    using coordinate_type = typename mesh_type::coordinate_type;
    using scalar_type = typename Basis::scalar_type;
    using point_type = typename mesh_type::point_type;
    static const size_t immersion_dimension = Basis::immersion_dimension;
    static const size_t basis_dimension = Basis::basis_dimension;

    static_assert(Basis::tensor_order > 0);
    static const size_t tensor_order = Basis::tensor_order-1;
    using value_type = typename tensor<scalar_type, basis_dimension, tensor_order>::value_type;
    using normal_type = Eigen::Matrix<coordinate_type, immersion_dimension, 1>;
    normal_type n;

    dot_evaluator() = delete;
    dot_evaluator(const dot_evaluator&) = delete;

    dot_evaluator(const Basis& pb, const normal_type& pn)
        : basis(pb), n(pn)
    {}

    auto operator()(const point_type& pt) const {
        typename basis_traits<Basis>::category basis_category;
        return eval_dot(pt, basis_category);
    }

    template<typename Range>
    auto operator()(const point_type& pt, const Range& dr) const {
        typename basis_traits<Basis>::category basis_category;
        return eval_dot(pt, dr, basis_category);
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

private:
    auto eval_dot(const point_type& pt, vector_basis_tag) const {
        return (basis(pt)*n).eval();
    }

    template<typename Range>
    auto eval_dot(const point_type& pt, const Range& dr, vector_basis_tag) const {
        return (basis(pt, dr)*n).eval();
    }

    //auto eval_dot(const point_type& pt, matrix_basis_tag) {
    //}
};

template<typename Basis>
struct grad_evaluator
{
    const Basis& basis;

    using mesh_type = typename Basis::mesh_type;
    using coordinate_type = typename mesh_type::coordinate_type;
    using scalar_type = typename Basis::scalar_type;
    using point_type = typename mesh_type::point_type;
    static const size_t immersion_dimension = Basis::immersion_dimension;
    static const size_t basis_dimension = Basis::basis_dimension;
    static const size_t tensor_order = Basis::tensor_order+1;
    using value_type = typename tensor<scalar_type, basis_dimension, tensor_order>::value_type;
    using normal_type = Eigen::Matrix<coordinate_type, immersion_dimension, 1>;
    normal_type n;

    grad_evaluator() = delete;
    grad_evaluator(const grad_evaluator&) = delete;

    grad_evaluator(const Basis& pb)
        : basis(pb)
    {}

    auto operator()(const point_type& pt) const {
        return basis.grad(pt);
    }

    template<typename Range>
    auto operator()(const point_type& pt, const Range& dr) const {
        return basis.grad(pt, dr);
    }

    auto dot(const normal_type& n) const {
        return dot_evaluator(*this, n);
    }

    size_t size() const {
        return basis.size();
    }

    size_t degree() const {
        return basis.degree();
    }

    size_t integration_degree() const {
        if (basis.degree() > 0)
            return basis.degree()-1;

        return basis.degree();
    }
};

template<typename Basis>
struct div_evaluator
{
    const Basis& basis;
    using mesh_type = typename Basis::mesh_type;
    using coordinate_type = typename mesh_type::coordinate_type;
    using scalar_type = typename Basis::scalar_type;
    using point_type = typename mesh_type::point_type;
    static const size_t immersion_dimension = Basis::immersion_dimension;
    static const size_t basis_dimension = Basis::basis_dimension;

    static_assert(Basis::tensor_order > 0, "Can't compute the divergence on a scalar basis");
    static const size_t tensor_order = Basis::tensor_order-1;
    using value_type = typename tensor<scalar_type, basis_dimension, tensor_order>::value_type;
    using sum_type = Eigen::Matrix<coordinate_type, immersion_dimension, 1>;

    div_evaluator() = delete;
    div_evaluator(const div_evaluator&) = delete;

    div_evaluator(const Basis& pb)
        : basis(pb)
    {}

    auto operator()(const point_type& pt) const {
        typename basis_traits<Basis>::category basis_category;
        return eval_div(pt, basis_category);
    }

    template<typename Range>
    auto operator()(const point_type& pt, const Range& dr) {
        typename basis_traits<Basis>::category basis_category;
        return eval_div(pt, dr, basis_category);
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

private:
    auto eval_div(const point_type& pt, vector_basis_tag) {
        throw "this implementation is just to test";
        sum_type sum = sum_type::Ones();
        return (basis(pt)*sum).eval();
    }

    template<typename Range>
    auto eval_div(const point_type& pt, const Range& dr, vector_basis_tag) {
        throw "this implementation is just to test";
        sum_type sum = sum_type::Ones();
        return (basis(pt, dr)*sum).eval();
    }

    //auto eval_div(const point_type& pt, matrix_basis_tag) {
    //}
};

template<typename Basis, typename MatTens>
struct material_prod_evaluator
{
    const Basis&    basis;
    const MatTens&  mt;

    using mesh_type = typename Basis::mesh_type;
    using coordinate_type = typename mesh_type::coordinate_type;
    using scalar_type = typename Basis::scalar_type;
    using point_type = typename mesh_type::point_type;
    static const size_t immersion_dimension = Basis::immersion_dimension;
    static const size_t basis_dimension = Basis::basis_dimension;

    static_assert(Basis::tensor_order == 1);
    static const size_t tensor_order = Basis::tensor_order;
    using value_type = typename tensor<scalar_type, immersion_dimension, tensor_order>::value_type;
    using normal_type = Eigen::Matrix<coordinate_type, immersion_dimension, 1>;

    material_prod_evaluator() = delete;
    material_prod_evaluator(const material_prod_evaluator&) = delete;

    material_prod_evaluator(const Basis& pb, const MatTens& pmt)
        : basis(pb), mt(pmt)
    {}

    auto operator()(const point_type& pt) const {
        return eval(pt, mt);
    }

    auto dot(const normal_type& n) const {
        return dot_evaluator(*this, n);
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

private:
    auto eval(const point_type& pt, const variable_material_tensor<Basis>& pmt) const {
        return (basis(pt)*pmt(pt)).eval();
    }

    auto eval(const point_type& pt, const piecewise_material_tensor<Basis>& pmt) const {
        return (basis(pt)*pmt).eval();
    }
};

} //namespace priv

template<typename Basis>
auto grad(const Basis& basis) {
    return priv::grad_evaluator(basis);
}

template<typename Basis>
auto div(const Basis& basis) {
    return priv::div_evaluator(basis);
}

template<vector_basis Basis>
auto
operator*(const variable_material_tensor<Basis>& mt, const Basis& basis)
{
    return priv::material_prod_evaluator(basis, mt);
}

template<vector_basis Basis>
auto
operator*(piecewise_material_tensor<Basis>& mt, const Basis& basis)
{
    return priv::material_prod_evaluator(basis, mt);
}


namespace priv {
template<typename Basis>
typename Basis::value_type
eval(const typename Basis::point_type& pt,
     const Eigen::Matrix<typename Basis::scalar_type, Eigen::Dynamic, 1>& dofs,
     const Basis& basis,
     scalar_basis_tag)
{
    assert(dofs.rows() == basis.size());
    auto phi = basis(pt);
    return phi.transpose()*dofs;
}

template<typename Basis>
typename Basis::value_type
eval(const typename Basis::point_type& pt,
     const Eigen::Matrix<typename Basis::scalar_type, Eigen::Dynamic, 1>& dofs,
     const Basis& basis,
     vector_basis_tag)
{
    assert(dofs.rows() == basis.size());
    auto phi = basis(pt);
    return phi.transpose()*dofs;
}

} //namespace priv

template<typename Basis>
typename Basis::value_type
eval(const typename Basis::point_type& pt,
     const Eigen::Matrix<typename Basis::scalar_type, Eigen::Dynamic, 1>& dofs,
     const Basis& basis)

{
    typename basis_traits<Basis>::category basis_category;
    return priv::eval(pt, dofs, basis, basis_category);
}

template<basis Test>
using fun = std::function<typename Test::value_type(typename Test::point_type)>;

template<basis Basis, typename Mesh, typename Element>
auto
L2_project(const Mesh& msh, const Element& elem,
    const fun<Basis>& f, const Basis& basis)
{
    auto mass = integrate(msh, elem, basis, basis);
    auto rhs = integrate(msh, elem, f, basis);
    return mass.ldlt().solve(rhs).eval();
}

} //namespace disk::basis
