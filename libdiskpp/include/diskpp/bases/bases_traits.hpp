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

namespace disk::basis {

template<typename T>
concept basis = requires(T t) {
    T::immersion_dimension;
    T::basis_dimension;
    T::tensor_order;
    t.size();
    t.degree();
    t.integration_degree();
};

template<typename T>
concept vector_basis = basis<T> && requires {
    requires T::tensor_order == 1;
};

template<typename Trial, typename Test>
struct can_take_scalar_product
{};

template<basis Trial, basis Test>
struct can_take_scalar_product<Trial, Test>
    : std::bool_constant<Trial::tensor_order == Test::tensor_order> { };

template<size_t order>
struct basis_category_tag
{};

using scalar_basis_tag = basis_category_tag<0>;
using vector_basis_tag = basis_category_tag<1>;
using matrix_basis_tag = basis_category_tag<2>;

template<typename Basis>
struct basis_traits
{
    using category = basis_category_tag<Basis::tensor_order>;
};






template<typename T, size_t dimension, size_t order>
struct tensor {
    static_assert(sizeof(T) == -1, "invalid tensor dimension or order");
};

template<typename T, size_t DIM>
struct tensor<T, DIM, 0> {
    static const size_t order = 0;
    static const size_t dimension = DIM;
    using scalar_type = T;
    using value_type = scalar_type;
    using array_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
};

template<typename T, size_t DIM>
struct tensor<T, DIM, 1> {
    static const size_t order = 1;
    static const size_t dimension = DIM;
    using scalar_type = T;
    using value_type = Eigen::Matrix<T, dimension, 1>;
    using array_type = Eigen::Matrix<T, Eigen::Dynamic, dimension>;
};

template<typename T>
struct tensor<T, 0, 0> {
    static const size_t order = 0;
    static const size_t dimension = 0;
    using scalar_type = T;
    using value_type = scalar_type;
    using array_type = scalar_type;
};

template<typename T>
struct tensor<T, 0, 1> {
    static const size_t order = 1;
    static const size_t dimension = 0;
    using scalar_type = T;
    using value_type = Eigen::Matrix<T, 1, 1>;
    using array_type = Eigen::Matrix<T, 1, 1>;
};

} // namespace disk::basis