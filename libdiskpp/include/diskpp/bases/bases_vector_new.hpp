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

#include "diskpp/common/eigen.hpp"
#include "diskpp/bases/bases_traits.hpp"
#include "diskpp/quadratures/quadrature_point.hpp"

namespace disk::basis {

/**********************************************************************************************
 * Cells, 2D
 */
template<template<typename, size_t, typename> class Mesh,
    typename CoordT, typename Storage, typename ScalT>
class vector_monomial<Mesh<CoordT, 2, Storage>, typename Mesh<CoordT, 2, Storage>::cell_type, ScalT>
{
public:
    using mesh_type = Mesh<CoordT, 2, Storage>;
    using cell_type = typename mesh_type::cell_type;
    using element_type = cell_type;
    using coordinate_type = typename mesh_type::coordinate_type;
    using scalar_type = ScalT;
    using point_type = typename mesh_type::point_type;
    static const size_t immersion_dimension = mesh_type::dimension;
    static const size_t basis_dimension = mesh_type::dimension;
    static const size_t tensor_order = 1; /* vector */
    using value_type = typename tensor<scalar_type, basis_dimension, tensor_order>::value_type;
    using value_array_type = typename tensor<scalar_type, basis_dimension, tensor_order>::array_type;
    using gradient_array_type = typename tensor<scalar_type, basis_dimension, tensor_order+1>::array_type;

private:
    using tr_mat_type = Eigen::Matrix<coordinate_type, immersion_dimension, immersion_dimension>;
    using vec_type = Eigen::Matrix<coordinate_type, immersion_dimension, 1>;

    point_type          bar_;
    size_t              degree_;
    size_t              size_;
    tr_mat_type         tr_;
    vec_type            v0_, v1_;
};

} // namespace disk::basis