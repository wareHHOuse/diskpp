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





template<typename Mesh, typename Element, typename ScalT>
class scalar_monomial;




template<template<typename, size_t, typename> class Mesh,
    typename CoordT, typename Storage, typename ScalT>
class scalar_monomial<Mesh<CoordT, 1, Storage>, typename Mesh<CoordT, 1, Storage>::cell_type, ScalT>
{
public:
    using mesh_type = Mesh<CoordT, 1, Storage>;
    using cell_type = typename mesh_type::cell_type;
    using point_type = typename mesh_type::point_type;
    static const size_t dimension = 1;
    using coordinate_type = typename mesh_type::coordinate_type;
    using scalar_type = ScalT;
    using basis_category = scalar_basis;
    using function_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
    using gradient_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;

private:
    point_type      pa_, pb_;
    point_type      bar_;
    coordinate_type h_;
    size_t          degree_;
    size_t          size_;

    point_type phys2ref(const point_type& pt) {
        auto bv = (pb_ - bar_).to_vector();
        auto dv = (pt - bar_).to_vector();
        auto nbv = h_/2.;
        return point_type( bv.dot(dv)/(nbv*nbv) );
    }

public:
    scalar_monomial(const mesh_type& msh, const cell_type& cl, size_t degree)
        : degree_(degree), size_( size_of_degree(degree) )
    {
        bar_ = barycenter(msh, cl);
        auto pts = points(msh, cl);
        assert(pts.size() == 2);
        pa_ = pts[0];
        pb_ = pts[1];
        h_ = distance(pa_, pb_);
    }

    function_type operator()(const point_type& pt) {
        function_type ret = function_type::Zero(size_);
        auto ep = phys2ref(pt);

        ret(0) = 1.0;
        for (size_t k = 1; k <= degree_; k++)
            ret(k) = ret(k-1)*ep.x();
        
        return ret;
    }

    gradient_type grad(const point_type& pt) {
        gradient_type ret = gradient_type::Zero(size_);
        auto ep = phys2ref(pt);
        
        ret(0) = 0.;
        ret(1) = 2./h_;
        for (size_t k = 2; k <= degree_; k++)
            ret(k) = ret(k-1)*k*ep.x()/(k-1);

        return ret;
    }

    size_t size() const {
        return size_;
    }

    static size_t size_of_degree(size_t k) {
        return k+1;
    }

    size_t degree() const {
        return degree_;
    }

    size_t integration_degree() const {
        return degree_;
    }
};


template<template<typename, size_t, typename> class Mesh,
    typename CoordT, typename Storage, typename ScalT>
class scalar_monomial<Mesh<CoordT, 1, Storage>, typename Mesh<CoordT, 1, Storage>::face_type, ScalT>
{
public:
    using mesh_type = Mesh<CoordT, 1, Storage>;
    using cell_type = typename mesh_type::cell_type;
    using point_type = typename mesh_type::point_type;
    static const size_t dimension = 1;
    using coordinate_type = typename mesh_type::coordinate_type;
    using scalar_type = ScalT;
    using basis_category = scalar_basis;
    using function_type = scalar_type;
    using gradient_type = scalar_type;

private:
    size_t          degree_;
    size_t          size_;

public:
    scalar_monomial(const mesh_type& msh, const cell_type& cl, size_t degree)
        : degree_(degree), size_( size_of_degree(degree) )
    {}

    function_type operator()(const point_type& pt) {
        return 1.;
    }

    gradient_type grad(const point_type& pt) {
        return 0.;
    }

    size_t size() const {
        return size_;
    }

    static size_t size_of_degree(size_t k) {
        return 1;
    }

    size_t degree() const {
        return degree_;
    }

    size_t integration_degree() const {
        return 0;
    }
};


template<typename Mesh, typename Element, typename ScalT = typename Mesh::coordinate_type>
auto scaled_monomial_basis(const Mesh& msh, const Element& elem, size_t degree)
{
    return scalar_monomial<Mesh, Element, ScalT>(msh, elem, degree);
}




#if 0

template<typename Mesh, typename ScalT = typename Mesh::coordinate_type>
class scalar_monomial<Mesh, typename Mesh::cell_type, ScalT> {

    using coordinate_type = typename Mesh::coordinate_type;
    using scalar_type = ScalT;
    static const size_t DIM = Mesh::dimension;

    scalar_monomial_basis<coordinate_type, scalar_type, DIM, DIM> basis;

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

#endif


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