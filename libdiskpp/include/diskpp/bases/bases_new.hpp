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

template<typename Mesh, typename Element, typename ScalT>
class scalar_monomial;

/**********************************************************************************************
 * Cells, 1D
 */
template<template<typename, size_t, typename> class Mesh,
    typename CoordT, typename Storage, typename ScalT>
class scalar_monomial<Mesh<CoordT, 1, Storage>, typename Mesh<CoordT, 1, Storage>::cell_type, ScalT>
{
public:
    using mesh_type = Mesh<CoordT, 1, Storage>;
    using cell_type = typename mesh_type::cell_type;
    using element_type = cell_type;
    using coordinate_type = typename mesh_type::coordinate_type;
    using scalar_type = ScalT;
    using point_type = typename mesh_type::point_type;
    static const size_t immersion_dimension = mesh_type::dimension;
    static const size_t basis_dimension = mesh_type::dimension;
    static const size_t tensor_order = 0; /* scalar */
    using value_type = typename tensor<scalar_type, basis_dimension, tensor_order>::value_type;
    using value_array_type = typename tensor<scalar_type, basis_dimension, tensor_order>::array_type;
    using gradient_array_type = typename tensor<scalar_type, basis_dimension, tensor_order+1>::array_type;

private:
    point_type      pa_, pb_;
    point_type      bar_;
    coordinate_type h_;
    size_t          degree_;
    size_t          size_;

    point_type phys2ref(const point_type& pt) const {
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

    value_array_type operator()(const point_type& pt) const {
        value_array_type ret = value_array_type::Zero(size_);
        auto ep = phys2ref(pt);

        ret(0) = 1.0;
        for (size_t k = 1; k <= degree_; k++)
            ret(k) = ret(k-1)*ep.x();
        
        return ret;
    }

    gradient_array_type grad(const point_type& pt) const {
        gradient_array_type ret = gradient_array_type::Zero(size_);
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

/**********************************************************************************************
 * Faces, 1D
 */
template<template<typename, size_t, typename> class Mesh,
    typename CoordT, typename Storage, typename ScalT>
class scalar_monomial<Mesh<CoordT, 1, Storage>, typename Mesh<CoordT, 1, Storage>::face_type, ScalT>
{
public:
    using mesh_type = Mesh<CoordT, 1, Storage>;
    using face_type = typename mesh_type::face_type;
    using element_type = face_type;
    using coordinate_type = typename mesh_type::coordinate_type;
    using scalar_type = ScalT;
    using point_type = typename mesh_type::point_type;
    static const size_t immersion_dimension = mesh_type::dimension;
    static const size_t basis_dimension = mesh_type::dimension - 1;
    static const size_t tensor_order = 0; /* scalar */
    using value_type = typename tensor<scalar_type, basis_dimension, tensor_order>::value_type;
    using value_array_type = typename tensor<scalar_type, basis_dimension, tensor_order>::array_type;
    using gradient_array_type = typename tensor<scalar_type, basis_dimension, tensor_order+1>::array_type;

private:
    size_t          degree_;
    size_t          size_;

public:
    scalar_monomial(const mesh_type& msh, const face_type& fc, size_t degree)
        : degree_(degree), size_( size_of_degree(degree) )
    {}

    value_array_type operator()(const point_type& pt) const {
        return 1.;
    }

    gradient_array_type grad(const point_type& pt) const {
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


/**********************************************************************************************
 * Cells, 2D
 */
template<template<typename, size_t, typename> class Mesh,
    typename CoordT, typename Storage, typename ScalT>
class scalar_monomial<Mesh<CoordT, 2, Storage>, typename Mesh<CoordT, 2, Storage>::cell_type, ScalT>
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
    static const size_t tensor_order = 0; /* scalar */
    using value_type = typename tensor<scalar_type, basis_dimension, tensor_order>::value_type;
    using value_array_type = typename tensor<scalar_type, basis_dimension, tensor_order>::array_type;
    using gradient_array_type = typename tensor<scalar_type, basis_dimension, tensor_order+1>::array_type;

private:
    using tr_mat_type = Eigen::Matrix<coordinate_type, immersion_dimension, immersion_dimension>;
    using vec_type = Eigen::Matrix<coordinate_type, immersion_dimension, 1>;

    point_type      bar_;
    size_t          degree_;
    size_t          size_;
    tr_mat_type     tr_;
    vec_type        v0_, v1_;

    void
    compute_reference_frame(const mesh_type& msh, const cell_type& cl)
    {
        auto pts = points(msh, cl);
        assert(pts.size() >= 3);
        v0_ = (pts[1] - pts[0]).to_vector();
        v1_ = (pts[2] - pts[1]).to_vector();
        v1_ = v1_ - v1_.dot(v0_)*v0_/(v0_.dot(v0_)); // Gram-Schmidt
        tr_.col(0) = v0_;
        tr_.col(1) = v1_;
        tr_ = tr_.inverse().eval();
    }

    point_type phys2ref(const point_type& pt) const {
        auto v = tr_*(pt - bar_).to_vector();
        return point_type(v(0), v(1));
    }

public:
    scalar_monomial(const mesh_type& msh, const cell_type& cl, size_t degree)
        : degree_(degree), size_( size_of_degree(degree) )
    {
        bar_ = barycenter(msh, cl);
        compute_reference_frame(msh, cl);
    }

    value_array_type operator()(const point_type& pt) const {
        value_array_type powx = value_array_type::Zero(size_);
        value_array_type powy = value_array_type::Zero(size_);
        auto ep = phys2ref(pt);

        powx(0) = 1.0;
        powy(0) = 1.0;
        for (size_t k = 1; k <= degree_; k++) {
            powx(k) = powx(k-1) * ep.x();
            powy(k) = powy(k-1) * ep.y();
        }

        size_t pos = 0;
        value_array_type ret = value_array_type::Zero(size_);
        for (size_t k = 0; k <= degree_; k++)
            for (size_t i = 0; i <= k; i++)
                ret(pos++) = powx(k-i) * powy(i);
        assert(pos == size_);

        return ret;
    }

    gradient_array_type grad(const point_type& pt) const {
        value_array_type powx = value_array_type::Zero(size_);
        value_array_type powy = value_array_type::Zero(size_);

        auto ep = phys2ref(pt);
        
        powx(0) = 1.0;
        powy(0) = 1.0;
        for (size_t k = 1; k <= degree_; k++) {
            powx(k) = powx(k-1) * ep.x();
            powy(k) = powy(k-1) * ep.y();
        }

        gradient_array_type ret = gradient_array_type::Zero(size_, immersion_dimension);
        size_t pos = 0;
        for (size_t k = 0; k <= degree_; k++) {
            for (size_t i = 0; i <= k; i++) {
                const auto ex = k-i;
                const auto ey = i;
                const auto px = powx(ex);
                const auto py = powy(ey);
                const auto dx = (ex == 0) ? 0 : ex * powx(ex-1);
                const auto dy = (ey == 0) ? 0 : ey * powy(ey-1);
                vec_type tmp_grad;
                tmp_grad(0) = dx*py;
                tmp_grad(1) = px*dy;
                ret.block(pos, 0, 1, immersion_dimension) = tmp_grad.transpose()*tr_;
                pos++;
            }
        }
        assert(pos == size_);

        return ret;
    }

    size_t size() const {
        return size_;
    }

    static size_t size_of_degree(size_t k) {
        return ((k+2)*(k+1))/2;
    }

    size_t degree() const {
        return degree_;
    }

    size_t integration_degree() const {
        return degree_;
    }
};


/**********************************************************************************************
 * Faces, 2D
 */
template<template<typename, size_t, typename> class Mesh,
    typename CoordT, typename Storage, typename ScalT>
class scalar_monomial<Mesh<CoordT, 2, Storage>, typename Mesh<CoordT, 2, Storage>::face_type, ScalT>
{
public:
    using mesh_type = Mesh<CoordT, 2, Storage>;
    using face_type = typename mesh_type::face_type;
    using element_type = face_type;
    using coordinate_type = typename mesh_type::coordinate_type;
    using scalar_type = ScalT;
    using point_type = typename mesh_type::point_type;
    static const size_t immersion_dimension = mesh_type::dimension;
    static const size_t basis_dimension = mesh_type::dimension - 1;
    static const size_t tensor_order = 0; /* scalar */
    using value_type = typename tensor<scalar_type, basis_dimension, tensor_order>::value_type;
    using value_array_type = typename tensor<scalar_type, basis_dimension, tensor_order>::array_type;
    using gradient_array_type = typename tensor<scalar_type, basis_dimension, tensor_order+1>::array_type;

private:
    point_type  pa_, pb_;
    point_type  bar_;
    CoordT      h_;
    size_t      degree_;
    size_t      size_;

    coordinate_type phys2ref(const point_type& pt) {
        auto bv = (pb_ - bar_).to_vector();
        auto dv = (pt - bar_).to_vector();
        auto nbv = h_/2.;
        return rs_point_type( bv.dot(dv)/(nbv*nbv) );
    }

public:
    scalar_monomial(const mesh_type& msh, const face_type& fc, size_t degree)
        : degree_(degree), size_( size_of_degree(degree) )
    {
        auto pts = points(msh, fc);
        assert(pts.size() == 2);
        pa_ = pts[0];
        pb_ = pts[1];
        bar_ = (pa_ + pb_)/2.;
        h_ = distance(pa_, pb_);
    }

    value_array_type operator()(const point_type& pt) const {
        value_array_type ret = value_array_type::Zero(size_);
        auto ep = phys2ref(pt);

        ret(0) = 1.0;
        for (size_t k = 1; k <= degree_; k++)
            ret(k) = ret(k-1)*ep;
        
        return ret;
    }

    gradient_array_type grad(const point_type& pt) const {
        throw std::logic_error("not implemented yet");
        value_array_type ret = value_array_type::Zero(size_, immersion_dimension);
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
        return 0;
    }
};




template<typename Mesh, typename Element, typename ScalT = typename Mesh::coordinate_type>
auto scaled_monomial_basis(const Mesh& msh, const Element& elem, size_t degree)
{
    return scalar_monomial<Mesh, Element, ScalT>(msh, elem, degree);
}

template<typename Basis>
class L2_projector
{
    using mesh_type = typename Basis::mesh_type;
    using element_type = typename Basis::element_type;
    using coordinate_type = typename Basis::coordinate_type;
    using scalar_type = typename Basis::scalar_type;
    using rhs_func_type = std::function<typename Basis::value_type(typename Basis::point_type)>;
    static const size_t immersion_dimension = Basis::immersion_dimension;
    
    const Basis& basis_;

    using matrix_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;

    matrix_type mass_;
    std::vector<quadrature_point<coordinate_type, immersion_dimension>> qps_;

    Eigen::LDLT<matrix_type> ldlt_;

public:
    L2_projector(const mesh_type& msh, const element_type& elem, const Basis& pb)
        : basis_(pb)
    {
        auto degree = basis_.integration_degree()*2+1;
        mass_ = matrix_type::Zero(basis_.size(), basis_.size());
        qps_ = integrate(msh, elem, degree);
        for (auto& qp : qps_) {
            auto phi = basis_(qp.point());
            mass_ += qp.weight() * phi * phi.transpose();
        }
        ldlt_.compute(mass_);
    }

    vector_type operator()(const rhs_func_type& f) const {
        auto degree = basis_.integration_degree()*2+1;
        vector_type rhs = vector_type::Zero(basis_.size());
        for (auto& qp : qps_) {
            auto phi = basis_(qp.point());
            /* This must be fixed, for now it works only with
             * scalar-valued and vector-valued functions */
            rhs += qp.weight() * phi * f(qp.point());
        }

        vector_type dofs = ldlt_.solve(rhs);
        return dofs;
    }

    template<typename OtherBasis>
    typename std::enable_if< can_take_scalar_product<Basis,OtherBasis>::value, matrix_type >::type
    basis(const OtherBasis& other) const
    {
        auto degree = basis_.integration_degree() + other.integration_degree();
        matrix_type rhs = matrix_type::Zero(basis_.size(), other.size());
        for (auto& qp : qps_) {
            auto test = basis_(qp.point());
            auto trial = other(qp.point());
            /* This must be fixed, for now it works only with
             * scalar-valued and vector-valued functions */
            rhs += qp.weight() * test * trial.transpose();
        }

        matrix_type oper = ldlt_.solve(rhs);
        return oper;
    }
};

template<typename Mesh, typename Element, typename IntegrandFunction>
auto integrate(const Mesh& msh, const Element& elem, size_t degree, IntegrandFunction f)
{
    using return_type = decltype( f(typename Mesh::point_type{}) );
    
    return_type ret{};

    auto qps = integrate(msh, elem, degree);
    for (auto& qp : qps) {
        ret += qp.weight() * f(qp.point());
    };

    return ret;
}

template<typename Trial, typename Test>
class L2_scalar_product
{
    static_assert(can_take_scalar_product<Trial, Test>::value);
    static_assert(Trial::immersion_dimension == Test::immersion_dimension);
    static_assert(std::is_same<typename Trial::coordinate_type, typename Test::coordinate_type>::value);
    using coordinate_type = typename Test::coordinate_type;
    static const size_t immersion_dimension = Test::immersion_dimension;
    using scalar_type = decltype(typename Trial::scalar_value{} * typename Test::scalar_value{});
    
    std::vector<quadrature_point<coordinate_type, immersion_dimension>> qps_;

    const Trial& trial_;
    const Test& test_;

    using matrix_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;

    matrix_type mat_;

    template<typename Mesh, typename Element>
    L2_scalar_product(const Mesh& msh, const Element& elem, const Trial& trial, const Test& test)
        : trial_(trial), test_(test)
    {
        auto degree = trial_.integration_degree() + test_.integration_degree();
        mat_ = matrix_type::Zero(test_.size(), trial_.size());
        qps_ = integrate(msh, elem, degree);
        for (auto& qp : qps_) {
            /* This must be fixed, for now it works only with
             * scalar-valued and vector-valued functions */
            mat_ += qp.weight() * test_(qp.point()) * trial(qp.point()).transpose();
        }
    }
};


} // namespace disk::basis