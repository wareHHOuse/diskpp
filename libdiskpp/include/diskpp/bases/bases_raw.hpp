#pragma once

#include "diskpp/mesh/point.hpp"
#include "diskpp/common/eigen.hpp"

namespace disk {
namespace basis {

/* Raw scaled monomial scalar basis */
template<typename CoordT, typename ScalT, size_t DIM>
class scaled_monomial_basis;

/* Raw scaled monomial scalar basis: 1D specialization */
template<typename CoordT, typename ScalT>
class scaled_monomial_basis<CoordT, ScalT, 1> {
    using scalar_type = ScalT;
    using point_type = point<CoordT,1>;

    using function_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
    using gradient_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;

    point_type  center_;
    CoordT      h_;
    size_t      degree_;
    size_t      size_;

    size_t size_from_degree(size_t k) {
        return k+1;
    }

public:
    scaled_monomial_basis()
        : center_({0.0}), h_(2.0), degree_(1), size_(2)
    {}
    scaled_monomial_basis(const point_type& center, const CoordT& h, size_t degree)
        : center_(center), h_(h), degree_(degree), size_( size_from_degree(degree) )
    {}

    function_type operator()(const point_type& pt) {
        function_type ret = function_type::Zero(size_);
        auto ep = 0.5*(pt - center_)/h_;

        ret(0) = 1.0;
        for (size_t k = 1; k <= degree_; k++)
            ret(k) = ret(k-1)*ep.x();
        
        return ret;
    }

    gradient_type grad(const point_type& pt) {
        gradient_type ret = gradient_type::Zero(size_);
        auto ep = 0.5*(pt - center_)/h_;
        
        ret(0) = 0;
        ret(1) = 0.5/h_;
        for (size_t k = 2; k <= degree_; k++)
            ret(k) = ret(k-1)*k*ep.x()/(k-1);

        return ret;
    }

    size_t size() const {
        return size_;
    }

    size_t degree() const {
        return degree_;
    }
};

/* Raw scaled monomial scalar basis: 2D specialization */
template<typename CoordT, typename ScalT>
class scaled_monomial_basis<CoordT, ScalT, 2> {

    scaled_monomial_basis() {}
    scaled_monomial_basis(const CoordT& ph) {}
    scaled_monomial_basis(const CoordT& phx, const CoordT& phy) {}

    ScalT operator()(const point<CoordT,2>& pt) {}
};

/* Raw scaled monomial scalar basis: 3D specialization */
template<typename CoordT, typename ScalT>
class scaled_monomial_basis<CoordT, ScalT, 3> {

    scaled_monomial_basis() {}
    scaled_monomial_basis(const CoordT& ph) {}
    scaled_monomial_basis(const CoordT& phx, const CoordT& phy, const CoordT& phz) {}

    ScalT operator()(const point<CoordT,3>& pt) {}
};

} // namespace basis
} // namespace disk