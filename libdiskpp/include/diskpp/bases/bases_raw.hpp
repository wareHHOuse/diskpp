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

#include "diskpp/mesh/point.hpp"
#include "diskpp/common/eigen.hpp"

namespace disk {
namespace basis {

/* Raw scaled monomial scalar basis */
template<typename CoordT, typename ScalT, size_t IDIM, size_t DIM>
class scaled_monomial_basis
{
    static_assert(IDIM >= DIM, "Immersion space dimension cannot be smaller than basis space dimension");
};

/* Raw scaled monomial scalar basis: 1D specialization */
template<typename CoordT, typename ScalT, size_t IDIM>
class scaled_monomial_basis<CoordT, ScalT, IDIM, 1> {
    using scalar_type = ScalT;
    using rs_point_type = point<CoordT, 1>;     // reference space
    using is_point_type = point<CoordT, IDIM>;  // immersion space
    using function_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
    using gradient_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, IDIM>;

    is_point_type   pa_, pb_;
    is_point_type   bar_;
    CoordT          h_;
    size_t          degree_;
    size_t          size_;

    rs_point_type phys2ref(const is_point_type& pt) {
        auto bv = (pb_ - bar_).to_vector();
        auto dv = (pt - bar_).to_vector();
        auto nbv = h_/2.;
        return rs_point_type( bv.dot(dv)/(nbv*nbv) );
    }

public:
    scaled_monomial_basis(const is_point_type& pa, const is_point_type& pb, size_t degree)
        : pa_(pa), pb_(pb), bar_( midpoint(pa, pb) ), h_( distance(pa, pb) ),
          degree_(degree), size_( size_from_degree(degree) )
    {}

    function_type operator()(const is_point_type& pt) {
        function_type ret = function_type::Zero(size_);
        auto ep = phys2ref(pt);

        ret(0) = 1.0;
        for (size_t k = 1; k <= degree_; k++)
            ret(k) = ret(k-1)*ep.x();
        
        return ret;
    }

    gradient_type grad(const is_point_type& pt) {
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

    static size_t size_from_degree(size_t k) {
        return k+1;
    }

    size_t degree() const {
        return degree_;
    }
};

/* Raw scaled monomial scalar basis: 2D specialization */
template<typename CoordT, typename ScalT, size_t IDIM>
class scaled_monomial_basis<CoordT, ScalT, IDIM, 2> {
    using scalar_type = ScalT;
    using rs_point_type = point<CoordT, 2>;     // reference space
    using is_point_type = point<CoordT, IDIM>;  // immersion space
    using is_vector_type = Eigen::Matrix<CoordT, IDIM, 1>;
    using function_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
    using gradient_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, IDIM>;

    is_vector_type  v0_, v1_;
    is_point_type   bar_;
    size_t          degree_;
    CoordT          h2_;
    size_t          size_;
    Eigen::Matrix<CoordT, IDIM, 2> tr_;

    rs_point_type phys2ref(const is_point_type& pt) {
        auto v = (pt - bar_).to_vector();   /* origin is at (0)^IDIM */
        auto xi = v.dot(v0_);               /* v is in the plane spanned by v0_, v1_ */
        auto eta = v.dot(v1_);
        return rs_point_type(xi, eta);
    }

public:
    /* This expects an orthogonal basis s.t. |v0| = |v1| = h/2 centered on the element.
     * h is the _diameter_ of the element */
    scaled_monomial_basis(const is_vector_type& v0, const is_vector_type& v1,
        const is_point_type& bar, size_t degree)
        : v0_(v0), v1_(v1), bar_(bar), degree_(degree), h2_(v0.norm()),
          size_( size_from_degree(degree) )
    {
        /* orthogonal */
        //assert( v0.dot(v1) < std::numeric_limits<CoordT>::epsilon );
        /* same length */
        //assert( std::abs(v0.norm()-v1.norm()) < std::numeric_limits<CoordT>::epsilon );
        tr_.col(0) = v0;
        tr_.col(1) = v1;
        std::cout << tr_ << std::endl;
    }

    function_type operator()(const is_point_type& pt) {
        function_type powx = function_type::Zero(size_);
        function_type powy = function_type::Zero(size_);
        auto ep = phys2ref(pt);

        powx(0) = 1.0;
        powy(0) = 1.0;
        for (size_t k = 1; k <= degree_; k++) {
            powx(k) = powx(k-1) * ep.x();
            powy(k) = powy(k-1) * ep.y();
        }

        size_t pos = 0;
        function_type ret = function_type::Zero(size_);
        for (size_t k = 0; k <= degree_; k++)
            for (size_t i = 0; i <= k; i++)
                ret(pos++) = powx(k-i) * powy(i);
        assert(pos == size_);

        return ret;
    }

    gradient_type grad(const is_point_type& pt) {
        function_type powx = function_type::Zero(size_);
        function_type powy = function_type::Zero(size_);

        auto ep = phys2ref(pt);
        
        powx(0) = 1.0;
        powy(0) = 1.0;
        for (size_t k = 1; k <= degree_; k++) {
            powx(k) = powx(k-1) * ep.x();
            powy(k) = powy(k-1) * ep.y();
        }

        gradient_type ret = gradient_type::Zero(size_, IDIM);
        size_t pos = 0;
        for (size_t k = 0; k <= degree_; k++) {
            for (size_t i = 0; i <= k; i++) {
                const auto ex = k-i;
                const auto ey = i;
                const auto px = powx(ex);
                const auto py = powy(ey);
                const auto dx = (ex == 0) ? 0 : ex * powx(ex-1);
                const auto dy = (ey == 0) ? 0 : ey * powy(ey-1);
                Eigen::Matrix<ScalT, 2, 1> grad_ref;
                grad_ref(0) = dx*py;
                grad_ref(1) = px*dy;
                ret.block(pos,0,1,IDIM) = (tr_ * grad_ref).transpose();
                pos++;
            }
        }
        assert(pos == size_);

        return ret;
    }

    size_t size() const {
        return size_;
    }

    static size_t size_from_degree(size_t k) {
        return ((k+2)*(k+1))/2;
    }

    size_t degree() const {
        return degree_;
    }
};

/* Raw scaled monomial scalar basis: 3D specialization */
template<typename CoordT, typename ScalT, size_t IDIM>
class scaled_monomial_basis<CoordT, ScalT, IDIM, 3> {
    using scalar_type = ScalT;
    using rs_point_type = point<CoordT, 3>;     // reference space
    using is_point_type = point<CoordT, IDIM>;  // immersion space
    using function_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
    using gradient_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, IDIM>;

    size_t size_;
    size_t degree_;

    ScalT operator()(const point<CoordT,3>& pt) {}

        size_t size() const {
        return size_;
    }

    static size_t size_from_degree(size_t k) {
        return ((k+3)*(k+2)*(k+1))/6;
    }

    size_t degree() const {
        return degree_;
    }
};

} // namespace basis
} // namespace disk