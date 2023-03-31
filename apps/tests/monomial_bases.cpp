/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

/* Test for 1D raw monomial bases: interpolate a sine, measure the error
 * on the interpolation and its gradient. */

#include "diskpp/bases/bases_raw.hpp"
#include "diskpp/quadratures/quad_raw_gauss.hpp"

using namespace disk;
using namespace disk::quadrature;
using namespace disk::basis;

int main(void)
{
    using T = double;
    size_t degree = 8;

    point<T,1> center({0});
    T a = 0;
    T b = M_PI;
    T h = b-a;

    scaled_monomial_basis<T, T, 1> basis(center, h, degree);

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mass =
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(basis.size(), basis.size());

    Eigen::Matrix<T, Eigen::Dynamic, 1> rhs =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(basis.size());

    auto qps = gauss_legendre(2*degree+2, a, b);
    for (auto& qp : qps) {
        auto x = qp.point().x();
        auto phi = basis(x);
        mass += qp.weight() * phi * phi.transpose();
        rhs += qp.weight() * sin(x) * phi;
    }

    Eigen::Matrix<T, Eigen::Dynamic, 1> proj = mass.ldlt().solve(rhs);
    
    T fun_error = 0.0;
    T grad_error = 0.0;
    for (auto& qp : qps) {
        auto x = qp.point().x();
        auto phi = basis(x);
        auto dphi = basis.grad(x);
        auto diff_fun = sin(x) - proj.dot(phi);
        auto diff_grad = cos(x) - proj.dot(dphi);
        fun_error += qp.weight() * diff_fun * diff_fun;
        grad_error += qp.weight() * diff_grad * diff_grad;
    }

    auto fe = std::sqrt(fun_error);
    auto ge = std::sqrt(grad_error);

    const double DBL_FE_MAX_ERR = 7.57e-11;
    const double DBL_GE_MAX_ERR = 1.56e-06;

    if (fe > DBL_FE_MAX_ERR or ge > DBL_GE_MAX_ERR) {
        std::cout << __FILE__ << ": Test FAILED: fe = " << fe << ", ge = " << ge << std::endl;
        return 1; 
    }

    return 0;
}