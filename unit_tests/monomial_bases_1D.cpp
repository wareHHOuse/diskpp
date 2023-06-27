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

#include <iostream>
#include <fstream>

#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/bases/bases_new.hpp"
#include "diskpp/bases/bases_operations.hpp"
#include "diskpp/quadratures/quadratures_raw.hpp"

using namespace disk;
using namespace disk::quadrature;
using namespace disk::basis;

int main(void)
{
    using T = double;
    size_t degree = 8;

    T a = 0.0;
    T b = M_PI;

    generic_mesh<T,1> msh;
    make_single_element_mesh(msh, a, b);

    auto cl = *msh.cells_begin();

    auto basis = scaled_monomial_basis(msh, cl, degree);

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mass =
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(basis.size(), basis.size());

    Eigen::Matrix<T, Eigen::Dynamic, 1> rhs =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(basis.size());

    auto qps = gauss_legendre(2*degree+2, a, b);
    for (auto& qp : qps) {
        T x = qp.point().x();
        auto phi = basis(x);
        mass += qp.weight() * phi * phi.transpose();
        rhs += qp.weight() * std::sin(x) * phi;
    }

    auto sol_deriv = [](const point<T,1>& pt) {
        Eigen::Matrix<T, 1, 1> ret;
        ret(0,0) = cos(pt.x());
        return ret;
    };

    Eigen::Matrix<T, Eigen::Dynamic, 1> dofs = mass.ldlt().solve(rhs);

    T fun_error = 0.0, fun_norm = 0.0;
    T grad_error = 0.0, grad_norm = 0.0;
    for (auto& qp : qps) {
        auto fnum = eval(qp.point(), dofs, basis);
        auto fval = std::sin(qp.point().x());
        auto diff_fun = fval - fnum;
        fun_error += qp.weight() * diff_fun * diff_fun;
        fun_norm += qp.weight() * fval * fval;

        auto gnum = eval(qp.point(), dofs, grad(basis));
        auto gval = sol_deriv(qp.point());
        auto diff_grad = gval - gnum;
        grad_error += qp.weight() * diff_grad * diff_grad;
        grad_norm += qp.weight() * gval.dot(gval);
    }

    auto fe = 100.0*std::sqrt(fun_error/fun_norm);
    auto ge = 100.0*std::sqrt(grad_error/grad_norm);

    bool fepass = false, gepass = false;
    if (fe < 1.96e-12) fepass = true;
    if (ge < 8.79e-5) gepass = true;

    auto passfail = [](bool pass) {
        return pass ? "[PASS]" : "[FAIL]";
    };

    std::cout << __FILE__ << std::endl;
    std::cout << "  Function relative error = " << fe << "% " << passfail(fepass) << std::endl;
    std::cout << "  Gradient relative error = " << ge << "% " << passfail(gepass) << std::endl;

    return not (fepass and gepass);
}