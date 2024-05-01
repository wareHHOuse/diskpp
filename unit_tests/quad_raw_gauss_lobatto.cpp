#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>

#include "diskpp/quadratures/bits/quad_raw_gauss_lobatto.hpp"
#include "sgr.hpp"

int main(void)
{
    size_t deg_min = 0;
    size_t deg_max = 35;

    double a = -1;
    double b = 2;

    int fails = 0;

    for (size_t k = deg_min; k <= deg_max; k++) {

        auto qps = disk::quadrature::gauss_lobatto(k, a, b);
        double num_int = 0.0;
        for (auto& qp : qps)
            num_int += qp.weight() * std::pow(qp.point().x(), k);

        double ana_int = (std::pow(b, k+1) - std::pow(a, k+1))/(k+1);

        double abserr = std::abs(num_int - ana_int);
        double relerr = abserr / std::abs(ana_int);

        printf("degree = %lu, num = %g, ana = %g, relerr = %g\n",
            k, num_int, ana_int, relerr);
    
        if (relerr > 1e-14)
            fails++;
    }

    return fails;
}