#include <iostream>
#include <iomanip>

#include "diskpp/quadratures/quad_raw_triangle_gauss.hpp"
#include "diskpp/quadratures/quad_raw_dunavant.hpp"
#include "sgr.hpp"

using namespace disk;
using namespace disk::quadrature;
using namespace sgr;

int main(void)
{
    using T = double;

    T x0 = 0.5, x1 = 1.0;
    T y0 = 0.5, y1 = 1.0;

    auto int_ana_fun = [&](int m, int n) {
        auto x1m1 = std::pow(x1, m+1);
        auto x0m1 = std::pow(x0, m+1);
        auto y1n1 = std::pow(y1, n+1);
        auto y0n1 = std::pow(y0, n+1);
        auto a = (x1m1 - x0m1)/(m+1);
        auto b = (y1n1 - y0n1)/(n+1);
        return a*b;
    };

    point<T,2> p0(x0,y0);
    point<T,2> p1(x1,y0);
    point<T,2> p2(x0,y1);
    point<T,2> p3(x1,y1);

    auto ifun = [](const point<T,2>& pt, int m, int n) {
        return std::pow(pt.x(), m)*std::pow(pt.y(), n);
    };

    const int max_order = QUAD_RAW_TRIANGLE_GAUSS_MAX_ORDER;
    for (int m = 0; m <= max_order; m++)
    {
        for (int n = 0; n <= max_order; n++)
        {
            if (m+n > max_order)
                continue;
            
            auto order = m+n;
            T int_num = 0.0;
            auto qps1 = triangle_gauss(order, p0, p1, p2);
            for (auto& qp : qps1)
                int_num += qp.weight() * ifun(qp.point(), m, n);
            auto qps2 = triangle_gauss(order, p1, p2, p3);
            for (auto& qp : qps2)
                int_num += qp.weight() * ifun(qp.point(), m, n);
            
            auto int_ana = int_ana_fun(m,n);
            auto relerr = 100*std::abs(int_num - int_ana)/std::abs(int_ana);
            if (relerr > 1e-11)
                std::cout << redfg;
            else
                std::cout << greenfg;
            
            std::cout << std::setw(3) << m << " ";
            std::cout << std::setw(3) << n << "    ";
            std::cout << std::setw(15) << relerr << "%    ";
            std::cout << std::setw(10) << int_num << "    ";
            std::cout << std::setw(10) << int_ana << nofg << std::endl;
        }
    }

    return 0;
}