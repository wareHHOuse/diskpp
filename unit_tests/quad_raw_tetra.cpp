#include <iostream>
#include <iomanip>

#include "diskpp/quadratures/bits/quad_raw_tetra.hpp"
#include "sgr.hpp"

using namespace disk;
using namespace disk::quadrature;
using namespace sgr;

#define THRESH 2e-15

using T = double;

template<typename Function>
bool test_tetra_quadrature(Function& quad, int max_order)
{
    T x0 = 0.5, x1 = 1.0;
    T y0 = 0.5, y1 = 1.0;
    T z0 = 0.5, z1 = 1.0;

    auto int_ana_fun = [&](int m, int n, int p) {
        auto x1m1 = std::pow(x1, m+1);
        auto x0m1 = std::pow(x0, m+1);
        auto y1n1 = std::pow(y1, n+1);
        auto y0n1 = std::pow(y0, n+1);
        auto z1p1 = std::pow(z1, p+1);
        auto z0p1 = std::pow(z0, p+1);
        auto a = (x1m1 - x0m1)/(m+1);
        auto b = (y1n1 - y0n1)/(n+1);
        auto c = (z1p1 - z0p1)/(p+1);
        return a*b*c;
    };

    point<T,3> p0(x0,y0,z0);
    point<T,3> p1(x1,y0,z0);
    point<T,3> p2(x1,y1,z0);
    point<T,3> p3(x0,y1,z0);
    point<T,3> p4(x0,y0,z1);
    point<T,3> p5(x1,y0,z1);
    point<T,3> p6(x1,y1,z1);
    point<T,3> p7(x0,y1,z1);

    auto ifun = [](const point<T,3>& pt, int m, int n, int p) {
        return std::pow(pt.x(), m)*std::pow(pt.y(), n)*std::pow(pt.z(), p);
    };

    bool fail = false;

    for (int m = 0; m <= max_order; m++)
    {
        for (int n = 0; n <= max_order; n++)
        {
            for (int p = 0; p <= max_order; p++)
            {
                if (m+n+p > max_order)
                    continue;
                
                auto order = m+n+p;
                T int_num = 0.0;
                auto qps1 = quad(order, p0, p1, p3, p7);
                for (auto& qp : qps1)
                    int_num += qp.weight() * ifun(qp.point(), m, n, p);
                auto qps2 = quad(order, p0, p1, p4, p7);
                for (auto& qp : qps2)
                    int_num += qp.weight() * ifun(qp.point(), m, n, p);
                auto qps3 = quad(order, p1, p2, p3, p7);
                for (auto& qp : qps3)
                    int_num += qp.weight() * ifun(qp.point(), m, n, p);
                auto qps4 = quad(order, p1, p2, p6, p7);
                for (auto& qp : qps4)
                    int_num += qp.weight() * ifun(qp.point(), m, n, p);
                auto qps5 = quad(order, p1, p4, p5, p7);
                for (auto& qp : qps5)
                    int_num += qp.weight() * ifun(qp.point(), m, n, p);
                auto qps6 = quad(order, p1, p5, p6, p7);
                for (auto& qp : qps6)
                    int_num += qp.weight() * ifun(qp.point(), m, n, p);
                
                auto int_ana = int_ana_fun(m,n,p);
                auto err = std::abs(int_num - int_ana);
                if (err > THRESH) {
                    std::cout << redfg << "[FAIL] ";
                    std::cout << std::setw(3) << m << " ";
                    std::cout << std::setw(3) << n << " ";
                    std::cout << std::setw(3) << p << " ";
                    std::cout << std::setw(15) << err << "     ";
                    std::cout << std::setw(10) << int_num << "    ";
                    std::cout << std::setw(10) << int_ana << nofg << std::endl;
                    fail = true;
                }
            }
        }
    }

    return fail;
}

int main(void)
{
    bool fail = false;

    using pt_t = const point<T,3>&;

    auto q_arbq = [](size_t degree, pt_t p0, pt_t p1, pt_t p2, pt_t p3) {
        return arbq(degree, p0, p1, p2, p3);
    };

    fail |= test_tetra_quadrature(q_arbq, 15);

    auto q_gm = [](size_t degree, pt_t p0, pt_t p1, pt_t p2, pt_t p3) {
        return grundmann_moeller(degree, p0, p1, p2, p3);
    };

    fail |= test_tetra_quadrature(q_gm, 15);

    return fail;
}