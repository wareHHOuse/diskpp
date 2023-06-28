#include <iostream>
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/quadratures/quadratures.hpp"

int main(void)
{
    using T = double;
    disk::generic_mesh<T,2> msh;
    auto mesher = make_fvca5_hex_mesher(msh);
    mesher.make_level(0);

    auto ifun = [](const disk::point<T,2>& pt, int m, int n) {
        return std::pow(pt.x(), m)*std::pow(pt.y(), n);
    };

    for (size_t m = 0; m < 5; m++)
    {
        for (size_t n = 0; n < 5; n++)
        {
            double int_val = 0.0;
            for (auto& cl : msh)
            {
                auto tris = disk::triangulate_nonconvex_polygon(msh, cl);
                for (auto &tri : tris)
                {
                    auto qps = disk::quadrature::dunavant(m+n, tri.p0, tri.p1, tri.p2);
                    for (auto& qp : qps)
                        int_val += qp.weight() * ifun(qp.point(), m, n);
                }
            }

            double int_ana = 1.0/((m+1)*(n+1));

            double relerr = 100.0*std::abs(int_ana-int_val)/std::abs(int_ana);
            if (relerr > 1e-12) {
                std::cout << m << " " << n << ": " << relerr << std::endl;
                return 1;
            }
        }
    }

    return 0;
}