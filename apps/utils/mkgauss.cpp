#include <iostream>

#include "diskpp/quadratures/quad_bones.hpp"

using namespace disk;

int main(void)
{
    const int order_max = 40;
    
    for (int order = 0; order < order_max; order+=2)
    {
        auto qps = golub_welsch<double>(order);

        std::cout << std::endl;
        std::cout << "static gauss_point gauss_rule_" << 1+order/2 << "[] = {" << std::endl;
        for (size_t i = 0; i < qps.size(); i++) {
            auto& qp = qps[i];
            auto p = ( std::abs(qp.first.x()) < 1e-16 ) ? 0 : qp.first.x();
            auto w = qp.second;

            
            std::cout << "    { ";
            std::cout << std::setprecision(15) << std::setw(20) << p << ", ";
            std::cout << std::setprecision(15) << std::setw(20) << w;
            std::cout << " }" << (( i < qps.size()-1) ? "," : "") << std::endl;
        }
        std::cout << "};" << std::endl;
    }
}