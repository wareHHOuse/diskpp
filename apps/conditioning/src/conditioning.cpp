#include <iostream>

#include "diskpp/mesh/point.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/bases/bases_new.hpp"
#include "diskpp/bases/bases_operations.hpp"

#include "sgr.hpp"

using namespace sgr;
using namespace disk::basis;

template<typename T>
auto cond(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A)
{
    using MT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    Eigen::JacobiSVD<MT> svd(A);
    auto lmax = svd.singularValues()(0);
    auto lmin = svd.singularValues()(svd.singularValues().size()-1);
    return lmax/lmin; 
}

int main(void)
{
    using T = double;
    using mesh_type = disk::generic_mesh<T,1>;
    using cell_type = typename mesh_type::cell_type;
    using face_type = typename mesh_type::face_type;
    using ScalT = double;
    using cell_monomial_basis = scalar_monomial<mesh_type, cell_type, ScalT>;
    

    disk::point<T,1>  p0{0.0};
    disk::point<T,1>  p1{1.0};

    mesh_type msh;
    make_single_element_mesh(msh, p0.x(), p1.x());

    std::cout << Bgreenfg << "Scale factor = 2.0" << reset << std::endl;
    for (size_t degree = 1; degree < 6; degree++)
    {
        std::cout << Bon << "  Degree " << degree << Boff << std::endl; 
        for (auto& cl : msh)
        {
            cell_monomial_basis phi(msh, cl, degree);

            std::cout << "  phi(-1):  " << phi(p0).transpose() << std::endl;
            std::cout << "  grad(-1): " << grad(phi)(p0).transpose() << std::endl;

            auto mass = integrate(msh, cl, phi, phi);
            std::cout << yellowfg << "    Mass cond: " << cond(mass) << ", ";

            auto stiff = integrate(msh, cl, grad(phi), grad(phi));
            auto n = phi.size() - 1;
            auto stiff_nc = stiff.bottomRightCorner(n,n).eval();
            std::cout << "stiff cond: " << cond(stiff_nc) << nofg << std::endl;
        }
    }

    std::cout << Bgreenfg << "Scale factor = 1.0" << reset << std::endl;
    for (size_t degree = 1; degree < 6; degree++)
    {
        std::cout << Bon << "  Degree " << degree << Boff << std::endl;  
        for (auto& cl : msh)
        {
            cell_monomial_basis phi(msh, cl, degree);
            phi.scalefactor(1.0);

            std::cout << "  phi(-1):  " << phi(p0).transpose() << std::endl;
            std::cout << "  grad(-1): " << grad(phi)(p0).transpose() << std::endl;

            auto mass = integrate(msh, cl, phi, phi);
            std::cout << yellowfg << "    Mass cond: " << cond(mass) << ", ";

            auto stiff = integrate(msh, cl, grad(phi), grad(phi));
            auto n = phi.size() - 1;
            auto stiff_nc = stiff.bottomRightCorner(n,n).eval();
            std::cout << "stiff cond: " << cond(stiff_nc) << nofg << std::endl;
        }
    }

    return 0;

}