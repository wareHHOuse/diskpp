#include <iostream>

#include "diskpp/mesh/point.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/bases/bases_new.hpp"
#include "diskpp/bases/bases_operations.hpp"
#include "diskpp/bases/bases.hpp"

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

template<typename Mesh>
void test_conditioning(const Mesh& msh, double scalefactor,
    const typename Mesh::point_type& tp, rescaling_strategy rs)
{
    using mesh_type = Mesh;
    using cell_type = typename mesh_type::cell_type;
    using face_type = typename mesh_type::face_type;
    using ScalT = double;
    using cell_monomial_basis = scalar_monomial<mesh_type, cell_type, ScalT>;
    
    std::cout << Bgreenfg << "Scale factor = " << scalefactor << reset << std::endl;
    for (size_t degree = 1; degree < 3; degree++)
    {
        std::cout << Bon << "  Degree " << degree << Boff << std::endl;  
        for (auto& cl : msh)
        {
            cell_monomial_basis phi(msh, cl, degree, rs);
            phi.scalefactor(scalefactor);

            //std::cout << "  phi:  " << phi(tp).transpose() << std::endl;
            //std::cout << "  grad: " << grad(phi)(tp).transpose() << std::endl;

            auto mass = integrate(msh, cl, phi, phi);
            std::cout << yellowfg << "    Mass cond: " << cond(mass) << ", ";

            auto stiff = integrate(msh, cl, grad(phi), grad(phi));
            auto n = phi.size() - 1;
            auto stiff_nc = stiff.bottomRightCorner(n,n).eval();
            std::cout << "stiff cond: " << cond(stiff_nc) << nofg << std::endl;

            std::cout << mass << std::endl;
            std::cout << "  Diameter " << diameter(msh, cl) << std::endl;
        }
    }
}

/*
template<typename Mesh>
void cmm(const Mesh& msh)
{
    auto degree = 2;
    for (auto& cl : msh) {
        const auto cb = make_scalar_monomial_basis(msh, cl, degree);
        disk::dynamic_matrix<double> H2 =
            disk::dynamic_matrix<double>::Zero(cb.size(), cb.size());
        const auto qps = integrate(msh, cl, 2 * degree);
        for (auto& qp : qps) {
            const auto phi = cb.eval_functions(qp.point());
            H2 += qp.weight() * phi * phi.transpose();
        }
        std::cout << "****** H2 ******" << std::endl;
        std::cout << H2 << std::endl;
        std::cout << yellowfg << "  H2 cond: " << cond(H2) << std::endl;
        std::cout << "  Diameter " << diameter(msh, cl) << std::endl;
    }
}
*/

int main(int argc, char **argv)
{
    using T = double;

    disk::generic_mesh<T,1> msh_1D;
    make_single_element_mesh(msh_1D, 0.0, 1.0);

    test_conditioning(msh_1D, 2.0, {1.0}, rescaling_strategy::none);
    test_conditioning(msh_1D, 1.0, {1.0}, rescaling_strategy::none);

    disk::generic_mesh<T,2> msh_2D;
    make_single_element_mesh(msh_2D, 1.0, 5);
    auto tr2D = [](const disk::point<T,2>& pt) -> disk::point<T,2> {
        return {pt.x(), 0.1*pt.y()};
    };
    msh_2D.transform(tr2D);

    std::cout << " ----------- NONE " << std::endl;
    test_conditioning(msh_2D, 2.0, {1.0, 0.0}, rescaling_strategy::none);
    std::cout << " ----------- INERTIAL " << std::endl;
    test_conditioning(msh_2D, 2.0, {1.0, 0.0}, rescaling_strategy::inertial);
    std::cout << " ----------- G-S " << std::endl;
    test_conditioning(msh_2D, 2.0, {1.0, 0.0}, rescaling_strategy::gram_schmidt);

    if (argc > 1) {
        disk::generic_mesh<T,2> msh_2D_fromfile;
        load_single_element_csv(msh_2D_fromfile, argv[1]);
        std::cout << " ----------- NONE " << std::endl;
        test_conditioning(msh_2D_fromfile, 2.0, {1.0, 0.0}, rescaling_strategy::none);
        std::cout << " ----------- INERTIAL " << std::endl;
        test_conditioning(msh_2D_fromfile, 2.0, {1.0, 0.0}, rescaling_strategy::inertial);
    }

    return 0;
}