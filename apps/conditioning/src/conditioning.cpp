#include <iostream>

#include "diskpp/mesh/point.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/bases/bases_new.hpp"
#include "diskpp/bases/bases_operations.hpp"
#include "diskpp/bases/bases.hpp"
#include "diskpp/output/silo.hpp"
#include "diskpp/common/timecounter.hpp"

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
    for (size_t degree = 1; degree < 4; degree++)
    {
        std::cout << Bon << "  Degree " << degree << Boff << std::endl;  
        for (auto& cl : msh)
        {
            std::cout << diameters(msh, cl).transpose() << std::endl;
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

            std::cout << stiff << std::endl;
            //std::cout << "  Diameter " << diameter(msh, cl) << std::endl;
        }
    }
}

#if 0
int main2(int argc, char **argv)
{
    using T = double;

    disk::simplicial_mesh<T,3> msh;
    auto mesher = make_simple_mesher(msh);
    for (size_t i = 0; i < 5; i++)
        mesher.refine();

    std::cout << msh.cells_size() << std::endl;
    
    
    using mesh_type = disk::simplicial_mesh<T,3>;
    using cell_type = typename mesh_type::cell_type;
    using face_type = typename mesh_type::face_type;
    using ScalT = double;
    using cell_monomial_basis = scalar_monomial<mesh_type, cell_type, ScalT>;
    
    timecounter tc;
    
    disk::dynamic_matrix<T> V = disk::dynamic_matrix<T>::Random(56,74);
    Eigen::DiagonalMatrix<T, Eigen::Dynamic> A(74);

    tc.tic();
    for (auto& cl : msh) {
        //cell_monomial_basis phi(msh, cl, 5);
        disk::dynamic_matrix<T> mass = V*A*V.transpose();
    }
    std::cout << tc.toc() << std::endl;

    tc.tic();
    for (auto& cl : msh) {
        cell_monomial_basis phi(msh, cl, 5);
        disk::dynamic_matrix<T> mass = integrate(msh, cl, phi, phi);
    }
    std::cout << tc.toc() << std::endl;
    
    tc.tic();
    for (auto& cl : msh) {
        cell_monomial_basis phi(msh, cl, 5);
        disk::dynamic_matrix<T> stiffness = integrate(msh, cl, grad(phi), grad(phi));
    }
    std::cout << tc.toc() << std::endl;

    return 0;
}
#endif

int main(int argc, char **argv)
{
    using T = double;

    disk::generic_mesh<T,1> msh_1D;
    make_single_element_mesh(msh_1D, 0.0, 1.0);

    //test_conditioning(msh_1D, 2.0, {1.0}, rescaling_strategy::none);
    //test_conditioning(msh_1D, 1.0, {1.0}, rescaling_strategy::none);

    disk::generic_mesh<T,2> msh_2D;
    make_single_element_mesh(msh_2D, 0.4, 5);
    double cf = 0.1;
    double theta = M_PI/4;
    double c = std::cos(theta);
    double s = std::sin(theta);
    auto tr2D = [&](const disk::point<T,2>& pt) -> disk::point<T,2> {
        double xn = pt.x()*c - pt.y()*s*cf;
        double yn = pt.x()*s + pt.y()*c*cf;
        return {xn, yn};
        //return {pt.x(), 0.1*pt.y()};
    };
    msh_2D.transform(tr2D);

    //disk::silo_database db;
    //db.create("elem.silo");
    //db.add_mesh(msh_2D, "mesh");

    //std::cout << " ----------- NONE " << std::endl;
    test_conditioning(msh_2D, 2.0, {1.0, 0.0}, rescaling_strategy::none);
    //std::cout << " ----------- INERTIAL " << std::endl;
    //test_conditioning(msh_2D, 2.0, {1.0, 0.0}, rescaling_strategy::inertial);
    //std::cout << " ----------- G-S " << std::endl;
    //test_conditioning(msh_2D, 2.0, {1.0, 0.0}, rescaling_strategy::gram_schmidt);

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
