#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/hho"
#include "diskpp/bases/bases_utils.hpp"
#include "diskpp/output/silo.hpp"

#include "sgr.hpp"

using namespace disk;
using namespace sgr;

template<typename Mesh>
auto
test_consistency(Mesh& msh, size_t degree, size_t increment)
{
    using T = typename Mesh::coordinate_type;

    hho_degree_info hdi;
    hdi.cell_degree(degree+1);
    hdi.face_degree(degree);

    auto delta = (hdi.cell_degree() > hdi.face_degree()) ? 1 : 2;

    hdi.reconstruction_degree(hdi.cell_degree()+delta+increment);

    for (auto& cl : msh)
    {
        auto cd = hdi.cell_degree();
        auto fd = hdi.face_degree();
        auto rd = (delta == 1) ? hdi.reconstruction_degree() : cd+1;
        auto hd = (delta == 1) ? 0 : increment;

        auto dofsT = ((cd+2)*(cd+1))/2;
        auto dofsF = fd+1;
        auto dofsR = ((rd+2)*(rd+1))/2;
        auto dofsH = 2*hd;
        auto n = faces(msh,cl).size();

        auto totfrom = dofsT + n*dofsF;
        auto totto = dofsR + dofsH;

        std::cout << Bon << "HHO(" << cd << "," << fd << ") -> (" << rd << "," << hd;
        std::cout << ")" << reset << std::endl; 

        std::cout << "dofsT = " << dofsT << ", dofsF = " << dofsF;
        std::cout << ", total = " << totfrom << ", ";
        if (totto >= totfrom)
            std::cout << greenfg << "rec dofs = " << totto;
        else
            std::cout << redfg << "rec dofs = " << totto;
        
        std::cout << reset << std::endl;

        auto poly = [&](const point<T,2>& pt) {
            auto x = pt.x();
            auto y = pt.y();

            auto ret = 0.0;

            auto deg = degree+1;

            for (size_t k = 0; k <= deg; k++)
                for (size_t i = 0; i <= k; i++)
                    ret += std::pow(x,k-i)*std::pow(y,i);

            return ret;
        };

        auto cb = make_scalar_monomial_basis(msh, cl, degree+1);
        Eigen::Matrix<T, Eigen::Dynamic, 1> poly_k1 =
            disk::project_function(msh, cl, cb, poly);

        auto hb = make_scalar_harmonic_top_basis(msh, cl, hdi.reconstruction_degree());
        hb.maximum_polynomial_degree(hdi.cell_degree()+2);

        //auto hb = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
        
        Eigen::Matrix<T, Eigen::Dynamic, 1> poly_h =
            disk::project_function(msh, cl, hb, poly);

        Eigen::Matrix<T, Eigen::Dynamic, 1> poly_hho =
            project_function(msh, cl, hdi, poly, 2);


        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> GR;

        Eigen::Matrix<T,2,2> Id = Eigen::Matrix<T,2,2>::Identity();
        if (delta == 1) {
            auto oper = make_shl_face_proj(msh, cl, hdi, Id);
            A = oper.second;
            GR = oper.first;
        } else {
            auto oper = make_sfl(msh, cl, hdi, Id);
            A = oper.second;
            GR = oper.first;
        }


        Eigen::Matrix<T, Eigen::Dynamic, 1> poly_rec =
            GR * poly_hho;

        Eigen::Matrix<T, Eigen::Dynamic, 1> diff =
            poly_rec.head( cb.size()-1 ) - poly_k1.tail( cb.size()-1 );

        std::cout << greenfg << "DIFFERENCE NORM: " << diff.norm() << nofg << std::endl;
        std::cout << "Reconstructed:   " << poly_rec.transpose() << std::endl;
        std::cout << "Expected:        " << poly_h.tail(poly_h.rows()-1).transpose() << std::endl;
        std::cout << "k+1 degree proj: " << poly_k1.tail(poly_k1.rows()-1).transpose() << std::endl;

        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> eigs = A.eigenvalues();
        Eigen::Matrix<T, Eigen::Dynamic, 1> eigsa = eigs.cwiseAbs();
        
        std::cout << Cyanfg << "Eigenvalues:" << reset << std::endl;

        int zero_eigs = 0;
        T tol = 1e-10;
        T min_eig = std::max(eigsa[0], tol);
        for (size_t i = 0; i < eigsa.rows(); i++)
        {
            if ( eigsa[i] < tol ) {
                zero_eigs++;
                std::cout << Bcyanfg << eigsa[i] << " " << reset;
            }
            else {
                min_eig = std::min(min_eig, eigsa[i]);
                std::cout << eigsa[i] << " ";
            }
        }

        std::cout << Bcyanfg << "[" << zero_eigs << " null eigs]" << reset << std::endl;
        std::cout << "Smallest nonzero eig: " << min_eig << std::endl;

    }
}


int main(int argc, char **argv)
{
    using T = double;

    size_t degree = 1;
    size_t increment = 0;

    if (argc > 1)
        degree = std::stoi(argv[1]);

    if (argc > 2)
        increment = std::stoi(argv[2]);

    disk::generic_mesh<T,2> msh_gen;
    load_single_element_csv(msh_gen, argv[3]);
    test_consistency(msh_gen, degree, increment);

    /*
    triangular_mesh<T> msh;
    auto mesher = make_simple_mesher(msh);
    //mesher.refine();
    test_consistency(msh);

    for (auto& cl : msh) {
        std::cout << cl << std::endl;
        auto pts = points(msh, cl);
        for (auto& pt : pts) {
            std::cout << pt << " ";
        }
        std::cout << std::endl;
    }
    */

    /*
    for (size_t i = 0; i < 3; i++)
    {
        auto div = std::pow(2, i);
        std::cout << Yellowfg << "Triangle, div = " << div << reset << std::endl;
        disk:triangular_mesh<T> msh_tri;
        make_single_element_mesh(msh_tri, {0.0,0.0}, {1.0/div,0.0}, {0.7/div,0.5/div});
        //make_single_element_mesh(msh_tri, {1.0, 0.0}, {-0.5, 0.866025}, {-0.5, -0.866025});
        test_consistency(msh_tri, degree, increment);
    }
    */

    /*    
    for (size_t i = 0; i < 3; i++)
    {
        auto div = std::pow(2, i);
        std::cout << yellowfg << "Quad, div = " << div << nofg << std::endl;
        disk::cartesian_mesh<T,2> msh_quad;
        disk::make_single_element_mesh(msh_quad, {0,0}, 1./div, 1./div);
        test_consistency(msh_quad, degdeg);
    }
    */

    /*
    for (size_t i = 3; i < 10; i++) {
        std::cout << Yellowfg << " **** Faces = " << i << " ****" << nofg << std::endl;
        disk::generic_mesh<T,2> msh_gen;
        disk::make_single_element_mesh(msh_gen, 1.0, i);
        test_consistency(msh_gen, degree, increment);
    }
    */


    return 0;
}
