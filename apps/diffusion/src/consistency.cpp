/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#define PRINT_RANKS_AND_OTHER_STUFF
#include "diskpp/loaders/loader.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/methods/hho"
#include "diskpp/bases/bases_utils.hpp"

#include "sgr.hpp"

#include <unistd.h>

using namespace disk;
using namespace sgr;

enum class hho_variant {
    mixed_order_low,
    equal_order,
    mixed_order_high
};

template<typename Mesh>
auto
test_consistency(Mesh& msh, size_t degree, size_t increment, hho_variant hv)
{
    using T = typename Mesh::coordinate_type;

    hho_degree_info hdi;
    if (hv == hho_variant::mixed_order_low)
        hdi.cell_degree(degree-1);
    if (hv == hho_variant::equal_order)
        hdi.cell_degree(degree);
    if (hv == hho_variant::mixed_order_high)
        hdi.cell_degree(degree+1);
    hdi.face_degree(degree);
    hdi.reconstruction_degree(hdi.cell_degree()+2+increment);

    auto cd = hdi.cell_degree();
    auto fd = hdi.face_degree();
    auto rd = hdi.cell_degree()+2;
    auto hd = increment;

    auto dofsT = ((cd+2)*(cd+1))/2;
    auto dofsF = fd+1;
    auto dofsR = ((rd+2)*(rd+1))/2;
    auto dofsH = 2*hd;

    std::cout << Bwhitefg << "HHO(" << cd << "," << fd << "). Reconstruction: ";
    std::cout << rd+increment << " (poly: " << cd+2 << ", harmonic increment: ";
    std::cout << increment << ")" << reset << std::endl;

    for (auto& cl : msh)
    {
        auto num_faces = faces(msh,cl).size();

        /* Compute space dimensions and print them */
        auto totfrom = dofsT + num_faces*dofsF;
        auto totto = dofsR + dofsH;

        std::cout << "dofsT = " << dofsT << ", dofsF = " << dofsF;
        std::cout << ", total = " << totfrom << ", ";
        if (totto >= totfrom)
            std::cout << greenfg << "rec dofs = " << totto;
        else
            std::cout << redfg << "rec dofs = " << totto;

        std::cout << " (poly dofs: " << dofsR << ", harmonic dofs: " << dofsH <<")";
        std::cout << reset << std::endl;

        /* This is the polynomial with all the monomials up to some degree */
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
        
        /* Polynomial of degree k+1 to check consistency */
        Eigen::Matrix<T, Eigen::Dynamic, 1> poly_h = disk::project_function(msh, cl, hb, poly);

        /* These are the coefficients of the HHO function (vT, vF1, ..., vFn) */
        Eigen::Matrix<T, Eigen::Dynamic, 1> poly_hho = project_function(msh, cl, hdi, poly, 2);


        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> GR;

        Eigen::Matrix<T,2,2> Id = Eigen::Matrix<T,2,2>::Identity();
        if (hv == hho_variant::mixed_order_high) {
            /* Build the standard reconstruction + projection on the cells */
            auto oper = make_shl_face_proj_harmonic(msh, cl, hdi, Id);
            A = oper.second;
            GR = oper.first;
            std::cout << "+GR rows: " << GR.rows() << std::endl;
        } else {
            /* Build the nested reconstruction */
            auto oper = make_sfl(msh, cl, hdi, Id);
            A = oper.second;
            GR = oper.first;
            std::cout << "-GR rows: " << GR.rows() << std::endl;
        }

        Eigen::Matrix<T, Eigen::Dynamic, 1> poly_rec =
            GR * poly_hho;

        Eigen::Matrix<T, Eigen::Dynamic, 1> diff =
            poly_rec.head( cb.size()-1 ) - poly_k1.tail( cb.size()-1 );

        std::cout << Greenfg << "DIFFERENCE NORM: " << diff.norm() << nofg << std::endl;
        std::cout << Bwhitefg << "Reconstructed: " << reset << poly_rec.transpose() << std::endl;
        std::cout << Bwhitefg << "Expected:      " << reset << poly_h.tail(poly_h.rows()-1).transpose() << std::endl;

        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> eigs = A.eigenvalues();
        Eigen::Matrix<T, Eigen::Dynamic, 1> eigsa = eigs.cwiseAbs();
        
        std::cout << Cyanfg << "Eigenvalues:" << reset << std::endl;

        int zero_eigs = 0;
        T tol = 1e-10;
        T min_eig = eigsa[0];
        T max_eig = eigsa[0];
        for (size_t i = 1; i < eigsa.rows(); i++) {
            min_eig = std::min(min_eig, eigsa[i]);
            max_eig = std::max(max_eig, eigsa[i]);
        }

        T min2_eig = max_eig;
        for (size_t i = 0; i < eigsa.rows(); i++)
        {
            if ( eigsa[i] < tol ) {
                zero_eigs++;
                std::cout << Bcyanfg << eigsa[i] << " " << reset;
            }
            else {
                std::cout << eigsa[i] << " ";
            }

            if (eigsa[i] > min_eig)
                min2_eig = std::min(eigsa[i], min2_eig);
        }

        std::cout << Cyanfg << "[" << zero_eigs << " null eigs]" << reset << std::endl;
        std::cout << "Smallest nonzero eig: " << Magentafg << min2_eig << reset << std::endl;

    }
}


int main(int argc, char **argv)
{
    using T = double;

    int degree;
    int increment;
    hho_variant variant = hho_variant::equal_order;
    std::string mesh_filename;

    int ch;
    while ( (ch = getopt(argc, argv, "k:i:v:m:")) != -1 )
    {
        switch(ch)
        {
            case 'k':
                degree = std::stoi(optarg);
                break;

            case 'i':
                increment = std::stoi(optarg);
                break;

            case 'v':
                if (std::string(optarg) == "low")
                    variant = hho_variant::mixed_order_low;
                else if (std::string(optarg) == "equal")
                    variant = hho_variant::equal_order;
                else if (std::string(optarg) == "high")
                    variant = hho_variant::mixed_order_high;
                else {
                    std::cout << "Unknown variant '" << optarg << "'" << std::endl;
                    return 1;
                } 
                break;

            case 'm':
                mesh_filename = optarg;
                break;

            case '?':
            default:
                std::cout << "Invalid option" << std::endl;
                return 1;
        }
    }

    if ( mesh_filename.length() > 0 )
    {
        disk::generic_mesh<T,2> msh_gen;
        load_single_element_csv(msh_gen, mesh_filename);
        test_consistency(msh_gen, degree, increment, variant);
        return 0;
    }


    for (size_t i = 3; i < 11; i++) {
        std::cout << BYellowfg << " **** Faces = " << i << " ****" << reset << std::endl;
        disk::generic_mesh<T,2> msh_gen;
        disk::make_single_element_mesh(msh_gen, 1.0, i);
        test_consistency(msh_gen, degree, increment, variant);
    }    

    return 0;
}
