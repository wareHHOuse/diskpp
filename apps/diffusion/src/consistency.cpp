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
#include "diskpp/loaders/loader_utils.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/methods/hho"
#include "diskpp/bases/bases_utils.hpp"
#include "diskpp/output/silo.hpp"

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
    std::ofstream ofs("eigdata.dat");

    using T = typename Mesh::coordinate_type;
    using point_type = typename Mesh::point_type;

    const size_t cr_delta = 2;

    hho_degree_info hdi;
    if (hv == hho_variant::mixed_order_low)
        hdi.cell_degree(degree-1);
    if (hv == hho_variant::equal_order)
        hdi.cell_degree(degree);
    if (hv == hho_variant::mixed_order_high)
        hdi.cell_degree(degree+1);
    hdi.face_degree(degree);
    hdi.reconstruction_degree(hdi.cell_degree()+cr_delta+increment);

    auto cd = hdi.cell_degree();
    auto fd = hdi.face_degree();
    auto rd = hdi.cell_degree()+cr_delta;
    auto hd = increment;

    auto dofsT = scalar_basis_size(cd, Mesh::dimension);
    auto dofsF = scalar_basis_size(fd, Mesh::dimension-1);
    auto dofsR = scalar_basis_size(rd, Mesh::dimension);
    auto dofsH = harmonic_basis_size(hd+rd, Mesh::dimension) - harmonic_basis_size(rd, Mesh::dimension);

    std::cout << Bwhitefg << "HHO(" << cd << "," << fd << "). Reconstruction: ";
    std::cout << rd+increment << " (poly: " << rd << ", harmonic increment: ";
    std::cout << increment << ")" << reset << std::endl;

    double Tmin_eig = 99999999.0;
    double Tmax_eig = 0.0;

    std::vector<double> var_hv;
    std::vector<double> var_eigv;

    for (auto& cl : msh)
    {
        auto num_faces = faces(msh,cl).size();

        /* Compute space dimensions and print them */
        auto totfrom = dofsT + num_faces*dofsF;
        auto totto = dofsR + dofsH;

        std::cout << "dofsT = " << dofsT << ", dofsF = " << dofsF;
        std::cout << ", total = " << totfrom << ", ";
        if (totto >= totfrom)
            std::cout << Greenfg << "rec dofs = " << totto;
        else
            std::cout << redfg << "rec dofs = " << totto;

        std::cout << " (poly dofs: " << dofsR << ", harmonic dofs: " << dofsH <<")";
        std::cout << reset << std::endl;

        /* This is the polynomial with all the monomials up to some degree */
        auto poly = [&](const point_type& pt) {
            auto ret = 0.0;
            auto deg = degree+1;
            for (size_t k = 0; k <= deg; k++)
                for (size_t i = 0; i <= k; i++)
                    ret += std::pow(pt.x(), k-i)*std::pow(pt.y(), i);

            return ret;
        };

        auto cb = make_scalar_monomial_basis(msh, cl, degree+1);
        Eigen::Matrix<T, Eigen::Dynamic, 1> poly_k1 =
            disk::project_function(msh, cl, cb, poly);

        auto hb = make_scalar_harmonic_top_basis(msh, cl, hdi.reconstruction_degree());
        hb.maximum_polynomial_degree(rd);

        disk::dynamic_vector<T> p = cb.eval_functions(barycenter(msh, cl));
        disk::dynamic_vector<T> q = hb.eval_functions(barycenter(msh, cl)).head(cb.size());
        std::cout << "DIFF pq: " << (p-q).norm() << std::endl;

        disk::dynamic_matrix<T> Gp = cb.eval_gradients(barycenter(msh, cl));
        disk::dynamic_matrix<T> Gq = hb.eval_gradients(barycenter(msh, cl)).block(0,0,cb.size(),Mesh::dimension);
        std::cout << "DIFF GpGq: " << (Gp-Gq).norm() << std::endl;

        /* Polynomial of degree k+1 to check consistency */
        Eigen::Matrix<T, Eigen::Dynamic, 1> poly_h = disk::project_function(msh, cl, hb, poly);

        /* These are the coefficients of the HHO function (vT, vF1, ..., vFn) */
        Eigen::Matrix<T, Eigen::Dynamic, 1> poly_hho = project_function(msh, cl, hdi, poly, 2);


        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> GR;

        Eigen::Matrix<T,Mesh::dimension,Mesh::dimension> Id =
            Eigen::Matrix<T,Mesh::dimension,Mesh::dimension>::Identity();
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

        assert(GR.rows() == totto-1);

        Eigen::Matrix<T, Eigen::Dynamic, 1> poly_rec = GR * poly_hho;

        Eigen::Matrix<T, Eigen::Dynamic, 1> diff =
            poly_rec.head( cb.size()-1 ) - poly_k1.tail( cb.size()-1 );

        if (diff.norm() > 1e-12) std::cout << BRedfg; else std::cout << Greenfg;
        std::cout << "DIFFERENCE NORM: " << diff.norm() << nofg << std::endl;
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
                std::cout << Cyanfg << eigsa[i] << " " << reset;
            }
            else {
                std::cout << eigsa[i] << " ";
            }

            if (eigsa[i] > min_eig)
                min2_eig = std::min(eigsa[i], min2_eig);
        }

        std::cout << Cyanfg << "[" << zero_eigs << " null eigs]" << reset << std::endl;
        if (zero_eigs == 1) {
            auto h = diameter(msh, cl);
            if (Mesh::dimension == 2)
                std::cout << "Smallest nonzero eig: " << BMagentafg << min2_eig << reset << std::endl;
            if (Mesh::dimension == 3) {
                
                std::cout << "Smallest nonzero eig: " << BMagentafg << min2_eig;
                std::cout << " (h = " << h << ", corrected by h: " << min2_eig/h << ")" << reset << std::endl;
            }

            Tmin_eig = std::min(Tmin_eig, min2_eig/h);
            Tmax_eig = std::max(Tmax_eig, min2_eig/h);

            var_hv.push_back(h);
            var_eigv.push_back(min2_eig);

            ofs << h << " " << min2_eig << " " << min2_eig/h << " " << measure(msh, cl) << " ";
            
            
            auto pts = points(msh, cl);
            for (auto& pt : pts)
                ofs << pt << " ";
        
            ofs << "bar: " << barycenter(msh, cl) << " ";

            ofs << std::endl;

            continue;
        }

        std::cout << "NOT STABLE" << std::endl;
        break;
    }

    std::cout << Tmin_eig << " " << Tmax_eig << std::endl;

    disk::silo_database db;
    db.create("consist.silo");
    db.add_mesh(msh, "mesh");
    db.add_variable("mesh", "h", var_hv, disk::zonal_variable_t);
    db.add_variable("mesh", "eig", var_eigv, disk::zonal_variable_t);
}


int main(int argc, char **argv)
{
    using T = double;

    int degree = 0;
    int increment = 0;
    int refs = 0;
    double radius = 1.0;
    hho_variant variant = hho_variant::equal_order;
    std::string mesh_filename;

    int ch;
    while ( (ch = getopt(argc, argv, "k:i:v:m:r:z:")) != -1 )
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

            case 'r':
                radius = std::stod(optarg);
                if (radius <= 0)
                    radius = 1.0;
                break;

            case 'z':
                refs = std::stod(optarg);
                if (refs < 0)
                    refs = 0;
                break;

            case '?':
            default:
                std::cout << "Invalid option" << std::endl;
                return 1;
        }
    }

    if ( mesh_filename.length() > 0 )
    {
        /* Single element CSV 2D */
        if (std::regex_match(mesh_filename, std::regex(".*\\.csv$") ))
        {
            std::cout << "Guessed mesh format: CSV 2D" << std::endl;
            disk::generic_mesh<T,2> msh_gen;
            load_single_element_csv(msh_gen, mesh_filename);
            test_consistency(msh_gen, degree, increment, variant);
            std::cout << std::endl;
            return 0;
        }

        /* FVCA5 2D */
        if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
        {
            std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
            disk::generic_mesh<T,2> msh;
            disk::load_mesh_fvca5_2d<T>(mesh_filename.c_str(), msh);
            test_consistency(msh, degree, increment, variant);
            return 0;
        }

        /* FVCA6 3D */
        if (std::regex_match(mesh_filename, std::regex(".*\\.msh$") ))
        {
            std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
            disk::generic_mesh<T,3> msh;
            disk::load_mesh_fvca6_3d<T>(mesh_filename.c_str(), msh);
            test_consistency(msh, degree, increment, variant);
            return 0;
        }

        if (std::regex_match(mesh_filename, std::regex(".*\\.geo3s$") ))
        {
            std::cout << "Guessed mesh format: GMSH simplicial 3D" << std::endl;
            disk::simplicial_mesh<T,3> msh;
            disk::gmsh_geometry_loader< disk::simplicial_mesh<T,3> > loader;
            loader.read_mesh(mesh_filename);
            loader.populate_mesh(msh);
            test_consistency(msh, degree, increment, variant);
            return 0;
        }
    }

    /*
    for (size_t i = 3; i < 11; i++) {
        std::cout << BYellowfg << " **** Faces = " << i << " ****" << reset << std::endl;
        disk::generic_mesh<T,2> msh_gen;
        disk::make_single_element_mesh(msh_gen, radius, i);
        test_consistency(msh_gen, degree, increment, variant);
        std::cout << std::endl;
    }
    */

    {
        using mesh_type = disk::simplicial_mesh<T,3>;
        mesh_type msh;
        auto mesher = disk::make_simple_mesher(msh);
        for (size_t i = 0; i < refs; i++)
            mesher.refine();
        test_consistency(msh, degree, increment, variant);
    }

    return 0;
}
