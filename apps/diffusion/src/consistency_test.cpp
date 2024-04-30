/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2024
 * matteo.cicuttin@polito.it
 *
 * Karol Cascavita (C) 2024
 * karol.cascavita@polito.it

 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#define PRINT_RANKS_AND_OTHER_STUFF
#include "diskpp/loaders/loader.hpp"
#include "diskpp/loaders/loader_utils.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/methods/hho"
#include "diskpp/bases/bases_utils.hpp"

#include "sgr.hpp"

#include <unistd.h>
#include <numeric>

using namespace disk;
using namespace sgr;

enum class hho_variant {
    mixed_order_low,
    equal_order,
    mixed_order_high
};

std::string dir = "/home/karol/Documents/UNIVERSITA/POLITO/STABFREE/stabfree_hho/test_coercivita/polygons/";

std::ostream& operator <<(std::ostream& os, const hho_variant& o)
{
    switch(o)
    {
        case (hho_variant::mixed_order_low):
            os << "low-mixed";
            return os;
        case (hho_variant::equal_order):
            os << "equal";
            return os;
        case (hho_variant::mixed_order_high):
            os << "high-mixed";
            return os;
        default:
            throw std::invalid_argument("ERROR: hho_variant not found.");
    }
}


using sparse_matrix_t = Eigen::SparseMatrix<double>; 
using triplet_t = Eigen::Triplet<double>;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

template<typename Mesh>
auto
test_consistency(Mesh& msh, size_t degree, size_t increment, hho_variant hv)
{
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
            //std::cout << "+GR rows: " << GR.rows() << std::endl;
        } else {
            /* Build the nested reconstruction */
            auto oper = make_sfl(msh, cl, hdi, Id);
            A = oper.second;
            GR = oper.first;
            //std::cout << "-GR rows: " << GR.rows() << std::endl;
        }

        Eigen::Matrix<T, Eigen::Dynamic, 1> poly_rec =
            GR * poly_hho;

        Eigen::Matrix<T, Eigen::Dynamic, 1> diff =
            poly_rec.head( cb.size()-1 ) - poly_k1.tail( cb.size()-1 );

        if (diff.norm() > 1e-12) std::cout << BRedfg; else std::cout << Greenfg;
        //std::cout << "DIFFERENCE NORM: " << diff.norm() << nofg << std::endl;
        //std::cout << Bwhitefg << "Reconstructed: " << reset << poly_rec.transpose() << std::endl;
        //std::cout << Bwhitefg << "Expected:      " << reset << poly_h.tail(poly_h.rows()-1).transpose() << std::endl;

        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> eigs = A.eigenvalues();
        Eigen::Matrix<T, Eigen::Dynamic, 1> eigsa = eigs.cwiseAbs();
        
        //std::cout << Cyanfg << "Eigenvalues:" << reset << std::endl;

        int zero_eigs = 0;
        T tol = 1e-10;
        T min_eig = eigsa[0];
        T max_eig = eigsa[0];
        for (size_t i = 1; i < eigsa.rows(); i++) {
            min_eig = std::min(min_eig, eigsa[i]);
            max_eig = std::max(max_eig, eigsa[i]);
        }

        T min2_eig = max_eig;Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
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

        //std::cout << Cyanfg << "[" << zero_eigs << " null eigs]" << reset << std::endl;
        if (zero_eigs == 1)
        {
            std::cout << std::endl;
            std::cout << "Smallest nonzero eig: " << BMagentafg << min2_eig << reset << std::endl;
            return std::make_pair(true, min2_eig);
        }
 
        break;
    }
    return std::make_pair(false, 0.0);
}


enum class polygon_type
{
    REGULAR,
    IRREGULAR,
    IRREGULAR_ALIGNED,
    ALIGNED,
    CONCAVE,
    SINGLE,
};

std::ostream& operator <<(std::ostream& os, const polygon_type& o)
{
    switch(o)
    {
        case (polygon_type::REGULAR):
            os << "Regular";
            return os;
        case (polygon_type::IRREGULAR):
            os << "Irregular";
            return os;
        case (polygon_type::IRREGULAR_ALIGNED):
            os << "IrregularAligned";
            return os;
        case (polygon_type::ALIGNED):
            os << "Aligned";
            return os;
        case (polygon_type::SINGLE):
            os << "Single";
            return os;
        case (polygon_type::CONCAVE):
            os << "Concave";
            return os;
        default:
            throw std::invalid_argument("ERROR: polygon type not found.");
    }
}


enum class generated_by
{
    CSV,
    DISKPP,
};


bool test_polygons_first_try(size_t degree, double radius,
                    const hho_variant& variant, 
                    polygon_type polygon)
{
    size_t max_increment = 40; //This can be changed.
    generated_by g;
    std::vector<int> poly_faces(10);

    switch(polygon)
    {
        case(polygon_type::REGULAR):
        case(polygon_type::IRREGULAR):
        case(polygon_type::IRREGULAR_ALIGNED):
            g = generated_by::CSV;
            std::iota(poly_faces.begin(), poly_faces.end(), 3); 
            break;
        case(polygon_type::ALIGNED):    
            g = generated_by::CSV;
            std::iota(poly_faces.begin(), poly_faces.end(), 5); 
            break;
        case(polygon_type::CONCAVE):
            g = generated_by::CSV;
            poly_faces = std::vector<int>({4,5,6,7,9,11,13,15,17,19});
            break;
        case(polygon_type::SINGLE):
            g = generated_by::DISKPP;
            std::iota(poly_faces.begin(), poly_faces.end(), 3); 
            break;
        default:
            throw std::invalid_argument("POLYGON TYPE not found");        
    }

    matrix_t mat = matrix_t::Zero(10,3);

    for (size_t f = 0; f < 10; f++) 
    {
        size_t num_faces = poly_faces[f];     
        std::cout << BYellowfg << " **** Faces = " << num_faces << " ****" << reset << std::endl;
        disk::generic_mesh<double,2> msh_gen;
        switch(g)
        {
            case(generated_by::DISKPP):
                disk::make_single_element_mesh(msh_gen, radius, num_faces);
                break;
            case(generated_by::CSV): 
            {
                std::stringstream mesh_filename;
                mesh_filename << dir << "Test"<< polygon << "/Csv/Domain_" << num_faces << ".csv"; 
                std::cout << mesh_filename.str() << std::endl;
                load_single_element_csv(msh_gen, mesh_filename.str());
                break;
            }   
            default:
                throw std::invalid_argument("");
        }

        size_t increment;
        for( increment = 0; increment <= max_increment; increment++)   
        {
            std::cout<<"------------------------------------------------------------"<<std::endl;
            std::cout<<"* Increment..."<< increment <<std::endl;
            std::cout<<"------------------------------------------------------------"<<std::endl;

            auto [pass, min2_eig] = test_consistency(msh_gen, degree, increment, variant);
            double val = (pass)? min2_eig : -1;

            if(pass || increment == max_increment)
            {
                std::cout << std::endl;
                mat(f, 0) = num_faces;
                mat(f, 1) = val;
                mat(f, 2) = increment;
                break;
            }
        }
    }

    std::cout<<"======================================================="<<std::endl;
    std::cout<<"======================================================="<<std::endl;

    std::stringstream name;
    name << "matrix_resume_"<< polygon << "_" << variant << "_k_" << degree << ".txt";
    std::ofstream ofs(name.str());               
            
    for (int i = 0; i < mat.rows(); ++i)
    {
        ofs << std::setprecision(2) << mat(i,0) << " ,"
            << std::setprecision(4) << std::scientific << mat(i,1) << " ,"
            << std::setprecision(3) << mat(i,2) << std::endl ;
    }
    ofs.close();            
    std::cout<<"======================================================="<<std::endl;
    std::cout<<"======================================================="<<std::endl; 

    return true;    
}


bool test_polygons( size_t degree, double radius,
                    const hho_variant& variant, 
                    polygon_type polygon)
{
    generated_by g;
    std::vector<int> poly_faces(10);    
    
    switch(polygon)
    {
        case(polygon_type::REGULAR):
        case(polygon_type::IRREGULAR):
        case(polygon_type::IRREGULAR_ALIGNED):
            g = generated_by::CSV;
            std::iota(poly_faces.begin(), poly_faces.end(), 3); 
            break;
        case(polygon_type::ALIGNED):    
            g = generated_by::CSV;
            std::iota(poly_faces.begin(), poly_faces.end(), 5); 
            break;
        case(polygon_type::CONCAVE):
            g = generated_by::CSV;
            poly_faces = std::vector<int>({4,5,6,7,9,11,13,15,17,19});
            break;
        case(polygon_type::SINGLE):
            g = generated_by::DISKPP;
            std::iota(poly_faces.begin(), poly_faces.end(), 3); 
            break;
        default:
            throw std::invalid_argument("POLYGON TYPE not found");        
    }
    size_t max_faces = 10;

    std::vector<triplet_t> triplets;

    size_t global_increment = 40;
    for(size_t increment = 0; increment <=40 ; increment++)
    {
        bool pass_all = true;
        std::cout<<"------------------------------------------------------------"<<std::endl;
        std::cout<<"* Increment..."<< increment <<std::endl;
        std::cout<<"------------------------------------------------------------"<<std::endl;

        for (size_t f = 0; f < 10; f++) 
        {
            std::cout << BYellowfg << " **** Faces = " << poly_faces[f] << " ****" << reset << std::endl;
            disk::generic_mesh<double,2> msh_gen;
            switch(g)
            {
                case(generated_by::DISKPP):
                    disk::make_single_element_mesh(msh_gen, radius, poly_faces[f]);
                    break;
                case(generated_by::CSV): 
                {
                    std::stringstream mesh_filename;
                    mesh_filename << dir << "Test"<< polygon << "/Csv/Domain_" << poly_faces[f] << ".csv"; 
                    std::cout << mesh_filename.str() << std::endl;
                    load_single_element_csv(msh_gen, mesh_filename.str());
                    break;
                }   
                default:
                    throw std::invalid_argument("");
            }

            auto [pass, min2_eig] = test_consistency(msh_gen, degree, increment, variant);
            double val = (pass)? min2_eig : -1;
            pass_all = pass_all && pass; 

            triplets.push_back(triplet_t(poly_faces[f], increment, val));

        }

        if(pass_all)
        {
            global_increment = increment;
            break;
        }


    }

    sparse_matrix_t smat(11, global_increment+1);    
    smat.setFromTriplets(triplets.begin(), triplets.end()); 

    matrix_t mat = matrix_t(smat);

    std::stringstream name;
    name << "matrix_"<< polygon << "_" << variant << "_k_" << degree << ".txt";
    std::ofstream ofs(name.str());               
           
    for (int i = 0; i < mat.rows(); ++i)
    {
        ofs << poly_faces[i]<<",";
        for (int j = 0; j < mat.cols(); ++j)
            ofs << std::setprecision(4) <<  mat(i,j) << ((j==(mat.cols()-1))? "": ",");
        ofs<< "" << std::endl ;
    }
    ofs.close();            

    return true;
}

# if 0
void test_dispp_varying_degree(double radius , const hho_variant& variant)
{
    for(size_t degree = 0; degree <= 2; degree++)
        test_dispp_polygons(degree,radius, variant);
}

void test_dispp_high_order(double radius)
{
    hho_variant variant = hho_variant::mixed_order_high;
    test_dispp_varying_degree(radius, variant);
    return;
}

void test_dispp_low_order(double radius)
{
    hho_variant variant = hho_variant::mixed_order_low;
    test_dispp_varying_degree(radius, variant);
    return;
}

void test_dispp_equal_order(double radius)
{
    hho_variant variant = hho_variant::equal_order;
    test_dispp_varying_degree(radius, variant);
    return;
}
#endif

int main(int argc, char **argv)
{
    using T = double;

    int degree = 0;
    int increment = 0;
    double radius = 1.0;
    hho_variant variant = hho_variant::equal_order;
    std::string mesh_filename;
    polygon_type polygon = polygon_type::SINGLE;
    bool run_all = false;
    bool resume = false;

    int ch;
    while ( (ch = getopt(argc, argv, "ark:i:v:p:")) != -1 )
    {
        switch(ch)
        {
            case 'a':
                run_all =true;
                break;
            case 'r':
                resume = true;
                break;
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
            case 'p':
                if (std::string(optarg) == "single")
                    polygon = polygon_type::SINGLE;
                else if (std::string(optarg) == "regular")
                    polygon = polygon_type::REGULAR;
                else if (std::string(optarg) == "irregular")
                    polygon = polygon_type::IRREGULAR;
                else if (std::string(optarg) == "irregular_aligned")
                    polygon = polygon_type::IRREGULAR_ALIGNED;
                else if (std::string(optarg) == "aligned")
                    polygon = polygon_type::ALIGNED;
                else if (std::string(optarg) == "concave")
                    polygon = polygon_type::CONCAVE;
                else {
                    std::cout << "Unknown polygon type '" << optarg << "'" << std::endl;
                    return 1;
                } 
                break;
            case '?':
            default:
                std::cout << "Invalid option" << std::endl;
                return 1;
        }
    }

    if(resume)
    {
        if(run_all)
        {
            size_t start = (variant == hho_variant::mixed_order_low)? 1 : 0;  
            for(size_t degree = start; degree <= 2; degree++)
                test_polygons_first_try(degree, radius, variant, polygon);
            std::cout << "Running resume test with " << variant << " order and polygon "
                << polygon << std::endl;  
            return 0;
        }

        test_polygons_first_try(degree, radius, variant, polygon);
        std::cout << "Running test with K = "<<  degree <<" and with " << variant << " order and polygon "
                << polygon << std::endl;  

        return 0;
    }

    if(run_all)
    {
        size_t start = (variant == hho_variant::mixed_order_low)? 1 : 0;  
        for(size_t degree = start; degree <= 2; degree++)
            test_polygons(degree, radius, variant, polygon);
        std::cout << "Running test with " << variant << " order and polygon "
            << polygon << std::endl;  
        return 0;
    }

    test_polygons(degree, radius, variant, polygon);
    std::cout << "Running test with K = "<<  degree <<" and with " << variant << " order and polygon "
                << polygon << std::endl;  


    return 0;
}
