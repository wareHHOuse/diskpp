#include <iostream>
#include <iomanip>
#include <regex>

#include <sstream>
#include <string>

#include <unistd.h>

#include "bases/bases.hpp"
#include "quadratures/quadratures.hpp"
#include "methods/hho"
#include "methods/implementation_hho/curl.hpp"
#include "core/loaders/loader.hpp"
#include "output/silo.hpp"
#include "solvers/solver.hpp"

template<typename Mesh>
void rectest(Mesh& msh, size_t order)
{
    using T = typename Mesh::coordinate_type;
    using point_type = typename Mesh::point_type;

    disk::hho_degree_info chdi( { .rd = (size_t) order+1,
                                  .cd = (size_t) order,
                                  .fd = (size_t) order } );

    auto fun = [&](const point_type& pt) -> Matrix<T, 3, 1> {
        Matrix<T, 3, 1> ret;
        ret(0) = 0.0;
        ret(1) = std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
        ret(2) = std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());

        //ret(0) = M_PI * M_PI * std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y())*std::sin(M_PI*pt.z());
        //ret(1) = M_PI * M_PI * std::cos(M_PI*pt.x())*std::cos(M_PI*pt.y())*std::sin(M_PI*pt.z());
        //ret(2) = M_PI * M_PI * std::cos(M_PI*pt.x())*std::sin(M_PI*pt.y())*std::cos(M_PI*pt.z());

        return ret;
    };

    std::vector<T> data_ex, data_ey, data_ez;

    for (auto& cl : msh)
    {
        auto cb = make_vector_monomial_basis(msh, cl, chdi.cell_degree());
        auto rb = make_vector_monomial_basis(msh, cl, chdi.reconstruction_degree());

        auto CR = disk::curl_reconstruction(msh, cl, chdi);

        Matrix<T, Dynamic, 1> pfun = disk::project_tangent(msh, cl, chdi, fun, 1);
        
        Matrix<T, Dynamic, 1> rfun = Matrix<T, Dynamic, 1>::Zero(rb.size());
        rfun.segment(3, rb.size()-3) = CR.first * pfun;

        auto bar = barycenter(msh, cl);
        auto h = diameter(msh, cl) / 10.;
        auto bar_h = point_type{bar.x() + h, bar.y() + h, bar.z() + h};
        auto Rphi = rb.eval_functions(bar_h);

        Matrix<T,3,1> Rval = disk::eval(rfun, Rphi);
        data_ex.push_back( Rval(0) );
        data_ey.push_back( Rval(1) );
        data_ez.push_back( Rval(2) );
    }

    disk::silo_database silo_db;
    silo_db.create("maxwell.silo");
    silo_db.add_mesh(msh, "mesh");

    disk::silo_zonal_variable<T> ex("Ru_x", data_ex);
    silo_db.add_variable("mesh", ex);

    disk::silo_zonal_variable<T> ey("Ru_y", data_ey);
    silo_db.add_variable("mesh", ey);

    disk::silo_zonal_variable<T> ez("Ru_z", data_ez);
    silo_db.add_variable("mesh", ez);

    silo_db.add_expression("Ru", "{Ru_x, Ru_y, Ru_z}", DB_VARTYPE_VECTOR);

    silo_db.close();
}

int main(int argc, char **argv)
{
    using T = double;

    char *      mesh_filename = nullptr;
    size_t      degree;
    int ch;
    while ( (ch = getopt(argc, argv, "m:k:")) != -1 )
    {
        switch(ch)
        {

            case 'm':
                mesh_filename = optarg;
                break;

            case 'k':
                degree = std::stoi(optarg);
                break;
            
            case '?':
            default:
                std::cout << "Invalid option" << std::endl;
                return 1;
        }
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<T>(mesh_filename);

        return 0;
    }

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        auto msh = disk::load_fvca5_2d_mesh<T>(mesh_filename);

        return 0;
    }

    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        auto msh = disk::load_netgen_3d_mesh<T>(mesh_filename);

        rectest(msh, degree);

        return 0;
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        auto msh = disk::load_cartesian_3d_mesh<T>(mesh_filename);
        
        rectest(msh, degree);
        
        return 0;
    }

    /* FVCA6 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
        disk::generic_mesh<T,3> msh;
        
        rectest(msh, degree);
        
        return 0;
    }
}