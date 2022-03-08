#include <iostream>
#include <iomanip>
#include <regex>
#include <fstream>
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

    auto fs = [&](const point_type& pt) -> T {
        return std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
    };

    std::vector<T> data_ex, data_ey, data_ez;

    T min = 9999999;
    T max = 0;
    T err = 0;

    std::ofstream errfile("error.dat");

    for (auto& cl : msh)
    {
        //std::cout << "**********" << std::endl;
        Matrix<T, Dynamic, 1> pfun = disk::project_function(msh, cl, chdi, fs, 1);
        auto GR = make_scalar_hho_laplacian(msh, cl, chdi);
        auto ST = make_scalar_hdg_stabilization(msh, cl, chdi);

        auto cb = make_scalar_monomial_basis(msh, cl, chdi.cell_degree());
        auto cbs = cb.size();

        auto ofs = cbs;
        auto fcs = faces(msh, cl);
        for (auto& fc : fcs)
        {
            //std::cout << "   +++++" << std::endl;
            auto fb = make_scalar_monomial_basis(msh, fc, chdi.face_degree());
            auto fbs = fb.size();

            Matrix<T, Dynamic, Dynamic> mass = Matrix<T, Dynamic, Dynamic>::Zero(fbs, fbs);
            Matrix<T, Dynamic, Dynamic> trace  = Matrix<T, Dynamic, Dynamic>::Zero(fbs, cbs);

            const auto qps = integrate(msh, fc, 2 * chdi.face_degree());
            for (auto& qp : qps)
            {
                const auto c_phi = cb.eval_functions(qp.point());
                const auto f_phi = fb.eval_functions(qp.point());

                mass += qp.weight() * (f_phi * f_phi.transpose());
                trace += qp.weight() * (f_phi * c_phi.transpose());
            }

            Matrix<T, Dynamic, 1> vT_F = mass.ldlt().solve(trace*pfun.segment(0,cbs));

            Matrix<T, Dynamic, 1> vF = pfun.segment(ofs, fb.size());
            Matrix<T, Dynamic, 1> diff = (vT_F - vF);
            /*
            std::cout << "   vT|F   = " << vT_F.transpose() << std::endl;
            std::cout << "   vF     = " << vF.transpose() << std::endl;
            std::cout << "   diff   = " << diff.transpose() << std::endl;
            std::cout << "   1/hF   = " << 1./diameter(msh, fc) << std::endl;
            std::cout << "   dot    = " << diff.dot(diff) << std::endl;
            std::cout << "   dotM   = " << diff.dot(mass*diff) << std::endl;
            */
            ofs += fb.size();
        }

        auto vSv = pfun.transpose().dot(ST*pfun);
        err += vSv;
        min = std::min(vSv, min);
        max = std::max(vSv, max);

        //std::cout << " v'Sv = " << vSv << " " << std::sqrt(vSv) << std::endl;

        JacobiSVD<MatrixXd> svd(GR.second + ST);
        double cond = svd.singularValues()(0) 
            / svd.singularValues()(svd.singularValues().size()-1);

        //std::cout << cond << std::endl;

        /*
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
        */

        auto bar = barycenter(msh, cl);
        errfile << bar.x() << " " << bar.y() << " " << vSv << std::endl;
    }

    std::cout << min << " " << max << " " << std::sqrt(err) << std::endl;

    /*
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
    */
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

        rectest(msh, degree);

        return 0;
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        auto msh = disk::load_cartesian_2d_mesh<T>(mesh_filename);
        
        rectest(msh, degree);
        
        return 0;
    }

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        auto msh = disk::load_fvca5_2d_mesh<T>(mesh_filename);

        rectest(msh, degree);

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
        auto msh = disk::load_fvca6_3d_mesh<T>(mesh_filename);
        rectest(msh, degree);
        
        return 0;
    }
}