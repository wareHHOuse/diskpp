/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *  
 * Matteo Cicuttin (C) 2020
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group  
 */

#include <iostream>
#include <iomanip>
#include <regex>

#include <sstream>
#include <string>

#include <unistd.h>

#include "bases/bases.hpp"
#include "quadratures/quadratures.hpp"
#include "loaders/loader.hpp"
#include "output/silo.hpp"
#include "solvers/solver.hpp"
#include "methods/dg"

#include "compinfo.h"

#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include "paramloader.hpp"

namespace disk {

/* Temporarly copied from curl.hpp */
template<typename T>
Matrix<T, Dynamic, 3>
vcross(const static_vector<T, 3>& normal, const Matrix<T, Dynamic, 3>& field)
{
    Matrix<T, Dynamic, 3>   ret;
    ret = Matrix<T, Dynamic, 3>::Zero(field.rows(), field.cols());

    for (size_t i = 0; i < field.rows(); i++)
    {
        static_vector<T, 3> f;
        f(0) = field(i,0);
        f(1) = field(i,1);
        f(2) = field(i,2);

        auto r = normal.cross( f );
        ret(i,0) = r(0);
        ret(i,1) = r(1);
        ret(i,2) = r(2);
    }

    return ret;
}

/* Temporarly copied from curl.hpp */
template<typename T>
Matrix<T, Dynamic, 3>
vcross(const Matrix<T, Dynamic, 3>& field, const static_vector<T, 3>& normal)
{
    Matrix<T, Dynamic, 3>   ret;
    ret = Matrix<T, Dynamic, 3>::Zero(field.rows(), field.cols());

    for (size_t i = 0; i < field.rows(); i++)
    {
        static_vector<T, 3> f;
        f(0) = field(i,0);
        f(1) = field(i,1);
        f(2) = field(i,2);

        auto r = f.cross( normal );
        ret(i,0) = r(0);
        ret(i,1) = r(1);
        ret(i,2) = r(2);
    }

    return ret;
}

}

/***************************************************************************/
/* RHS definition */
template<typename Mesh>
struct rhs_functor;

/*
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        return 2.0 * M_PI * M_PI * sin_px * sin_py;
    }
};
*/

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    Matrix<T,3,1> operator()(const point_type& pt) const
    {
        Matrix<T, 3, 1> ret;
        //ret(0) = M_PI * M_PI * std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y())*std::sin(M_PI*pt.z());
        //ret(1) = M_PI * M_PI * std::cos(M_PI*pt.x())*std::cos(M_PI*pt.y())*std::sin(M_PI*pt.z());
        //ret(2) = M_PI * M_PI * std::cos(M_PI*pt.x())*std::sin(M_PI*pt.y())*std::cos(M_PI*pt.z());
        ret(0) = 0.0;
        ret(1) = 0.0;
        ret(2) = M_PI*M_PI * std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
        return ret;
    }
};

template<typename Mesh>
auto make_rhs_function(const Mesh& msh)
{
    return rhs_functor<Mesh>();
}

/***************************************************************************/
/* Expected solution definition */
template<typename Mesh>
struct solution_functor;

/*
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        return sin_px * sin_py;
    }
};
*/

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    Matrix<T,3,1> operator()(const point_type& pt) const
    {
        Matrix<T, 3, 1> ret;
        ret(0) = std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y())*std::sin(M_PI*pt.z());
        ret(1) = 0.0;
        ret(2) = 0.0;
        return ret;
    }
};

template<typename Mesh>
auto make_solution_function(const Mesh& msh)
{
    return solution_functor<Mesh>();
}


template<typename Mesh>
void
run_maxwell_solver(Mesh& msh, size_t degree, const typename Mesh::coordinate_type eta,
    const std::string& cfg_fn)
{
    parameter_loader<double> pl;
    bool ok = pl.load(cfg_fn);
    if (!ok)
        return;

    double omega = 2*M_PI*pl.frequency();

    auto cvf = connectivity_via_faces(msh);
    using coordinate_type = typename Mesh::coordinate_type;
    using scalar_type = std::complex<double>;

    typedef Matrix<scalar_type, Dynamic, Dynamic> matrix_type;
    typedef Matrix<scalar_type, Dynamic, 1>       vector_type;

    auto cbs = disk::vector_basis_size(degree, Mesh::dimension, Mesh::dimension);
    auto assm = disk::make_discontinuous_galerkin_assembler<scalar_type>(msh, cbs);


    for (auto& tcl : msh)
    {
        auto di = msh.domain_info(tcl);

        auto eps = pl.epsilon( di.tag() );
        auto mu = pl.mu( di.tag() );
        auto Z = std::sqrt(mu/eps);

        auto tbasis = disk::make_vector_monomial_basis(msh, tcl, degree);
        auto qps = disk::integrate(msh, tcl, 2*degree);
        
        matrix_type M = matrix_type::Zero(tbasis.size(), tbasis.size());
        matrix_type K = matrix_type::Zero(tbasis.size(), tbasis.size());
        vector_type loc_rhs = vector_type::Zero(tbasis.size());
        for (auto& qp : qps)
        {
            auto phi = tbasis.eval_functions(qp.point());
            auto cphi = tbasis.eval_curls2(qp.point());
            
            M += qp.weight() * phi * phi.transpose();
            K += qp.weight() * cphi * cphi.transpose();
            //loc_rhs += qp.weight() * phi * f(qp.point());
        }

        auto fcs = faces(msh, tcl);
        for (auto& fc : fcs)
        {
            auto bi = msh.boundary_info(fc);

            auto f = [&](const disk::point<double,3>& pt) -> Matrix<std::complex<double>,3,1> {
                Matrix<std::complex<double>,3,1> ret;
                auto src = pl.plane_wave_source(bi.tag(), pt);
                ret(0) = std::complex<double>(src.Ex_re, src.Ex_im);
                ret(1) = std::complex<double>(src.Ey_re, src.Ey_im);
                ret(2) = std::complex<double>(src.Ez_re, src.Ez_im);
                return ret;
            };

            matrix_type Att = matrix_type::Zero(tbasis.size(), tbasis.size());
            matrix_type Atn = matrix_type::Zero(tbasis.size(), tbasis.size());
            
            auto nv = cvf.neighbour_via(msh, tcl, fc);
            auto ncl = nv.first;
            auto nbasis = disk::make_vector_monomial_basis(msh, ncl, degree);
            assert(tbasis.size() == nbasis.size());
            
            auto n     = normal(msh, tcl, fc);
            auto eta_l = eta / diameter(msh, fc);
            auto f_qps = disk::integrate(msh, fc, 2*degree);

            for (auto& fqp : f_qps)
            {
                auto tphi       = tbasis.eval_functions(fqp.point());
                auto tcphi      = tbasis.eval_curls2(fqp.point());
                auto n_x_tphi   = disk::vcross(n, tphi);
                auto n_x_tphi_x_n = disk::vcross(n_x_tphi, n);
                
                if (nv.second)
                {   /* NOT on a boundary */
                    Att += + fqp.weight() * eta_l * n_x_tphi * n_x_tphi.transpose();
                    Att += - fqp.weight() * 0.5 * n_x_tphi * tcphi.transpose();
                    Att += - fqp.weight() * 0.5 * tcphi * n_x_tphi.transpose();
                }
                else
                {
                   /* On a boundary*/
                    if ( not pl.is_magnetic_like(bi.tag()) and not pl.is_impedance_like(bi.tag()) )
                    {
                        Att += + fqp.weight() * eta_l * n_x_tphi * n_x_tphi.transpose();
                        Att += - fqp.weight() * tcphi * n_x_tphi.transpose();
                        Att += - fqp.weight() * n_x_tphi * tcphi.transpose();
                    }

                    
                    //loc_rhs -= fqp.weight() * tcphi_x_n;
                    //loc_rhs += fqp.weight() * eta_l * tphi;

                    std::complex<double> jomega(0, omega);

                    if ( pl.is_impedance_like(bi.tag()) )
                        Att -= fqp.weight() * (mu*jomega/Z) * n_x_tphi_x_n * n_x_tphi_x_n.transpose();

                    if ( pl.is_plane_wave(bi.tag()) )
                        loc_rhs += 2*fqp.weight()*(mu*jomega/Z)*n_x_tphi_x_n*f(fqp.point());

                    continue;
                }
                
                auto nphi       = nbasis.eval_functions(fqp.point());
                auto ncphi      = nbasis.eval_curls2(fqp.point());
                auto n_x_nphi   = disk::vcross(n, nphi);
                
                Atn += - fqp.weight() * eta_l * n_x_tphi * n_x_nphi.transpose();
                Atn += - fqp.weight() * 0.5 * n_x_tphi * ncphi.transpose();
                Atn += + fqp.weight() * 0.5 * tcphi * n_x_nphi.transpose();
            }
            
            assm.assemble(msh, tcl, tcl, Att);
            if (nv.second)
                assm.assemble(msh, tcl, ncl, Atn);
        }

        matrix_type LC = K - (eps*omega)*(omega*mu)*M;

        assm.assemble(msh, tcl, LC, loc_rhs);
    }

    assm.finalize();

    //disk::dump_sparse_matrix(assm.LHS, "dg.txt");

    std::cout << "Mesh has " << msh.cells_size() << " elements." << std::endl;
    std::cout << "System has " << assm.LHS.rows() << " unknowns and ";
    std::cout << assm.LHS.nonZeros() << " nonzeros." << std::endl;

    disk::dynamic_vector<scalar_type> sol = disk::dynamic_vector<scalar_type>::Zero(assm.syssz);

    std::cout << "Running pardiso" << std::endl;
    disk::solvers::pardiso_params<scalar_type> pparams;
    pparams.report_factorization_Mflops = true;
    mkl_pardiso(pparams, assm.LHS, assm.RHS, sol);

    /*
    Eigen::GMRES<Eigen::SparseMatrix<T>, Eigen::IdentityPreconditioner> gmres;
    gmres.compute(assm.LHS);
    sol = gmres.solve(assm.RHS);
    */

    std::vector<scalar_type> data_ux, data_uy, data_uz, data_div;

    //T err = 0.0;
    size_t cell_i = 0;
    //T err_divergence = 0.0;
    for (auto& cl : msh)
    {
        auto cb = make_vector_monomial_basis(msh, cl, degree);
        Matrix<scalar_type, Dynamic, 1> lsol = sol.segment(cell_i*cb.size(), cb.size());

        //scalar_type divt = 0.0;
        //auto qps = integrate(msh, cl, degree);
        //for (auto& qp : qps)
        //{
        //    auto divphi = cb.eval_divergences(qp.point());
        //    divt += qp.weight() * (asol-lsol).dot(divphi);
        //}
        //data_div.push_back( divt );
        //err_divergence += divt;

        Matrix<scalar_type,1,3> acc = Matrix<scalar_type,1,3>::Zero();
        auto bar = barycenter(msh, cl);
        for (size_t i = 0; i < cb.size(); i++)
        {
            auto phi = cb.eval_functions(bar);
            acc += lsol(i)*phi.row(i);
        }

        data_ux.push_back( acc(0) );
        data_uy.push_back( acc(1) );
        data_uz.push_back( acc(2) );

        cell_i++;
    }

    //std::cout << "h = " << disk::average_diameter(msh) << " ";
    //std::cout << "err = " << std::sqrt(err) << std::endl;
    //std::cout << "Error divergence: " << err_divergence << std::endl;

    disk::silo_database silo_db;
    silo_db.create("maxwell_sip_dg.silo");
    silo_db.add_mesh(msh, "mesh");

    disk::silo_zonal_variable<scalar_type> ux("ux", data_ux);
    silo_db.add_variable("mesh", ux);

    disk::silo_zonal_variable<scalar_type> uy("uy", data_uy);
    silo_db.add_variable("mesh", uy);

    disk::silo_zonal_variable<scalar_type> uz("uz", data_uz);
    silo_db.add_variable("mesh", uz);

    silo_db.add_expression("u", "{ux, uy, uz}", DB_VARTYPE_VECTOR);
    silo_db.add_expression("mag_u", "magnitude(u)", DB_VARTYPE_SCALAR);

    //disk::silo_zonal_variable<scalar_type> div("div", data_div);
    //silo_db.add_variable("mesh", div);

    silo_db.close();

    return;
}

template<typename Mesh>
void
run_maxwell_eigenvalue_solver(Mesh& msh, size_t degree, const typename Mesh::coordinate_type eta)
{   
    auto cvf = connectivity_via_faces(msh);
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;
    
    auto f = make_rhs_function(msh);

    auto cbs = disk::vector_basis_size(degree, Mesh::dimension, Mesh::dimension);
    auto assm = make_discontinuous_galerkin_eigenvalue_assembler(msh, cbs);

    for (auto& tcl : msh)
    {
        auto tbasis = disk::make_vector_monomial_basis(msh, tcl, degree);
        auto qps = disk::integrate(msh, tcl, 2*degree);
        
        matrix_type M = matrix_type::Zero(tbasis.size(), tbasis.size());
        matrix_type K = matrix_type::Zero(tbasis.size(), tbasis.size());
        vector_type loc_rhs = vector_type::Zero(tbasis.size());
        for (auto& qp : qps)
        {
            auto phi = tbasis.eval_functions(qp.point());
            auto cphi = tbasis.eval_curls2(qp.point());
            
            M += qp.weight() * phi * phi.transpose();
            K += qp.weight() * cphi * cphi.transpose();
        }

        assm.assemble(msh, tcl, tcl, K);
        assm.assemble(msh, tcl, M);

        auto fcs = faces(msh, tcl);
        for (auto& fc : fcs)
        {
            matrix_type Att = matrix_type::Zero(tbasis.size(), tbasis.size());
            matrix_type Atn = matrix_type::Zero(tbasis.size(), tbasis.size());
            
            auto nv = cvf.neighbour_via(msh, tcl, fc);
            auto ncl = nv.first;
            auto nbasis = disk::make_vector_monomial_basis(msh, ncl, degree);
            assert(tbasis.size() == nbasis.size());
            
            auto n     = normal(msh, tcl, fc);
            auto eta_l = eta / diameter(msh, fc);
            auto f_qps = disk::integrate(msh, fc, 2*degree);
            
            for (auto& fqp : f_qps)
            {
                auto tphi       = tbasis.eval_functions(fqp.point());
                auto tcphi      = tbasis.eval_curls2(fqp.point());
                auto n_x_tphi   = disk::vcross(n, tphi);
                
                if (nv.second)
                {   /* NOT on a boundary */
                    Att += + fqp.weight() * eta_l * n_x_tphi * n_x_tphi.transpose();
                    Att += - fqp.weight() * 0.5 * n_x_tphi * tcphi.transpose();
                    Att += - fqp.weight() * 0.5 * tcphi * n_x_tphi.transpose();
                }
                else
                {   /* On a boundary*/
                    Att += + fqp.weight() * eta_l * n_x_tphi * n_x_tphi.transpose();
                    Att += - fqp.weight() * n_x_tphi * tcphi.transpose();
                    Att += - fqp.weight() * tcphi * n_x_tphi.transpose();
                    continue;
                }
                
                auto nphi       = nbasis.eval_functions(fqp.point());
                auto ncphi      = nbasis.eval_curls2(fqp.point());
                auto n_x_nphi   = disk::vcross(n, nphi);
                
                Atn += - fqp.weight() * eta_l * n_x_tphi * n_x_nphi.transpose();
                Atn += - fqp.weight() * 0.5 * n_x_tphi * ncphi.transpose();
                Atn += + fqp.weight() * 0.5 * tcphi * n_x_nphi.transpose();
            }
            
            assm.assemble(msh, tcl, tcl, Att);
            if (nv.second)
                assm.assemble(msh, tcl, ncl, Atn);
        }
        
    }

    assm.finalize();


    disk::feast_eigensolver_params<T> fep;

    fep.verbose = true;
    fep.tolerance = 8;
    fep.min_eigval = 1;
    fep.max_eigval = 100;
    fep.subspace_size = 50;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    dg_eigvecs;
    Eigen::Matrix<T, Eigen::Dynamic, 1>                 dg_eigvals;

    //disk::dump_sparse_matrix(assm.gK, "gK.txt");
    //disk::dump_sparse_matrix(assm.gM, "gM.txt");

    generalized_eigenvalue_solver(fep, assm.gK, assm.gM, dg_eigvecs, dg_eigvals);

    for (size_t i = 0; i < fep.eigvals_found; i++)
    {
        std::vector<T> data_ux, data_uy, data_uz;

        std::stringstream fn;
        fn << "eigenfun_" << i << ".silo";

        disk::silo_database silo_db;
        silo_db.create(fn.str().c_str());
        silo_db.add_mesh(msh, "mesh");

        size_t cell_i = 0;
        for (auto& cl : msh)
        {
            auto cb = make_vector_monomial_basis(msh, cl, degree);
            Matrix<T, Dynamic, 1> lsol = dg_eigvecs.block(cell_i*cb.size(), i, cb.size(), 1);

            Matrix<T,1,3> acc = Matrix<T,1,3>::Zero();
            auto bar = barycenter(msh, cl);
            for (size_t i = 0; i < cb.size(); i++)
            {
                auto phi = cb.eval_functions(bar);
                acc += lsol(i)*phi.row(i);
            }

            data_ux.push_back( acc(0) );
            data_uy.push_back( acc(1) );
            data_uz.push_back( acc(2) );

            cell_i++;
        }

        disk::silo_zonal_variable<T> ux("ux", data_ux);
        silo_db.add_variable("mesh", ux);

        disk::silo_zonal_variable<T> uy("uy", data_uy);
        silo_db.add_variable("mesh", uy);

        disk::silo_zonal_variable<T> uz("uz", data_uz);
        silo_db.add_variable("mesh", uz);
    
        silo_db.add_expression("u", "{ux, uy, uz}", DB_VARTYPE_VECTOR);
        silo_db.add_expression("mag_u", "magnitude(u)", DB_VARTYPE_SCALAR);
    }


    for (size_t i = 0; i < fep.eigvals_found; i++)
    {
        std::cout << std::setprecision(8) << dg_eigvals(i) << " -> ";
        std::cout << std::sqrt( dg_eigvals(i) ) << std::endl;
    }
}

#if 0
void autotest_convergence(size_t order_min, size_t order_max)
{
    using T = double;
    
    using Mesh = disk::simplicial_mesh<T,3>;

    std::ofstream ofs( "dg_convergence_convt.txt" );

    double sp[] = { 50, 50, 50, 50 };

    for (size_t i = 1; i < 5; i++)
    {
        Mesh msh;
        std::stringstream ss;
        ss << "../../../meshes/3D_tetras/netgen/convt0" << i << ".mesh";

        disk::load_mesh_netgen(ss.str().c_str(), msh);
        auto diam = disk::average_diameter(msh);

        ofs << diam << " ";

        for (size_t order = order_min; order <= order_max; order++)
        {
            auto ci = run_maxwell_solver(msh, order, sp[order]);
            ofs << ci.l2_error_e << " " << ci.nrg_error << " " << ci.mflops << " ";
            ofs << ci.dofs << " " << ci.nonzeros << " ";
        }

        ofs << std::endl;
    }
}
#endif

int main(int argc, char **argv)
{
    rusage_monitor rm;

    using T = double;

    T           stab_param = 1.0;
    size_t      degree = 1;
    bool        solve_eigvals = false;
    char *      mesh_filename = nullptr;

    int ch;
    while ( (ch = getopt(argc, argv, "a:k:m:e")) != -1 )
    {
        switch(ch)
        {
            case 'a':
                stab_param = std::stod(optarg);
                break;

            case 'e':
                solve_eigvals = true;
                break;

            case 'k':
                degree = std::stoi(optarg);
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

    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        auto msh = disk::load_netgen_3d_mesh<T>(mesh_filename);
        run_maxwell_solver(msh, degree, stab_param, "params.lua");
        return 0;
    }

    /* GMSH */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo$") ))
    {
        std::cout << "Guessed mesh format: GMSH" << std::endl;
        disk::simplicial_mesh<T,3> msh;
        disk::gmsh_geometry_loader< disk::simplicial_mesh<T,3> > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        /*
        for (auto itor = msh.points_begin(); itor != msh.points_end(); itor++)
        {
            auto pt = *itor;
            auto newpt = disk::point<T,3>(pt.x()/1000, pt.y()/1000, pt.z()/1000);
            *itor = newpt;
        }
        */

        run_maxwell_solver(msh, degree, stab_param, "params.lua");

        return 0;
    }
}

