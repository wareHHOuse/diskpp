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

#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

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
        ret(0) = 0.0;
        ret(1) = 0.0;
        ret(2) = M_PI * M_PI * std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
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
        ret(0) = 0.0;//M_PI*std::sin(M_PI*pt.x());
        ret(1) = 0.0;
        ret(2) = std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
        return ret;
    }
};

template<typename Mesh>
auto make_solution_function(const Mesh& msh)
{
    return solution_functor<Mesh>();
}


template<typename Mesh>
class discontinuous_galerkin_assembler
{
    size_t cbs;

    using mesh_type = Mesh;
    using cell_type = typename mesh_type::cell_type;
    using scal_type = typename mesh_type::coordinate_type;
    using matrix_type = Matrix<scal_type, Dynamic, Dynamic>;
    using vector_type = Matrix<scal_type, Dynamic, 1>;
    using triplet_type = Triplet<scal_type>;
    using T = scal_type;

    std::vector<triplet_type> triplets;

public:
    SparseMatrix<T>         LHS;
    Matrix<T, Dynamic, 1>   RHS;

    size_t                  syssz;

    discontinuous_galerkin_assembler(const Mesh& msh, size_t pcbs)
    {
        cbs = pcbs;
        syssz = cbs * msh.cells_size();

        LHS = SparseMatrix<T>(syssz, syssz);
        RHS = Matrix<T, Dynamic, 1>::Zero(syssz);
    }

    void
    assemble(const mesh_type& msh, const cell_type& clp, const cell_type& clm,
             const matrix_type& A)
    {
        auto ofs_i = cbs * offset(msh, clp);
        auto ofs_j = cbs * offset(msh, clm);

        for (size_t i = 0; i < cbs; i++)
            for (size_t j = 0; j < cbs; j++)
                triplets.push_back( triplet_type(ofs_i+i, ofs_j+j, A(i,j)) );
    }

    void
    assemble(const mesh_type& msh, const cell_type& cl, const matrix_type& A,
             const vector_type& b)
    {
        auto ofs = cbs * offset(msh, cl);

        for (size_t i = 0; i < cbs; i++)
            for (size_t j = 0; j < cbs; j++)
                triplets.push_back( triplet_type(ofs+i, ofs+j, A(i,j)) );

        RHS.segment(ofs, cbs) += b;
    }

    void
    finalize()
    {
        LHS.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();
    }
};

template<typename Mesh>
auto
make_discontinuous_galerkin_assembler(const Mesh& msh, size_t cbs)
{
    return discontinuous_galerkin_assembler<Mesh>(msh, cbs);
}

template<typename Mesh>
void
run_maxwell_solver(Mesh& msh, size_t degree, const typename Mesh::coordinate_type eta)
{   
    auto cvf = connectivity_via_faces(msh);
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;
    
    auto f = make_rhs_function(msh);

    auto cbs = disk::vector_basis_size(degree, Mesh::dimension, Mesh::dimension);
    auto assm = make_discontinuous_galerkin_assembler(msh, cbs);

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
            loc_rhs += qp.weight() * phi * f(qp.point());
        }

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
                    
                    //loc_rhs -= fqp.weight() * tcphi_x_n;
                    //loc_rhs += fqp.weight() * eta_l * tphi;
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
        
        auto omega = M_PI;
        matrix_type LC = K - omega*omega*M;
        assm.assemble(msh, tcl, LC, loc_rhs);
    }

    assm.finalize();

    //disk::dump_sparse_matrix(assm.LHS, "dg.txt");

    std::cout << "Mesh has " << msh.cells_size() << " elements." << std::endl;
    std::cout << "System has " << assm.LHS.rows() << " unknowns and ";
    std::cout << assm.LHS.nonZeros() << " nonzeros." << std::endl;

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.syssz);

    std::cout << "Running pardiso" << std::endl;
    disk::solvers::pardiso_params<T> pparams;
    pparams.report_factorization_Mflops = true;
    mkl_pardiso(pparams, assm.LHS, assm.RHS, sol);

    /*
    Eigen::GMRES<Eigen::SparseMatrix<T>, Eigen::IdentityPreconditioner> gmres;
    gmres.compute(assm.LHS);
    sol = gmres.solve(assm.RHS);
    */

    std::vector<T> data_ux, data_uy, data_uz;

    auto sol_fun = make_solution_function(msh);

    T err = 0.0; size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto cb = make_vector_monomial_basis(msh, cl, degree);

        Matrix<T, Dynamic, Dynamic> MMe = disk::make_mass_matrix(msh, cl, cb, 1);
        Matrix<T, Dynamic, 1> arhs = disk::make_rhs(msh, cl, cb, sol_fun, 1);
        Matrix<T, Dynamic, 1> asol = MMe.llt().solve(arhs);
        Matrix<T, Dynamic, 1> lsol = sol.segment(cell_i*cb.size(), cb.size());

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

        Matrix<T, Dynamic, 1> diff = lsol - asol;

        err += diff.dot(MMe*diff);

        cell_i++;
    }

    std::cout << "h = " << disk::average_diameter(msh) << " ";
    std::cout << "err = " << std::sqrt(err) << std::endl;

    disk::silo_database silo_db;
    silo_db.create("maxwell_sip_dg.silo");
    silo_db.add_mesh(msh, "mesh");

    disk::silo_zonal_variable<T> ux("ux", data_ux);
    silo_db.add_variable("mesh", ux);

    disk::silo_zonal_variable<T> uy("uy", data_uy);
    silo_db.add_variable("mesh", uy);

    disk::silo_zonal_variable<T> uz("uz", data_uz);
    silo_db.add_variable("mesh", uz);

    silo_db.add_expression("u", "{ux, uy, uz}", DB_VARTYPE_VECTOR);
    silo_db.add_expression("mag_u", "magnitude(u)", DB_VARTYPE_SCALAR);

    silo_db.close();
}

int main(int argc, char **argv)
{
    using T = double;

    T           stab_param = 1.0;
    size_t      degree = 1;
    char *      mesh_filename = nullptr;

    int ch;
    while ( (ch = getopt(argc, argv, "a:k:m:")) != -1 )
    {
        switch(ch)
        {
            case 'a':
                stab_param = std::stod(optarg);
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
        run_maxwell_solver(msh, degree, stab_param);
        return 0;
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        auto msh = disk::load_cartesian_3d_mesh<T>(mesh_filename);
        run_maxwell_solver(msh, degree, stab_param);
        return 0;
    }
}

