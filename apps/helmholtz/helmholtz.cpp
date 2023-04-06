#include <iostream>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include "diskpp/common/colormanip.h"

#include "diskpp/methods/hho"
#include "diskpp/methods/implementation_hho/curl.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/output/silo.hpp"
#include "diskpp/solvers/solver.hpp"

/***************************************************************************/
/* RHS definition */
template<typename Mesh>
struct rhs_functor;

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
        return M_PI * M_PI * sin_px * sin_py;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        auto sin_pz = std::sin(M_PI * pt.z());
        return 2 * M_PI * M_PI * sin_px * sin_py * sin_pz;
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

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        auto sin_pz = std::sin(M_PI * pt.z());
        return sin_px * sin_py * sin_pz;
    }
};

template<typename Mesh>
auto make_solution_function(const Mesh& msh)
{
    return solution_functor<Mesh>();
}

using namespace disk;

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_mass_oper(const Mesh&                       msh,
                      const typename Mesh::cell_type&   cl,
                      const hho_degree_info&            cell_infos)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto facdeg = cell_infos.face_degree();
    const auto celdeg = cell_infos.cell_degree();

    const auto cb = make_scalar_monomial_basis(msh, cl, celdeg);
    const auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);
    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);

    const auto fcs = faces(msh, cl);
    const auto num_faces_dofs = fcs.size() * fbs;

    matrix_type mass = matrix_type::Zero(cbs + num_faces_dofs, cbs + num_faces_dofs);

    const auto qps = integrate(msh, cl, 2*celdeg);
    for (auto& qp : qps)
    {
        const auto phi = cb.eval_functions(qp.point());
        mass.block(0,0,cbs,cbs) += qp.weight() * phi * phi.transpose();
    }

    return mass;
}

template<typename Mesh>
void
run_hho_helmholtz_solver(Mesh& msh, size_t degree, const typename Mesh::coordinate_type& alpha)
{
    using T = typename Mesh::coordinate_type;

    hho_degree_info hdi(degree);

    auto rhs_fun = make_rhs_function(msh);
    auto sol_fun = make_solution_function(msh);

    auto assembler = make_diffusion_assembler(msh, hdi);

    T omega = M_PI;


    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
        auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);
        auto mm     = make_scalar_mass_oper(msh, cl, hdi);
        auto rhs    = make_rhs(msh, cl, cb, rhs_fun);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = gr.second + alpha*stab - omega*omega*mm;
        auto sc     = make_scalar_static_condensation(msh, cl, hdi, A, rhs);
        assembler.assemble(msh, cl, sc.first, sc.second, sol_fun);
    }

    assembler.finalize();

    size_t systsz = assembler.LHS.rows();
    size_t nnz = assembler.LHS.nonZeros();

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(systsz);

    disk::solvers::pardiso_params<T> pparams;
    pparams.report_factorization_Mflops = true;
    mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

    T l2_error = 0.0;
    T nrg_error = 0.0;

    std::vector<T> data_u;

    for (auto& cl : msh)
    {
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
        auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);
        auto mm     = make_scalar_mass_oper(msh, cl, hdi);
        auto rhs    = make_rhs(msh, cl, cb, rhs_fun);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = gr.second + alpha*stab - omega*omega*mm;

        Eigen::Matrix<T, Eigen::Dynamic, 1> locsol =
            assembler.take_local_data(msh, cl, sol, sol_fun);

        Eigen::Matrix<T, Eigen::Dynamic, 1> fullsol = make_scalar_static_decondensation(msh, cl, hdi, A, rhs, locsol);

        Eigen::Matrix<T, Eigen::Dynamic, 1> realsol = project_function(msh, cl, hdi, sol_fun, 2);

        data_u.push_back( fullsol(0) );

        auto diff = realsol - fullsol;
        l2_error += diff.dot(mm*diff);
        nrg_error += diff.dot(A*diff);

    }

    std::cout << "h = " << disk::average_diameter(msh) << std::endl;
    std::cout << "L2 error:     " << std::sqrt(l2_error) << std::endl; 
    std::cout << "Energy error: " << std::sqrt(nrg_error) << std::endl; 

    disk::silo_database silo_db;
    silo_db.create("helmholtz.silo");
    silo_db.add_mesh(msh, "mesh");

    disk::silo_zonal_variable<T> u("u", data_u);
    silo_db.add_variable("mesh", u);

    silo_db.close();
}

int main(int argc, char **argv)
{
    using T = double;

    T           stab_param = 1.0;
    bool        solve_eigvals = false;
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

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<T>(mesh_filename);

        run_hho_helmholtz_solver(msh, degree, stab_param);

        return 0;
    }

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        auto msh = disk::load_fvca5_2d_mesh<T>(mesh_filename);

        run_hho_helmholtz_solver(msh, degree, stab_param);

        return 0;
    }

    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        auto msh = disk::load_netgen_3d_mesh<T>(mesh_filename);

        run_hho_helmholtz_solver(msh, degree, stab_param);

        return 0;
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        auto msh = disk::load_cartesian_3d_mesh<T>(mesh_filename);
        
        run_hho_helmholtz_solver(msh, degree, stab_param);
        
        return 0;
    }

    /* FVCA6 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
        disk::generic_mesh<T,3> msh;
        
        disk::load_mesh_fvca6_3d<T>(mesh_filename, msh);
        
        run_hho_helmholtz_solver(msh, degree, stab_param);
        
        return 0;
    }
}