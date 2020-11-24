/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet (C) 2018                      nicolas.pignet@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */

#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <unistd.h>

#include <map>

#include "colormanip.h"

#include "geometry/geometry.hpp"
#include "loaders/loader.hpp"
#include "methods/hho"
#include "solvers/solver.hpp"

#include "../tests/common.hpp"

/***************************************************************************/
/* RHS definition */
template<typename Mesh>
struct rhs_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor<Mesh<T, 2, Storage>>
{
    typedef Mesh<T, 2, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;

    scalar_type
    operator()(const point_type& pt) const
    {
        const auto sin_px = std::sin(M_PI * pt.x());
        const auto sin_py = std::sin(M_PI * pt.y());
        const auto cos_px = std::cos(M_PI * pt.x());
        const auto cos_py = std::cos(M_PI * pt.y());

        return M_PI * M_PI * (2.0 * sin_px * sin_py - cos_px * cos_py);
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor<Mesh<T, 3, Storage>>
{
    typedef Mesh<T, 3, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;

    scalar_type
    operator()(const point_type& pt) const
    {
        const auto sin_px = std::sin(M_PI * pt.x());
        const auto sin_py = std::sin(M_PI * pt.y());
        const auto sin_pz = std::sin(M_PI * pt.z());
        const auto cos_px = std::cos(M_PI * pt.x());
        const auto cos_py = std::cos(M_PI * pt.y());
        const auto cos_pz = std::cos(M_PI * pt.z());

        return M_PI * M_PI *
               (3.0 * sin_px * sin_py * sin_pz -
                2.0 * (0.5 * cos_px * cos_py * sin_pz + 0.4 * sin_px * cos_py * cos_pz));
    }
};

template<typename Mesh>
auto
make_rhs_function(const Mesh& msh)
{
    return rhs_functor<Mesh>();
}

/***************************************************************************/
/* Expected solution definition */
template<typename Mesh>
struct solution_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_functor<Mesh<T, 2, Storage>>
{
    typedef Mesh<T, 2, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;

    scalar_type
    operator()(const point_type& pt) const
    {
        const auto sin_px = std::sin(M_PI * pt.x());
        const auto sin_py = std::sin(M_PI * pt.y());
        return sin_px * sin_py;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_functor<Mesh<T, 3, Storage>>
{
    typedef Mesh<T, 3, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;

    scalar_type
    operator()(const point_type& pt) const
    {
        const auto sin_px = std::sin(M_PI * pt.x());
        const auto sin_py = std::sin(M_PI * pt.y());
        const auto sin_pz = std::sin(M_PI * pt.z());
        return sin_px * sin_py * sin_pz;
    }
};

template<typename Mesh>
auto
make_solution_function(const Mesh& msh)
{
    return solution_functor<Mesh>();
}

/***************************************************************************/
/* Diffusion tensor definition */
template<typename Mesh>
struct diffusion_tensor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct diffusion_tensor<Mesh<T, 2, Storage>>
{
    typedef Mesh<T, 2, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;

    typedef disk::static_matrix<scalar_type, mesh_type::dimension, mesh_type::dimension> tensor_type;

    tensor_type
    operator()(const point_type& pt) const
    {
        tensor_type ret = tensor_type::Identity();
        ret(0, 1)       = 0.5;
        ret(1, 0)       = 0.5;
        return ret;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct diffusion_tensor<Mesh<T, 3, Storage>>
{
    typedef Mesh<T, 3, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;

    typedef disk::static_matrix<scalar_type, mesh_type::dimension, mesh_type::dimension> tensor_type;

    tensor_type
    operator()(const point_type& pt) const
    {
        tensor_type ret = tensor_type::Identity();
        ret(0, 1)       = 0.5;
        ret(1, 0)       = 0.5;
        ret(1, 2)       = 0.4;
        ret(2, 1)       = 0.4;
        return ret;
    }
};

template<typename Mesh>
auto
make_diffusion_tensor(const Mesh& msh)
{
    return diffusion_tensor<Mesh>();
}

using namespace disk;

void print_error(const size_t& degree, const double diam, const double H1_error)
{
    std::cout << "Degree: " << degree << std::endl;
    std::cout << "h-diameter: " << diam << std::endl;
    std::cout << "H1-error: " << H1_error << std::endl;
}

/* Solve anisotropic diffusion with HHO method on simplicial meshes using Raviart-Thomas */

template<typename Mesh>
typename Mesh::coordinate_type
run_hho_variable_diffusion_solver_RT(const Mesh& msh, const size_t degree, bool print = false)
{
    using T = typename Mesh::coordinate_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1> vector_type;

    hho_degree_info hdi(degree, degree, degree + 1);

    const auto rhs_fun          = make_rhs_function(msh);
    const auto sol_fun          = make_solution_function(msh);
    const auto diffusion_tensor = make_diffusion_tensor(msh);

    scalar_boundary_conditions<Mesh> bnd(msh);
    bnd.addDirichletEverywhere(sol_fun);

    auto assembler = make_scalar_primal_hho_assembler(msh, hdi, bnd);

    for (auto& cl : msh)
    {
        const auto cb  = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        const auto gr  = make_vector_hho_gradrec_RT(msh, cl, hdi, diffusion_tensor);
        const auto rhs = make_rhs(msh, cl, cb, rhs_fun, 2);

        const auto sc = make_scalar_static_condensation(msh, cl, hdi, gr.second, rhs);

        assembler.assemble(msh, cl, bnd, sc.first, sc.second);
    }

    assembler.impose_neumann_boundary_conditions(msh, bnd);

    assembler.finalize();

    size_t systsz = assembler.LHS.rows();
    size_t nnz    = assembler.LHS.nonZeros();

    //  std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
    //  std::cout << "systsz: " << systsz << std::endl;
    //  std::cout << "nnz: " << nnz << std::endl;

    vector_type sol = vector_type::Zero(systsz);

    disk::solvers::pardiso_params<T> pparams;
    pparams.report_factorization_Mflops = false;
    mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

    T error = 0.0;

    for (auto& cl : msh)
    {
        const auto cb  = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        const auto gr  = make_vector_hho_gradrec_RT(msh, cl, hdi, diffusion_tensor);
        const auto rhs = make_rhs(msh, cl, cb, rhs_fun);

        vector_type locsol = assembler.take_local_solution(msh, cl, bnd, sol);

        vector_type sol = make_scalar_static_decondensation(msh, cl, hdi, gr.second, rhs, locsol);

        vector_type realsol = project_function(msh, cl, hdi, sol_fun, 2);

        const auto diff = realsol - sol;
        error += diff.dot(gr.second * diff);
    }

    if(print)
    {
        print_error(degree, average_diameter(msh), error);
    }

    return std::sqrt(error);
}

template<typename Mesh>
typename Mesh::coordinate_type
run_hho_variable_diffusion_solver(const Mesh& msh, const size_t degree, bool print = false)
{
    using T = typename Mesh::coordinate_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              vector_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_type;

    hho_degree_info hdi(degree + 1, degree, degree);

    const auto rhs_fun          = make_rhs_function(msh);
    const auto sol_fun          = make_solution_function(msh);
    const auto diffusion_tensor = make_diffusion_tensor(msh);

    scalar_boundary_conditions<Mesh> bnd(msh);
    bnd.addDirichletEverywhere(sol_fun);

    auto assembler = make_scalar_primal_hho_assembler(msh, hdi, bnd);

    for (auto& cl : msh)
    {
        const auto gr   = make_vector_hho_gradrec(msh, cl, hdi, diffusion_tensor);
        const auto stab = make_scalar_hdg_stabilization(msh, cl, hdi);

        matrix_type A = gr.second + stab;

        const auto cb  = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        const auto rhs = make_rhs(msh, cl, cb, rhs_fun, 2);

        const auto sc = make_scalar_static_condensation(msh, cl, hdi, A, rhs);

        assembler.assemble(msh, cl, bnd, sc.first, sc.second);
    }

    assembler.impose_neumann_boundary_conditions(msh, bnd);

    assembler.finalize();

    size_t systsz = assembler.LHS.rows();
    size_t nnz    = assembler.LHS.nonZeros();

    //  std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
    //  std::cout << "systsz: " << systsz << std::endl;
    //  std::cout << "nnz: " << nnz << std::endl;

    vector_type sol = vector_type::Zero(systsz);

    disk::solvers::pardiso_params<T> pparams;
    pparams.report_factorization_Mflops = false;
    mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

    T error = 0.0;

    for (auto& cl : msh)
    {
        const auto gr   = make_vector_hho_gradrec(msh, cl, hdi, diffusion_tensor);
        const auto stab = make_scalar_hdg_stabilization(msh, cl, hdi);

        matrix_type A = gr.second + stab;

        const auto cb  = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        const auto rhs = make_rhs(msh, cl, cb, rhs_fun);

        vector_type locsol = assembler.take_local_solution(msh, cl, bnd, sol);

        vector_type sol = make_scalar_static_decondensation(msh, cl, hdi, A, rhs, locsol);

        vector_type realsol = project_function(msh, cl, hdi, sol_fun, 2);

        const auto diff = realsol - sol;
        error += diff.dot(gr.second * diff);
    }

    if (print)
    {
        print_error(degree, average_diameter(msh), error);
    }

    return std::sqrt(error);
}

template<typename Mesh>
class test_functor
{
  public:
    /* Expect k+1 convergence (hho stabilization, energy norm) */
    typename Mesh::coordinate_type
    operator()(const Mesh& msh, size_t degree) const
    {
        return run_hho_variable_diffusion_solver(msh, degree, false);
    }

    size_t
    expected_rate(size_t k)
    {
        return k + 1;
    }
};

template<typename Mesh>
class test_functor_RT
{
  public:
    /* Expect k+1 convergence (energy norm) */
    typename Mesh::coordinate_type
    operator()(const Mesh& msh, size_t degree) const
    {
        return run_hho_variable_diffusion_solver_RT(msh, degree, false);
    }

    size_t
    expected_rate(size_t k)
    {
        return k + 1;
    }
};

void
print_info()
{
    std::cout << "Arguments" << std::endl;
    std::cout << "-m : use the specified mesh file" << std::endl;
    std::cout << "-k : degree of the HHO method" << std::endl;
    std::cout << "-r : use Raiart-Thomas gradient" << std::endl;
    std::cout << "default: convergence test" << std::endl;
}

int
main(int argc, char** argv)
{
    using RealType = double;

    char* mesh_filename = nullptr;

    bool   use_mesh = false;
    bool   use_RT   = false;
    size_t degree   = 1;
    int    ch;

    while ((ch = getopt(argc, argv, "k:m:r")) != -1)
    {
        switch (ch)
        {
            case 'k': degree = std::stoi(optarg); break;

            case 'm':
                mesh_filename = optarg;
                use_mesh      = true;
                break;

            case 'r': use_RT = true; break;

            default:
                std::cout << "wrong arguments" << std::endl;
                print_info();
                exit(1);
        }
    }

    if (use_mesh)
    {
        /* FVCA5 2D */
        if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$")))
        {
            std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
            disk::generic_mesh<RealType, 2> msh;
            disk::load_mesh_fvca5_2d(mesh_filename, msh);

            if (use_RT)
            {
                run_hho_variable_diffusion_solver_RT(msh, degree, true);
            }
            else
            {
                run_hho_variable_diffusion_solver(msh, degree, true);
            }

            return 0;
        }

        /* Netgen 2D */
        if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$")))
        {
            std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
            disk::simplicial_mesh<RealType, 2> msh;
            disk::load_mesh_netgen(mesh_filename, msh);

            if (use_RT)
            {
                run_hho_variable_diffusion_solver_RT(msh, degree, true);
            }
            else
            {
                run_hho_variable_diffusion_solver(msh, degree, true);
            }

            return 0;
        }

        /* DiSk++ cartesian 2D */
        if (std::regex_match(mesh_filename, std::regex(".*\\.quad$")))
        {
            std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
            disk::cartesian_mesh<RealType, 2> msh;
            disk::load_mesh_diskpp_cartesian(mesh_filename, msh);

            if (use_RT)
            {
                run_hho_variable_diffusion_solver_RT(msh, degree, true);
            }
            else
            {
                run_hho_variable_diffusion_solver(msh, degree, true);
            }
            return 0;
        }

        /* Netgen 3D */
        if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$")))
        {
            std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
            disk::simplicial_mesh<RealType, 3> msh;
            disk::load_mesh_netgen(mesh_filename, msh);

            if (use_RT)
            {
                run_hho_variable_diffusion_solver_RT(msh, degree, true);
            }
            else
            {
                run_hho_variable_diffusion_solver(msh, degree, true);
            }

            return 0;
        }

        /* DiSk++ cartesian 3D */
        if (std::regex_match(mesh_filename, std::regex(".*\\.hex$")))
        {
            std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
            disk::cartesian_mesh<RealType, 3> msh;
            disk::load_mesh_diskpp_cartesian(mesh_filename, msh);

            if (use_RT)
            {
                run_hho_variable_diffusion_solver_RT(msh, degree, true);
            }
            else
            {
                run_hho_variable_diffusion_solver(msh, degree, true);
            }
            return 0;
        }

        /* FVCA6 3D */
        if (std::regex_match(mesh_filename, std::regex(".*\\.msh$")))
        {
            std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
            disk::generic_mesh<RealType, 3> msh;
            disk::load_mesh_fvca6_3d(mesh_filename, msh);

            if (use_RT)
            {
                run_hho_variable_diffusion_solver_RT(msh, degree, true);
            }
            else
            {
                run_hho_variable_diffusion_solver(msh, degree, true);
            }

            return 0;
        }
    }
    else
    {
        if (use_RT)
        {
            tester_simplicial<test_functor_RT> tstr;
            tstr.run();
        }
        else
        {
            tester<test_functor> tstr;
            tstr.run();
        }
    }

    return 0;
}