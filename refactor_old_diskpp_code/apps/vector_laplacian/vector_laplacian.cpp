/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
 * Nicolas Pignet  (C) 2018                     nicolas.pignet@enpc.fr
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

#include "boundary_conditions/boundary_conditions.hpp"
#include "geometry/geometry.hpp"
#include "loaders/loader.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "vector_laplacian_solver.hpp"

struct run_params
{
    size_t degree;
    int    l;
    bool   verbose;
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_vector_laplacian_solver(const Mesh<T, 2, Storage>& msh, const run_params& rp, const Parameters material_data)
{
    typedef Mesh<T, 2, Storage>                            mesh_type;
    typedef disk::static_vector<T, 2>                     result_type;
    typedef disk::static_matrix<T, 2, 2>                  result_grad_type;
    typedef disk::vector_boundary_conditions<mesh_type>    Bnd_type;

    timecounter tc;
    tc.tic();

    auto load = [material_data](const disk::point<T, 2>& p) -> result_type {
        T fx = 2. * material_data.lambda * M_PI * M_PI * sin(M_PI * p.x()) * sin(M_PI * p.y());
        T fy = 2. * material_data.lambda * M_PI * M_PI * cos(M_PI * p.x()) * cos(M_PI * p.y());

        return result_type{fx, fy};
    };

    auto solution = [material_data](const disk::point<T, 2>& p) -> result_type {
        T fx = sin(M_PI * p.x()) * sin(M_PI * p.y());
        T fy = cos(M_PI * p.x()) * cos(M_PI * p.y());

        return result_type{fx, fy};
    };

    auto gradient = [material_data](const disk::point<T, 2>& p) -> result_grad_type {
        result_grad_type grad = result_grad_type::Zero();

        grad(0, 0) = cos(M_PI * p.x()) * sin(M_PI * p.y());
        grad(0, 1) = sin(M_PI * p.x()) * cos(M_PI * p.y());
        grad(1, 0) = -sin(M_PI * p.x()) * cos(M_PI * p.y());
        grad(1, 1) = -cos(M_PI * p.x()) * sin(M_PI * p.y());
        return M_PI * grad;
    };

    Bnd_type bnd(msh);
    bnd.addDirichletEverywhere(solution);

    vector_laplacian_solver<mesh_type> le(msh, bnd, material_data, rp.degree, rp.l);
    le.verbose(rp.verbose);

    assembly_info assembling_info = le.assemble(load);

    solver_info solve_info = le.solve();

    postprocess_info post_info = le.postprocess(load);

    tc.toc();

    if (le.verbose())
    {
        std::cout << std::endl;
        std::cout << "************************************************************" << std::endl;
        std::cout << "** Time to solve the problem " << tc.to_double() << " sec" << std::endl;
        std::cout << "**** Assembly time: " << assembling_info.time_assembly << " sec" << std::endl;
        std::cout << "****** Gradient reconstruction: " << assembling_info.time_gradrec << " sec" << std::endl;
        std::cout << "****** Divergence reconstruction: " << assembling_info.time_divrec << " sec" << std::endl;
        std::cout << "****** Stabilisation: " << assembling_info.time_stab << " sec" << std::endl;
        std::cout << "****** Static condensation: " << assembling_info.time_statcond << " sec" << std::endl;
        std::cout << "**** Solver time: " << solve_info.time_solver << " sec" << std::endl;
        std::cout << "**** Postprocess time: " << post_info.time_postprocess << " sec" << std::endl;
        std::cout << "***********************************************************" << std::endl;
    }

    std::cout << "Discetisation h: " << disk::average_diameter(msh) << std::endl;
    std::cout << "L2 error: " << le.compute_l2_displacement_error(solution) << std::endl;
    std::cout << "L2 grad error: " << le.compute_l2_gradient_error(gradient) << std::endl;

    le.compute_continuous_solution("sol2d.msh");
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_vector_laplacian_solver(const Mesh<T, 3, Storage>& msh, const run_params& rp, const Parameters material_data)
{
    typedef Mesh<T, 3, Storage>                            mesh_type;
    typedef disk::static_vector<T, 3>                     result_type;
    typedef disk::static_matrix<T, 3, 3>                  result_grad_type;
    typedef disk::vector_boundary_conditions<mesh_type>    Bnd_type;

    timecounter tc;
    tc.tic();

    auto load = [material_data](const disk::point<T, 3>& p) -> auto
    {
        T fx = 2. * material_data.lambda * M_PI * M_PI * cos(M_PI * p.x()) * sin(M_PI * p.y());

        T fy = 2. * material_data.lambda * M_PI * M_PI * cos(M_PI * p.y()) * sin(M_PI * p.z());

        T fz = 2. * material_data.lambda * M_PI * M_PI * cos(M_PI * p.z()) * sin(M_PI * p.x());

        return result_type{fx, fy, fz};
    };

    auto solution = [material_data](const disk::point<T, 3>& p) -> auto
    {
        T fx = cos(M_PI * p.x()) * sin(M_PI * p.y());
        T fy = cos(M_PI * p.y()) * sin(M_PI * p.z());
        T fz = cos(M_PI * p.z()) * sin(M_PI * p.x());
        return result_type{fx, fy, fz};
    };

    auto gradient = [material_data](const disk::point<T, 3>& p) -> result_grad_type {
        result_grad_type grad = result_grad_type::Zero();

        grad(0, 0) = -sin(M_PI * p.x()) * sin(M_PI * p.y());
        grad(0, 1) = cos(M_PI * p.x()) * cos(M_PI * p.y());
        grad(1, 1) = -sin(M_PI * p.y()) * sin(M_PI * p.z());
        grad(1, 2) = cos(M_PI * p.y()) * cos(M_PI * p.z());
        grad(2, 2) = -sin(M_PI * p.z()) * sin(M_PI * p.x());
        grad(2, 0) = cos(M_PI * p.z()) * cos(M_PI * p.x());
        return M_PI * grad;
    };

    Bnd_type bnd(msh);
    bnd.addDirichletEverywhere(solution);

    vector_laplacian_solver<mesh_type> le(msh, bnd, material_data, rp.degree, rp.l);
    le.verbose(rp.verbose);

    assembly_info assembling_info = le.assemble(load);

    solver_info solve_info = le.solve();

    postprocess_info post_info = le.postprocess(load);

    tc.toc();

    if (le.verbose())
    {
        std::cout << std::endl;
        std::cout << "***********************************************************" << std::endl;
        std::cout << "** Time to solve the problem " << tc.to_double() << " sec" << std::endl;
        std::cout << "**** Assembly time: " << assembling_info.time_assembly << " sec" << std::endl;
        std::cout << "****** Gradient reconstruction: " << assembling_info.time_gradrec << " sec" << std::endl;
        std::cout << "****** Divergence reconstruction: " << assembling_info.time_divrec << " sec" << std::endl;
        std::cout << "****** Stabilisation: " << assembling_info.time_stab << " sec" << std::endl;
        std::cout << "****** Static condensation: " << assembling_info.time_statcond << " sec" << std::endl;
        std::cout << "**** Solver time: " << solve_info.time_solver << " sec" << std::endl;
        std::cout << "**** Postprocess time: " << post_info.time_postprocess << " sec" << std::endl;
        std::cout << "***********************************************************" << std::endl;
    }

    std::cout << "Discetisation h: " << disk::average_diameter(msh) << std::endl;
    std::cout << "L2 error: " << le.compute_l2_displacement_error(solution) << std::endl;
    std::cout << "L2 grad error: " << le.compute_l2_gradient_error(gradient) << std::endl;

    le.compute_continuous_solution("sol3d.msh");
}

int
main(int argc, char** argv)
{
    using RealType = double;

    char* mesh_filename = nullptr;
    int   degree        = 1;
    int   l             = 0;

    run_params rp;
    rp.degree  = 1;
    rp.l       = 0;
    rp.verbose = true;

    // Parameters
    Parameters material_data;

    material_data.lambda = 1;

    int ch;

    while ((ch = getopt(argc, argv, "k:l:v:")) != -1)
    {
        switch (ch)
        {
            case 'k':
                degree = atoi(optarg);
                if (degree < 0)
                {
                    std::cout << "Degree must be positive. Falling back to 1." << std::endl;
                    degree = 1;
                }
                rp.degree = degree;
                break;

            case 'l':
                rp.l = atoi(optarg);
                if (rp.l < -1 or rp.l > 1)
                {
                    std::cout << "l can be -1, 0 or 1. Falling back to 0." << std::endl;
                    rp.l = 0;
                }
                break;

            case 'v': rp.verbose = true; break;

            case '?':
            default: std::cout << "wrong arguments" << std::endl; exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    mesh_filename = argv[0];

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$")))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        auto msh = disk::load_fvca5_2d_mesh<RealType>(mesh_filename);
        run_vector_laplacian_solver(msh, rp, material_data);
        return 0;
    }

    /* Medit 2d*/
    if (std::regex_match(mesh_filename, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        auto msh = disk::load_medit_2d_mesh<RealType>(mesh_filename);
        run_vector_laplacian_solver(msh, rp, material_data);
        return 0;
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$")))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);
        run_vector_laplacian_solver(msh, rp, material_data);
        return 0;
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$")))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        auto msh = disk::load_cartesian_2d_mesh<RealType>(mesh_filename);
        run_vector_laplacian_solver(msh, rp, material_data);
        return 0;
    }

    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$")))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        auto msh = disk::load_netgen_3d_mesh<RealType>(mesh_filename);
        run_vector_laplacian_solver(msh, rp, material_data);
        return 0;
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$")))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        auto msh = disk::load_cartesian_3d_mesh<RealType>(mesh_filename);
        run_vector_laplacian_solver(msh, rp, material_data);
        return 0;
    }

    std::cout << "Unkwnon mesh format" << std::endl;
    return 0;
}
