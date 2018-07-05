/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
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

#include <map>

#include "Informations.hpp"
#include "Parameters.hpp"
#include "loaders/loader.hpp"
#include "mechanics/BoundaryConditions.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "plasticity_solver.hpp"

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_plasticity_solver(const Mesh<T, 2, Storage>&        msh,
                      const ParamRun<T>&                rp,
                      const NLE::MaterialParameters<T>& material_data)
{
    typedef Mesh<T, 2, Storage>                            mesh_type;
    typedef static_vector<T, 2>                            result_type;
    typedef static_matrix<T, 2, 2>                         result_grad_type;
    typedef disk::mechanics::BoundaryConditions<mesh_type> Bnd_type;

    auto load = [material_data](const point<T, 2>& p) -> result_type { return result_type{0, 0}; };

    auto solution = [material_data](const point<T, 2>& p) -> result_type { return result_type{0, 0}; };

    auto gradient = [material_data](const point<T, 2>& p) -> result_grad_type { return result_grad_type::Zero(); };

    Bnd_type bnd(msh);

    // Plaque with quadrilaterals
    // auto zero = [material_data](const point<T, 2>& p) -> result_type { return result_type{0.0, 0}; };

    // auto trac = [material_data](const point<T, 2>& p) -> result_type { return result_type{0.0, 5}; };

    // bnd.addDirichletBC(disk::mechanics::DX, 4, zero);
    // bnd.addDirichletBC(disk::mechanics::DX, 14, zero);
    // bnd.addDirichletBC(disk::mechanics::DY, 23, zero);
    // bnd.addDirichletBC(disk::mechanics::DY, 16, trac);

    // Cook with quadrilaterals
    auto zero = [material_data](const point<T, 2>& p) -> result_type { return result_type{0.0, 0}; };

    auto trac = [material_data](const point<T, 2>& p) -> result_type { return result_type{0.0, 0.1125}; };
    auto depltest = [material_data](const point<T, 2>& p) -> result_type { return result_type{0.0, 10}; };

    bnd.addDirichletBC(disk::mechanics::CLAMPED, 3, zero);
    //bnd.addNeumannBC(disk::mechanics::NEUMANN, 8, trac);
    bnd.addDirichletBC(disk::mechanics::DY, 8, depltest);

    plasticity_solver<mesh_type> nl(msh, bnd, rp, material_data);

    if (nl.verbose())
    {
        std::cout << "Solving the problem ..." << '\n';
    }

    SolverInfo solve_info = nl.compute(load);

    if (nl.verbose())
    {
        std::cout << " " << std::endl;
        std::cout << "------------------------------------------------------- " << std::endl;
        std::cout << "Summaring: " << std::endl;
        std::cout << "Total Newton's iterations: " << solve_info.m_iter << " in " << solve_info.m_time_step
                  << " load increments" << std::endl;
        std::cout << "Total time to solve the problem: " << solve_info.m_time_solver << " sec" << std::endl;
        std::cout << "**** Assembly time: " << solve_info.m_newton_info.m_assembly_info.m_time_assembly << " sec"
                  << std::endl;
        std::cout << "****** Gradient reconstruction: " << solve_info.m_newton_info.m_assembly_info.m_time_gradrec
                  << " sec" << std::endl;
        std::cout << "****** Stabilisation: " << solve_info.m_newton_info.m_assembly_info.m_time_stab << " sec"
                  << std::endl;
        std::cout << "****** Elementary computation: " << solve_info.m_newton_info.m_assembly_info.m_time_elem << " sec"
                  << std::endl;
        std::cout << "       *** Behavior computation: " << solve_info.m_newton_info.m_assembly_info.m_time_law
                  << " sec" << std::endl;
        std::cout << "****** Static condensation: " << solve_info.m_newton_info.m_assembly_info.m_time_statcond
                  << " sec" << std::endl;
        std::cout << "**** Postprocess time: " << solve_info.m_newton_info.m_assembly_info.m_time_postpro << " sec"
                  << std::endl;
        std::cout << "**** Solver time: " << solve_info.m_newton_info.m_solve_info.m_time_solve << " sec" << std::endl;
        std::cout << "------------------------------------------------------- " << std::endl;
        std::cout << " " << std::endl;
    }

    if (nl.test_convergence())
    {
        std::cout << "average diameter h: " << average_diameter(msh) << std::endl;
        std::cout << "energy mechanic: " << nl.energy_mechanic() << std::endl;
    }

    nl.compute_discontinuous_displacement("depl2D_disc.msh");
    nl.compute_continuous_displacement("depl2D_cont.msh");
    nl.compute_discontinuous_stress("stress2D_disc.msh");
    nl.compute_continuous_stress("stress2D_cont.msh");
    nl.compute_stress_GP("stress2D_GP.msh");
    nl.compute_discontinuous_equivalent_plastic_strain("p2D_disc.msh");
    nl.compute_continuous_equivalent_plastic_strain("p2D_cont.msh");
    nl.compute_equivalent_plastic_strain_GP("p2D_GP.msh");
    nl.compute_is_plastic_GP("state2D_GP.msh");
    nl.compute_is_plastic_continuous("state2D_cont.msh");
    nl.compute_continuous_deformed("deformed2D_cont.msh");
    nl.compute_discontinuous_deformed("deformed2D_disc.msh");
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_plasticity_solver(const Mesh<T, 3, Storage>&        msh,
                      const ParamRun<T>&                rp,
                      const NLE::MaterialParameters<T>& material_data)
{
    typedef Mesh<T, 3, Storage>                            mesh_type;
    typedef static_vector<T, 3>                            result_type;
    typedef static_matrix<T, 3, 3>                         result_grad_type;
    typedef disk::mechanics::BoundaryConditions<mesh_type> Bnd_type;

    auto load = [material_data](const point<T, 3>& p) -> result_type { return result_type{0, 0, 0}; };

    auto solution = [material_data](const point<T, 3>& p) -> result_type { return result_type{0, 0, 0}; };

    auto gradient = [material_data](const point<T, 3>& p) -> result_grad_type { return result_grad_type::Zero(); };

    Bnd_type bnd(msh);
    // bnd.addDirichletEverywhere(solution);

    // // Sphere

    auto zero = [material_data](const point<T, 3>& p) -> result_type { return result_type{0.0, 0.0, 0.0}; };

    auto pres = [material_data](const point<T, 3>& p) -> result_type {
        result_type er = result_type::Zero();

        er(0) = p.x();
        er(1) = p.y();
        er(2) = p.z();

        er /= er.norm();

        return er;
    };

    bnd.addDirichletBC(disk::mechanics::DX, 3, zero);
    bnd.addDirichletBC(disk::mechanics::DY, 13, zero);
    bnd.addDirichletBC(disk::mechanics::DZ, 24, zero);
    bnd.addNeumannBC(disk::mechanics::NEUMANN, 27, pres);

    // Cube + pres
    //    auto zero = [material_data](const point<T, 3>& p) -> result_type {
    //       return result_type{0.0, 0.0, 0.0};
    //    };

    //    auto pres = [material_data](const point<T, 3>& p) -> result_type {
    //       return result_type{0.0, 0.0, -356};
    //    };

    //    bnd.addDirichletBC(disk::mechanics::CLAMPED, 69, zero);
    //    bnd.addDirichletBC(disk::mechanics::CLAMPED, 55, zero);
    //    bnd.addDirichletBC(disk::mechanics::CLAMPED, 31, zero);
    //    bnd.addDirichletBC(disk::mechanics::CLAMPED, 96, zero);
    //    bnd.addDirichletBC(disk::mechanics::DY, 62, zero);
    //    bnd.addDirichletBC(disk::mechanics::DY, 14, zero);
    //    bnd.addDirichletBC(disk::mechanics::DX, 4, zero);
    //    bnd.addDirichletBC(disk::mechanics::DX, 38, zero);
    //    bnd.addNeumannBC(disk::mechanics::NEUMANN, 21, pres);

    // Beam 3D

    // auto zero = [material_data](const point<T, 3>& p) -> result_type { return result_type{0.0, 0.0, 0.0}; };

    // auto pres = [material_data](const point<T, 3>& p) -> result_type { return result_type{0.0, 0.0, -3000.0}; };

    // bnd.addDirichletBC(disk::mechanics::CLAMPED, 3, zero);
    // bnd.addNeumannBC(disk::mechanics::NEUMANN, 13, pres);


    // Thin Plate 3D

    // auto zero = [material_data](const point<T, 3>& p) -> result_type { return result_type{0.0, 0.0, 0.0}; };

    // auto pres = [material_data](const point<T, 3>& p) -> result_type { return result_type{0.0, 0.0, -2.0E-1}; };

    // bnd.addDirichletBC(disk::mechanics::CLAMPED, 3, zero);
    // bnd.addDirichletBC(disk::mechanics::CLAMPED, 13, zero);
    // bnd.addDirichletBC(disk::mechanics::CLAMPED, 23, zero);
    // bnd.addDirichletBC(disk::mechanics::CLAMPED, 27, zero);
    // bnd.addNeumannBC(disk::mechanics::NEUMANN, 33, pres);

    // Solver

    plasticity_solver<mesh_type> nl(msh, bnd, rp, material_data);

    if (nl.verbose())
    {
        std::cout << "Solving the problem ..." << '\n';
    }

    SolverInfo solve_info = nl.compute(load);

    if (nl.verbose())
    {
        std::cout << " " << std::endl;
        std::cout << "------------------------------------------------------- " << std::endl;
        std::cout << "Summaring: " << std::endl;
        std::cout << "Total Newton's iterations: " << solve_info.m_iter << " in " << solve_info.m_time_step
                  << " load increments" << std::endl;
        std::cout << "Total time to solve the problem: " << solve_info.m_time_solver << " sec" << std::endl;
        std::cout << "**** Assembly time: " << solve_info.m_newton_info.m_assembly_info.m_time_assembly << " sec"
                  << std::endl;
        std::cout << "****** Gradient reconstruction: " << solve_info.m_newton_info.m_assembly_info.m_time_gradrec
                  << " sec" << std::endl;
        std::cout << "****** Stabilisation: " << solve_info.m_newton_info.m_assembly_info.m_time_stab << " sec"
                  << std::endl;
        std::cout << "****** Elementary computation: " << solve_info.m_newton_info.m_assembly_info.m_time_elem << " sec"
                  << std::endl;
        std::cout << "       *** Behavior computation: " << solve_info.m_newton_info.m_assembly_info.m_time_law
                  << " sec" << std::endl;
        std::cout << "****** Static condensation: " << solve_info.m_newton_info.m_assembly_info.m_time_statcond
                  << " sec" << std::endl;
        std::cout << "**** Postprocess time: " << solve_info.m_newton_info.m_assembly_info.m_time_postpro << " sec"
                  << std::endl;
        std::cout << "**** Solver time: " << solve_info.m_newton_info.m_solve_info.m_time_solve << " sec" << std::endl;
        std::cout << "------------------------------------------------------- " << std::endl;
        std::cout << " " << std::endl;
    }

    if (nl.test_convergence())
    {
        std::cout << "average diameter h: " << average_diameter(msh) << std::endl;
        // std::cout << "l2 error: " << nl.compute_l2_error(solution) << std::endl;
    }

    nl.compute_discontinuous_displacement("depl3D_disc.msh");
    nl.compute_continuous_displacement("depl3D_cont.msh");
    nl.compute_discontinuous_stress("stress3D_disc.msh");
    nl.compute_continuous_stress("stress3D_cont.msh");
    nl.compute_stress_GP("stress3D_GP.msh");
    nl.compute_discontinuous_equivalent_plastic_strain("p3D_disc.msh");
    nl.compute_continuous_equivalent_plastic_strain("p3D_cont.msh");
    nl.compute_equivalent_plastic_strain_GP("p3D_GP.msh");
    nl.compute_is_plastic_GP("state3D_GP.msh");
    nl.compute_is_plastic_continuous("state3D_cont.msh");
    nl.compute_continuous_deformed("deformed3D_cont.msh");
    nl.compute_discontinuous_deformed("deformed3D_disc.msh");

    nl.compute_sphere("resu_sphere.dat");
}

int
main(int argc, char** argv)
{
    using RealType = double;

    char* mesh_filename = nullptr;

    ParamRun<RealType> rp;

    // Elasticity Parameters
    NLE::MaterialParameters<RealType> material_data;

    //    material_data.mu       = 82E4;
    //    material_data.lambda   = 11E5;
    //    RealType E  = 200000.;
    //    RealType nu = 0.3;
    //    RealType ET = 40000;

    //    material_data.mu     = material_data.converttomu(E, nu);
    //    material_data.lambda = material_data.converttolambda(E, nu);

    //    material_data.K = 20000;
    //    material_data.H = material_data.converttoH(E, ET, material_data.K);

    //    material_data.sigma_y0 = 400;

    // // Cook Parameters (mm, GPa, kN)
    // RealType E  = 70;
    // RealType nu = 0.4999;

    // material_data.mu     = material_data.converttomu(E, nu);
    // material_data.lambda = material_data.converttolambda(E, nu);

    // material_data.K = 0.0;
    // material_data.H = 0.135;

    // material_data.sigma_y0 = 0.243;

    // Sphere Parameters (mm, GPa, kN)
    RealType E  = 210000;
    RealType nu = 0.3;
    RealType ET = 0;

    material_data.mu     = material_data.converttomu(E, nu);
    material_data.lambda = material_data.converttolambda(E, nu);

    material_data.K = 0;
    material_data.H = 0;

    material_data.sigma_y0 = 240;

    // Plaque Parameters(mm, GPa, kN)
    // RealType E  = 70;
    // RealType nu = 0.3;

    // material_data.mu     = material_data.converttomu(E, nu);
    // material_data.lambda = material_data.converttolambda(E, nu);

    // material_data.K = 0;
    // material_data.H = 10;

    // material_data.sigma_y0 = 0.8;

    // Cube Parameters (mm, MPa, N)
    //    RealType E  = 200000;
    //    RealType ET = 0;
    //    RealType nu = 0.499;

    //    material_data.mu     = material_data.converttomu(E, nu);
    //    material_data.lambda = material_data.converttolambda(E, nu);

    //    material_data.K = 0;
    //    material_data.H = material_data.converttoH(E, ET, material_data.K);

    //    material_data.sigma_y0 = 150;

    // Test aster Parameters (mm, GPa, kN)
    // RealType E  = 70;
    // RealType nu = 0.4999;

    // material_data.mu     = material_data.converttomu(E, nu);
    // material_data.lambda = material_data.converttolambda(E, nu);

    // material_data.K = 0.0;
    // material_data.H = 70;

    // material_data.sigma_y0 = 0.5;


    // Test Beam (elastic)
    // RealType E  = 30E6;
    // RealType nu = 0.3;

    // material_data.mu     = material_data.converttomu(E, nu);
    // material_data.lambda = material_data.converttolambda(E, nu);

    // material_data.K = 0.0;
    // material_data.H = 30E6;

    // material_data.sigma_y0 = 10E12;

    //Test Thin Plate (elastic)
    // RealType E  = 250;
    // RealType nu = 0.3;

    // material_data.mu     = material_data.converttomu(E, nu);
    // material_data.lambda = material_data.converttolambda(E, nu);

    // material_data.K = 0.0;
    // material_data.H = 500;

    // material_data.sigma_y0 = 10E20;

    int ch;

    while ((ch = getopt(argc, argv, "r:")) != -1)
    {
        switch (ch)
        {
            case 'r':
                if (!rp.readParameters(optarg))
                    exit(1);
                break;
            default: std::cout << "wrong arguments" << std::endl; exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    if (argc == 0)
    {
        std::cout << "Error" << std::endl;
        return 0;
    }

    mesh_filename = argv[0];

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$")))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        auto msh = disk::load_fvca5_2d_mesh<RealType>(mesh_filename);
        run_plasticity_solver(msh, rp, material_data);
        return 0;
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$")))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);
        run_plasticity_solver(msh, rp, material_data);
        return 0;
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$"))) {
       std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
       auto msh = disk::load_cartesian_2d_mesh<RealType>(mesh_filename);
       run_plasticity_solver(msh, rp, material_data);
       return 0;
    }

    /* Medit 2d*/
    if (std::regex_match(mesh_filename, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        auto msh = disk::load_medit_2d_mesh<RealType>(mesh_filename);
        run_plasticity_solver(msh, rp, material_data);
        return 0;
    }

    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$"))) {
       std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
       auto msh = disk::load_netgen_3d_mesh<RealType>(mesh_filename);
       run_plasticity_solver(msh, rp, material_data);
       return 0;
    }

    /* Medit 3d*/
    if (std::regex_match(mesh_filename, std::regex(".*\\.medit3d$"))) {
       std::cout << "Guessed mesh format: Medit format" << std::endl;
       auto msh = disk::load_medit_3d_mesh<RealType>(mesh_filename);
       run_plasticity_solver(msh, rp, material_data);
       return 0;
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$"))) {
       std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
       auto msh = disk::load_cartesian_3d_mesh<RealType>(mesh_filename);
       run_plasticity_solver(msh, rp, material_data);
       return 0;
    }

    /* FVCA6 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.msh$"))) {
       std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
       auto msh = disk::load_fvca6_3d_mesh<RealType>(mesh_filename);
       run_plasticity_solver(msh, rp, material_data);
       return 0;
    }
}
