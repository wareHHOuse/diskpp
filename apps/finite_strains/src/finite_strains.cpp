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

#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <unistd.h>

#include <map>

#include "Informations.hpp"
#include "Parameters.hpp"
#include "diskpp/mechanics/behaviors/laws/materialData.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/boundary_conditions/boundary_conditions.hpp"

#include "diskpp/common/timecounter.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

#include "finite_strains_solver.hpp"

template<typename T>
void
readCurve(const std::string filename, disk::MaterialData<T>& material_data)
{
    std::ifstream file(filename, std::ios::in);

    if (file)
    {
        // L'ouverture s'est bien passée, on peut donc lire

        T p, Rp;

        while (!file.eof()) // Tant qu'on n'est pas à la fin, on lit
        {
            file >> p >> Rp;
            if (!file.eof())
            {
                material_data.addCurvePoint(p, Rp);
            }
        }
        material_data.checkRpCurve();
        file.close();
    }
    else
    {
        std::cout << "ERREUR: Imposible to read the file" << std::endl;
    }
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_finite_strains_solver(const Mesh<T, 2, Storage>&   msh,
                          const ParamRun<T>&           rp,
                          const disk::MaterialData<T>& material_data)
{
    typedef Mesh<T, 2, Storage>                            mesh_type;
    typedef disk::static_vector<T, 2>                      result_type;
    typedef disk::static_matrix<T, 2, 2>                   result_grad_type;
    typedef disk::vector_boundary_conditions<mesh_type>    Bnd_type;

    auto load = [material_data](const disk::point<T, 2>& p) -> result_type { return result_type{0, 0}; };

    auto solution = [material_data](const disk::point<T, 2>& p) -> result_type { return result_type{0, 0}; };

    auto gradient = [material_data](const disk::point<T, 2>& p) -> result_grad_type { return result_grad_type::Zero(); };

    Bnd_type bnd(msh);


    // Cook with quadrilaterals
    auto zero = [material_data](const disk::point<T, 2>& p) -> result_type { return result_type{0.0, 0}; };

    auto trac = [material_data](const disk::point<T, 2>& p) -> result_type { return result_type{0.0, 0.3125}; };

    bnd.addDirichletBC(disk::CLAMPED, 3, zero);
    bnd.addNeumannBC(disk::NEUMANN, 8, trac);

    finite_strains_solver<mesh_type> nl(msh, bnd, rp, material_data);

    if (nl.verbose())
    {
        std::cout << "Solving the problem ..." << '\n';
    }

    SolverInfo solve_info = nl.compute(load);

    if (nl.verbose())
    {
        solve_info.printInfo();
    }

    if (nl.test_convergence())
    {
        std::cout << "average diameter h: " << average_diameter(msh) << std::endl;
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
run_finite_strains_solver(const Mesh<T, 3, Storage>&   msh,
                          const ParamRun<T>&           rp,
                          const disk::MaterialData<T>& material_data)
{
    typedef Mesh<T, 3, Storage>                            mesh_type;
    typedef disk::static_vector<T, 3>                      result_type;
    typedef disk::static_matrix<T, 3, 3>                   result_grad_type;
    typedef disk::vector_boundary_conditions<mesh_type>    Bnd_type;

    auto load = [material_data](const disk::point<T, 3>& p) -> result_type { return result_type{0, 0, 0}; };

    auto solution = [material_data](const disk::point<T, 3>& p) -> result_type { return result_type{0, 0, 0}; };

    auto gradient = [material_data](const disk::point<T, 3>& p) -> result_grad_type { return result_grad_type::Zero(); };

    Bnd_type bnd(msh);
    // bnd.addDirichletEverywhere(solution);

    // // Sphere

    auto zero = [material_data](const disk::point<T, 3>& p) -> result_type { return result_type{0.0, 0.0, 0.0}; };

    auto pres = [material_data](const disk::point<T, 3>& p) -> result_type {
        result_type er = result_type::Zero();

        er(0) = p.x();
        er(1) = p.y();
        er(2) = p.z();

        er /= er.norm();

        return 3*er;
    };

    auto deplr = [material_data](const disk::point<T, 3>& p) -> result_type {
        result_type er = result_type::Zero();

        er(0) = p.x();
        er(1) = p.y();
        er(2) = p.z();

        er /= er.norm();

        return 0.157*er;
    };

    // Sphere Hpp
    // bnd.addDirichletBC(disk::DX, 3, zero);
    // bnd.addDirichletBC(disk::DY, 13, zero);
    // bnd.addDirichletBC(disk::DZ, 24, zero);
    // bnd.addNeumannBC(disk::NEUMANN, 27, pres);

    // // Sphere GDEF
    bnd.addDirichletBC(disk::DX, 12, zero);
    bnd.addDirichletBC(disk::DY, 24, zero);
    bnd.addDirichletBC(disk::DZ, 19, zero);
    bnd.addDirichletBC(disk::DIRICHLET, 27, deplr);
    //bnd.addNeumannBC(disk::NEUMANN, 27, pres);

    // Cylindre GDEF
    // auto depl = [material_data](const disk::point<T, 3>& p) -> result_type { return result_type{0.0, 0.0, 1.0}; };
    // bnd.addDirichletBC(disk::DZ, 125, zero);
    // bnd.addDirichletBC(disk::DZ, 50, zero);
    // bnd.addDirichletBC(disk::DZ, 96, zero);
    // bnd.addDirichletBC(disk::DX, 113, zero);
    // bnd.addDirichletBC(disk::DX, 91, zero);
    // bnd.addDirichletBC(disk::DX, 74, zero);
    // bnd.addDirichletBC(disk::DX, 120, zero);
    // bnd.addDirichletBC(disk::DY, 31, zero);
    // bnd.addDirichletBC(disk::DY, 108, zero);
    // bnd.addDirichletBC(disk::DY, 128, zero);
    // bnd.addDirichletBC(disk::DY, 55, zero);
    // bnd.addDirichletBC(disk::DZ, 103, depl);
    // bnd.addDirichletBC(disk::DZ, 14, depl);
    // bnd.addDirichletBC(disk::DZ, 69, depl);

    // Cube + pres
    //    auto zero = [material_data](const disk::point<T, 3>& p) -> result_type {
    //       return result_type{0.0, 0.0, 0.0};
    //    };

    //    auto pres = [material_data](const disk::point<T, 3>& p) -> result_type {
    //       return result_type{0.0, 0.0, -356};
    //    };

    //    bnd.addDirichletBC(disk::CLAMPED, 69, zero);
    //    bnd.addDirichletBC(disk::CLAMPED, 55, zero);
    //    bnd.addDirichletBC(disk::CLAMPED, 31, zero);
    //    bnd.addDirichletBC(disk::CLAMPED, 96, zero);
    //    bnd.addDirichletBC(disk::DY, 62, zero);
    //    bnd.addDirichletBC(disk::DY, 14, zero);
    //    bnd.addDirichletBC(disk::DX, 4, zero);
    //    bnd.addDirichletBC(disk::DX, 38, zero);
    //    bnd.addNeumannBC(disk::NEUMANN, 21, pres);


    // auto trac = [material_data](const disk::point<T, 3>& p) -> result_type {
    //     return result_type { 0.0, 0.3125,  0.0}; };

    // bnd.addDirichletBC(disk::CLAMPED, 3, zero);
    // bnd.addNeumannBC(disk::NEUMANN, 20, trac);

    // Solver

    finite_strains_solver<mesh_type> nl(msh, bnd, rp, material_data);

    if (nl.verbose())
    {
        std::cout << "Solving the problem ..." << '\n';
    }

    SolverInfo solve_info = nl.compute(load);

    if (nl.verbose())
    {
        solve_info.printInfo();
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

    const RealType MPa = 10E6;
    const RealType GPa = 10E9;

    // Elasticity Parameters
    disk::MaterialData<RealType> material_data;

    // // Cook Parameters (mm, MPa, kN)
    // RealType E  = 206.9;
    // RealType nu = 0.29;

    // material_data.setMu(E, nu);
    // material_data.setLambda(E, nu);

    // readCurve("VEM2_2d.dat", material_data);

    // material_data.setK(0.0);
    // material_data.setH(E, 0.13, 0.0);

    // material_data.setSigma_y0(0.450);

    // Old cook
    // RealType E  = 70;
    // RealType nu = 0.4999;

    // material_data.setMu(E, nu);
    // material_data.setLambda(E, nu);
    // material_data.setK(0.0);
    // material_data.setH(0.135);
    // material_data.setSigma_y0(0.243);

    // Sphere Parameters (mm, MPa, kN)
    RealType E = 28.95;
    RealType nu = 0.3;
    RealType ET = 0;

    material_data.setMu(E, nu);
    material_data.setLambda(E, nu);
    material_data.setK(0);
    material_data.setH(0.0);
    material_data.setSigma_y0(6);

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

    // /* FVCA5 2D */
    // if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$")))
    // {
    //     std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
    //     auto msh = disk::load_fvca5_2d_mesh<RealType>(mesh_filename);
    //     run_finite_strains_solver(msh, rp, material_data);
    //     return 0;
    // }

    // /* Netgen 2D */
    // if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$")))
    // {
    //     std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
    //     auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);
    //     run_finite_strains_solver(msh, rp, material_data);
    //     return 0;
    // }

    // /* DiSk++ cartesian 2D */
    // if (std::regex_match(mesh_filename, std::regex(".*\\.quad$"))) {
    //    std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
    //    auto msh = disk::load_cartesian_2d_mesh<RealType>(mesh_filename);
    //    run_finite_strains_solver(msh, rp, material_data);
    //    return 0;
    // }

    /* Medit 2d*/
    if (std::regex_match(mesh_filename, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        auto msh = disk::load_medit_2d_mesh<RealType>(mesh_filename);
        run_finite_strains_solver(msh, rp, material_data);
        return 0;
    }

    // /* Netgen 3D */
    // if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$")))
    // {
    //     std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
    //     auto msh = disk::load_netgen_3d_mesh<RealType>(mesh_filename);
    //     run_finite_strains_solver(msh, rp, material_data);
    //     return 0;
    // }
    //
    /* Medit 3d*/
    if (std::regex_match(mesh_filename, std::regex(".*\\.medit3d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        auto msh = disk::load_medit_3d_mesh<RealType>(mesh_filename);
        run_finite_strains_solver(msh, rp, material_data);
        return 0;
    }
    //
    // /* DiSk++ cartesian 3D */
    // if (std::regex_match(mesh_filename, std::regex(".*\\.hex$"))) {
    //    std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
    //    auto msh = disk::load_cartesian_3d_mesh<RealType>(mesh_filename);
    //    run_finite_strains_solver(msh, rp, material_data);
    //    return 0;
    //}

    // /* FVCA6 3D */
    // if (std::regex_match(mesh_filename, std::regex(".*\\.msh$"))) {
    //    std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
    //    auto msh = disk::load_fvca6_3d_mesh<RealType>(mesh_filename);
    //    run_finite_strains_solver(msh, rp, material_data);
    //    return 0;
    // }
}
