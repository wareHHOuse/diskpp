/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
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
#include "boundary_conditions/boundary_conditions.hpp"
#include "loaders/loader.hpp"
#include "mechanics/behaviors/laws/materialData.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "tresca_solver.hpp"

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_tresca_solver(const Mesh<T, 2, Storage>& msh, const ParamRun<T>& rp, const disk::MaterialData<T>& material_data)
{
    typedef Mesh<T, 2, Storage>                         mesh_type;
    typedef static_vector<T, 2>                         result_type;
    typedef static_matrix<T, 2, 2>                      result_grad_type;
    typedef disk::vector_boundary_conditions<mesh_type> Bnd_type;

    auto load = [material_data](const disk::point<T, 2>& p) -> result_type { return result_type{0, 0}; };

    auto solution = [material_data](const disk::point<T, 2>& p) -> result_type { return result_type{0, 0}; };

    Bnd_type bnd(msh);
    // bnd.addDirichletEverywhere(solution);

    auto s = [rp](const disk::point<T, 2>& p) -> T { return rp.m_threshold; };

    // Bostan2
    auto zero = [material_data](const disk::point<T, 2>& p) -> result_type { return result_type{0.0, 0}; };
    auto neum = [material_data](const disk::point<T, 2>& p) -> result_type { return result_type{400.0, 0}; };

    bnd.addContactBC(disk::SIGNORINI_FACE, 6, s);
    bnd.addDirichletBC(disk::DIRICHLET, 10, zero);
    // bnd.addDirichletBC(disk::DX, 6, zero);
    bnd.addNeumannBC(disk::NEUMANN, 3, neum);

    // // Hertz
    // auto depl = [material_data](const disk::point<T, 2>& p) -> result_type { return result_type{0.0, -5.0}; };
    // auto zero = [material_data](const disk::point<T, 2>& p) -> result_type { return result_type{0.0, 0.0}; };

    // bnd.addContactBC(disk::SIGNORINI_FACE, 26, s);
    // bnd.addContactBC(disk::SIGNORINI_FACE, 16, s);

    // bnd.addDirichletBC(disk::DX, 7, zero);
    // bnd.addDirichletBC(disk::DX, 14, zero);
    // bnd.addDirichletBC(disk::DY, 4, depl);
    // bnd.addDirichletBC(disk::DY, 19, depl);

    // solver

    tresca_solver<mesh_type> nl(msh, bnd, rp, material_data);

    if (nl.verbose())
    {
        std::cout << "Solving the problem ..." << '\n';
    }

    SolverInfo solve_info = nl.compute(load);

    if (nl.verbose())
    {
        solve_info.printInfo();
    }

    // nl.compute_discontinuous_displacement("depl2D_disc.msh");
    nl.compute_continuous_displacement("depl2D_cont.msh");
    // nl.compute_discontinuous_stress("stress2D_disc.msh");
    nl.compute_continuous_stress("stress2D_cont.msh");
    nl.compute_stress_GP("stress2D_GP.msh");
    nl.compute_equivalent_plastic_strain_GP("plastic2D_p_GP_cont.msh");
    nl.compute_is_plastic_GP("isplastic2D_GP_cont.msh");
    nl.compute_continuous_deformed("deformed2D_cont.msh");
    // nl.compute_discontinuous_deformed("deformed2D_disc.msh");
    nl.contact_quantities("contact.dat");
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_tresca_solver(const Mesh<T, 3, Storage>& msh, const ParamRun<T>& rp, const disk::MaterialData<T>& material_data)
{
    typedef Mesh<T, 3, Storage>                         mesh_type;
    typedef static_vector<T, 3>                         result_type;
    typedef static_matrix<T, 3, 3>                      result_grad_type;
    typedef disk::vector_boundary_conditions<mesh_type> Bnd_type;

    auto load = [material_data](const disk::point<T, 3>& p) -> result_type { return result_type{0, 0, 0}; };

    auto solution = [material_data](const disk::point<T, 3>& p) -> result_type { return result_type{0, 0, 0}; };

    Bnd_type bnd(msh);
    // bnd.addDirichletEverywhere(solution);

    auto s = [rp](const disk::point<T, 3>& p) -> T { return rp.m_threshold; };

    // // Bostan2
    // auto zero = [material_data](const disk::point<T, 3>& p) -> result_type { return result_type{0.0, 0.0, 0.0}; };
    // auto neum = [material_data](const disk::point<T, 3>& p) -> result_type { return result_type{400.0, 0.0, 0.0}; };

    // bnd.addContactBC(disk::SIGNORINI_FACE, 23);
    // bnd.addDirichletBC(disk::CLAMPED, 27, zero);
    // bnd.addDirichletBC(disk::DZ, 31, zero);
    // bnd.addDirichletBC(disk::DZ, 33, zero);
    // bnd.addNeumannBC(disk::NEUMANN, 3, neum);



    // // Hertz
    auto depl= [material_data](const disk::point<T, 3>& p) -> result_type { return result_type{0.0, 0.0, -2.0}; };
    auto zero = [material_data](const disk::point<T, 3>& p) -> result_type { return result_type{0.0, 0.0, 0.0}; };

    bnd.addContactBC(disk::SIGNORINI_FACE, 31, s);
    bnd.addDirichletBC(disk::DX, 19, zero);
    bnd.addDirichletBC(disk::DX, 37, zero);
    bnd.addDirichletBC(disk::DY, 27, zero);
    bnd.addDirichletBC(disk::DY, 40, zero);
    bnd.addDirichletBC(disk::DZ, 14, depl);

    // Solver

    tresca_solver<mesh_type> nl(msh, bnd, rp, material_data);

    if (nl.verbose())
    {
        std::cout << "Solving the problem ..." << '\n';
    }

    SolverInfo solve_info = nl.compute(load);

    if (nl.verbose())
    {
        solve_info.printInfo();
    }

    nl.compute_discontinuous_displacement("depl3D_disc.msh");
    nl.compute_continuous_displacement("depl3D_cont.msh");
    // nl.compute_discontinuous_stress("stress3D_disc.msh");
    // nl.compute_continuous_stress("stress3D_cont.msh");
    nl.compute_stress_GP("stress3D_GP.msh");
    nl.compute_equivalent_plastic_strain_GP("plastic_p_GP_cont.msh");
    nl.compute_is_plastic_GP("isplastic_GP_cont.msh");
    // nl.compute_continuous_deformed("deformed3D_cont.msh");
    // nl.compute_discontinuous_deformed("deformed3D_disc.msh");
    nl.contact_quantities("contact.dat");
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

    // Hertz elasticity
    RealType E  = 70;
    RealType nu = 0.3;

    material_data.setMu(E, nu);
    material_data.setLambda(E, nu);

    // Hertz plasticity
        // RealType E  = 70;
        // RealType nu = 0.3;

        // material_data.setMu(E, nu);
        // material_data.setLambda(E, nu);

        // material_data.setK(0.0);
        // material_data.setH(7.0);
        // material_data.setSigma_y0(2.5);

        // material_data.setMu(1.0);
        // material_data.setLambda(1000.0);

        // Bostan
        // material_data.setMu(1000, 0.3);
        // material_data.setLambda(1000, 0.3);

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

    /* Medit 2d*/
    if (std::regex_match(mesh_filename, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        auto msh = disk::load_medit_2d_mesh<RealType>(mesh_filename);
        run_tresca_solver(msh, rp, material_data);
        return 0;
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$")))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);
        run_tresca_solver(msh, rp, material_data);
        return 0;
    }

    /* Medit 3d*/
    if (std::regex_match(mesh_filename, std::regex(".*\\.medit3d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        auto msh = disk::load_medit_3d_mesh<RealType>(mesh_filename);
        run_tresca_solver(msh, rp, material_data);
        return 0;
    }
}
