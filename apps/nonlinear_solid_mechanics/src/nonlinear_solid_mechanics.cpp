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
 * Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework.
 * M. Abbas, A. Ern, N. Pignet.
 * International Journal of Numerical Methods in Engineering (2019)
 * 120(3), 303-327
 * DOI: 10.1002/nme.6137
 */

#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <unistd.h>

#include "diskpp/boundary_conditions/boundary_conditions.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/mechanics/NewtonSolver/NewtonSolver.hpp"
#include "diskpp/mechanics/behaviors/laws/behaviorlaws.hpp"

#include "diskpp/common/timecounter.hpp"

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_nl_solid_mechanics_solver(const Mesh<T, 2, Storage>&      msh,
                              const NewtonSolverParameter<T>& rp,
                              const disk::MaterialData<T>&    material_data)
{
    typedef Mesh<T, 2, Storage>                         mesh_type;
    typedef disk::static_vector<T, 2>                   result_type;
    typedef disk::static_matrix<T, 2, 2>                result_grad_type;
    typedef disk::vector_boundary_conditions<mesh_type> Bnd_type;

    auto load = [material_data](const disk::point<T, 2>& p) -> result_type { return result_type{0, 0}; };

    auto solution = [material_data](const disk::point<T, 2>& p) -> result_type { return result_type{0, 0}; };

    Bnd_type bnd(msh);

    // Cook with quadrilaterals
    auto zero = [material_data](const disk::point<T, 2>& p) -> result_type { return result_type{0.0, 0}; };

    auto trac = [material_data](const disk::point<T, 2>& p) -> result_type { return result_type{0.0, 0.3125}; };

    bnd.addDirichletBC(disk::CLAMPED, 3, zero);
    bnd.addNeumannBC(disk::NEUMANN, 8, trac);

    disk::mechanics::NewtonSolver<mesh_type> nl(msh, bnd, rp);

#ifdef HAVE_MGIS
    // To use a law developped with Mfront
    const auto hypo = mgis::behaviour::Hypothesis::PLANESTRAIN;
    const std::string filename = "src/libBehaviour.dylib";
    // nl.addBehavior(filename, "IsotropicLinearHardeningPlasticity", hypo);
    nl.addBehavior(filename, "LogarithmicStrainPlasticity", hypo);
#else
    // To use a native law from DiSk++
    nl.addBehavior(disk::DeformationMeasure::LOGARITHMIC_DEF, disk::LawType::LINEAR_HARDENING);
#endif

    nl.addMaterialData(material_data);

    nl.initial_guess(zero);

    if (nl.verbose())
    {
        std::cout << "Solving the problem ..." << '\n';
    }

    SolverInfo solve_info = nl.compute(load);

    if (nl.verbose())
    {
        solve_info.printInfo();
    }

    if (nl.convergence())
    {
        std::cout << "average diameter h: " << average_diameter(msh) << std::endl;
    }
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_nl_solid_mechanics_solver(const Mesh<T, 3, Storage>&      msh,
                              const NewtonSolverParameter<T>& rp,
                              const disk::MaterialData<T>&    material_data)
{
    typedef Mesh<T, 3, Storage>                         mesh_type;
    typedef disk::static_vector<T, 3>                   result_type;
    typedef disk::static_matrix<T, 3, 3>                result_grad_type;
    typedef disk::vector_boundary_conditions<mesh_type> Bnd_type;

    Bnd_type bnd(msh);

    auto load = [material_data](const disk::point<T, 3>& p) -> result_type { return result_type{0, 0, 0}; };

    auto solution = [material_data](const disk::point<T, 3>& p) -> result_type { return result_type{0, 0, 0}; };

    auto zero = [material_data](const disk::point<T, 3>& p) -> result_type { return result_type{0.0, 0.0, 0.0}; };

    auto pres = [material_data](const disk::point<T, 3>& p) -> result_type {
        result_type er = result_type::Zero();

        er(0) = p.x();
        er(1) = p.y();
        er(2) = p.z();

        er /= er.norm();

        return 3 * er;
    };

    auto deplr = [material_data](const disk::point<T, 3>& p) -> result_type {
        result_type er = result_type::Zero();

        er(0) = p.x();
        er(1) = p.y();
        er(2) = p.z();

        er /= er.norm();

        return 0.157 * er;
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

    disk::mechanics::NewtonSolver<mesh_type> nl(msh, bnd, rp);

#ifdef HAVE_MGIS
    // To use a law developped with Mfront
    const auto        hypo     = mgis::behaviour::Hypothesis::TRIDIMENSIONAL;
    const std::string filename = "src/libBehaviour.dylib";
    nl.addBehavior(filename, "LogarithmicStrainPlasticity", hypo);
#else
    // To use a native law from DiSk++
    nl.addBehavior(disk::DeformationMeasure::SMALL_DEF, disk::LawType::LINEAR_HARDENING);
#endif

    nl.addMaterialData(material_data);

    nl.initial_guess(zero);

    if (nl.verbose())
    {
        std::cout << "Solving the problem ..." << '\n';
    }

    SolverInfo solve_info = nl.compute(load);

    if (nl.verbose())
    {
        solve_info.printInfo();
    }

    if (nl.convergence())
    {
        std::cout << "average diameter h: " << average_diameter(msh) << std::endl;
    }
}

int
main(int argc, char** argv)
{
    using RealType = double;

    char* mesh_filename = nullptr;

    NewtonSolverParameter<RealType> rp;

    const RealType MPa = 10E6;
    const RealType GPa = 10E9;

    // Elasticity Parameters
    disk::MaterialData<RealType> material_data;

    // // Cook Parameters HPP (mm, MPa, kN)
    RealType E  = 70;
    RealType nu = 0.4999;

    material_data.setMu(E, nu);
    material_data.setLambda(E, nu);

    material_data.setK(0.0);
    material_data.setH(0.135);

    material_data.setSigma_y0(0.243);

    material_data.addMfrontParameter("YoungModulus", material_data.getE());
    material_data.addMfrontParameter("PoissonRatio", material_data.getNu());
    material_data.addMfrontParameter("HardeningSlope", material_data.getH());
    material_data.addMfrontParameter("YieldStrength", material_data.getSigma_y0());

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
    // RealType E  = 28.95;
    // RealType nu = 0.3;
    // RealType ET = 0;

    // material_data.setMu(E, nu);
    // material_data.setLambda(E, nu);
    // material_data.setK(0);
    // material_data.setH(0.0);
    // material_data.setSigma_y0(6);

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

    /* Poly 2d*/
    if (std::regex_match(mesh_filename, std::regex(".*\\.poly2d$")))
    {
        std::cout << "Guessed mesh format: Poly2D format" << std::endl;
        disk::generic_mesh<RealType, 2> msh;
        disk::load_mesh_poly2d<RealType>(mesh_filename, msh);
        run_nl_solid_mechanics_solver(msh, rp, material_data);
        return 0;
    }

    /* Poly 3d*/
    if (std::regex_match(mesh_filename, std::regex(".*\\.poly3d$")))
    {
        std::cout << "Guessed mesh format: Poly3D format" << std::endl;
        disk::generic_mesh<RealType, 3> msh;
        disk::load_mesh_poly3d<RealType>(mesh_filename, msh);
        run_nl_solid_mechanics_solver(msh, rp, material_data);
        return 0;
    }

    /* Medit 2d*/
    if (std::regex_match(mesh_filename, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        auto msh = disk::load_medit_2d_mesh<RealType>(mesh_filename);
        run_nl_solid_mechanics_solver(msh, rp, material_data);
        return 0;
    }

    /* Medit 3d*/
    if (std::regex_match(mesh_filename, std::regex(".*\\.medit3d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        auto msh = disk::load_medit_3d_mesh<RealType>(mesh_filename);
        run_nl_solid_mechanics_solver(msh, rp, material_data);
        return 0;
    }
}
