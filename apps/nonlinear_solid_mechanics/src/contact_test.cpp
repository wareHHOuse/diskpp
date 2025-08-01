/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2020, 2024               nicolas.pignet@enpc.fr
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

#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <unistd.h>

#include "diskpp/common/colormanip.h"

#include "diskpp/boundary_conditions/boundary_conditions.hpp"
#include "diskpp/common/timecounter.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/mechanics/NewtonSolver/NewtonSolver.hpp"
#include "diskpp/mechanics/behaviors/laws/behaviorlaws.hpp"

struct error_type
{
    int    degree;
    int    nb_dof;
    double h;
    double error_L2;
    double error_H1;
};

struct run_params
{
    int  degree;
    int  l;
    bool verbose;
    bool cell_based;
};

void
usage()
{
    std::cout << "Usage: %s <options> <filename> " << std::endl;
    std::cout << "    -2: test 2D mesh (default)" << std::endl;
    std::cout << "    -3: test 3D mesh" << std::endl;
    std::cout << "    -c: cell_based method" << std::endl;
    std::cout << "    -k: face degree (>=1)" << std::endl;
    std::cout << "    -g: gamma_0 parameter for Nitsche method" << std::endl;
    std::cout << "    -t: theta parameter for Nitsche method" << std::endl;
    std::cout << "    -v: verbose" << std::endl;
    std::cout << "    -t: test robustness" << std::endl;
}

template<typename Mesh>
void
renumber_boundaries_2d(Mesh& msh)
{
    using T      = typename Mesh::coordinate_type;
    auto storage = msh.backend_storage();

    auto is_close_to = [](const T val, const T ref) -> bool
    {
        T eps = 1E-7;
        return std::abs(ref - val) < eps;
    };

    /*           -----------------
     *           !       3       !
     *           !               !
     *      4    !               ! 2
     *           !               !
     *           !               !
     *           -----------------
     *                   1
     * */

    for (size_t face_i = 0; face_i < msh.faces_size(); face_i++)
    {
        auto fc = *std::next(msh.faces_begin(), face_i);
        if (storage->boundary_info.at(face_i).is_boundary())
        {
            const auto bar = barycenter(msh, fc);
            if (is_close_to(bar.y(), T(0)))
            {
                storage->boundary_info.at(face_i).id(1);
            }
            else if (is_close_to(bar.x(), T(1)))
            {
                storage->boundary_info.at(face_i).id(2);
            }
            else if (is_close_to(bar.y(), T(1)))
            {
                storage->boundary_info.at(face_i).id(3);
            }
            else if (is_close_to(bar.x(), T(0)))
            {
                storage->boundary_info.at(face_i).id(4);
            }
            else
            {
                throw std::invalid_argument("dont find the boundaries");
            }
        }
    }
}

template<typename Mesh>
void
renumber_boundaries_3d(Mesh& msh)
{
    using T      = typename Mesh::coordinate_type;
    auto storage = msh.backend_storage();

    auto is_close_to = [](const T val, const T ref) -> bool
    {
        T eps = 1E-7;
        return std::abs(ref - val) < eps;
    };

    for (size_t face_i = 0; face_i < msh.faces_size(); face_i++)
    {
        auto fc = *std::next(msh.faces_begin(), face_i);
        if (storage->boundary_info.at(face_i).is_boundary())
        {
            const auto bar = barycenter(msh, fc);
            if (is_close_to(bar.z(), T(0)))
            {
                storage->boundary_info.at(face_i).id(1);
            }
            else if (is_close_to(bar.z(), T(1)))
            {
                storage->boundary_info.at(face_i).id(2);
            }
            else if (is_close_to(bar.x(), T(0)))
            {
                storage->boundary_info.at(face_i).id(2);
            }
            else if (is_close_to(bar.x(), T(1)))
            {
                storage->boundary_info.at(face_i).id(2);
            }
            else if (is_close_to(bar.y(), T(0)))
            {
                storage->boundary_info.at(face_i).id(2);
            }
            else if (is_close_to(bar.y(), T(1)))
            {
                storage->boundary_info.at(face_i).id(2);
            }
            else
            {
                throw std::invalid_argument("dont find the boundaries");
            }
        }
    }
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
error_type
run_tresca_solver(Mesh<T, 2, Storage>&            msh,
                  const NewtonSolverParameter<T>& rp,
                  const disk::MaterialData<T>&    material_data)
{
    typedef Mesh<T, 2, Storage>                         mesh_type;
    typedef disk::static_vector<T, 2>                   result_type;
    typedef disk::static_matrix<T, 2, 2>                grad_type;
    typedef disk::vector_boundary_conditions<mesh_type> Bnd_type;

    timecounter tc;
    tc.tic();

    renumber_boundaries_2d(msh);

    auto load = [material_data](const disk::point<T, 2>& p, const T& time) -> result_type
    {
        T y      = p.y();
        T x      = p.x();
        T exy    = std::exp(x * y);
        T mu     = material_data.getMu();
        T lambda = material_data.getLambda();
        T coeff  = 0.5 * lambda * y * y - x * x * (0.5 * lambda + 1);

        T fx = -mu * (lambda * y + x * coeff) + y * (lambda + mu * (lambda + 2)) * (x * y + 2);
        T fy = -lambda * x * (mu - 1) * (x * y + 2) + mu * (2 * x * (0.5 * lambda + 1) - y * coeff);

        return -time * result_type{fx, fy} * exy / (3.0 * (1.0 + material_data.getLambda()));
    };

    auto solution = [material_data](const disk::point<T, 2>& p, const T time) -> result_type
    {
        T y     = p.y();
        T x     = p.x();
        T exy   = std::exp(x * y);
        T coeff = 1.0 / (1.0 + material_data.getLambda());

        return time*result_type{x * exy * (1.0 + coeff), y * exy * (-1.0 + coeff)} / 6.0;
    };

    auto sol = [solution](const disk::point<T, 2>& p) -> auto
        { return solution(p, 1.0); };

    auto s = [rp, material_data](const disk::point<T, 2>& p) -> T
    {
        T y      = p.y();
        T x      = p.x();
        T mu     = material_data.getMu();
        T lambda = material_data.getLambda();

        return mu * x * x * (0.5 * lambda + 1.0) / (3 * lambda + 3.0);
    };

    Bnd_type bnd(msh);
    bnd.addContactBC(disk::SIGNORINI_FACE, 1, s);
    // bnd.addDirichletBC(disk::DIRICHLET, 1, solution);
    bnd.addDirichletBC(disk::DIRICHLET, 2, solution);
    bnd.addDirichletBC(disk::DIRICHLET, 3, solution);
    bnd.addDirichletBC(disk::DIRICHLET, 4, solution);

    disk::mechanics::NewtonSolver<mesh_type> nl(msh, bnd, rp);

    nl.addBehavior(disk::DeformationMeasure::SMALL_DEF, disk::LawType::ELASTIC);
    nl.addMaterialData(material_data);

    nl.initial_guess(sol);

    if (nl.verbose())
    {
        std::cout << "Solving the problem ..." << '\n';
    }

    SolverInfo solve_info = nl.compute(load);

    if (nl.verbose())
    {
        solve_info.printInfo();
    }

    error_type error;
    error.h        = average_diameter(msh);
    error.degree   = rp.m_face_degree;
    error.nb_dof   = nl.numberOfDofs();
    error.error_L2 = nl.compute_l2_displacement_error(sol);
    error.error_H1 = nl.compute_H1_error(sol);

    // nl.compute_stress_GP("stress2D_GP_test.msh");
    // nl.compute_continuous_displacement("depl2D_cont_test.msh");

    return error;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
error_type
run_tresca_solver(const Mesh<T, 3, Storage>&      msh,
                  const NewtonSolverParameter<T>& rp,
                  const disk::MaterialData<T>&    material_data)
{
    typedef Mesh<T, 3, Storage>                         mesh_type;
    typedef disk::static_vector<T, 3>                   result_type;
    typedef disk::static_matrix<T, 3, 3>                grad_type;
    typedef disk::vector_boundary_conditions<mesh_type> Bnd_type;

    timecounter tc;
    tc.tic();

    renumber_boundaries_3d(msh);

    auto load = [material_data](const disk::point<T, 3>& p, const T& time) -> result_type
    {
        T z      = p.z();
        T x      = p.x();
        T exz    = std::exp(x * z);
        T mu     = material_data.getMu();
        T lambda = material_data.getLambda();
        T coeff  = 0.5 * lambda * z * z - x * x * (0.5 * lambda + 1);

        T fx = -mu * (lambda * z + x * coeff) + z * (lambda + mu * (lambda + 2)) * (x * z + 2);
        T fz = -lambda * x * (mu - 1) * (x * z + 2) + mu * (2 * x * (0.5 * lambda + 1) - z * coeff);

        return -time * result_type{fx, 0.0, fz} * exz / (3.0 * (1.0 + material_data.getLambda()));
    };

    auto solution = [material_data](const disk::point<T, 3>& p, const T& time) -> result_type
    {
        T z     = p.z();
        T x     = p.x();
        T exz   = std::exp(x * z);
        T coeff = 1.0 / (1.0 + material_data.getLambda());

        return time*result_type{x * exz * (1.0 + coeff), 0.0, z * exz * (-1.0 + coeff)} / 6.0;
    };

    auto sol = [solution](const disk::point<T, 3>& p) -> auto
        { return solution(p, 1.0); };

    auto s = [rp, material_data](const disk::point<T, 3>& p) -> T
    {
        T x      = p.x();
        T mu     = material_data.getMu();
        T lambda = material_data.getLambda();

        return mu * x * x * (0.5 * lambda + 1.0) / (3 * lambda + 3.0);
    };

    Bnd_type bnd(msh);
    bnd.addContactBC(disk::SIGNORINI_FACE, 1, s);
    bnd.addDirichletBC(disk::DIRICHLET, 2, solution);

    disk::mechanics::NewtonSolver<mesh_type> nl(msh, bnd, rp);

    nl.addBehavior(disk::DeformationMeasure::SMALL_DEF, disk::LawType::ELASTIC);
    nl.addMaterialData(material_data);

    nl.initial_guess(sol);

    SolverInfo solve_info = nl.compute(load);

    if (nl.verbose())
    {
        solve_info.printInfo();
    }

    error_type error;
    error.h        = average_diameter(msh);
    error.degree   = rp.m_face_degree;
    error.nb_dof   = nl.numberOfDofs();
    error.error_L2 = nl.compute_l2_displacement_error(sol);
    error.error_H1 = nl.compute_H1_error(sol);

    return error;
}

void
printResults(const std::vector<error_type>& error)
{
    if (error.size() > 0)
    {
        std::ios::fmtflags f(std::cout.flags());
        std::cout.precision(4);
        std::cout.setf(std::iostream::scientific, std::iostream::floatfield);

        std::cout << "Convergence test for k = " << error[0].degree << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << "| Size mesh  | Displacement | Convergence | Displacem  | Convergence |    Total   |" << std::endl;
        std::cout << "|    h       |   L2 error   |     rate    |  H1 error  |     rate    | faces DOF  |" << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;

        std::string s_dof = " " + std::to_string(error[0].nb_dof) + "                  ";
        s_dof.resize(10);

        std::cout << "| " << error[0].h << " |  " << error[0].error_L2 << "  | "
                  << "     -     "
                  << " | " << error[0].error_H1 << " | "
                  << "     -     "
                  << " | " << s_dof << " |" << std::endl;

        for (int i = 1; i < error.size(); i++)
        {
            s_dof = " " + std::to_string(error[i].nb_dof) + "                  ";
            s_dof.resize(10);
            double rate_depl =
              (log10(error[i - 1].error_L2) - log10(error[i].error_L2)) / (log10(error[i - 1].h) - log10(error[i].h));
            double rate_stress =
              (log10(error[i - 1].error_H1) - log10(error[i].error_H1)) / (log10(error[i - 1].h) - log10(error[i].h));

            std::cout << "| " << error[i].h << " |  " << error[i].error_L2 << "  |  " << rate_depl << " | "
                      << error[i].error_H1 << " |  " << rate_stress << " | " << s_dof << " |" << std::endl;
        }

        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << "  " << std::endl;
        std::cout.flags(f);
    }
    else
        std::cout << "The file error is empty" << std::endl;
}

void
printResults2(const std::vector<error_type>& error)
{
    if (error.size() > 0)
    {
        std::ios::fmtflags f(std::cout.flags());
        std::cout.precision(4);
        std::cout.setf(std::iostream::scientific, std::iostream::floatfield);

        std::cout << "Robustesse test for k = " << error[0].degree << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << "| Lambda coef| Displacement | Difference  | Displacem  | Difference  |    Total   |" << std::endl;
        std::cout << "|    nu       |   L2 error   |             |  H1 error  |             | faces DOF  |"
                  << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;

        std::string s_dof = " " + std::to_string(error[0].nb_dof) + "                  ";
        s_dof.resize(10);

        std::cout << "| " << error[0].h << " |  " << error[0].error_L2 << "  | "
                  << "     -     "
                  << " | " << error[0].error_H1 << " | "
                  << "     -     "
                  << " | " << s_dof << " |" << std::endl;

        for (int i = 1; i < error.size(); i++)
        {
            s_dof = " " + std::to_string(error[i].nb_dof) + "                  ";
            s_dof.resize(10);
            double rate_depl   = error[i].error_L2 - error[i - 1].error_L2;
            double rate_stress = error[i].error_H1 - error[i - 1].error_H1;

            std::cout << "| " << error[i].h << " |  " << error[i].error_L2 << "  |  " << rate_depl << " | "
                      << error[i].error_H1 << " |  " << rate_stress << " | " << s_dof << " |" << std::endl;
        }

        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << "  " << std::endl;
        std::cout.flags(f);
    }
    else
        std::cout << "The file error is empty" << std::endl;
}

template<typename T>
void
test_triangles_fvca5(const NewtonSolverParameter<T>& rp, const disk::MaterialData<T>& material_data)
{
    int runs = 4;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_1.typ1");
    paths.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_2.typ1");
    paths.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_3.typ1");
    paths.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_4.typ1");
    paths.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_5.typ1");

    std::vector<error_type> error_sumup;

    for (int i = 0; i < runs; i++)
    {
        disk::generic_mesh<T, 2> msh;
        disk::load_mesh_fvca5_2d<T>(paths[i].c_str(), msh);
        error_sumup.push_back(run_tresca_solver(msh, rp, material_data));
    }
    printResults(error_sumup);
}

template<typename T>
void
test_triangles_netgen(const NewtonSolverParameter<T>& rp, const disk::MaterialData<T>& material_data)
{
    int runs = 4;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri01.mesh2d");
    paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri02.mesh2d");
    paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri03.mesh2d");
    paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri04.mesh2d");
    paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri05.mesh2d");

    std::vector<error_type> error_sumup;

    for (int i = 0; i < runs; i++)
    {
        disk::simplicial_mesh<T, 2> msh;
        disk::load_mesh_netgen(paths[i].c_str(), msh);
        error_sumup.push_back(run_tresca_solver(msh, rp, material_data));
    }
    printResults(error_sumup);
}

template<typename T>
void
test_hexagons(const NewtonSolverParameter<T>& rp, const disk::MaterialData<T>& material_data)
{
    int runs = 5;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_1.typ1");
    paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_2.typ1");
    paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_3.typ1");
    paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_4.typ1");
    paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_5.typ1");

    std::vector<error_type> error_sumup;

    for (int i = 0; i < runs; i++)
    {
        disk::generic_mesh<T, 2> msh;
        disk::load_mesh_fvca5_2d<T>(paths[i].c_str(), msh);
        error_sumup.push_back(run_tresca_solver(msh, rp, material_data));
    }
    printResults(error_sumup);
}

template<typename T>
void
test_kershaws(const NewtonSolverParameter<T>& rp, const disk::MaterialData<T>& material_data)
{
    int runs = 5;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_1.typ1");
    paths.push_back("../../../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_2.typ1");
    paths.push_back("../../../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_3.typ1");
    paths.push_back("../../../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_4.typ1");
    paths.push_back("../../../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_5.typ1");

    std::vector<error_type> error_sumup;

    for (int i = 0; i < runs; i++)
    {
        disk::generic_mesh<T, 2> msh;
        disk::load_mesh_fvca5_2d<T>(paths[i].c_str(), msh);
        error_sumup.push_back(run_tresca_solver(msh, rp, material_data));
    }
    printResults(error_sumup);
}

template<typename T>
void
test_quads_fvca5(const NewtonSolverParameter<T>& rp, const disk::MaterialData<T>& material_data)
{
    int runs = 5;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_1.typ1");
    paths.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_2.typ1");
    paths.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_3.typ1");
    paths.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_4.typ1");
    paths.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_5.typ1");

    std::vector<error_type> error_sumup;

    for (int i = 0; i < runs; i++)
    {
        disk::generic_mesh<T, 2> msh;
        disk::load_mesh_fvca5_2d<T>(paths[i].c_str(), msh);
        error_sumup.push_back(run_tresca_solver(msh, rp, material_data));
    }
    printResults(error_sumup);
}

template<typename T>
void
test_quads_diskpp(const NewtonSolverParameter<T>& rp, const disk::MaterialData<T>& material_data)
{
    int runs = 5;

    std::vector<std::string> paths;
    // paths.push_back("../../../diskpp/meshes/2D_quads/diskpp/testmesh-1-1.quad");
    paths.push_back("../../../diskpp/meshes/2D_quads/diskpp/testmesh-2-2.quad");
    paths.push_back("../../../diskpp/meshes/2D_quads/diskpp/testmesh-4-4.quad");
    paths.push_back("../../../diskpp/meshes/2D_quads/diskpp/testmesh-8-8.quad");
    paths.push_back("../../../diskpp/meshes/2D_quads/diskpp/testmesh-16-16.quad");
    paths.push_back("../../../diskpp/meshes/2D_quads/diskpp/testmesh-32-32.quad");
    paths.push_back("../../../diskpp/meshes/2D_quads/diskpp/testmesh-256-256.quad");

    std::vector<error_type> error_sumup;

    for (int i = 0; i < runs; i++)
    {
        disk::cartesian_mesh<T, 2> msh;
        disk::load_mesh_diskpp_cartesian(paths[i].c_str(), msh);
        error_sumup.push_back(run_tresca_solver(msh, rp, material_data));
    }
    printResults(error_sumup);
}

template<typename T>
void
test_hexahedra_diskpp(const NewtonSolverParameter<T>& rp, const disk::MaterialData<T>& material_data)
{
    int runs = 4;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-2-2-2.hex");
    paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-4-4-4.hex");
    paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-8-8-8.hex");
    paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-16-16-16.hex");
    paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-32-32-32.hex");

    std::vector<error_type> error_sumup;

    for (int i = 0; i < runs; i++)
    {
        disk::cartesian_mesh<T, 3> msh;
        disk::load_mesh_diskpp_cartesian(paths[i].c_str(), msh);
        error_sumup.push_back(run_tresca_solver(msh, rp, material_data));
    }
    printResults(error_sumup);
}

template<typename T>
void
test_hexahedra_fvca6(const NewtonSolverParameter<T>& rp, const disk::MaterialData<T>& material_data)
{
    int runs = 4;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_2x2x2.msh");
    paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_4x4x4.msh");
    paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_8x8x8.msh");
    paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_16x16x16.msh");
    paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_32x32x32.hex");

    std::vector<error_type> error_sumup;

    for (int i = 0; i < runs; i++)
    {
        disk::generic_mesh<T, 3> msh;
        disk::load_mesh_fvca6_3d<T>(paths[i].c_str(), msh);
        error_sumup.push_back(run_tresca_solver(msh, rp, material_data));
    }
    printResults(error_sumup);
}

template<typename T>
void
test_tetrahedra_netgen(const NewtonSolverParameter<T>& rp, const disk::MaterialData<T>& material_data)
{
    int runs = 4;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/3D_tetras/netgen/fvca6_tet0.mesh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/netgen/fvca6_tet1.mesh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/netgen/fvca6_tet2.mesh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/netgen/fvca6_tet3.mesh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/netgen/fvca6_tet4.mesh");

    std::vector<error_type> error_sumup;

    for (int i = 0; i < runs; i++)
    {
        disk::simplicial_mesh<T, 3> msh;
        disk::load_mesh_netgen(paths[i].c_str(), msh);
        error_sumup.push_back(run_tresca_solver(msh, rp, material_data));
    }
    printResults(error_sumup);
}

template<typename T>
void
test_polyhedra_fvca6(const NewtonSolverParameter<T>& rp, const disk::MaterialData<T>& material_data)
{
    int runs = 3;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_10.msh");
    paths.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_20.msh");
    paths.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_30.msh");
    paths.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_40.msh");

    std::vector<error_type> error_sumup;

    for (int i = 0; i < runs; i++)
    {
        disk::generic_mesh<T, 3> msh;
        disk::load_mesh_fvca6_3d<T>(paths[i].c_str(), msh);
        error_sumup.push_back(run_tresca_solver(msh, rp, material_data));
    }
    printResults(error_sumup);
}

template<typename T>
void
test_tetrahedra_fvca6(const NewtonSolverParameter<T>& rp, const disk::MaterialData<T>& material_data)
{
    int runs = 4;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/3D_tetras/fvca6/tet.0.msh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/fvca6/tet.1.msh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/fvca6/tet.2.msh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/fvca6/tet.3.msh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/fvca6/tet.4.msh");

    std::vector<error_type> error_sumup;

    for (int i = 0; i < runs; i++)
    {
        disk::generic_mesh<T, 3> msh;
        disk::load_mesh_fvca6_3d<T>(paths[i].c_str(), msh);
        error_sumup.push_back(run_tresca_solver(msh, rp, material_data));
    }
    printResults(error_sumup);
}

template<typename T>
void
test_triangles_fvca5_robust(const NewtonSolverParameter<T>& rp, disk::MaterialData<T>& material_data)
{

    std::vector<error_type> error_sumup;

    for (int i = 0; i < 5; i++)
    {
        material_data.setLambda(std::pow(10.0, i));
        std::string              path = "../../../diskpp/meshes/2D_triangles/fvca5/mesh1_3.typ1";
        disk::generic_mesh<T, 2> msh;
        disk::load_mesh_fvca5_2d<T>(path.c_str(), msh);
        error_sumup.push_back(run_tresca_solver(msh, rp, material_data));
        error_sumup[i].h = material_data.getLambda();
    }
    printResults2(error_sumup);
}

template<typename T>
void
test_quads_diskpp_robust(const NewtonSolverParameter<T>& rp, disk::MaterialData<T>& material_data)
{

    std::vector<error_type> error_sumup;

    for (int i = 0; i < 5; i++)
    {
        material_data.setLambda(std::pow(10.0, i));
        std::string                path = "../../../diskpp/meshes/2D_quads/diskpp/testmesh-16-16.quad";
        disk::cartesian_mesh<T, 2> msh;
        disk::load_mesh_diskpp_cartesian(path.c_str(), msh);
        error_sumup.push_back(run_tresca_solver(msh, rp, material_data));
        error_sumup[i].h = material_data.getLambda();
    }
    printResults2(error_sumup);
}

template<typename T>
void
test_hexagons_robust(const NewtonSolverParameter<T>& rp, disk::MaterialData<T>& material_data)
{

    std::vector<error_type> error_sumup;

    for (int i = 0; i < 5; i++)
    {
        material_data.setLambda(std::pow(10.0, i));
        std::string              path = "../../../diskpp/meshes/2D_hex/fvca5/hexagonal_4.typ1";
        disk::generic_mesh<T, 2> msh;
        disk::load_mesh_fvca5_2d<T>(path.c_str(), msh);
        error_sumup.push_back(run_tresca_solver(msh, rp, material_data));
        error_sumup[i].h = material_data.getLambda();
    }
    printResults2(error_sumup);
}

template<typename T>
void
test_kershaws_robust(const NewtonSolverParameter<T>& rp, disk::MaterialData<T>& material_data)
{

    std::vector<error_type> error_sumup;

    for (int i = 0; i < 5; i++)
    {
        material_data.setLambda(std::pow(10.0, i));
        std::string              path = "../../../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_3.typ1";
        disk::generic_mesh<T, 2> msh;
        disk::load_mesh_fvca5_2d<T>(path.c_str(), msh);
        error_sumup.push_back(run_tresca_solver(msh, rp, material_data));
        error_sumup[i].h = material_data.getLambda();
    }
    printResults2(error_sumup);
}

template<typename T>
void
test_tetrahedra_fvca6_robust(const NewtonSolverParameter<T>& rp, disk::MaterialData<T>& material_data)
{

    std::vector<error_type> error_sumup;

    for (int i = 0; i < 5; i++)
    {
        material_data.setLambda(std::pow(10.0, i));
        std::string              path = "../../../diskpp/meshes/3D_tetras/fvca6/tet.2.msh";
        disk::generic_mesh<T, 3> msh;
        disk::load_mesh_fvca6_3d<T>(path.c_str(), msh);
        error_sumup.push_back(run_tresca_solver(msh, rp, material_data));
        error_sumup[i].h = material_data.getLambda();
    }
    printResults2(error_sumup);
}

template<typename T>
void
test_hexahedra_fvca6_robust(const NewtonSolverParameter<T>& rp, disk::MaterialData<T>& material_data)
{

    std::vector<error_type> error_sumup;

    for (int i = 0; i < 5; i++)
    {
        material_data.setLambda(std::pow(10.0, i));
        std::string              path = "../../../diskpp/meshes/3D_hexa/fvca6/hexa_8x8x8.msh";
        disk::generic_mesh<T, 3> msh;
        disk::load_mesh_fvca6_3d<T>(path.c_str(), msh);
        error_sumup.push_back(run_tresca_solver(msh, rp, material_data));
        error_sumup[i].h = material_data.getLambda();
    }
    printResults2(error_sumup);
}

template<typename T>
void
test_polyhedra_fvca6_robust(const NewtonSolverParameter<T>& rp, disk::MaterialData<T>& material_data)
{

    std::vector<error_type> error_sumup;

    for (int i = 0; i < 5; i++)
    {
        material_data.setLambda(std::pow(10.0, i));
        std::string              path = "../../../diskpp/meshes/3D_general/fvca6/dbls_20.msh";
        disk::generic_mesh<T, 3> msh;
        disk::load_mesh_fvca6_3d<T>(path.c_str(), msh);
        error_sumup.push_back(run_tresca_solver(msh, rp, material_data));
        error_sumup[i].h = material_data.getLambda();
    }
    printResults2(error_sumup);
}

int
main(int argc, char** argv)
{
    using RealType = double;

    int degree = 1;
    int dim    = 2;

    // Elasticity Parameters
    disk::MaterialData<RealType> material_data;

    material_data.setMu(2.0);
    material_data.setLambda(2.E4);

    std::cout << "Material coefficients:  E = " << material_data.getE() << " , nu = " << material_data.getNu()
              << std::endl;

    // Solver parameters
    NewtonSolverParameter<RealType> rp;
    rp.m_precomputation    = true;
    rp.m_stab_type         = HHO;
    rp.m_epsilon           = 2.0E-6;
    rp.m_time_step.front() = std::make_pair(1.0, 1);
    rp.m_sublevel          = 2;
    rp.m_beta              = 2 * material_data.getMu();
    rp.m_gamma_0           = 2 * material_data.getMu();
    rp.m_frot_type         = TRESCA;
    rp.m_threshold         = 0;
    rp.m_theta             = 0;
    rp.m_iter_max          = 20;
    rp.m_sublevel          = 0;

    int  ch;
    bool robust = false;

    while ((ch = getopt(argc, argv, "23g:k:t:vr")) != -1)
    {
        switch (ch)
        {
            case '2': dim = 2; break;
            case '3':
                dim = 3;
                break;

                // case 'g': rp.m_gamma_0 = atof(optarg); break;

            case 'k':
                degree = atoi(optarg);
                if (degree < 0)
                {
                    std::cout << "Degree must be positive. Falling back to 1." << std::endl;
                    degree = 1;
                }
                rp.m_face_degree = degree;
                rp.m_cell_degree = degree;
                rp.m_grad_degree = degree;
                break;

            case 't': rp.m_theta = atof(optarg); break;

            case 'v': rp.m_verbose = true; break;

            case 'r': robust = true; break;

            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                usage();
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    timecounter tc;

    if (rp.m_theta != -1.0)
    {
        rp.m_gamma_0 *= (rp.m_face_degree + 1.0) * (rp.m_face_degree + dim);
    }

    std::cout << " Test convergence rates for: " << std::endl;
    std::cout << " ** Face_Degree = " << rp.m_face_degree << std::endl;
    std::cout << " " << std::endl;

    if (dim == 3)
    {
        if (robust)
        {
            tc.tic();
            std::cout << "-Tetrahedras fvca6:" << std::endl;
            test_tetrahedra_fvca6_robust<RealType>(rp, material_data);
            tc.toc();
            std::cout << "Time to test robustess: " << tc.elapsed() << std::endl;
            std::cout << " " << std::endl;

            tc.tic();
            std::cout << "-Hexahedras fvca6:" << std::endl;
            test_hexahedra_fvca6_robust<RealType>(rp, material_data);
            tc.toc();
            std::cout << "Time to test robustess: " << tc.elapsed() << std::endl;
            std::cout << " " << std::endl;

            tc.tic();
            std::cout << "-Polyhedra:" << std::endl;
            test_polyhedra_fvca6_robust<RealType>(rp, material_data);
            tc.toc();
            std::cout << "Time to test robustess: " << tc.elapsed() << std::endl;
            std::cout << " " << std::endl;
        }
        else
        {
            tc.tic();
            std::cout << "-Tetrahedras fvca6:" << std::endl;
            test_tetrahedra_fvca6<RealType>(rp, material_data);
            tc.toc();
            std::cout << "Time to test convergence rates: " << tc.elapsed() << std::endl;
            std::cout << " " << std::endl;

            tc.tic();
            std::cout << "-Tetrahedras netgen:" << std::endl;
            test_tetrahedra_netgen<RealType>(rp, material_data);
            tc.toc();
            std::cout << "Time to test convergence rates: " << tc.elapsed() << std::endl;
            std::cout << " " << std::endl;

            tc.tic();
            std::cout << "-Hexahedras fvca6:" << std::endl;
            test_hexahedra_fvca6<RealType>(rp, material_data);
            tc.toc();
            std::cout << "Time to test convergence rates: " << tc.elapsed() << std::endl;
            std::cout << " " << std::endl;

            tc.tic();
            std::cout << "-Hexahedras diskpp:" << std::endl;
            test_hexahedra_diskpp<RealType>(rp, material_data);
            tc.toc();
            std::cout << "Time to test convergence rates: " << tc.elapsed() << std::endl;
            std::cout << " " << std::endl;

            tc.tic();
            std::cout << "-Polyhedra:" << std::endl;
            test_polyhedra_fvca6<RealType>(rp, material_data);
            tc.toc();
            std::cout << "Time to test convergence rates: " << tc.elapsed() << std::endl;
            std::cout << " " << std::endl;
        }
    }
    else if (dim == 2)
    {
        if (robust)
        {
            tc.tic();
            std::cout << "-Triangles fvca5:" << std::endl;
            test_triangles_fvca5_robust<RealType>(rp, material_data);
            tc.toc();
            std::cout << "Time to test robustess: " << tc.elapsed() << std::endl;
            std::cout << " " << std::endl;

            tc.tic();
            std::cout << "-Quadrangles diskpp:" << std::endl;
            test_quads_diskpp_robust<RealType>(rp, material_data);
            tc.toc();
            std::cout << "Time to test robustess: " << tc.elapsed() << std::endl;
            std::cout << " " << std::endl;

            tc.tic();
            std::cout << "-Hexagons:" << std::endl;
            test_hexagons_robust<RealType>(rp, material_data);
            tc.toc();
            std::cout << "Time to test robustess: " << tc.elapsed() << std::endl;
            std::cout << " " << std::endl;

            tc.tic();
            std::cout << "-Kershaws:" << std::endl;
            test_kershaws_robust<RealType>(rp, material_data);
            tc.toc();
            std::cout << "Time to test robustess: " << tc.elapsed() << std::endl;
            std::cout << " " << std::endl;
        }
        else
        {
            tc.tic();
            std::cout << "-Triangles fvca5:" << std::endl;
            test_triangles_fvca5<RealType>(rp, material_data);
            tc.toc();
            std::cout << "Time to test convergence rates: " << tc.elapsed() << std::endl;
            std::cout << " " << std::endl;

            tc.tic();
            std::cout << "-Triangles netgen:" << std::endl;
            test_triangles_netgen<RealType>(rp, material_data);
            tc.toc();
            std::cout << "Time to test convergence rates: " << tc.elapsed() << std::endl;
            std::cout << " " << std::endl;

            tc.tic();
            std::cout << "-Quadrangles fvca5:" << std::endl;
            test_quads_fvca5<RealType>(rp, material_data);
            tc.toc();
            std::cout << "Time to test convergence rates: " << tc.elapsed() << std::endl;
            std::cout << " " << std::endl;

            tc.tic();
            std::cout << "-Quadrangles diskpp:" << std::endl;
            test_quads_diskpp<RealType>(rp, material_data);
            tc.toc();
            std::cout << "Time to test convergence rates: " << tc.elapsed() << std::endl;
            std::cout << " " << std::endl;

            tc.tic();
            std::cout << "-Hexagons:" << std::endl;
            test_hexagons<RealType>(rp, material_data);
            tc.toc();
            std::cout << "Time to test convergence rates: " << tc.elapsed() << std::endl;
            std::cout << " " << std::endl;

            tc.tic();
            std::cout << "-Kershaws:" << std::endl;
            test_kershaws<RealType>(rp, material_data);
            tc.toc();
            std::cout << "Time to test convergence rates: " << tc.elapsed() << std::endl;
            std::cout << " " << std::endl;
        }
    }
}