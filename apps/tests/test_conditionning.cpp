/*
 *       /\        Matteo Cicuttin (C) 2016-2020
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2020                     nicolas.pignet@enpc.fr
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

#include <unistd.h>

#include "bases/bases.hpp"

#include "loaders/loader.hpp"
#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"
#include "output/postMesh.hpp"

#include "common.hpp"

  template<typename T>
  T
  compute_conditionning(const Matrix<T, Dynamic, Dynamic>& mat)
{
    using Mat             = Matrix<T, Dynamic, Dynamic>;
    const auto ev = SelfAdjointEigenSolver<Mat>(mat).eigenvalues();
    return ev(ev.size() - 1) / ev(0);
}

template<typename Mesh, typename Elem, typename Basis>
typename Mesh::coordinate_type
compute_mass_conditionning(const Mesh& msh, const Elem& elem, const Basis& basis)
{
    const auto mass = make_mass_matrix(msh, elem, basis);

    // std::cout << mass << std::endl;

    return compute_conditionning(mass);
}

template<typename Mesh>
void
conditionning(const Mesh& msh, const std::string& filename)
{

    std::ofstream ofs(filename);

    ofs << "\\begin{tabular} { | c | c | c | c | c | }" << std::endl;
    ofs << "\\hline" << std::endl;

    ofs << "degree &  CSC  &  CSC-rot  &   CLG  &   CLG-rot \\\\" << std::endl;
    ofs << "\\hline" << std::endl;

    for (size_t degree = 0; degree <= 6; degree++)
    {
        for(auto& cl : msh)
        {
            auto cb_SC     = make_scalar_monomial_basis(msh, cl, degree, false);
            auto cb_SC_rot = make_scalar_monomial_basis(msh, cl, degree, true);
            auto cb_LG     = make_scalar_legendre_basis(msh, cl, degree, false);
            auto cb_LG_rot = make_scalar_legendre_basis(msh, cl, degree, true);

            ofs << degree << " & " << compute_mass_conditionning(msh, cl, cb_SC) << " & "
                << compute_mass_conditionning(msh, cl, cb_SC_rot) << " & " << compute_mass_conditionning(msh, cl, cb_LG)
                << " & " << compute_mass_conditionning(msh, cl, cb_LG_rot) << " \\\\" << std::endl;

            ofs << "\\hline" << std::endl;

            // auto fcs = faces(msh, cl);
            // for(auto& fc : fcs)
            // {
            //     auto fb_SC = make_scalar_monomial_basis(msh, fc, degree, true);
            //     auto fb_SC_rot = make_scalar_monomial_basis(msh, fc, degree, true);

            //     std::cout << degree << "F " << compute_mass_conditionning(msh, fc, fb_SC) << " "
            //               << compute_mass_conditionning(msh, fc, fb_SC_rot) << std::endl;
            // }
        }
    }

    ofs << "\\end{tabular}";

    ofs.close();

    disk::PostMesh<Mesh> post_mesh = disk::PostMesh<Mesh>(msh);
    gmsh::Gmesh          gmsh      = disk::convertMesh(post_mesh);

    gmsh.writeGmesh("mesh.msh", 2);
}


int
main(int argc, char** argv)
{
    using T = double;

    char* mesh_filename = nullptr;
    char* filename      = nullptr;

    int ch;
    while ((ch = getopt(argc, argv, "f:m:")) != -1)
    {
        switch (ch)
        {
            case 'f': filename = optarg; break;

            case 'm': mesh_filename = optarg; break;

            case '?':
            default: std::cout << "Invalid option" << std::endl; return 1;
        }
    }

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$")))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        auto msh = disk::load_fvca5_2d_mesh<T>(mesh_filename);
        conditionning(msh, filename);
        return 0;
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$")))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<T>(mesh_filename);
        conditionning(msh, filename);
        return 0;
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$")))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        auto msh = disk::load_cartesian_2d_mesh<T>(mesh_filename);
        conditionning(msh, filename);
        return 0;
    }

    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$")))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        auto msh = disk::load_netgen_3d_mesh<T>(mesh_filename);
        conditionning(msh, filename);
        return 0;
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$")))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        auto msh = disk::load_cartesian_3d_mesh<T>(mesh_filename);
        conditionning(msh, filename);
        return 0;
    }

    /* FVCA6 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.msh$")))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
        disk::generic_mesh<T, 3> msh;

        disk::load_mesh_fvca6_3d<T>(mesh_filename, msh);
        conditionning(msh, filename);

        return 0;
    }
}