/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
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

#include <iostream>
#include <regex>

#include "loaders/loader.hpp"

#include "core/bases/bases.hpp"
#include "quadratures/quadratures.hpp"


template<typename MeshType>
void
process_mesh(const MeshType& msh)
{
    typedef MeshType                            mesh_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;
    typedef typename mesh_type::scalar_type     scalar_type;
    typedef dynamic_matrix<scalar_type>         matrix_type;

    const size_t degree = 1;

    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, degree);

        std::cout << "Diameter: " << diameter(msh, cl) << std::endl;
        std::cout << "Measure:  " << measure(msh, cl) << std::endl;

        matrix_type stiff_mat = matrix_type::Zero(cb.size(), cb.size());

        const auto qps = disk::integrate(msh, cl, 2 * degree);
        for (auto& qp : qps)
        {
            const auto dphi = cb.eval_gradients(qp.point());
            stiff_mat += qp.weight() * dphi * dphi.transpose();
        }

        std::cout << "Stiffness matrix for " << cl << std::endl;
        std::cout << stiff_mat << std::endl;
    }
}


int main(int argc, char **argv)
{
    using RealType = double;

    char    *filename       = nullptr;
    int     elems_1d        = 8;
    int ch;


    if (argc == 1)
    {
        // std::cout << "Mesh format: 1D uniform" << std::endl;

        // typedef disk::generic_mesh<RealType, 1>  mesh_type;

        // mesh_type msh;
        // disk::uniform_mesh_loader<RealType, 1> loader(0,1,elems_1d);
        // loader.populate_mesh(msh);

        // process_mesh(msh);

        return 0;
    }

    filename = argv[1];

    if (std::regex_match(filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;

        typedef disk::generic_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::fvca5_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        process_mesh(msh);

        dump_to_matlab(msh, "test.m");
    }

    if (std::regex_match(filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;

        typedef disk::simplicial_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::netgen_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        process_mesh(msh);
    }

    if (std::regex_match(filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;

        typedef disk::generic_mesh<RealType, 3>   mesh_type;

        mesh_type msh;
        disk::fvca6_mesh_loader<RealType, 3> loader;


        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }

        loader.populate_mesh(msh);

        process_mesh(msh);
    }

    if (std::regex_match(filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;

        typedef disk::simplicial_mesh<RealType, 3>   mesh_type;

        mesh_type msh;
        disk::netgen_mesh_loader<RealType, 3> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        process_mesh(msh);
    }

    if (std::regex_match(filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: Cartesian 2D" << std::endl;

        typedef disk::cartesian_mesh<RealType, 2>   mesh_type;

        mesh_type msh;
        disk::cartesian_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        process_mesh(msh);
    }

    if (std::regex_match(filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: Cartesian 3D" << std::endl;

        typedef disk::cartesian_mesh<RealType, 3>   mesh_type;

        mesh_type msh;
        disk::cartesian_mesh_loader<RealType, 3> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        process_mesh(msh);
    }

    return 0;
}
