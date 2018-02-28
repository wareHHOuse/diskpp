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
#include <iomanip>
#include <regex>

#include <unistd.h>

#include "revolution/bases"
#include "revolution/quadratures"
#include "revolution/methods/hho"

#include "core/loaders/loader.hpp"


template<typename Mesh>
typename Mesh::scalar_type
run_test(const Mesh& msh, size_t degree)
{
    typedef Mesh mesh_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;


    auto f = [&](const point_type& pt) -> scalar_type {
        return std::sin(2. * M_PI * pt.x()) * std::sin(2. * M_PI * pt.y());
    };

    typename revolution::hho_degree_info hdi(degree);

    scalar_type error = 0.0;
    for (auto& cl : msh)
    {
        Matrix<scalar_type, Dynamic, 1> proj = revolution::project_function(msh, cl, hdi, f);
        auto gr = revolution::make_hho_scalar_laplacian(msh, cl, hdi);

        size_t rec_size = revolution::scalar_basis_size(hdi.reconstruction_degree(), Mesh::dimension);

        Matrix<scalar_type, Dynamic, 1> reconstr = Matrix<scalar_type, Dynamic, 1>::Zero(rec_size);
        reconstr.tail(rec_size-1) = gr.first * proj;
        reconstr(0) = proj(0);

        auto cb = revolution::make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
        Matrix<scalar_type, Dynamic, Dynamic> mass = revolution::make_mass_matrix(msh, cl, cb);
        Matrix<scalar_type, Dynamic, 1> rhs = revolution::make_rhs(msh, cl, cb, f);
        Matrix<scalar_type, Dynamic, 1> exp_reconstr = mass.llt().solve(rhs);

        Matrix<scalar_type, Dynamic, 1> diff = reconstr - exp_reconstr;

        error += diff.dot(mass*diff);
    }

    return std::sqrt(error);
}


int main(void)
{
    using T = double;

    std::vector<std::string> meshfiles;
    
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_1.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_2.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_3.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_4.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_5.typ1");
    
    /*
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_1.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_2.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_3.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_4.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_5.typ1");
    */

    for (size_t k = 0; k < 6; k++)
    {
        std::cout << "DEGREE " << k << std::endl;

        std::vector<T> mesh_hs;
        std::vector<T> l2_errors;

        for (size_t i = 0; i < meshfiles.size(); i++)
        {
            typedef disk::generic_mesh<T, 2>  mesh_type;

            mesh_type msh;
            disk::fvca5_mesh_loader<T, 2> loader;
            if (!loader.read_mesh(meshfiles.at(i)))
            {
                std::cout << "Problem loading mesh." << std::endl;
                continue;
            }
            loader.populate_mesh(msh);

            auto error = run_test(msh, k);
            mesh_hs.push_back( disk::mesh_h(msh) );
            l2_errors.push_back(error);
        }

        for (size_t i = 0; i < mesh_hs.size(); i++)
        {
            if (i == 0)
            {
                std::cout << "    ";
                std::cout << std::scientific << std::setprecision(5) << mesh_hs.at(i) << "    ";
                std::cout << std::scientific << std::setprecision(5) << l2_errors.at(i);
                std::cout << "     -- " << std::endl;
            }
            else
            {
                auto rate = std::log( l2_errors.at(i)/l2_errors.at(i-1) ) /
                            std::log( mesh_hs.at(i)/mesh_hs.at(i-1) );
                std::cout << "    ";
                std::cout << std::scientific << std::setprecision(5) << mesh_hs.at(i) << "    ";
                std::cout << std::scientific << std::setprecision(5) << l2_errors.at(i) << "    ";
                std::cout << std::defaultfloat << std::setprecision(3) << rate << std::endl;
            }
        }
    }

    return 0;
}
