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
 * Karol Cascavita (C) 2018                     klcascavitam@unal.edu.co
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

#include "bases/bases.hpp"
#include "quadratures/quadratures.hpp"
#include "methods/hho"

#include "core/loaders/loader.hpp"


template<typename Mesh>
void
test_bases(const Mesh& msh)
{

    typedef Mesh mesh_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;
    typedef typename mesh_type::coordinate_type scalar_type;

    typedef dynamic_matrix<scalar_type>     matrix_type;

    using point_type = typename mesh_type::point_type;

    size_t degree = 1;

    auto f1 = [](const point_type& pt) -> scalar_type {
        return std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
    };

    auto f2 = [](const point_type& pt) -> scalar_type {
        return M_PI * (std::cos(M_PI*pt.x()) + std::cos(M_PI*pt.y()));
    };


    auto v1 = [](const point_type& pt) -> Matrix<scalar_type, 2, 1> {
        Matrix<scalar_type, 2, 1> ret;
        ret(0) = pt.x();
        ret(1) = 0;
        return ret;
    };

    auto v2 = [](const point_type& pt) -> Matrix<scalar_type, 2, 1> {
        Matrix<scalar_type, 2, 1> ret;
        ret(0) = pt.y();
        ret(1) = 0;
        return ret;
    };

    auto v3 = [](const point_type& pt) -> Matrix<scalar_type, 2, 1> {
        Matrix<scalar_type, 2, 1> ret;
        ret(0) = std::sin(M_PI*pt.x());
        ret(1) = std::sin(M_PI*pt.y());
        return ret;
    };

    typename disk::hho_degree_info hdi(degree);

    for(auto cl : msh)
    {
        //revolution::scaled_monomial_vector_basis<mesh_type, cell_type> vec_cell_basis(msh, cl, degree);

        //auto cell_basis = revolution::make_vector_monomial_basis(msh, cl, degree);

        //matrix_type mass = revolution::make_mass_matrix(msh, cl, cell_basis);

        //std::cout << "Cell mass matrix for " << cl << std::endl;
        //std::cout << mass << std::endl;


        //matrix_type stiff = revolution::make_stiffness_matrix(msh, cl, cell_basis);

        //std::cout << "Cell stiffness matrix for " << cl << std::endl;
        //std::cout << stiff << std::endl;


        /*
        auto fcs = faces(msh, cl);

        for(auto fc : fcs)
        {
            revolution::scaled_monomial_vector_basis<mesh_type, face_type> vec_face_basis(msh, fc, degree);

            matrix_type mass = revolution::make_mass_matrix(msh, fc, cell_basis);

            std::cout << "Face mass matrix for " << fc << std::endl;
            std::cout << mass << std::endl;

        }

        project_function(msh, cl, hdi, f1);

        make_hho_scalar_laplacian(msh, cl, hdi);
        */

    }

    scalar_type intval = 0.0;
    scalar_type stab_error = 0.0;
    for (auto cl : msh)
    {
        //auto cb = revolution::make_scalar_monomial_basis(msh, cl, degree);
        //Matrix<scalar_type, Dynamic, Dynamic> mass = revolution::make_mass_matrix(msh, cl, cb);
        //Matrix<scalar_type, Dynamic, 1> rhs1 = revolution::make_rhs(msh, cl, cb, v1);
        //Matrix<scalar_type, Dynamic, 1> rhs2 = revolution::make_rhs(msh, cl, cb, v2);

        //auto massInv = mass.llt();

        //Matrix<scalar_type, Dynamic, 1> dofs1 = massInv.solve(rhs1);
        //Matrix<scalar_type, Dynamic, 1> dofs2 = massInv.solve(rhs2);

        //intval += dofs1.dot(mass*dofs2);

        Matrix<scalar_type, Dynamic, 1> p = project_function(msh, cl, hdi, v3);

        auto gr = make_hho_vector_laplacian(msh, cl, hdi);

        Matrix<scalar_type, Dynamic, Dynamic> stab;
        stab = make_hho_vector_stabilization(msh, cl, gr.first, hdi);


        auto dr = make_hho_divergence_reconstruction(msh, cl, hdi);
        Matrix<scalar_type, Dynamic, 1> d = dr.first * p;

        stab_error += p.dot(stab*p);
    }

    std::cout << intval << std::endl;
    std::cout << std::sqrt(stab_error) << std::endl;
}

int main(int argc, char **argv)
{
    using RealType = double;

    char    *filename       = nullptr;
    int ch;


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

        test_bases(msh);
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

        test_bases(msh);
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

        test_bases(msh);
    }


    return 0;
}
