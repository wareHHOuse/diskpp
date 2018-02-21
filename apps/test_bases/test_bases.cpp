/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
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

#include "core/bases/bases_final.hpp"
#include "core/loaders/loader.hpp"
#include "core/quadratures/quadratures.hpp"

template<typename Mesh>
void
test_bases(const Mesh& msh)
{
    
    typedef Mesh mesh_type;
    typedef typename mesh_type::cell  cell_type;
    typedef typename mesh_type::face  face_type;
    typedef typename mesh_type::scalar_type  scalar_type;    

    typedef revolution::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    typedef revolution::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

    typedef disk::quadrature<mesh_type, cell_type>    cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>    face_quadrature_type;

    typedef dynamic_matrix<scalar_type> matrix_type;

    size_t degree = 1;

    cell_quadrature_type    cell_quad(2 * degree);
    face_quadrature_type    face_quad(2 * degree);
    
    for(auto cl : msh)
    {
        cell_basis_type          cell_basis(msh, cl, degree);

        auto qps = cell_quad.integrate(msh, cl);

        matrix_type mass = matrix_type::Zero(cell_basis.size(), cell_basis.size());
        
        for (auto& qp : qps)
        {
            auto phi = cell_basis.eval_functions(msh, cl, qp.point());
            mass += qp.weight() * phi * phi.transpose();
        }

        std::cout << "Cell mass matrix for " << cl << std::endl;
        std::cout << mass << std::endl;

        matrix_type stiff = matrix_type::Zero(cell_basis.size(), cell_basis.size());

        for (auto& qp : qps)
        {
            auto dphi = cell_basis.eval_gradients(msh, cl, qp.point());
            stiff += qp.weight() * dphi * dphi.transpose();
        }

        auto fcs = faces(msh, cl);
        
        for(auto fc : fcs)        
        {
            face_basis_type          face_basis(msh, fc, degree);
            matrix_type mass = matrix_type::Zero(face_basis.size(), face_basis.size());
    
            auto qps = face_quad.integrate(msh, fc);
            for (auto& qp : qps)
            {
                auto phi = face_basis.eval_basis(msh, fc, qp.point());
                mass += qp.weight() * phi * phi.transpose();
            }

            std::cout << "Face mass matrix for " << fc << std::endl;
            std::cout << mass << std::endl;
        
        }
            
    }
   dump_to_matlab(msh, "test.m");
   
   
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
