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

#include "core/quadratures/quad_from_proton.hpp"
#include "loaders/loader.hpp"


template<typename MeshType>
void
process_mesh(const MeshType& msh)
{
    size_t degree = 2;

    using mesh_type = MeshType;
    using point_type = typename mesh_type::point_type;
    using coordinate_type = typename mesh_type::coordinate_type;

    auto f1 = [](const point_type& pt) -> coordinate_type {
        return pt.x();
    };

    auto f2 = [](const point_type& pt) -> coordinate_type {
        return pt.y();
    };

    auto f3 = [](const point_type& pt) -> coordinate_type {
        return pt.x() * pt.y();
    };

    coordinate_type intval_f1 = 0.0;
    coordinate_type intval_f2 = 0.0;
    coordinate_type intval_f3 = 0.0;

    for (auto& cl : msh)
    {
        auto qps = revolution::integrate(msh, cl, degree);

        for (auto& qp : qps)
        {
            intval_f1 += qp.weight() * f1(qp.point());
            intval_f2 += qp.weight() * f2(qp.point());
            intval_f3 += qp.weight() * f3(qp.point());
        }

        /*
        auto fcs = faces(msh, cl);
        for (auto& fc : fcs)
        {
            revolution::integrate(msh, fc, degree);
        }
        */
    }

    std::cout << "f(x,y) = x   : " << intval_f1 << std::endl;
    std::cout << "f(x,y) = y   : " << intval_f2 << std::endl;
    std::cout << "f(x,y) = x*y : " << intval_f3 << std::endl;
}


int main(int argc, char **argv)
{
    using RealType = double;

    char    *filename       = nullptr;
    int     elems_1d        = 8;
    int ch;

    if (argc != 2)
    {
        std::cout << "Filename plz" << std::endl;
        return 0;
    }

/*
    if (argc == 1)
    {
        std::cout << "Mesh format: 1D uniform" << std::endl;

        typedef disk::generic_mesh<RealType, 1>  mesh_type;

        mesh_type msh;
        disk::uniform_mesh_loader<RealType, 1> loader(0,1,elems_1d);
        loader.populate_mesh(msh);

        process_mesh(msh);

        return 0;
    }
*/
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
/*
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
    */
    return 0;
}
