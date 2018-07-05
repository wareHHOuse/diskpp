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

#pragma once

#include <iomanip>

#include "core/loaders/loader.hpp"

const size_t MIN_TEST_DEGREE = 0;
const size_t MAX_TEST_DEGREE = 3;

/*****************************************************************************************/
template<typename Mesh>
struct scalar_testing_function;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct scalar_testing_function< Mesh<T,2,Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        return std::sin(M_PI * pt.x()) * std::sin(M_PI * pt.y());
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct scalar_testing_function< Mesh<T,3,Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        return std::sin(M_PI * pt.x()) * std::sin(M_PI * pt.y()) * std::sin(M_PI * pt.z());
    }
};

template<typename Mesh>
auto make_scalar_testing_data(const Mesh& msh)
{
    return scalar_testing_function<Mesh>();
}

/*****************************************************************************************/
template<typename Mesh>
struct vector_testing_function;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct vector_testing_function< Mesh<T,2,Storage> >
{
	typedef Mesh<T,2,Storage> 				mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;
    typedef Matrix<scalar_type, 2, 1> 		ret_type;

    ret_type operator()(const point_type& pt) const
    {
    	ret_type ret;
        ret(0) = std::sin(2. * M_PI * pt.x());
        ret(1) = std::sin(2. * M_PI * pt.y());
        return ret;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct vector_testing_function< Mesh<T,3,Storage> >
{
	typedef Mesh<T,3,Storage> 				mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;
    typedef Matrix<scalar_type, 3, 1> 		ret_type;

    ret_type operator()(const point_type& pt) const
    {
    	ret_type ret;
        ret(0) = std::sin(2. * M_PI * pt.x());
        ret(1) = std::sin(2. * M_PI * pt.y());
        ret(2) = std::sin(2. * M_PI * pt.z());
        return ret;
    }
};

template<typename Mesh>
auto make_vector_testing_data(const Mesh& msh)
{
	return vector_testing_function<Mesh>();
}

/*****************************************************************************************/
template<typename Mesh>
struct vector_testing_function_div;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct vector_testing_function_div< Mesh<T,2,Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;
    typedef Matrix<scalar_type, 2, 1>       ret_type;

    scalar_type operator()(const point_type& pt) const
    {
        return 2. * M_PI * std::cos(2. * M_PI * pt.x()) +
               2. * M_PI * std::cos(2. * M_PI * pt.y());
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct vector_testing_function_div< Mesh<T,3,Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;
    typedef Matrix<scalar_type, 3, 1>       ret_type;

    scalar_type operator()(const point_type& pt) const
    {
        return 2. * M_PI * std::cos(2. * M_PI * pt.x()) +
               2. * M_PI * std::cos(2. * M_PI * pt.y()) +
               2. * M_PI * std::cos(2. * M_PI * pt.z());
    }
};

template<typename Mesh>
auto make_vector_testing_data_div(const Mesh& msh)
{
    return vector_testing_function_div<Mesh>();
}

/*****************************************************************************************/

template<typename T>
std::vector< disk::generic_mesh<T, 2> >
get_triangle_generic_meshes(void)
{
	std::vector<std::string> meshfiles;
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_1.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_2.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_3.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_4.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_5.typ1");

    typedef disk::generic_mesh<T, 2>  mesh_type;

    std::vector< mesh_type > ret;
    for (size_t i = 0; i < meshfiles.size(); i++)
    {
        mesh_type msh;
        disk::fvca5_mesh_loader<T, 2> loader;

        if (!loader.read_mesh(meshfiles.at(i)))
        {
            std::cout << "Problem loading mesh." << std::endl;
            continue;
        }
        loader.populate_mesh(msh);

        ret.push_back(msh);
    }

    return ret;
}

template<typename T>
std::vector< disk::simplicial_mesh<T, 2> >
get_triangle_netgen_meshes(void)
{
	std::vector<std::string> meshfiles;
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri01.mesh2d");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri02.mesh2d");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri03.mesh2d");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri04.mesh2d");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri05.mesh2d");


    typedef disk::simplicial_mesh<T, 2>  mesh_type;

    std::vector< mesh_type > ret;
    for (size_t i = 0; i < meshfiles.size(); i++)
    {
        mesh_type msh;
        disk::netgen_mesh_loader<T, 2> loader;

        if (!loader.read_mesh(meshfiles.at(i)))
        {
            std::cout << "Problem loading mesh." << std::endl;
            continue;
        }
        loader.populate_mesh(msh);

        ret.push_back(msh);
    }

    return ret;
}

template<typename T>
std::vector< disk::simplicial_mesh<T, 3> >
get_tetrahedra_netgen_meshes(void)
{
    std::vector<std::string> meshfiles;
    meshfiles.push_back("../../../diskpp/meshes/3D_tetras/netgen/cube1.mesh");
    meshfiles.push_back("../../../diskpp/meshes/3D_tetras/netgen/cube2.mesh");
    meshfiles.push_back("../../../diskpp/meshes/3D_tetras/netgen/cube3.mesh");
    meshfiles.push_back("../../../diskpp/meshes/3D_tetras/netgen/cube4.mesh");
    meshfiles.push_back("../../../diskpp/meshes/3D_tetras/netgen/cube5.mesh");


    typedef disk::simplicial_mesh<T, 3>  mesh_type;

    std::vector< mesh_type > ret;
    for (size_t i = 0; i < meshfiles.size(); i++)
    {
        mesh_type msh;
        disk::netgen_mesh_loader<T, 3> loader;

        if (!loader.read_mesh(meshfiles.at(i)))
        {
            std::cout << "Problem loading mesh." << std::endl;
            continue;
        }
        loader.populate_mesh(msh);

        ret.push_back(msh);
    }

    return ret;
}

template<typename T>
std::vector< disk::cartesian_mesh<T, 3> >
get_cartesian_diskpp_meshes(void)
{
    std::vector<std::string> meshfiles;
    meshfiles.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-2-2-2.hex");
    meshfiles.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-4-4-4.hex");
    meshfiles.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-8-8-8.hex");
    meshfiles.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-16-16-16.hex");
    meshfiles.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-32-32-32.hex");

    typedef disk::cartesian_mesh<T, 3>  mesh_type;

    std::vector< mesh_type > ret;
    for (size_t i = 0; i < meshfiles.size(); i++)
    {
        mesh_type msh;
        disk::cartesian_mesh_loader<T, 3> loader;

        if (!loader.read_mesh(meshfiles.at(i)))
        {
            std::cout << "Problem loading mesh." << std::endl;
            continue;
        }
        loader.populate_mesh(msh);

        ret.push_back(msh);
    }

    return ret;
}

template<typename T>
std::vector< disk::generic_mesh<T, 3> >
get_generic_fvca6_meshes(void)
{
    std::vector<std::string> meshfiles;
    meshfiles.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_10.msh");
    meshfiles.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_20.msh");
    meshfiles.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_30.msh");
    meshfiles.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_40.msh");

    typedef disk::generic_mesh<T, 3>  mesh_type;

    std::vector< mesh_type > ret;
    for (size_t i = 0; i < meshfiles.size(); i++)
    {
        mesh_type msh;
        disk::fvca6_mesh_loader<T, 3> loader;

        if (!loader.read_mesh(meshfiles.at(i)))
        {
            std::cout << "Problem loading mesh." << std::endl;
            continue;
        }
        loader.populate_mesh(msh);

        ret.push_back(msh);
    }

    return ret;
}

template<typename T>
std::vector< disk::generic_mesh<T, 2> >
get_quad_generic_meshes(void)
{
	std::vector<std::string> meshfiles;
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_1.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_2.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_3.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_4.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_5.typ1");

    typedef disk::generic_mesh<T, 2>  mesh_type;

    std::vector< mesh_type > ret;
    for (size_t i = 0; i < meshfiles.size(); i++)
    {
        mesh_type msh;
        disk::fvca5_mesh_loader<T, 2> loader;

        if (!loader.read_mesh(meshfiles.at(i)))
        {
            std::cout << "Problem loading mesh." << std::endl;
            continue;
        }
        loader.populate_mesh(msh);

        /*
        // ADD A RANDOM TRANSFORM HERE
        auto tr = [](const typename mesh_type::point_type& pt) -> auto {

            auto px = -1 * ( 1-pt.x() ) + 1 * pt.x();
            auto py = -1 * ( 1-pt.y() ) + 1 * pt.y();
            return typename mesh_type::point_type({px, py});
        };

        msh.transform(tr);
        */

        ret.push_back(msh);
    }

    return ret;
}

template<typename Mesh, typename Function>
void
do_testing(std::vector<Mesh>& meshes, const Function& run_test)
{
	using T = typename Mesh::scalar_type;

	for (size_t k = MIN_TEST_DEGREE; k <= MAX_TEST_DEGREE; k++)
    {
        std::cout << "DEGREE " << k << std::endl;

        std::vector<T> mesh_hs;
        std::vector<T> l2_errors;

        for(auto& msh : meshes)
        {
            auto error = run_test(msh, k);
            mesh_hs.push_back( disk::average_diameter(msh) );
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
}
