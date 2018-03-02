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

#include "core/loaders/loader.hpp"

const size_t MIN_TEST_DEGREE = 0;
const size_t MAX_TEST_DEGREE = 5;

/*****************************************************************************************/
template<typename Mesh>
struct scalar_testing_function;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct scalar_testing_function< Mesh<T,2,Storage> >
{
	typedef Mesh<T,2,Storage> 				mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
    	return std::sin(2. * M_PI * pt.x()) * std::sin(2. * M_PI * pt.y());
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

template<typename Mesh>
auto make_vector_testing_data(const Mesh& msh)
{
	return vector_testing_function<Mesh>();
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


    std::vector< disk::generic_mesh<T, 2> > ret;
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


    std::vector< disk::generic_mesh<T, 2> > ret;
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
}

