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

#pragma once

#include <vector>

#include "common/eigen.hpp"

using namespace Eigen;

namespace revolution
{

/* Compute the size of a vector basis of degree k in dimension d. */
size_t vector_basis_size(size_t k, size_t sd, size_t vd)
{
    size_t num = 1;
    size_t den = 1;

    for (size_t i = 1; i <= sd; i++)
    {
        num *= k + i;
        den *= i;
    }

    return vd * (num/den);
}

/* Generic template for bases. */
template<typename MeshType, typename Element>
struct scaled_monomial_vector_basis
{
    static_assert(sizeof(MeshType) == -1, "scaled_monomial_vector_basis: not suitable for the requested kind of mesh");
    static_assert(sizeof(Element) == -1, "scaled_monomial_vector_basis: not suitable for the requested kind of element");
};

/* Basis 'factory'. */
template<typename MeshType, typename ElementType>
auto make_vector_monomial_basis(const MeshType& msh, const ElementType& elem, size_t degree)
{
    return scaled_monomial_vector_basis<MeshType, ElementType>(msh, elem, degree);
}

/* Specialization for 2D meshes, cells */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_vector_basis<Mesh<T,2,Storage>, typename Mesh<T,2,Storage>::cell>
{

public:
    typedef Mesh<T,2,Storage>                       mesh_type;
    typedef typename mesh_type::scalar_type         scalar_type;
    typedef typename mesh_type::cell                cell_type;
    typedef typename mesh_type::point_type          point_type;
    typedef Matrix<scalar_type, 2, 2>               gradient_type;
    typedef Matrix<scalar_type, Dynamic, 2>         function_type;

private:
    point_type          cell_bar;
    scalar_type         cell_h;
    size_t              basis_degree, basis_size;

    typedef scaled_monomial_scalar_basis<mesh_type, cell_type> scalar_basis_type;
    scalar_basis_type  scalar_basis;

public:
    scaled_monomial_vector_basis(const mesh_type& msh, const cell_type& cl, size_t degree)
    : scalar_basis(msh, cl, degree)
    {
        cell_bar        = barycenter(msh, cl);
        cell_h          = diameter(msh, cl);
        basis_degree    = degree;
        basis_size      = vector_basis_size(degree, 2, 2);
    }

    auto
    eval_functions(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size);

        auto phi = scalar_basis.eval_functions(pt);
        auto dim = mesh_type::dimension;

        for(size_t j = 0; j < dim; j++)
            for(size_t i = 0; i < scalar_basis.size(); i++ )
                ret(j, dim * i + j) = phi(i);

        return ret;
    }

    auto
    eval_gradients(const point_type& pt) const
    {
        std::vector<gradient_type> ret;
        ret.reserve(basis_size);

        auto dphi = scalar_basis.eval_functions(pt);

        size_t j = 0;
        for(size_t i = 0; i < scalar_basis.size(); i++)
        {
            auto dphi_i = dphi.row(i);
            gradient_type sg;

            sg = gradient_type::Zero();
            sg.row(0) = dphi_i;
            ret.push_back(sg);

            sg = gradient_type::Zero();
            sg.row(1) = dphi_i;
            ret.push_back(sg);
        }
        return ret;
    }

    size_t size() const
    {
        return basis_size;
    }

    size_t degree() const
    {
        return basis_degree;
    }
};

/* Specialization for 2D meshes, faces */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_vector_basis<Mesh<T,2,Storage>, typename Mesh<T,2,Storage>::face>
{

public:
    typedef Mesh<T,2,Storage>                       mesh_type;
    typedef typename mesh_type::scalar_type         scalar_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::face                face_type;
    typedef Matrix<scalar_type, Dynamic, 2>         function_type;

private:
    point_type           face_bar;
    point_type           base;
    scalar_type          face_h;
    size_t basis_degree, basis_size;

    typedef scaled_monomial_scalar_basis<mesh_type, face_type> scalar_basis_type;
    scalar_basis_type  scalar_basis;


public:
    scaled_monomial_vector_basis(const mesh_type& msh, const face_type& fc, size_t degree)
        : scalar_basis(msh, fc, degree)
    {
        face_bar        = barycenter(msh, fc);
        face_h          = diameter(msh, fc);
        basis_degree    = degree;
        basis_size      = vector_basis_size( degree, 1, 2);
    }

    auto
    eval_functions(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size);

        auto phi = scalar_basis.eval_functions(pt);
        auto dim = mesh_type::dimension;

        for(size_t j = 0; j < dim; j++)
            for(size_t i = 0; i < scalar_basis.size(); i++ )
                ret(j, dim * i + j) = phi(i);

        assert( 2 * scalar_basis == basis_size);
        return ret;
    }

    size_t size() const
    {
        return basis_size;
    }

    size_t degree() const
    {
        return basis_degree;
    }
};

}
