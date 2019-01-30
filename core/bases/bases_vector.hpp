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
 * Karol Cascavita (C) 2018                     karol.cascavita@enpc.fr
 * Nicolas Pignet  (C) 2018                     nicolas.pignet@enpc.fr
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

namespace disk
{

/* Compute the size of a vector basis of degree k in dimension d. */
size_t
vector_basis_size(size_t k, size_t sd, size_t vd)
{
    size_t num = 1;
    size_t den = 1;

    for (size_t i = 1; i <= sd; i++)
    {
        num *= k + i;
        den *= i;
    }

    return vd * (num / den);
}

/* Generic template for bases. */
template<typename MeshType, typename Element>
struct scaled_monomial_vector_basis
{
    static_assert(sizeof(MeshType) == -1, "scaled_monomial_vector_basis: not suitable for the requested kind of mesh");
    static_assert(sizeof(Element) == -1,
                  "scaled_monomial_vector_basis: not suitable for the requested kind of element");
};

/* Basis 'factory'. */
template<typename MeshType, typename ElementType>
auto
make_vector_monomial_basis(const MeshType& msh, const ElementType& elem, size_t degree)
{
    return scaled_monomial_vector_basis<MeshType, ElementType>(msh, elem, degree);
}

/* Specialization for 3D meshes, cells */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_vector_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::cell>
{

  public:
    typedef Mesh<T, 3, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::point_type      point_type;
    typedef Matrix<scalar_type, 3, 3>           gradient_type;
    typedef Matrix<scalar_type, Dynamic, 3>     function_type;
    typedef Matrix<scalar_type, Dynamic, 1>     divergence_type;

  private:
    size_t basis_degree, basis_size;

    typedef scaled_monomial_scalar_basis<mesh_type, cell_type> scalar_basis_type;
    scalar_basis_type                                          scalar_basis;

  public:
    scaled_monomial_vector_basis(const mesh_type& msh, const cell_type& cl, size_t degree) :
      scalar_basis(msh, cl, degree)
    {
        basis_degree = degree;
        basis_size   = vector_basis_size(degree, 3, 3);
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size, 3);

        auto phi = scalar_basis.eval_functions(pt);

        for (size_t i = 0; i < scalar_basis.size(); i++)
        {
            ret(3 * i, 0)     = phi(i);
            ret(3 * i + 1, 1) = phi(i);
            ret(3 * i + 2, 2) = phi(i);
        }

        return ret;
    }

    eigen_compatible_stdvector<gradient_type>
    eval_gradients(const point_type& pt) const
    {
        eigen_compatible_stdvector<gradient_type> ret;
        ret.reserve(basis_size);

        function_type dphi = scalar_basis.eval_gradients(pt);

        for (size_t i = 0; i < scalar_basis.size(); i++)
        {
            const Matrix<scalar_type, 1, 3> dphi_i = dphi.row(i);
            gradient_type                   g;

            g        = gradient_type::Zero();
            g.row(0) = dphi_i;
            ret.push_back(g);

            g        = gradient_type::Zero();
            g.row(1) = dphi_i;
            ret.push_back(g);

            g        = gradient_type::Zero();
            g.row(2) = dphi_i;
            ret.push_back(g);
        }
        assert(ret.size() == basis_size);

        return ret;
    }

    eigen_compatible_stdvector<gradient_type>
    eval_sgradients(const point_type& pt) const
    {
        eigen_compatible_stdvector<gradient_type> ret;
        ret.reserve(basis_size);

        function_type dphi = scalar_basis.eval_gradients(pt);

        for (size_t i = 0; i < scalar_basis.size(); i++)
        {
            const Matrix<scalar_type, 1, 3> dphi_i = dphi.row(i);
            gradient_type                   g;

            g        = gradient_type::Zero();
            g.row(0) = dphi_i;
            ret.push_back((g + g.transpose()) / scalar_type(2));

            g        = gradient_type::Zero();
            g.row(1) = dphi_i;
            ret.push_back((g + g.transpose()) / scalar_type(2));

            g        = gradient_type::Zero();
            g.row(2) = dphi_i;
            ret.push_back((g + g.transpose()) / scalar_type(2));
        }
        assert(ret.size() == basis_size);

        return ret;
    }

    function_type
    eval_curls(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size, 3);

        const function_type dphi = scalar_basis.eval_gradients(pt);

        size_t j = 0;
        for (size_t i = 0; i < scalar_basis.size(); i++)
        {
            const Matrix<scalar_type, 1, 3> dphi_i = dphi.row(i);
            // row 1
            ret(j, 0) = dphi_i(1);
            ret(j, 1) = dphi_i(2);
            j++;
            // row 2
            ret(j, 0) = -dphi_i(0);
            ret(j, 2) = dphi_i(2);
            j++;
            // row 2
            ret(j, 1) = -dphi_i(0);
            ret(j, 2) = -dphi_i(1);
            j++;
        }
        assert(j == basis_size);
        return ret;
    }

    divergence_type
    eval_divergences(const point_type& pt) const
    {
        divergence_type ret = divergence_type::Zero(basis_size);

        const function_type dphi = scalar_basis.eval_gradients(pt);

        for (size_t i = 0; i < scalar_basis.size(); i++)
        {
            ret(3 * i)     = dphi(i, 0);
            ret(3 * i + 1) = dphi(i, 1);
            ret(3 * i + 2) = dphi(i, 2);
        }

        return ret;
    }

    size_t
    size() const
    {
        return basis_size;
    }

    size_t
    degree() const
    {
        return basis_degree;
    }
};

/* Specialization for 3D meshes, faces */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_vector_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::face>
{

  public:
    typedef Mesh<T, 3, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;
    typedef typename mesh_type::face            face_type;
    typedef Matrix<scalar_type, Dynamic, 3>     function_type;

  private:
    size_t basis_degree, basis_size;

    typedef scaled_monomial_scalar_basis<mesh_type, face_type> scalar_basis_type;
    scalar_basis_type                                          scalar_basis;

  public:
    scaled_monomial_vector_basis(const mesh_type& msh, const face_type& fc, size_t degree) :
      scalar_basis(msh, fc, degree)
    {
        basis_degree = degree;
        basis_size   = vector_basis_size(degree, 2, 3);
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size, 3);

        const auto phi = scalar_basis.eval_functions(pt);

        for (size_t i = 0; i < scalar_basis.size(); i++)
        {
            ret(3 * i, 0)     = phi(i);
            ret(3 * i + 1, 1) = phi(i);
            ret(3 * i + 2, 2) = phi(i);
        }

        assert(3 * scalar_basis.size() == basis_size);
        return ret;
    }

    size_t
    size() const
    {
        return basis_size;
    }

    size_t
    degree() const
    {
        return basis_degree;
    }
};

/* Specialization for 2D meshes, cells */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_vector_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::cell>
{

  public:
    typedef Mesh<T, 2, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::point_type      point_type;
    typedef Matrix<scalar_type, 2, 2>           gradient_type;
    typedef Matrix<scalar_type, Dynamic, 2>     function_type;
    typedef Matrix<scalar_type, Dynamic, 1>     divergence_type;

  private:
    size_t basis_degree, basis_size;

    typedef scaled_monomial_scalar_basis<mesh_type, cell_type> scalar_basis_type;
    scalar_basis_type                                          scalar_basis;

  public:
    scaled_monomial_vector_basis(const mesh_type& msh, const cell_type& cl, size_t degree) :
      scalar_basis(msh, cl, degree)
    {
        basis_degree = degree;
        basis_size   = vector_basis_size(degree, 2, 2);
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size, 2);

        const auto phi = scalar_basis.eval_functions(pt);

        for (size_t i = 0; i < scalar_basis.size(); i++)
        {
            ret(2 * i, 0)     = phi(i);
            ret(2 * i + 1, 1) = phi(i);
        }

        return ret;
    }

    eigen_compatible_stdvector<gradient_type>
    eval_gradients(const point_type& pt) const
    {
        eigen_compatible_stdvector<gradient_type> ret;
        ret.reserve(basis_size);

        const function_type dphi = scalar_basis.eval_gradients(pt);

        for (size_t i = 0; i < scalar_basis.size(); i++)
        {
            const Matrix<scalar_type, 1, 2> dphi_i = dphi.row(i);
            gradient_type                   g;

            g        = gradient_type::Zero();
            g.row(0) = dphi_i;
            ret.push_back(g);

            g        = gradient_type::Zero();
            g.row(1) = dphi_i;
            ret.push_back(g);
        }
        assert(ret.size() == basis_size);

        return ret;
    }

    eigen_compatible_stdvector<gradient_type>
    eval_sgradients(const point_type& pt) const
    {
        eigen_compatible_stdvector<gradient_type> ret;
        ret.reserve(basis_size);

        const function_type dphi = scalar_basis.eval_gradients(pt);

        for (size_t i = 0; i < scalar_basis.size(); i++)
        {
            const Matrix<scalar_type, 1, 2> dphi_i = dphi.row(i);
            gradient_type                   g;

            g        = gradient_type::Zero();
            g.row(0) = dphi_i;
            ret.push_back(0.5 * (g + g.transpose()));

            g        = gradient_type::Zero();
            g.row(1) = dphi_i;
            ret.push_back(0.5 * (g + g.transpose()));
        }
        assert(ret.size() == basis_size);

        return ret;
    }

    Matrix<scalar_type, Dynamic, 1>
    eval_curls(const point_type& pt) const
    {
        Matrix<scalar_type, Dynamic, 1> ret = Matrix<scalar_type, Dynamic, 1>::Zero(basis_size);

        const function_type dphi = scalar_basis.eval_gradients(pt);

        size_t j = 0;
        for (size_t i = 0; i < scalar_basis.size(); i++)
        {
            Matrix<scalar_type, 1, 2> dphi_i = dphi.row(i);

            ret(j++) = dphi_i(1);
            ret(j++) = -dphi_i(0);
        }
        return ret;
    }

    divergence_type
    eval_divergences(const point_type& pt) const
    {
        divergence_type ret = divergence_type::Zero(basis_size);

        const function_type dphi = scalar_basis.eval_gradients(pt);

        for (size_t i = 0; i < scalar_basis.size(); i++)
        {
            ret(2 * i)     = dphi(i, 0);
            ret(2 * i + 1) = dphi(i, 1);
        }

        return ret;
    }

    size_t
    size() const
    {
        return basis_size;
    }

    size_t
    degree() const
    {
        return basis_degree;
    }
};

/* Specialization for 2D meshes, faces */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_vector_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::face>
{

  public:
    typedef Mesh<T, 2, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;
    typedef typename mesh_type::face            face_type;
    typedef Matrix<scalar_type, Dynamic, 2>     function_type;

  private:
    size_t basis_degree, basis_size;

    typedef scaled_monomial_scalar_basis<mesh_type, face_type> scalar_basis_type;
    scalar_basis_type                                          scalar_basis;

  public:
    scaled_monomial_vector_basis(const mesh_type& msh, const face_type& fc, size_t degree) :
      scalar_basis(msh, fc, degree)
    {
        basis_degree = degree;
        basis_size   = vector_basis_size(degree, 1, 2);
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size, 2);

        const auto phi = scalar_basis.eval_functions(pt);

        for (size_t i = 0; i < scalar_basis.size(); i++)
        {
            ret(2 * i, 0)     = phi(i);
            ret(2 * i + 1, 1) = phi(i);
        }

        assert(2 * scalar_basis.size() == basis_size);
        return ret;
    }

    size_t
    size() const
    {
        return basis_size;
    }

    size_t
    degree() const
    {
        return basis_degree;
    }
};

///////////////////////////////////////////////////
//// Raviart-Thomas elements on simplicial ////////
///////////////////////////////////////////////////

/* Compute the size of a vector basis of degree k in dimension d. */
size_t
vector_basis_size_RT(const size_t k, const size_t sd, const size_t vd)
{
    if(k == 0)
        throw std::invalid_argument("Raviart-Thomas basis: degree has to be > 0");

    if (sd == 3 && vd == 3)
    {
        return (k + 1) * (k + 2) * (k + 4) / 2;
    }
    else if (sd == 2 && vd == 2)
    {
        return (k + 1) * (k + 3);
    }

    throw std::invalid_argument("Raviart-Thomas basis: unknown case");

    return 0;
}

// RT^{k}(T; R^d) = P^{k-1}(T;R^d) + x * P^{k-1,H}(T;R)

/* Generic template for RT bases. */
template<typename MeshType, typename Element>
struct scaled_monomial_vector_basis_RT
{
    static_assert(sizeof(MeshType) == -1, "scaled_monomial_vector_basis_RT: not suitable for the requested kind of mesh");
    static_assert(sizeof(Element) == -1,
                  "scaled_monomial_vector_basis_RT: not suitable for the requested kind of element");
};

/* RT Basis 'factory'. */
template<typename MeshType, typename ElementType>
auto
make_vector_monomial_basis_RT(const MeshType& msh, const ElementType& elem, const size_t degree)
{
    return scaled_monomial_vector_basis_RT<MeshType, ElementType>(msh, elem, degree);
}

/* Specialization for 3D meshes, cells */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_vector_basis_RT<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::cell>
{

  public:
    typedef Mesh<T, 3, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::point_type      point_type;
    typedef Matrix<scalar_type, Dynamic, 3>     function_type;
    typedef Matrix<scalar_type, Dynamic, 1>     divergence_type;

  private:
    size_t basis_degree, basis_size;

    typedef scaled_monomial_scalar_basis<mesh_type, cell_type> scalar_basis_type;
    scalar_basis_type                                          scalar_basis;

    typedef scaled_monomial_vector_basis<mesh_type, cell_type> vector_basis_type;
    vector_basis_type                                        vector_basis;

  public:
    scaled_monomial_vector_basis_RT(const mesh_type& msh, const cell_type& cl, const size_t degree) :
      basis_degree(degree), basis_size(vector_basis_size_RT(degree, 3, 3)), scalar_basis(msh, cl, degree - 1),
      vector_basis(msh, cl, degree - 1)
    {
        if (degree <= 0)
            throw std::invalid_argument("Raviart-Thomas basis: degree has to be > 0");

        if (points(msh, cl).size() != 4)
            throw std::invalid_argument("Raviart-Thomas basis: available only on tetrahedron");
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size, 3);

        ret.block(0,0, vector_basis.size(), 3) = vector_basis.eval_functions(pt);

        const auto sphi = scalar_basis.eval_functions(pt);
        size_t     deg    = 0;
        if(basis_degree > 2)
            deg = basis_degree - 2;
        const auto beg    = scalar_basis_size(deg, 3);
        const auto offset = vector_basis.size();

        // compute x P^(k-1)_H (monomial of degree exactly k - 1)
        for (size_t i = beg; i < scalar_basis.size(); i++)
        {
            ret(offset + i - beg, 0)     = pt.x() * sphi(i);
            ret(offset + i - beg, 1)     = pt.y() * sphi(i);
            ret(offset + i - beg, 2)     = pt.z() * sphi(i);
        }

        return ret;
    }

    divergence_type
    eval_divergences(const point_type& pt) const
    {
        divergence_type ret = divergence_type::Zero(basis_size);

        ret.head(vector_basis.size()) = vector_basis.eval_divergences(pt);

        const auto sphi   = scalar_basis.eval_functions(pt);
        const auto sdphi  = scalar_basis.eval_gradients(pt);
        size_t     deg    = 0;
        if (basis_degree > 2)
            deg = basis_degree - 2;
        const auto beg    = scalar_basis_size(deg, 3);
        const auto offset = vector_basis.size();

        /// compute P^(k-1)_H + x.grad(P^(k-1)_H) (monomial of degree exactly k - 1)
        for (size_t i = beg; i < scalar_basis.size(); i++)
        {
            ret(offset + i - beg) = 3 * sphi(i) + pt.x() * sdphi(i, 0) + pt.y() * sdphi(i, 1) + pt.z() * sdphi(i, 2);
        }

        return ret;
    }

    size_t
    size() const
    {
        return basis_size;
    }

    size_t
    degree() const
    {
        return basis_degree;
    }
};

/* Specialization for 2D meshes, cells */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_vector_basis_RT<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::cell>
{

  public:
    typedef Mesh<T, 2, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::point_type      point_type;
    typedef Matrix<scalar_type, Dynamic, 2>     function_type;
    typedef Matrix<scalar_type, Dynamic, 1>     divergence_type;

  private:
    size_t basis_degree, basis_size;

    typedef scaled_monomial_scalar_basis<mesh_type, cell_type> scalar_basis_type;
    scalar_basis_type                                          scalar_basis;

    typedef scaled_monomial_vector_basis<mesh_type, cell_type> vector_basis_type;
    vector_basis_type                                        vector_basis;

  public:
    scaled_monomial_vector_basis_RT(const mesh_type& msh, const cell_type& cl, const size_t degree) :
      basis_degree(degree), basis_size(vector_basis_size_RT(degree, 2, 2)), scalar_basis(msh, cl, degree - 1),
      vector_basis(msh, cl, degree - 1)
    {
        if (degree <= 0)
            throw std::invalid_argument("Raviart-Thomas basis: degree has to be > 0");

        if (points(msh, cl).size() != 3)
            throw std::invalid_argument("Raviart-Thomas basis: available only on triangles");
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size, 2);

        ret.block(0, 0, vector_basis.size(), 2) = vector_basis.eval_functions(pt);

        const auto sphi   = scalar_basis.eval_functions(pt);
        size_t     deg    = 0;
        if (basis_degree > 2)
            deg = basis_degree - 2;
        const auto beg    = scalar_basis_size(deg, 2);
        const auto offset = vector_basis.size();

        // compute x P^(k-1)_H (monomial of degree exactly k - 1)
        for (size_t i = beg; i < scalar_basis.size(); i++)
        {
            ret(offset + i - beg, 0) = pt.x() * sphi(i);
            ret(offset + i - beg, 1) = pt.y() * sphi(i);
        }

        return ret;
    }

    divergence_type
    eval_divergences(const point_type& pt) const
    {
        divergence_type ret = divergence_type::Zero(basis_size);

        ret.head(vector_basis.size()) = vector_basis.eval_divergences(pt);

        const auto sphi   = scalar_basis.eval_functions(pt);
        const auto sdphi  = scalar_basis.eval_gradients(pt);
        size_t     deg    = 0;
        if (basis_degree > 2)
            deg = basis_degree - 2;
        const auto beg    = scalar_basis_size(deg, 2);
        const auto offset = vector_basis.size();

        /// compute P^(k-1)_H + x.grad(P^(k-1)_H) (monomial of degree exactly k - 1)
        for (size_t i = beg; i < scalar_basis.size(); i++)
        {
            ret(offset + i - beg) = 2 * sphi(i) + pt.x() * sdphi(i, 0) + pt.y() * sdphi(i, 1);
        }

        return ret;
    }

    size_t
    size() const
    {
        return basis_size;
    }

    size_t
    degree() const
    {
        return basis_degree;
    }
};


}
