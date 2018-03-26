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

namespace revolution {

/* Compute the size of a matrix basis of degree k in dimension d. */
size_t
matrix_basis_size(size_t k, size_t sd, size_t md)
{
   size_t num = 1;
   size_t den = 1;

   for (size_t i = 1; i <= sd; i++) {
      num *= k + i;
      den *= i;
   }

   return md * md * (num / den);
}

/* Generic template for bases. */
template<typename MeshType, typename Element>
struct scaled_monomial_matrix_basis
{
   static_assert(sizeof(MeshType) == -1,
                 "scaled_monomial_matrix_basis: not suitable for the requested kind of mesh");
   static_assert(sizeof(Element) == -1,
                 "scaled_monomial_matrix_basis: not suitable for the requested kind of element");
};

/* Basis 'factory'. */
template<typename MeshType, typename ElementType>
auto
make_matrix_monomial_basis(const MeshType& msh, const ElementType& elem, size_t degree)
{
   return scaled_monomial_matrix_basis<MeshType, ElementType>(msh, elem, degree);
}

/* Specialization for 3D meshes, cells */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_matrix_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::cell>
{

 public:
   typedef Mesh<T, 3, Storage>              mesh_type;
   typedef typename mesh_type::scalar_type  scalar_type;
   typedef typename mesh_type::cell         cell_type;
   typedef typename mesh_type::point_type   point_type;
   typedef static_matrix<scalar_type, 3, 3> function_type;

 private:
   size_t      basis_degree, basis_size;

   typedef scaled_monomial_scalar_basis<mesh_type, cell_type> scalar_basis_type;
   scalar_basis_type                                          scalar_basis;

 public:
   scaled_monomial_matrix_basis(const mesh_type& msh, const cell_type& cl, size_t degree) :
     scalar_basis(msh, cl, degree)
   {
      basis_degree = degree;
      basis_size   = matrix_basis_size(degree, 3, 3);
   }

   eigen_compatible_stdvector<function_type>
   eval_functions(const point_type& pt) const
   {
      eigen_compatible_stdvector<function_type> ret;
      ret.reserve(basis_size);

      const auto phi = scalar_basis.eval_functions(pt);

      for (size_t k = 0; k < scalar_basis.size(); k++) {
         function_type fc;

         for (size_t j = 0; j < 3; j++) {
            for (size_t i = 0; i < 3; i++) {
               fc       = function_type::Zero();
               fc(i, j) = phi(k);
               ret.push_back(fc);
            }
         }
      }

      assert(ret.size() == basis_size);

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
class scaled_monomial_matrix_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::face>
{

 public:
   typedef Mesh<T, 3, Storage>              mesh_type;
   typedef typename mesh_type::scalar_type  scalar_type;
   typedef typename mesh_type::point_type   point_type;
   typedef typename mesh_type::face         face_type;
   typedef static_matrix<scalar_type, 3, 3> function_type;

 private:
   size_t      basis_degree, basis_size;

   typedef scaled_monomial_scalar_basis<mesh_type, face_type> scalar_basis_type;
   scalar_basis_type                                          scalar_basis;

 public:
   scaled_monomial_matrix_basis(const mesh_type& msh, const face_type& fc, size_t degree) :
     scalar_basis(msh, fc, degree)
   {
      basis_degree = degree;
      basis_size   = matrix_basis_size(degree, 2, 3);
   }

   eigen_compatible_stdvector<function_type>
   eval_functions(const point_type& pt) const
   {
      eigen_compatible_stdvector<function_type> ret;
      ret.reserve(basis_size);

      const auto phi = scalar_basis.eval_functions(pt);

      for (size_t k = 0; k < scalar_basis.size(); k++) {
         function_type fc;

         for (size_t j = 0; j < 3; j++) {
            for (size_t i = 0; i < 3; i++) {
               fc       = function_type::Zero();
               fc(i, j) = phi(k);
               ret.push_back(fc);
            }
         }
      }

      assert(ret.size() == basis_size);
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
class scaled_monomial_matrix_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::cell>
{

 public:
   typedef Mesh<T, 2, Storage>              mesh_type;
   typedef typename mesh_type::scalar_type  scalar_type;
   typedef typename mesh_type::cell         cell_type;
   typedef typename mesh_type::point_type   point_type;
   typedef static_matrix<scalar_type, 2, 2> function_type;

 private:
   size_t      basis_degree, basis_size;

   typedef scaled_monomial_scalar_basis<mesh_type, cell_type> scalar_basis_type;
   scalar_basis_type                                          scalar_basis;

 public:
   scaled_monomial_matrix_basis(const mesh_type& msh, const cell_type& cl, size_t degree) :
     scalar_basis(msh, cl, degree)
   {
      basis_degree = degree;
      basis_size   = matrix_basis_size(degree, 2, 2);
   }

   eigen_compatible_stdvector<function_type>
   eval_functions(const point_type& pt) const
   {
      eigen_compatible_stdvector<function_type> ret;
      ret.reserve(basis_size);

      const auto phi = scalar_basis.eval_functions(pt);

      for (size_t k = 0; k < scalar_basis.size(); k++) {
         function_type fc;

         for (size_t j = 0; j < 2; j++) {
            for (size_t i = 0; i < 2; i++) {
               fc       = function_type::Zero();
               fc(i, j) = phi(k);
               ret.push_back(fc);
            }
         }
      }

      assert(ret.size() == basis_size);

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
class scaled_monomial_matrix_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::face>
{

 public:
   typedef Mesh<T, 2, Storage>              mesh_type;
   typedef typename mesh_type::scalar_type  scalar_type;
   typedef typename mesh_type::point_type   point_type;
   typedef typename mesh_type::face         face_type;
   typedef static_matrix<scalar_type, 2, 2> function_type;

 private:
   size_t      basis_degree, basis_size;

   typedef scaled_monomial_scalar_basis<mesh_type, face_type> scalar_basis_type;
   scalar_basis_type                                          scalar_basis;

 public:
   scaled_monomial_matrix_basis(const mesh_type& msh, const face_type& fc, size_t degree) :
     scalar_basis(msh, fc, degree)
   {
      basis_degree = degree;
      basis_size   = matrix_basis_size(degree, 1, 2);
   }

   eigen_compatible_stdvector<function_type>
   eval_functions(const point_type& pt) const
   {
      eigen_compatible_stdvector<function_type> ret;
      ret.reserve(basis_size);

      const auto phi = scalar_basis.eval_functions(pt);

      for (size_t k = 0; k < scalar_basis.size(); k++) {
         function_type fc;

         for (size_t j = 0; j < 2; j++) {
            for (size_t i = 0; i < 2; i++) {
               fc       = function_type::Zero();
               fc(i, j) = phi(k);
               ret.push_back(fc);
            }
         }
      }

      assert(ret.size() == basis_size);
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

/* Compute the size of a matrix basis of degree k in dimension d. */
size_t
sym_matrix_basis_size(size_t k, size_t sd, size_t vd)
{
    size_t num = 1;
    size_t den = 1;

    for (size_t i = 1; i <= sd; i++)
    {
        num *= k + i;
        den *= i;
    }

    size_t md = 0;

    switch (vd)
    {
        case 3:
            md = 6;
            break;

        case 2:
            md = 3;
            break;

        default:
            std::logic_error("Expected 3 >= dim > 1");
    }

    return md * (num / den);
}

/* Generic template for bases. */
template<typename MeshType, typename Element>
struct scaled_monomial_sym_matrix_basis
{
    static_assert(sizeof(MeshType) == -1,   "scaled_monomial_sym_matrix_basis: not "
                                            "suitable for the requested kind of mesh");
    static_assert( sizeof(Element) == -1,   "scaled_monomial_sym_matrix_basis: not "
                                            "suitable for the requested kind of element");
};

/* Basis 'factory'. */
template<typename MeshType, typename ElementType>
auto
make_sym_matrix_monomial_basis(const MeshType& msh, const ElementType& elem, size_t degree)
{
    return scaled_monomial_sym_matrix_basis<MeshType, ElementType>(msh, elem, degree);
}

/* Specialization for 3D meshes, cells */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_sym_matrix_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::cell>
{

public:
    typedef Mesh<T, 3, Storage>              mesh_type;
    typedef typename mesh_type::scalar_type  scalar_type;
    typedef typename mesh_type::cell         cell_type;
    typedef typename mesh_type::point_type   point_type;
    typedef static_matrix<scalar_type, 3, 3> function_type;

private:
    size_t      basis_degree, basis_size;

    typedef scaled_monomial_scalar_basis<mesh_type, cell_type> scalar_basis_type;
    scalar_basis_type                                          scalar_basis;

public:
    scaled_monomial_sym_matrix_basis(const mesh_type& msh, const cell_type& cl, size_t degree)
        : scalar_basis(msh, cl, degree)
    {
        basis_degree = degree;
        basis_size   = sym_matrix_basis_size(degree, 3, 3);
    }

    eigen_compatible_stdvector<function_type>
    eval_functions(const point_type& pt) const
    {
        eigen_compatible_stdvector<function_type> ret;
        ret.reserve(basis_size);

        auto phi = scalar_basis.eval_functions(pt);

        for (size_t k = 0; k < scalar_basis.size(); k++)
        {
            function_type fc;

            for (size_t j = 0; j < 3; j++)
            {
                for (size_t i = 0; i < j; i++)
                {
                    fc       = function_type::Zero();
                    fc(i, j) = phi(k);
                    fc(j, i) = phi(k);
                    ret.push_back(fc);
                }
            
                fc       = function_type::Zero();
                fc(j, j) = phi(k);
                ret.push_back(fc);
            }
        }

        assert(ret.size() == basis_size);
        return ret;
    }

    size_t
    size() const {
        return basis_size;
    }

    size_t
    degree() const {
        return basis_degree;
    }
};

/* Specialization for 3D meshes, faces */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_sym_matrix_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::face>
{

public:
    typedef Mesh<T, 3, Storage>              mesh_type;
    typedef typename mesh_type::scalar_type  scalar_type;
    typedef typename mesh_type::point_type   point_type;
    typedef typename mesh_type::face         face_type;
    typedef static_matrix<scalar_type, 3, 3> function_type;

private:
    size_t      basis_degree, basis_size;

    typedef scaled_monomial_scalar_basis<mesh_type, face_type> scalar_basis_type;
    scalar_basis_type                                          scalar_basis;

public:
    scaled_monomial_sym_matrix_basis(const mesh_type& msh, const face_type& fc, size_t degree)
        : scalar_basis(msh, fc, degree)
    {
        basis_degree = degree;
        basis_size   = sym_matrix_basis_size(degree, 2, 3);
    }

    eigen_compatible_stdvector<function_type>
    eval_functions(const point_type& pt) const
    {
        eigen_compatible_stdvector<function_type> ret;
        ret.reserve(basis_size);

        auto phi = scalar_basis.eval_functions(pt);

        for (size_t k = 0; k < scalar_basis.size(); k++)
        {
            function_type fc;

            for (size_t j = 0; j < 3; j++)
            {
                for (size_t i = 0; i < j; i++)
                {
                    fc       = function_type::Zero();
                    fc(i, j) = phi(k);
                    fc(j, i) = phi(k);
                    ret.push_back(fc);
                }
                fc       = function_type::Zero();
                fc(j, j) = phi(k);
                ret.push_back(fc);
            }
        }

        assert(ret.size() == basis_size);
        return ret;
    }
   
    size_t
    size() const {
        return basis_size;
    }

    size_t degree() const {
        return basis_degree;
    }
};

/* Specialization for 2D meshes, cells */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_sym_matrix_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::cell>
{
public:
    typedef Mesh<T, 2, Storage>              mesh_type;
    typedef typename mesh_type::scalar_type  scalar_type;
    typedef typename mesh_type::cell         cell_type;
    typedef typename mesh_type::point_type   point_type;
    typedef static_matrix<scalar_type, 2, 2> function_type;

private:
    size_t      basis_degree, basis_size;

    typedef scaled_monomial_scalar_basis<mesh_type, cell_type> scalar_basis_type;
    scalar_basis_type                                          scalar_basis;

 public:
    scaled_monomial_sym_matrix_basis(const mesh_type& msh, const cell_type& cl, size_t degree)
        : scalar_basis(msh, cl, degree)
    {
        basis_degree = degree;
        basis_size   = sym_matrix_basis_size(degree, 2, 2);
    }

    eigen_compatible_stdvector<function_type>
    eval_functions(const point_type& pt) const
    {
        eigen_compatible_stdvector<function_type> ret;
        ret.reserve(basis_size);

        auto phi = scalar_basis.eval_functions(pt);

        for (size_t k = 0; k < scalar_basis.size(); k++)
        {
            function_type fc;

            for (size_t j = 0; j < 2; j++)
            {
                for (size_t i = 0; i < j; i++)
                {
                    fc       = function_type::Zero();
                    fc(i, j) = phi(k);
                    fc(j, i) = phi(k);
                    ret.push_back(fc);
                }
            
                fc       = function_type::Zero();
                fc(j, j) = phi(k);
                ret.push_back(fc);
            }
        }

        assert(ret.size() == basis_size);
        return ret;
    }

    size_t
    size() const {
        return basis_size;
    }

    size_t degree() const {
        return basis_degree;
    }
};

/* Specialization for 2D meshes, faces */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_sym_matrix_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::face>
{
public:
    typedef Mesh<T, 2, Storage>              mesh_type;
    typedef typename mesh_type::scalar_type  scalar_type;
    typedef typename mesh_type::point_type   point_type;
    typedef typename mesh_type::face         face_type;
    typedef static_matrix<scalar_type, 2, 2> function_type;

private:
    size_t      basis_degree, basis_size;

    typedef scaled_monomial_scalar_basis<mesh_type, face_type> scalar_basis_type;
    scalar_basis_type                                          scalar_basis;

public:
    scaled_monomial_sym_matrix_basis(const mesh_type& msh, const face_type& fc, size_t degree)
        : scalar_basis(msh, fc, degree)
    {
        basis_degree = degree;
        basis_size   = sym_matrix_basis_size(degree, 1, 2);
    }

    eigen_compatible_stdvector<function_type>
    eval_functions(const point_type& pt) const
    {
        eigen_compatible_stdvector<function_type> ret;
        ret.reserve(basis_size);

        const phi = scalar_basis.eval_functions(pt);

        for (size_t k = 0; k < scalar_basis.size(); k++)
        {
            function_type fc;

            for (size_t j = 0; j < 2; j++)
            {
                for (size_t i = 0; i < j; i++)
                {
                    fc       = function_type::Zero();
                    fc(i, j) = phi(k);
                    fc(j, i) = phi(k);
                    ret.push_back(fc);
                }
                fc       = function_type::Zero();
                fc(j, j) = phi(k);
                ret.push_back(fc);
            }
        }

        assert(ret.size() == basis_size);
        return ret;
    }

    size_t
    size() const {
      return basis_size;
    }

    size_t
    degree() const {
        return basis_degree;
    }
};

} // namespace revolution
