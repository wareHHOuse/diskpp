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

#include "diskpp/bases/bases_all.hpp"
#include "diskpp/common/eigen.hpp"
#include "diskpp/mesh/mesh.hpp"

using namespace Eigen;

namespace disk
{

/* Compute the size of a scalar basis of degree k in dimension d. */
inline size_t
scalar_basis_size(size_t k, size_t d)
{
    size_t num = 1;
    size_t den = 1;

    for (size_t i = 1; i <= d; i++)
    {
        num *= k + i;
        den *= i;
    }

    return num / den;
}

/* Generic template for bases. */
template<typename MeshType, typename Element, typename ScalarType>
struct scaled_monomial_scalar_basis
{
    static_assert(sizeof(MeshType) == -1, "scaled_monomial_scalar_basis: not suitable for the requested kind of mesh");
    static_assert(sizeof(Element) == -1,
                  "scaled_monomial_scalar_basis: not suitable for the requested kind of element");
};

/* Basis 'factory'. */
#ifndef USE_LEGENDRE
template<typename MeshType, typename ElementType, typename ScalarType = typename MeshType::coordinate_type>
auto
make_scalar_monomial_basis(const MeshType&    msh,
                           const ElementType& elem,
                           size_t             degree,
                           bool               use_inertia_axes = USE_INERTIA_AXES)
{
    return scaled_monomial_scalar_basis<MeshType, ElementType, ScalarType>(msh, elem, degree, use_inertia_axes);
}
#endif

/***************************************************************************************************/
/***************************************************************************************************/
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename ScalarType>
class scaled_monomial_scalar_basis<Mesh<T, 1, Storage>, typename Mesh<T, 1, Storage>::cell, ScalarType>
{
public:
    typedef Mesh<T, 1, Storage>                     mesh_type;
    typedef ScalarType                              scalar_type;
    typedef typename mesh_type::coordinate_type     coordinate_type;

    typedef typename mesh_type::cell                cell_type;
    typedef typename mesh_type::point_type          point_type;
    typedef Matrix<scalar_type, Dynamic, 1>         function_type;
    typedef Matrix<scalar_type, Dynamic, 1>         gradient_type;

private:
    size_t          basis_degree, basis_size;
    point_type      bar;
    coordinate_type h;

public:
    scaled_monomial_scalar_basis(const mesh_type& msh, const cell_type& cl, size_t degree, bool)
    {
        basis_degree = degree;
        basis_size = degree + 1;
        bar = barycenter(msh, cl);
        h = diameter(msh, cl);
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size);

        const auto bp = 2*(pt - bar)/h;

        ret(0) = 1;
        for (size_t k = 1; k < basis_size; k++)
            ret(k) = ret(k-1)*bp.x();

        return ret;
    }

    gradient_type
    eval_gradients(const point_type& pt) const
    {
        gradient_type ret = gradient_type::Zero(basis_size);

        const auto bp = 2*(pt - bar)/h;

        for (size_t k = 1; k < basis_size; k++)
            ret(k) = (2*k/h)*iexp_pow(bp.x(), k-1);

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

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename ScalarType>
class scaled_monomial_scalar_basis<Mesh<T, 1, Storage>, typename Mesh<T, 1, Storage>::face, ScalarType>
{
public:
    typedef Mesh<T, 1, Storage>                     mesh_type;
    typedef ScalarType                              scalar_type;
    typedef typename mesh_type::coordinate_type     coordinate_type;

    typedef typename mesh_type::face                face_type;
    typedef typename mesh_type::point_type          point_type;
    typedef Matrix<scalar_type, Dynamic, 1>         function_type;

private:
    size_t          basis_degree, basis_size;

public:
    scaled_monomial_scalar_basis(const mesh_type& msh, const face_type& cl, size_t degree, bool)
    {
        basis_degree = degree;
        basis_size = 1;
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size);
        ret(0) = 1;
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

/* Specialization for 2D meshes, cells */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename ScalarType>
class scaled_monomial_scalar_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::cell, ScalarType>
  : public scaled_monomial_abstract_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::cell, ScalarType>
{

  public:
    typedef Mesh<T, 2, Storage>                 mesh_type;
    typedef ScalarType                          scalar_type;
    typedef typename mesh_type::coordinate_type coordinate_type;
    ;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::point_type  point_type;
    typedef Matrix<scalar_type, Dynamic, 2> gradient_type;
    typedef Matrix<scalar_type, Dynamic, 1> function_type;

    using base = scaled_monomial_abstract_basis<mesh_type, cell_type, scalar_type>;

  private:
    size_t basis_degree, basis_size;

  public:
    scaled_monomial_scalar_basis(const mesh_type& msh,
                                 const cell_type& cl,
                                 size_t           degree,
                                 bool             use_inertia_axes = false) :
      base(msh, cl, use_inertia_axes)
    {
        basis_degree = degree;
        basis_size   = scalar_basis_size(degree, 2);
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size);

        const auto bp = this->scaling_point(pt);

        size_t pos = 0;
        for (size_t k = 0; k <= basis_degree; k++)
        {
            for (size_t i = 0; i <= k; i++)
            {
                const auto pow_x = k - i;
                const auto pow_y = i;

                const auto px = iexp_pow(bp.x(), pow_x);
                const auto py = iexp_pow(bp.y(), pow_y);

                ret(pos++) = px * py;
            }
        }

        assert(pos == basis_size);

        return ret;
    }

    gradient_type
    eval_gradients(const point_type& pt) const
    {
        gradient_type ret = gradient_type::Zero(basis_size, 2);

        const auto bp = this->scaling_point(pt);

        const auto bx = bp.x();
        const auto by = bp.y();

        const auto                        passage = this->passage_new2old();
        static_vector<coordinate_type, 2> grad_new, grad;

        size_t pos = 0;
        for (size_t k = 0; k <= basis_degree; k++)
        {
            for (size_t i = 0; i <= k; i++)
            {
                const auto pow_x = k - i;
                const auto pow_y = i;

                const auto px = iexp_pow(bx, pow_x);
                const auto py = iexp_pow(by, pow_y);
                const auto dx = (pow_x == 0) ? 0 : pow_x * iexp_pow(bx, pow_x - 1);
                const auto dy = (pow_y == 0) ? 0 : pow_y * iexp_pow(by, pow_y - 1);

                grad_new(0) = dx * py;
                grad_new(1) = px * dy;
                grad        = passage * grad_new;

                ret(pos, 0) = grad(0);
                ret(pos, 1) = grad(1);

                pos++;
            }
        }

        assert(pos == basis_size);

        return ret;
    }

    gradient_type
    eval_curls2(const point_type& pt) const
    {
        gradient_type ret = gradient_type::Zero(basis_size, 2);

        const auto dphi = eval_gradients(pt);

        for (size_t i = 0; i < basis_size; i++)
        {
            ret(i, 0) = dphi(i, 1);
            ret(i, 1) = -dphi(i, 0);
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
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename ScalarType>
class scaled_monomial_scalar_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::face, ScalarType>
  : public scaled_monomial_abstract_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::face, ScalarType>
{

  public:
    typedef Mesh<T, 2, Storage>                 mesh_type;
    typedef ScalarType                          scalar_type;
    typedef typename mesh_type::coordinate_type coordinate_type;
    typedef typename mesh_type::point_type      point_type;
    typedef typename mesh_type::face            face_type;
    typedef Matrix<scalar_type, Dynamic, 1>     function_type;

    using base = scaled_monomial_abstract_basis<mesh_type, face_type, scalar_type>;

  private:
    size_t basis_degree, basis_size;

  public:
    scaled_monomial_scalar_basis(const mesh_type& msh,
                                 const face_type& fc,
                                 size_t           degree,
                                 bool             use_inertia_axes = false) :
      base(msh, fc, use_inertia_axes)
    {
        basis_degree = degree;
        basis_size   = degree + 1;
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size);

        const auto ep = this->scaling_point(pt);

        for (size_t i = 0; i <= basis_degree; i++)
        {
            ret(i) = iexp_pow(ep.x(), i);
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

/***************************************************************************************************/
/***************************************************************************************************/

/* Specialization for 3D meshes, cells */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename ScalarType>
class scaled_monomial_scalar_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::cell, ScalarType>
  : public scaled_monomial_abstract_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::cell, ScalarType>
{

  public:
    typedef Mesh<T, 3, Storage>                 mesh_type;
    typedef ScalarType                          scalar_type;
    typedef typename mesh_type::coordinate_type coordinate_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::point_type      point_type;
    typedef Matrix<scalar_type, Dynamic, 3>     gradient_type;
    typedef Matrix<scalar_type, Dynamic, 1>     function_type;

    using base = scaled_monomial_abstract_basis<mesh_type, cell_type, scalar_type>;

  private:
    size_t basis_degree, basis_size;

  public:
    scaled_monomial_scalar_basis(const mesh_type& msh,
                                 const cell_type& cl,
                                 size_t           degree,
                                 bool             use_inertia_axes = false) :
      base(msh, cl, use_inertia_axes)
    {
        basis_degree = degree;
        basis_size   = scalar_basis_size(degree, 3);
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size);

        const auto bp = this->scaling_point(pt);

        size_t pos = 0;
        for (size_t k = 0; k <= basis_degree; k++)
        {
            for (size_t pow_x = 0; pow_x <= k; pow_x++)
            {
                for (size_t pow_y = 0, pow_z = k - pow_x; pow_y <= k - pow_x; pow_y++, pow_z--)
                {
                    const auto px = iexp_pow(bp.x(), pow_x);
                    const auto py = iexp_pow(bp.y(), pow_y);
                    const auto pz = iexp_pow(bp.z(), pow_z);

                    ret(pos++) = px * py * pz;
                }
            }
        }
        assert(pos == basis_size);

        return ret;
    }

    gradient_type
    eval_gradients(const point_type& pt) const
    {
        gradient_type ret = gradient_type::Zero(basis_size, 3);

        const auto scale_fact = this->scaling_factor();
        const auto bp         = this->scaling_point(pt);

        const auto bx = bp.x();
        const auto by = bp.y();
        const auto bz = bp.z();

        const auto                        passage = this->passage_new2old();
        static_vector<coordinate_type, 3> grad_new, grad;

        size_t pos = 0;
        for (size_t k = 0; k <= basis_degree; k++)
        {
            for (size_t pow_x = 0; pow_x <= k; pow_x++)
            {
                for (size_t pow_y = 0, pow_z = k - pow_x; pow_y <= k - pow_x; pow_y++, pow_z--)
                {
                    const auto px = iexp_pow(bx, pow_x);
                    const auto py = iexp_pow(by, pow_y);
                    const auto pz = iexp_pow(bz, pow_z);
                    const auto dx = (pow_x == 0) ? 0 : pow_x * iexp_pow(bx, pow_x - 1);
                    const auto dy = (pow_y == 0) ? 0 : pow_y * iexp_pow(by, pow_y - 1);
                    const auto dz = (pow_z == 0) ? 0 : pow_z * iexp_pow(bz, pow_z - 1);

                    grad_new(0) = dx * py * pz;
                    grad_new(1) = px * dy * pz;
                    grad_new(2) = px * py * dz;
                    grad        = passage * grad_new;

                    ret(pos, 0) = grad(0);
                    ret(pos, 1) = grad(1);
                    ret(pos, 2) = grad(2);

                    pos++;
                }
            }
        }

        assert(pos == basis_size);

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
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename ScalarType>
class scaled_monomial_scalar_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::face, ScalarType>
  : public scaled_monomial_abstract_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::face, ScalarType>
{

  public:
    typedef Mesh<T, 3, Storage>                 mesh_type;
    typedef ScalarType                          scalar_type;
    typedef typename mesh_type::coordinate_type coordinate_type;
    typedef typename mesh_type::point_type      point_type;
    typedef typename mesh_type::face            face_type;
    typedef Matrix<scalar_type, Dynamic, 1>     function_type;

    using base = scaled_monomial_abstract_basis<mesh_type, face_type, scalar_type>;

  private:
    size_t          basis_degree, basis_size;

  public:
    scaled_monomial_scalar_basis(const mesh_type& msh,
                                 const face_type& fc,
                                 size_t           degree,
                                 bool             use_inertia_axes = false) :
      base(msh, fc, use_inertia_axes)
    {
        basis_degree = degree;
        basis_size   = scalar_basis_size(degree, 2);
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        const auto ep = this->scaling_point(pt);

        function_type ret = function_type::Zero(basis_size);
        size_t        pos = 0;
        for (size_t k = 0; k <= basis_degree; k++)
        {
            for (size_t i = 0; i <= k; i++)
            {
                const auto pow_x = k - i;
                const auto pow_y = i;
                const auto px    = iexp_pow(ep.x(), pow_x);
                const auto py    = iexp_pow(ep.y(), pow_y);
                ret(pos++)       = px * py;
            }
        }

        assert(pos == basis_size);

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

/* ====================================================================================================== */
/* ====================================================================================================== */
/*            LEGENDRE BASIS                                                                              */
/* ====================================================================================================== */
/* ====================================================================================================== */

// Compute Legendre basis on [-1, 1] -> orthonormal
template<typename T>
class legendre_1D
{
  public:
    typedef Matrix<T, Dynamic, 1> function_type;
    typedef Matrix<T, Dynamic, 1> gradient_type;

  private:
    size_t basis_degree, basis_size;

    std::array<T, 11>
    eval_monomials(const T& x) const
    {
        std::array<T, 11> pows;
        pows[0] = 1;
        for (size_t i = 1; i <= basis_degree; i++)
            pows[i] = x * pows[i - 1];

        return pows;
    }

    T
    eval_poly(const std::array<T, 11>& pows, const size_t degree) const
    {
        T val;
        switch (degree)
        {
            case 0: val = 1; break;
            case 1: val = pows[1]; break;
            case 2: val = (3 * pows[2] - 1) / 2; break;
            case 3: val = (5 * pows[3] - 3 * pows[1]) / 2; break;
            case 4: val = (35 * pows[4] - 30 * pows[2] + 3) / 8; break;
            case 5: val = (63 * pows[5] - 70 * pows[3] + 15 * pows[1]) / 8; break;
            case 6: val = (231 * pows[6] - 315 * pows[4] + 105 * pows[2] - 5) / 16; break;
            case 7: val = (429 * pows[7] - 693 * pows[5] + 315 * pows[3] - 35 * pows[1]) / 16; break;
            case 8: val = (6435 * pows[8] - 12012 * pows[6] + 6930 * pows[4] - 1260 * pows[2] + 35) / 128; break;
            case 9:
                val = (12155 * pows[9] - 25740 * pows[7] + 18018 * pows[5] - 4620 * pows[3] + 315 * pows[1]) / 128;
                break;
            case 10:
                val =
                  (46189 * pows[10] - 109395 * pows[8] + 90090 * pows[6] - 30030 * pows[4] + 3465 * pows[2] - 63) / 256;
                break;
        }
        return val / sqrt(T(2) / T(2 * degree + 1));
    }

    T
    eval_deriv(const std::array<T, 11>& pows, const size_t degree) const
    {
        T val;
        switch (degree)
        {
            case 0: val = 0; break;
            case 1: val = 1; break;
            case 2: val = (6 * pows[1]) / 2; break;
            case 3: val = (15 * pows[2] - 3) / 2; break;
            case 4: val = (140 * pows[3] - 60 * pows[1]) / 8; break;
            case 5: val = (315 * pows[4] - 210 * pows[2] + 15) / 8; break;
            case 6: val = (1386 * pows[5] - 1260 * pows[3] + 210 * pows[1]) / 16; break;
            case 7: val = (3003 * pows[6] - 3465 * pows[4] + 945 * pows[2] - 35) / 16; break;
            case 8: val = (51480 * pows[7] - 72012 * pows[5] + 27720 * pows[3] - 2520 * pows[1]) / 128; break;
            case 9: val = (109395 * pows[8] - 180180 * pows[6] + 90090 * pows[4] - 13860 * pows[2] + 315) / 128; break;
            case 10:
                val =
                  (461890 * pows[9] - 875160 * pows[7] + 450450 * pows[5] - 120120 * pows[3] + 6930 * pows[1]) / 256;
                break;
        }
        return val / sqrt(T(2) / T(2 * degree + 1));
    }

  public:
    legendre_1D() : basis_degree(10), basis_size(11){};

    legendre_1D(const size_t degree)
    {
        if (degree > 10)
            throw std::invalid_argument("Sorry, I don't have a Legendre basis of order > 10.");

        basis_degree = degree;
        basis_size   = degree + 1;
    }

    function_type
    eval_functions(const T& x) const
    {
        function_type ret = function_type::Zero(basis_size);

        const auto pows = eval_monomials(x);

        for (size_t i = 0; i <= basis_degree; i++)
        {
            ret(i) = eval_poly(pows, i);
        }

        return ret;
    }

    gradient_type
    eval_gradients(const T& x) const
    {
        gradient_type ret = gradient_type::Zero(basis_size);

        const auto pows = eval_monomials(x);

        for (size_t i = 0; i <= basis_degree; i++)
        {
            ret(i) = eval_deriv(pows, i);
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

// Compute Legendre basis on the bounding box of the element

template<typename MeshType, typename Element, typename ScalarType>
struct scaled_legendre_scalar_basis
{
    static_assert(sizeof(MeshType) == -1, "scaled_monomial_scalar_basis: not suitable for the requested kind of mesh");
    static_assert(sizeof(Element) == -1,
                  "scaled_monomial_scalar_basis: not suitable for the requested kind of element");
};

/* Basis 'factory'. */
template<typename MeshType, typename ElementType, typename ScalarType = typename MeshType::coordinate_type>
auto
make_scalar_legendre_basis(const MeshType&    msh,
                           const ElementType& elem,
                           size_t             degree,
                           bool               use_inertia_axes = USE_INERTIA_AXES)
{
    return scaled_legendre_scalar_basis<MeshType, ElementType, ScalarType>(msh, elem, degree, use_inertia_axes);
}

/* Specialization for 2D meshes, cells */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename ScalarType>
class scaled_legendre_scalar_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::cell, ScalarType>
  : public scaled_monomial_abstract_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::cell, ScalarType>
{

  public:
    typedef Mesh<T, 2, Storage>                 mesh_type;
    typedef ScalarType                          scalar_type;
    typedef typename mesh_type::coordinate_type coordinate_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::point_type      point_type;
    typedef Matrix<scalar_type, Dynamic, 2>     gradient_type;
    typedef Matrix<scalar_type, Dynamic, 1>     function_type;

    using base = scaled_monomial_abstract_basis<mesh_type, cell_type, scalar_type>;

  private:
    size_t                   basis_degree, basis_size;
    legendre_1D<scalar_type> leg_1D;

  public:
    scaled_legendre_scalar_basis(const mesh_type& msh,
                                 const cell_type& cl,
                                 size_t           degree,
                                 bool             use_inertia_axes = false) :
      base(msh, cl, use_inertia_axes)
    {
        basis_degree = degree;
        basis_size   = scalar_basis_size(degree, 2);
        leg_1D       = legendre_1D<scalar_type>(degree);
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size);

        const auto bp = this->scaling_point(pt);

        const auto poly_x = leg_1D.eval_functions(bp.x());
        const auto poly_y = leg_1D.eval_functions(bp.y());

        const auto scale_fact = this->scaling_factor();
        const auto scaling    = sqrt(scale_fact(0) * scale_fact(1));

        size_t pos = 0;
        for (size_t k = 0; k <= basis_degree; k++)
        {
            for (size_t i = 0; i <= k; i++)
            {
                const auto pow_x = k - i;
                const auto pow_y = i;

                const auto px = poly_x(pow_x);
                const auto py = poly_y(pow_y);

                ret(pos++) = px * py * scaling;
            }
        }

        assert(pos == basis_size);

        return ret;
    }

    gradient_type
    eval_gradients(const point_type& pt) const
    {
        gradient_type ret = gradient_type::Zero(basis_size, 2);

        const auto scale_fact = this->scaling_factor();
        const auto scaling    = sqrt(scale_fact(0) * scale_fact(1));

        const auto bp = this->scaling_point(pt);

        const auto                        passage = this->passage_new2old();
        static_vector<coordinate_type, 2> grad_new, grad;

        const auto poly_x = leg_1D.eval_functions(bp.x());
        const auto poly_y = leg_1D.eval_functions(bp.y());

        const auto dpoly_x = leg_1D.eval_gradients(bp.x());
        const auto dpoly_y = leg_1D.eval_gradients(bp.y());

        size_t pos = 0;
        for (size_t k = 0; k <= basis_degree; k++)
        {
            for (size_t i = 0; i <= k; i++)
            {
                const auto pow_x = k - i;
                const auto pow_y = i;

                const auto px = poly_x(pow_x);
                const auto py = poly_y(pow_y);
                const auto dx = dpoly_x(pow_x);
                const auto dy = dpoly_y(pow_y);

                grad_new(0) = dx * py;
                grad_new(1) = px * dy;
                grad        = scaling * passage * grad_new;

                ret(pos, 0) = grad(0);
                ret(pos, 1) = grad(1);

                pos++;
            }
        }

        assert(pos == basis_size);

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
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename ScalarType>
class scaled_legendre_scalar_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::face, ScalarType>
  : public scaled_monomial_abstract_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::face, ScalarType>
{

  public:
    typedef Mesh<T, 2, Storage>             mesh_type;
    typedef ScalarType                      scalar_type;
    typedef typename mesh_type::point_type  point_type;
    typedef typename mesh_type::face        face_type;
    typedef Matrix<scalar_type, Dynamic, 1> function_type;

    using base = scaled_monomial_abstract_basis<mesh_type, face_type, scalar_type>;

  private:
    size_t                   basis_degree, basis_size;
    legendre_1D<scalar_type> leg_x;

  public:
    scaled_legendre_scalar_basis(const mesh_type& msh,
                                 const face_type& fc,
                                 size_t           degree,
                                 bool             use_inertia_axes = false) :
      base(msh, fc, use_inertia_axes)
    {
        basis_degree = degree;
        basis_size   = degree + 1;
        leg_x        = legendre_1D<scalar_type>(degree);
    }

    bool
    is_orthonormal() const
    {
        return true;
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        const auto        ep      = this->scaling_point(pt);
        const scalar_type scaling = this->scaling_factor();

        return leg_x.eval_functions(ep.x()) * scaling;
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

/***************************************************************************************************/
/***************************************************************************************************/

/* Specialization for 3D meshes, cells */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename ScalarType>
class scaled_legendre_scalar_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::cell, ScalarType>
  : public scaled_monomial_abstract_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::cell, ScalarType>

{

  public:
    typedef Mesh<T, 3, Storage>                 mesh_type;
    typedef ScalarType                          scalar_type;
    typedef typename mesh_type::coordinate_type coordinate_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::point_type      point_type;
    typedef Matrix<scalar_type, Dynamic, 3>     gradient_type;
    typedef Matrix<scalar_type, Dynamic, 1>     function_type;

    using base = scaled_monomial_abstract_basis<mesh_type, cell_type, scalar_type>;

  private:
    size_t                   basis_degree, basis_size;
    legendre_1D<scalar_type> leg_1D;

  public:
    scaled_legendre_scalar_basis(const mesh_type& msh,
                                 const cell_type& cl,
                                 size_t           degree,
                                 bool             use_inertia_axes = false) :
      base(msh, cl, use_inertia_axes)
    {
        basis_degree = degree;
        basis_size   = scalar_basis_size(degree, 3);
        leg_1D       = legendre_1D<scalar_type>(degree);
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size);

        const auto scale_fact = this->scaling_factor();
        const auto bp         = this->scaling_point(pt);

        const auto poly_x = leg_1D.eval_functions(bp.x());
        const auto poly_y = leg_1D.eval_functions(bp.y());
        const auto poly_z = leg_1D.eval_functions(bp.z());

        const auto scaling = sqrt(scale_fact(0) * scale_fact(1) * scale_fact(2));

        size_t pos = 0;
        for (size_t k = 0; k <= basis_degree; k++)
        {
            for (size_t pow_x = 0; pow_x <= k; pow_x++)
            {
                for (size_t pow_y = 0, pow_z = k - pow_x; pow_y <= k - pow_x; pow_y++, pow_z--)
                {
                    const auto px = poly_x(pow_x);
                    const auto py = poly_y(pow_y);
                    const auto pz = poly_z(pow_z);

                    ret(pos++) = px * py * pz * scaling;
                }
            }
        }
        assert(pos == basis_size);

        return ret;
    }

    gradient_type
    eval_gradients(const point_type& pt) const
    {
        gradient_type ret = gradient_type::Zero(basis_size, 3);

        const auto scale_fact = this->scaling_factor();
        const auto bp         = this->scaling_point(pt);

        const auto poly_x = leg_1D.eval_functions(bp.x());
        const auto poly_y = leg_1D.eval_functions(bp.y());
        const auto poly_z = leg_1D.eval_functions(bp.z());

        const auto dpoly_x = leg_1D.eval_gradients(bp.x());
        const auto dpoly_y = leg_1D.eval_gradients(bp.y());
        const auto dpoly_z = leg_1D.eval_gradients(bp.z());

        const auto                        passage = this->passage_new2old();
        static_vector<coordinate_type, 3> grad_new, grad;

        const auto scaling = sqrt(scale_fact(0) * scale_fact(1) * scale_fact(2));

        size_t pos = 0;
        for (size_t k = 0; k <= basis_degree; k++)
        {
            for (size_t pow_x = 0; pow_x <= k; pow_x++)
            {
                for (size_t pow_y = 0, pow_z = k - pow_x; pow_y <= k - pow_x; pow_y++, pow_z--)
                {
                    const auto px = poly_x(pow_x);
                    const auto py = poly_y(pow_y);
                    const auto pz = poly_z(pow_z);
                    const auto dx = dpoly_x(pow_x);
                    const auto dy = dpoly_y(pow_y);
                    const auto dz = dpoly_z(pow_z);

                    grad_new(0) = dx * py * pz;
                    grad_new(1) = px * dy * pz;
                    grad_new(2) = px * py * dz;
                    grad        = scaling * passage * grad_new;

                    ret(pos, 0) = grad(0);
                    ret(pos, 1) = grad(1);
                    ret(pos, 2) = grad(2);

                    pos++;
                }
            }
        }

        assert(pos == basis_size);

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
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename ScalarType>
class scaled_legendre_scalar_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::face, ScalarType>
  : public scaled_monomial_abstract_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::face, ScalarType>
{

  public:
    typedef Mesh<T, 3, Storage>                 mesh_type;
    typedef ScalarType                          scalar_type;
    typedef typename mesh_type::coordinate_type coordinate_type;
    typedef typename mesh_type::point_type      point_type;
    typedef typename mesh_type::face            face_type;
    typedef Matrix<scalar_type, Dynamic, 1>     function_type;

    using base = scaled_monomial_abstract_basis<mesh_type, face_type, scalar_type>;

  private:
    size_t                   basis_degree, basis_size;
    legendre_1D<scalar_type> leg_1D;

  public:
    scaled_legendre_scalar_basis(const mesh_type& msh,
                                 const face_type& fc,
                                 size_t           degree,
                                 bool             use_inertia_axes = false) :
      base(msh, fc, use_inertia_axes)
    {
        basis_degree = degree;
        basis_size   = scalar_basis_size(degree, 2);
        leg_1D       = legendre_1D<scalar_type>(degree);
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        const auto ep = this->scaling_point(pt);

        const auto poly_x = leg_1D.eval_functions(ep.x());
        const auto poly_y = leg_1D.eval_functions(ep.y());

        const auto scaling_fact = this->scaling_factor();

        const auto scaling = sqrt(scaling_fact(0) * scaling_fact(1));

        function_type ret = function_type::Zero(basis_size);
        size_t        pos = 0;
        for (size_t k = 0; k <= basis_degree; k++)
        {
            for (size_t i = 0; i <= k; i++)
            {
                const auto pow_x = k - i;
                const auto pow_y = i;
                const auto px    = poly_x(pow_x);
                const auto py    = poly_y(pow_y);
                ret(pos++)       = scaling * px * py;
            }
        }

        assert(pos == basis_size);

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

#ifdef USE_LEGENDRE
template<typename MeshType, typename ElementType, typename ScalarType = typename MeshType::coordinate_type>
auto
make_scalar_monomial_basis(const MeshType&    msh,
                           const ElementType& elem,
                           size_t             degree,
                           bool               use_inertia_axes = USE_INERTIA_AXES)
{
    return make_scalar_legendre_basis<MeshType, ElementType, ScalarType>(msh, elem, degree, use_inertia_axes);
}
#endif

/*
template<typename Mesh, typename Element, typename Basis>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
compute_averages(const Mesh& msh, const Element& elem, const Basis& basis)
{
    using T = typename Mesh::coordinate_type;
    auto meas = measure(msh, elem);

    Matrix<T, Dynamic, 1> avgs = Matrix<T, Dynamic, 1>::Zero(basis.size());

    auto qps = integrate(msh, elem, basis.degree());
    for (auto& qp : qps)
        avgs += qp.weight() * basis.eval_functions(qp.point());

    avgs(0) = 0;

    return avgs * (1/meas);
}
*/

//#define USE_H2DT

inline size_t
harmonic_basis_size(size_t k, size_t d)
{
    if (d == 2)
        return 2*k+1;

    if (d == 3) {
#ifdef USE_H2DT
        return (k == 0) ? 1 : 4+6*(k-1);
#else
        return (k+1)*(k+1);
#endif
    }

    throw std::invalid_argument("not yet implemented");
}

#if 0
/* Generic template for bases. */
template<typename MeshType, typename Element, typename ScalarType>
struct scaled_harmonic_scalar_basis
{
    static_assert(sizeof(MeshType) == -1, "scaled_harmonic_scalar_basis: not suitable for the requested kind of mesh");
    static_assert(sizeof(Element) == -1,
                  "scaled_harmonic_scalar_basis: not suitable for the requested kind of element");
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename ScalarType>
class scaled_harmonic_scalar_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::cell, ScalarType>
    : public scaled_monomial_abstract_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::cell, ScalarType>
{

public:
    typedef Mesh<T, 2, Storage>                 mesh_type;
    typedef ScalarType                          scalar_type;
    typedef typename mesh_type::coordinate_type coordinate_type;

    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::point_type      point_type;
    typedef Matrix<scalar_type, Dynamic, 2>     gradient_type;
    typedef Matrix<scalar_type, Dynamic, 1>     function_type;

    using base = scaled_monomial_abstract_basis<mesh_type, cell_type, scalar_type>;

private:
    size_t basis_degree, basis_size;

public:
    scaled_harmonic_scalar_basis(const mesh_type& msh, const cell_type& cl, size_t degree)
        : base(msh, cl, false)
    {
        basis_degree = degree;
        basis_size   = 1+2*basis_degree;
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        const auto bp = this->scaling_point(pt);
        function_type ret = function_type::Zero(basis_size);
        ret(0) = 1.0;

        if ( basis_degree > 0 ) {
            ret(1) = bp.x();
            ret(2) = bp.y();
        }

        for (size_t i = 2; i <= basis_degree; i++) {
            auto Pi = ret(2*i-3);
            auto Qi = ret(2*i-2);
            ret(2*i-1) = bp.x()*Pi - bp.y()*Qi;
            ret(2*i) = bp.y()*Pi + bp.x()*Qi;
        }

        return ret;
    }

    gradient_type
    eval_gradients(const point_type& pt) const
    {
        const auto bp = this->scaling_point(pt);
        gradient_type ret = gradient_type::Zero(basis_size, 2);
        ret(0,0) = 0.0; ret(0,1) = 0.0;

        if ( basis_degree > 0 ) {
            ret(1,0) = 1.0; ret(1,1) = 0.0;
            ret(2,0) = 0.0; ret(2,1) = 1.0;

            if (basis_degree == 1)
                return ret;
        }

        auto funcs = eval_functions(pt);

        for (size_t i = 2; i <= basis_degree; i++) {
            auto iPi = i*funcs(2*i-3);
            auto iQi = i*funcs(2*i-2);
            ret(2*i-1, 0) = iPi; ret(2*i-1, 1) = -iQi;
            ret(2*i, 0) = iQi; ret(2*i, 1) = iPi;
        }

        return ret;
    }

    size_t size() const {
        return basis_size;
    }
};


template<typename Mesh, typename ElemType, typename ScalarType = typename Mesh::coordinate_type>
auto
make_scalar_harmonic_basis(const Mesh& msh, const ElemType& elem,size_t degree)
{
    return scaled_harmonic_scalar_basis<Mesh, ElemType, ScalarType>(msh, elem, degree);
}
#endif





template<typename MeshType, typename Element, typename ScalarType>
struct scaled_harmonic_top_scalar_basis
{
    static_assert(sizeof(MeshType) == -1, "scaled_harmonic_top_scalar_basis: not suitable for the requested kind of mesh");
    static_assert(sizeof(Element) == -1,
                  "scaled_harmonic_top_scalar_basis: not suitable for the requested kind of element");
};


template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename ScalarType>
class scaled_harmonic_top_scalar_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::cell, ScalarType>
    : public scaled_monomial_abstract_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::cell, ScalarType>
{

public:
    typedef Mesh<T, 2, Storage>                 mesh_type;
    typedef ScalarType                          scalar_type;
    typedef typename mesh_type::coordinate_type coordinate_type;

    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::point_type      point_type;
    typedef Matrix<scalar_type, Dynamic, 2>     gradient_type;
    typedef Matrix<scalar_type, Dynamic, 1>     function_type;

    using base = scaled_monomial_abstract_basis<mesh_type, cell_type, scalar_type>;

private:
    size_t full_basis_degree, full_basis_size;
    size_t polynomial_max_degree, polynomial_part_size;
    size_t harmonic_min_degree, harmonic_part_size;
    point_type bar;
    T diam;

    function_type
    eval_harmonic(const point_type& bp) const
    {
        auto hbs = harmonic_basis_size(full_basis_degree, 2);
        function_type ret = function_type::Zero(hbs);
        ret(0) = 1.0;

        if ( full_basis_degree > 0 ) {
            ret(1) = bp.x();
            ret(2) = bp.y();
        }

        for (size_t i = 2; i <= full_basis_degree; i++) {
            auto Pi = ret(2*i-3);
            auto Qi = ret(2*i-2);
            ret(2*i-1) = bp.x()*Pi - bp.y()*Qi;
            ret(2*i) = bp.y()*Pi + bp.x()*Qi;
        }

        return ret;
    }

    gradient_type
    eval_harmonic_gradients(const point_type& bp) const
    {
        auto hbs = harmonic_basis_size(full_basis_degree, 2);
        gradient_type ret = gradient_type::Zero(hbs, 2);
        ret(0,0) = 0.0; ret(0,1) = 0.0;

        //const auto p = this->passage_new2old();
        auto sx = 1.0;//p(0,0);
        auto sy = 1.0;//p(1,1);

        if ( full_basis_degree > 0 ) {
            ret(1,0) = sx; ret(1,1) = 0.0;
            ret(2,0) = 0.0; ret(2,1) = sy;

            if (full_basis_degree == 1)
                return ret;
        }

        auto funcs = eval_harmonic(bp);

        for (size_t i = 2; i <= full_basis_degree; i++) {
            auto iPi = i*funcs(2*i-3);
            auto iQi = i*funcs(2*i-2);
            ret(2*i-1, 0) = sx*iPi; ret(2*i-1, 1) = -sy*iQi;
            ret(2*i, 0) = sx*iQi; ret(2*i, 1) = sy*iPi;
        }

        return ret;
    }

public:
    scaled_harmonic_top_scalar_basis(const mesh_type& msh, const cell_type& cl, size_t fbd)
        : base(msh, cl, false)
    {
        full_basis_degree       = fbd;
        full_basis_size         = scalar_basis_size(full_basis_degree, 2);
        polynomial_max_degree   = full_basis_degree;
        polynomial_part_size    = full_basis_size;
        harmonic_min_degree     = full_basis_degree+1;
        harmonic_part_size      = 0;
        bar                     = barycenter(msh, cl);
        diam                    = diameter(msh, cl);
    }

    void
    maximum_polynomial_degree(size_t pd)
    {
        if (pd > full_basis_degree)
            return;
        
        polynomial_max_degree = pd;
        polynomial_part_size = scalar_basis_size(polynomial_max_degree, 2);
        harmonic_min_degree = polynomial_max_degree+1;
        assert(harmonic_min_degree > 0);
        size_t hs_notused = harmonic_basis_size(harmonic_min_degree-1, 2);
        size_t hs_used = harmonic_basis_size(full_basis_degree, 2);
        assert(hs_notused <= hs_used);
        harmonic_part_size = hs_used - hs_notused;
        full_basis_size = polynomial_part_size + harmonic_part_size;
    }

    size_t
    maximum_polynomial_degree(void) const
    {
        return polynomial_max_degree;
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        const auto bp = this->scaling_point(pt);
        function_type ret = function_type::Zero(full_basis_size);

        size_t pos = 0;
        for (size_t k = 0; k <= polynomial_max_degree; k++) {
            for (size_t i = 0; i <= k; i++) {
                const auto pow_x = k - i;
                const auto pow_y = i;
                const auto px = iexp_pow(bp.x(), pow_x);
                const auto py = iexp_pow(bp.y(), pow_y);
                ret(pos++) = px * py;
            }
        }
        assert(pos == polynomial_part_size);
        if (harmonic_part_size > 0) {
            assert(polynomial_max_degree < full_basis_degree);
            auto hep = 2.0*(pt-bar)/diam; /* beware of scaling: must be the same in all dirs */
            function_type hp = eval_harmonic(hep);
            ret.tail(harmonic_part_size) = hp.tail(harmonic_part_size);
        }
        return ret;
    }

    gradient_type
    eval_gradients(const point_type& pt) const
    {
        gradient_type ret = gradient_type::Zero(full_basis_size, 2);
        const auto bp = this->scaling_point(pt);
        const auto bx = bp.x();
        const auto by = bp.y();

        const auto p = this->passage_new2old();
        size_t pos = 0;
        for (size_t k = 0; k <= polynomial_max_degree; k++) {
            for (size_t i = 0; i <= k; i++) {
                const auto pow_x = k - i;
                const auto pow_y = i;
                const auto px = iexp_pow(bx, pow_x);
                const auto py = iexp_pow(by, pow_y);
                const auto dx = (pow_x == 0) ? 0 : pow_x * iexp_pow(bx, pow_x - 1);
                const auto dy = (pow_y == 0) ? 0 : pow_y * iexp_pow(by, pow_y - 1);
                ret(pos, 0) = p(0,0) * dx * py;
                ret(pos, 1) = p(1,1) * px * dy;
                pos++;
            }
        }
        assert(pos == polynomial_part_size);
        if (harmonic_part_size > 0) {
            assert(polynomial_max_degree < full_basis_degree);
            auto hep = 2.0*(pt-bar)/diam; /* beware of scaling: must be the same in all dirs */
            gradient_type hg = eval_harmonic_gradients(hep);
            ret.block(ret.rows()-harmonic_part_size, 0, harmonic_part_size, 2) =
                hg.block(hg.rows()-harmonic_part_size, 0, harmonic_part_size, 2);
        }
        return ret;
    }

    size_t size() const {
        return full_basis_size;
    }

    size_t poly_size() const {
        return polynomial_part_size;
    }

    size_t harmonic_size() const {
        return harmonic_part_size;
    }

    size_t degree() const {
        return full_basis_degree;
    }
};

#ifdef USE_H2DT

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename ScalarType>
class scaled_harmonic_top_scalar_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::cell, ScalarType>
    : public scaled_monomial_abstract_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::cell, ScalarType>
{

public:
    typedef Mesh<T, 3, Storage>                 mesh_type;
    typedef ScalarType                          scalar_type;
    typedef typename mesh_type::coordinate_type coordinate_type;

    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::point_type      point_type;
    typedef Matrix<scalar_type, Dynamic, 3>     gradient_type;
    typedef Matrix<scalar_type, Dynamic, 1>     function_type;

    using base = scaled_monomial_abstract_basis<mesh_type, cell_type, scalar_type>;

private:
    size_t full_basis_degree, full_basis_size;
    size_t polynomial_max_degree, polynomial_part_size;
    size_t harmonic_min_degree, harmonic_part_size;
    point_type bar;
    T diam;

    function_type
    eval_harmonic(const point_type& bp) const
    {
        auto hbs = harmonic_basis_size(full_basis_degree, 3);
        function_type ret = function_type::Zero(hbs);
        ret(0) = 1.0;

        if ( full_basis_degree > 0 ) {
            ret(1) = bp.x();
            ret(2) = bp.y();
            ret(3) = bp.z();
        }

        if ( full_basis_degree > 1){//} && harmonic_min_degree < 2 ) {
            auto x = bp.x(); auto y = bp.y(); auto z = bp.z();
            ret(4) = x*x - y*y;
            ret(5) = 2*x*y;
            ret(6) = x*x - z*z;
            ret(7) = 2*x*z;
            ret(8) = y*y - z*z;
            ret(9) = 2*y*z;
        }

        for (size_t i = 3; i <= full_basis_degree; i++) {
            auto x = bp.x(); auto y = bp.y(); auto z = bp.z();
            auto iPi = ret(6*i-14);
            auto iQi = ret(6*i-13);
            auto jPi = ret(6*i-12);
            auto jQi = ret(6*i-11);
            auto kPi = ret(6*i-10);
            auto kQi = ret(6*i-9);
            ret(6*i-8) = x*iPi - y*iQi;
            ret(6*i-7) = y*iPi + x*iQi;
            ret(6*i-6) = x*jPi - z*jQi;
            ret(6*i-5) = z*jPi + x*jQi;
            ret(6*i-4) = y*kPi - z*kQi;
            ret(6*i-3) = z*kPi + y*kQi;
        }

        return ret;
    }

    gradient_type
    eval_harmonic_gradients(const point_type& bp) const
    {
        auto hbs = harmonic_basis_size(full_basis_degree, 3);
        gradient_type ret = gradient_type::Zero(hbs, 3);
        ret(0,0) = 0.0; ret(0,1) = 0.0; ret(0,2) = 0.0;

        auto sx = 1.0;
        auto sy = 1.0;
        auto sz = 1.0;

        if ( full_basis_degree > 0 ) {
            ret(1,0) = sx; ret(1,1) = 0.0; ret(1,2) = 0.0;
            ret(2,0) = 0.0; ret(2,1) = sy; ret(2,2) = 0.0;
            ret(3,0) = 0.0; ret(3,1) = 0.0; ret(3,2) = sz;
        }

        if (full_basis_degree > 1){//} && harmonic_min_degree < 2) {
            ret(4,0) = sx*2.0*bp.x(); ret(4,1) = -sy*2.0*bp.y(); ret(4,2) = 0.0;
            ret(5,0) = sx*2.0*bp.y(); ret(5,1) =  sy*2.0*bp.x(); ret(5,2) = 0.0;
            ret(6,0) = sx*2.0*bp.x(); ret(6,1) = 0.0; ret(6,2) = -sz*2.0*bp.z();
            ret(7,0) = sx*2.0*bp.z(); ret(7,1) = 0.0; ret(7,2) =  sz*2.0*bp.x();
            ret(8,0) = 0.0; ret(8,1) = sy*2.0*bp.y(); ret(8,2) = -sz*2.0*bp.z();
            ret(9,0) = 0.0; ret(9,1) = sy*2.0*bp.z(); ret(9,2) =  sz*2.0*bp.y();
        }

        auto funcs = eval_harmonic(bp);

        for (size_t i = 3; i <= full_basis_degree; i++) {
            auto i_iPi = i*funcs(6*i-14);
            auto i_iQi = i*funcs(6*i-13);
            auto j_iPi = i*funcs(6*i-12);
            auto j_iQi = i*funcs(6*i-11);
            auto k_iPi = i*funcs(6*i-10);
            auto k_iQi = i*funcs(6*i-9);
            ret(6*i-8, 0) = sx*i_iPi; ret(6*i-8, 1) = -sy*i_iQi; ret(6*i-8, 2) = 0.0;
            ret(6*i-7, 0) = sx*i_iQi; ret(6*i-7, 1) =  sy*i_iPi; ret(6*i-7, 2) = 0.0;
            ret(6*i-6, 0) = sx*j_iPi; ret(6*i-6, 1) = 0.0; ret(6*i-6, 2) = -sz*j_iQi; 
            ret(6*i-5, 0) = sx*j_iQi; ret(6*i-5, 1) = 0.0; ret(6*i-5, 2) =  sz*j_iPi;
            ret(6*i-4, 0) = 0.0; ret(6*i-4, 1) = sy*k_iPi; ret(6*i-4, 2) = -sz*k_iQi; 
            ret(6*i-3, 0) = 0.0; ret(6*i-3, 1) = sy*k_iQi; ret(6*i-3, 2) =  sz*k_iPi;  
        }

        return ret;
    }

public:
    scaled_harmonic_top_scalar_basis(const mesh_type& msh, const cell_type& cl, size_t fbd)
        : base(msh, cl, false)
    {
        full_basis_degree       = fbd;
        full_basis_size         = scalar_basis_size(full_basis_degree, 3);
        polynomial_max_degree   = full_basis_degree;
        polynomial_part_size    = full_basis_size;
        harmonic_min_degree     = full_basis_degree+1;
        harmonic_part_size      = 0;
        bar                     = barycenter(msh, cl);
        diam                    = diameter(msh, cl);
    }

    void
    maximum_polynomial_degree(size_t pd)
    {
        if (pd > full_basis_degree)
            return;
        
        polynomial_max_degree = pd;
        polynomial_part_size = scalar_basis_size(polynomial_max_degree, 3);
        harmonic_min_degree = polynomial_max_degree+1;
        assert(harmonic_min_degree > 0);
        size_t hs_notused = harmonic_basis_size(harmonic_min_degree-1, 3);
        size_t hs_used = harmonic_basis_size(full_basis_degree, 3);
        assert(hs_notused <= hs_used);
        harmonic_part_size = hs_used - hs_notused;
        full_basis_size = polynomial_part_size + harmonic_part_size;

        std::cout << "HMD: " << harmonic_min_degree << std::endl;
    }

    size_t
    maximum_polynomial_degree(void) const
    {
        return polynomial_max_degree;
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        const auto bp = this->scaling_point(pt);
        function_type ret = function_type::Zero(full_basis_size);

        size_t pos = 0;
        for (size_t k = 0; k <= polynomial_max_degree; k++) {
            for (size_t pow_x = 0; pow_x <= k; pow_x++) {
                for (size_t pow_y = 0, pow_z = k - pow_x; pow_y <= k - pow_x; pow_y++, pow_z--) {
                    const auto px = iexp_pow(bp.x(), pow_x);
                    const auto py = iexp_pow(bp.y(), pow_y);
                    const auto pz = iexp_pow(bp.z(), pow_z);
                    ret(pos++) = px * py * pz;
                }
            }
        }
        assert(pos == polynomial_part_size);

        if (harmonic_part_size > 0) {
            assert(polynomial_max_degree < full_basis_degree);
            auto hep = 2.0*(pt-bar)/diam; /* beware of scaling: must be the same in all dirs */
            function_type hp = eval_harmonic(hep);
            ret.tail(harmonic_part_size) = hp.tail(harmonic_part_size);
        }
        return ret;
    }

    gradient_type
    eval_gradients(const point_type& pt) const
    {
        gradient_type ret = gradient_type::Zero(full_basis_size, 3);
        const auto bp = this->scaling_point(pt);

        const auto bx = bp.x();
        const auto by = bp.y();
        const auto bz = bp.z();

        const auto                        passage = this->passage_new2old();
        static_vector<coordinate_type, 3> grad_new, grad;

        size_t pos = 0;
        for (size_t k = 0; k <= polynomial_max_degree; k++) {
            for (size_t pow_x = 0; pow_x <= k; pow_x++) {
                for (size_t pow_y = 0, pow_z = k - pow_x; pow_y <= k - pow_x; pow_y++, pow_z--) {
                    const auto px = iexp_pow(bx, pow_x);
                    const auto py = iexp_pow(by, pow_y);
                    const auto pz = iexp_pow(bz, pow_z);
                    const auto dx = (pow_x == 0) ? 0 : pow_x * iexp_pow(bx, pow_x - 1);
                    const auto dy = (pow_y == 0) ? 0 : pow_y * iexp_pow(by, pow_y - 1);
                    const auto dz = (pow_z == 0) ? 0 : pow_z * iexp_pow(bz, pow_z - 1);

                    grad_new(0) = dx * py * pz;
                    grad_new(1) = px * dy * pz;
                    grad_new(2) = px * py * dz;
                    grad        = passage * grad_new;

                    ret(pos, 0) = grad(0);
                    ret(pos, 1) = grad(1);
                    ret(pos, 2) = grad(2);

                    pos++;
                }
            }
        }

        assert(pos == polynomial_part_size);

        if (harmonic_part_size > 0) {
            assert(polynomial_max_degree < full_basis_degree);
            auto hep = 2.0*(pt-bar)/diam; /* beware of scaling: must be the same in all dirs */
            gradient_type hg = eval_harmonic_gradients(hep);
            ret.block(ret.rows()-harmonic_part_size, 0, harmonic_part_size, 3) =
                hg.block(hg.rows()-harmonic_part_size, 0, harmonic_part_size, 3);
        }
        return ret;
    }

    size_t size() const {
        return full_basis_size;
    }

    size_t poly_size() const {
        return polynomial_part_size;
    }

    size_t harmonic_size() const {
        return harmonic_part_size;
    }

    size_t degree() const {
        return full_basis_degree;
    }
};

#else /* USE_H2DT */

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename ScalarType>
class scaled_harmonic_top_scalar_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::cell, ScalarType>
    : public scaled_monomial_abstract_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::cell, ScalarType>
{

public:
    typedef Mesh<T, 3, Storage>                 mesh_type;
    typedef ScalarType                          scalar_type;
    typedef typename mesh_type::coordinate_type coordinate_type;

    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::point_type      point_type;
    typedef Matrix<scalar_type, Dynamic, 3>     gradient_type;
    typedef Matrix<scalar_type, Dynamic, 1>     function_type;
    typedef Matrix<scalar_type, Dynamic, Dynamic>     cs_matrix_type;

    using base = scaled_monomial_abstract_basis<mesh_type, cell_type, scalar_type>;

private:
    size_t full_basis_degree, full_basis_size;
    size_t polynomial_max_degree, polynomial_part_size;
    size_t harmonic_min_degree, harmonic_part_size;
    point_type bar;
    T diam;

    auto eval_cs(const point_type& bp) const
    {
        auto hmd = full_basis_degree;
        cs_matrix_type Cs = cs_matrix_type::Zero(hmd+1, hmd+1);
        cs_matrix_type Ss = cs_matrix_type::Zero(hmd+1, hmd+1);
        Cs(0,0) = 1.0;

        if (hmd > 0) {
            Cs(1,0) = bp.x();
            Cs(1,1) = bp.y();
            Ss(1,1) = bp.z();
        }

        auto R2 = bp.x()*bp.x() + bp.y()*bp.y() + bp.z()*bp.z();
        auto r2 = bp.y()*bp.y() + bp.z()*bp.z();

        for(int n = 0; n <= hmd; n++) {
            for (int m = 0; m <= n; m++) {
                if (n >= 2 && m < n) {
                    auto multn = 2*(m+1);
                    auto multd = (n+m+1)*r2;
                    auto mult = multn/multd;
                    auto c1 = (R2*Cs(n-1,m)-bp.x()*Cs(n,m))*bp.y();
                    auto c2 = (R2*Ss(n-1,m)-bp.x()*Ss(n,m))*bp.z();
                    Cs(n,m+1) = mult*(c1-c2);
                    auto s1 = (R2*Ss(n-1,m)-bp.x()*Ss(n,m))*bp.y();
                    auto s2 = (R2*Cs(n-1,m)-bp.x()*Cs(n,m))*bp.z();
                    Ss(n,m+1) = mult*(s1+s2);
                }
                if (n >= 1 and m < 2 and n < hmd) {
                    auto mult = 1.0/(n+m+1.0);
                    Cs(n+1,m) = mult*( (2*n+1)*bp.x()*Cs(n,m) - (n-m)*R2*Cs(n-1, m) );
                    Ss(n+1,m) = mult*( (2*n+1)*bp.x()*Ss(n,m) - (n-m)*R2*Ss(n-1, m) );
                }
            }
        }

        //std::cout << " *** Cs ***" << std::endl;
        //std::cout << Cs << std::endl;
        //std::cout << " *** Ss ***" << std::endl;
        //std::cout << Ss << std::endl;

        return std::pair(Cs, Ss);
    }

    function_type
    eval_harmonic(const point_type& bp) const
    {
        auto hmd = full_basis_degree;
        auto bsize = harmonic_basis_size(full_basis_degree, 3);
        auto [Cs, Ss] = eval_cs(bp);
        function_type ret = function_type::Zero(bsize);
        assert(bsize>0);
        ret(0) = 1.0;
        size_t ofs = 1;
        for(size_t n = 1; n <= hmd; n++) {
            for (size_t m = 0; m <= n; m++) {
                ret(ofs++) = Cs(n,m);
            }
            
            for (size_t m = 1; m <= n; m++) {
                ret(ofs++) = Ss(n,m);
            }
        }
        assert(ofs == bsize);
        return ret;
    }

    gradient_type
    eval_harmonic_gradients(const point_type& bp) const
    {
        auto hmd = full_basis_degree;
        auto bsize = harmonic_basis_size(full_basis_degree, 3);
        auto [Cs, Ss] = eval_cs(bp);
        gradient_type ret = gradient_type::Zero(bsize, 3);
        
        auto r2 = bp.y()*bp.y() + bp.z()*bp.z();

        size_t ofs = 1;
        for(size_t n = 1; n <= hmd; n++) {
            for (size_t m = 0; m <= n; m++) {
                ret(ofs, 0) = (n-m)*Cs(n-1,m);
                ret(ofs, 1) = (bp.y()/r2)*(n*Cs(n,m) - (n-m)*bp.x()*Cs(n-1,m)) + m*Ss(n,m)*bp.z()/r2;
                ret(ofs, 2) = (bp.z()/r2)*(n*Cs(n,m) - (n-m)*bp.x()*Cs(n-1,m)) - m*Ss(n,m)*bp.y()/r2;
                ofs++;
            }

            for (size_t m = 1; m <= n; m++) {
                ret(ofs, 0) = (n-m)*Ss(n-1,m);
                ret(ofs, 1) = (bp.y()/r2)*(n*Ss(n,m) - (n-m)*bp.x()*Ss(n-1,m)) - m*Cs(n,m)*bp.z()/r2;
                ret(ofs, 2) = (bp.z()/r2)*(n*Ss(n,m) - (n-m)*bp.x()*Ss(n-1,m)) + m*Cs(n,m)*bp.y()/r2;
                ofs++;
            }
        }
        assert(ofs == bsize);

        return ret;
    }

public:
    scaled_harmonic_top_scalar_basis(const mesh_type& msh, const cell_type& cl, size_t fbd)
        : base(msh, cl, false)
    {
        full_basis_degree       = fbd;
        full_basis_size         = scalar_basis_size(full_basis_degree, 3);
        polynomial_max_degree   = full_basis_degree;
        polynomial_part_size    = full_basis_size;
        harmonic_min_degree     = full_basis_degree+1;
        harmonic_part_size      = 0;
        bar                     = barycenter(msh, cl);
        diam                    = diameter(msh, cl);
    }

    void
    maximum_polynomial_degree(size_t pd)
    {
        if (pd > full_basis_degree)
            return;
        
        polynomial_max_degree = pd;
        polynomial_part_size = scalar_basis_size(polynomial_max_degree, 3);
        harmonic_min_degree = polynomial_max_degree+1;
        assert(harmonic_min_degree > 0);
        size_t hs_notused = harmonic_basis_size(harmonic_min_degree-1, 3);
        size_t hs_used = harmonic_basis_size(full_basis_degree, 3);
        assert(hs_notused <= hs_used);
        harmonic_part_size = hs_used - hs_notused;
        full_basis_size = polynomial_part_size + harmonic_part_size;
    }

    size_t
    maximum_polynomial_degree(void) const
    {
        return polynomial_max_degree;
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        const auto bp = this->scaling_point(pt);
        function_type ret = function_type::Zero(full_basis_size);

        size_t pos = 0;
        for (size_t k = 0; k <= polynomial_max_degree; k++) {
            for (size_t pow_x = 0; pow_x <= k; pow_x++) {
                for (size_t pow_y = 0, pow_z = k - pow_x; pow_y <= k - pow_x; pow_y++, pow_z--) {
                    const auto px = iexp_pow(bp.x(), pow_x);
                    const auto py = iexp_pow(bp.y(), pow_y);
                    const auto pz = iexp_pow(bp.z(), pow_z);
                    ret(pos++) = px * py * pz;
                }
            }
        }
        assert(pos == polynomial_part_size);

        if (harmonic_part_size > 0) {
            assert(polynomial_max_degree < full_basis_degree);
            auto hep = 2.0*(pt-bar)/diam; /* beware of scaling: must be the same in all dirs */
            function_type hp = eval_harmonic(hep);
            ret.tail(harmonic_part_size) = hp.tail(harmonic_part_size);
        }
        return ret;
    }

    gradient_type
    eval_gradients(const point_type& pt) const
    {
        gradient_type ret = gradient_type::Zero(full_basis_size, 3);
        const auto bp = this->scaling_point(pt);

        const auto bx = bp.x();
        const auto by = bp.y();
        const auto bz = bp.z();

        const auto                        passage = this->passage_new2old();
        static_vector<coordinate_type, 3> grad_new, grad;

        size_t pos = 0;
        for (size_t k = 0; k <= polynomial_max_degree; k++) {
            for (size_t pow_x = 0; pow_x <= k; pow_x++) {
                for (size_t pow_y = 0, pow_z = k - pow_x; pow_y <= k - pow_x; pow_y++, pow_z--) {
                    const auto px = iexp_pow(bx, pow_x);
                    const auto py = iexp_pow(by, pow_y);
                    const auto pz = iexp_pow(bz, pow_z);
                    const auto dx = (pow_x == 0) ? 0 : pow_x * iexp_pow(bx, pow_x - 1);
                    const auto dy = (pow_y == 0) ? 0 : pow_y * iexp_pow(by, pow_y - 1);
                    const auto dz = (pow_z == 0) ? 0 : pow_z * iexp_pow(bz, pow_z - 1);

                    grad_new(0) = dx * py * pz;
                    grad_new(1) = px * dy * pz;
                    grad_new(2) = px * py * dz;
                    grad        = passage * grad_new;

                    ret(pos, 0) = grad(0);
                    ret(pos, 1) = grad(1);
                    ret(pos, 2) = grad(2);

                    pos++;
                }
            }
        }

        assert(pos == polynomial_part_size);

        if (harmonic_part_size > 0) {
            assert(polynomial_max_degree < full_basis_degree);
            auto hep = 2.0*(pt-bar)/diam; /* beware of scaling: must be the same in all dirs */
            gradient_type hg = eval_harmonic_gradients(hep);
            ret.block(ret.rows()-harmonic_part_size, 0, harmonic_part_size, 3) =
                hg.block(hg.rows()-harmonic_part_size, 0, harmonic_part_size, 3);
        }
        return ret;
    }

    size_t size() const {
        return full_basis_size;
    }

    size_t poly_size() const {
        return polynomial_part_size;
    }

    size_t harmonic_size() const {
        return harmonic_part_size;
    }

    size_t degree() const {
        return full_basis_degree;
    }
};

#endif /* USE_H2DT */

template<typename Mesh, typename ElemType, typename ScalarType = typename Mesh::coordinate_type>
auto
make_scalar_harmonic_top_basis(const Mesh& msh, const ElemType& elem,size_t degree)
{
    return scaled_harmonic_top_scalar_basis<Mesh, ElemType, ScalarType>(msh, elem, degree);
}


} // namespace disk


