/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#pragma once

#include <vector>

#include "geometry/geometry.hpp"
#include "common/eigen.hpp"
#include "bases/bases_bones.hpp"
#include "bases/bases_templates.hpp"
#include "bases/monomial_generator.hpp"


namespace disk {

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_scalar_basis<Mesh<T,2,Storage>, typename Mesh<T,2,Storage>::cell>
    : public priv::scaled_monomial_basis_base<2>
{
    typedef Mesh<T,2,Storage>                       mesh_type;
    typedef typename mesh_type::cell                cell_type;
    typedef typename mesh_type::point_type          point_type;
    typedef priv::scaled_monomial_basis_base<2>     base;

public:
    typedef T                                       function_value_type;
    typedef static_vector<T,2>                      gradient_value_type;

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    std::vector<function_value_type>
    eval_functions(const mesh_type& msh, const cell_type& cl, const point_type& pt) const
    {
        auto bar = barycenter(msh, cl);
        auto h = measure(msh, cl);

        auto ep = (pt - bar)/h;

        std::vector<function_value_type> ret;
        ret.reserve( this->size() );

        for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
        {
            auto m = *itor;
            auto vx = iexp_pow(ep.x(), m[0]);
            auto vy = iexp_pow(ep.y(), m[1]);
            ret.push_back( vx * vy );
        }

        return ret;
    }

    std::vector<gradient_value_type>
    eval_gradients(const mesh_type& msh, const cell_type& cl, const point_type& pt) const
    {
        auto bar = barycenter(msh, cl);
        auto h = measure(msh, cl);

        auto ep = (pt - bar)/h;

        std::vector<gradient_value_type> ret;
        ret.reserve( this->size() );

        for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
        {
            auto m = *itor;
            gradient_value_type grad;

            auto px = iexp_pow(ep.x(), m[0]);
            auto py = iexp_pow(ep.y(), m[1]);

            auto dx = (m[0] == 0) ? 0 : (m[0]/h)*iexp_pow(ep.x(), m[0]-1);
            auto dy = (m[1] == 0) ? 0 : (m[1]/h)*iexp_pow(ep.y(), m[1]-1);

            grad(0) = dx * py;
            grad(1) = dy * px;

            ret.push_back(grad);
        }

        return ret;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_scalar_basis<Mesh<T,2,Storage>, typename Mesh<T,2,Storage>::face>
    : public priv::scaled_monomial_basis_base<1>
{
    typedef Mesh<T,2,Storage>                           mesh_type;
    typedef typename mesh_type::point_type              point_type;
    typedef typename mesh_type::face                    face_type;
    typedef priv::scaled_monomial_basis_base<1>         base;

public:
    typedef T                                           function_value_type;

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    std::vector<function_value_type>
    eval_functions(const mesh_type& msh, const face_type& fc, const point_type& pt) const
    {
        auto pts = points(msh, fc);
        auto bar = barycenter(msh, fc);
        auto h = measure(msh, fc);
        auto v = (pts[1] - pts[0]).to_vector();
        auto t = (pt - bar).to_vector();
        //v = v/v.norm();
        T dot = v.dot(t);
        auto ep = point<T, 1>({dot/(h*h)});

        std::vector<function_value_type> ret;
        ret.reserve( this->size() );

        for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
        {
            auto m = *itor;
            auto vx = iexp_pow(ep.x(), m[0]);
            ret.push_back( vx );
        }

        return ret;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_scalar_basis<Mesh<T,1,Storage>, typename Mesh<T,1,Storage>::cell>
    : public priv::scaled_monomial_basis_base<1>
{
    typedef Mesh<T,1,Storage>                       mesh_type;
    typedef typename mesh_type::cell                cell_type;
    typedef typename mesh_type::point_type          point_type;
    typedef priv::scaled_monomial_basis_base<1>     base;

public:
    typedef T                                       function_value_type;
    typedef T                                       gradient_value_type;

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    std::vector<function_value_type>
    eval_functions(const mesh_type& msh, const cell_type& cl, const point_type& pt) const
    {
        auto bar = barycenter(msh, cl);
        auto h = measure(msh, cl);

        auto ep = (pt - bar)/h;

        std::vector<function_value_type> ret;
        ret.reserve( this->size() );

        for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
        {
            auto m = *itor;
            auto vx = iexp_pow(ep.x(), m[0]);
            ret.push_back( vx );
        }

        return ret;
    }

    std::vector<gradient_value_type>
    eval_gradients(const mesh_type& msh, const cell_type& cl, const point_type& pt) const
    {
        auto bar = barycenter(msh, cl);
        auto h = measure(msh, cl);

        auto ep = (pt - bar)/h;

        std::vector<gradient_value_type> ret;
        ret.reserve( this->size() );

        for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
        {
            auto m = *itor;

            auto dx = (m[0] == 0) ? 0 : (m[0]/h)*iexp_pow(ep.x(), m[0]-1);

            ret.push_back(dx);
        }

        return ret;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_scalar_basis<Mesh<T,1,Storage>, typename Mesh<T,1,Storage>::face>
    : public priv::scaled_monomial_basis_base<0>
{
    typedef Mesh<T,1,Storage>                           mesh_type;
    typedef typename mesh_type::point_type              point_type;
    typedef typename mesh_type::face                    face_type;
    typedef priv::scaled_monomial_basis_base<0>         base;

public:
    typedef T                                           function_value_type;

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    std::vector<function_value_type>
    eval_functions(const mesh_type& msh, const face_type& fc, const point_type& pt) const
    {
        return std::vector<function_value_type>(1,1);
    }
};


template<template<typename, size_t, typename> class Mesh,
                  typename T, typename Storage>
std::vector<point<T,2>>
make_test_points(const Mesh<T,2,Storage>& msh,
                 const typename Mesh<T,2,Storage>::cell& cl)
{
    std::vector<point<T,2>> test_points;
    auto pts = points(msh, cl);
    auto bar = barycenter(msh, cl);

    test_points.insert(test_points.begin(), pts.begin(), pts.end()); //vertices
    test_points.push_back(bar); //barycenter

    //split in triangles and take barycenter
    for (size_t i = 0; i < pts.size(); i++)
        test_points.push_back( (pts[i] + pts[(i+1) % pts.size()] + bar) / 3. );

    return test_points;
}

template<template<typename, size_t, typename> class Mesh,
                  typename T, typename Storage>
std::vector<point<T,2>>
make_test_points(const Mesh<T,2,Storage>& msh,
                 const typename Mesh<T,2,Storage>::face& fc)
{
    std::vector<point<T,2>> test_points(5);
    auto pts = points(msh, fc);
    assert(pts.size() == 2); /* we are in 2D and we are dealing with a face,
                                which is an edge */

    test_points[0] = pts[0];
    test_points[1] = pts[1];
    test_points[2] = (pts[0] + pts[1]) / 2.;
    test_points[3] = (test_points[2] + pts[0]) / 2.;
    test_points[4] = (test_points[2] + pts[1]) / 2.;

    return test_points;
}

template<template<typename, size_t, typename> class Mesh,
                  typename T, typename Storage>
std::vector<point<T,1>>
make_test_points(const Mesh<T,1,Storage>& msh,
                 const typename Mesh<T,1,Storage>::cell& cl)
{
    std::vector<point<T,1>> test_points;
    auto pts = points(msh, cl);
    auto bar = barycenter(msh, cl);

    test_points.insert(test_points.begin(), pts.begin(), pts.end()); //vertices
    test_points.push_back(bar); //barycenter
    test_points.push_back((pts[0] + bar) / 2.);
    test_points.push_back((pts[1] + bar) / 2.);

    return test_points;
}

} // namespace disk
