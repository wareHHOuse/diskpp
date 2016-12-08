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


#include "common/eigen.hpp"

namespace disk {

template<typename T>
class scaled_monomial_scalar_basis<simplicial_mesh<T,3>, typename simplicial_mesh<T,3>::cell>
    : public priv::monomial_basis_bones<3>
{
    typedef disk::simplicial_mesh<T, 3>             mesh_type;
    typedef typename mesh_type::scalar_type         scalar_type;
    typedef typename mesh_type::cell                cell_type;
    typedef priv::monomial_basis_bones<3>           base;

public:
    typedef T                                       function_value_type;
    typedef static_vector<T,3>                      gradient_value_type;

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    eval_functions(const mesh_type& msh, const cell_type& cl,
                   const point<T,3>& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? max_degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto bar = barycenter(msh, cl);
        auto h = diameter(msh, cl);

        auto ep = (pt - bar)/h;

        Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> ret;
        ret.resize( eval_range.size(), 1 );

#ifdef POWER_CACHE
        power_cache<scalar_type, 3> pow(ep, this->max_degree()+1);
#endif
        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());
        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;

#ifdef POWER_CACHE
            auto vx = pow.x(m[0]);
            auto vy = pow.y(m[1]);
            auto vz = pow.z(m[2]);
#else
            auto vx = iexp_pow(ep.x(), m[0]);
            auto vy = iexp_pow(ep.y(), m[1]);
            auto vz = iexp_pow(ep.z(), m[2]);
#endif

            ret(i++) = vx * vy * vz;
        }

        return ret;
    }

    Eigen::Matrix<scalar_type, Eigen::Dynamic, 3>
    eval_gradients(const mesh_type& msh, const cell_type& cl,
                   const point<T,3>& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? max_degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto bar = barycenter(msh, cl);
        auto h = diameter(msh, cl);
        auto ih = 1./h;

        auto ep = (pt - bar)/h;

        Eigen::Matrix<scalar_type, Eigen::Dynamic, 3> ret;
        ret.resize( eval_range.size(), 3 );

#ifdef POWER_CACHE
        power_cache<scalar_type, 3> pow(ep, this->max_degree()+1);
#endif

        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());
        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;
            gradient_value_type grad;

#ifdef POWER_CACHE
            auto px = pow.x(m[0]);
            auto py = pow.y(m[1]);
            auto pz = pow.z(m[2]);

            auto dx = (m[0] == 0) ? 0 : m[0]*ih*pow.x(m[0]-1);
            auto dy = (m[1] == 0) ? 0 : m[1]*ih*pow.y(m[1]-1);
            auto dz = (m[2] == 0) ? 0 : m[2]*ih*pow.z(m[2]-1);
#else
            auto px = iexp_pow(ep.x(), m[0]);
            auto py = iexp_pow(ep.y(), m[1]);
            auto pz = iexp_pow(ep.z(), m[2]);

            auto dx = (m[0] == 0) ? 0 : m[0]*ih*iexp_pow(ep.x(), m[0]-1);
            auto dy = (m[1] == 0) ? 0 : m[1]*ih*iexp_pow(ep.y(), m[1]-1);
            auto dz = (m[2] == 0) ? 0 : m[2]*ih*iexp_pow(ep.z(), m[2]-1);
#endif

            ret(i,0) = dx * py * pz;
            ret(i,1) = dy * px * pz;
            ret(i,2) = dz * px * py;

            i++;
        }

        return ret;
    }
};

template<typename T>
class scaled_monomial_scalar_basis<simplicial_mesh<T,3>, typename simplicial_mesh<T,3>::face>
    : public priv::monomial_basis_bones<2>
{
    typedef disk::simplicial_mesh<T, 3>     mesh_type;
    typedef typename mesh_type::scalar_type         scalar_type;
    typedef typename mesh_type::face        face_type;

    typedef priv::monomial_basis_bones<2>   base;

public:
    typedef T                               function_value_type;

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    eval_functions(const mesh_type& msh, const face_type& fc,
                   const point<T,3>& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? max_degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto ep = map_to_reference(msh, fc, pt);

        Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> ret;
        ret.resize( eval_range.size(), 1 );

#ifdef POWER_CACHE
        power_cache<scalar_type, 2> pow(ep, this->max_degree()+1);
#endif

        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());
        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;

#ifdef POWER_CACHE
            auto vx = pow.x(m[0]);
            auto vy = pow.y(m[1]);
#else
            auto vx = iexp_pow(ep.x(), m[0]);
            auto vy = iexp_pow(ep.y(), m[1]);
#endif

            ret(i++) = vx * vy;
        }

        return ret;
    }
};


template<typename T>
class scaled_monomial_vector_sg_basis<simplicial_mesh<T,3>, typename simplicial_mesh<T,3>::cell>
    : public priv::monomial_basis_bones<3,3>
{
    typedef disk::simplicial_mesh<T, 3>             mesh_type;
    typedef typename mesh_type::cell                cell_type;
    typedef priv::monomial_basis_bones<3,3>         base;

public:
    typedef static_vector<T,3>              function_value_type;
    typedef static_matrix<T,3,3>            gradient_value_type;

    scaled_monomial_vector_sg_basis()
        : base(1)
    {}

    scaled_monomial_vector_sg_basis(size_t degree)
        : base(degree)
    {}

    std::vector<function_value_type>
    eval_functions(const mesh_type& msh, const cell_type& cl, const point<T,3>& pt) const
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
            auto vz = iexp_pow(ep.z(), m[2]);
            auto val = vx * vy * vz;
            ret.push_back( static_vector<T,3>({val,   0,   0}) );
            ret.push_back( static_vector<T,3>({  0, val,   0}) );
            ret.push_back( static_vector<T,3>({  0,   0, val}) );
        }

        return ret;
    }

    std::vector<gradient_value_type>
    eval_gradients(const mesh_type& msh, const cell_type& cl, const point<T,3>& pt) const
    {
        auto bar = barycenter(msh, cl);
        auto h = measure(msh, cl);

        auto ep = (pt - bar)/h;

        std::vector<gradient_value_type> ret;
        ret.reserve( this->size() );

        for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
        {
            auto m = *itor;

            auto px = iexp_pow(ep.x(), m[0]);
            auto py = iexp_pow(ep.y(), m[1]);
            auto pz = iexp_pow(ep.z(), m[2]);

            auto dx = (m[0] == 0) ? 0 : (m[0]/h)*iexp_pow(ep.x(), m[0]-1);
            auto dy = (m[1] == 0) ? 0 : (m[1]/h)*iexp_pow(ep.y(), m[1]-1);
            auto dz = (m[2] == 0) ? 0 : (m[2]/h)*iexp_pow(ep.z(), m[2]-1);

            gradient_value_type sg;
            sg = gradient_value_type::Zero();
            sg(0,0) = dx * py * pz;
            sg(0,1) = px * dy * pz;
            sg(0,2) = px * py * dz;
            ret.push_back(0.5*(sg + sg.transpose()));

            sg = gradient_value_type::Zero();
            sg(1,0) = dx * py * pz;
            sg(1,1) = px * dy * pz;
            sg(1,2) = px * py * dz;
            ret.push_back(0.5*(sg + sg.transpose()));

            sg = gradient_value_type::Zero();
            sg(2,0) = dx * py * pz;
            sg(2,1) = px * dy * pz;
            sg(2,2) = px * py * dz;
            ret.push_back(0.5*(sg + sg.transpose()));
        }

        return ret;
    }
};

template<typename T>
static_vector<T,3>
map_to_physical(const simplicial_mesh<T,3>& msh,
                const simplicial_element<3,1>& face,
                const static_vector<T,2>& pm)
{
    auto pts = points(msh, face);
    auto ret = pts[0] + (pts[1]-pts[0])*pm(0) + (pts[2]-pts[0])*pm(1);
    return ret.to_vector();
}

template<typename T>
class scaled_monomial_vector_sg_basis<simplicial_mesh<T,3>, typename simplicial_mesh<T,3>::face>
    : public priv::monomial_basis_bones<2,3>

{
    typedef disk::simplicial_mesh<T, 3>         mesh_type;
    typedef typename mesh_type::face            face_type;
    typedef priv::monomial_basis_bones<2,3>     base;

public:
    typedef static_vector<T,3>                      function_value_type;

    scaled_monomial_vector_sg_basis()
        : base(1)
    {}

    scaled_monomial_vector_sg_basis(size_t degree)
        : base(degree)
    {}

    std::vector<function_value_type>
    eval_functions(const mesh_type& msh, const face_type& fc, const point<T,3>& pt) const
    {
        //auto bar = barycenter(msh, fc);
        //auto h = measure(msh, fc);

        auto ep = map_to_reference(msh, fc, pt);

        std::vector<function_value_type> ret;
        ret.reserve( this->size() );

        for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
        {
            auto m = *itor;
            auto vx = iexp_pow(ep.x(), m[0]);
            auto vy = iexp_pow(ep.y(), m[1]);

            auto val = vx * vy;

            ret.push_back( static_vector<T,3>({val,   0,   0}) );
            ret.push_back( static_vector<T,3>({  0, val,   0}) );
            ret.push_back( static_vector<T,3>({  0,   0, val}) );
        }

        return ret;
    }
};



template<typename T>
std::vector<point<T,3>>
make_test_points(const simplicial_mesh<T,3>& msh,
                 const typename simplicial_mesh<T,3>::cell& cl)
{
    std::vector<point<T,3>> test_points(9);
    auto pts = points(msh, cl);
    auto bar = barycenter(msh, cl);

    test_points[0] = pts[0];
    test_points[1] = pts[1];
    test_points[2] = pts[2];
    test_points[3] = pts[3];
    test_points[4] = bar;
    test_points[5] = (pts[0] + pts[1] + pts[2] + bar)/4.;
    test_points[6] = (pts[0] + pts[1] + pts[3] + bar)/4.;
    test_points[7] = (pts[1] + pts[2] + pts[3] + bar)/4.;
    test_points[8] = (pts[0] + pts[2] + pts[3] + bar)/4.;

    return test_points;
}

template<typename T>
std::vector<point<T,3>>
make_test_points(const simplicial_mesh<T,3>& msh,
                 const typename simplicial_mesh<T,3>::face& fc)
{
    std::vector<point<T,3>> test_points(7);
    auto pts = points(msh, fc);
    auto bar = barycenter(msh, fc);

    test_points[0] = pts[0];
    test_points[1] = pts[1];
    test_points[2] = pts[2];
    test_points[3] = bar;
    test_points[4] = (pts[0] + pts[1] + bar)/3.;
    test_points[5] = (pts[1] + pts[2] + bar)/3.;
    test_points[6] = (pts[0] + pts[2] + bar)/3.;

    return test_points;
}

template<typename T>
std::vector<point<T,2>>
make_test_points(const simplicial_mesh<T,2>& msh,
                 const typename simplicial_mesh<T,2>::cell& cl,
                 size_t levels)
{
    std::vector<point<T,2>> test_points;
    auto pts = points(msh, cl);

    auto d0 = (pts[0] - pts[2]) / levels;
    auto d1 = (pts[1] - pts[2]) / levels;

    test_points.push_back( pts[2] );
    for (size_t i = 1; i < levels; i++)
    {
        auto p0 = pts[2] + i * d0;
        auto p1 = pts[2] + i * d1;

        auto d2 = (p1 - p0) / i;

        for (size_t j = 0; j <= i; j++)
            test_points.push_back( p0 + j * d2 );
    }

    //std::cout << "Test points: " << test_points.size() << std::endl;

    return test_points;
}


} // namespace disk
