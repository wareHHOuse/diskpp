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

#include "quadratures/quad_bones.hpp"
#include "geometry/geometry_generic.hpp"

namespace disk {

template<typename T>
class quadrature<generic_mesh<T,2>, typename generic_mesh<T,2>::cell>
{
    size_t                                          m_order;
    std::vector<std::pair<point<T,2>, T>>           m_quadrature_data;

public:
    typedef generic_mesh<T,2>                       mesh_type;
    typedef typename mesh_type::cell                cell_type;
    typedef quadrature_point<T,2>                   quadpoint_type;
    typedef typename mesh_type::point_type          point_type;
    typedef T                                       weight_type;

    quadrature()
        : m_order(1)
    {
        m_quadrature_data = triangle_quadrature(1);
    }

    quadrature(size_t order)
        : m_order(order)
    {
        m_quadrature_data = triangle_quadrature(m_order);
    }

    std::vector<quadpoint_type>
    integrate(const mesh_type& msh, const cell_type& cl) const
    {
        auto pts        = points(msh, cl);
        auto c_center   = barycenter(msh, cl);

        std::vector<quadpoint_type> ret;;

        /* OPTIMIZE THE TRIANGULAR CASE!! */

        /* Break the cell in triangles, compute the transformation matrix and
         * map quadrature data in the physical space. Edges of the triangle as
         * column vectors in the transformation matrix. */
        ret.resize( m_quadrature_data.size() * pts.size() );
        for (size_t i = 0; i < pts.size(); i++)
        {
            auto pt1 = pts[i];
            auto pt2 = pts[(i+1)%pts.size()];
            auto col1 = pt1 - c_center;
            auto col2 = pt2 - c_center;

            /* Compute the area of the sub-triangle */
            auto tm = (col1.x()*col2.y() - col2.x()*col1.y())/2.;

            auto tr = [&](const std::pair<point<T,2>, T>& qd) -> auto {
                auto point = col1*qd.first.x() + col2*qd.first.y() + c_center;
                auto weight = qd.second * std::abs(tm);
                return make_qp(point, weight);
            };

            auto retbegin = ret.begin();
            std::advance(retbegin, m_quadrature_data.size()*i);

            std::transform(m_quadrature_data.begin(), m_quadrature_data.end(),
                           retbegin, tr);
        }



        return ret;
    }
};

template<typename T>
class quadrature<generic_mesh<T,2>, typename generic_mesh<T,2>::face>
{
    size_t                                          m_order;
    std::vector<std::pair<point<T,1>, T>>           m_quadrature_data;

public:
    typedef generic_mesh<T,2>                       mesh_type;
    typedef typename mesh_type::face                face_type;
    typedef quadrature_point<T,2>                   quadpoint_type;
    typedef typename mesh_type::point_type          point_type;
    typedef T                                       weight_type;

    quadrature()
        : m_order(1)
    {
        m_quadrature_data = edge_quadrature<T>(1);
    }

    quadrature(size_t order)
        : m_order(order)
    {
        m_quadrature_data = edge_quadrature<T>(m_order);
    }

    std::vector<quadpoint_type>
    integrate(const mesh_type& msh, const face_type& fc) const
    {
        auto meas = measure(msh, fc);
        auto pts = points(msh, fc);
        auto tr = [&](const std::pair<point<T,1>, T>& qd) -> auto {
            auto point = (pts[1] - pts[0])*qd.first.x() + pts[0];
            auto weight = qd.second * meas;
            return make_qp(point, weight);
        };

        std::vector<quadpoint_type> ret(m_quadrature_data.size());
        std::transform(m_quadrature_data.begin(), m_quadrature_data.end(),
                       ret.begin(), tr);

        return ret;
    }
};

template<typename T>
class quadrature<generic_mesh<T,1>, typename generic_mesh<T,1>::cell>
{
    size_t                                          m_order;
    std::vector<std::pair<point<T,1>, T>>           m_quadrature_data;

public:
    typedef generic_mesh<T,1>                       mesh_type;
    typedef typename mesh_type::cell                cell_type;
    typedef quadrature_point<T,1>                   quadpoint_type;
    typedef typename mesh_type::point_type          point_type;
    typedef T                                       weight_type;

    quadrature()
        : m_order(1)
    {
        m_quadrature_data = edge_quadrature<T>(1);
    }

    quadrature(size_t order)
        : m_order(order)
    {
        m_quadrature_data = edge_quadrature<T>(m_order);
    }

    std::vector<quadpoint_type>
    integrate(const mesh_type& msh, const cell_type& cl) const
    {
        auto pts        = points(msh, cl);
        auto c_center   = barycenter(msh, cl);
        auto meas       = measure(msh, cl);

        assert(pts.size() == 2);

        std::vector<quadpoint_type> ret;;
        ret.resize( m_quadrature_data.size() );


        auto tr = [&](const std::pair<point<T,1>, T>& qd) -> auto {
            auto point = (pts[1] - pts[0]) * qd.first.x() + pts[0];
            auto weight = qd.second * meas;
            return make_qp(point, weight);
        };

        auto retbegin = ret.begin();

        std::transform(m_quadrature_data.begin(), m_quadrature_data.end(),
                       retbegin, tr);

        return ret;
    }
};

template<typename T>
class quadrature<generic_mesh<T,1>, typename generic_mesh<T,1>::face>
{

public:
    typedef generic_mesh<T,1>                       mesh_type;
    typedef typename mesh_type::face                face_type;
    typedef quadrature_point<T,1>                   quadpoint_type;
    typedef typename mesh_type::point_type          point_type;
    typedef T                                       weight_type;

    quadrature()
    {}

    quadrature(size_t)
    {}

    std::vector<quadpoint_type>
    integrate(const mesh_type& msh, const face_type& fc) const
    {
        return std::vector<quadpoint_type>( 1, make_qp(barycenter(msh, fc), 1.) );
    }
};

} // namespace disk
