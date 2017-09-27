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

#ifndef _QUADRATURES_HPP_WAS_INCLUDED_
    #error "You must NOT include this file. Include quadratures.hpp"
#endif

#ifndef _QUAD_GENERIC_HPP_
#define _QUAD_GENERIC_HPP_

#include "raw_simplices.hpp"

namespace disk {

template<typename T>
class quadrature<generic_mesh<T,3>, typename generic_mesh<T,3>::cell>
{
    size_t                                          m_order;
    std::vector<std::pair<point<T,3>, T>>           m_quadrature_data;

public:
    typedef generic_mesh<T,3>                       mesh_type;
    typedef typename mesh_type::cell                cell_type;
    typedef quadrature_point<T,3>                   quadpoint_type;
    typedef typename mesh_type::point_type          point_type;
    typedef T                                       weight_type;
    typedef T                                       scalar_type;

    quadrature()
        : m_order(1)
    {
        m_quadrature_data = tetrahedron_quadrature(1);
    }

    quadrature(size_t order)
        : m_order(order)
    {
        m_quadrature_data = tetrahedron_quadrature(m_order);
    }

    std::vector<quadpoint_type>
    integrate(const mesh_type& msh, const cell_type& cl) const
    {
        auto rss = split_in_raw_tetrahedra(msh, cl);

        std::vector<quadpoint_type> ret;
        for (auto& rs : rss)
        {
            auto meas = measure(rs);
            for (auto& qd : m_quadrature_data)
            {
                auto point = map_from_reference(qd.first, rs);
                auto weight = qd.second * meas;
                ret.push_back( make_qp(point, weight) );
            }
        }

        return ret;
    }

};

template<typename T>
class quadrature<generic_mesh<T,3>, typename generic_mesh<T,3>::face>
{
    size_t                                          m_order;
    std::vector<std::pair<point<T,2>, T>>           m_quadrature_data;

public:
    typedef generic_mesh<T,3>                       mesh_type;
    typedef typename mesh_type::face                face_type;
    typedef quadrature_point<T,3>                   quadpoint_type;
    typedef typename mesh_type::point_type          point_type;
    typedef T                                       weight_type;
    typedef T                                       scalar_type;

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
    integrate(const mesh_type& msh, const face_type& fc) const
    {
        auto rss = split_in_raw_triangles(msh, fc);

        std::vector<quadpoint_type> ret;
        for (auto& rs : rss)
        {
            auto meas = measure(rs);
            for (auto& qd : m_quadrature_data)
            {
                auto point = map_from_reference(qd.first, rs);
                auto weight = qd.second * meas;
                ret.push_back( make_qp(point, weight) );
            }
        }

        return ret;
    }

};

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

private:

    template<typename PtA>
    std::vector<quadpoint_type>
    integrate_triangle(const mesh_type& msh, const cell_type& cl,
                       const PtA& pts) const
    {
        std::vector<quadpoint_type> ret;

        ret.resize( m_quadrature_data.size() );

        auto col1 = pts[1] - pts[0];
        auto col2 = pts[2] - pts[0];

        /* Compute the area of the sub-triangle */
        auto tm = (col1.x()*col2.y() - col2.x()*col1.y())/2.;

        auto tr = [&](const std::pair<point<T,2>, T>& qd) -> auto {
            auto point = col1*qd.first.x() + col2*qd.first.y() + pts[0];
            auto weight = qd.second * std::abs(tm);
            return make_qp(point, weight);
        };

        auto retbegin = ret.begin();

        std::transform(m_quadrature_data.begin(), m_quadrature_data.end(),
                       retbegin, tr);

        return ret;
    }


    template<typename PtA>
    std::vector<quadpoint_type>
    integrate_quad(const mesh_type& msh, const cell_type& cl,
                       const PtA& pts) const
    {
        std::vector<quadpoint_type> ret;

        ret.resize( m_quadrature_data.size() * 2 );
        for (size_t i = 1; i < 3; i++)
        {
            auto pt1 = pts[i];
            auto pt2 = pts[i+1];
            auto col1 = pt1 - pts[0];
            auto col2 = pt2 - pts[0];

            /* Compute the area of the sub-triangle */
            auto tm = (col1.x()*col2.y() - col2.x()*col1.y())/2.;

            auto tr = [&](const std::pair<point<T,2>, T>& qd) -> auto {
                auto point = col1*qd.first.x() + col2*qd.first.y() + pts[0];
                auto weight = qd.second * std::abs(tm);
                return make_qp(point, weight);
            };

            auto retbegin = ret.begin();
            std::advance(retbegin, m_quadrature_data.size()*(i-1));

            std::transform(m_quadrature_data.begin(), m_quadrature_data.end(),
                           retbegin, tr);
        }

        return ret;
    }


//#define OPTIMAL_TRIANGLE_NUMBER

    /* The 'optimal triangle number' version gives almost the same results
     * of the other version. In bigger meshes there are some advantages in
     * assembly time. The problem with this version is that it could generate
     * triangles with a very bad aspect ratio and at the moment I don't know
     * if and how this can affect computations, so I leave it turned off.
     */
#ifdef OPTIMAL_TRIANGLE_NUMBER
    template<typename PtA>
    std::vector<quadpoint_type>
    integrate_other(const mesh_type& msh, const cell_type& cl,
                    const PtA& pts) const
    {
        std::vector<quadpoint_type> ret;

        /* Break the cell in triangles, compute the transformation matrix and
         * map quadrature data in the physical space. Edges of the triangle as
         * column vectors in the transformation matrix. */
        ret.resize( m_quadrature_data.size() * pts.size()-2 );
        for (size_t i = 1; i < pts.size()-1; i++)
        {
            auto pt1 = pts[i];
            auto pt2 = pts[i+1];
            auto col1 = pt1 - pts[0];
            auto col2 = pt2 - pts[0];

            /* Compute the area of the sub-triangle */
            auto tm = (col1.x()*col2.y() - col2.x()*col1.y())/2.;

            auto tr = [&](const std::pair<point<T,2>, T>& qd) -> auto {
                auto point = col1*qd.first.x() + col2*qd.first.y() + pts[0];
                auto weight = qd.second * std::abs(tm);
                return make_qp(point, weight);
            };

            auto retbegin = ret.begin();
            std::advance(retbegin, m_quadrature_data.size()*(i-1));

            std::transform(m_quadrature_data.begin(), m_quadrature_data.end(),
                           retbegin, tr);
        }

        return ret;
    }
#else
    template<typename PtA>
    std::vector<quadpoint_type>
    integrate_other(const mesh_type& msh, const cell_type& cl,
                    const PtA& pts) const
    {
        auto c_center   = barycenter(msh, cl);

        std::vector<quadpoint_type> ret;;

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
#endif

public:
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

        switch(pts.size())
        {
            case 3:
                return integrate_triangle(msh, cl, pts);

            case 4:
                return integrate_quad(msh, cl, pts);

            default:
                return integrate_other(msh, cl, pts);
        }

        throw std::logic_error("Shouldn't have arrived here");
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

#endif /* _QUAD_GENERIC_HPP_ */
