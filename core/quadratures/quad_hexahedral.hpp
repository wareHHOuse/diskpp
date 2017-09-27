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

#ifndef _QUAD_HEXAHEDRAL_HPP_
#define _QUAD_HEXAHEDRAL_HPP_


namespace disk {

template<typename T>
class quadrature<cartesian_mesh<T,3>, typename cartesian_mesh<T,3>::cell>
{
    std::vector<std::pair<point<T,1>, T>> m_quadrature_data;

public:
    typedef cartesian_mesh<T,3>                     mesh_type;
    typedef typename mesh_type::cell                cell_type;
    typedef quadrature_point<T,3>                   quadpoint_type;
    typedef typename mesh_type::point_type          point_type;
    typedef T                                       weight_type;

private:

    size_t m_order;

public:
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

        assert(pts.size() == 8);

        size_t qd_size = m_quadrature_data.size();
        std::vector<std::pair<point<T,3>, T>> qd;
        qd.reserve(qd_size * qd_size * qd_size);

        for (size_t i = 0; i < qd_size; i++)
        {
            auto px = m_quadrature_data[i].first;
            auto wx = m_quadrature_data[i].second;

            for (size_t j = 0; j < qd_size; j++)
            {
                auto py = m_quadrature_data[j].first;
                auto wy = m_quadrature_data[j].second;

                for (size_t k = 0; k < qd_size; k++)
                {
                    auto pz = m_quadrature_data[k].first;
                    auto wz = m_quadrature_data[k].second;

                    auto pt = point<T,3>({px[0], py[0], pz[0]});
                    auto wt = wx * wy * wz;

                    qd.push_back( std::make_pair(pt, wt) );
                }
            }
        }

        auto tr = [&](const std::pair<point<T,3>, T>& qd) -> auto {

            auto point = (pts[1] - pts[0]) * qd.first.x() +
                         (pts[2] - pts[0]) * qd.first.y() +
                         (pts[4] - pts[0]) * qd.first.z() +
                          pts[0];

            auto weight = qd.second * meas;
            return make_qp(point, weight);
        };

        std::vector<quadpoint_type> ret;
        ret.resize(qd_size * qd_size * qd_size);
        std::transform(qd.begin(), qd.end(), ret.begin(), tr);

        return ret;
    }
};

template<typename T>
class quadrature<cartesian_mesh<T,3>, typename cartesian_mesh<T,3>::face>
{
    std::vector<std::pair<point<T,1>, T>> m_quadrature_data;

public:
    typedef cartesian_mesh<T,3>                     mesh_type;
    typedef typename mesh_type::face                face_type;
    typedef quadrature_point<T,3>                   quadpoint_type;
    typedef T                                       weight_type;

private:
    size_t m_order;

public:
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
        auto pts        = points(msh, fc);
        auto meas       = measure(msh, fc);

        assert(pts.size() == 4);

        size_t qd_size = m_quadrature_data.size();
        std::vector<std::pair<point<T,2>, T>> qd;
        qd.reserve(qd_size * qd_size);

        for (size_t i = 0; i < qd_size; i++)
        {
            auto px = m_quadrature_data[i].first;
            auto wx = m_quadrature_data[i].second;

            for (size_t j = 0; j < qd_size; j++)
            {
                auto py = m_quadrature_data[j].first;
                auto wy = m_quadrature_data[j].second;

                auto pt = point<T,2>({px[0], py[0]});
                auto wt = wx * wy;

                qd.push_back( std::make_pair(pt, wt) );

            }
        }

        auto tr = [&](const std::pair<point<T,2>, T>& qd) -> auto {

            auto point = (pts[1] - pts[0]) * qd.first.x() +
                         (pts[3] - pts[0]) * qd.first.y() +
                          pts[0];

            auto weight = qd.second * meas;
            return make_qp(point, weight);
        };

        std::vector<quadpoint_type> ret;
        ret.resize(qd_size * qd_size);
        std::transform(qd.begin(), qd.end(), ret.begin(), tr);

        return ret;
    }
};








template<typename T>
class quadrature<cartesian_mesh<T,2>, typename cartesian_mesh<T,2>::cell>
{
    std::vector<std::pair<point<T,1>, T>> m_quadrature_data;

public:
    typedef cartesian_mesh<T,2>                     mesh_type;
    typedef typename mesh_type::cell                cell_type;
    typedef quadrature_point<T,2>                   quadpoint_type;
    typedef typename mesh_type::point_type          point_type;
    typedef T                                       weight_type;

private:

    size_t m_order;

public:
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

        assert(pts.size() == 4);

        size_t qd_size = m_quadrature_data.size();
        std::vector<std::pair<point<T,2>, T>> qd;
        qd.reserve(qd_size * qd_size);

        for (size_t i = 0; i < qd_size; i++)
        {
            auto px = m_quadrature_data[i].first;
            auto wx = m_quadrature_data[i].second;

            for (size_t j = 0; j < qd_size; j++)
            {
                auto py = m_quadrature_data[j].first;
                auto wy = m_quadrature_data[j].second;

                auto pt = point<T,2>({px[0], py[0]});
                auto wt = wx * wy;

                qd.push_back( std::make_pair(pt, wt) );

            }
        }

        auto tr = [&](const std::pair<point<T,2>, T>& qd) -> auto {

            auto point = (pts[1] - pts[0]) * qd.first.x() +
                         (pts[2] - pts[0]) * qd.first.y() +
                          pts[0];

            auto weight = qd.second * meas;
            return make_qp(point, weight);
        };

        std::vector<quadpoint_type> ret;
        ret.resize(qd_size * qd_size);
        std::transform(qd.begin(), qd.end(), ret.begin(), tr);

        return ret;
    }
};

template<typename T>
class quadrature<cartesian_mesh<T,2>, typename cartesian_mesh<T,2>::face>
{
    std::vector<std::pair<point<T,1>, T>> m_quadrature_data;

public:
    typedef cartesian_mesh<T,2>                     mesh_type;
    typedef typename mesh_type::face                face_type;
    typedef quadrature_point<T,2>                   quadpoint_type;
    typedef T                                       weight_type;

private:
    size_t m_order;

public:
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
        auto pts        = points(msh, fc);
        auto meas       = measure(msh, fc);

        assert(pts.size() == 2);

        auto tr = [&](const std::pair<point<T,1>, T>& qd) -> auto {

            auto point = (pts[1] - pts[0]) * qd.first.x() + pts[0];
            auto weight = qd.second * meas;
            return make_qp(point, weight);
        };

        std::vector<quadpoint_type> ret;
        ret.resize(m_quadrature_data.size());
        std::transform(m_quadrature_data.begin(), m_quadrature_data.end(),
                       ret.begin(), tr);

        return ret;
    }
};









} // namespace disk


#endif /* _QUAD_GENERIC_HPP_ */
