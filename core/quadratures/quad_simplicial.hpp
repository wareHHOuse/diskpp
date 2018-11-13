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

#ifndef _QUADRATURES_HPP_WAS_INCLUDED_
#error "You must NOT include this file. Include quadratures.hpp"
#endif

#ifndef _QUAD_SIMPLICIAL_HPP_
#define _QUAD_SIMPLICIAL_HPP_


namespace disk {

template<typename T>
std::vector<disk::quadrature_point<T, 2>>
integrate(const disk::simplicial_mesh<T, 2>& msh,
          const typename disk::simplicial_mesh<T, 2>::cell& cl,
          const size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, cl);
    }

    const auto pts = points(msh, cl);
    return priv::integrate_triangle<T>(degree, pts);
}

template<typename T>
std::vector<disk::quadrature_point<T, 2>>
integrate(const disk::simplicial_mesh<T, 2>& msh,
          const typename disk::simplicial_mesh<T, 2>::face& fc,
          const size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, fc);
    }

    return priv::integrate_2D_face(msh, fc, degree);
}

namespace priv
{

template<typename T>
point<T, 3>
map_to_physical(const disk::simplicial_mesh<T, 3>&                msh,
                const typename disk::simplicial_mesh<T, 3>::face& face,
                const point<T, 2>&                                pm)
{
    const auto pts = points(msh, face);
    return pts[0] + (pts[1] - pts[0]) * pm.x() + (pts[2] - pts[0]) * pm.y();
}

template<typename T>
point<T, 3>
map_to_physical(const disk::simplicial_mesh<T, 3>&                msh,
                const typename disk::simplicial_mesh<T, 3>::cell& cell,
                const point<T, 3>&                                p)
{
    const auto pts = points(msh, cell);

    const auto pp = (pts[1] - pts[0]) * p.x() +
                    (pts[2] - pts[0]) * p.y() +
                    (pts[3] - pts[0]) * p.z() + pts[0];

    return pp;
}

} // namespace priv


template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate(const disk::simplicial_mesh<T, 3>& msh,
          const typename disk::simplicial_mesh<T, 3>::face& fc,
          size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, fc);
    }

    const auto m_quadrature_data = disk::triangle_quadrature(degree);

    const auto pts = points(msh, fc);
    assert(pts.size() == 3);
    const auto meas = measure(msh, fc);

    const auto col1 = pts[1] - pts[0];
    const auto col2 = pts[2] - pts[0];

    auto tr = [&](const std::pair<point<T, 2>, T>& qd) -> auto
    {
        const auto point  = col1 * qd.first.x() + col2 * qd.first.y() + pts[0];
        const auto weight = qd.second * meas;
        return disk::make_qp(point, weight);
    };

    std::vector<disk::quadrature_point<T, 3>> ret(m_quadrature_data.size());
    std::transform(m_quadrature_data.begin(), m_quadrature_data.end(), ret.begin(), tr);

    return ret;
}

template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate(const disk::simplicial_mesh<T, 3>& msh,
          const typename disk::simplicial_mesh<T, 3>::cell& cl,
          size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, cl);
    }

    const auto m_quadrature_data = disk::tetrahedron_quadrature(degree);

    const auto pts = points(msh, cl);
    assert(pts.size() == 4);
    const auto meas = measure(msh, cl);

    const auto col1 = pts[1] - pts[0];
    const auto col2 = pts[2] - pts[0];
    const auto col3 = pts[3] - pts[0];

    auto tr = [&](const std::pair<point<T, 3>, T>& qd) -> auto
    {
        const auto point  = col1 * qd.first.x() +
                            col2 * qd.first.y() +
                            col3 * qd.first.z() + pts[0];
        const auto weight = qd.second * meas;
        return disk::make_qp(point, weight);
    };

    std::vector<disk::quadrature_point<T, 3>> ret(m_quadrature_data.size());
    std::transform(m_quadrature_data.begin(), m_quadrature_data.end(), ret.begin(), tr);

    return ret;
}

} // namespace disk

#endif /* _QUAD_SIMPLICIAL_HPP_ */
