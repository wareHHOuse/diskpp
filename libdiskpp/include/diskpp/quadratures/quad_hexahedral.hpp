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

#ifndef _QUAD_HEXAHEDRAL_HPP_
#define _QUAD_HEXAHEDRAL_HPP_

#include "diskpp/quadratures/quad_raw_gauss.hpp"
namespace disk {

template<typename T>
std::vector<disk::quadrature_point<T, 2>>
integrate(const disk::cartesian_mesh<T, 2>& msh, const typename disk::cartesian_mesh<T, 2>::cell& cl, const size_t degree)
{
    const auto pts = points(msh, cl);
    return disk::quadrature::tensorized_gauss_legendre(degree, pts);
}

template<typename T>
std::vector<disk::quadrature_point<T, 2>>
integrate(const disk::cartesian_mesh<T, 2>& msh, const typename disk::cartesian_mesh<T, 2>::face& fc, const size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, fc);
    }

    return priv::integrate_2D_face(msh, fc, degree);
}

namespace priv {

template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate_hexahedron_tens(const size_t degree, const std::array<point<T, 3>, 8>& pts)
{
    assert(degree >= 0);
    const auto qps = disk::edge_quadrature<T>(degree);

    //std::cout << "deg:" << degree << " nb : " << 3 * qps.size() << std::endl;

    std::vector<disk::quadrature_point<T, 3>> ret;
    ret.reserve(3 * qps.size());

    auto P = [&pts](T xi, T eta, T zeta) -> T {
        const T val = pts[0].x() * (1 - xi) * (1 - eta) * (1 - zeta) + pts[1].x() * (1 + xi) * (1 - eta) * (1 - zeta) +
                      pts[2].x() * (1 + xi) * (1 + eta) * (1 - zeta) + pts[3].x() * (1 - xi) * (1 + eta) * (1 - zeta) +
                      pts[4].x() * (1 - xi) * (1 - eta) * (1 + zeta) + pts[5].x() * (1 + xi) * (1 - eta) * (1 + zeta) +
                      pts[6].x() * (1 + xi) * (1 + eta) * (1 + zeta) + pts[7].x() * (1 - xi) * (1 + eta) * (1 + zeta);
        return T(0.125) * val;
    };

    auto Q = [&pts](T xi, T eta, T zeta) -> T {
        const T val = pts[0].y() * (1 - xi) * (1 - eta) * (1 - zeta) + pts[1].y() * (1 + xi) * (1 - eta) * (1 - zeta) +
                      pts[2].y() * (1 + xi) * (1 + eta) * (1 - zeta) + pts[3].y() * (1 - xi) * (1 + eta) * (1 - zeta) +
                      pts[4].y() * (1 - xi) * (1 - eta) * (1 + zeta) + pts[5].y() * (1 + xi) * (1 - eta) * (1 + zeta) +
                      pts[6].y() * (1 + xi) * (1 + eta) * (1 + zeta) + pts[7].y() * (1 - xi) * (1 + eta) * (1 + zeta);
        return T(0.125) * val;
    };

    auto R = [&pts](T xi, T eta, T zeta) -> T {
        const T val = pts[0].z() * (1 - xi) * (1 - eta) * (1 - zeta) + pts[1].z() * (1 + xi) * (1 - eta) * (1 - zeta) +
                      pts[2].z() * (1 + xi) * (1 + eta) * (1 - zeta) + pts[3].z() * (1 - xi) * (1 + eta) * (1 - zeta) +
                      pts[4].z() * (1 - xi) * (1 - eta) * (1 + zeta) + pts[5].z() * (1 + xi) * (1 - eta) * (1 + zeta) +
                      pts[6].z() * (1 + xi) * (1 + eta) * (1 + zeta) + pts[7].z() * (1 - xi) * (1 + eta) * (1 + zeta);
        return T(0.125) * val;
    };

    auto J = [&pts](T xi, T eta, T zeta) -> T {
        static_matrix<T, 3, 3> Jac = static_matrix<T, 3, 3>::Zero();
        Jac(0, 0) = - pts[0].x() * (1 - eta) * (1 - zeta) + pts[1].x() * (1 - eta) * (1 - zeta)
                    + pts[2].x() * (1 + eta) * (1 - zeta) - pts[3].x() * (1 + eta) * (1 - zeta)
                    - pts[4].x() * (1 - eta) * (1 + zeta) + pts[5].x() * (1 - eta) * (1 + zeta)
                    + pts[6].x() * (1 + eta) * (1 + zeta) - pts[7].x() * (1 + eta) * (1 + zeta);

        Jac(0, 1) = - pts[0].y() * (1 - eta) * (1 - zeta) + pts[1].y() * (1 - eta) * (1 - zeta)
                    + pts[2].y() * (1 + eta) * (1 - zeta) - pts[3].y() * (1 + eta) * (1 - zeta)
                    - pts[4].y() * (1 - eta) * (1 + zeta) + pts[5].y() * (1 - eta) * (1 + zeta)
                    + pts[6].y() * (1 + eta) * (1 + zeta) - pts[7].y() * (1 + eta) * (1 + zeta);

        Jac(0, 2) = - pts[0].z() * (1 - eta) * (1 - zeta) + pts[1].z() * (1 - eta) * (1 - zeta)
                    + pts[2].z() * (1 + eta) * (1 - zeta) - pts[3].z() * (1 + eta) * (1 - zeta)
                    - pts[4].z() * (1 - eta) * (1 + zeta) + pts[5].z() * (1 - eta) * (1 + zeta)
                    + pts[6].z() * (1 + eta) * (1 + zeta) - pts[7].z() * (1 + eta) * (1 + zeta);

        Jac(1, 0) = - pts[0].x() * (1 - xi) * (1 - zeta) - pts[1].x() * (1 + xi) * (1 - zeta)
                    + pts[2].x() * (1 + xi) * (1 - zeta) + pts[3].x() * (1 - xi) * (1 - zeta)
                    - pts[4].x() * (1 - xi) * (1 + zeta) - pts[5].x() * (1 + xi) * (1 + zeta)
                    + pts[6].x() * (1 + xi) * (1 + zeta) + pts[7].x() * (1 - xi) * (1 + zeta);

        Jac(1, 1) = - pts[0].y() * (1 - xi) * (1 - zeta) - pts[1].y() * (1 + xi) * (1 - zeta)
                    + pts[2].y() * (1 + xi) * (1 - zeta) + pts[3].y() * (1 - xi) * (1 - zeta)
                    - pts[4].y() * (1 - xi) * (1 + zeta) - pts[5].y() * (1 + xi) * (1 + zeta)
                    + pts[6].y() * (1 + xi) * (1 + zeta) + pts[7].y() * (1 - xi) * (1 + zeta);

        Jac(1, 2) = - pts[0].z() * (1 - xi) * (1 - zeta) - pts[1].z() * (1 + xi) * (1 - zeta)
                    + pts[2].z() * (1 + xi) * (1 - zeta) + pts[3].z() * (1 - xi) * (1 - zeta)
                    - pts[4].z() * (1 - xi) * (1 + zeta) - pts[5].z() * (1 + xi) * (1 + zeta)
                    + pts[6].z() * (1 + xi) * (1 + zeta) + pts[7].z() * (1 - xi) * (1 + zeta);

        Jac(2, 0) = - pts[0].x() * (1 - xi) * (1 - eta) - pts[1].x() * (1 + xi) * (1 - eta)
                    - pts[2].x() * (1 + xi) * (1 + eta) - pts[3].x() * (1 - xi) * (1 + eta)
                    + pts[4].x() * (1 - xi) * (1 - eta) + pts[5].x() * (1 + xi) * (1 - eta)
                    + pts[6].x() * (1 + xi) * (1 + eta) + pts[7].x() * (1 - xi) * (1 + eta);

        Jac(2, 1) = - pts[0].y() * (1 - xi) * (1 - eta) - pts[1].y() * (1 + xi) * (1 - eta)
                    - pts[2].y() * (1 + xi) * (1 + eta) - pts[3].y() * (1 - xi) * (1 + eta)
                    + pts[4].y() * (1 - xi) * (1 - eta) + pts[5].y() * (1 + xi) * (1 - eta)
                    + pts[6].y() * (1 + xi) * (1 + eta) + pts[7].y() * (1 - xi) * (1 + eta);

        Jac(2, 2) = - pts[0].z() * (1 - xi) * (1 - eta) - pts[1].z() * (1 + xi) * (1 - eta)
                    - pts[2].z() * (1 + xi) * (1 + eta) - pts[3].z() * (1 - xi) * (1 + eta)
                    + pts[4].z() * (1 - xi) * (1 - eta) + pts[5].z() * (1 + xi) * (1 - eta)
                    + pts[6].z() * (1 + xi) * (1 + eta) + pts[7].z() * (1 - xi) * (1 + eta);

        Jac *= T(0.125);
        return std::abs(Jac.determinant());
    };

    for (auto jtor = qps.begin(); jtor != qps.end(); jtor++)
    {
        for (auto itor = qps.begin(); itor != qps.end(); itor++)
        {
            for (auto ktor = qps.begin(); ktor != qps.end(); ktor++)
            {
                const auto qp_x = *itor;
                const auto qp_y = *jtor;
                const auto qp_z = *ktor;

                const auto xi   = qp_x.first.x();
                const auto eta  = qp_y.first.x();
                const auto zeta = qp_z.first.x();

                const auto px = P(xi, eta, zeta);
                const auto py = Q(xi, eta, zeta);
                const auto pz = R(xi, eta, zeta);

                // std::cout << xi << " " << px << std::endl;
                // std::cout << eta << " " << py << std::endl;
                // std::cout << zeta << " " << pz << std::endl;
                // std::cout << J(xi, eta, zeta) << std::endl;

                const auto w = qp_x.second * qp_y.second * qp_z.second * J(xi, eta, zeta);

                ret.push_back(disk::make_qp(point<T, 3>({px, py, pz}), w));
            }
        }
    }

    return ret;
}

} // end priv

template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate(const disk::cartesian_mesh<T, 3>& msh, const typename disk::cartesian_mesh<T, 3>::cell& cl, const size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, cl);
    }

    const auto pts = points(msh, cl);
    // transform arry to vector
    std::array<point<T, 3>, 8> ptsv{pts[0], pts[1], pts[3], pts[2], pts[4], pts[5], pts[7], pts[6]};
    //std::cout << "pts : " << pts[0]<< pts[1]<< pts[3]  << pts[2]<< pts[4]<< pts[5]<< pts[7]<< pts[6] << std::endl;
    return priv::integrate_hexahedron_tens(degree, ptsv);
}

template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate(const disk::cartesian_mesh<T, 3>& msh, const typename disk::cartesian_mesh<T, 3>::face& fc, const size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, fc);
    }

    const auto pts               = points(msh, fc);
    const auto meas              = measure(msh, fc);
    const auto bar               = barycenter(msh, fc);
    const auto m_quadrature_data = disk::quadrangle_quadrature<T>(degree);

    std::vector<disk::quadrature_point<T, 3>> ret(m_quadrature_data.size());

    const auto v0 = (pts[1] - pts[0]) / T(2);
    const auto v1 = (pts[3] - pts[0]) / T(2);
    const T    meas4 = meas / T(4);

    for (auto& qp : m_quadrature_data)
    {
        const auto point  = v0 * qp.first.x() + v1 * qp.first.y() + bar;
        const auto weight = qp.second * meas4;
        ret.push_back(disk::make_qp(point, weight));
    }

    return ret;
}

} // namespace disk


#endif /* _QUAD_GENERIC_HPP_ */
