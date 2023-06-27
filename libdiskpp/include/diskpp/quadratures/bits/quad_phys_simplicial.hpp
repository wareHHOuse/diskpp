/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */
 
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
 * Nicolas Pignet  (C) 2018, 2019               nicolas.pignet@enpc.fr
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

#include "diskpp/common/simplicial_formula.hpp"
#include "quad_raw_gauss.hpp"
#include "quad_raw_dunavant.hpp"

namespace disk
{

/**
 * @brief Return quadrature points on the interior of a physical 2D triangle
 *
 * @param msh Reference to the mesh
 * @param cl Refererence to the cell
 * @param degree Order of the quadrature
 * @return Quadrature points
 */
template<typename T>
std::vector<disk::quadrature_point<T, 2>>
integrate(const disk::simplicial_mesh<T, 2>&                msh,
          const typename disk::simplicial_mesh<T, 2>::cell& cl,
          const size_t                                      degree)
{
    const auto pts = points(msh, cl);
    assert(pts.size() == 3);
    return disk::quadrature::dunavant(degree, pts);
}

/**
 * @brief Return quadrature points on a face of a physical 2D triangle
 *
 * @param msh Reference to the mesh
 * @param fc Refererence to the face
 * @param degree Order of the quadrature
 * @return Quadrature points
 */
template<typename T>
std::vector<disk::quadrature_point<T, 2>>
integrate(const disk::simplicial_mesh<T, 2>&                msh,
          const typename disk::simplicial_mesh<T, 2>::face& fc,
          const size_t                                      degree)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 2);
    return disk::quadrature::gauss_legendre(degree, pts[0], pts[1]);
}

namespace priv
{

/**
 * @brief Map a point from the reference triangle to the physical triangle
 *
 * @tparam T scalar type
 * @param msh (simplicial) mesh
 * @param face specified 3D-face (triangle)
 * @param pm point in the reference triangle
 * @return point<T, 3> point in the physical triangle
 */
template<typename T>
point<T, 3>
map_to_physical(const disk::simplicial_mesh<T, 3>&                msh,
                const typename disk::simplicial_mesh<T, 3>::face& face,
                const point<T, 2>&                                pm)
{
    const auto pts = points(msh, face);
    // Compute the integration basis
    const auto [p0, v0, v1] = integration_basis(pts[0], pts[1], pts[2]);
    return p0 + v0 * pm.x() + v1 * pm.y();
}

/**
 * @brief Map a point from the reference tetrahedra to the physical tetrahedra
 *
 * @tparam T scalar type
 * @param msh (simplicial) mesh
 * @param cl specified 3D-cell (tetrahedra)
 * @param p point in the reference tetrahedra
 * @return point<T, 3> point in the physical tetrahedra
 */
template<typename T>
point<T, 3>
map_to_physical(const disk::simplicial_mesh<T, 3>&                msh,
                const typename disk::simplicial_mesh<T, 3>::cell& cl,
                const point<T, 3>&                                p)
{
    const auto pts = points(msh, cl);

    // Compute the integration basis
    const auto [p0, v0, v1, v2] = integration_basis(pts[0], pts[1], pts[2], pts[3]);

    return p0 + v0 * p.x() + v1 * p.y() + v2 * p.z();
}

} // namespace priv

/**
 * @brief Compute a quadrature of order "degree" in the physical space of the specified 3D-face (triangle)
 *
 * @tparam T scalar type
 * @param msh (simplicial) mesh
 * @param fc specified 3D-face (triangle)
 * @param degree order of the quadrature
 * @return std::vector<disk::quadrature_point<T, 2>> quadrature
 */
template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate(const disk::simplicial_mesh<T, 3>& msh, const typename disk::simplicial_mesh<T, 3>::face& fc, size_t degree)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 3);
    return disk::quadrature::dunavant(degree, pts[0], pts[1], pts[2]);
}

/**
 * @brief Compute a quadrature of order "degree" in the physical space of the specified 3D-cell (tetrahedra)
 *
 * @tparam T scalar type
 * @param msh (simplicial) mesh
 * @param cl specified 3D-cell (tetrahedra)
 * @param degree order of the quadrature
 * @return std::vector<disk::quadrature_point<T, 3>> quadrature
 */
template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate(const disk::simplicial_mesh<T, 3>& msh, const typename disk::simplicial_mesh<T, 3>::cell& cl, size_t degree)
{
    if (degree == 0)
    {
        return priv::integrate_degree0(msh, cl);
    }

    const auto m_quadrature_data = disk::tetrahedron_quadrature(degree);

    const auto pts = points(msh, cl);
    assert(pts.size() == 4);
    const auto meas = measure(msh, cl);

    // Compute the integration basis
    const auto [p0, v0, v1, v2] = integration_basis(pts[0], pts[1], pts[2], pts[3]);

    auto tr = [ p0 = p0, v0 = v0, v1 = v1, v2 = v2, meas ](const std::pair<point<T, 3>, T>& qd) -> auto
    {
        const auto point  = p0 + v0 * qd.first.x() + v1 * qd.first.y() + v2 * qd.first.z();
        const auto weight = qd.second * meas;
        return disk::make_qp(point, weight);
    };

    std::vector<disk::quadrature_point<T, 3>> ret(m_quadrature_data.size());
    std::transform(m_quadrature_data.begin(), m_quadrature_data.end(), ret.begin(), tr);

    return ret;
}

} // namespace disk

