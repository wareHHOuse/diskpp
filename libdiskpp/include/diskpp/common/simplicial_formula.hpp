/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018, 2019
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet (C) 2019                      nicolas.pignet@enpc.fr
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

#include <algorithm>
#include <tuple>

#include "diskpp/common/eigen.hpp"
#include "diskpp/mesh/point.hpp"

namespace disk
{

    /**
     * @brief Compute the area of a triangle by using the Kahan formula,
     * which minimize the round-off error
     *https://inria.hal.science/hal-00790071/document
     *
     * @param p0 first point of the triangle
     * @param p1 second point of the triangle
     * @param p2 third point of the triangle
     * @return T area
     */
    template <typename T, size_t N>
    T area_triangle_kahan(const point<T, N> &p0, const point<T, N> &p1, const point<T, N> &p2)
    {
        const T l10 = (p1 - p0).to_vector().norm();
        const T l20 = (p2 - p0).to_vector().norm();
        const T l21 = (p2 - p1).to_vector().norm();

        std::array<T, 3> length = {l10, l20, l21};

        std::sort(length.begin(), length.end());

        const T a = length[2];
        const T b = length[1];
        const T c = length[0];

        // Minimal assumption on geometry
        if (a > (b + c))
        {
            throw std::runtime_error("Area is negative");
        }

        return T(0.25) * std::sqrt((((a + (b + c)) * (a + (b - c))) * (c + (a - b))) * (c - (a - b)));
    }

/**
 * @brief Compute the integration basis in order to map a point from the reference to physcal frame (for a triangle)
 *
 * We search the two vectors generated form the edges of the triangle which are the most orthogonal.
 *
 * @tparam T scalar type
 * @tparam N dimension
 * @param p0 first point of the triangle
 * @param p1 second point of the triangle
 * @param p2 third point of the triangle
 * @return std::tuple<point<T, N>, point<T, N>, point<T, N>> the first point is the origin of the basis, the second and the third point are the basis
 */

template<typename T, size_t N>
std::tuple<point<T, N>, point<T, N>, point<T, N>>
integration_basis(const point<T, N>& p0, const point<T, N>& p1, const point<T, N>& p2)
{
    const std::array<point<T, N>, 3> pts   = {p0, p1, p2};
    size_t                           node  = 0;
    T                                pscal = T(1);

    for (size_t i = 0; i < 3; i++)
    {
        const auto v0 = (pts[(i + 1) % 3] - pts[i]).to_vector();
        const auto v1 = (pts[(i + 2) % 3] - pts[i]).to_vector();

        const auto v0n = v0 / v0.norm();
        const auto v1n = v1 / v1.norm();

        if (v0n.dot(v1n) < pscal) // we want the maximum angle;
        {
            node  = i;
            pscal = v0n.dot(v1n);
            //std::cout << "node: " << node << ", pscal: " << pscal << std::endl;
        }
    }

    const point<T, N> pbasis = pts[node];
    const point<T, N> b0     = pts[(node + 1) % 3] - pts[node];
    const point<T, N> b1     = pts[(node + 2) % 3] - pts[node];

    return std::make_tuple(pbasis, b0, b1);
}

/**
 * @brief Compute the volume of a tetrahedron by using the Kahan formula,
 * which minimize the round-off error
 * https://people.eecs.berkeley.edu/~wkahan/VtetLang.pdf
 *
 * @tparam T scalar type
 * @param p0 first point of the tetrahedron
 * @param p1 second point of the tetrahedron
 * @param p2 third point of the tetrahedron
 * @param p3 fourth point of the tetrahedron
 *
 * @return T volume
 */
template<typename T>
T
volume_tetrahedron_kahan(const point<T, 3>& p0, const point<T, 3>& p1, const point<T, 3>& p2, const point<T, 3>& p3)
{
    // facial difference = (u-v+w)
    auto fd = [](const T &u, const T &v, const T &w)
    { return (std::max(u, w) - v) + std::min(u, w); };

    // lengths of the edges
    const auto l10 = (p1 - p0).to_vector().norm();
    const auto l20 = (p2 - p0).to_vector().norm();
    const auto l30 = (p3 - p0).to_vector().norm();

    const auto l32 = (p3 - p2).to_vector().norm();
    const auto l31 = (p3 - p1).to_vector().norm();
    const auto l21 = (p2 - p1).to_vector().norm();

    // sort edges
    std::array<std::pair<T, T>, 3> length = {
        std::make_pair(l10, l32),
        std::make_pair(l20, l31),
        std::make_pair(l30, l21)};

    std::sort(length.begin(), length.end(), [&](const std::pair<T, T> &va, const std::pair<T, T> &vb)
              { return (va.first + va.second < vb.first + vb.second); });

    const auto [u, U] = length[2];
    const auto [v, V] = length[1];
    const auto [w, W] = length[0];

    // Accurate products
    const auto X = fd(w, U, v) * (U + v + w);
    const auto Y = fd(u, V, w) * (V + w + u);
    const auto Z = fd(v, W, u) * (W + u + v);

    const auto x = fd(U, v, w) * fd(v, w, U);
    const auto y = fd(V, w, u) * fd(w, u, V);
    const auto z = fd(W, u, v) * fd(u, v, W);

    // elementary factors
    const auto xi = std::sqrt(x * Y * Z);
    const auto eta = std::sqrt(y * Z * X);
    const auto zeta = std::sqrt(z * X * Y);
    const auto lambda = std::sqrt(x * y * z);

    const auto det = T(192.) * u * v * w;

    return std::sqrt((xi + eta + zeta - lambda) * (lambda + xi + eta - zeta) * (eta + zeta + lambda - xi) * (zeta + lambda + xi - eta)) / det;
}

/**
 * @brief Compute the integration basis in order to map a point from the reference to physcal frame (for a tetrahedron)
 *
 *
 * @tparam T scalar type
 * @tparam N dimension
 * @param p0 first point of the tetrahedron
 * @param p1 second point of the tetrahedron
 * @param p2 third point of the tetrahedron
 * @param p3 fourth point of the tetrahedron
 *
 * @return std::tuple<point<T, N>, point<T, N>, point<T, N>> the first point is the origin of the basis, the second and
 * the third point are the basis
 */
template<typename T>
std::tuple<point<T, 3>, point<T, 3>, point<T, 3>, point<T, 3>>
integration_basis(const point<T, 3>& p0, const point<T, 3>& p1, const point<T, 3>& p2, const point<T, 3> p3)
{
    return std::make_tuple(p0, p1 - p0, p2 - p0, p3 - p0);
}
}