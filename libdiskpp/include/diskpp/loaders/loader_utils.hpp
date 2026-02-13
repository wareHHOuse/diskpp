/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2024
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */
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

#pragma once

#include <thread>
#include <vector>

#include "diskpp/loaders/strtot.hpp"

namespace disk
{

template<typename T>
std::optional<size_t>
boundary_number(const Eigen::Matrix<T, 2, 1>& n)
{
    /* Give a number to the boundaries of the unit
     * cube [0,1]^3 depending on the normal */

    Eigen::Matrix<T, 2, 1> x  = {1.0, 0.0};
    auto                   xn = n.dot(x);
    if (xn > -1.01 && xn < -0.99) // Plane x = 0
        return 3;

    if (xn > 0.99 && xn < 1.01) // Plane x = 1
        return 1;

    Eigen::Matrix<T, 2, 1> y  = {0.0, 1.0};
    auto                   yn = n.dot(y);
    if (yn > -1.01 && yn < -0.99) // Plane y = 0
        return 0;

    if (yn > 0.99 && yn < 1.01) // Plane y = 1
        return 2;

    return {};
}

template<typename T>
std::optional<size_t>
boundary_number(const Eigen::Matrix<T, 3, 1>& n)
{
    /* Give a number to the boundaries of the unit
     * cube [0,1]^3 depending on the normal */

    Eigen::Matrix<T, 3, 1> x  = {1.0, 0.0, 0.0};
    auto                   xn = n.dot(x);
    if (xn > -1.01 && xn < -0.99) // Plane x = 0
        return 2;

    if (xn > 0.99 && xn < 1.01) // Plane x = 1
        return 5;

    Eigen::Matrix<T, 3, 1> y  = {0.0, 1.0, 0.0};
    auto                   yn = n.dot(y);
    if (yn > -1.01 && yn < -0.99) // Plane y = 0
        return 1;

    if (yn > 0.99 && yn < 1.01) // Plane y = 1
        return 4;

    Eigen::Matrix<T, 3, 1> z  = {0.0, 0.0, 1.0};
    auto                   zn = n.dot(z);
    if (zn > -1.01 && zn < -0.99) // Plane z = 0
        return 0;

    if (zn > 0.99 && zn < 1.01) // Plane z = 1
        return 3;

    return {};
}

/* Some meshes from standard benchmarks (FVCA5, FVCA6) do not carry
 * boundary info. As these meshes usually represent the unit cube, try
 * to renumber the boundaries using the normals. */
template<typename Mesh>
bool
renumber_hypercube_boundaries(Mesh& msh)
{
    auto cvf     = connectivity_via_faces(msh);
    auto storage = msh.backend_storage();

    for (const auto& cl : msh)
    {
        const auto fcs = faces(msh, cl);
        for (const auto& fc : fcs)
        {
            auto neigh = cvf.neighbour_via(msh, cl, fc);
            if (neigh)
                continue; /* Not an external boundary */

            auto bn = boundary_number(normal(msh, cl, fc));
            if (!bn)
                return false; /* Can't renumber, unexpected normal. */

            auto& bi = storage->boundary_info[offset(msh, fc)];
            bi.is_boundary(true);
            bi.is_internal(false);
            bi.id(bn.value());
        }
    }

    return true;
}

inline bool
expect(std::ifstream& ifs, const std::string& str)
{
    std::string keyword;
    ifs >> keyword;
    if (keyword != str)
    {
        std::cout << "Expected keyword \"" << str << "\"" << std::endl;
        std::cout << "Found \"" << keyword << "\"" << std::endl;
        return false;
    }

    return true;
}

namespace priv
{

template<typename T>
void
sort_uniq(std::vector<T>& v)
{
    std::sort(v.begin(), v.end());
    auto uniq_iter = std::unique(v.begin(), v.end());
    v.erase(uniq_iter, v.end());
}

template<typename T>
point<T, 2>
read_2d_point_line(const char* str, char** endptr, T scalefactor)
{
    T t1, t2;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);

    return point<T, 2>{t1 * scalefactor, t2 * scalefactor};
}

template<typename T>
point<T, 3>
read_3d_point_line(const char* str, char** endptr, T scalefactor)
{
    T t1, t2, t3;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);

    return point<T, 3>{t1 * scalefactor, t2 * scalefactor, t3 * scalefactor};
}

} // namespace priv

} // namespace disk

#define THREADED
#ifdef THREADED
#define THREAD(name, body) std::thread name([&] {  body })
#define WAIT_THREAD(name) name.join()
#else
#define THREAD(name, body)                                                                                             \
    {                                                                                                                  \
        body                                                                                                           \
    }
#define WAIT_THREAD(name)
#endif
