/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
 * Karol Cascavita (C) 2018                     klcascavitam@unal.edu.co
 * Nicolas Pignet (C) 2018                      nicolas.pignet@enpc.fr
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

#include "common/eigen.hpp"

namespace disk
{

namespace mechanics
{

// G: gradient

// Compute F = G + I

template<typename T, int DIM>
static_matrix<T, DIM, DIM>
convertGtoF(const static_matrix<T, DIM, DIM>& Gradient)
{
    return Gradient + static_matrix<T, DIM, DIM>::Identity();
}

template<typename T, int DIM>
static_matrix<T, DIM, DIM>
convertFtoG(const static_matrix<T, DIM, DIM>& F)
{
    return F - static_matrix<T, DIM, DIM>::Identity();
}

// Compute C = F^T * F

template<typename T, int DIM>
static_matrix<T, DIM, DIM>
convertFtoCauchyGreenRight(const static_matrix<T, DIM, DIM>& F)
{
    return F.transpose() * F;
}

template<typename T, int DIM>
static_matrix<T, DIM, DIM>
convertGtoCauchyGreenRight(const static_matrix<T, DIM, DIM>& G)
{
    const auto F = convertGtoF(G);
    return convertFtoCauchyGreenRight(F);
}

// Compute b = F * F^t

template<typename T, int DIM>
static_matrix<T, DIM, DIM>
convertFtoCauchyGreenLeft(const static_matrix<T, DIM, DIM>& F)
{
    return F * F.transpose();
}

template<typename T, int DIM>
static_matrix<T, DIM, DIM>
convertGtoCauchyGreenLeft(const static_matrix<T, DIM, DIM>& G)
{
    const auto F = convertGtoF(G);
    return convertFtoCauchyGreenLeft(F);
}

// Compute E = 1/2 *( C - I)

template<typename T, int DIM>
static_matrix<T, DIM, DIM>
convertCauchyGreenRighttoGreenLagrange(const static_matrix<T, DIM, DIM>& CauchyGreenRight)
{
    return T(0.5) * (CauchyGreenRight - static_matrix<T, DIM, DIM>::Identity());
}

template<typename T, int DIM>
static_matrix<T, DIM, DIM>
convertFtoGreenLagrange(const static_matrix<T, DIM, DIM>& F)
{
    const auto CauchyGreenRight = convertFtoCauchyGreenRight(F);
    return convertCauchyGreenRighttoGreenLagrange(CauchyGreenRight);
}

} // end mechanics

} // end disk