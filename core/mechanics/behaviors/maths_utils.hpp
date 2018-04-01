/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
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

#pragma once

#include "mechanics/behaviors/maths_tensor.hpp"
#include "common/eigen.hpp"

namespace disk {

// deviatoric part of mat
template<typename T, int DIM>
static_matrix<T, DIM, DIM>
deviator(const static_matrix<T, DIM, DIM> mat)
{
   return mat - mat.trace() / T(DIM) * static_matrix<T, DIM, DIM>::Identity();
}

// spheric part of mat
template<typename T, int DIM>
static_matrix<T, DIM, DIM>
spheric(const static_matrix<T, DIM, DIM> mat)
{
   return mat.trace() / T(DIM) * static_matrix<T, DIM, DIM>::Identity();
}

// compute the first derivate of the Frobenius-norm of the deviatoric part of a matrix
template<typename T, int DIM>
static_matrix<T, DIM, DIM>
computeNormFroFirstDerivate(const static_matrix<T, DIM, DIM>& dev)
{
   const T normFrodev    = dev.norm();
   const T normFrodev3_2 = normFrodev * std::sqrt(normFrodev);

   auto mat = dev;
   mat.diagonal() *= (1 - 1.0 / T(DIM));

   return -mat / normFrodev3_2;
}

} // end disk