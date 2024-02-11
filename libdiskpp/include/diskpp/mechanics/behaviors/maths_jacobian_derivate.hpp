/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
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

#include "diskpp/common/eigen.hpp"
#include "mechanis/bevaviors/maths_tensor.hpp"

namespace disk
{

// the fourth order tensor is stored like a Matrix
// | A1111  A1112  A1211  A1212 |
// | A1121  A1122  A1221  A1222 |
// | A2111  A2112  A2211  A2212 |
// | A2121  A2122  A2221  A2222 |

// Compute \partial J / \partial M

template<typename T, int DIM>
void
JacobianFirstDerivate(const static_matrix<T, DIM, DIM>& M)
{
    static_assert((DIM == 2 || DIM == 3), "Can not compute jacobian derivate for this dimension");
}

template<typename T>
static_matrix<T, 2, 2>
JacobianFirstDerivate(const static_matrix<T, 2, 2>& M)
{
    static_matrix<T, 2, 2> ret = static_matrix<T, 2, 2>::Zero();

    ret(0, 0) = M(1, 1);
    ret(0, 1) = -M(1, 0);

    ret(1, 0) = -M(0, 1);
    ret(1, 1) = M(0, 0);

    return ret;
}

template<typename T>
static_matrix<T, 3, 3>
JacobianFirstDerivate(const static_matrix<T, 3, 3>& M)
{
    static_matrix<T, 3, 3> ret = static_matrix<T, 3, 3>::Zero();

    ret(0, 0) = M(1, 1) * M(2, 2) - M(2, 1) * M(1, 2);
    ret(0, 1) = M(1, 2) * M(2, 0) - M(2, 2) * M(1, 0);
    ret(0, 2) = M(1, 0) * M(2, 1) - M(1, 1) * M(2, 0);

    ret(1, 0) = M(0, 2) * M(2, 1) - M(2, 2) * M(0, 1);
    ret(1, 1) = M(0, 0) * M(2, 2) - M(0, 2) * M(2, 0);
    ret(1, 2) = M(2, 0) * M(0, 1) - M(0, 0) * M(2, 1);

    ret(2, 0) = M(0, 1) * M(1, 2) - M(1, 1) * M(0, 2);
    ret(2, 1) = M(0, 2) * M(1, 0) - M(0, 0) * M(1, 2);
    ret(2, 2) = M(0, 0) * M(1, 1) - M(0, 1) * M(1, 0);

    return ret;
}

// Compute \partial^2 J / \partial M^2

template<typename T, int DIM>
void
JacobianSecondDerivate(const static_matrix<T, DIM, DIM>& M)
{
    static_assert((DIM == 2 || DIM == 3), "Can not compute jacobian derivate for this dimension");
}

template<typename T>
static_tensor<T, 2>
JacobianSecondDerivate(const static_matrix<T, 2, 2>& M)
{
    static_tensor<T, 2> ret = static_tensor<T, 2>::Zero();
    T                   one = T(1);

    ret(1, 1) = one;
    ret(1, 2) = -one;
    ret(2, 1) = -one;
    ret(2, 2) = one;

    return ret;
}

template<typename T>
static_tensor<T, 3>
JacobianSecondDerivate(const static_matrix<T, 3, 3>& M)
{
    static_tensor<T, 3> ret = static_tensor<T, 3>::Zero();

    ret(1, 1) = M(2, 2);
    ret(1, 2) = -M(2, 1);
    ret(1, 3) = -M(2, 2);
    ret(2, 1) = -M(1, 2);
    ret(2, 2) = M(1, 1);
    ret(2, 3) = M(1, 2);
    ret(3, 1) = -M(2, 2);
    ret(3, 2) = M(2, 1);
    ret(3, 3) = M(2, 2);

    ret(1, 5) = M(2, 0);
    ret(1, 6) = M(2, 1);
    ret(1, 7) = -M(2, 0);
    ret(2, 5) = -M(1, 0);
    ret(2, 6) = -M(1, 1);
    ret(2, 7) = M(1, 0);
    ret(3, 5) = -M(2, 0);
    ret(3, 6) = -M(2, 1);
    ret(3, 7) = M(2, 0);

    ret(5, 1) = M(0, 2);
    ret(5, 2) = -M(0, 1);
    ret(5, 3) = -M(0, 2);
    ret(6, 1) = M(1, 2);
    ret(6, 2) = -M(1, 1);
    ret(6, 3) = -M(2, 1);
    ret(7, 1) = -M(0, 2);
    ret(7, 2) = M(0, 1);
    ret(7, 3) = M(0, 2);

    ret(5, 5) = M(0, 0);
    ret(5, 6) = M(0, 1);
    ret(5, 7) = -M(0, 0);
    ret(6, 5) = M(1, 0);
    ret(6, 6) = M(1, 1);
    ret(6, 7) = -M(1, 0);
    ret(7, 5) = -M(0, 0);
    ret(7, 6) = -M(0, 1);
    ret(7, 7) = M(0, 0);

    return ret;
}
}