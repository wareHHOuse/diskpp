/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#pragma once

#include "common/eigen.hpp"
#include "quadratures/quadratures.hpp"
#include "bases/bases.hpp"

namespace disk {

template<typename T>
T
mm_prod(const T& a, const T& b)
{
    return a*b;
}

template<typename T, int N>
T
mm_prod(const static_vector<T,N>& a, const static_vector<T,N>& b)
{
    return a.dot(b);
}

template<typename T, int N>
static_vector<T,N>
mm_prod(const static_vector<T,N>& a, const T& b)
{
    return a * b;
}

template<typename T, int N>
static_vector<T,N>
mm_prod(const T& a, const static_vector<T,N>& b)
{
    return a * b;
}

template<typename T, int N>
static_vector<T,N>
mm_prod(const static_matrix<T,N,N>& a, const static_vector<T,N>& b)
{
    return a * b;
}

template<typename T, int N>
T
mm_prod(const static_matrix<T,N,N>& a, const static_matrix<T,N,N>& b)
{
    return a.cwiseProduct(b).sum();
}

} // namespace disk
