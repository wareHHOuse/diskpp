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
 * Karol Cascavita (C) 2018                     klcascavitam@unal.edu.co
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

#include "common/eigen.hpp"
#include "quadratures/quadratures.hpp"
#include "bases/bases.hpp"

namespace disk {

namespace priv
{
template<typename T, int N>
Matrix<T, Dynamic, Dynamic>
outer_product(const Matrix<T, Dynamic, N>& a, const Matrix<T, Dynamic, N>& b)
{
    return b * a.transpose();
}

template<typename T, int N>
Matrix<T, Dynamic, N>
outer_product(const eigen_compatible_stdvector<Matrix<T, N, N>>& a, const Matrix<T, N, 1>& b)
{
    Matrix<T, Dynamic, N> ret(a.size(), N);
    for (size_t i = 0; i < a.size(); i++)
    {
        Matrix<T, N, 1> t = a[i] * b;
        ret.row(i)        = t.transpose();
    }

    return ret;
}

template<typename T, int N>
Matrix<T, Dynamic, Dynamic>
outer_product(const eigen_compatible_stdvector<Matrix<T, N, N>>& a,
              const eigen_compatible_stdvector<Matrix<T, N, N>>& b)
{
    Matrix<T, Dynamic, Dynamic> ret(a.size(), b.size());

    for (size_t i = 0; i < a.size(); i++)
        for (size_t j = 0; j < b.size(); j++)
            ret(i, j) = a[i].cwiseProduct(b[j]).sum();

    return ret;
}

template<typename T, int N>
Matrix<T, Dynamic, Dynamic>
outer_product(const eigen_compatible_stdvector<Matrix<T, N, N>>& a, const Matrix<T, N, N>& b)
{
    Matrix<T, Dynamic, 1> ret(a.size());

    for (size_t i = 0; i < a.size(); i++)
        ret(i) = a[i].cwiseProduct(b).sum();

    return ret;
}

template<typename T>
Matrix<T, Dynamic, 1>
inner_product(const T& a, const Matrix<T, Dynamic, 1>& b)
{
    return a * b;
}

template<typename T, int N>
Matrix<T, N, N>
inner_product(const T& a, const Matrix<T, N, N>& b)
{
    return a * b;
}

template<typename T, int N>
Matrix<T, N, N>
inner_product(const Matrix<T, N, N>& b, const T& a)
{
    return a * b;
}

template<typename T, int N>
Matrix<T, N, 1>
inner_product(const Matrix<T, N, 1>& b, const T& a)
{
    return a * b;
}

template<typename T, int N>
Matrix<T, N, 1>
inner_product(const T& a, const Matrix<T, N, 1>& b)
{
    return a * b;
}

template<typename T, int N>
T
inner_product(const Matrix<T, N, 1>& a, const Matrix<T, N, 1>& b)
{
    return a.dot(b);
}

template<typename T, int N>
T
inner_product(const Matrix<T, N, N>& b, const Matrix<T, N, N>& a)
{
    return a.cwiseProduct(b).sum();
}

template<typename T, int N>
Matrix<T, Dynamic, 1>
inner_product(const Matrix<T, N, 1>& a, const Matrix<T, Dynamic, N>& b)
{
    return b * a;
}

} // priv

template<typename T>
size_t howmany_dofs(const T& basis)
{
    return basis.size();
}

template<typename T>
size_t howmany_dofs(const T& basis, size_t min_degree, size_t max_degree)
{
    throw std::invalid_argument("error not used anymore");
    return 0;
}

} // namespace disk
