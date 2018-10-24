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

#include "common/eigen.hpp"
#include "mechanics/behaviors/maths_tensor.hpp"

namespace disk
{

// the fourth order tensor is stored like a Matrix
// | A1111  A1112  A1211  A1212 |
// | A1121  A1122  A1221  A1222 |
// | A2111  A2112  A2211  A2212 |
// | A2121  A2122  A2221  A2222 |

// Convert Matrix in Coulum vector

template<typename T, int DIM>
void
converttovector(const static_matrix<T, DIM, DIM>& mat)
{
    static_assert((DIM == 2 || DIM == 3), "Can not compute conversion for this dimension");
}

template<typename T>
static_vector<T, 4>
converttovector(const static_matrix<T, 2, 2>& mat)
{
    return static_vector<T, 4>{mat(0, 0), mat(1, 0), mat(0, 1), mat(1, 1)};
}

template<typename T>
static_vector<T, 9>
converttovector(const static_matrix<T, 3, 3>& mat)
{
    return static_vector<T, 9>{
      mat(0, 0), mat(1, 0), mat(2, 0), mat(0, 1), mat(1, 1), mat(2, 1), mat(0, 2), mat(1, 2), mat(2, 2)};
}

// Convert vector in matrix

template<typename T, int DIM>
void
converttomatrix(const static_vector<T, DIM>& vec)
{
    static_assert((DIM == 4 || DIM == 9), "Can not compute conversion for this dimension");
}

template<typename T>
static_matrix<T, 2, 2>
converttomatrix(const static_vector<T, 4>& vec)
{
    static_matrix<T, 2, 2> mat;

    mat(0, 0) = vec(0);
    mat(1, 0) = vec(1);

    mat(0, 1) = vec(2);
    mat(1, 1) = vec(3);

    return mat;
}

template<typename T>
static_matrix<T, 3, 3>
converttomatrix(const static_vector<T, 9>& vec)
{
    static_matrix<T, 3, 3> mat;

    mat(1, 0) = vec(0);
    mat(2, 0) = vec(1);
    mat(3, 0) = vec(2);

    mat(1, 1) = vec(3);
    mat(2, 1) = vec(4);
    mat(3, 1) = vec(5);

    mat(1, 2) = vec(6);
    mat(2, 2) = vec(7);
    mat(3, 2) = vec(8);

    return mat;
}

template<typename T>
static_matrix<T, 3, 3>
convertMatrix3D(const static_matrix<T, 3, 3>& mat)
{
    return mat;
}

template<typename T>
static_matrix<T, 3, 3>
convertMatrix3D(const static_matrix<T, 2, 2>& mat)
{
    static_matrix<T, 3, 3> ret = static_matrix<T, 3, 3>::Zero();
    ret.block(0, 0, 2, 2)      = mat;

    return ret;
}

template<typename T>
static_matrix<T, 3, 3>
convertMatrix3DwithOne(const static_matrix<T, 3, 3>& mat)
{
    return mat;
}

template<typename T>
static_matrix<T, 3, 3>
convertMatrix3DwithOne(const static_matrix<T, 2, 2>& mat)
{
    static_matrix<T, 3, 3> ret = static_matrix<T, 3, 3>::Identity();
    ret.block(0, 0, 2, 2)      = mat;

    return ret;
}

template<typename T, int DIM>
static_matrix<T, DIM, DIM>
convertMatrix(const static_matrix<T, 3, 3>& mat)
{
    static_assert((DIM == 2 || DIM == 3), "Can not compute conversion for this dimension");

    return mat.block(0, 0, DIM, DIM);
}

template<typename T, int DIM>
static_tensor<T, DIM>
convertTensor(const static_tensor<T, 3>& tens)
{
    static_assert((DIM == 2 || DIM == 3), "Can not compute conversion for this dimension");

    if (DIM == 2)
    {
        static_tensor<T, 3> ret = static_tensor<T, 3>::Zero();

        ret.block(0, 0, 2, 2) = tens.block(0, 0, 2, 2);
        ret.block(0, 2, 2, 2) = tens.block(0, 3, 2, 2);
        ret.block(2, 0, 2, 2) = tens.block(3, 0, 2, 2);
        ret.block(2, 2, 2, 2) = tens.block(3, 3, 2, 2);

        return ret.block(0, 0, DIM * DIM, DIM * DIM);
    }

    return tens.block(0, 0, DIM * DIM, DIM * DIM);
}

template<typename T, int DIM>
static_matrix<T, 6, 6>
convertTensorNotationMangel(const static_tensor<T, DIM>& tens)
{
    static_assert((DIM == 2 || DIM == 3), "Can not compute conversion for this dimension");

    static_matrix<T, 6, 6> ret = static_matrix<T, 6, 6>::Zero();

    ret(0, 0) = coeff<T, DIM>(tens, 0, 0, 0, 0);
    ret(0, 1) = coeff<T, DIM>(tens, 0, 0, 1, 1);
    ret(1, 0) = coeff<T, DIM>(tens, 1, 1, 0, 0);
    ret(1, 1) = coeff<T, DIM>(tens, 1, 1, 1, 1);

    ret(0, 3) = coeff<T, DIM>(tens, 0, 0, 0, 1) * sqrt(T(2));
    ret(1, 3) = coeff<T, DIM>(tens, 1, 1, 0, 1) * sqrt(T(2));
    ret(3, 0) = coeff<T, DIM>(tens, 0, 1, 0, 0) * sqrt(T(2));
    ret(3, 1) = coeff<T, DIM>(tens, 0, 1, 1, 1) * sqrt(T(2));
    ret(3, 3) = coeff<T, DIM>(tens, 0, 1, 0, 1) * T(2);

    if (DIM == 3)
    {
        ret(0, 2) = coeff<T, DIM>(tens, 0, 0, 2, 2);
        ret(1, 2) = coeff<T, DIM>(tens, 1, 1, 2, 2);
        ret(2, 2) = coeff<T, DIM>(tens, 2, 2, 2, 2);
        ret(2, 0) = coeff<T, DIM>(tens, 2, 2, 0, 0);
        ret(2, 1) = coeff<T, DIM>(tens, 2, 2, 1, 1);

        ret(2, 3) = coeff<T, DIM>(tens, 2, 2, 0, 1) * sqrt(T(2));
        ret(3, 2) = coeff<T, DIM>(tens, 0, 1, 2, 2) * sqrt(T(2));

        ret(0, 5) = coeff<T, DIM>(tens, 0, 0, 1, 2) * sqrt(T(2));
        ret(0, 4) = coeff<T, DIM>(tens, 0, 0, 0, 2) * sqrt(T(2));
        ret(1, 5) = coeff<T, DIM>(tens, 1, 1, 1, 2) * sqrt(T(2));
        ret(1, 4) = coeff<T, DIM>(tens, 1, 1, 0, 2) * sqrt(T(2));
        ret(2, 5) = coeff<T, DIM>(tens, 2, 2, 1, 2) * sqrt(T(2));
        ret(2, 4) = coeff<T, DIM>(tens, 2, 2, 0, 2) * sqrt(T(2));

        ret(4, 0) = coeff<T, DIM>(tens, 0, 2, 0, 0) * sqrt(T(2));
        ret(4, 1) = coeff<T, DIM>(tens, 0, 2, 1, 1) * sqrt(T(2));
        ret(4, 2) = coeff<T, DIM>(tens, 0, 2, 2, 2) * sqrt(T(2));
        ret(5, 0) = coeff<T, DIM>(tens, 1, 2, 0, 0) * sqrt(T(2));
        ret(5, 1) = coeff<T, DIM>(tens, 1, 2, 1, 1) * sqrt(T(2));
        ret(5, 2) = coeff<T, DIM>(tens, 1, 2, 2, 2) * sqrt(T(2));

        ret(3, 5) = coeff<T, DIM>(tens, 0, 1, 1, 2) * T(2);
        ret(3, 4) = coeff<T, DIM>(tens, 0, 1, 0, 2) * T(2);
        ret(4, 3) = coeff<T, DIM>(tens, 0, 2, 0, 1) * T(2);
        ret(4, 5) = coeff<T, DIM>(tens, 0, 2, 1, 2) * T(2);
        ret(4, 4) = coeff<T, DIM>(tens, 0, 2, 0, 2) * T(2);
        ret(5, 3) = coeff<T, DIM>(tens, 1, 2, 0, 1) * T(2);
        ret(5, 5) = coeff<T, DIM>(tens, 1, 2, 1, 2) * T(2);
        ret(5, 4) = coeff<T, DIM>(tens, 1, 2, 0, 2) * T(2);
    }

    return ret;
}

template<typename T>
static_tensor<T, 3>
convertCtoA(const static_tensor<T, 3>& C, const static_matrix<T, 3, 3>& PK2, const static_matrix<T, 3, 3>& F)
{
    static_tensor<T, 3> A = static_tensor<T, 3>::Zero();

    for (int i = 0; i < 3; i++)
    {
        for (int J = 0; J < 3; J++)
        {
            for (int k = 0; k < 3; k++)
            {
                for (int L = 0; L < 3; L++)
                {
                    T sum = T(0);
                    for (int M = 0; M < 3; M++)
                    {
                        T sum1 = T(0);
                        for (int N = 0; N < 3; N++)
                        {
                            sum1 += coeff<T, 3>(C, M, J, N, L) * F(k, N);
                        }
                        sum += sum1 * F(i, M);
                    }

                    if (i == k)
                    {
                        sum += PK2(J, L);
                    }

                    coeff<T, 3>(A, i, J, k, L, sum);
                    coeff<T, 3>(A, k, L, i, J, sum);
                }
            }
        }
    }

    return A;
}
}
