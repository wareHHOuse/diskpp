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

namespace disk
{

// the fourth order tensor is stored like a Matrix
// | A1111  A1112  A1211  A1212 |
// | A1121  A1122  A1221  A1222 |
// | A2111  A2112  A2211  A2212 |
// | A2121  A2122  A2221  A2222 |

// Put it in Line
template<typename T, int DIM2> // DIM2 = DIM*DIM
static_matrix<T, DIM2, DIM2>
changeFormatRowTensor(const static_matrix<T, DIM2, DIM2>& tens)
{
    int DIM = 0;
    if (DIM2 == 1)
        DIM = 1;
    else if (DIM2 == 4)
        DIM = 2;
    else if (DIM2 == 9)
        DIM = 3;
    else
        assert(false);

    static_matrix<T, DIM2, DIM2> ret = static_matrix<T, DIM2, DIM2>::Zero();

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            const int row = i * DIM + j;
            for (int k = 0; k < DIM; k++)
            {
                for (int l = 0; l < DIM; l++)
                {
                    ret(row, k * DIM + l) = tens(i * DIM + k, j * DIM + l);
                }
            }
        }
    }

    //    for (size_t i = 0; i < DIM; i++)
    //       for (size_t j = 0; j < DIM; j++)
    //          ret.block(i*DIM + j, 0, 1, DIM2) = converttovector( tens.block(i*DIM, j*DIM, DIM,
    //          DIM);

    return ret;
}

template<typename T, int DIM2>
static_matrix<T, DIM2, DIM2>
changeFormatColTensor(const static_matrix<T, DIM2, DIM2>& tens)
{
    return changeFormatRowTensor(tens).transpose();
}

// Product Tensor - Matrix

template<typename T, int DIM>
static_matrix<T, DIM, DIM>
tm_prod(const static_tensor<T, DIM>& tens, const static_matrix<T, DIM, DIM>& mat)
{
    static_matrix<T, DIM, DIM> ret;

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            ret(i, j) = (tens.block(i * DIM, j * DIM, DIM, DIM).cwiseProduct(mat)).sum();
        }
    }
    return ret;
}

// optimization mat(row,col) neq 0, 0 else
template<typename T, int DIM>
static_matrix<T, DIM, DIM>
tm_prod(const static_tensor<T, DIM>& tens, const static_matrix<T, DIM, DIM>& mat, const int row, const int col)
{
    static_matrix<T, DIM, DIM> ret;
    ret.setConstant(mat(row, col));

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            ret(i, j) *= tens(i * DIM + row, j * DIM + col);
        }
    }

    return ret;
}

template<typename T>
T
tm_prod(const T& tens, const T& mat)
{
    return tens * mat;
}

// Compute Kronecker product

template<typename T, int M, int N, int P, int Q>
void
computeKroneckerProduct(const static_matrix<T, M, N>& A, const static_matrix<T, P, Q>& B)
{
    static_assert((M == N && N == P && P == Q), "Kronecker product : Not yet develloped");
}

// T_ijkl = A_ij B_kl

template<typename T, int DIM>
static_tensor<T, DIM>
computeKroneckerProduct(const static_matrix<T, DIM, DIM>& A, const static_matrix<T, DIM, DIM>& B)
{
    static_tensor<T, DIM> ret = static_tensor<T, DIM>::Zero();

    for (int j = 0; j < DIM; j++)
    {
        for (int i = 0; i < DIM; i++)
        {
            ret.block(i * DIM, j * DIM, DIM, DIM) = A(i, j) * B;
        }
    }
    return ret;
}

template<typename T, int DIM>
static_matrix<T, DIM, DIM>
computeKroneckerProduct(const static_vector<T, DIM>& A, const static_vector<T, DIM>& B)
{
    return A * B.transpose();
}

// contracted product
template<typename T, int DIM>
T
computeContractedProduct(const static_vector<T, DIM>& A, const static_vector<T, DIM>& B)
{
    return A.dot(B);
}

template<typename T, int DIM>
T
computeContractedProduct(const static_matrix<T, DIM, DIM>& A, const static_matrix<T, DIM, DIM>& B)
{
    return A.cwiseProduct(B).sum();
}

template<typename T, int DIM>
static_matrix<T, DIM, DIM>
computeContractedProduct(const static_tensor<T, DIM>& Tens, const static_matrix<T, DIM, DIM>& B)
{
    return tm_prod(Tens, B);
}

// T_ijkl = A_ik B_jl

template<typename T, int DIM>
static_tensor<T, DIM>
computeProductSup(const static_matrix<T, DIM, DIM>& A, const static_matrix<T, DIM, DIM>& B)
{
    static_tensor<T, DIM> ret = static_tensor<T, DIM>::Zero();

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                for (int l = 0; l < DIM; l++)
                {
                    ret(i * DIM + k, j * DIM + l) = A(i, k) * B(j, l);
                }
            }
        }
    }

    return ret;
}

// T_ijkl = A_il B_jk

template<typename T, int DIM>
static_tensor<T, DIM>
computeProductInf(const static_matrix<T, DIM, DIM>& A, const static_matrix<T, DIM, DIM>& B)
{
    static_tensor<T, DIM> ret = static_tensor<T, DIM>::Zero();

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                for (int l = 0; l < DIM; l++)
                {
                    ret(i * DIM + k, j * DIM + l) = A(i, l) * B(j, k);
                }
            }
        }
    }
    return ret;
}

template<typename T, int DIM>
static_tensor<T, DIM>
compute_IdentityTensor()
{
    static_tensor<T, DIM> ret = static_tensor<T, DIM>::Zero();
    T                     one = T{1};

    if (DIM == 1)
        ret(0, 0) = one; // I1111
    else if (DIM == 2)
    {
        ret(0, 0) = one; // I1111
        ret(0, 3) = one; // I1212
        ret(3, 0) = one; // I2121
        ret(3, 3) = one; // I2222
    }
    else if (DIM == 3)
    {
        ret(0, 0) = one; // I1111
        ret(0, 4) = one; // I1212
        ret(0, 8) = one; // I1313
        ret(4, 0) = one; // I2121
        ret(4, 4) = one; // I2222
        ret(4, 8) = one; // I2323
        ret(8, 0) = one; // I3131
        ret(8, 4) = one; // I3232
        ret(8, 8) = one; // I3333
    }
    else
        static_assert((DIM == 1 || DIM == 2 || DIM == 3), "Wrong dimension only 2 and 3");

    return ret;
}

template<typename T, int DIM>
static_tensor<T, DIM>
compute_IdentitySymTensor()
{
    static_tensor<T, DIM> ret  = static_tensor<T, DIM>::Zero();
    T                     one  = T{1};
    T                     half = one / T{2};

    if (DIM == 1)
        ret(0, 0) = one; // I1111
    else if (DIM == 2)
    {
        ret(0, 0) = one;  // I1111
        ret(0, 3) = half; // I1212
        ret(1, 2) = half; // I1221
        ret(2, 1) = half; // I2112
        ret(3, 0) = half; // I2121
        ret(3, 3) = one;  // I2222
    }
    else if (DIM == 3)
    {
        ret(0, 0) = one;  // I1111
        ret(0, 4) = half; // I1212
        ret(0, 8) = half; // I1331
        ret(1, 3) = half; // I1221
        ret(2, 6) = half; // I1331
        ret(3, 1) = half; // I2112
        ret(4, 0) = half; // I2121
        ret(4, 4) = one;  // I2222
        ret(4, 8) = half; // I2323
        ret(5, 7) = half; // I2332
        ret(6, 2) = half; // I3113
        ret(8, 0) = half; // I3131
        ret(7, 5) = half; // I3223
        ret(8, 4) = half; // I3232
        ret(8, 8) = one;  // I333
    }
    else
        static_assert((DIM == 1 || DIM == 2 || DIM == 3), "Wrong dimension only 2 and 3");

    return ret;
}

template<typename T, int DIM>
static_tensor<T, DIM>
compute_IxI()
{
    return static_tensor<T, DIM>::Identity();
}

template<typename T>
static_tensor<T, 3>
transpose(const static_tensor<T, 3>& tens)
{
    static_tensor<T, 3> ret;

    // block 11
    ret(0, 0) = tens(0, 0);
    ret(0, 1) = tens(0, 3);
    ret(0, 2) = tens(0, 6);
    ret(1, 0) = tens(3, 0);
    ret(1, 1) = tens(3, 3);
    ret(1, 2) = tens(3, 6);
    ret(2, 0) = tens(6, 0);
    ret(2, 1) = tens(6, 3);
    ret(2, 2) = tens(6, 6);

    // block 12
    ret(0, 3) = tens(0, 1);
    ret(0, 4) = tens(0, 4);
    ret(0, 5) = tens(0, 7);
    ret(1, 3) = tens(3, 1);
    ret(1, 4) = tens(3, 4);
    ret(1, 5) = tens(3, 7);
    ret(2, 3) = tens(6, 1);
    ret(2, 4) = tens(6, 4);
    ret(2, 5) = tens(6, 7);

    // block 13
    ret(0, 6) = tens(0, 2);
    ret(0, 7) = tens(0, 5);
    ret(0, 8) = tens(0, 8);
    ret(1, 6) = tens(3, 2);
    ret(1, 7) = tens(3, 5);
    ret(1, 8) = tens(3, 8);
    ret(2, 6) = tens(6, 2);
    ret(2, 7) = tens(6, 5);
    ret(2, 8) = tens(6, 0);

    // block 21
    ret(3, 0) = tens(1, 0);
    ret(3, 1) = tens(1, 3);
    ret(3, 2) = tens(1, 6);
    ret(4, 0) = tens(4, 0);
    ret(4, 1) = tens(4, 3);
    ret(4, 2) = tens(4, 6);
    ret(5, 0) = tens(7, 0);
    ret(5, 1) = tens(7, 3);
    ret(5, 2) = tens(7, 6);

    // block 22
    ret(3, 3) = tens(1, 1);
    ret(3, 4) = tens(1, 4);
    ret(3, 5) = tens(1, 7);
    ret(4, 3) = tens(4, 1);
    ret(4, 4) = tens(4, 4);
    ret(4, 5) = tens(4, 7);
    ret(5, 3) = tens(7, 1);
    ret(5, 4) = tens(7, 4);
    ret(5, 5) = tens(7, 7);

    // block 23
    ret(3, 6) = tens(1, 2);
    ret(3, 7) = tens(1, 5);
    ret(3, 8) = tens(1, 8);
    ret(4, 6) = tens(4, 2);
    ret(4, 7) = tens(4, 5);
    ret(4, 8) = tens(4, 8);
    ret(5, 6) = tens(7, 2);
    ret(5, 7) = tens(7, 5);
    ret(5, 8) = tens(7, 8);

    // block 31
    ret(6, 0) = tens(2, 0);
    ret(6, 1) = tens(2, 3);
    ret(6, 2) = tens(2, 6);
    ret(7, 0) = tens(5, 0);
    ret(7, 1) = tens(5, 3);
    ret(7, 2) = tens(5, 6);
    ret(8, 0) = tens(8, 0);
    ret(8, 1) = tens(8, 3);
    ret(8, 2) = tens(8, 6);

    // block 32
    ret(6, 3) = tens(2, 1);
    ret(6, 4) = tens(2, 4);
    ret(6, 5) = tens(2, 7);
    ret(7, 3) = tens(5, 1);
    ret(7, 4) = tens(5, 4);
    ret(7, 5) = tens(5, 7);
    ret(8, 3) = tens(8, 1);
    ret(8, 4) = tens(8, 4);
    ret(8, 5) = tens(8, 7);

    // block 33
    ret(6, 6) = tens(2, 2);
    ret(6, 7) = tens(2, 5);
    ret(6, 8) = tens(2, 8);
    ret(7, 6) = tens(5, 2);
    ret(7, 7) = tens(5, 5);
    ret(7, 8) = tens(5, 8);
    ret(8, 6) = tens(8, 2);
    ret(8, 7) = tens(8, 5);
    ret(8, 8) = tens(8, 8);

    return ret;
}

template<typename T>
static_tensor<T, 2>
transpose(const static_tensor<T, 2>& tens)
{
    static_tensor<T, 2> ret;

    // block 11
    ret(0, 0) = tens(0, 0);
    ret(0, 1) = tens(0, 2);
    ret(1, 0) = tens(2, 0);
    ret(1, 1) = tens(2, 2);

    // block 12
    ret(0, 2) = tens(0, 1);
    ret(0, 3) = tens(0, 3);
    ret(1, 2) = tens(2, 1);
    ret(1, 3) = tens(2, 3);

    // block 21
    ret(2, 0) = tens(1, 0);
    ret(2, 1) = tens(1, 2);
    ret(3, 0) = tens(3, 0);
    ret(3, 1) = tens(3, 2);

    // block 22
    ret(2, 2) = tens(1, 1);
    ret(2, 3) = tens(1, 3);
    ret(3, 2) = tens(3, 1);
    ret(3, 3) = tens(3, 3);

    return ret;
}

template<typename T, int DIM>
static_tensor<T, DIM>
symetric_part(const static_tensor<T, DIM>& tens)
{
    return (tens + transpose(tens)) / T(2);
}
}
