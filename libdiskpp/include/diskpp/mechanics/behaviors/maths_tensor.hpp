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

namespace disk
{

// See https://trilinos.org/docs/r11.14/packages/intrepid/doc/html/Intrepid__MiniTensor__Tensor4_8t_8h_source.html
// for details about operations

// the fourth order tensor is stored like a Matrix
// | A1111  A1112  A1211  A1212 |
// | A1121  A1122  A1221  A1222 |
// | A2111  A2112  A2211  A2212 |
// | A2121  A2122  A2221  A2222 |

// return Aijkl
template<typename T, int DIM>
T
coeff(const static_tensor<T, DIM>& A, const int i, const int j, const int k, const int l)
{
    if (i < 0 || j < 0 || k < 0 || l < 0 || i >= DIM || j >= DIM || k >= DIM || l >= DIM)
        throw std::invalid_argument("Invalid coefficient");

    return A(i * DIM + k, j * DIM + l);
}

template<typename T, int DIM>
void
coeff(static_tensor<T, DIM>& A, const int i, const int j, const int k, const int l, const T val)
{
    if (i < 0 || j < 0 || k < 0 || l < 0 || i >= DIM || j >= DIM || k >= DIM || l >= DIM)
        throw std::invalid_argument("Invalid coefficient");

    A(i * DIM + k, j * DIM + l) = val;
}

// Product Tensor - Matrix
// aij = Aijkl Bkl
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

// aij = Bkl Aklij
template<typename T, int DIM>
static_matrix<T, DIM, DIM>
tm_prod(const static_matrix<T, DIM, DIM>& mat, const static_tensor<T, DIM>& tens)
{
    static_matrix<T, DIM, DIM> ret = static_matrix<T, DIM, DIM>::Zero();

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                for (int l = 0; l < DIM; l++)
                {
                    ret(i, j) += mat(k, l) * coeff<T, DIM>(tens, k, l, i, j);
                }
            }
        }
    }

    return ret;
}

// aij = Aijkl Bkl
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
            ret(i, j) *= coeff<T, DIM>(tens, i, j, row, col);
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

// T_ijkl = A_ij B_kl

template<typename T, int DIM>
static_tensor<T, DIM>
Kronecker(const static_matrix<T, DIM, DIM>& A, const static_matrix<T, DIM, DIM>& B)
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
Kronecker(const static_vector<T, DIM>& A, const static_vector<T, DIM>& B)
{
    return A * B.transpose();
}

// contracted product
template<typename T, int DIM>
T
InnerProduct(const static_vector<T, DIM>& A, const static_vector<T, DIM>& B)
{
    return A.dot(B);
}

template<typename T, int DIM>
T
InnerProduct(const static_matrix<T, DIM, DIM>& A, const static_matrix<T, DIM, DIM>& B)
{
    return A.cwiseProduct(B).sum();
}

template<typename T, int DIM>
static_matrix<T, DIM, DIM>
ContractedProduct(const static_tensor<T, DIM>& Tens, const static_matrix<T, DIM, DIM>& B)
{
    return tm_prod(Tens, B);
}

template<typename T, int DIM>
static_matrix<T, DIM, DIM>
ContractedProduct(const static_matrix<T, DIM, DIM>& B, const static_tensor<T, DIM>& Tens)
{
    return tm_prod(B, Tens);
}

// Cijkl = Aijpq Bpqkl
template<typename T, int DIM>
static_tensor<T, DIM>
ContractedProduct(const static_tensor<T, DIM>& A, const static_tensor<T, DIM>& B)
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
                    for (int p = 0; p < DIM; p++)
                    {
                        for (int q = 0; q < DIM; q++)
                        {
                            ret(i * DIM + k, j * DIM + l) +=
                              coeff<T, DIM>(A, i, j, p, q) * coeff<T, DIM>(B, p, q, k, l);
                        }
                    }
                }
            }
        }
    }

    return ret;
}

// T_ijkl = A_ik B_jl

template<typename T, int DIM>
static_tensor<T, DIM>
ProductSup(const static_matrix<T, DIM, DIM>& A, const static_matrix<T, DIM, DIM>& B)
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
ProductInf(const static_matrix<T, DIM, DIM>& A, const static_matrix<T, DIM, DIM>& B)
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
IdentityTensor4()
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
IdentitySymTensor4()
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
IxI()
{
    return static_tensor<T, DIM>::Identity();
}

//(A^T)ijkl = Aklij
template<typename T, int DIM>
static_tensor<T, DIM>
transpose(const static_tensor<T, DIM>& tens)
{
    static_tensor<T, DIM> ret;

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                for (int l = 0; l < DIM; l++)
                {
                    ret(i * DIM + k, j * DIM + l) = coeff<T, DIM>(tens, k, l, i, j);
                }
            }
        }
    }

    return ret;
}

template<typename T, int DIM>
static_tensor<T, DIM>
Odot(const static_matrix<T, DIM, DIM>& A, const static_matrix<T, DIM, DIM>& B)
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
                    T val = (A(i, k) * B(j, l) + A(i, l) * B(j, k)) / T(2);
                    coeff<T, DIM>(ret, i, j, k, l, val);
                }
            }
        }
    }

    return ret;
}
}
