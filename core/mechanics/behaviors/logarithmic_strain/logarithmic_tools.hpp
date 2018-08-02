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
#include <array>
#include <tuple>

/* For details see the paper:
 *   Anisotropic additive plasticity in the logaritmic strain: modular
 *   kinematic formulation and implementation based on incremental minimization
 *   principles for standard materials
 *   C. Miehe, 3. Apel, M. Lambrecht
 *   Comput. Methods Appl. Mech. Engrg. (2002)
 */

namespace disk
{

namespace mechanics
{

// compute eigenvalues and eigenvectors
// ! Be careful Mat hes to be symmetric
template<typename T, int N>
std::pair<static_vector<T, N>, static_matrix<T, N, N>>
compute_eigenvalues(const static_matrix<T, N, N>& Mat)
{
    typedef static_matrix<T, N, N> matrix_type;

    SelfAdjointEigenSolver<matrix_type> es(Mat);

    return std::make_pair(es.eigenvalues(), es.eigenvectors());
}

// compute ei
template<typename T>
static_vector<T, 3>
compute_ei(const static_vector<T, 3>& lambda_i)
{
    static_vector<T, 3> ei;

    for (int i = 0; i < 3; i++)
    {
        ei(i) = std::log(lambda_i(i)) / T(2);
    }
    return ei;
}

// compute Elog Elog = P *log(D) * P^T
template<typename T>
static_matrix<T, 3, 3>
compute_Elog(const static_vector<T, 3>& lambda_i, const static_matrix<T, 3, 3>& P)
{
    const static_vector<T, 3>         ei = compute_ei(lambda_i);
    const Eigen::DiagonalMatrix<T, 3> D  = Eigen::DiagonalMatrix<T, 3>(ei[0], ei[1], ei[2]);

    return P * D * P.transpose();
}

// compute coefficient for logarithmic strain

// compute Ni
template<typename T>
std::array<static_vector<T, 3>, 3>
compute_Ni(const static_matrix<T, 3, 3>& evec)
{
    std::array<static_vector<T, 3>, 3> Ni;

    for (int i = 0; i < 3; i++)
    {
        Ni[i] = evec.col(i);
    }

    return Ni;
}

// compute ni = F * Ni
template<typename T>
std::array<static_vector<T, 3>, 3>
compute_ni(const static_matrix<T, 3, 3>& F, const std::array<static_vector<T, 3>, 3>& Ni)
{
    std::array<static_vector<T, 3>, 3> ni;

    for (int i = 0; i < 3; i++)
    {
        ni[i] = F * Ni[i];
    }

    return ni;
}

// compute di
template<typename T>
static_vector<T, 3>
compute_di(const static_vector<T, 3>& lambda_i)
{
    static_vector<T, 3> di;

    for (int i = 0; i < 3; i++)
    {
        di(i) = T(1) / lambda_i(i);
    }
    return di;
}

// compute fi
template<typename T>
static_vector<T, 3>
compute_fi(const static_vector<T, 3>& lambda_i)
{
    static_vector<T, 3> fi;

    for (int i = 0; i < 3; i++)
    {
        fi(i) = -T(2) / (lambda_i(i) * lambda_i(i));
    }
    return fi;
}

// compute zeta
template<typename T>
static_matrix<T, 3, 3>
compute_zeta(const static_matrix<T, 3, 3>& stress_T, const std::array<static_vector<T, 3>, 3>& Ni)
{
    static_matrix<T, 3, 3> zeta = static_matrix<T, 3, 3>::Zero();

    for (int j = 0; j < 3; j++)
    {
        for (int i = 0; i < 3; i++)
        {
            zeta(i, j) = computeInnerProduct(stress_T, computeKroneckerProduct(Ni[i], Ni[j]));
        }
    }
    return zeta;
}

// compute operator tensor
// compute Mij
template<typename T>
static_matrix<T, 3, 3>
compute_Mij(const std::array<static_vector<T, 3>, 3>& Ni,
            const std::array<static_vector<T, 3>, 3>& ni,
            const int                                 i,
            const int                                 j)
{
    static_matrix<T, 3, 3> Mij = static_matrix<T, 3, 3>::Zero();

    for (int a = 0; a < 3; a++)
    {
        for (int b = 0; b < 3; b++)
        {
            Mij(a, b) = ni[i](a) * Ni[j](b) + ni[j](a) * Ni[i](b);
        }
    }
    return Mij;
}

enum EigenCase
{
    THREE_EQUAL,
    TWO_EQUAL01,
    TWO_EQUAL02,
    TWO_EQUAL12,
    THREE_DIFF
};

namespace priv
{
template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
almost_equal(T x, T y, int ulp = 2)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x - y) <= std::numeric_limits<T>::epsilon() * std::abs(x + y) * ulp
           // unless the result is subnormal
           || std::abs(x - y) < std::numeric_limits<T>::min();
}

template<typename T>
int
selectCase(const static_vector<T, 3>& lambda_i)
{
    if (almost_equal(lambda_i(0), lambda_i(1)))
    {
        if (almost_equal(lambda_i(1), lambda_i(2)))
        {
            return THREE_EQUAL;
        }
        else
        {
            return TWO_EQUAL01;
        }
    }
    else
    {
        if (almost_equal(lambda_i(1), lambda_i(2)))
        {
            return TWO_EQUAL12;
        }
        else if (almost_equal(lambda_i(0), lambda_i(2)))
        {
            return TWO_EQUAL02;
        }
        else
        {
            return THREE_DIFF;
        }
    }

    return -1;
}

// lambda_0 != lambda_1 != lambda_2 && lambda_0 < lambda_1 < lambda_2
template<typename T>
std::tuple<static_matrix<T, 3, 3>, static_matrix<T, 3, 3>, T>
compute_three_diff(const static_vector<T, 3>& lambda_i, const static_vector<T, 3>& ei, const static_vector<T, 3>& di)
{
    T                      eta   = T(0);
    static_matrix<T, 3, 3> theta = static_matrix<T, 3, 3>::Zero();
    static_matrix<T, 3, 3> xi    = static_matrix<T, 3, 3>::Zero();

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (j != i)
            {
                const T lilj = lambda_i(i) - lambda_i(j);
                theta(i, j)  = (ei(i) - ei(j)) / lilj;
                xi(i, j)     = (theta(i, j) - di(j) / T(2)) / lilj;

                for (int k = 0; k < 3; k++)
                {
                    if (k != i && k != j)
                    {
                        const T lilk = lambda_i(i) - lambda_i(k);
                        eta += ei(i) / (T(2) * lilj * lilk);
                    }
                }
            }
        }
    }

    return std::make_tuple(theta, xi, eta);
}

// lambda_0 == lambda_1 == lambda_2
template<typename T>
std::tuple<static_matrix<T, 3, 3>, static_matrix<T, 3, 3>, T>
compute_three_equal(const static_vector<T, 3>& di, const static_vector<T, 3>& fi)
{
    const T d = di(0) / T(2);
    const T f = fi(0) / T(8);

    const T                eta   = f;
    static_matrix<T, 3, 3> theta = static_matrix<T, 3, 3>::Zero();
    static_matrix<T, 3, 3> xi    = static_matrix<T, 3, 3>::Zero();

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (j != i)
            {
                theta(i, j) = d;
                xi(i, j)    = f;
            }
        }
    }

    return std::make_tuple(theta, xi, eta);
}

// lambda_0 == lambda_1 != lambda_2
template<typename T>
std::tuple<static_matrix<T, 3, 3>, static_matrix<T, 3, 3>, T>
compute_two_equal(const static_vector<T, 3>& lambda_i,
                  const static_vector<T, 3>& ei,
                  const static_vector<T, 3>& di,
                  const static_vector<T, 3>& fi,
                  const int                  CASE)
{
    T                      eta   = T(0);
    static_matrix<T, 3, 3> theta = static_matrix<T, 3, 3>::Zero();
    static_matrix<T, 3, 3> xi    = static_matrix<T, 3, 3>::Zero();

    switch (CASE)
    {
        case TWO_EQUAL01:
        {
            theta(0, 1) = di(0) / T(2);
            theta(1, 0) = theta(0, 1);

            theta(0, 2) = (ei(0) - ei(2)) / (lambda_i(0) - lambda_i(2));
            theta(2, 0) = theta(0, 2);

            theta(1, 2) = (ei(1) - ei(2)) / (lambda_i(1) - lambda_i(2));
            theta(2, 1) = theta(1, 2);

            xi(0, 1) = fi(0) / T(8);
            xi(1, 0) = xi(0, 1);

            xi(0, 2) = (theta(0, 2) - di(2) / T(2)) / (lambda_i(0) - lambda_i(2));
            xi(2, 0) = (theta(2, 0) - di(0) / T(2)) / (lambda_i(2) - lambda_i(0));

            xi(1, 2) = (theta(1, 2) - di(2) / T(2)) / (lambda_i(1) - lambda_i(2));
            xi(2, 1) = (theta(2, 1) - di(1) / T(2)) / (lambda_i(2) - lambda_i(1));

            eta = xi(2, 0);

            break;
        }
        case TWO_EQUAL02:
        {
            theta(0, 1) = (ei(0) - ei(1)) / (lambda_i(0) - lambda_i(1));
            theta(1, 0) = theta(0, 1);

            theta(0, 2) = di(0) / T(2);
            theta(2, 0) = theta(0, 2);

            theta(1, 2) = (ei(1) - ei(2)) / (lambda_i(1) - lambda_i(2));
            theta(2, 1) = theta(1, 2);

            xi(0, 1) = (theta(0, 1) - di(1) / T(2)) / (lambda_i(0) - lambda_i(1));
            xi(1, 0) = (theta(1, 0) - di(0) / T(2)) / (lambda_i(1) - lambda_i(0));

            xi(0, 2) = fi(0) / T(8);
            xi(2, 0) = xi(0, 2);

            xi(1, 2) = (theta(1, 2) - di(2) / T(2)) / (lambda_i(1) - lambda_i(2));
            xi(2, 1) = (theta(2, 1) - di(1) / T(2)) / (lambda_i(2) - lambda_i(1));

            eta = xi(0, 1);
            break;
        }
        case TWO_EQUAL12:
        {
            theta(0, 1) = (ei(0) - ei(1)) / (lambda_i(0) - lambda_i(1));
            theta(1, 0) = theta(0, 1);

            theta(0, 2) = (ei(0) - ei(2)) / (lambda_i(0) - lambda_i(2));
            theta(2, 0) = theta(0, 2);

            theta(1, 2) = di(0) / T(2);
            theta(2, 1) = theta(1, 2);

            xi(0, 1) = (theta(0, 1) - di(1) / T(2)) / (lambda_i(0) - lambda_i(1));
            xi(1, 0) = (theta(1, 0) - di(0) / T(2)) / (lambda_i(1) - lambda_i(0));

            xi(0, 2) = (theta(0, 2) - di(2) / T(2)) / (lambda_i(0) - lambda_i(2));
            xi(2, 0) = (theta(2, 0) - di(0) / T(2)) / (lambda_i(2) - lambda_i(0));

            xi(1, 2) = fi(0) / T(8);
            xi(2, 1) = xi(1, 2);

            eta = xi(1, 2);
            break;
        }
        default: std::invalid_argument("LogarithmicStrain: Wrong case");
    }

    return std::make_tuple(theta, xi, eta);
}
}

template<typename T>
std::tuple<static_matrix<T, 3, 3>, static_matrix<T, 3, 3>, T>
compute_coefficient(const static_vector<T, 3>& lambda_i,
                    const static_vector<T, 3>& ei,
                    const static_vector<T, 3>& di,
                    const static_vector<T, 3>& fi)
{
    const int evcase = priv::selectCase(lambda_i);
    // std::cout << "CASE: " << evcase << std::endl;
    switch (evcase)
    {
        case THREE_DIFF: return priv::compute_three_diff(lambda_i, ei, di); break;
        case TWO_EQUAL01: return priv::compute_two_equal(lambda_i, ei, di, fi, TWO_EQUAL01); break;
        case TWO_EQUAL02: return priv::compute_two_equal(lambda_i, ei, di, fi, TWO_EQUAL02); break;
        case TWO_EQUAL12: return priv::compute_two_equal(lambda_i, ei, di, fi, TWO_EQUAL12); break;
        case THREE_EQUAL: return priv::compute_three_equal(di, fi); break;
        default: std::invalid_argument("LogarithmicStrain: Wrong case");
    }

    return std::make_tuple(static_matrix<T, 3, 3>::Zero(), static_matrix<T, 3, 3>::Zero(), T(0));
}

// compute projector tensor
template<typename T>
std::pair<static_tensor<T, 3>, static_tensor<T, 3>>
compute_projector(const static_matrix<T, 3, 3>& F,
                  const static_matrix<T, 3, 3>& stress_T,
                  const static_vector<T, 3>&    lambda_i,
                  const static_matrix<T, 3, 3>& evec,
                  bool                          compute_TL = true)
{
    typedef static_tensor<T, 3>    tensor_type;
    typedef static_matrix<T, 3, 3> matrix_type;
    typedef static_vector<T, 3>    vector_type;

    // projector tensor
    tensor_type P  = tensor_type::Zero();
    tensor_type TL = tensor_type::Zero();

    // compute Normal
    const std::array<vector_type, 3> Ni = compute_Ni(evec);
    const std::array<vector_type, 3> ni = compute_ni(F, Ni);

    // compute quantites
    const vector_type ei   = compute_ei(lambda_i);
    const vector_type di   = compute_di(lambda_i);
    const vector_type fi   = compute_fi(lambda_i);
    const matrix_type zeta = compute_zeta(stress_T, Ni);

    const auto        coefficient = compute_coefficient(lambda_i, ei, di, fi);
    const matrix_type theta       = std::get<0>(coefficient);
    const matrix_type xi          = std::get<1>(coefficient);
    const T           eta         = std::get<2>(coefficient);

    // std::cout << "T" << std::endl;
    // std::cout << stress_T << std::endl;
    // std::cout << "F" << std::endl;
    // std::cout << F << std::endl;

    // std::cout << "ei" << std::endl;
    // std::cout << ei.transpose() << std::endl;
    // std::cout << "di" << std::endl;
    // std::cout << di.transpose() << std::endl;
    // std::cout << "fi" << std::endl;
    // std::cout << fi.transpose() << std::endl;
    // std::cout << "zeta" << std::endl;
    // std::cout << zeta << std::endl;
    // std::cout << "eta" << std::endl;
    // std::cout << eta << std::endl;
    // std::cout << "theta" << std::endl;
    // std::cout << theta << std::endl;
    // std::cout << "xi" << std::endl;
    // std::cout << xi << std::endl;

    for (int i = 0; i < 3; i++)
    {
        const matrix_type Mii      = compute_Mij(Ni, ni, i, i);
        const matrix_type di2_NiNi = di(i) / T(2) * computeKroneckerProduct(Ni[i], Ni[i]);

        P += computeKroneckerProduct(di2_NiNi, Mii);

        if (compute_TL)
        {
            const T fiZeta = fi(i) / T(4) * zeta(i, i);
            TL += fiZeta * computeKroneckerProduct(Mii, Mii);
        }

        for (int j = 0; j < 3; j++)
        {
            if (j != i)
            {
                const matrix_type Mij        = compute_Mij(Ni, ni, i, j);
                const matrix_type theta_NiNj = theta(i, j) * computeKroneckerProduct(Ni[i], Ni[j]);

                P += computeKroneckerProduct(theta_NiNj, Mij);

                // for TL
                const matrix_type zeta_Mjj    = zeta(i, j) * compute_Mij(Ni, ni, j, j);
                const tensor_type zeta_MijMjj = computeKroneckerProduct(Mij, zeta_Mjj);
                const tensor_type zeta_MjjMij = computeKroneckerProduct(zeta_Mjj, Mij);
                const tensor_type MijMij      = computeKroneckerProduct(Mij, Mij);

                if (compute_TL)
                {
                    const T twoxi = T(2) * xi(i, j);
                    TL += twoxi * ((zeta_MijMjj + zeta_MjjMij) + zeta(j, j) * MijMij);

                    for (int k = 0; k < 3; k++)
                    {
                        if (k != i && k != j)
                        {
                            const T           etaZeta     = T(2) * eta * zeta(i, j);
                            const matrix_type etaZeta_Mik = etaZeta * compute_Mij(Ni, ni, i, k);
                            const matrix_type Mjk         = compute_Mij(Ni, ni, j, k);

                            TL += computeKroneckerProduct(etaZeta_Mik, Mjk);
                        }
                    }
                }
            }
        }
    }

    if (compute_TL)
    {
        // last term partie symmetric
        const matrix_type inv_F   = F.inverse();
        const matrix_type invF_TP = inv_F * computeContractedProduct<T, 3>(stress_T, P);
        const matrix_type Id      = matrix_type::Identity();
        TL += symetric_part<T, 3>(computeKroneckerProduct(invF_TP, Id));
    }
    // std::cout << "P" << std::endl;
    // std::cout << P << std::endl;
    // std::cout << "TL" << std::endl;
    // std::cout << TL << std::endl;

    return std::make_pair(P, TL);
}

// compute projector tensor
template<typename T>
std::pair<static_tensor<T, 3>, static_tensor<T, 3>>
compute_projector_PK2(const static_matrix<T, 3, 3>& F,
                      const static_matrix<T, 3, 3>& stress_T,
                      const static_vector<T, 3>&    lambda_i,
                      const static_matrix<T, 3, 3>& evec)
{
    typedef static_tensor<T, 3>    tensor_type;
    typedef static_matrix<T, 3, 3> matrix_type;
    typedef static_vector<T, 3>    vector_type;

    // projector tensor
    tensor_type P  = tensor_type::Zero();
    tensor_type TL = tensor_type::Zero();

    // compute Normal
    const std::array<vector_type, 3> Ni = compute_Ni(evec);

    // compute quantites
    const vector_type ei   = compute_ei(lambda_i);
    const vector_type di   = compute_di(lambda_i);
    const vector_type fi   = compute_fi(lambda_i);
    const matrix_type zeta = compute_zeta(stress_T, Ni);

    const auto        coefficient = compute_coefficient(lambda_i, ei, di, fi);
    const matrix_type theta       = std::get<0>(coefficient);
    const matrix_type xi          = std::get<1>(coefficient);
    const T           eta         = std::get<2>(coefficient);

    // std::cout << "T" << std::endl;
    // std::cout << stress_T << std::endl;
    // std::cout << "F" << std::endl;
    // std::cout << F << std::endl;

    // std::cout << "ei" << std::endl;
    // std::cout << ei.transpose() << std::endl;
    // std::cout << "zeta" << std::endl;
    // std::cout << zeta << std::endl;

    for (int i = 0; i < 3; i++)
    {
        const matrix_type Mii  = compute_Mij(Ni, Ni, i, i);
        const matrix_type NiNi = computeKroneckerProduct(Ni[i], Ni[i]);

        P += di(i) / T(2) * computeKroneckerProduct(NiNi, Mii);
        TL += fi(i) / T(4) * zeta(i, i) * computeKroneckerProduct(Mii, Mii);

        for (int j = 0; j < 3; j++)
        {
            if (j != i)
            {
                const matrix_type Mij  = compute_Mij(Ni, Ni, i, j);
                const matrix_type NiNj = computeKroneckerProduct(Ni[i], Ni[j]);

                P += theta(i, j) * computeKroneckerProduct(NiNj, Mij);

                // for TL
                const matrix_type Mjj    = compute_Mij(Ni, Ni, j, j);
                const tensor_type MijMjj = computeKroneckerProduct(Mij, Mjj);
                const tensor_type MjjMij = computeKroneckerProduct(Mjj, Mij);
                const tensor_type MijMij = computeKroneckerProduct(Mij, Mij);

                TL += T(2) * xi(i, j) * (zeta(i, j) * (MijMjj + MjjMij) + zeta(j, j) * MijMij);

                for (int k = 0; k < 3; k++)
                {
                    if (k != i && k != j)
                    {
                        const matrix_type Mik = compute_Mij(Ni, Ni, i, k);
                        const matrix_type Mjk = compute_Mij(Ni, Ni, j, k);

                        TL += T(2) * eta * zeta(i, j) * computeKroneckerProduct(Mik, Mjk);
                    }
                }
            }
        }
    }

    // std::cout << "P" << std::endl;
    // std::cout << P << std::endl;
    // std::cout << "TL" << std::endl;
    // std::cout << TL << std::endl;

    return std::make_pair(P, TL);
}
}
}