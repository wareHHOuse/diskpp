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

#include "../../config.h"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wall"
#pragma clang diagnostic ignored "-Wshadow"

#ifdef HAVE_INTEL_MKL
    /* Don't use MKL! It makes everything slower! */
    //#define EIGEN_USE_MKL_ALL
    #include <Eigen/PardisoSupport>
#endif

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <unsupported/Eigen/SparseExtra>

#pragma clang diagnostic pop

template<typename T>
using dynamic_matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
using dynamic_vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename T, size_t M, size_t N>
using static_matrix = Eigen::Matrix<T, M, N>;

template<typename T, size_t N>
using static_vector = Eigen::Matrix<T, N, 1>;

template<typename T, size_t M, size_t N>
using material_tensor = static_matrix<T, M, N>;

template<typename T>
using sparse_matrix = Eigen::SparseMatrix<T>;

template<typename T>
using triplet = Eigen::Triplet<T>;

template<typename T>
static_vector<T, 3>
cross(const static_vector<T, 2>& v1, const static_vector<T, 2>& v2)
{
    static_vector<T, 3> ret;

    ret(0) = T(0);
    ret(1) = T(0);
    ret(2) = v1(0)*v2(1) - v1(1)*v2(0);

    return ret;
}
