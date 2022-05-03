/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
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

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wall"
#pragma clang diagnostic ignored "-Wshadow"
#pragma clang diagnostic ignored "-Weverything"

#ifdef HAVE_INTEL_MKL
/* Don't use MKL! It makes everything slower! */
//#define EIGEN_USE_MKL_ALL
// Fix for eigen version > 3.3.7
#if !defined(EIGEN_USING_STD)
#define EIGEN_USING_STD(X) using std::X
#endif
#include <Eigen/PardisoSupport>
#endif

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <unsupported/Eigen/SparseExtra>

#include <Eigen/StdVector>

#include <iomanip>

#pragma clang diagnostic pop

  namespace disk{

template<typename T>
using dynamic_matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
using dynamic_vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename T, int M, int N>
using static_matrix = Eigen::Matrix<T, M, N>;

template<typename T, int N>
using static_vector = Eigen::Matrix<T, N, 1>;

// to mimic fourth order tensor
template<typename T, int N>
using static_tensor = Eigen::Matrix<T, N * N, N * N>;

template<typename T, size_t M, size_t N>
using material_tensor = static_matrix<T, M, N>;

template<typename T>
using sparse_matrix = Eigen::SparseMatrix<T>;

template<typename T>
using triplet = Eigen::Triplet<T>;

template<typename T>
using eigen_compatible_stdvector = std::vector<T, Eigen::aligned_allocator<T>>;

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

template<typename T>
class eigen_sparse_raw
{
    bool m_one_based, m_already_restored;
public:

    T *     data;
    int *   ia;
    int *   ja;
    int     nnz;
    int     n;

    template<int Options, typename Index>
    eigen_sparse_raw(Eigen::SparseMatrix<T, Options, Index>& A, bool one_based = false)
        : m_one_based(one_based), m_already_restored(false)
    {
        A.makeCompressed();

        data    = A.valuePtr();
        ia      = A.outerIndexPtr();    //colmaj size(ia) == N+1
        ja      = A.innerIndexPtr();    //colmaj size(ja) == size(data) == nnz
        nnz     = A.nonZeros();
        n       = A.rows();

        if (one_based)
        {
            for (size_t i = 0; i < n+1; i++)
                ia[i] += 1;

            for (size_t i = 0; i < nnz; i++)
                ja[i] += 1;
        }
    }

    void restore(void)
    {
        if ( !m_one_based || m_already_restored )
            return;

        for (size_t i = 0; i < n+1; i++)
            ia[i] -= 1;

        for (size_t i = 0; i < nnz; i++)
            ja[i] -= 1;

        m_already_restored = true;
    }

    void show(void)
    {
        std::cout << "A: ";
        for (size_t i = 0; i < nnz; i++)
            std::cout << data[i] << " ";
        std::cout << std::endl;

        std::cout << "ja: ";
        for (size_t i = 0; i < nnz; i++)
            std::cout << ja[i] << " ";
        std::cout << std::endl;

        std::cout << "ia: ";
        for (size_t i = 0; i < n+1; i++)
            std::cout << ia[i] << " ";
        std::cout << std::endl;
    }

    ~eigen_sparse_raw()
    {
        if ( m_one_based && !m_already_restored )
            restore();
    }
};

void dump_sparse_matrix(Eigen::SparseMatrix<double>& M, const std::string& filename)
{
    std::ofstream ofs(filename);

    for (int k=0; k < M.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(M,k); it; ++it)
            ofs << it.row() << " " << it.col() << " " << std::setprecision(20) << it.value() << std::endl;

    ofs.close();
}

void dump_sparse_matrix_with_header(Eigen::SparseMatrix<double>& M, const std::string& filename)
{
    std::ofstream ofs(filename);

    ofs << M.rows() << " " << M.cols() << " " << M.nonZeros() << std::endl;

    for (int k=0; k < M.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(M,k); it; ++it)
            ofs << it.row() << " " << it.col() << " " << std::setprecision(20) << it.value() << std::endl;

    ofs.close();
}

} // end disk
