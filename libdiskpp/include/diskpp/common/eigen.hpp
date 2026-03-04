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
#include <Eigen/SparseCholesky>
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

template<typename T>
void dump_sparse_matrix(Eigen::SparseMatrix<T>& M, const std::string& filename)
{
    std::ofstream ofs(filename);

    for (int k=0; k < M.outerSize(); ++k)
        for (typename Eigen::SparseMatrix<T>::InnerIterator it(M,k); it; ++it)
            ofs << it.row() << " " << it.col() << " " << std::setprecision(20) << it.value() << std::endl;

    ofs.close();
}

template<typename T>
void dump_sparse_matrix_with_header(Eigen::SparseMatrix<T>& M, const std::string& filename)
{
    std::ofstream ofs(filename);

    ofs << M.rows() << " " << M.cols() << " " << M.nonZeros() << std::endl;

    for (int k=0; k < M.outerSize(); ++k)
        for (typename Eigen::SparseMatrix<T>::InnerIterator it(M,k); it; ++it)
            ofs << it.row() << " " << it.col() << " " << std::setprecision(20) << it.value() << std::endl;

    ofs.close();
}

} //namespace disk



#ifdef HAVE_HDF5

#include <highfive/highfive.hpp>
#include <highfive/H5Easy.hpp>

/* To load a sparse matrix saved with this in Matlab,
 * just use the following commands:
 *
 * nrows = h5read("matrix.h5", "/sparsematrix/nrows");
 * ncols = h5read("matrix.h5", "/sparsematrix/ncols");
 * is = h5read("matrix.h5", "/sparsematrix/is");
 * js = h5read("matrix.h5", "/sparsematrix/js");
 * vals = h5read("matrix.h5", "/sparsematrix/vals");
 * A = sparse(is+1, js+1, vals, nrows, ncols);
 */

namespace disk {

template<typename T>
bool save_to_hdf5(Eigen::SparseMatrix<T>& M, const std::string& filename)
{
    M.makeCompressed();
    int nrows = M.rows();
    int ncols = M.cols();
    std::vector<int>    is;
    std::vector<int>    js;
    std::vector<T>      vals;

    for (int k=0; k < M.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<T>::InnerIterator it(M,k); it; ++it) {
            is.push_back( it.row() );
            js.push_back( it.col() );
            vals.push_back( it.value() );
        }
    }

    HighFive::File file(filename, HighFive::File::Truncate);
    file.createDataSet("/sparsematrix/nrows", nrows);
    file.createDataSet("/sparsematrix/ncols", ncols);
    file.createDataSet("/sparsematrix/is", is);
    file.createDataSet("/sparsematrix/js", js);
    file.createDataSet("/sparsematrix/vals", vals);

    return true;
}

template<typename T>
bool load_from_hdf5(Eigen::SparseMatrix<T>& M, const std::string& filename)
{
    using trip_t = Eigen::Triplet<T>;
    std::vector<trip_t> trips;

    HighFive::File file(filename, HighFive::File::ReadOnly);

    auto nrows_ds = file.getDataSet("/sparsematrix/nrows");
    auto nrows = nrows_ds.read<int>();
    auto ncols_ds = file.getDataSet("/sparsematrix/ncols");
    auto ncols = ncols_ds.read<int>();
    auto is_ds = file.getDataSet("/sparsematrix/is");
    auto is = is_ds.read<std::vector<int>>();
    auto js_ds = file.getDataSet("/sparsematrix/js");
    auto js = js_ds.read<std::vector<int>>();
    auto vals_ds = file.getDataSet("/sparsematrix/vals");
    auto vals = vals_ds.read<std::vector<T>>();

    bool sizes_ok = (is.size() == js.size()) and (js.size() == vals.size());
    if (not sizes_ok) {
        return false;
    } 

    M.resize(nrows, ncols);

    for (size_t i = 0; i < vals.size(); i++) {
        trips.push_back( {is[i], js[i], vals[i]} );
    }

    M.setFromTriplets(trips.begin(), trips.end());

    return true;
}

/* To read dense objects from matlab, just use
 * 
 * Z = h5read("dense.h5", "/densematrix");
 * 
 * BEWARE THAT MATLAB LOADS THE TRANSPOSE OF WHAT YOU SAVE 
 */

template<typename T, int nrows, int ncols>
bool save_to_hdf5(Eigen::Matrix<T, nrows, ncols>& M, const std::string& filename)
{
    H5Easy::File file(filename, H5Easy::File::Truncate);
    H5Easy::dump(file, "/densematrix", M);
    return true;
}

template<typename T, int nrows, int ncols>
bool load_from_hdf5(Eigen::Matrix<T, nrows, ncols>& M, const std::string& filename)
{
    using mtype = Eigen::Matrix<T, nrows, ncols>;
    H5Easy::File file(filename, H5Easy::File::ReadOnly);
    M = file.getDataSet("/densematrix").read<mtype>();
    return true;
}

} // namespace disk

#endif


