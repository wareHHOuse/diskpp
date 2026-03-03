/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2026
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#pragma once

#ifdef HAVE_MUMPS
#include "external/mumps.hpp"
#endif

#include "diskpp/common/eigen.hpp"

namespace disk::solvers {

enum class direct_solver {
    autosel,
#ifdef HAVE_MUMPS
    mumps,
#endif
    sparselu,
#ifdef HAVE_PARDISO
    pardiso
#endif
};

enum class direct_solver_status {
    ok,
    failure,
};

template<typename T, int _Options, typename _Index, _Index nrhs>
direct_solver_status
sparse_lu(Eigen::SparseMatrix<T, _Options, _Index>& A,
    Eigen::Matrix<T, Eigen::Dynamic, nrhs>& b,
    Eigen::Matrix<T, Eigen::Dynamic, nrhs>& x,
    direct_solver solver = direct_solver::autosel)
{
    using spmat_t = Eigen::SparseMatrix<T, _Options, _Index>;
    switch (solver) {
    /* Beware here: if we have mumps, autosel falls back to mumps, otherwise
     * it falls back to SparseLU (which is the only we are guaranteed to have)
     */

    case direct_solver::autosel:

#ifdef HAVE_MUMPS
    case direct_solver::mumps: {
        x = mumps_lu(A, b);
        /* check needed */
    } break;
#endif

    case direct_solver::sparselu: {
        Eigen::SparseLU<spmat_t> LU(A);
        x = LU.solve(b);
        if (LU.info() != Eigen::Success) {
            return direct_solver_status::failure;
        }
    } break;

#ifdef HAVE_PARDISO
    case direct_solver::pardiso: {
        Eigen::PardisoLU<spmat_t> LU(A);
        x = LU.solve(b);
        if (LU.info() != Eigen::Success) {
            return direct_solver_status::failure;
        }
    } break;
#endif

    default:
        std::cerr << "sparse_lu(): inexistent solver selected" << std::endl;
        return  direct_solver_status::failure;
    } // switch

   return  direct_solver_status::ok;
}

template<typename T, int _Options, typename _Index, _Index nrhs>
direct_solver_status
sparse_ldlt(Eigen::SparseMatrix<T, _Options, _Index>& A,
    Eigen::Matrix<T, Eigen::Dynamic, nrhs>& b,
    Eigen::Matrix<T, Eigen::Dynamic, nrhs>& x,
    direct_solver solver = direct_solver::autosel)
{
    using spmat_t = Eigen::SparseMatrix<T, _Options, _Index>;
    switch (solver) {
    /* Beware here: if we have mumps, autosel falls back to mumps, otherwise
     * it falls back to SparseLU (which is the only we are guaranteed to have)
     */

    case direct_solver::autosel:

#ifdef HAVE_MUMPS
    case direct_solver::mumps: {
        x = mumps_ldlt(A, b);
        /* check needed */
    } break;
#endif

    case direct_solver::sparselu: {
        Eigen::SimplicialLDLT<spmat_t> LDLT(A);
        x = LDLT.solve(b);
        if (LDLT.info() != Eigen::Success) {
            return direct_solver_status::failure;
        }
    } break;

#ifdef HAVE_PARDISO
    case direct_solver::pardiso: {
        Eigen::PardisoLDLT<spmat_t> LDLT(A);
        x = LDLT.solve(b);
        if (LDLT.info() != Eigen::Success) {
            return direct_solver_status::failure;
        }
    } break;
#endif

    default:
        std::cerr << "sparse_ldlt(): inexistent solver selected" << std::endl;
        return  direct_solver_status::failure;
    } // switch

   return  direct_solver_status::ok;
}

};