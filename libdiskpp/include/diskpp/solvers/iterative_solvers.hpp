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

#include "matrix_free_solvers/solver_cg.hpp"
#include "matrix_free_solvers/solver_bicgstab.hpp"
#include "matrix_free_solvers/solver_minres.hpp"
#include "matrix_free_solvers/solver_qmr.hpp"

namespace disk::solvers {

template<typename T>
iterative_solver_status
cg(const iterative_solver_params& params, const priv::spmat<T>& A,
    const priv::vec<T>& b, priv::vec<T>& x)
{
    auto operA = operator_from_matrix(A);
    auto id = identity{};

    return cg_mf(params, operA, b, x, id);
}

template<typename T, typename Precond>
iterative_solver_status
cg(const iterative_solver_params& params, const priv::spmat<T>& A,
    const priv::vec<T>& b, priv::vec<T>& x, Precond M)
{
    auto operA = operator_from_matrix(A);

    return cg_mf(params, operA, b, x, M);
}

} // namespace disk::solvers