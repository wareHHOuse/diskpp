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

#include <concepts>
#include <utility>
#include <complex>

#include "diskpp/common/eigen.hpp"

namespace disk::solvers {

struct iterative_solver_params {
    int         max_iter = 100;
    double      max_residual = 1e5;
    double      tol = 1e-8;
    bool        verbose = false;   
};

enum class iterative_solver_status {
    undefined,
    converged,
    hit_max_iter,
    diverged,
    breakdown,
};

namespace priv {
template<typename T>
using vec = disk::dynamic_vector<T>;

template<typename T>
using mat = disk::dynamic_matrix<T>;

template<typename T>
using spmat = disk::sparse_matrix<T>;
}

template<typename F, typename T>
concept oper =
    std::invocable<F, priv::vec<T>> &&
    std::same_as<std::invoke_result_t<F, priv::vec<T>>, priv::vec<T>>;

template<typename T>
struct remove_complex {
    using type = T;
};

template<typename T>
struct remove_complex<std::complex<T>> {
    using type = T;
};

template<typename T>
using remove_complex_t = typename remove_complex<T>::type;

template<typename T>
struct operator_from_matrix {
    const priv::mat<T>&  M;

    operator_from_matrix(const priv::mat<T>& pM) : M(pM) {}

    priv::vec<T> operator()(const priv::vec<T>& v) const {
        return M*v;
    }
};
struct identity {
    template<typename T>
    T&& operator()(T&& v) const noexcept {
        return std::forward<T>(v);
    }
};

}