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

#include <iostream>
#include <vector>
#include <optional>
#include <utility>

#include "diskpp/common/eigen.hpp"
#include "diskpp/solvers/defs.hpp"

namespace disk::solvers {

template<typename T>
std::optional<std::pair<priv::vec<T>, T>>
powiter(const priv::spmat<T>& A,
    remove_complex_t<T> tol = 1e-10,
    size_t maxiter = 1000)
{
    assert(A.rows() == A.cols());
    using dv = priv::vec<T>;
    dv x = dv::Ones(A.cols());
    x = x/x.norm();
    T lambda = 0;
    size_t i = 0;
    for (; i < maxiter; i++) {
        dv new_x = A*x;
        T new_lambda = x.dot(new_x);
        x = new_x/new_x.norm();
        auto err = std::abs(lambda-new_lambda);
        if ( err <= tol*std::abs(new_lambda) ) {
            std::cout << std::endl;
            return std::pair{x, new_lambda};
        }
        std::cout << "\rPower iteration " << i << ": " << err/std::abs(new_lambda) << std::flush;
        lambda = new_lambda;
    }
    std::cout << std::endl;
    if (i == maxiter) {
        return {};
    }

    return std::pair{x, lambda};
}

template<typename T>
std::optional<std::pair<priv::vec<T>, T>>
inv_powiter(const priv::spmat<T>& A,
    remove_complex_t<T> mu,
    remove_complex_t<T> tol = 1e-10,
    size_t maxiter = 1000)
{
    using sm = priv::spmat<T>;
    using dv = priv::vec<T>;
    assert(A.rows() == A.cols());
    sm muI(A.rows(), A.cols());
    for (int i = 0; i < A.rows(); i++)
        muI.insert(i,i) = mu;
    sm M = A - muI;
    mumps_solver<T> s;
    s.factorize(M);
    
    dv x = dv::Ones(M.cols());
    x = x/x.norm();
    T lambda = mu;
    size_t i = 0;
    for (; i < maxiter; i++) {
        dv new_x = s.solve(x);
        T new_lambda = mu + 1.0/x.dot(new_x);
        x = new_x/new_x.norm();
        auto err = std::abs(lambda-new_lambda);
        if ( err <= tol*std::abs(new_lambda) ) {
            std::cout << std::endl;
            return std::pair{x, new_lambda};
        }
        std::cout << "\rInverse power iteration " << i << ": " << err/std::abs(new_lambda) << std::flush;
        lambda = new_lambda;
    }
    std::cout << std::endl;
    if (i == maxiter) {
        return {};
    }

    return std::pair{x, lambda};
}

template<typename T>
std::optional<std::pair<priv::vec<T>, T>>
inv_powiter(const priv::spmat<T>& A,
    const priv::spmat<T>& B,
    remove_complex_t<T> mu,
    remove_complex_t<T> tol = 1e-10,
    size_t maxiter = 1000)
{
    assert(A.rows() == A.cols());
    assert(B.rows() == B.cols());
    assert(A.rows() == B.rows());

    using sm = priv::spmat<T>;
    using dv = priv::vec<T>;
    
    sm M = A - mu*B;
    mumps_solver<T> s;
    s.factorize(M);
    dv x = dv::Ones(M.cols());
    x = x/x.norm();
    T lambda = mu;
    size_t i = 0;
    for (; i < maxiter; i++) {
        dv Bx = B*x;
        dv new_x = s.solve(Bx);
        T new_lambda = mu + 1.0/x.dot(new_x);
        x = new_x/new_x.norm();
        auto err = std::abs(lambda-new_lambda);
        if ( err <= tol*std::abs(new_lambda) ) {
            std::cout << std::endl;
            return std::pair{x, new_lambda};
        }
        std::cout << "\rInverse power iteration " << i << ": " << err/std::abs(new_lambda) << std::flush;
        lambda = new_lambda;
    }
    std::cout << std::endl;
    if (i == maxiter) {
        return {};
    }

    return std::pair{x, lambda};
}

} // namespace disk::solvers