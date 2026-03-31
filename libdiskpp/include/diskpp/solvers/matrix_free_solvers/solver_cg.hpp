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

#include "diskpp/solvers/defs.hpp"

namespace disk::solvers {


template<typename T, typename Oper, typename Precond>
iterative_solver_status
cg_mf(const iterative_solver_params& params, Oper A,
    const priv::vec<T>& b, priv::vec<T>& x, Precond precond)
{
    using dv = priv::vec<T>;
    using real = remove_complex<T>;

    assert( b.size() == x.size() );

    double  nr, nr0;
    T  alpha, beta, rho;
    auto N = b.size();
    dv d(N), r(N), r0(N), y(N), Pr(N), Pr1(N);

    r0 = r = b - A(x);
    d = precond(r0);
    nr = nr0 = r.norm();

    for (size_t iter = 0; iter < params.max_iter; iter++)
    {
        if (params.verbose) {
            std::cout << "\rCG iteration " << iter << std::flush;
        }

        double rr = nr/nr0;
        if (rr < params.tol) {
            return iterative_solver_status::converged;
        }

        if (rr > params.max_residual) {
            return iterative_solver_status::diverged;
        }

        y     = A(d);
        Pr    = precond(r);
        rho   = r.dot(Pr);
        alpha = rho / d.dot(y);
        x     = x + alpha * d;
        r     = r - alpha * y;
        Pr1   = precond(r);
        beta  = r.dot(Pr1) / rho;
        d     = Pr1 + beta * d;

        nr = r.norm();
    }
    
    if (params.verbose) {
        std::cout << std::endl;
    }

    return iterative_solver_status::hit_max_iter;
}

} // namespace disk::solvers