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

/* Matrix-free BiCGStab. */
template<typename T, typename Oper, typename Precond>
iterative_solver_status
bicgstab_mf(const iterative_solver_params& params, Oper A,
    const priv::vec<T>& b, priv::vec<T>& x, Precond precond)
{
    using dv = priv::vec<T>;
    using real = remove_complex_t<T>;

    assert( b.size() == x.size() );

    dv r = b - A(x);
    dv rhat = b;
    T rho = rhat.dot(r);
    dv p = r;

    for (size_t i = 0; i < params.max_iter; i++) {
        dv y = precond(p);
        dv v = A(y);
        T rhat_dot_v = rhat.dot(v);
        if (rhat_dot_v == T{0.0}) {
            if (params.verbose) {
                std::cout << "BiCGStab breakdown at iteration " << i;
                std::cout << std::endl;
            }
            return iterative_solver_status::breakdown;
        }
        T alpha = rho/rhat_dot_v;
        dv h = x + alpha*y;
        dv s = r - alpha*v;
        dv z = precond(s);
        dv t = A(z);
        dv Mt = precond(t);
        T Mt_dot_Mt = Mt.dot(Mt);
        T omega = Mt.dot(z)/Mt_dot_Mt;
        x = h + omega*z;
        r = s - omega*t;

        real nr = r.norm();
        if ( nr < params.tol ) {
            if (params.verbose) {
                std::cout << "BiCGStab converged with residual ";
                std::cout << nr << " at iteration " << i << std::endl;
            }
            return iterative_solver_status::converged;
        }

        if ( nr > params.max_residual ) {
            if (params.verbose) {
                std::cout << "BiCGStab DIVERGED with residual ";
                std::cout << nr << " at iteration " << i << std::endl;
            }
            return iterative_solver_status::diverged;
        }

        T rho_prev = rho;
        rho = rhat.dot(r);
        T beta = (rho/rho_prev)*(alpha/omega);
        p = r + beta*(p - omega*v);
    }

    if (params.verbose) {
        std::cout << "BiCGStab reached max iterations with residual ";
        std::cout << r.norm() << std::endl;
    }

    return iterative_solver_status::hit_max_iter;
}

} // namespace disk::solvers