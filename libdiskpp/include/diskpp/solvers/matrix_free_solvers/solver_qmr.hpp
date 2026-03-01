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
tfqmr_mf(const iterative_solver_params& params, Oper A,
    const priv::vec<T>& b, priv::vec<T>& x, Precond precond)
{
    using dv = priv::vec<T>;

    assert( b.size() == x.size() );

    dv r = b - A(x); 
    if (r.norm() < params.tol) {
        return iterative_solver_status::converged;
    }

    dv r_tld = r;

    dv u = r;
    dv p = r;

    dv v = precond(A(p));

    double rho = r_tld.dot(r);
    double alpha = 0.0;
    double omega = 1.0;

    dv d = dv::Zero( b.size() );
    dv w = r;

    double tau = r.norm();
    double theta = 0.0;
    double eta = 0.0;

    for (int i = 0; i < params.max_iter; i++)
    {
        double sigma = r_tld.dot(v);
        if (std::abs(sigma) < 1e-30) {
            if (params.verbose) {
                std::cout << "QMR breakdown at iteration " << i;
                std::cout << std::endl;
            }
            return iterative_solver_status::breakdown;
        }
        
        alpha = rho / sigma;
        dv q = u - alpha * v;
        dv u_new = u + q;
        dv y = precond(u_new);
        x += alpha * y;
        w -= alpha * precond(A(u_new));

        double rho_new = r_tld.dot(w);
        if (std::abs(rho) < 1e-30) {
            if (params.verbose) {
                std::cout << "QMR breakdown at iteration " << i;
                std::cout << std::endl;
            }
            return iterative_solver_status::breakdown;
        }

        double beta = rho_new / rho;
        rho = rho_new;

        u = w + beta * q;
        p = u + beta * (q + beta * p);

        v = precond(A(p));

        if (w.norm() < params.tol) {
            if (params.verbose) {
                std::cout << "QMR converged with residual ";
                std::cout << w.norm() << " at iteration " << i << std::endl;
            }
            return iterative_solver_status::converged;
        }

        if (w.norm() > params.max_residual) {
            if (params.verbose) {
                std::cout << "QMR DIVERGED with residual ";
                std::cout << w.norm() << " at iteration " << i << std::endl;
            }
            return iterative_solver_status::diverged;
        }
    
    }

    if (params.verbose) {
        std::cout << "QMR reached max iterations with residual ";
        std::cout << w.norm() << std::endl;
    }

    return iterative_solver_status::hit_max_iter;
}

} // namespace disk::solvers