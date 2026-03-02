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
minres_mf(const iterative_solver_params& params, Oper A,
    const priv::vec<T>& b, priv::vec<T>& x, Precond precond)
{
    using dv = priv::vec<T>;

    assert( b.size() == x.size() );
    const int n = b.size();

    dv r0 = b - A(x);
    if (r0.norm() < params.tol) {
        return iterative_solver_status::converged;
    }

    dv v_old = dv::Zero(n);
    dv v = r0;
    dv z = precond(v);

    T beta1 = std::sqrt(v.dot(z));
    if (beta1 == T{0.0}) {
        return iterative_solver_status::breakdown;
    }

    v /= beta1;
    z /= beta1;

    dv w_old = dv::Zero(n);
    dv w = dv::Zero(n);

    T beta = 0.0;
    T c_old = 1.0, s_old = 0.0;
    T c = 1.0, s = 0.0;
    T eta = beta1;

    for (int k = 0; k < params.max_iter; ++k)
    {
        dv Av = A(v);
        dv Az = precond(Av);

        T alpha = v.dot(Az);

        dv r = Az - alpha * v - beta * v_old;

        dv z_new = precond(r);
        T beta_new = std::sqrt( std::real(r.dot(z_new)) );

        v_old = v;
        v = r;
        if (beta_new > 0) {
            v /= beta_new;
        }

        T delta = c_old * alpha - c_old * s_old * beta;
        T gamma = std::sqrt(delta * delta + beta_new * beta_new);

        if (gamma == T{0.0}) {
            return iterative_solver_status::breakdown;
        }

        c = delta / gamma;
        s = beta_new / gamma;

        T epsilon = s_old * alpha + c_old * beta;
        T zeta = c * eta;
        eta = -s * eta;

        dv w_new = (v_old - epsilon * w_old - delta * w) / gamma;
        x += zeta * w_new;

        if (std::abs(eta) < params.tol) {
            return iterative_solver_status::converged;
        }

        if (std::abs(eta) > params.max_residual) {
            return iterative_solver_status::diverged;
        }

        w_old = w;
        w = w_new;

        beta = beta_new;
        c_old = c;
        s_old = s;
    }

    return iterative_solver_status::hit_max_iter;
}

} // namespace disk::solvers