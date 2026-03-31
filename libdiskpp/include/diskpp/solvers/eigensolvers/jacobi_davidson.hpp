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

#include "diskpp/common/eigen.hpp"
#include "diskpp/solvers/matrix_free_solvers.hpp"

namespace disk::solvers {


struct bjd_params {
    int     max_outer_iters         = 100;
    int     max_inner_iters         = 30;
    double  outer_tol               = 1e-8;
    double  inner_tol               = 1e-3;
    double  ev_tol                  = 1e-8;
    double  expand_thresh           = 1e-12;
    int     max_subspace_growth     = 5;
    int     block_size              = 10;
    bool    verbose                 = false;
};

enum class bdj_status {
    converged,
    hit_max_iter,
};

template<typename T>
void orthonormalize(priv::mat<T>& V, const priv::spmat<T>& B)
{
    // Modified G-S B-orthonormalization
    for (int i = 0; i < V.cols(); i++) {
        double Bnorm = std::sqrt(V.col(i).dot(B*V.col(i)));
        V.col(i) /= Bnorm;
        for (int j = i+1; j < V.cols(); j++) {
            double proj = V.col(i).dot(B*V.col(j));
            V.col(j) = V.col(j) - proj * V.col(i);
        }
    }
}

/* This is a block Jacobi-Davidson solver to find the smallest
 * eigenvalues of Ax = λBx. A and B are assumed SPD, the inner
 * solver is QMR.
 * The implementation is matrix-free and is suited for HHO
 * when the static condensation is done on the cells (and
 * therefore must be done globally and not locally).
 * If you don't need a matrix-free solver you should probably
 * use FEAST.
 */
template<typename T>
bdj_status
block_jacobi_davidson(const bjd_params& params,
    auto applyA, const priv::spmat<T>& B,
    priv::mat<T>& eigenvectors,
    priv::vec<T>& eigenvalues)
{
    const int N = B.rows();

    using dv = priv::vec<T>;
    using dm = priv::mat<T>;

    dm V = dm::Random(N, params.block_size);
    orthonormalize(V, B);

    // Keep track of converged eigenpairs and block
    // them during iteration to save operations
    std::vector<bool> conv(params.block_size, false);

    disk::solvers::iterative_solver_params innerp;
    innerp.max_iter = params.max_inner_iters;

    double prev_r = 1;

    for (int outer = 0; outer < params.max_outer_iters; outer++) {
        const int M = V.cols();
        dm TA = V.transpose() * applyA(V);
        /* Since we orthonormalized V, the matrix V'*B*V is the
         * identity, so we don't need to solve a generalized
         * problem here */
        Eigen::SelfAdjointEigenSolver<dm> es(TA);
        dv theta = es.eigenvalues();
        dm Y = es.eigenvectors();
        dm X = V * Y.leftCols(params.block_size); //smallest eigs

        /* If the subspace gets too big we restart */
        if ( V.cols() > params.block_size * params.max_subspace_growth ) {
            V = X.leftCols(params.block_size);
        }

        double maxnorm = 0.0;
    
        std::cout << "outer " << outer << ": " << std::flush; 

        for (int i_ev = 0; i_ev < params.block_size; i_ev++) {
            if (conv[i_ev]) {
                continue;
            }
    
            dv xi = X.col(i_ev);
            T lambda = theta(i_ev);
            dv r = applyA(xi) - lambda * (B * xi);
            if (r.norm() < params.ev_tol * std::abs(lambda)) {
                conv[i_ev] = true;
                continue;
            }

            maxnorm = std::max(r.norm(), maxnorm);

            /* This is the LHS for the correction equation */
            auto corrLHS = [&](const dv& s) -> dv {
                dv Bs = B*s;
                dv v1 = s - xi * xi.dot(Bs);
                dv v2 = applyA(v1) - lambda*B*v1;
                return v2 - B*(xi * xi.dot(v2));
            };

            /* Call iterative solver to solve the correction equation */
            
            innerp.tol = 0.1*r.norm();
            innerp.verbose = false;
            dv s = dv::Zero( r.rows() );
            using id = disk::solvers::identity;
            auto status = bicgstab_mf(innerp, corrLHS, (-r).eval(), s, id{});

            switch (status)
            {
                case iterative_solver_status::converged: std::cout << 'C'; break;
                case iterative_solver_status::diverged: std::cout << 'D'; break;
                case iterative_solver_status::hit_max_iter: std::cout << 'I'; break;
            }
            std::cout << std::flush;

            /* Add the solution to the current search space */
            if (s.norm() > params.expand_thresh) {
                V.conservativeResize(N, V.cols() + 1);
                V.col(V.cols() - 1) = s;
            }
        }
        //std::cout << maxnorm << " " << prev_r << std::endl;

        if ( (maxnorm/prev_r) > 1.1 ) {
            innerp.max_iter *= 2;
        }
        prev_r = maxnorm;

        std::cout << " maxnorm = " << maxnorm << ", subspace size = " << V.cols() << " " << innerp.max_iter << "\n";

        /* Check if all the eigenpairs have converged. In that case copy
         * them in the output parameters and return, we're finished.
         */
        if ( std::all_of(conv.begin(), conv.end(), [](bool v) { return v; }) ) {
            std::cout << "All eigenpairs have converged" << std::endl;
            eigenvalues.resize(params.block_size);
            eigenvalues = theta.segment(0, params.block_size);
            eigenvectors = X.leftCols(params.block_size);
            return bdj_status::converged;
        }

        /* Otherwise we orthonormalize the whole V and we go for another
         * round. It is not necessary to go through the whole V, but I am
         * lazy and for now (read: forever) it's OK like this. */
        orthonormalize(V, B);
    }

    return bdj_status::hit_max_iter;
}

} // namespace disk::solvers