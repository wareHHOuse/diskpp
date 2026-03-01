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

template<typename T, oper<T> SolveB>
void block_jacobi_davidson(const bjd_params& params,
    auto applyA, const priv::spmat<T>& B, SolveB solveB,
    priv::vec<T>& eigenvalues,
    priv::mat<T>& eigenvectors)
{
    const int N = B.rows();

    using dv = priv::vec<T>;
    using dm = priv::mat<T>;

    dm V = dm::Random(N, params.block_size);
    orthonormalize(V, B);

    // Keeps track of converged eigenpairs
    std::vector<bool> conv(params.block_size, false);

    for (int outer = 0; outer < params.max_outer_iters; outer++) {
        std::cout << "outer: " << outer << std::endl;
        
        const int M = V.cols();
        dm TA = V.transpose() * applyA(V);
        //dm TB = V.transpose() * B * V;
        //dm TB = dm::Identity(TA.rows(), TA.cols());
        Eigen::SelfAdjointEigenSolver<dm> ges(TA);
        dv theta = ges.eigenvalues();
        dm Y = ges.eigenvectors();
        dm X = V * Y.leftCols(params.block_size);

        if ( V.cols() > params.block_size * params.max_subspace_growth ) {
            V = X.leftCols(params.block_size);
        }

        for (int i_ev = 0; i_ev < params.block_size; i_ev++) {
            if (conv[i_ev]) {
                continue; /* This pair has already converged. */
            }

            //std::cout << "  Eigenvalue " << i_ev << std::endl;
    
            dv xi = X.col(i_ev);
            T lambda = theta(i_ev);
            dv r = applyA(xi) - lambda * (B * xi);

            //std::cout << "   - norm " << i_ev << ": " << r.norm() << std::endl;
            if (r.norm() < params.ev_tol * std::abs(lambda)) {
                //std::cout << "  locked" << std::endl;
                conv[i_ev] = true;
                continue;
            }

            auto JDproj = [&](const dv& s) -> dv {
                dv Bs = B*s;
                dv v1 = s - xi * xi.dot(Bs); // (I - u*u'*B)*s
                dv v2 = applyA(v1) - lambda*B*v1; //(A - lambda*B)*v1
                return v2 - B*(xi * xi.dot(v2)); //(I - B*x*x')*v2
            };

            disk::solvers::iterative_solver_params tfqmr_p;
            tfqmr_p.max_iter = params.max_inner_iters;
            tfqmr_p.tol = 0.1*r.norm();
            tfqmr_p.verbose = false;
            dv s = dv::Zero( r.rows() );
            tfqmr_mf(tfqmr_p, JDproj, (-r).eval(), s, disk::solvers::identity{});
            
            if (s.norm() > params.expand_thresh) {
                //std::cout << "  - Expanding subspace. Norm: " << s.norm();
                //std::cout << ", size: " << V.cols()+1 << std::endl;
                V.conservativeResize(N, V.cols() + 1);
                V.col(V.cols() - 1) = s;
            }
        }

        orthonormalize(V, B);

        // Check if all the eigenpairs have converged
        if ( std::all_of(conv.begin(), conv.end(), [](bool v) { return v; }) ) {
            std::cout << "All eigenpairs have converged" << std::endl;
            break;
        }
    }

    // Final Rayleigh-Ritz to compute eigenvalues/eigenvectors
    int M = V.cols();
    dm TA = V.transpose() * applyA(V);
    //dm TB = V.transpose() * B * V;
    Eigen::SelfAdjointEigenSolver<dm> ges(TA);

    eigenvalues.resize(params.block_size);
    eigenvectors = V * ges.eigenvectors().leftCols(params.block_size);
    for (int i = 0; i < params.block_size; i++)
        eigenvalues(i) = ges.eigenvalues()(i);
}

} // namespace disk::solvers