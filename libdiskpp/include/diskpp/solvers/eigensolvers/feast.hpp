/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */
/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */
/*
 * This source file is part of EMT, the ElectroMagneticTool.
 *
 * Copyright (C) 2013-2015, Matteo Cicuttin - matteo.cicuttin@uniud.it
 * Department of Electrical Engineering, University of Udine
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Udine nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR(s) ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(s) BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include <iostream>

#include "diskpp/solvers/defs.hpp"

namespace disk::solvers {

enum class feast_inner_solver {
    eigen_sparselu,
    mumps
};

enum class feast_status {
    success,
    did_not_converge,
    invalid_input,
    inner_solver_problem,
    subspace_too_small
};

std::ostream& operator<<(std::ostream& os, const feast_status& fs)
{
    switch (fs) {
        case feast_status::success:
            os << "success";
            break;
        
        case feast_status::did_not_converge:
            os << "did_not_converge";
            break;

        case feast_status::invalid_input:
            os << "invalid_input";
            break;

        case feast_status::inner_solver_problem:
            os << "inner_solver_problem";
            break;

        case feast_status::subspace_too_small:
            os << "subspace_too_small";
            break;
    }

    return os;
}

template<typename T>
struct feast_eigensolver_params
{
    bool    verbose;
    int     tolerance;
    T       min_eigval;
    T       max_eigval;
    int     subspace_size;
    int     eigvals_found;
    int     feast_info;
    size_t  max_iter;
    feast_inner_solver fis = feast_inner_solver::eigen_sparselu;
};

static double quadrature_xs[] = {
    0.183434642495649, -0.183434642495649,
    0.525532409916328, -0.525532409916328,
    0.796666477413626, -0.796666477413626,
    0.960289856497536, -0.960289856497536
};

static double quadrature_ws[] = {
    0.362683783378361, 0.362683783378361,
    0.313706645877887, 0.313706645877887,
    0.222381034453374, 0.222381034453374,
    0.101228536290376, 0.101228536290376
};

/* This function implements the algorithm presented in
 *   "A Density Matrix-based Algorithm for Solving Eigenvalue Problems"
 *   by E. Polizzi. arXiv:0901.2665v1.
 *
 * I wrote this because the FEAST implemented in the Intel MKL crashes with
 * the following error
 *
 *   Intel MKL Extended Eigensolvers ERROR: Problem from Inner Linear System Solver
 *   ==>INFO code =: -2
 *
 * already on very small matrices (less that 10000x10000). INFO code -2 from PARDISO
 * means out of memory, however there are no memory problems at all. It looks like
 * that it has something to do with the pivoting strategy used in PARDISO. As a matter
 * of fact, also this implementation fails if instead of using Eigen::SparseLU or
 * mumps_lu() you use Eigen::PardisoLU.
 *
 * This implementation is slightly less robust than the MKL one (it requires you to
 * get a rough idea about how FEAST works to tune correctly the input parameters
 * and get fast convergence). 
 *
 * The convergence test is done on the maximum residual, see:
 *   "FEAST as a subspace iteration eigensolver accelerated by approximate
 *   spectral projection" by Tang & Polizzi, arXiv:1302:0432v4.
 *
 * This is momentarily implemented only for double. I promise I'll make it generic.
 */

feast_status
feast(const feast_eigensolver_params<double>& params,
    const priv::spmat<double>& A,
    const priv::spmat<double>& B,
    priv::mat<double>& eigvecs,
    priv::vec<double>& eigvals)
{
    using complex = std::complex<double>;
    using rdv = priv::vec<double>;
    using rdm = priv::mat<double>;
    using cdm = priv::mat<complex>;
    using csm = priv::spmat<complex>;
    using eigsolver = Eigen::SelfAdjointEigenSolver<rdm>;
    using generalized_eigsolver = Eigen::GeneralizedSelfAdjointEigenSolver<rdm>;
    int eigvals_found = 0;
    const double *xs = quadrature_xs;
    const double *omegas = quadrature_ws;

    bool A_square = A.rows() == A.cols();
    bool B_square = B.rows() == B.cols();
    bool A_B_same_size = (A.rows() == B.rows()) and (A.cols() == B.cols());

    if ( (not A_square) or (not B_square) or (not A_B_same_size) ) {
        if (params.verbose) {
            std::cout << "FEAST: invalid input" << std::endl;
        }
        return feast_status::invalid_input;
    }

    auto N = A.rows();
    auto M0 = params.subspace_size;
    auto r = (params.max_eigval - params.min_eigval)/2.0;

    rdm Y = rdm::Random(N, M0);
    double trace_prev = 0.0;
    for (size_t iter = 0; iter < params.max_iter; iter++) {        
        cdm Yc = Y;
        rdm Q = rdm::Zero(N, M0);
        
        /* Subspace projection */
        for (size_t e = 0; e < 8; e++) {
            auto x_e = xs[e];
            auto theta_e = -(M_PI/2.)*(x_e-1.);
            auto mid = (params.max_eigval + params.min_eigval)/2.0;
            auto Z_e = mid + r*std::exp( std::complex<double>(0.0, theta_e) );
            csm lhs = Z_e*B - std::complex<double>(1.0, 0.0)*A;
            
            cdm Qe;
            switch(params.fis) {
                case feast_inner_solver::eigen_sparselu: {
                    Eigen::SparseLU<csm> solver;
                    solver.compute(lhs);
                    if(solver.info() != Eigen::Success) {
                        std::cout << "FEAST: SparseLU failed" << std::endl;
                        return feast_status::inner_solver_problem;
                    }
                    Qe = solver.solve(Yc);
                } break;

                case feast_inner_solver::mumps: {
#ifdef HAVE_MUMPS
                    Qe = mumps_lu(lhs, Yc);
#else
                    throw std::runtime_error("Mumps is not installed");
#endif /* HAVE_MUMPS */
                }
            }

            cdm T = r * std::exp(std::complex<double>(0.0, theta_e)) * Qe;
            auto omega_e = omegas[e];
            Q = Q - (omega_e/2.0)*T.real();
        }

        /* Form matrices for reduced problem */
        rdm Aq = Q.transpose() * A * Q;
        rdm Bq = Q.transpose() * B * Q;

        /* Estimate the number of eigenvalues and check subspace size */
        if (iter == 1) {
            eigsolver Bq_es(Bq, Eigen::EigenvaluesOnly);
            rdv Bq_eigvals = Bq_es.eigenvalues();
            const double thresh = 0.25; /* See Tang & Polizzi 2014. */

            size_t Bq_eigvals_abovetr = 0;
            double Bq_min_eigval = Bq_eigvals(0);
            for (size_t i = 0; i < Bq_eigvals.size(); i++) {
                Bq_min_eigval = std::min(Bq_min_eigval, Bq_eigvals(i));
                if ( Bq_eigvals(i) >= thresh )
                    Bq_eigvals_abovetr++;
            }

            if (Bq_min_eigval >= thresh) {
                if (params.verbose) {
                    std::cout << "FEAST: subspace too small" << std::endl;
                }   
                
                return feast_status::subspace_too_small;
            }
        }
        
        /* Solve the reduced eigenvalue problem */
        generalized_eigsolver es(Aq,Bq);
        rdv epsilon = es.eigenvalues();
        rdm Phi = es.eigenvectors();

        /* Compute eigenvectors of original problem */
        rdm X = Q * Phi;

        size_t found_eigs = 0;
        for (size_t m = 0; m < epsilon.size(); m++)
            if (epsilon(m) >= params.min_eigval and epsilon(m) <= params.max_eigval)
                found_eigs++;

        eigvecs = rdm::Zero(X.rows(), found_eigs);
        eigvals = rdv::Zero(found_eigs);

        size_t m_out = 0;
        for (size_t m = 0; m < epsilon.size(); m++) {
            if (epsilon(m) >= params.min_eigval and epsilon(m) <= params.max_eigval) {
                eigvecs.col(m_out) = X.col(m);
                eigvals(m_out) = epsilon(m);
                m_out++;
            }
        }
        
        /* Compute residual and eigenvalue trace to check convergence */
        auto rmax = std::abs(std::max(params.min_eigval, params.max_eigval));
        rdm BX = B * eigvecs;
        rdm residual = A * eigvecs - BX * eigvals.asDiagonal();
        
        double maxres = 0.0;
        for (size_t i = 0; i < residual.cols(); i++) {
            auto num = residual.col(i).norm();
            auto den = rmax*BX.col(i).norm();
            maxres = std::max(maxres, num/den);
        }

        double trace = eigvals.sum();
        double trerr = std::abs(trace - trace_prev)/rmax;
        trace_prev = trace;
        if (params.verbose) {
            std::cout << "Iteration: " << iter+1 << "/" << params.max_iter << ", ";
            std::cout << "Eigvals found: " << found_eigs << ", ";
            std::cout << "Trace: " << trace << ", Trace error: " << trerr;
            std::cout << ", maximum residual: " << maxres << std::endl;
        }

        if (maxres < std::pow(10,-params.tolerance)) {
            eigvals_found = eigvals.size();
            if (params.verbose) {
                std::cout << "FEAST: success" << std::endl;
            }
            return feast_status::success;
        }

        Y = B * X;
    }

    return feast_status::did_not_converge;
}

} // namespace disk::solvers