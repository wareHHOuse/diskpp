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
#include <vector>
#include <complex>
#include <array>
#include <type_traits>
#include <fstream>

#include "diskpp/common/eigen.hpp"
#include "diskpp/common/colormanip.h"
#include "diskpp/common/timecounter.hpp"

/****************************************************************************/
// feast prototypes
extern "C" {

void feastinit(int *);

void sfeast_scsrev(const char& uplo, const int& n,
                   const float * a, const int * ia, const int * ja,
                   int * fpm, float& epsout, int& loop,
                   const float& emin, const float& emax, int& m0,
                   float * e, float * x, int& m, float * res, int& info);

void dfeast_scsrev(const char& uplo, const int& n,
                   const double * a, const int * ia, const int * ja,
                   int * fpm, double& epsout, int& loop,
                   const double& emin, const double& emax, int& m0,
                   double * e, double * x, int& m, double * res, int& info);

void cfeast_hcsrev(const char& uplo, const int& n,
                   const std::complex<float> * a, const int * ia, const int * ja,
                   int * fpm, float& epsout, int& loop,
                   const float& emin, const float& emax, int& m0,
                   float * e, std::complex<float> * x, int& m,
                   float * res, int& info);

void zfeast_hcsrev(const char& uplo, const int& n,
                   const std::complex<double> * a, const int * ia, const int * ja,
                   int * fpm, double& epsout, int& loop,
                   const double& emin, const double& emax, int& m0,
                   double * e, std::complex<double> * x, int& m,
                   double * res, int& info);

void sfeast_scsrgv(const char& uplo, const int& n,
                   const float * a, const int * ia, const int * ja,
                   const float * b, const int * ib, const int * jb,
                   int * fpm, float& epsout, int& loop,
                   const float& emin, const float& emax, int& m0,
                   float * e, float * x, int& m, float * res, int& info);

void dfeast_scsrgv(const char& uplo, const int& n,
                   const double * a, const int * ia, const int * ja,
                   const double * b, const int * ib, const int * jb,
                   int * fpm, double& epsout, int& loop,
                   const double& emin, const double& emax, int& m0,
                   double * e, double * x, int& m, double * res, int& info);

void cfeast_hcsrgv(const char& uplo, const int& n,
                   const std::complex<float> * a, const int * ia, const int * ja,
                   const std::complex<float> * b, const int * ib, const int * jb,
                   int * fpm, float& epsout, int& loop,
                   const float& emin, const float& emax, int& m0,
                   float * e, std::complex<float> * x, int& m,
                   float * res, int& info);

void zfeast_hcsrgv(const char& uplo, const int& n,
                   const std::complex<double> * a, const int * ia, const int * ja,
                   const std::complex<double> * b, const int * ib, const int * jb,
                   int * fpm, double& epsout, int& loop,
                   const double& emin, const double& emax, int& m0,
                   double * e, std::complex<double> * x, int& m,
                   double * res, int& info);

} // extern "C"

#define FEASTPARM_LEN   128

namespace disk{

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
    T       min_eigval, max_eigval;
    int     subspace_size;
    int     eigvals_found;
    int     feast_info;
    size_t  max_iter;
    feast_inner_solver fis;
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
feast(feast_eigensolver_params<double>& params,
    const Eigen::SparseMatrix<double>& A,
    const Eigen::SparseMatrix<double>& B,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& eigvecs,
    Eigen::Matrix<double, Eigen::Dynamic, 1>& eigvals)
{
    using rdv = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using rdm = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    using cdm = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;
    using csm = Eigen::SparseMatrix<std::complex<double>>;
    using eigsolver = Eigen::SelfAdjointEigenSolver<rdm>;
    using generalized_eigsolver = Eigen::GeneralizedSelfAdjointEigenSolver<rdm>;

    const double *xs = quadrature_xs;
    const double *omegas = quadrature_ws;

    bool A_square = A.rows() == A.cols();
    bool B_square = B.rows() == B.cols();
    bool A_B_same_size = (A.rows() == B.rows()) and (A.cols() == B.cols());

    if ( (not A_square) or (not B_square) or (not A_B_same_size) )
        return feast_status::invalid_input;

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
                    Qe = mumps_lu(lhs, Yc);
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

            if (Bq_min_eigval >= thresh)
                return feast_status::subspace_too_small;
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
            params.eigvals_found =  eigvals.size();
            return feast_status::success;
        }

        Y = B * X;
    }

    return feast_status::did_not_converge;
}

template<int _Options, typename _Index>
int
generalized_eigenvalue_solver(feast_eigensolver_params<double>& params,
                              Eigen::SparseMatrix<double, _Options, _Index>& L,
                              Eigen::SparseMatrix<double, _Options, _Index>& R,
                              Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& eigvecs,
                              Eigen::Matrix<double, Eigen::Dynamic, 1>& eigvals)
{
    if ( L.rows() != R.rows() && L.cols() != R.cols() && L.rows() == L.cols() )
    {
        std::cout << "Two square matrices of the same size are needed." << std::endl;
        return false;
    }

    std::array<int, FEASTPARM_LEN>    fpm;
    feastinit(fpm.data());

    auto Lc = eigen_sparse_raw(L, true);
    auto Rc = eigen_sparse_raw(R, true);
    int N = Lc.n;

    if (params.verbose)
        fpm.at(0) = 1;

    if (params.tolerance > 0 && params.tolerance < 16)
        fpm.at(2) = params.tolerance;
    else
    {
        std::cout << "Invalid tolerance." << std::endl;
        return false;
    }

    if (params.subspace_size < 1 || params.subspace_size > N)
    {
        std::cout << "Invalid subspace size." << std::endl;
        return false;
    }

    double  eps;
    int     loop;

    Eigen::Matrix<double, Eigen::Dynamic, 1> res;
    res     = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(params.subspace_size);
    eigvals = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(params.subspace_size);
    eigvecs = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(N, params.subspace_size);

    dfeast_scsrgv('F', N, Lc.data, Lc.ia, Lc.ja, Rc.data, Rc.ia, Rc.ja,
                       fpm.data(), eps, loop,
                       params.min_eigval,
                       params.max_eigval,
                       params.subspace_size,
                       eigvals.data(),
                       eigvecs.data(),
                       params.eigvals_found,
                       res.data(),
                       params.feast_info);

    return params.feast_info;
}

template<typename T>
bool
setup_feast(sol::state& lua, feast_eigensolver_params<T>& fep)
{
    fep.verbose     = lua["solver"]["feast"]["verbose"].get_or(false);
    fep.tolerance   = lua["solver"]["feast"]["tolerance"].get_or(9);

    auto min_ev = lua["solver"]["feast"]["min_eigval"];
    if (!min_ev.valid())
    {
        std::cout << "solver.feast.min_eigval not set." << std::endl;
        return false;
    }
    fep.min_eigval = min_ev;

    auto max_ev = lua["solver"]["feast"]["max_eigval"];
    if (!max_ev.valid())
    {
        std::cout << "solver.feast.max_eigval not set." << std::endl;
        return false;
    }
    fep.max_eigval = max_ev;

    auto subsp = lua["solver"]["feast"]["subspace_size"];
    if (!subsp.valid())
    {
        std::cout << "solver.feast.subspace_size not set." << std::endl;
        return false;
    }
    fep.subspace_size = subsp;

    return true;
}



/*
bool test_eigenvalue_solver(void)
{
    int N = 256;

    int aN = N-2;
    double h = 1./N;
    Eigen::SparseMatrix<double> K(aN, aN);
    Eigen::SparseMatrix<double> M(aN, aN);

    typedef Eigen::Triplet<double> triplet_type;

    std::vector<triplet_type> tK, tM;

    for (size_t i = 0; i < aN; i++)
    {
        if (i > 0)
        {
            tK.push_back( triplet_type(i-1, i, -1.0) );
            tM.push_back( triplet_type(i-1, i,  1.0/6.0) );
        }

        tK.push_back( triplet_type(i, i, 2.0) );
        tM.push_back( triplet_type(i, i, 4.0/6.0) );

        if (i < aN-1)
        {
            tK.push_back( triplet_type(i+1, i, -1.0) );
            tM.push_back( triplet_type(i+1, i,  1.0/6.0) );
        }
    }

    K.setFromTriplets(tK.begin(), tK.end());
    M.setFromTriplets(tM.begin(), tM.end());

    K *= 1./(h*h);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigvecs;
    Eigen::Matrix<double, Eigen::Dynamic, 1> eigvals;
    
    disk::feast_eigensolver_params<double> fep;
    
    fep.verbose = true;
    fep.tolerance = 9;
    fep.min_eigval = 1;
    fep.max_eigval = 30;
    fep.subspace_size = 10;

    generalized_eigenvalue_solver(fep, K, M, eigvecs, eigvals);

    std::cout << eigvals << std::endl;

    return true;
}

*/

} // end disk
