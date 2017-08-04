/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016-2017 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

/*
 * Copyright (C) 2013-2016, Matteo Cicuttin - matteo.cicuttin@uniud.it
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




#include <Eigen/Dense>
#include <Eigen/Sparse>

#define CALL_FORTRAN(function) function##_

namespace disk { namespace solvers {

namespace agmg_priv {

extern "C" void CALL_FORTRAN(sagmg)(int&, float *, int *, int *,
                                    float *, float *,
                                    int&, int&, int&, int&, double&);

extern "C" void CALL_FORTRAN(dagmg)(int&, double *, int *, int *,
                                    double *, double *,
                                    int&, int&, int&, int&, double&);


/* Standard mandates that std::complex must be compatible with C99 format:
 * array T[2] with re [0] and im [1]. This seems to be compatible with fortran,
 * so it is assumed that these two are safe.
 */
extern "C" void CALL_FORTRAN(cagmg)(int&, std::complex<float> *, int *, int *,
                                    std::complex<float> *, std::complex<float> *,
                                    int&, int&, int&, int&, double&);

extern "C" void CALL_FORTRAN(zagmg)(int&, std::complex<double> *, int *, int *,
                                    std::complex<double> *, std::complex<double> *,
                                    int&, int&, int&, int&, double&);


template<typename T>
void
    call_agmg(int& N, T *a, int *ja, int *ia, T *f, T *x,
              int& ijob, int& iprint, int& nrest, int& iter, double& tol);

template<>
void
call_agmg<float>(int& N, float *a, int *ja, int *ia, float *f, float *x,
                 int& ijob, int& iprint, int& nrest, int& iter, double& tol)
{
    CALL_FORTRAN(sagmg)(N, a, ja, ia, f, x, ijob, iprint, nrest, iter, tol);
}

template<>
void
call_agmg<double>(int& N, double *a, int *ja, int *ia, double *f, double *x,
                  int& ijob, int& iprint, int& nrest, int& iter, double& tol)
{
    CALL_FORTRAN(dagmg)(N, a, ja, ia, f, x, ijob, iprint, nrest, iter, tol);
}

template<>
void
call_agmg<std::complex<float>>(int& N, std::complex<float> *a, int *ja,
                               int *ia, std::complex<float> *f,
                               std::complex<float> *x, int& ijob,
                               int& iprint, int& nrest, int& iter, double& tol)
{
    CALL_FORTRAN(cagmg)(N, a, ja, ia, f, x, ijob, iprint, nrest, iter, tol);
}

template<>
void
call_agmg<std::complex<double>>(int& N, std::complex<double> *a, int *ja,
                                int *ia, std::complex<double> *f,
                                std::complex<double> *x, int& ijob,
                                int& iprint, int& nrest, int& iter, double& tol)
{
    CALL_FORTRAN(zagmg)(N, a, ja, ia, f, x, ijob, iprint, nrest, iter, tol);
}

} // namespace agmg_priv


template<typename T>
struct agmg_solver_config
{
    int         ijob;
    int         iprint;
    int         num_restarts;
    int         max_iters;
    double      tolerance;

    agmg_solver_config() : ijob(0),
                           iprint(6),
                           num_restarts(1),
                           max_iters(10000),
                           tolerance(1e-8)
    {}
};

template<typename T>
class agmg_solver
{
    agmg_solver_config<T>   config;

public:
    agmg_solver()
    {}

    agmg_solver(const agmg_solver_config<T>& conf) : config(conf) {}

    template<int _Options, typename _Index>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    solve(Eigen::SparseMatrix<T, _Options, _Index>& A,
          Eigen::Matrix<T, Eigen::Dynamic, 1>& b)
    {
        if ( A.rows() != A.cols() )
            throw std::invalid_argument("Only square matrices");

        A.makeCompressed();

        int     N       = A.rows();
        T *     data    = A.valuePtr();
        int *   ia      = A.outerIndexPtr();
        int     js      = A.nonZeros();
        int *   ja      = A.innerIndexPtr();
        int     is      = A.innerSize();
        T *     f       = b.data();

        /* Convert to one-based */
        for (size_t i = 0; i < is+1; i++)
            ia[i] += 1;

        for (size_t i = 0; i < js; i++)
            ja[i] += 1;

        assert(data != nullptr);
        assert(ia != nullptr);
        assert(ja != nullptr);

        Eigen::Matrix<T, Eigen::Dynamic, 1> ret;
        ret.resize(N);
        T * x = ret.data();

        if (_Options == Eigen::ColMajor)
        {
            // User asked to use the transpose, but backend is CSC: don't ask
            // AGMG to compute with the transpose.
            if (config.ijob > 100)
                config.ijob -= 100;

            // User did not ask to use the transpose, but backend is CSC: ask
            // AGMG to NOT compute with the transpose.
            if (config.ijob < 100)
                config.ijob += 100;
        }

        agmg_priv::call_agmg<T>(N, data, ja, ia, f, x, config.ijob, config.iprint,
                                config.num_restarts, config.max_iters, config.tolerance);

        /* Convert back to zero-based */
        for (size_t i = 0; i < is+1; i++)
            ia[i] -= 1;

        for (size_t i = 0; i < js; i++)
            ja[i] -= 1;

        return ret;
    }
};


template<typename T>
bool
agmg_multigrid_solver(sol::state& lua,
                      Eigen::SparseMatrix<T>& A,
                      Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
                      Eigen::Matrix<T, Eigen::Dynamic, 1>& x)
{
    agmg_solver<T> agmg;

    x = agmg.solve(A, b);

    return true;
}

}} // namespaces