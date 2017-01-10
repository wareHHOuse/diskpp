/*
 * Copyright (C) 2013-2017, Matteo Cicuttin - matteo.cicuttin@uniud.it
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

#include <complex>

#include <smumps_c.h>
#include <dmumps_c.h>
#include <cmumps_c.h>
#include <zmumps_c.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace mumps_priv {

template<typename T>
struct mumps_types;

template<>
struct mumps_types<float> {
    typedef SMUMPS_STRUC_C          MUMPS_STRUC_C;
    typedef float                   mumps_elem_t;
};

template<>
struct mumps_types<double> {
    typedef DMUMPS_STRUC_C          MUMPS_STRUC_C;
    typedef double                  mumps_elem_t;
};

template<>
struct mumps_types<std::complex<float>> {
    typedef CMUMPS_STRUC_C          MUMPS_STRUC_C;
    typedef mumps_complex           mumps_elem_t;
};

template<>
struct mumps_types<std::complex<double>> {
    typedef ZMUMPS_STRUC_C          MUMPS_STRUC_C;
    typedef mumps_double_complex    mumps_elem_t;
};

template<typename From>
typename mumps_types<From>::mumps_elem_t *
mumps_cast_from(From *from)
{
    return reinterpret_cast<typename mumps_types<From>::mumps_elem_t *>(from);
}

template<typename T>
void
call_mumps(T *);

template<>
void
call_mumps(typename mumps_types<float>::MUMPS_STRUC_C *id)
{
    smumps_c(id);
}

template<>
void
call_mumps(typename mumps_types<double>::MUMPS_STRUC_C *id)
{
    dmumps_c(id);
}

template<>
void
call_mumps(typename mumps_types<std::complex<float>>::MUMPS_STRUC_C *id)
{
    cmumps_c(id);
}

template<>
void
call_mumps(typename mumps_types<std::complex<double>>::MUMPS_STRUC_C *id)
{
    zmumps_c(id);
}

} // namespace mumps_priv

template<typename T>
class mumps_solver
{

public:
    mumps_solver()
    {}

    template<int _Options, typename _Index>
    Eigen::Matrix<T, 1, Eigen::Dynamic>
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

        std::vector<int>    ii( A.nonZeros() );

        size_t vpos = 0;
        for (int i = 0; i < is; i++)
        {
            int begin = ia[i];
            int end = ia[i+1];

            for (size_t j = begin; j < end; j++)
                ii[j] = i+1;
        }

        for (size_t i = 0; i < js; i++)
            ja[i] += 1;

        Eigen::Matrix<T, Eigen::Dynamic, 1> ret = b;
        T * x = ret.data();

        typename mumps_priv::mumps_types<T>::MUMPS_STRUC_C   id;

        id.job = -1;
        id.par = 1;
        id.sym = 0;
        id.comm_fortran = -987654; /* <- This is typical Fortran madness: the
                                    magic number -987654 is really required
                                    to inform the mumps code to use the
                                    MPI communicator MPI_COMM_WORLD
                                    */

        mumps_priv::call_mumps(&id);

        id.a = mumps_priv::mumps_cast_from<T>(data);
        id.irn = ii.data();
        id.jcn = ja;
        id.n = A.rows();
        id.nz = A.nonZeros();
        id.rhs = mumps_priv::mumps_cast_from<T>(ret.data());


        //id.icntl[0] = -1; id.icntl[1] = -1; id.icntl[2] = -1; id.icntl[3] = 0;
        id.icntl[0] = 1; id.icntl[1] = 1; id.icntl[2] = 6; id.icntl[3] = 3;

        id.job = 6;
        mumps_priv::call_mumps(&id);

        id.job = -2;
        mumps_priv::call_mumps(&id);

        for (size_t i = 0; i < js; i++)
            ja[i] -= 1;

        return ret;
    }
};
