/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *  
 * Matteo Cicuttin (C) 2020,2021
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group  
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

#ifdef HAVE_MUMPS

#include <complex>

#include <smumps_c.h>
#include <dmumps_c.h>
#include <cmumps_c.h>
#include <zmumps_c.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#define MUMPS_JOB_INIT              -1
#define MUMPS_JOB_TERMINATE         -2
#define MUMPS_JOB_ANALYZE            1
#define MUMPS_JOB_FACTORIZE          2
#define MUMPS_JOB_SOLVE              3
#define MUMPS_JOB_ANALYZE_FACTORIZE  4
#define MUMPS_JOB_FACTORIZE_SOLVE    5
#define MUMPS_JOB_EVERYTHING         6

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

void
call_mumps(typename mumps_types<float>::MUMPS_STRUC_C *id)
{
    smumps_c(id);
}

void
call_mumps(typename mumps_types<double>::MUMPS_STRUC_C *id)
{
    dmumps_c(id);
}

void
call_mumps(typename mumps_types<std::complex<float>>::MUMPS_STRUC_C *id)
{
    cmumps_c(id);
}

void
call_mumps(typename mumps_types<std::complex<double>>::MUMPS_STRUC_C *id)
{
    zmumps_c(id);
}

} // namespace mumps_priv

template<typename T>
class mumps_solver
{
    using MUMPS_STRUC = typename mumps_priv::mumps_types<T>::MUMPS_STRUC_C;
    
    MUMPS_STRUC         id;
    size_t              solver_Mflops;
    std::vector<int>    Aii, Aja;

public:
    mumps_solver()
    {
        id.job = MUMPS_JOB_INIT;
        id.par = 1;
        id.sym = 0;
        id.comm_fortran = -987654; /* <- This is typical Fortran madness: the
                                    magic number -987654 is really required
                                    to inform the mumps code to use the
                                    MPI communicator MPI_COMM_WORLD
                                    */

        mumps_priv::call_mumps(&id);
    }
    
    ~mumps_solver()
    {
        id.job = MUMPS_JOB_TERMINATE;
        mumps_priv::call_mumps(&id);
    }
    
    template<int _Options, typename _Index>
    void
    factorize(Eigen::SparseMatrix<T, _Options, _Index>& A)
    {
        if ( A.rows() != A.cols() )
            throw std::invalid_argument("Only square matrices");

        static_assert( !(_Options & Eigen::RowMajor), "CSR not supported yet.");
        
        A.makeCompressed();

        int     N       = A.rows();
        T *     data    = A.valuePtr();
        int *   ia      = A.outerIndexPtr();
        int     js      = A.nonZeros();
        int *   ja      = A.innerIndexPtr();
        int     is      = A.innerSize();

        Aii.resize( A.nonZeros() );

        /* Convert CSC to COO */
        for (int i = 0; i < is; i++)
        {
            int begin = ia[i];
            int end = ia[i+1];

            for (size_t j = begin; j < end; j++)
                Aii[j] = i+1;
        }

        Aja.resize( A.nonZeros() );
        for (size_t i = 0; i < js; i++)
            Aja[i] = ja[i] + 1;
        
        id.a = mumps_priv::mumps_cast_from<T>(data);
        id.irn = Aii.data();
        id.jcn = Aja.data();
        id.n = A.rows();
        id.nz = A.nonZeros();
        
        id.icntl[0]= -1;//6     // Suppress error output
        id.icntl[1]= -1;//6     // Suppress diagnostic output
        id.icntl[2]= -1;//6;    // Suppress global output
        id.icntl[3]= 2;         // Loglevel

        id.job = MUMPS_JOB_ANALYZE_FACTORIZE;
        mumps_priv::call_mumps(&id);
        solver_Mflops = (size_t)((id.rinfog[0] + id.rinfog[1])/1e6);
    }

    Eigen::Matrix<T, Eigen::Dynamic, 1>
    solve(Eigen::Matrix<T, Eigen::Dynamic, 1>& b)
    {
        Eigen::Matrix<T, Eigen::Dynamic, 1> ret = b;
        T* x = ret.data();
        
        id.rhs = mumps_priv::mumps_cast_from<T>(ret.data());
        
        id.job = MUMPS_JOB_SOLVE;
        mumps_priv::call_mumps(&id);

        return ret;
    }

    size_t get_Mflops()
    {
        return solver_Mflops;
    }
};

template<typename T, int _Options, typename _Index>
Eigen::Matrix<T, Eigen::Dynamic, 1>
mumps(Eigen::SparseMatrix<T, _Options, _Index>& A, Eigen::Matrix<T, Eigen::Dynamic, 1>& b)
{
    using MUMPS_STRUC = typename mumps_priv::mumps_types<T>::MUMPS_STRUC_C;
    
    MUMPS_STRUC         id;
    
    id.job = MUMPS_JOB_INIT;
    id.par = 1;
    id.sym = 0;
    id.comm_fortran = -987654; /* <- This is typical Fortran madness: the
                                magic number -987654 is really required
                                to inform the mumps code to use the
                                MPI communicator MPI_COMM_WORLD
                                */

    mumps_priv::call_mumps(&id);
    
    if ( A.rows() != A.cols() )
        throw std::invalid_argument("Only square matrices");

    static_assert( !(_Options & Eigen::RowMajor), "CSR not supported yet.");
    
    A.makeCompressed();

    int     N       = A.rows();
    T *     data    = A.valuePtr();
    int *   ia      = A.outerIndexPtr();
    int     js      = A.nonZeros();
    int *   ja      = A.innerIndexPtr();
    int     is      = A.innerSize();

    std::vector<int>    Aii;
    Aii.resize( A.nonZeros() );

    /* Convert CSC to COO */
    for (int i = 0; i < is; i++)
    {
        int begin = ia[i];
        int end = ia[i+1];

        for (size_t j = begin; j < end; j++)
            Aii[j] = i+1;
    }

    for (size_t i = 0; i < js; i++)
        ja[i] += 1.;
    
    id.a = mumps_priv::mumps_cast_from<T>(data);
    id.irn = Aii.data();
    id.jcn = ja;
    id.n = A.rows();
    id.nz = A.nonZeros();

    id.icntl[0]= -1;//6     // Suppress error output
    id.icntl[1]= -1;//6     // Suppress diagnostic output
    id.icntl[2]= -1;//6;    // Suppress global output
    id.icntl[3]= 2;         // Loglevel
    
    Eigen::Matrix<T, Eigen::Dynamic, 1> ret = b;
    T* x = ret.data();
    
    id.rhs = mumps_priv::mumps_cast_from<T>(ret.data());

    id.job = MUMPS_JOB_EVERYTHING;
    mumps_priv::call_mumps(&id);
    
    for (size_t i = 0; i < js; i++)
        ja[i] -= 1.;
    
    id.job = MUMPS_JOB_TERMINATE;
    mumps_priv::call_mumps(&id);
    
    return ret;
}

#endif /* HAVE_MUMPS */
