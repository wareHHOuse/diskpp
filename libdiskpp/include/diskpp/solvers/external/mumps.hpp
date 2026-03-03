/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *  
 * Matteo Cicuttin (C) 2020-2022
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

#ifdef HAVE_MUMPSS
#include <smumps_c.h>
#endif

#ifdef HAVE_MUMPSD
#include <dmumps_c.h>
#endif

#ifdef HAVE_MUMPSC
#include <cmumps_c.h>
#endif

#ifdef HAVE_MUMPSZ
#include <zmumps_c.h>
#endif

#define MUMPS_JOB_INIT              -1
#define MUMPS_JOB_TERMINATE         -2
#define MUMPS_JOB_ANALYZE            1
#define MUMPS_JOB_FACTORIZE          2
#define MUMPS_JOB_SOLVE              3
#define MUMPS_JOB_ANALYZE_FACTORIZE  4
#define MUMPS_JOB_FACTORIZE_SOLVE    5
#define MUMPS_JOB_EVERYTHING         6

#define MUMPS_OUTPUT_ERROR           1
#define MUMPS_OUTPUT_DIAG            2
#define MUMPS_OUTPUT_GLOBAL          4

/* This is typical Fortran madness: the
 * magic number -987654 is really required
 * to inform the mumps code to use the
 * MPI communicator MPI_COMM_WORLD
 */
#define FORTRAN_MADNESS_MAGIC       -987654

namespace mumps_priv {

template<typename T>
struct mumps_types;

#ifdef HAVE_MUMPSS
template<>
struct mumps_types<float> {
    typedef SMUMPS_STRUC_C          MUMPS_STRUC_C;
    typedef float                   mumps_elem_t;
};

void
call_mumps(typename mumps_types<float>::MUMPS_STRUC_C *id)
{
    smumps_c(id);
}
#endif

#ifdef HAVE_MUMPSD
template<>
struct mumps_types<double> {
    typedef DMUMPS_STRUC_C          MUMPS_STRUC_C;
    typedef double                  mumps_elem_t;
};

void
call_mumps(typename mumps_types<double>::MUMPS_STRUC_C *id)
{
    dmumps_c(id);
}
#endif

#ifdef HAVE_MUMPSC
template<>
struct mumps_types<std::complex<float>> {
    typedef CMUMPS_STRUC_C          MUMPS_STRUC_C;
    typedef mumps_complex           mumps_elem_t;
};

void
call_mumps(typename mumps_types<std::complex<float>>::MUMPS_STRUC_C *id)
{
    cmumps_c(id);
}
#endif

#ifdef HAVE_MUMPSZ
template<>
struct mumps_types<std::complex<double>> {
    typedef ZMUMPS_STRUC_C          MUMPS_STRUC_C;
    typedef mumps_double_complex    mumps_elem_t;
};

void
call_mumps(typename mumps_types<std::complex<double>>::MUMPS_STRUC_C *id)
{
    zmumps_c(id);
}
#endif

template<typename From>
typename mumps_types<From>::mumps_elem_t *
mumps_cast_from(From *from)
{
    return reinterpret_cast<typename mumps_types<From>::mumps_elem_t *>(from);
}

} // namespace mumps_priv

template<typename MUMPS_STRUC>
struct mumps_state
{
    MUMPS_STRUC         id;
    size_t              mflops;
    /* From the MUMPS documentation it is not clear which
     * should be the lifetime of Aii and Aja */
    std::vector<int>    Aii, Aja;
};

template<typename T>
class mumps_solver
{
    using MUMPS_STRUC = typename mumps_priv::mumps_types<T>::MUMPS_STRUC_C;
    using mstate = mumps_state<MUMPS_STRUC>;
    using mstate_uptr = std::unique_ptr<mstate>;

    mstate_uptr state;
    int         symmetric_flag;
    int         parallel_flag;

public:
    mumps_solver()
    {
        state = std::make_unique<mstate>();
        state->id.job = MUMPS_JOB_INIT;
        state->id.par = 1;
        state->id.sym = 0;
        state->id.comm_fortran = FORTRAN_MADNESS_MAGIC;
        
        mumps_priv::call_mumps(&state->id);
        
        state->id.icntl[0]= -1;//6;    // Suppress error output
        state->id.icntl[1]= -1;//6;    // Suppress diagnostic output
        state->id.icntl[2]= -1;//6;    // Suppress global output
        state->id.icntl[3]= 2;         // Loglevel
    }
    
    ~mumps_solver()
    {
        state->id.job = MUMPS_JOB_TERMINATE;
        mumps_priv::call_mumps(&state->id);
    }
    
    template<int _Options, typename _Index>
    void
    factorize(Eigen::SparseMatrix<T, _Options, _Index>& A) const
    {
        static_assert( !(_Options & Eigen::RowMajor), "CSR not supported yet.");
        if ( A.rows() != A.cols() )
            throw std::invalid_argument("Only square matrices");

        A.makeCompressed();

        int     N       = A.rows();
        T *     data    = A.valuePtr();
        int *   ia      = A.outerIndexPtr();
        int     js      = A.nonZeros();
        int *   ja      = A.innerIndexPtr();
        int     is      = A.innerSize();

        state->Aii.resize( A.nonZeros() );

        /* Convert CSC to COO */
        for (int i = 0; i < is; i++) {
            int begin = ia[i];
            int end = ia[i+1];

            for (size_t j = begin; j < end; j++)
                state->Aii[j] = i+1;
        }

        state->Aja.resize( A.nonZeros() );
        for (size_t i = 0; i < js; i++)
            state->Aja[i] = ja[i] + 1;
        
        state->id.a = mumps_priv::mumps_cast_from<T>(data);
        state->id.irn = state->Aii.data();
        state->id.jcn = state->Aja.data();
        state->id.n = A.rows();
        state->id.nz = A.nonZeros();

        state->id.job = MUMPS_JOB_ANALYZE_FACTORIZE;
        mumps_priv::call_mumps(&state->id);
        state->mflops = (size_t)((state->id.rinfog[0] + state->id.rinfog[1])/1e6);
    }

    /* Pseudo-compatibility with Eigen solver interface */
    template<int _Options, typename _Index>
    void
    compute(Eigen::SparseMatrix<T, _Options, _Index>& A) const
    {
        factorize(A);
    }

    template<int nrhs>
    Eigen::Matrix<T, Eigen::Dynamic, nrhs>
    solve(Eigen::Matrix<T, Eigen::Dynamic, nrhs>& b) const
    {
        Eigen::Matrix<T, Eigen::Dynamic, nrhs> ret = b;
        T* x = ret.data();

        state->id.nrhs = (nrhs > 0) ? nrhs : b.cols();
        state->id.rhs = mumps_priv::mumps_cast_from<T>(ret.data());
        
        state->id.job = MUMPS_JOB_SOLVE;
        mumps_priv::call_mumps(&state->id);

        return ret;
    }
    
    void set_output(int oflags)
    {
        state->id.icntl[0]= -1;//6     // Suppress error output
        state->id.icntl[1]= -1;//6     // Suppress diagnostic output
        state->id.icntl[2]= -1;//6;    // Suppress global output
        
        if (oflags & MUMPS_OUTPUT_ERROR)
            state->id.icntl[0] = 6;
        
        if (oflags & MUMPS_OUTPUT_DIAG)
            state->id.icntl[1] = 6;
        
        if (oflags & MUMPS_OUTPUT_GLOBAL)
            state->id.icntl[2] = 6;
    }
    
    int symmetric() const
    {
        return state->id.sym;
    }
    
    void symmetric(int flag)
    {
        state->id.sym = flag;
    }
    
    int parallel() const
    {
        return state->id.par;
    }
    
    void parallel(int flag)
    {
        state->id.par = flag;
    }

    size_t get_Mflops() const
    {
        return state->mflops;
    }
};

namespace mumps_priv {

template<bool symmetric, typename T, int _Options, typename _Index, _Index nrhs>
Eigen::Matrix<T, Eigen::Dynamic, nrhs>
mumps(Eigen::SparseMatrix<T, _Options, _Index>& A, Eigen::Matrix<T, Eigen::Dynamic, nrhs>& b)
{
    using MUMPS_STRUC = typename mumps_priv::mumps_types<T>::MUMPS_STRUC_C;
    
    MUMPS_STRUC         id;
    
    id.job = MUMPS_JOB_INIT;
    id.par = 1;

    if constexpr (symmetric)
        id.sym = 1;
    else
        id.sym = 0;

    id.comm_fortran = FORTRAN_MADNESS_MAGIC;

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

    id.icntl[0]= -1;//6;    // Suppress error output
    id.icntl[1]= -1;//6;    // Suppress diagnostic output
    id.icntl[2]= -1;//6;    // Suppress global output
    id.icntl[3]= 2;         // Loglevel
    
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ret = b;
    T* x = ret.data();
    
    id.nrhs = (nrhs > 0) ? nrhs : b.cols();
    id.lrhs = b.rows();
    id.rhs = mumps_priv::mumps_cast_from<T>(ret.data());

    id.job = MUMPS_JOB_EVERYTHING;
    mumps_priv::call_mumps(&id);
    
    for (size_t i = 0; i < js; i++)
        ja[i] -= 1.;
    
    id.job = MUMPS_JOB_TERMINATE;
    mumps_priv::call_mumps(&id);
    
    return ret;
}

} //namespace mumps_priv

template<typename T, int _Options, typename _Index, _Index nrhs>
Eigen::Matrix<T, Eigen::Dynamic, nrhs>
mumps_ldlt(Eigen::SparseMatrix<T, _Options, _Index>& A, Eigen::Matrix<T, Eigen::Dynamic, nrhs>& b)
{
    return mumps_priv::mumps<true>(A,b);
}

template<typename T, int _Options, typename _Index, _Index nrhs>
Eigen::Matrix<T, Eigen::Dynamic, nrhs>
mumps_lu(Eigen::SparseMatrix<T, _Options, _Index>& A, Eigen::Matrix<T, Eigen::Dynamic, nrhs>& b)
{
    return mumps_priv::mumps<false>(A,b);
}

#endif /* HAVE_MUMPS */
