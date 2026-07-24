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

#ifdef HAVE_MUMPS
#include "external/mumps.hpp"
#endif

#include "diskpp/common/eigen.hpp"

namespace disk::solvers {

enum class direct_solver {
    autosel,
#ifdef HAVE_MUMPS
    mumps,
#endif
    sparselu,
#ifdef HAVE_PARDISO
    pardiso
#endif
};

enum class direct_solver_status {
    ok,
    failure,
};

template<typename T, int _Options, typename _Index, _Index nrhs>
direct_solver_status
sparse_lu(Eigen::SparseMatrix<T, _Options, _Index>& A,
    Eigen::Matrix<T, Eigen::Dynamic, nrhs>& b,
    Eigen::Matrix<T, Eigen::Dynamic, nrhs>& x,
    direct_solver solver = direct_solver::autosel)
{
    using spmat_t = Eigen::SparseMatrix<T, _Options, _Index>;
    switch (solver) {
    /* Beware here: if we have mumps, autosel falls back to mumps, otherwise
     * it falls back to SparseLU (which is the only we are guaranteed to have)
     */

    case direct_solver::autosel:

#ifdef HAVE_MUMPS
    case direct_solver::mumps: {
        mumps_solver<T> solver;
        solver.compute(A);
        x = solver.solve(b);
        if (solver.failure()) {
            return direct_solver_status::failure;
        }
    } break;
#endif

    case direct_solver::sparselu: {
        Eigen::SparseLU<spmat_t> LU(A);
        x = LU.solve(b);
        if (LU.info() != Eigen::Success) {
            return direct_solver_status::failure;
        }
    } break;

#ifdef HAVE_PARDISO
    case direct_solver::pardiso: {
        Eigen::PardisoLU<spmat_t> LU(A);
        x = LU.solve(b);
        if (LU.info() != Eigen::Success) {
            return direct_solver_status::failure;
        }
    } break;
#endif

    default:
        std::cerr << "sparse_lu(): inexistent solver selected" << std::endl;
        return  direct_solver_status::failure;
    } // switch

   return  direct_solver_status::ok;
}

template<typename T, int _Options, typename _Index, _Index nrhs>
direct_solver_status
sparse_ldlt(Eigen::SparseMatrix<T, _Options, _Index>& A,
    Eigen::Matrix<T, Eigen::Dynamic, nrhs>& b,
    Eigen::Matrix<T, Eigen::Dynamic, nrhs>& x,
    direct_solver solver = direct_solver::autosel)
{
    using spmat_t = Eigen::SparseMatrix<T, _Options, _Index>;
    switch (solver) {
    /* Beware here: if we have mumps, autosel falls back to mumps, otherwise
     * it falls back to SparseLU (which is the only we are guaranteed to have)
     */

    case direct_solver::autosel:

#ifdef HAVE_MUMPS
    case direct_solver::mumps: {
        mumps_solver<T> solver;
        solver.symmetric(true);
        solver.compute(A);
        x = solver.solve(b);
        if (solver.failure()) {
            return direct_solver_status::failure;
        }
    } break;
#endif

    case direct_solver::sparselu: {
        Eigen::SimplicialLDLT<spmat_t> LDLT(A);
        x = LDLT.solve(b);
        if (LDLT.info() != Eigen::Success) {
            return direct_solver_status::failure;
        }
    } break;

#ifdef HAVE_PARDISO
    case direct_solver::pardiso: {
        Eigen::PardisoLDLT<spmat_t> LDLT(A);
        x = LDLT.solve(b);
        if (LDLT.info() != Eigen::Success) {
            return direct_solver_status::failure;
        }
    } break;
#endif

    default:
        std::cerr << "sparse_ldlt(): inexistent solver selected" << std::endl;
        return  direct_solver_status::failure;
    } // switch

   return  direct_solver_status::ok;
}

template < typename T, int _Options = Eigen::ColMajor, typename _Index = int >
class sparse_solver {
  private:
    using spmat_t = Eigen::SparseMatrix< T, _Options, _Index >;
    template < _Index nrhs >
    using dense_matrix_t = Eigen::Matrix< T, Eigen::Dynamic, nrhs >;
    direct_solver _solver_type;
    bool _isFacto;

#ifdef HAVE_MUMPS
    mumps_solver< T > _mumps_solver;
#endif
    Eigen::SparseLU< spmat_t > _sparselu_solver;

#ifdef HAVE_PARDISO
    Eigen::PardisoLU< spmat_t > _pardisolu_solver;
#endif

  public:
    explicit sparse_solver( direct_solver solver = direct_solver::autosel )
        : _solver_type( solver ), _isFacto( false ) {
        /*
         * Automatic solver selection:
         *
         * 1. MUMPS, if available;
         * 2. PARDISO, if available;
         * 3. Eigen::SparseLU otherwise.
         */
        if ( _solver_type == direct_solver::autosel ) {
#if defined( HAVE_MUMPS )
            _solver_type = direct_solver::mumps;
#elif defined( HAVE_PARDISO )
            _solver_type = direct_solver::pardiso;
#else
            _solver_type = direct_solver::sparselu;
#endif
        }

        /*
         * Check that the explicitly selected solver is available in the
         * current build.
         */
        switch ( _solver_type ) {
        case direct_solver::mumps:
#ifndef HAVE_MUMPS
            throw std::invalid_argument( "sparse_solver: MUMPS is not available in this build" );
#endif
            break;
        case direct_solver::pardiso:
#ifndef HAVE_PARDISO
            throw std::invalid_argument( "sparse_solver: PARDISO is not available in this build" );
#endif
            break;
        case direct_solver::sparselu:
            // Eigen::SparseLU is always available.
            break;
        case direct_solver::autosel:
            /*
             * This should be unreachable because autosel has already been
             * resolved above.
             */
            throw std::logic_error( "sparse_solver: automatic solver selection failed" );
        default:
            throw std::invalid_argument( "sparse_solver: unknown direct solver" );
        }
    }

    [[nodiscard]]
    direct_solver solver_type() const noexcept {
        return _solver_type;
    }

    [[nodiscard]]
    bool is_factorized() const noexcept {
        return _isFacto;
    }

    void reset() noexcept { _isFacto = false; }

    direct_solver_status factorize( spmat_t &A ) {
        _isFacto = false;
        switch ( _solver_type ) {
#ifdef HAVE_MUMPS
        case direct_solver::mumps:
            _mumps_solver.compute( A );
            if ( _mumps_solver.failure() ) {
                return direct_solver_status::failure;
            }
            break;
#endif
        case direct_solver::sparselu:
            _sparselu_solver.compute( A );
            if ( _sparselu_solver.info() != Eigen::Success ) {
                return direct_solver_status::failure;
            }
            break;

#ifdef HAVE_PARDISO
        case direct_solver::pardiso:
            _pardisolu_solver.compute( A );
            if ( _pardisolu_solver.info() != Eigen::Success ) {
                return direct_solver_status::failure;
            }
            break;
#endif
        case direct_solver::autosel:
            throw std::logic_error( "sparse_solver::factorize: unresolved automatic selection" );
        default:
            throw std::invalid_argument( "sparse_solver::factorize: unknown solver" );
        }

        _isFacto = true;
        return direct_solver_status::ok;
    }

    template < _Index nrhs >
    direct_solver_status solve( const dense_matrix_t< nrhs > &b, dense_matrix_t< nrhs > &x ) {
        if ( !_isFacto ) {
            throw std::runtime_error( "Factorize matrix before solving." );
        }

        switch ( _solver_type ) {
#ifdef HAVE_MUMPS
        case direct_solver::mumps:
            x = _mumps_solver.solve( b );
            if ( _mumps_solver.failure() ) {
                return direct_solver_status::failure;
            }
            break;
#endif
        case direct_solver::sparselu:
            x = _sparselu_solver.solve( b );
            if ( _sparselu_solver.info() != Eigen::Success ) {
                return direct_solver_status::failure;
            }
            break;

#ifdef HAVE_PARDISO
        case direct_solver::pardiso:
            x = _pardisolu_solver.solve( b );
            if ( _pardisolu_solver.info() != Eigen::Success ) {
                return direct_solver_status::failure;
            }
            break;
#endif
        case direct_solver::autosel:
            throw std::logic_error( "sparse_solver::solve: unresolved automatic selection" );
        default:
            throw std::invalid_argument( "sparse_solver::solve: unknown solver" );
        }
        return direct_solver_status::ok;
    }

    template < _Index nrhs >
    [[nodiscard]]
    std::pair< direct_solver_status, dense_matrix_t< nrhs > >
    solve( const dense_matrix_t< nrhs > &b ) {
        dense_matrix_t< nrhs > x;
        const auto status = solve( b, x );
        return { status, std::move( x ) };
    }

    template < _Index nrhs >
    [[nodiscard]]
    std::pair< direct_solver_status, dense_matrix_t< nrhs > >
    solve( const spmat_t &A, const dense_matrix_t< nrhs > &b ) {
        const auto factorisation_status = factorize( A );
        if ( factorisation_status != direct_solver_status::ok ) {
            return { factorisation_status, dense_matrix_t< nrhs > {} };
        }
        return solve( b );
    }

    template < _Index nrhs >
    direct_solver_status solve( spmat_t &A, const dense_matrix_t< nrhs > &b,
                                dense_matrix_t< nrhs > &x ) {
        const auto factorisation_status = factorize( A );
        if ( factorisation_status != direct_solver_status::ok ) {
            return factorisation_status;
        }
        return solve( b, x );
    }
};
};