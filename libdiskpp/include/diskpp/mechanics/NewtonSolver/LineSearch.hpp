/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework.
 * M. Abbas, A. Ern, N. Pignet.
 * International Journal of Numerical Methods in Engineering (2019)
 * 120(3), 303-327
 * DOI: 10.1002/nme.6137
 */

#pragma once

#include <vector>

namespace disk {

namespace mechanics {

/**
 * @brief LineSearch
 *
 */

// Référence of Aitken and Anderson algorithm.
// Isabelle Ramière, Thomas Helfer. Iterative residual-based vector methods to accelerate fixed
// point iterations.
// Computers & Mathematics with Applications, 2015, 70, pp.2210 - 2226.
// 10.1016/j.camwa.2015.08.025. cea-01403292

template < typename T >
class ConvergenceAcceleration {

    typedef dynamic_vector< T > vector_type;
    typedef dynamic_matrix< T > matrix_type;

    int n_iter;

    enum kiter {
        k = 1,
        km = 0,
    };
    std::vector< vector_type > G, Xa;

  public:
    ConvergenceAcceleration() : n_iter( 0 ) {
        G.resize( 2 );
        Xa.resize( 2 );
    }

    vector_type aitken( const vector_type &G_k ) {
        // Also called crossed secand method - eq 44.
        if ( n_iter == 0 ) {
            n_iter++;
            Xa[km] = G_k;
            return G_k;
        } else if ( n_iter == 1 ) {
            n_iter++;
            G[km] = G_k;
            Xa[k] = G_k;
            return G_k;
        } else {
            n_iter++;
            G[k] = G_k;
            const auto dG_k = G[k] - G[km];

            const auto dX_k = G[k] - Xa[k];
            const auto dX_km = G[km] - Xa[km];
            const auto ddX = dX_k - dX_km;

            // compute acceleration
            const T wr = ddX.dot( dG_k ) / ddX.squaredNorm();

            // compute accelerted solution
            const vector_type Xa_kp = G_k - wr * dX_k;

            // update
            G[km] = G[k];

            Xa[km] = Xa[k];
            Xa[k] = Xa_kp;

            return Xa_kp;
        }
    }

    vector_type relaxation( const vector_type &X_kp, const T omega = 0.5 ) {
        if ( n_iter == 0 ) {
            n_iter++;
            Xa[k] = X_kp;
            return X_kp;
        } else {
            n_iter++;

            // compute accelerted solution
            const vector_type Xa_kp = ( 1.0 - omega ) * Xa[k] + omega * X_kp;

            // update
            Xa[k] = Xa_kp;

            return Xa_kp;
        }
    }

    template < typename Func >
    void secant( const Func &func, const double ALF = 0.1, const int MAXIT = 10 ) {

        const T TOLX = std::numeric_limits< T >::epsilon();
        T p0, p1, f1, f0, rho_0, rho_1;
        T rho, rho_neg, rho_pos, rho_opt, rho_new, rho_cur;
        T f, f_opt, f_cur;
        bool b_pos;

        // fixed paramters - from code_aster
        const T rho_min = 1e-2, rho_max = 10., rho_excl = 0.9e-2;
        const T parmul = 3.0;

        // Doc:
        // https://codeaster.pages.pleiade.edf.fr/doc/docaster/manuals/man_r/r5/r5.03.01/Recherche_lin_aire.html
        // METHODE="SECANT"

        // Compute residual.dot(increment) (and update solution)
        // auto _f = [&fvec, &func, &dx, &xold, &n, &x]( const T &rho ) {
        //     for ( int j = 0; j < n; j++ )
        //         G[j] = xold[j] + rho * dG[j];
        //     const auto norm = func( x );

        //     T f = 0.0;
        //     for ( int j = 0; j < n; j++ )
        //         f += fvec[j] * dG[j];
        //     return f;
        // };

        // project bound on admissible interval
        auto _proj = [rho_min, rho_max, rho_excl]( T &rho ) {
            const T rho_tmp = rho;
            if ( rho_tmp < rho_min ) {
                rho = rho_min;
            }
            if ( rho_tmp > rho_max ) {
                rho = rho_max;
            }
            if ( rho_tmp < 0.0 && rho_tmp >= -rho_excl ) {
                rho = -rho_excl;
            }
            if ( rho_tmp >= 0 && rho_tmp <= rho_excl ) {
                rho = rho_excl;
            }
        };

        // initial values
        const T f_old = func( 0.0 );
        const T f_cvg = ALF * std::abs( f_old );
        const T sens = ( f_old <= 0.0 ) ? 1.0 : -1.0;

        rho_opt = 1.0, rho = sens * 1.0;
        rho_neg = 0.0, rho_pos = std::numeric_limits< T >::signaling_NaN();
        f_opt = 10e100;

        rho_0 = 0.0, rho_1 = rho_0;
        f0 = sens * f_old, f1 = f0;

        b_pos = false;

        for ( int its = 0; its < MAXIT; its++ ) {
            // Compute new residual
            try {
                f = func( rho );
            } catch ( ... ) {
                break;
            }

            rho_cur = sens * rho;
            f_cur = sens * f;

            // Store value
            rho_0 = rho_1, f0 = f1;
            rho_1 = rho_cur, f1 = f_cur;

            // Update bounds
            if ( f_cur < 0.0 ) {
                rho_neg = rho_cur;
            } else {
                b_pos = true;
                rho_pos = rho_cur;
            }

            // Optimal solution until now ?
            if ( std::abs( f_cur ) < std::abs( f_opt ) ) {
                rho_opt = rho_cur;
                f_opt = f_cur;
                _proj( rho_opt );
            }

            // Search maximal bound
            if ( b_pos ) {
                if ( std::abs( f1 ) >= std::abs( f0 ) ) {
                    // f is not decreased - use dichotomie
                    rho_new = 0.5 * ( rho_neg + rho_pos );
                } else {
                    // linear interpolation
                    if ( std::abs( rho_1 - rho_0 ) > TOLX ) {
                        p1 = ( f1 - f0 ) / ( rho_1 - rho_0 );
                        p0 = f0 - p1 * rho_0;

                        if ( std::abs( p1 ) <= std::abs( f0 ) / ( rho_pos + rho_0 ) ) {
                            rho_new = 0.5 * ( rho_neg + rho_pos );
                        } else {
                            rho_new = -p0 / p1;
                        }
                    } else {
                        // failed
                        break;
                    }
                }
            } else {
                rho_new = parmul * rho_cur;
            }

            // minimal bound
            if ( rho_new < rho_neg ) {
                if ( b_pos ) {
                    rho_new = 0.5 * ( rho_neg + rho_pos );
                } else {
                    // failed
                    break;
                }
            }

            // maximal bound
            if ( b_pos && rho_new > rho_pos ) {
                rho_new = 0.5 * ( rho_neg + rho_pos );
            }

            // project bound
            _proj( rho_new );

            // update
            rho = sens * rho_new;

            // Test convergence ?
            if ( std::abs( f_opt ) <= f_cvg ) {
                break;
            }
        }

        /* Return optimal value */
        // std::cout << "rho_opt: " << rho_opt << std::endl;
        f = func( rho_opt, false );
    }

    vector_type anderson( const vector_type &G_k ) {
        // also called alternate secant method - eq.45
        // special version for M = 1
        if ( n_iter == 0 ) {
            n_iter++;
            Xa[km] = G_k;
            return G_k;
        } else if ( n_iter == 1 ) {
            n_iter++;
            G[km] = G_k;
            Xa[k] = G_k;
            return G_k;
        }
        n_iter++;

        G[k] = G_k;

        const auto dX_k = G[k] - Xa[k];
        const auto dX_km = G[km] - Xa[km];
        const auto ddX = dX_k - dX_km;

        // compute acceleration
        const T wr = ddX.dot( dX_k ) / ddX.squaredNorm();

        // compute accelerted solution
        const vector_type Xa_kp = ( 1.0 - wr ) * G[k] + wr * G[km];

        // update
        G[km] = G[k];

        Xa[km] = Xa[k];
        Xa[k] = Xa_kp;

        return Xa_kp;
    }

    vector_type anderson( const vector_type &G_k, const int M ) {
        // anderson method with M values.
        if ( M < 1 ) {
            return G_k;
        } else if ( M == 1 ) {
            return anderson( G_k );
        };
        if ( n_iter == 0 ) {
            Xa.resize( M + 1 );
            G.resize( M + 1 );
            n_iter++;
            Xa[km] = G_k;
            return G_k;
        } else if ( n_iter == 1 ) {
            n_iter++;
            G[km] = G_k;
            Xa[k] = G_k;
            return G_k;
        }

        const int n = G_k.size();
        const int m_k = std::min( M, n_iter - 1 );

        G[m_k] = G_k;

        matrix_type ddX( n, m_k );
        matrix_type dG( n, m_k );

        for ( int i = 1; i <= m_k; i++ ) {
            const auto dX_i = G[i] - Xa[i];
            const auto dX_im = G[i - 1] - Xa[i - 1];
            ddX.col( i - 1 ) = dX_i - dX_im;
            dG.col( i - 1 ) = G[i] - G[i - 1];
        }

        const vector_type dX_k = G[m_k] - Xa[m_k];

        // Solve least-square problem
        const FullPivHouseholderQR< matrix_type > qr( ddX );
        const vector_type gamma = qr.solve( dX_k );

        // std::cout << "gamma: " << gamma.rows() << ", " << gamma.cols() << std::endl;
        // std::cout << "dG: " << dG.rows() << ", " << dG.cols() << std::endl;

        const vector_type Xa_kp = G[m_k] - dG * gamma;

        // Update solution
        if ( m_k < M ) {
            Xa[n_iter] = Xa_kp;
        } else {
            for ( int i = 0; i < m_k; i++ ) {
                G[i] = G[i + 1];
                Xa[i] = Xa[i + 1];
            }
            Xa[m_k] = Xa_kp;
        }
        n_iter++;

        return Xa_kp;
    }
};
} // namespace mechanics
} // namespace disk