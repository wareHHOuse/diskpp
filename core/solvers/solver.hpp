/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
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

#include "contrib/sol2/sol.hpp"
#include "common/eigen.hpp"

namespace disk { namespace solvers {

void
init_lua(sol::state& lua)
{
    lua["solver"] = lua.create_table();
    lua["solver"]["cg"] = lua.create_table();
}

template<typename T>
bool
conjugated_gradient(sol::state& lua,
                    const Eigen::SparseMatrix<T>& A,
                    const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
                    Eigen::Matrix<T, Eigen::Dynamic, 1>& x)
{
    size_t                      N = A.cols();
    size_t                      iter = 0;
    T                           nr, nr0;
    T                           alpha, beta, rho;
    
    Eigen::Matrix<T, Eigen::Dynamic, 1> d(N), r(N), r0(N), y(N);
    x = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(N);
 

    r0 = d = r = b - A*x;
    nr = nr0 = r.norm();

    size_t max_iter = lua["solver"]["cg"]["max_iter"].get_or(1000);
    double rr_tol = lua["sover"]["cg"]["rr_tol"].get_or(1e-8);
    double rr_max = lua["solver"]["cg"]["rr_max"].get_or(20.0);
    
    //std::ofstream ofs("cocg_nopre_convergence.txt");
    
    while ( nr/nr0 > rr_tol && iter < max_iter && nr/nr0 < rr_max )
    {
        std::cout << "                                                 \r";
        std::cout << " -> Iteration " << iter << ", rr = ";
        std::cout << nr/nr0 << "\b\r";
        std::cout.flush();
        
        //ofs << nr/nr0 << std::endl;
        y = A*d;
        rho = r.dot(r);
        alpha = rho/d.dot(y);
        x = x + alpha * d;
        r = r - alpha * y;
        beta = r.dot(r)/rho;
        d = r + beta * d;
        
        nr = r.norm();
        iter++;
    }
    
    //ofs << nr/nr0 << std::endl;
    //ofs.close();
    
    std::cout << " -> Iteration " << iter << ", rr = " << nr/nr0 << std::endl;
    
    return true;
}

template<typename T>
bool
mkl_pardiso(sol::state& lua,
            const Eigen::SparseMatrix<T>& A,
            const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
            Eigen::Matrix<T, Eigen::Dynamic, 1>& x)
{
    Eigen::PardisoLU<Eigen::SparseMatrix<T>>  solver;
    //solver.pardisoParameterArray()[59] = 0; //out-of-core

    solver.analyzePattern(A);
    solver.factorize(A);
    x = solver.solve(b);

    return true;
}

template<typename T>
bool
linear_solver(sol::state& lua,
              const Eigen::SparseMatrix<T>& A,
              const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
              Eigen::Matrix<T, Eigen::Dynamic, 1>& x)
{
    std::string solver_type = lua["solver"]["solver_type"].get_or(std::string("cg"));

    if ( solver_type == "cg" )
        return conjugated_gradient(lua, A, b, x);
#ifdef HAVE_INTEL_MKL
    else if ( solver_type == "pardiso" )
        return mkl_pardiso(lua, A, b, x);
#endif
        
//        Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;


    return false;
}


} //namespace solvers
} //namespace disk