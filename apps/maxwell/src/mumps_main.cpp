#include <iostream>
#include <vector>

//clang++ -std=c++17 -I/usr/local//Cellar/brewsci-mumps/5.3.5/include/ -I/usr/local/include/eigen3 -L/usr/local//Cellar/brewsci-mumps/5.3.5/lib main.cpp -lsmumps -lcmumps -ldmumps -lzmumps -lmpiseq

#define HAVE_MUMPS
#include "mumps.hpp"

int main(void)
{
    using T = double;
    size_t ms = 10;
    
    using trip = Eigen::Triplet<T>;
    
    mumps_solver<T> solver;
    
    Eigen::SparseMatrix<T>              A(ms, ms);
    Eigen::Vector<T, Eigen::Dynamic>    b(ms), x(ms);
    std::vector<trip>                   triplets;
    
    for (size_t i = 0; i < ms; i++)
    {
        triplets.push_back( trip(i, i, i+1) );
        b(i) = 1;
    }
    
    A.setFromTriplets(triplets.begin(), triplets.end());
    
    solver.factorize(A);
    x = solver.solve(b);
    
    std::cout << x.transpose() << std::endl;
    
    x = mumps(A,b);
    
    std::cout << x.transpose() << std::endl;
    
    std::cout << A << std::endl;
    
    return 0;
}

