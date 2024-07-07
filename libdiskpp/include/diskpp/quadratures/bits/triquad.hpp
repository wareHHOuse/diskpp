#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <numeric> 

/*
This is a translation to C++ from the matlab code:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% triquad.m - Gaussian Quadrature for a triangular domain
%
% This scripts computes the N^2 nodes and weights for a triangle with
% vertices given by the 3x2 vector v. The nodes are produced by collapsing
% the square to a triangle. 
%
% Sample Usage: 
%
% >>[X,Y,Wx,Wy]=triquad(8,[0 0; 0 2; 2 1])
% >>f=@(x,y) exp(x+y);
% >>Q=Wx'*feval(f,X,Y)*Wy;
%
% Reference:  J.N. Lyness, Ronald Cools, A Survey of Numerical Cubature
%             over Triangles (1994)
%             http://citeseer.ist.psu.edu/lyness94survey.html
%
% Written by: Greg von Winckel
% Contact: gregvw(at)math(dot)unm(dot)edu
% http://math.unm.edu/~gregvw
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

namespace disk {
namespace quadrature {


using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;


matrix_t
compute_AB(size_t N)
{
    vector_t A = vector_t::Zero(N+1);
    A(0) =  1.0/3.0;
    for(size_t i = 1; i <= N; i++)
    {
        size_t nnk = 2 * i + 1;
        A[i] =  1.0 / double(nnk * (nnk + 2.0));
    }

    vector_t B = vector_t::Zero(N + 1);
    B(0) = 2.0;
    B(1) = 2.0/9.0;
    for(size_t n = 2; n <= N; n++)
    {   
        size_t nk = n + 1;
        size_t nnk = 2 * n + 1;
        size_t nnk2 = nnk * nnk;
    
        B[n] =  4.0 * (n * nk)* (n * nk) / (nnk2 * (nnk2 - 1.0));
    }
   
    matrix_t AB(N + 1, 2);
    AB.col(0) = A;
    AB.col(1) = B;

    return AB; 
}



matrix_t
compute_M(size_t N, const matrix_t& ab)
{
    vector_t s = ab.block(1, 1, N, 1).cwiseSqrt();

    matrix_t  mat = matrix_t::Zero(N, N);
    mat.diagonal() = ab.block(0,0,N,1);

    for(size_t n = 0; n < N-1; n++)
    {
        mat(1+n,n) =  s[n];
        mat(n,n+1) =  s[n];
    }

    return mat;
}



std::pair<std::pair<vector_t, matrix_t>, std::vector<size_t>>
compute_eigensolver(const matrix_t& mat)
{
    size_t N = mat.rows();
    Eigen::SelfAdjointEigenSolver<matrix_t> eigensolver(mat);
    vector_t eigen_vals = eigensolver.eigenvalues().real();
    matrix_t eigen_vecs = eigensolver.eigenvectors().real();

    std::vector<size_t> I (N);
    std::iota(I.begin(), I.end(), 0);
    std::sort(I.begin(), I.end(), 
               [&eigen_vals](size_t i0, size_t i1) 
                { return eigen_vals[i0] < eigen_vals[i1];}
              );

    vector_t sorted_evals = vector_t::Zero(N);
    for(int i = 0; i < sorted_evals.size(); ++i)
        sorted_evals[i] = eigen_vals[I[i]];

    auto E = std::make_pair(sorted_evals, eigen_vecs);
    return std::make_pair(E, I);
}


std::pair<vector_t, vector_t>
compute_y(size_t N)
{
    double eps = 2.2204e-16;
    vector_t ones = vector_t::Ones(N);
    vector_t y    = vector_t::Zero(N);
    for(size_t i = 0; i < N; i++)
        y[i] = std::cos( (2.0 * (N -1 - i) + 1.0) * M_PI/(2.0 * (N-1) + 2.0));

    matrix_t L = matrix_t::Zero(N, N+1);
    vector_t Lp = vector_t::Zero(N);

    L.col(0) = ones;    

    size_t iter = 0;
    double error = 1.0;
    while (error > eps)    
    {
        L.col(1) = y;   
        for (size_t k = 2; k < N+1; k++)
            L.col(k) = ((2.0 * k - 1.0) * y.cwiseProduct(L.col(k-1))
                                 - (k-1.0) * L.col(k-2))/double(k);

        vector_t dividend = L.col(N-1) - y.cwiseProduct(L.col(N));
        vector_t divisor = ones.array() - y.array().square();
        Lp = (N + 1)* dividend.cwiseQuotient(divisor);   
        
        vector_t diff = - (L.col(N)).cwiseQuotient(Lp); 
        y += diff; 
        error = diff.maxCoeff();
        iter = iter + 1;
    }


    return std::make_pair(y, Lp);
}



matrix_t
quad_weights(size_t N, const vector_t& y,
                const vector_t& Lp,
                const std::vector<size_t>& I,
                const matrix_t& AB,
                const matrix_t& Evecs,
                const matrix_t& vertices)
{
    // coords   | 1., 0., 0.| |x0 y0|
    //          |-1., 0., 1.| |x1 y1|
    //          | 0., 1.,-1.| |x2 y2|

    double det23 = std::abs((vertices(2,0) - vertices(0,0))
                          * (vertices(1,1) - vertices(2,1)) 
                          - (vertices(1,0) - vertices(2,0))
                          * (vertices(2,1) - vertices(0,1)));

    double coeff = double((N+1) * (N+1)) /double(N*N);

    vector_t ones = vector_t::Ones(N);
    vector_t divisor = (ones.array() - y.array().square()) * Lp.array().square();

    vector_t wx = vector_t::Zero(N); 
    for(size_t i = 0; i <N; i++)
        wx[i] = 0.25 * AB(0,1) * Evecs(0, I[i]) *  Evecs(0, I[i]);

    vector_t Wx = det23 * wx;  
    vector_t Wy = coeff * ones.array() /  divisor.array();
    matrix_t qw = Wx * Wy.transpose(); 

    return qw;
}



std::pair<matrix_t, matrix_t>
quad_points(size_t N, const vector_t& x, const vector_t& y,
            const matrix_t&  v)
{
    vector_t t  = 0.5 * (y + vector_t::Ones(N));
    
    matrix_t xx = matrix_t::Zero(N,N);
    matrix_t yy = matrix_t::Zero(N,N);

    for(size_t i = 0; i < x.size(); i++)
        xx.col(i) = x;
    
    for(size_t i = 0; i < x.size(); i++)
        yy.col(i) = t[i] * x;

    // coords   | 1., 0., 0.| |vx0 vy0| = |   vx0      vy0  |
    //          |-1., 0., 1.| |vx1 vy1|   |vx2-vx0   vy2-vy0|
    //          | 0., 1.,-1.| |vx2 vy2|   |vx1-vx2   vy1-vy2|

    matrix_t Xc = v(0,0) * matrix_t::Ones(N,N);
    matrix_t Yc = v(0,1) * matrix_t::Ones(N,N);

    Xc += (v(2,0)-v(0,0)) * xx + (v(1,0)-v(2,0)) * yy;
    Yc += (v(2,1)-v(0,1)) * xx + (v(1,1)-v(2,1)) * yy;  

    return std::make_pair(Xc, Yc);
}


std::tuple<vector_t, vector_t, vector_t>
triquad(size_t degree, const matrix_t& vertices)
{
    auto AB = compute_AB(degree);
    auto M  = compute_M( degree, AB);
    auto [E ,I] = compute_eigensolver(M);
    auto [Evals, Evecs] = E;

    auto x  = 0.5 * (Evals + vector_t::Ones(degree));
    auto [y, Lp] = compute_y(degree);

    auto qw  = quad_weights(degree, y, Lp, I, AB, Evecs, vertices);
    auto [X, Y] = quad_points(degree, x, y, vertices);

    return std::make_tuple( vector_t {X.reshaped()},
                            vector_t {Y.reshaped()},
                            vector_t {qw.reshaped()});
}


} //namespace quadrature
} //namespace disk