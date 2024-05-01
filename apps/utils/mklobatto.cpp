#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <Eigen/Dense>

using matrix_t = Eigen::MatrixXd;
using vector_t = Eigen::VectorXd;

void gen_lobatto(size_t npts)
{
    vector_t x = vector_t::Zero(npts);
    
    /* Initial guess for nodes: Chebyshev interpolation nodes */
    for (size_t i = 0; i < npts; i++)
        x(i) = -std::cos(M_PI*i/(npts-1));

    matrix_t P = matrix_t::Zero(npts, npts);

    double tol = 1;
    while (tol > 1e-15) {
        vector_t x_old = x;
        /* Generate Legendre polynomials */
        P.col(0) = vector_t::Ones(npts);
        P.col(1) = x;
        for (size_t j = 2; j < npts; j++) 
            P.col(j) = ((2*j-1)*x.cwiseProduct(P.col(j-1)) - (j-1)*P.col(j-2))/j;

        vector_t err = (x.cwiseProduct(P.col(npts-1)) - P.col(npts-2));
        err.array() /= (npts*P.col(npts-1)).array();
        x = x_old - err;
        tol = (x-x_old).lpNorm<Eigen::Infinity>();
        //std::cout << x_old.transpose() << std::endl;
        //std::cout << tol << std::endl;
    }
    //std::cout << x.transpose() << std::endl;
    vector_t w = 2.0*(((npts*(npts-1))*P.col(npts-1).array().square()).cwiseInverse());
    //std::cout << w.transpose() << std::endl;
    //std::cout << w.sum() << std::endl;

    for (size_t i = 0; i < npts; i++)
        printf("    { %16.16f, %16.16f }%c\n", x[i], w[i], (i < npts-1) ? ',' : ' ');
}

int main(void)
{
    for (size_t i = 2; i < 21; i++) {

        printf("struct lobatto_point lobatto_rule_%zu[] = {\n", i);
        gen_lobatto(i);
        printf("};\n");
    }
}