
#include <iostream>
#include <vector>
#include <array>

extern "C" {

void feastinit(int *);

void dfeast_scsrgv(const char& uplo, const int& n,
                   const double * a, const int * ia, const int * ja,
                   const double * b, const int * ib, const int * jb,
                   int * fpm, double& epsout, int& loop,
                   const double& emin, const double& emax, int& m0,
                   double * e, double * x, int& m, double * res, int& info);

}


int main(void)
{

    std::vector<double> A({ 2, -1, -1, 2, -1, -1, 2, -1, -1, 2});
    std::vector<double> B({ 4, 1, 1, 4, 1, 1, 4, 1, 1, 4});
    std::vector<int> ia({1,3,6,9,11});
    std::vector<int> ja({1,2,1,2,3,2,3,4,3,4});
    std::array<int, 128> fpm;

    for (auto& a : A)
        a /= 5;
    
    for (auto& b : B)
        b *= 5./6.;

    std::vector<double> eigvals, eigvecs, res;

    double epsout, emin = 0.1, emax = 100;
    int M0=3, M, loop, info;

    eigvals.resize(M0);
    eigvecs.resize(M0*4);
    res.resize(M0);

    feastinit(fpm.data());

    dfeast_scsrgv('F', 4, A.data(), ia.data(), ja.data(),
                          B.data(), ia.data(), ja.data(), fpm.data(), epsout,
                          loop, emin, emax, M0,
                          eigvals.data(), eigvecs.data(), M, res.data(), info);

    for (auto& ev : eigvals)
        std::cout << ev << " ";
    std::cout << std::endl;

    std::cout << M << " " << info << std::endl;

    return 0;
    
}