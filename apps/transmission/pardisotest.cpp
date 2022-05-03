#include <iostream>
#include <fstream>
#include <vector>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/hho"
#include "diskpp/output/silo.hpp"

#include "sol/sol.hpp"
#include "mumps.hpp"

#define PARDISO_IN_CORE                 0
#define PARDISO_OUT_OF_CORE_IF_NEEDED   1
#define PARDISO_OUT_OF_CORE_ALWAYS      2

template<typename T>
struct pardiso_params
{
    bool    report_factorization_Mflops;
    int     out_of_core; //0: IC, 1: IC, OOC if limits passed, 2: OOC
    int     mflops;

    pardiso_params() :
        report_factorization_Mflops(false), out_of_core(0), mflops(0)
    {}
};

template<typename T>
bool
mkl_pardiso(pardiso_params<T>& params,
            const Eigen::SparseMatrix<T>& A,
            const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
            Eigen::Matrix<T, Eigen::Dynamic, 1>& x)
{
    Eigen::PardisoLU<Eigen::SparseMatrix<T>>  solver;

    if (params.out_of_core >= 0 && params.out_of_core <= 2)
        solver.pardisoParameterArray()[59] = params.out_of_core;

    if (params.report_factorization_Mflops)
        solver.pardisoParameterArray()[18] = -1; //report flops

    solver.analyzePattern(A);
    if (solver.info() != Eigen::Success) {
       std::cerr << "ERROR: analyzePattern failed" << std::endl;
       return false;
    }

    solver.factorize(A);
    if (solver.info() != Eigen::Success) {
       std::cerr << "ERROR: Could not factorize the matrix" << std::endl;
       std::cerr << "Try to tweak MKL_PARDISO_OOC_MAX_CORE_SIZE" << std::endl;
       return false;
    }

    x = solver.solve(b);
    if (solver.info() != Eigen::Success) {
       std::cerr << "ERROR: Could not solve the linear system" << std::endl;
       return false;
    }

    if (params.report_factorization_Mflops)
    {
        int mflops = solver.pardisoParameterArray()[18];
        params.mflops = mflops;
        std::cout << "[PARDISO] Factorization Mflops: " << mflops << std::endl;
    }

    return true;
}

int main(void)
{
	std::cout << "PARDISOTEST START" << std::endl;
	using T = double;
	using triplet = Eigen::Triplet<T>;
    std::vector< triplet > triplets;

    disk::generic_mesh<T,2> msh;
    disk::fvca5_mesh_loader<T,2> loader;

	
	//disk::simplicial_mesh<T,2> msh2;
    //disk::gmsh_geometry_loader< disk::simplicial_mesh<T,2> > loader2;
	
	gmsh::initialize(0,nullptr);

	std::ifstream ifs("trmat.txt");

	size_t rows, cols, nnz;
	ifs >> rows >> cols >> nnz;

	std::cout << rows << " " << cols << " " << nnz << std::endl;

	while (!ifs.eof())
	{
		int i, j;
		double val;
		ifs >> i >> j >> val;

		triplets.push_back( triplet(i, j, val) );
	}

	Eigen::SparseMatrix<T> LHS(rows, cols);
	LHS.setFromTriplets( triplets.begin(), triplets.end() );
	triplets.clear();

	std::ifstream ifsv("trvec.txt");
	size_t sz;
	ifsv >> sz;

	disk::dynamic_vector<T> RHS = disk::dynamic_vector<T>::Zero(sz);
	for (size_t i = 0; i < sz; i++)
		ifsv >> RHS(i);


    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(rows);

	std::cout << "Running pardiso" << std::endl;
    pardiso_params<T> pparams;
    pparams.report_factorization_Mflops = true;
    pparams.out_of_core = PARDISO_OUT_OF_CORE_IF_NEEDED;
	//bool success = mkl_pardiso(pparams, LHS, RHS, sol);

    sol = mumps_lu(LHS, RHS);

	std::ofstream ofs("trsol.txt");
	ofs << sol.size() << std::endl;
	for (size_t i = 0; i < sol.size(); i++)
		ofs << sol(i) << std::endl;

	//std::cout << "Success: " << std::boolalpha << success << std::endl;

	std::cout << "PARDISOTEST END" << std::endl;

	return 0;
}
