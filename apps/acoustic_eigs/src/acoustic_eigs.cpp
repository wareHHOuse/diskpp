/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2024
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <cstddef>
#include <iostream>
#include <regex>
#include <set>
#include <string>
#include <filesystem>



#include "sol/sol.hpp"
#include "diskpp/common/eigen.hpp"
#include "mumps.hpp"
#include "diskpp/common/util.h"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/loaders/loader_gmsh.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/output/silo.hpp"

#include "diskpp/bases/bases.hpp"
#include "diskpp/bases/bases_new.hpp"
#include "diskpp/bases/bases_operations.hpp"
#include "diskpp/methods/dg"
#include "diskpp/solvers/feast.hpp"

#include "diskpp/methods/hho_assemblers.hpp"
#include "diskpp/methods/hho_slapl.hpp"
#include "diskpp/solvers/solver.hpp"

namespace priv {

/* it exists already in MC/hho+dd */
template<typename T>
struct toltype {
    using type = T;
};

template<typename T>
struct toltype<std::complex<T>> {
    using type = T;
}; 

template<typename T>
std::optional<std::pair<Eigen::Matrix<T, Eigen::Dynamic, 1>, T>>
inv_powiter(const Eigen::SparseMatrix<T>& A,
    const Eigen::SparseMatrix<T>& B,
    typename toltype<T>::type mu,
    typename toltype<T>::type tol = 1e-10,
    size_t maxiter = 1000)
{
    assert(A.rows() == A.cols());
    assert(B.rows() == B.cols());
    assert(A.rows() == B.rows());

    Eigen::SparseMatrix<T> M = A - mu*B;
    mumps_solver<T> s;
    s.factorize(M);
    using vec_t = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    vec_t x = vec_t::Ones(M.cols());
    x = x/x.norm();
    T lambda = mu;
    size_t i = 0;
    for (; i < maxiter; i++) {
        vec_t Bx = B*x;
        vec_t new_x = s.solve(Bx);
        T new_lambda = mu + 1.0/x.dot(new_x);
        x = new_x/new_x.norm();
        auto err = std::abs(lambda-new_lambda);
        if ( err <= tol*std::abs(new_lambda) ) {
            std::cout << std::endl;
            return std::pair{x, new_lambda};
        }
        std::cout << "\rInverse power iteration " << i << ": " << err/std::abs(new_lambda) << std::flush;
        lambda = new_lambda;
    }
    std::cout << std::endl;
    if (i == maxiter) {
        return {};
    }

    return std::pair{x, lambda};
}
}

namespace autogen {


using Vec = Eigen::VectorXd;
using Mat = Eigen::MatrixXd;

/*
 * B-orthonormalize a vector 's' against columns of 'V'
 * Uses Modified Gram-Schmidt with B-inner product
 * Ensures numerical stability and maintains B-orthogonality
 */
void B_orthonormalize(Vec& s, const Mat& V, const Mat& B) {
    for (int j = 0; j < V.cols(); ++j) {
        double proj = V.col(j).dot(B*s); // project onto existing basis vector
        s -= proj * V.col(j);
    }
    double norm = std::sqrt(s.dot(B*s));
    /*if (norm > 1e-12)*/ s /= norm; // normalize if not nearly zero
}

Vec cg(const std::function<Vec(const Vec&)>& J,
           const Vec& b,
           const std::function<Vec(const Vec&)>& precond,
           int max_iter = 50,
           double tol = 1e-6)
{
    double  nr, nr0;
    double  alpha, beta, rho;
    auto N = b.size();
    Vec d(N), r(N), r0(N), y(N), Pr(N), Pr1(N);
    Vec x = Vec::Zero(N);

    r0 = r = b - J(x);
    d = precond(r0);
    nr = nr0 = r.norm();

    for (size_t iter = 0; iter < max_iter; iter++)
    {
        double rr = nr/nr0;
        if (rr < tol) {
            std::cout << "  CG: converged with RR = " << rr << std::endl;
            return x;
        }

        y     = J(d);
        Pr    = precond(r);
        rho   = r.dot(Pr);
        alpha = rho / d.dot(y);
        x     = x + alpha * d;
        r     = r - alpha * y;
        Pr1   = precond(r);
        beta  = r.dot(Pr1) / rho;
        d     = Pr1 + beta * d;

        nr = r.norm();
    }

    std::cout << "  CG reached max_iter with RR = " << nr/nr0 << std::endl;
    return x;
}

#include <Eigen/Dense>
#include <iostream>
#include <cmath>

template <typename Function>
int minres(const Function A,
           const Vec& b,
           Vec& x,
           int maxIters = 1000,
           double tol = 1e-8)
{
    using Scalar = typename Vec::Scalar;
    using std::sqrt;

    int n = b.rows();
    Vec r = b - A(x);
    Vec v_old = Vec::Zero(n);
    Vec v = r;

    Scalar nr = r.norm();
    Scalar nr0 = nr;

    Scalar beta = r.norm();
    if (beta == 0.0) return 0;

    v /= beta;

    Scalar beta_old = 0;
    Scalar c_old = 1, s_old = 0;
    Scalar c = 1, s = 0;
    Scalar eta = beta;

    Vec w_old = Vec::Zero(n);
    Vec w = Vec::Zero(n);
    Vec w_new = Vec::Zero(n);

    for (int k = 0; k < maxIters; ++k)
    {
        Vec Av = A(v);
        Scalar alpha = v.dot(Av);
        Av -= alpha * v + beta_old * v_old;

        beta_old = Av.norm();
        v_old = v;

        if (beta_old > 0)
            v = Av / beta_old;
        else
            break;

        // Apply previous rotation
        Scalar delta = c * alpha - c_old * s * beta_old;
        Scalar gamma = sqrt(delta * delta + beta_old * beta_old);
        c = delta / gamma;
        s = beta_old / gamma;

        // Update direction
        w_new = (v_old - s_old * w_old - c_old * w) / gamma;

        x += c * eta * w_new;
        eta = -s * eta;

        auto rr = nr/nr0;

        if (rr < tol) {
            std::cout << "MINRES: converged with rr = " << rr << std::endl;
            return k + 1;
        }

        // Rotate variables
        w_old = w;
        w = w_new;
        c_old = c;
        s_old = s;

        nr = (b - A(x)).norm();
    }

    std::cout << "MINRES: reached max_iter with rr = " << nr/nr0 << std::endl;

    return maxIters;
}


/*
 * Robust block Jacobi-Davidson for generalized Hermitian eigenproblem:
 *      A x = lambda B x
 * Parameters:
 *  n             : dimension of problem
 *  k             : number of eigenpairs to compute
 *  applyA        : matrix-free function to compute A * x
 *  B             : explicit SPD matrix (mass matrix)
 *  solveB        : function to solve B * y = x (used for preconditioning, optional)
 *  eigenvalues   : output vector of computed eigenvalues (size k)
 *  eigenvectors  : output matrix of eigenvectors (size n x k)
 *  maxOuter      : maximum outer JD iterations
 *  subspaceLimit : subspace size factor for thick restart (maxSubspace = subspaceLimit * k)
 */
void block_jacobi_davidson(
    int n,
    int k,  
    std::function<Vec(const Vec&)> applyA,
    const Mat& B,
    std::function<Vec(const Vec&)> solveB,
    Vec& eigenvalues,
    Mat& eigenvectors,
    int maxOuter = 100,
    int subspaceLimit = 3)
{
    int blockSize = k; // initial block size = number of requested eigenpairs

    // --- 1. Initialize subspace V with random vectors ---
    Mat V = Mat::Random(n, blockSize);
    for (int i = 0; i < V.cols(); ++i) {
        Vec v = V.col(i);
        B_orthonormalize(v, V.leftCols(i), B); // B-orthonormalization
        V.col(i) = v;
    }

    // --- 2. Locking array for converged eigenvectors ---
    std::vector<bool> locked(k, false);

    // --- 3. Outer JD loop ---
    for (int outer = 0; outer < maxOuter; ++outer) {
        std::cout << "outer: " << outer << std::endl;
        int m = V.cols(); // current subspace size

        // --- 3a. Compute AV and BV matrices (apply A and B to subspace vectors) ---
        Mat AV(n, m), BV(n, m);
        for (int i = 0; i < m; ++i) {
            AV.col(i) = applyA(V.col(i));
            BV.col(i) = B * V.col(i);
        }
        // --- 3b. Projected Rayleigh-Ritz problem ---
        Mat TA = V.transpose() * AV; // projected A
        Mat TB = V.transpose() * BV; // projected B
        
        Eigen::GeneralizedSelfAdjointEigenSolver<Mat> ges(TA, TB);
        Vec theta = ges.eigenvalues();        // Ritz values
        Mat Y = ges.eigenvectors();           // Ritz vectors in subspace

        // --- 3c. Compute Ritz vectors in original space ---
        Mat X = V * Y.leftCols(k);

        bool allConverged = true;

        // --- 4. Loop over each Ritz pair ---
        for (int i = 0; i < k; ++i) {
            if (locked[i])
                continue;

            Vec xi = X.col(i);         // current Ritz vector
            double lambda = theta(i);  // corresponding Ritz value

            // --- 4a. Compute residual r = A x - lambda B x ---
            Vec r = applyA(xi) - lambda * (B * xi);

            // --- 4b. Check convergence ---
            std::cout << "  norm " << i << ": " << r.norm() << std::endl;
            if (r.norm() < 1e-5 * std::abs(lambda)) {
                locked[i] = true;
                continue;
            }
            allConverged = false;

            // --- 4c. Define Jacobi-Davidson projected operator ---
            auto JD_operator = [&](const Vec& s) -> Vec {
                Vec Bs = B*s;
                Vec v1 = s - xi * xi.dot(Bs); // (I - u*u'*B)*s
                Vec v2 = applyA(v1) - lambda*B*v1; //(A - lambda*B)*v1
                return v2 - B* (xi * xi.dot(v2)); //(I - B*x*x')*v2
            };

            auto id = [](const Vec& v) -> Vec { Vec w = v; return w; };

            // --- 4d. Solve correction equation approximately using MINRES ---
            Vec s = cg(JD_operator, -r, solveB, 10, 0.2);
            //Vec s = Vec::Zero(r.rows());
            //minres(JD_operator, -r, s, 10, 0.2);

            // --- 4e. B-orthonormalize new vector against subspace V ---
            B_orthonormalize(s, V, B);
            
            // --- 4f. Expand subspace if vector is nonzero ---
            if (s.norm() > 1e-12) {
                std::cout << "Expanding subspace. Norm: " << s.norm();
                std::cout << ", size: " << V.cols()+1 << std::endl;
                V.conservativeResize(n, V.cols() + 1);
                V.col(V.cols() - 1) = s;
            }
        }

        // --- 5. Stop if all eigenpairs converged ---
        if (allConverged)
            break;

        // --- 6. Thick restart if subspace grows too large ---
        if (V.cols() > subspaceLimit * k) {
            V = X.leftCols(k); // retain best k Ritz vectors
        }
    }

    // --- 7. Final Rayleigh-Ritz to compute eigenvalues/eigenvectors ---
    int m = V.cols();
    Mat AV(n, m), BV(n, m);
    for (int i = 0; i < m; ++i) {
        AV.col(i) = applyA(V.col(i));
        BV.col(i) = B * V.col(i);
    }
    Mat TA = V.transpose() * AV;
    Mat TB = V.transpose() * BV;
    Eigen::GeneralizedSelfAdjointEigenSolver<Mat> ges(TA, TB);

    eigenvalues.resize(k);
    eigenvectors = V * ges.eigenvectors().leftCols(k);
    for (int i = 0; i < k; ++i)
        eigenvalues(i) = ges.eigenvalues()(i);
}

}


namespace disk {

template<typename Mesh>
auto
acoustic_eigs_dg(Mesh& msh, size_t degree,
    const typename Mesh::coordinate_type eta,
    disk::silo_database& silo)
{   
    std::cout << "DG eigsolver" << std::endl;
    auto cvf = connectivity_via_faces(msh);
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    auto basis_rescaling = disk::basis::rescaling_strategy::inertial;

    auto cbs = disk::scalar_basis_size(degree, Mesh::dimension);
    auto assm = make_discontinuous_galerkin_eigenvalue_assembler(msh, cbs);
    
    timecounter tc;
    tc.tic();
    for (auto& tcl : msh)
    {
        auto tbasis = disk::basis::scaled_monomial_basis(msh, tcl, degree, basis_rescaling);
        
        matrix_type M = integrate(msh, tcl, tbasis, tbasis);
        matrix_type K = integrate(msh, tcl, grad(tbasis), grad(tbasis));

        assm.assemble(msh, tcl, tcl, K+M);
        assm.assemble(msh, tcl, M);

        auto fcs = faces(msh, tcl);
        for (auto& fc : fcs)
        {   
            auto n     = normal(msh, tcl, fc);
            auto eta_l = eta / diameter(msh, fc);
            
            auto nv = cvf.neighbour_via(msh, tcl, fc);
            if (nv) {
                matrix_type Att = matrix_type::Zero(tbasis.size(), tbasis.size());
                matrix_type Atn = matrix_type::Zero(tbasis.size(), tbasis.size());
                
                auto ncl = nv.value();
                auto nbasis = disk::basis::scaled_monomial_basis(msh, ncl, degree, basis_rescaling);
                assert(tbasis.size() == nbasis.size());

                Att += + eta_l * integrate(msh, fc, tbasis, tbasis);
                Att += - 0.5 * integrate(msh, fc, grad(tbasis).dot(n), tbasis);
                Att += - 0.5 * integrate(msh, fc, tbasis, grad(tbasis).dot(n));

                Atn += - eta_l * integrate(msh, fc, nbasis, tbasis);
                Atn += - 0.5 * integrate(msh, fc, grad(nbasis).dot(n), tbasis);
                Atn += + 0.5 * integrate(msh, fc, nbasis, grad(tbasis).dot(n));

                assm.assemble(msh, tcl, tcl, Att);
                assm.assemble(msh, tcl, ncl, Atn);
            }
            else {
                //matrix_type Att = matrix_type::Zero(tbasis.size(), tbasis.size());
                //Att += + eta_l * integrate(msh, fc, tbasis, tbasis);
                //Att += - integrate(msh, fc, grad(tbasis).dot(n), tbasis);
                //Att += - integrate(msh, fc, tbasis, grad(tbasis).dot(n));
                //assm.assemble(msh, tcl, tcl, Att);
            }   
        }
    }

    assm.finalize();
    std::cout << " Assembly time: " << tc.toc() << std::endl;

    std::cout << " Unknowns: " << assm.gK.rows() << " ";
    std::cout << " Nonzeros: " << assm.gK.nonZeros() << std::endl;

    disk::feast_eigensolver_params<T> fep;
    fep.subspace_size = 50;
    fep.min_eigval = 5;
    fep.max_eigval = 100;
    fep.verbose = true;
    fep.max_iter = 50;
    fep.tolerance = 8;
    fep.fis = feast_inner_solver::mumps;

    T pisq = M_PI * M_PI;

    Eigen::Matrix<T, Eigen::Dynamic, 1> eigvals_ref =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(12);
    eigvals_ref <<
        pisq+1, pisq+1, 2*pisq+1, 4*pisq+1, 4*pisq+1, 5*pisq+1,
        5*pisq+1, 8*pisq+1, 9*pisq+1, 9*pisq+1, 10*pisq+1, 10*pisq+1
    ;


    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigvecs;
    Eigen::Matrix<T, Eigen::Dynamic, 1> eigvals;

    std::cout << "Running FEAST" << std::endl;

    auto fs = disk::feast(fep, assm.gK, assm.gM, eigvecs, eigvals);

    std::cout << (eigvals - eigvals_ref).transpose() << std::endl;

    silo.add_mesh(msh, "mesh");
    for (size_t col = 0; col < eigvecs.cols(); col++) {
        Eigen::Matrix<T, Eigen::Dynamic, 1> eigvec = eigvecs.col(col);
        std::cout << eigvals(col) << std::endl;
        std::vector<T> u;
        for (auto& cl : msh) {
            auto ofs = cbs * offset(msh, cl);
            u.push_back(eigvecs(ofs, col));
        }
        
        std::string vname = "eigfun_" + std::to_string(col);
        silo.add_variable("mesh", vname, u, disk::zonal_variable_t);
    }

    auto opt_evs = ::priv::inv_powiter(assm.gK, assm.gM, 7.0);
    if (opt_evs) {
        auto [vec, val] = *opt_evs;
        std::cout << "with powiter: " << val << std::endl;
    }
}

template<typename Mesh>
void
acoustic_eigs_hho(const Mesh& msh, size_t degree, disk::silo_database& silo)
{
    std::cout << "HHO eigsolver" << std::endl;
    using namespace disk::basis;
    using namespace disk::hho::slapl;

    using mesh_type = Mesh;
    using T = typename hho_space<Mesh>::scalar_type;
    using cbasis_type = typename hho_space<Mesh>::cell_basis_type;
    using fbasis_type = typename hho_space<Mesh>::face_basis_type;

    degree_info di(degree);

    auto assm = disk::hho::eigenvalue_block_assembler<Mesh, cbasis_type, fbasis_type>(
        msh, di.cell, di.face
    );

    timecounter tc;
    tc.tic();
    for (auto& cl : msh)
    {
        auto [R, A] = local_operator(msh, cl, di);
        auto S = local_stabilization(msh, cl, di, R);
        disk::dynamic_matrix<T> lhs = A+S;

        auto phiT = hho_space<mesh_type>::cell_basis(msh, cl, di.cell);
        disk::dynamic_matrix<T> rhs = integrate(msh, cl, phiT, phiT);

        auto cbs = phiT.size();
        lhs.block(0,0,cbs,cbs) += rhs;

        assm.assemble(msh, cl, lhs, rhs);
    }
    assm.finalize();

    std::cout << " Assembly time: " << tc.toc() << std::endl;
    //std::cout << " Unknowns: " << assm.LHS.rows() << " ";
    //std::cout << " Nonzeros: " << assm.LHS.nonZeros() << std::endl;
    tc.tic();

    Eigen::PardisoLDLT< Eigen::SparseMatrix<T> > AFF_lu(assm.AFF);
    Eigen::PardisoLDLT< Eigen::SparseMatrix<T> > BTT_lu(assm.BTT);

    //Eigen::SparseMatrix<T> KTT = assm.ATT - assm.ATF*AFF_lu.solve(assm.AFT);

    //std::cout << KTT.rows() << " " << KTT.cols() << std::endl;
    //std::cout << " Nonzeros: " << KTT.nonZeros() << std::endl;

    /*
    disk::feast_eigensolver_params<T> fep;
    fep.subspace_size = 50;
    fep.min_eigval = 5;
    fep.max_eigval = 100;
    fep.verbose = true;
    fep.max_iter = 50;
    fep.tolerance = 8;
    fep.fis = feast_inner_solver::mumps;

    */

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigvecs;
    Eigen::Matrix<T, Eigen::Dynamic, 1> eigvals;


    auto apply_KTT = [&](const dynamic_vector<T>& v) -> dynamic_vector<T> {
        dynamic_vector<T> z = assm.AFT*v;
        dynamic_vector<T> w = AFF_lu.solve(z);
        return assm.ATT*v - assm.ATF*w;
    };

    auto solve_BTT = [&](const dynamic_vector<T>& v) -> dynamic_vector<T> {
        return BTT_lu.solve(v);
    };

    std::cout << "starting BDJ" << std::endl;
    autogen::block_jacobi_davidson(assm.BTT.rows(), 10, apply_KTT,
        assm.BTT, solve_BTT, eigvals, eigvecs);

    T pisq = M_PI * M_PI;

    Eigen::Matrix<T, Eigen::Dynamic, 1> eigvals_ref =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(12);
    eigvals_ref <<
        pisq+1, pisq+1, 2*pisq+1, 4*pisq+1, 4*pisq+1, 5*pisq+1,
        5*pisq+1, 8*pisq+1, 9*pisq+1, 9*pisq+1, 10*pisq+1, 10*pisq+1
    ;

    //std::cout << "Running FEAST" << std::endl;

    //auto fs = disk::feast(fep, KTT, assm.BTT, eigvecs, eigvals);

    //std::cout << (eigvals - eigvals_ref).transpose() << std::endl;

    silo.add_mesh(msh, "hmesh");
    for (size_t col = 0; col < eigvecs.cols(); col++) {
        Eigen::Matrix<T, Eigen::Dynamic, 1> eigvec = eigvecs.col(col);
        std::cout << eigvals(col) << std::endl;
        std::vector<T> u;
        for (auto& cl : msh) {
            auto ofs = cbasis_type::size_of_degree(di.cell) * offset(msh, cl);
            u.push_back(eigvecs(ofs, col));
        }
        
        std::string vname = "hho_eigfun_" + std::to_string(col);
        silo.add_variable("hmesh", vname, u, disk::zonal_variable_t);
    }

    //auto opt_evs = ::priv::inv_powiter(KTT, assm.BTT, 57.0);
    //if (opt_evs) {
    //    auto [vec, val] = *opt_evs;
    //    std::cout << "with powiter: " << val << std::endl;
    //}
}

}

int main(int argc, char **argv)
{
    std::string mesh_filename = argv[1];
    using T = double;
    using mesh_type = disk::simplicial_mesh<T,2>;
    mesh_type msh;
    disk::gmsh_geometry_loader< mesh_type > loader;
    loader.read_mesh(mesh_filename);
    loader.populate_mesh(msh);

    disk::silo_database db;
    db.create("acoustic_eigs.silo");
    //acoustic_eigs_dg(msh, 2, 10, db);

    acoustic_eigs_hho(msh, 1, db);


    return 0;
}