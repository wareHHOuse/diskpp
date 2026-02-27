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

struct bjd_params {
    int     max_outer_iters         = 100;
    int     max_inner_iters         = 30;
    double  outer_tol               = 1e-8;
    double  inner_tol               = 1e-3;
    double  ev_tol                  = 1e-8;
    double  expand_thresh           = 1e-12;
    int     max_subspace_growth     = 5;
    int     block_size              = 10;
    bool    verbose                 = false;
};

template<typename T>
void orthonormalize(disk::dynamic_matrix<T>& V, const disk::sparse_matrix<T>& B)
{
    // Modified G-S B-orthonormalization
    for (int i = 0; i < V.cols(); i++) {
        double Bnorm = std::sqrt(V.col(i).dot(B*V.col(i)));
        V.col(i) /= Bnorm;
        for (int j = i+1; j < V.cols(); j++) {
            double proj = V.col(i).dot(B*V.col(j));
            V.col(j) = V.col(j) - proj * V.col(i);
        }
    }
}

template<typename T>
void block_jacobi_davidson(const bjd_params& params,
    auto applyA, const disk::sparse_matrix<T>& B,
    auto solveB, disk::dynamic_vector<T>& eigenvalues,
    disk::dynamic_matrix<T>& eigenvectors)
{
    const int N = B.rows();
    
    using dv = disk::dynamic_vector<T>;
    using dm = disk::dynamic_matrix<T>;

    dm V = dm::Random(N, params.block_size);
    orthonormalize(V, B);

    // Keeps track of converged eigenpairs
    std::vector<bool> conv(params.block_size, false);

    for (int outer = 0; outer < params.max_outer_iters; outer++) {
        std::cout << "outer: " << outer << std::endl;
        
        const int M = V.cols();
        dm TA = V.transpose() * applyA(V);
        dm TB = V.transpose() * B * V;
        Eigen::GeneralizedSelfAdjointEigenSolver<dm> ges(TA, TB);
        dv theta = ges.eigenvalues();
        dm Y = ges.eigenvectors();
        dm X = V * Y.leftCols(params.block_size);

        for (int i_ev = 0; i_ev < params.block_size; i_ev++) {
            if (conv[i_ev]) {
                continue; /* This pair has already converged. */
            }

            dv xi = X.col(i_ev);
            T lambda = theta(i_ev);
            dv r = applyA(xi) - lambda * (B * xi);

            std::cout << "  norm " << i_ev << ": " << r.norm() << std::endl;
            if (r.norm() < params.ev_tol * std::abs(lambda)) {
                conv[i_ev] = true;
                continue;
            }

            auto JDproj = [&](const dv& s) -> dv {
                Vec Bs = B*s;
                Vec v1 = s - xi * xi.dot(Bs); // (I - u*u'*B)*s
                Vec v2 = applyA(v1) - lambda*B*v1; //(A - lambda*B)*v1
                return v2 - B*(xi * xi.dot(v2)); //(I - B*x*x')*v2
            };

            auto JDproj2 = [&](const dv& s) -> dv {
                Vec Bs = B * s;
                Vec Ps = s - xi * (xi.dot(Bs));
                Vec KPs = applyA(Ps) - lambda * (B * Ps);
                return KPs;
            };

            auto id = [](const dv& v) -> dv { return v; };

            dv s = cg(JDproj, -r, solveB, params.max_inner_iters,
                params.inner_tol);
            
            orthonormalize(V, B);
            
            if (s.norm() > params.expand_thresh) {
                std::cout << "Expanding subspace. Norm: " << s.norm();
                std::cout << ", size: " << V.cols()+1 << std::endl;
                V.conservativeResize(N, V.cols() + 1);
                V.col(V.cols() - 1) = s;
            }
        }

        // Check if all the eigenpairs have converged
        if ( std::all_of(conv.begin(), conv.end(), [](bool v) { return v; }) ) {
            std::cout << "All eigenpairs have converged" << std::endl;
            break;
        }

        // Restart if subspace grows too much
        if ( V.cols() > params.block_size * params.max_subspace_growth ) {
            V = X.leftCols(params.block_size);
        }
    }

    // Final Rayleigh-Ritz to compute eigenvalues/eigenvectors
    int M = V.cols();
    dm TA = V.transpose() * applyA(V);
    dm TB = V.transpose() * B * V;
    Eigen::GeneralizedSelfAdjointEigenSolver<dm> ges(TA, TB);

    eigenvalues.resize(params.block_size);
    eigenvectors = V * ges.eigenvectors().leftCols(params.block_size);
    for (int i = 0; i < params.block_size; i++)
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

    Eigen::SparseMatrix<T> KTT = assm.ATT - assm.ATF*AFF_lu.solve(assm.AFT);

    //std::cout << KTT.rows() << " " << KTT.cols() << std::endl;
    //std::cout << " Nonzeros: " << KTT.nonZeros() << std::endl;

    
    disk::feast_eigensolver_params<T> fep;
    fep.subspace_size = 50;
    fep.min_eigval = 0.1;
    fep.max_eigval = 100;
    fep.verbose = true;
    fep.max_iter = 50;
    fep.tolerance = 8;
    fep.fis = feast_inner_solver::mumps;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigvecs;
    Eigen::Matrix<T, Eigen::Dynamic, 1> eigvals;


    auto apply_KTT = [&](const dynamic_vector<T>& v) -> dynamic_vector<T> {
        dynamic_vector<T> z = assm.AFT*v;
        dynamic_vector<T> w = AFF_lu.solve(z);
        return assm.ATT*v - assm.ATF*w;
    };

    auto apply_A = [&]<int ncols>(
        const Eigen::Matrix<T, Eigen::Dynamic, ncols>& v) ->
            Eigen::Matrix<T, Eigen::Dynamic, ncols> {
        Eigen::Matrix<T, Eigen::Dynamic, ncols> z = assm.AFT*v;
        Eigen::Matrix<T, Eigen::Dynamic, ncols> w = AFF_lu.solve(z);
        return assm.ATT*v - assm.ATF*w;
    };

    auto solve_BTT = [&](const dynamic_vector<T>& v) -> dynamic_vector<T> {
        return BTT_lu.solve(v);
    };

    std::cout << "starting BDJ" << std::endl;
    autogen::bjd_params params;
    params.block_size = 13;
    autogen::block_jacobi_davidson(params, apply_A,
        assm.BTT, solve_BTT, eigvals, eigvecs);

    T pisq = M_PI * M_PI;

    Eigen::Matrix<T, Eigen::Dynamic, 1> eigvals_ref =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(12);
    eigvals_ref << 1,
        pisq+1, pisq+1, 2*pisq+1, 4*pisq+1, 4*pisq+1, 5*pisq+1,
        5*pisq+1, 8*pisq+1, 9*pisq+1, 9*pisq+1, 10*pisq+1, 10*pisq+1
    ;

    //std::cout << "Running FEAST" << std::endl;
    //auto fs = disk::feast(fep, KTT, assm.BTT, eigvecs, eigvals);

    std::cout << (eigvals - eigvals_ref).transpose() << std::endl;

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
    acoustic_eigs_dg(msh, 2, 10, db);

    acoustic_eigs_hho(msh, 1, db);


    return 0;
}