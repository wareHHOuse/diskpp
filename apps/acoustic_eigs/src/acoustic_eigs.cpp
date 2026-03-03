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
#include "diskpp/solvers/direct_solvers.hpp"
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

#include "diskpp/solvers/eigensolvers.hpp"

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

    disk::solvers::feast_eigensolver_params<T> fep;
    fep.subspace_size = 50;
    fep.min_eigval = 5;
    fep.max_eigval = 100;
    fep.verbose = true;
    fep.max_iter = 50;
    fep.tolerance = 8;
    fep.fis = disk::solvers::feast_inner_solver::mumps;

    T pisq = M_PI * M_PI;

    Eigen::Matrix<T, Eigen::Dynamic, 1> eigvals_ref =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(12);
    eigvals_ref <<
        pisq+1, pisq+1, 2*pisq+1, 4*pisq+1, 4*pisq+1, 5*pisq+1,
        5*pisq+1, 8*pisq+1, 9*pisq+1, 9*pisq+1, 10*pisq+1, 10*pisq+1
    ;


    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigvecs;
    Eigen::Matrix<T, Eigen::Dynamic, 1> eigvals;

    //std::cout << "Running FEAST" << std::endl;
    //auto fs = disk::feast(fep, assm.gK, assm.gM, eigvecs, eigvals);

    std::cout << "starting BDJ" << std::endl;
    disk::solvers::bjd_params params;
    params.block_size = 13;
    params.max_inner_iters = 30;
    params.inner_tol = 1e-4;
    auto opA = disk::solvers::operator_from_matrix(assm.gK);
    disk::solvers::block_jacobi_davidson(params, opA,
        assm.gM, eigvals, eigvecs);

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

    
    disk::solvers::feast_eigensolver_params<T> fep;
    fep.subspace_size = 50;
    fep.min_eigval = 0.1;
    fep.max_eigval = 100;
    fep.verbose = true;
    fep.max_iter = 50;
    fep.tolerance = 8;
    fep.fis = disk::solvers::feast_inner_solver::mumps;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigvecs;
    Eigen::Matrix<T, Eigen::Dynamic, 1> eigvals;


    //auto apply_KTT = [&](const dynamic_vector<T>& v) -> dynamic_vector<T> {
    //    dynamic_vector<T> z = assm.AFT*v;
    //    dynamic_vector<T> w = AFF_lu.solve(z);
    //    return assm.ATT*v - assm.ATF*w;
    //};

    auto apply_A = [&]<int ncols>(
        const Eigen::Matrix<T, Eigen::Dynamic, ncols>& v) ->
            Eigen::Matrix<T, Eigen::Dynamic, ncols> {
        Eigen::Matrix<T, Eigen::Dynamic, ncols> z = assm.AFT*v;
        return assm.ATT*v - assm.ATF*AFF_lu.solve(z);
    };

    auto solve_BTT = [&](const dynamic_vector<T>& v) -> dynamic_vector<T> {
        return BTT_lu.solve(v);
    };

    std::cout << "starting BDJ" << std::endl;
    disk::solvers::bjd_params params;
    params.block_size = 13;
    params.max_inner_iters = 30;
    params.inner_tol = 1e-4;
    disk::solvers::block_jacobi_davidson(params, apply_A,
        assm.BTT, eigvals, eigvecs);

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
    //acoustic_eigs_dg(msh, 2, 10, db);

    acoustic_eigs_hho(msh, 1, db);


    return 0;
}
