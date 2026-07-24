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
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <regex>
#include <set>
#include <string>
#include <unistd.h>

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
#include "diskpp/solvers/direct_solvers.hpp"

#include "diskpp/solvers/eigensolvers.hpp"



enum class eigsolver_type {
    feast_full,
    feast_mf,
    bjd_mf
};

enum class mesh_source {
    internal_tri,
    internal_quad,
    internal_hex,
    internal_tet,
    internal_brick,
    external
};

struct config {
    eigsolver_type  eigsolver = eigsolver_type::feast_full;
    size_t          order = 0;
    size_t          reflevels = 5;
    mesh_source     source = mesh_source::external;
    std::string     mesh_filename;
    std::string     silo_filename;
};

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
        assm.gM, eigvecs, eigvals);

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

template<typename T>
void solve_feast_dense(auto assm, disk::dynamic_matrix<T>& eigvecs,
    disk::dynamic_vector<T>& eigvals)
{
    timecounter tc;

    std::cout << "MUMPS factorization..." << std::flush;
    tc.tic();
    disk::solvers::mumps_solver<T> AFF_lu_mumps;
    AFF_lu_mumps.factorize(assm.AFF);
    std::cout << tc.toc() << " seconds\n";

    std::cout << "FEAST eigensolver (dense)" << std::endl;
    tc.tic();
    disk::solvers::feast_eigensolver_params<T> params;
    params.subspace_size = 20;
    params.max_iter = 10;
    params.min_eigval = 1;
    params.max_eigval = 50;
    params.verbose = true;
    params.tolerance = 7;

    std::cout << "Computing KTT\n";
    disk::dynamic_matrix<T> KTT =
        assm.ATT - assm.ATF*AFF_lu_mumps.solve(assm.AFT);
    std::cout << "Entering FEAST\n";
    disk::solvers::feast(params, KTT, assm.BTT, eigvecs, eigvals);

    std::cout << "Eigensolver time: " << tc.toc() << " seconds\n";
}

template<typename T>
void solve_feast_mf(auto assm, disk::dynamic_matrix<T>& eigvecs,
    disk::dynamic_vector<T>& eigvals)
{
    timecounter tc;

    Eigen::PardisoLDLT< Eigen::SparseMatrix<T> > AFF_lu(assm.AFF);

    auto apply_A = [&]<int ncols>(
        const Eigen::Matrix<T, Eigen::Dynamic, ncols>& v) ->
            Eigen::Matrix<T, Eigen::Dynamic, ncols> {
        Eigen::Matrix<T, Eigen::Dynamic, ncols> z = assm.AFT*v;
        return assm.ATT*v - assm.ATF*AFF_lu.solve(z);
    };

    std::cout << "FEAST eigensolver (matrix-free)" << std::endl;
    tc.tic();
    disk::solvers::feast_eigensolver_params<T> params;
    params.subspace_size = 20;
    params.max_iter = 10;
    params.min_eigval = 1;
    params.max_eigval = 50;
    params.verbose = true;
    params.tolerance = 7;
    disk::solvers::feast_mf(params, apply_A, assm.BTT, eigvecs, eigvals);

    std::cout << "Eigensolver time: " << tc.toc() << " seconds\n";
}


template<typename T>
void solve_bjd_mf(auto assm, disk::dynamic_matrix<T>& eigvecs,
    disk::dynamic_vector<T>& eigvals)
{
    timecounter tc;

    //Eigen::PardisoLDLT< Eigen::SparseMatrix<T> > AFF_lu(assm.AFF);
    disk::solvers::mumps_solver<T> AFF_lu;
    AFF_lu.symmetric(true);
    AFF_lu.factorize(assm.AFF);

    auto apply_A = [&]<int ncols>(
        const Eigen::Matrix<T, Eigen::Dynamic, ncols>& v) ->
            Eigen::Matrix<T, Eigen::Dynamic, ncols> {
        Eigen::Matrix<T, Eigen::Dynamic, ncols> z = assm.AFT*v;
        return assm.ATT*v - assm.ATF*AFF_lu.solve(z);
    };

    std::cout << "Block Jacobi-Davidson eigensolver" << std::endl;
    tc.tic();
    disk::solvers::bjd_params params;
    params.block_size = 10;
    params.max_outer_iters = 200;
    params.max_inner_iters = 5;
    params.max_subspace_growth = 10;
    params.inner_tol = 1e-4;
    params.verbose = true;
    disk::solvers::block_jacobi_davidson(params, apply_A,
        assm.BTT, eigvecs, eigvals);
    std::cout << "Eigensolver time: " << tc.toc() << " seconds\n";
}

#if 0
template<typename T>
void solve_spectra(auto assm, disk::dynamic_matrix<T>& eigvecs,
    disk::dynamic_vector<T>& eigvals)
{
    timecounter tc;

    Eigen::PardisoLDLT< Eigen::SparseMatrix<T> > AFF_lu(assm.AFF);
    Eigen::PardisoLDLT< Eigen::SparseMatrix<T> > BTT_lu(assm.BTT);

    auto apply_A = [&]<int ncols>(
        const Eigen::Matrix<T, Eigen::Dynamic, ncols>& v) ->
            Eigen::Matrix<T, Eigen::Dynamic, ncols> {
        Eigen::Matrix<T, Eigen::Dynamic, ncols> z = assm.AFT*v;
        return assm.ATT*v - assm.ATF*AFF_lu.solve(z);
    };

    auto solve_B = [&](const Eigen::Matrix<T, Eigen::Dynamic, ncols>& v) ->
            Eigen::Matrix<T, Eigen::Dynamic, ncols> {
        Eigen::Matrix<T, Eigen::Dynamic, ncols> z = assm.AFT*v;
        return assm.ATT*v - assm.ATF*AFF_lu.solve(z);
    };

    std::cout << "Block Jacobi-Davidson eigensolver" << std::endl;
    tc.tic();
    disk::solvers::bjd_params params;
    params.block_size = 10;
    params.max_outer_iters = 200;
    params.max_inner_iters = 5;
    params.max_subspace_growth = 10;
    params.inner_tol = 1e-4;
    params.verbose = true;
    disk::solvers::block_jacobi_davidson(params, apply_A,
        assm.BTT, eigvecs, eigvals);
    std::cout << "Eigensolver time: " << tc.toc() << " seconds\n";
}
#endif


template<typename Mesh>
void
acoustic_eigs_hho(const Mesh& msh, const config& cfg, disk::silo_database& silo)
{
    std::cout << "HHO eigsolver" << std::endl;
    using namespace disk::basis;
    using namespace disk::hho::slapl;

    using mesh_type = Mesh;
    using T = typename hho_space<Mesh>::scalar_type;
    using cbasis_type = typename hho_space<Mesh>::cell_basis_type;
    using fbasis_type = typename hho_space<Mesh>::face_basis_type;

    degree_info di(cfg.order);

    auto assm = disk::hho::eigenvalue_block_assembler<Mesh, cbasis_type, fbasis_type>(
        msh, di.cell, di.face
    );

    timecounter tc;

    /********* ASSEMBLY *********/
    std::cout << "Assembling matrices..." << std::flush;
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

    std::cout << "AFF: " << assm.AFF.rows() << " dofs, BTT: ";
    std::cout << assm.BTT.rows() << std::endl;

    std::cout << tc.toc() << " seconds\n";


    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigvecs;
    Eigen::Matrix<T, Eigen::Dynamic, 1> eigvals;

    if (cfg.eigsolver == eigsolver_type::feast_full) {
        solve_feast_dense(assm, eigvecs, eigvals);
    }

    if (cfg.eigsolver == eigsolver_type::feast_mf) {
        solve_feast_mf(assm, eigvecs, eigvals);
    }

    if (cfg.eigsolver == eigsolver_type::bjd_mf) {
        solve_bjd_mf(assm, eigvecs, eigvals);
    }

    for (auto& ev : eigvals) {
        ev -= 1.0;
    }

    std::cout << "Computed: " << std::setprecision(15) << eigvals.transpose() << std::endl;

    silo.add_mesh(msh, "hmesh");
    for (size_t col = 0; col < eigvecs.cols(); col++) {
        Eigen::Matrix<T, Eigen::Dynamic, 1> eigvec = eigvecs.col(col);
        //std::cout << eigvals(col) << std::endl;
        std::vector<T> u;
        for (auto& cl : msh) {
            auto ofs = cbasis_type::size_of_degree(di.cell) * offset(msh, cl);
            u.push_back(eigvecs(ofs, col));
        }

        std::string vname = "hho_eigfun_" + std::to_string(col);
        silo.add_variable("hmesh", vname, u, disk::zonal_variable_t);
    }
}

}

template<typename Mesh>
void
run_eigsolver(const Mesh& msh, const config& cfg)
{
    disk::silo_database db;
    db.create(cfg.silo_filename);

    acoustic_eigs_hho(msh, cfg, db);
}

int main(int argc, char **argv)
{
    resmon rm("main");

    using T = double;

    config cfg;

    int opt;
    while ((opt = getopt(argc, argv, "e:f:k:m:o:r:")) != -1) {
        switch (opt) {

        case 'e':
            if (std::string(optarg) == "feast_full") {
                cfg.eigsolver = eigsolver_type::feast_full;
            } else if (std::string(optarg) == "feast_mf") {
                cfg.eigsolver = eigsolver_type::feast_mf;
            } else if (std::string(optarg) == "bjd_mf") {
                cfg.eigsolver = eigsolver_type::bjd_mf;
            } else {
                std::cout << "invalid solver type" << std::endl;
            }
            break;

        case 'f':
            cfg.mesh_filename = optarg;
            break;

        case 'k':
            cfg.order = std::stoul(optarg);
            break;

        case 'm':
            if (std::string(optarg) == "tri") {
                cfg.source = mesh_source::internal_tri;
            } else if (std::string(optarg) == "quad") {
                cfg.source = mesh_source::internal_quad;
            } else if (std::string(optarg) == "hex") {
                cfg.source = mesh_source::internal_hex;
            } else if (std::string(optarg) == "tet") {
                cfg.source = mesh_source::internal_tet;
            } else if (std::string(optarg) == "brick") {
                cfg.source = mesh_source::internal_brick;
            } else {
                std::cout << "invalid mesh type" << std::endl;
            }
            break;

        case 'o':
            cfg.silo_filename = optarg;
            break;

        case 'r':
            cfg.reflevels = std::stoul(optarg);
            break;
        }
    }


    if (cfg.mesh_filename != "") {

        if (std::regex_match(cfg.mesh_filename, std::regex(".*\\.geo2s$") ))
        {
            std::cout << "Guessed mesh format: GMSH 2D simplicials" << std::endl;
            using mesh_type = disk::triangular_mesh<T>;
            mesh_type msh;
            disk::gmsh_geometry_loader< mesh_type > loader;
            loader.read_mesh(cfg.mesh_filename);
            loader.populate_mesh(msh);

            run_eigsolver(msh, cfg);
            return 0;
        }

        if (std::regex_match(cfg.mesh_filename, std::regex(".*\\.geo3s$") ))
        {
            std::cout << "Guessed mesh format: GMSH 3D simplicials" << std::endl;
            using mesh_type = disk::tetrahedral_mesh<T>;
            mesh_type msh;
            disk::gmsh_geometry_loader< mesh_type > loader;
            loader.read_mesh(cfg.mesh_filename);
            loader.populate_mesh(msh);

            run_eigsolver(msh, cfg);
            return 0;
        }
    }

    if (cfg.source == mesh_source::internal_tri) {
        using mesh_type = disk::simplicial_mesh<T, 2>;
        mesh_type msh;
        auto mesher = disk::make_simple_mesher(msh);
        mesher.refine();

        msh.transform( [&](const typename mesh_type::point_type& pt) {
            return typename mesh_type::point_type{pt.x(), 1.1*pt.y()};
        } );

        for (int i = 0; i < cfg.reflevels; i++) {
            mesher.refine();

            std::cout << ">>>>>>>> DIAM: " << disk::average_diameter(msh) << std::endl;
            run_eigsolver(msh, cfg);
        }
    }

    if (cfg.source == mesh_source::internal_tet) {
        using mesh_type = disk::simplicial_mesh<T, 3>;
        mesh_type msh;
        auto mesher = disk::make_simple_mesher(msh);
        mesher.refine();

        msh.transform( [&](const typename mesh_type::point_type& pt) {
            return typename mesh_type::point_type{pt.x(), 1.1*pt.y(), 1.2*pt.z()};
        } );

        for (int i = 0; i < cfg.reflevels; i++) {
            mesher.refine();

            std::cout << ">>>>>>>> DIAM: " << disk::average_diameter(msh) << std::endl;
            run_eigsolver(msh, cfg);
        }
    }

    if (cfg.source == mesh_source::internal_quad) {
        using mesh_type = disk::cartesian_mesh<T, 2>;
        mesh_type msh;
        auto mesher = disk::make_simple_mesher(msh);
        mesher.refine();

        msh.transform( [&](const typename mesh_type::point_type& pt) {
            return typename mesh_type::point_type{pt.x(), 1.1*pt.y()};
        } );

        for (int i = 0; i < cfg.reflevels; i++) {
            mesher.refine();

            std::cout << ">>>>>>>> DIAM: " << disk::average_diameter(msh) << std::endl;
            run_eigsolver(msh, cfg);
        }
    }

    /*
    if (cfg.source == mesh_source::internal_brick) {
        using mesh_type = disk::cartesian_mesh<T, 3>;
        mesh_type msh;
        auto mesher = disk::make_simple_mesher(msh);
        mesher.refine();

        msh.transform( [&](const typename mesh_type::point_type& pt) {
            return typename mesh_type::point_type{pt.x(), 1.1*pt.y(), 1.2*pt.z()};
        } );

        for (int i = 0; i < cfg.reflevels; i++) {
            mesher.refine();

            std::cout << ">>>>>>>> DIAM: " << disk::average_diameter(msh) << std::endl;
            run_eigsolver(msh, cfg);
        }
    }
    */

    if (cfg.source == mesh_source::internal_hex) {
        using mesh_type = disk::generic_mesh<T, 2>;
        mesh_type msh;
        auto mesher = disk::make_fvca5_hex_mesher(msh);

        for (int i = 0; i < cfg.reflevels; i++) {
            mesher.make_level(i);


            msh.transform( [&](const typename mesh_type::point_type& pt) {
                return typename mesh_type::point_type{pt.x(), 1.1*pt.y()};
            } );

            std::cout << ">>>>>>>> DIAM: " << disk::average_diameter(msh) << std::endl;
            run_eigsolver(msh, cfg);
        }
    }

    return 0;
}
