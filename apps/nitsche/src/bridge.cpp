/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2025
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

 #include <iostream>
 #include <algorithm>
 
 #include "diskpp/mesh/mesh.hpp"
 #include "diskpp/mesh/meshgen.hpp"
 #include "diskpp/loaders/loader.hpp"
 #include "diskpp/bases/bases.hpp"
 #include "diskpp/methods/hho"
 #include "diskpp/methods/implementation_hho/curl.hpp"
 #include "diskpp/methods/hho_slapl.hpp"
 #include "diskpp/methods/hho_assemblers.hpp"
 #include "mumps.hpp"
 #include "diskpp/common/timecounter.hpp"
 #include "diskpp/output/silo.hpp"
 #include "operators.hpp"
 #include "diskpp/solvers/solver.hpp"
 #include "asm.hpp"


template<typename Mesh>
void
set_dirichlet(const Mesh& msh, std::vector<bool>& bcs, size_t bnd)
{
    if (bcs.size() != msh.faces_size()) {
        bcs.resize( msh.faces_size() );
    }

    std::cout << "set dirichlet" << std::endl;

    size_t fcnum = 0;
    for (auto& fc : faces(msh)) {
        auto bi = msh.boundary_info(fc);
        if (bi.is_boundary() and bi.id() == bnd) {
            std::cout << "fcnum: " << fcnum << ", id: " << bi.id() << ", tag: " << bi.tag() << std::endl;
            bcs[fcnum] = true;
        }
        fcnum++;
    }   
}

int main(int argc, char **argv)
{
    using T = double;
    size_t degree = 1;
    hho_mode mode = hho_mode::nitsche;
    T mu = 10.0;
    T lambda = 1000.0;

    if (argc < 2) {
        std::cout << "args!" << std::endl;
        return 1;
    }

    std::string mesh_filename = argv[1];

    using mesh_type = disk::simplicial_mesh<T,3>;

    mesh_type msh;
        
    disk::gmsh_geometry_loader< mesh_type > loader;    
    loader.read_mesh(mesh_filename);
    loader.populate_mesh(msh);

    std::cout << "Cells: " << msh.cells_size() << std::endl;
    std::cout << "Faces: " << msh.faces_size() << std::endl;
    std::cout << "Boundary faces: " << msh.boundary_faces_size() << std::endl;

    const static size_t DIM = mesh_type::dimension;
    auto fbs = disk::vector_basis_size(degree, DIM-1, DIM);

    std::vector<bc> bcs(msh.faces_size());
    set_all_neumann(msh, bcs);
    set_boundary(msh, bcs, bc::dirichlet, 967);
    set_boundary(msh, bcs, bc::dirichlet, 994);

    std::vector<bool> dirfaces;
    dirfaces.resize( bcs.size() );

    auto df = [&](bc b) {
        if (mode == hho_mode::standard) {
            return b == bc::dirichlet;
        }
        
        return (b == bc::dirichlet) or (b == bc::neumann);
    };

    std::transform(bcs.begin(), bcs.end(), dirfaces.begin(), df);

    condensed_assembler assm(msh, fbs, dirfaces);

    using MT = disk::dynamic_matrix<T>;
    using VT = disk::dynamic_vector<T>;
    std::vector<std::pair<MT, VT>> lcs;

    timecounter tc;
    tc.tic();
    std::cout << "ASM: " << std::flush;
    for (auto& cl : msh) {
        auto [SGR, Asgr] = hho_mixedhigh_symlapl(msh, cl, degree, mode, bcs);
        auto [DR, Adr] = hho_mixedhigh_divrec(msh, cl, degree, mode, bcs);
        disk::hho_degree_info hdi(degree+1, degree);
        //auto S = disk::make_vector_hho_stabilization(msh, cl, SGR, hdi);
        auto S = vstab(msh, cl, degree, mode, bcs);

        MT lhs = 2*mu*Asgr + lambda*Adr + 2*mu*S;
        auto cb = disk::make_vector_monomial_basis(msh, cl, degree+1);
        auto cbs = cb.size();
        VT rhs = VT::Zero(lhs.rows());

        auto qps = disk::integrate(msh, cl, 2*degree+2);
        disk::static_vector<T,DIM> f =
            disk::static_vector<T,DIM>::Zero();
        f(0) = 0;
        f(1) = 0;
        f(2) = -0.001;

        for (auto& qp : qps) {
            auto phi = cb.eval_functions(qp.point());
            rhs.head(cbs) += qp.weight() * (phi * f);
        }

        lcs.push_back({lhs, rhs});

        auto [Lc, Rc] = disk::static_condensation(lhs, rhs, cbs);

        assm.assemble(msh, cl, Lc, Rc);
    }
    std::cout << tc.toc() << std::endl;

    assm.finalize();

    std::cout << "DoFs: " << assm.LHS.rows() << std::endl;
    std::cout << "NNZ: " << assm.LHS.nonZeros() << std::endl;

    tc.tic();
    std::cout << "MUMPS: " << std::flush;
    disk::dynamic_vector<T> sol = mumps_lu(assm.LHS, assm.RHS);
    //Eigen::SparseLU<Eigen::SparseMatrix<T>> solver(assm.LHS);
    //disk::dynamic_vector<T> sol = solver.solve(assm.RHS);
    /*
    disk::dynamic_vector<T> sol = assm.RHS;
    disk::solvers::conjugated_gradient_params<T> cgp;
    cgp.verbose = true;
    cgp.max_iter = assm.LHS.rows();
    disk::solvers::conjugated_gradient(cgp, assm.LHS, assm.RHS, sol);
    */

    Eigen::Matrix<T, Eigen::Dynamic, DIM> u_data =
        Eigen::Matrix<T, Eigen::Dynamic, DIM>::Zero(msh.cells_size(), DIM);

    std::cout << tc.toc() << std::endl;

    size_t cell_i = 0;
    for (auto& cl : msh)
    {
        const auto& [lhs, rhs] = lcs[cell_i];
        auto locsolF = assm.take_local_solution(msh, cl, sol);
        disk::dynamic_vector<T> locsol =
            disk::static_decondensation(lhs, rhs, locsolF);
        
        u_data(cell_i, 0) = locsol(0);
        u_data(cell_i, 1) = locsol(1);
        if constexpr (DIM == 3) {
            u_data(cell_i, 2) = locsol(2);
        }

        cell_i++;
    }

    disk::silo_database silo;
    silo.create("bridge.silo");
    silo.add_mesh(msh, "mesh");
    silo.add_variable("mesh", "u", u_data, disk::zonal_variable_t);

    return 0;
} 