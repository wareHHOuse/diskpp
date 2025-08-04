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

    size_t fcnum = 0;
    for (auto& fc : faces(msh)) {
        auto bi = msh.boundary_info(fc);
        if (bi.is_boundary() and bi.id() == bnd) {
            bcs[fcnum] = true;
        }
        fcnum++;
    }   
}

template<typename Mesher>
void make_mesh_level(Mesher& mesher, size_t level)
{
    for (size_t i = 0; i < level; i++) {
        mesher.refine();
    }
}

template<typename T>
void make_mesh_level(disk::fvca5_hex_mesher<disk::generic_mesh<T,2>>& mesher, size_t level)
{
    mesher.make_level(level);
}

template<typename Mesh>
auto get_mesher(Mesh& msh)
{
    return make_simple_mesher(msh);
}

template<typename T>
auto get_mesher(disk::generic_mesh<T,2>& msh)
{
    return make_fvca5_hex_mesher(msh);
}


int main(int argc, char **argv) {

    using T = double;

    int opt;
    std::string input_geometry;
    size_t degree = 1;
    size_t reflevel = 1;
    hho_mode mode = hho_mode::standard;
    T mu = 1.0;
    T lambda = 1.0;
    T fx = 0.0;
    T fy = 0.0;
    bool do_cook = false;

    while ((opt = getopt(argc, argv, "k:m:l:r:g:ncx:y:")) != -1) {
        switch (opt) {
        case 'k': /* method order */
            degree = std::max(1, std::stoi(optarg));
            break;
        case 'r': /* builtin mesh refinement levels */
            reflevel = std::max(1, std::stoi(optarg));
            break;
        case 'n': /* enable nitsche */
            mode = hho_mode::nitsche;
            break;
        case 'm':
            mu = std::stod(optarg);
            break;
        case 'l':
            lambda = std::stod(optarg);
            break;
        case 'c':
            do_cook = true;
            break;
        case 'x':
            fx = std::stod(optarg);
            break;
        case 'y':
            fy = std::stod(optarg);
            break;
        default:
            std::cout << "Invalid options." << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    std::cout << "mu = " << mu << ", lambda = " << lambda << std::endl; 

    using mesh_type = disk::simplicial_mesh<T,2>;
    //using mesh_type = disk::cartesian_mesh<T,2>;
    //using mesh_type = disk::generic_mesh<T,2>;
    //using mesh_type = disk::simplicial_mesh<T,3>;
    using point_type = typename mesh_type::point_type;

    mesh_type msh;
    auto mesher = get_mesher(msh);
    make_mesh_level(mesher, reflevel);

    if (do_cook) {
        auto tr = [](const point_type& pt) {
            auto y = 0.044*pt.x() + pt.y()*( 0.044*(1.0-pt.x()) + 0.016*pt.x() );
            return point_type{0.048*pt.x(), y};
        };
        msh.transform(tr);
    }
    

    const static size_t DIM = mesh_type::dimension;
    auto fbs = disk::vector_basis_size(degree, DIM-1, DIM);

    

    std::vector<bc> bcs;
    set_boundary(msh, bcs, bc::neumann, 0);
    set_boundary(msh, bcs, bc::neumann, 1);
    set_boundary(msh, bcs, bc::neumann, 2);
    set_boundary(msh, bcs, bc::dirichlet, 3);

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
        auto S = vstab(msh, cl, degree, mode, bcs);

        MT lhs = 2*mu*Asgr + lambda*Adr + 2*mu*S;
        auto cb = disk::make_vector_monomial_basis(msh, cl, degree+1);
        auto cbs = cb.size();
        VT rhs = VT::Zero(lhs.rows());

        auto qps = disk::integrate(msh, cl, 2*degree+2);
        disk::static_vector<T,DIM> f =
            disk::static_vector<T,DIM>::Zero();
        f(0) = fx;
        f(1) = fy;

        /*        
        f *= 10;
        for (auto& qp : qps) {
            auto phi = cb.eval_functions(qp.point());
            rhs.head(cbs) += qp.weight() * (phi * f);
        }
        */
        
        auto fcs = faces(msh, cl);
        for (size_t fcnum = 0; fcnum < fcs.size(); fcnum++) {
            const auto& fc = fcs[fcnum];
            auto fb = disk::make_vector_monomial_basis(msh, fc, degree);
            auto fcofs = cbs + fb.size()*fcnum;

            auto bi = msh.boundary_info(fc);
            if (bi.is_boundary() and bi.id() == 1) {
                if (mode == hho_mode::standard) {
                    auto fqps = disk::integrate(msh, fc, 2*degree);
                    for (auto& qp : fqps) {
                        auto phi = fb.eval_functions(qp.point());
                        rhs.segment(fcofs, fb.size()) += qp.weight() * (phi * f);
                    }
                } else {
                    auto fqps = disk::integrate(msh, fc, 2*degree+2);
                    for (auto& qp : fqps) {
                        auto phi = cb.eval_functions(qp.point());
                        rhs.head(cb.size()) += qp.weight() * (phi * f);
                    }
                }
            }
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
        auto cbs = disk::scalar_basis_size(degree+1, DIM);
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
    silo.create("linelast.silo");
    silo.add_mesh(msh, "mesh");
    silo.add_variable("mesh", "u", u_data, disk::zonal_variable_t);

    return 0;
}