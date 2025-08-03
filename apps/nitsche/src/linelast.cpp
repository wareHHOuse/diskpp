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
class condensed_assembler {
    using CoordT = typename Mesh::coordinate_type;
    using ScalT = CoordT;

    using mesh_type = Mesh;
    using cell_type = typename mesh_type::cell_type;
    using face_type = typename mesh_type::face_type;
    using trip_type = typename Eigen::Triplet<ScalT>;

    std::vector<size_t>             compress_table;
    std::vector<size_t>             expand_table;
    std::vector<bool>               dirichlet_faces;

    std::vector<trip_type>          triplets;

    size_t                          num_face_dofs;
    size_t                          num_all_faces;
    size_t                          num_dirichlet_faces;
    size_t                          num_other_faces;
    size_t                          system_size;

    bool is_dirichlet(const mesh_type& msh, const face_type& fc) const {
        if (dirichlet_faces.size() == 0) {
            auto bi = msh.boundary_info(fc);
            return bi.is_boundary() and not bi.is_internal();
        }

        const auto fc_id = msh.lookup(fc);
        assert(fc_id < dirichlet_faces.size());
        return dirichlet_faces[fc_id];
    }

public:
    using matrix_type = disk::dynamic_matrix<ScalT>;
    using vector_type = disk::dynamic_vector<ScalT>;

    Eigen::SparseMatrix<ScalT>      LHS;
    vector_type                     RHS;

    condensed_assembler(const mesh_type& msh, size_t nfd)
        : num_face_dofs(nfd) 
    {
        using namespace std::placeholders;
        num_all_faces = msh.faces_size();
        dirichlet_faces.resize( num_all_faces );

        size_t fnum = 0;
        for (auto& fc : faces(msh)) {
            auto bi = msh.boundary_info(fc);
            if (bi.is_boundary() and not bi.is_internal())
                dirichlet_faces[fnum] = true;
            fnum++;
        }
        
        num_dirichlet_faces = std::count(dirichlet_faces.begin(), dirichlet_faces.end(), true);
        num_other_faces = num_all_faces - num_dirichlet_faces;

        compress_table.resize(num_all_faces);
        expand_table.resize(num_other_faces);

        size_t compressed_offset = 0;
        for (size_t i = 0; i < num_all_faces; i++) {
            const auto fc = *std::next(msh.faces_begin(), i);
            if (!is_dirichlet(msh, fc)) {
                compress_table.at(i) = compressed_offset;
                expand_table.at(compressed_offset) = i;
                compressed_offset++;
            }
        }

        const auto fbs = num_face_dofs;
        system_size = fbs * num_other_faces;

        LHS = Eigen::SparseMatrix<ScalT>(system_size, system_size);
        RHS = vector_type::Zero(system_size);
    }

    condensed_assembler(const mesh_type& msh, size_t nfd, std::vector<bool> dirichlet)
        : num_face_dofs(nfd), dirichlet_faces(dirichlet)
    {
        using namespace std::placeholders;
        num_all_faces = msh.faces_size();
        num_dirichlet_faces = std::count(dirichlet_faces.begin(), dirichlet_faces.end(), true);
        num_other_faces = num_all_faces - num_dirichlet_faces;

        compress_table.resize(num_all_faces);
        expand_table.resize(num_other_faces);

        size_t compressed_offset = 0;
        for (size_t i = 0; i < num_all_faces; i++) {
            const auto fc = *std::next(msh.faces_begin(), i);
            if (!is_dirichlet(msh, fc)) {
                compress_table.at(i) = compressed_offset;
                expand_table.at(compressed_offset) = i;
                compressed_offset++;
            }
        }

        const auto fbs = num_face_dofs;
        system_size = fbs * num_other_faces;

        LHS = Eigen::SparseMatrix<ScalT>(system_size, system_size);
        RHS = vector_type::Zero(system_size);
    }

    void
    assemble(const mesh_type& msh, const cell_type& cl, const matrix_type& lhs,
             const vector_type& rhs)
    {
        const auto fbs = num_face_dofs;
        const auto fcs = faces(msh, cl);
        const auto fcs_id = faces_id(msh, cl);
        assert(fbs*fcs.size() == lhs.rows());
        assert(fbs*fcs.size() == lhs.cols());
        assert(fbs*fcs.size() == rhs.rows());

        auto bases = [&](size_t fi, size_t fj) {
            assert(fi < fcs_id.size() and fj < fcs_id.size());
            auto idi = fcs_id[fi]; assert(idi < compress_table.size());
            auto idj = fcs_id[fj]; assert(idj < compress_table.size());
            auto gi = compress_table[idi];
            auto gj = compress_table[idj]; 
            return std::tuple(fi*fbs, fj*fbs, gi*fbs, gj*fbs);
        };

        for (size_t lnum_i = 0; lnum_i < fcs.size(); lnum_i++) {
            if (dirichlet_faces[ fcs_id[lnum_i] ])
                continue;
            
            auto gnum_i = compress_table[ fcs_id[lnum_i] ];

            for (size_t lnum_j = 0; lnum_j < fcs.size(); lnum_j++) {
                if (dirichlet_faces[ fcs_id[lnum_j] ]) {
                    //RHS.segment(gnum_i*fbs, fbs) -=
                    //    lhs.block(lnum_i*fbs, lnum_j*fbs, fbs, fbs) * dirichlet_data.segment(lnum_i*fbs, fbs);
                }
                else {
                    auto gnum_j = compress_table[ fcs_id[lnum_j] ];
                    for (size_t i = 0; i < fbs; i++) {
                        auto li = lnum_i*fbs + i;
                        auto gi = gnum_i*fbs + i;
                        for (size_t j = 0; j < fbs; j++) {
                            auto lj = lnum_j*fbs + j;
                            auto gj = gnum_j*fbs + j;
                            triplets.push_back( {gi, gj, lhs(li, lj)} );
                        } // for j
                    } // for i
                } // else
            } // for lnum_j

            RHS.segment(gnum_i*fbs, fbs) += rhs.segment(lnum_i*fbs, fbs);
        } // for lnum_i
    }

    vector_type
    take_local_solution(const mesh_type& msh, const cell_type& cl, const vector_type& solution)
    {
        const auto fbs = num_face_dofs;
        const auto fcs = faces(msh, cl);
        const auto fcs_id = faces_id(msh, cl);
        const auto num_faces = fcs.size();
        vector_type ret = vector_type::Zero(num_faces * fbs);
        
        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            const auto fc = fcs[face_i];
            const bool is_dirichlet = dirichlet_faces[fcs_id[face_i]];
        
            if (is_dirichlet)
            {}
            else
            {
                const auto face_SOL_offset = compress_table.at(fcs_id[face_i])*fbs;
                ret.segment(face_i*fbs, fbs) = solution.segment(face_SOL_offset, fbs);
            }
        }

        return ret;
    }

    auto dofs() const {
        return system_size;
    }

    void
    finalize(void)
    {
        LHS.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();
    }
};




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

int main(int argc, char **argv) {

    using T = double;
    //using mesh_type = disk::simplicial_mesh<T,2>;
    //using mesh_type = disk::cartesian_mesh<T,2>;
    using mesh_type = disk::simplicial_mesh<T,3>;

    mesh_type msh;
    auto mesher = make_simple_mesher(msh);
    mesher.refine();
    mesher.refine();
    mesher.refine();
    mesher.refine();
    //mesher.refine();

    size_t degree = 2;

    const static size_t DIM = mesh_type::dimension;
    auto fbs = disk::vector_basis_size(degree, DIM-1, DIM);

    hho_mode mode = hho_mode::nitsche;

    std::vector<bc> bcs;
    set_boundary(msh, bcs, bc::dirichlet, 0);
    set_boundary(msh, bcs, bc::neumann, 1);
    set_boundary(msh, bcs, bc::neumann, 2);
    set_boundary(msh, bcs, bc::neumann, 3);

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

    
    auto mu = std::stod(argv[1]);
    auto lambda = std::stod(argv[2]);

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
        f(0) = 0.1;
        f(1) = 0.1;

        
        f *= 10;
        for (auto& qp : qps) {
            auto phi = cb.eval_functions(qp.point());
            rhs.head(cbs) += qp.weight() * (phi * f);
        }
        
        /*
        f *= 10;
        auto fcs = faces(msh, cl);
        for (size_t fcnum = 0; fcnum < fcs.size(); fcnum++) {
            const auto& fc = fcs[fcnum];
            auto fb = disk::make_vector_monomial_basis(msh, fc, degree);
            auto fcofs = cbs + fb.size()*fcnum;

            auto bi = msh.boundary_info(fc);
            if (bi.id() == 2) {
                auto fqps = disk::integrate(msh, fc, 2*degree);
                for (auto& qp : fqps) {
                    auto phi = fb.eval_functions(qp.point());
                    rhs.segment(fcofs, fb.size()) += qp.weight() * (phi * f);
                }
            }
        }
        */

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