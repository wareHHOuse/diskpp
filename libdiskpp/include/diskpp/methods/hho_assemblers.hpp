/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#pragma once

#include <functional>

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/bases/bases_new.hpp"
#include "diskpp/bases/bases_operations.hpp"
#include "diskpp/bases/bases_traits.hpp"

namespace disk::hho {

template<typename Mesh, disk::basis::basis FaceBasis>
class basic_condensed_assembler {
    using CoordT = typename Mesh::coordinate_type;
    using ScalT = typename FaceBasis::scalar_type;

    using mesh_type = Mesh;
    using cell_type = typename mesh_type::cell_type;
    using face_type = typename mesh_type::face_type;
    using trip_type = typename Eigen::Triplet<ScalT>;

    std::vector<size_t>             compress_table;
    std::vector<size_t>             expand_table;
    std::vector<bool>               dirichlet_faces;

    std::vector<trip_type>          triplets;

    size_t                          face_poly_degree;
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
    using matrix_type = dynamic_matrix<ScalT>;
    using vector_type = dynamic_vector<ScalT>;

    Eigen::SparseMatrix<ScalT>      LHS;
    vector_type                     RHS;

    basic_condensed_assembler(const mesh_type& msh, size_t deg)
        : face_poly_degree(deg) 
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

        const auto fbs = FaceBasis::size_of_degree(face_poly_degree);
        system_size = fbs * num_other_faces;

        LHS = Eigen::SparseMatrix<ScalT>(system_size, system_size);
        RHS = vector_type::Zero(system_size);
    }

    basic_condensed_assembler(const mesh_type& msh, size_t deg, std::vector<bool> dirichlet)
        : face_poly_degree(deg), dirichlet_faces(dirichlet)
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

        const auto fbs = FaceBasis::size_of_degree(face_poly_degree);
        system_size = fbs * num_other_faces;

        LHS = Eigen::SparseMatrix<ScalT>(system_size, system_size);
        RHS = vector_type::Zero(system_size);
    }

    void
    assemble(const mesh_type& msh, const cell_type& cl, const matrix_type& lhs,
             const vector_type& rhs)
    {
        const auto fbs = FaceBasis::size_of_degree(face_poly_degree);
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
        const auto fbs = FaceBasis::size_of_degree(face_poly_degree);
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

    auto dirichlet_faces_flags() const {
        return dirichlet_faces;
    }
};

template<typename Basis, typename T>
auto
schur(const dynamic_matrix<T>& L, const dynamic_vector<T>& R, const Basis& basis)
{
    assert(L.cols() == L.rows());
    assert( (L.cols() == R.rows()) or (R.rows() == basis.size()) );

    using cdm = const dynamic_matrix<T>;
    using cdv = const dynamic_vector<T>;
    auto ts = basis.size();
    auto fs = L.cols() - basis.size();
    cdm MTT = L.topLeftCorner(ts, ts);
    cdm MTF = L.topRightCorner(ts, fs);
    cdm MFT = L.bottomLeftCorner(fs, ts);
    cdm MFF = L.bottomRightCorner(fs, fs);

    const auto MTT_ldlt = MTT.ldlt();
    if (MTT_ldlt.info() != Eigen::Success)
        throw std::invalid_argument("Schur: MTT is not positive definite");

    cdv RT = R.head(ts);
    cdm A = MTT_ldlt.solve(MTF);
    cdv B = MTT_ldlt.solve(RT);
    cdm Lc = MFF - MFT*A;

    if ( L.cols() == R.rows() ) {
        cdv RF = R.tail(fs);
        cdv Rc = RF - MFT*B;
        return std::pair(Lc, Rc);
    }

    cdv Rc = -MFT*B;
    return std::pair(Lc, Rc);
}

template<typename Basis, typename T>
auto
deschur(const dynamic_matrix<T>& lhs, const dynamic_vector<T>& rhs,
    const dynamic_vector<T>& solF, const Basis& basis)
{
    assert(lhs.cols() == lhs.rows());
    assert( (lhs.cols() == rhs.rows()) or (rhs.rows() == basis.size()) );

    using cdm = const dynamic_matrix<T>;
    using cdv = const dynamic_vector<T>;
    using dv = dynamic_vector<T>;
    auto ts = basis.size();
    auto fs = lhs.cols() - basis.size();
    cdm MTT = lhs.topLeftCorner(ts, ts);
    cdm MTF = lhs.topRightCorner(ts, fs);

    const auto MTT_ldlt = MTT.ldlt();
    if (MTT_ldlt.info() != Eigen::Success)
        throw std::invalid_argument("Schur: MTT is not positive definite");

    cdv rhsT = rhs.head(ts);
    cdv solT = MTT.ldlt().solve(rhsT - MTF*solF);

    dv ret = cdv::Zero(lhs.rows());
    ret.head(ts) = solT;
    ret.tail(fs) = solF;

    return ret;
}

};

