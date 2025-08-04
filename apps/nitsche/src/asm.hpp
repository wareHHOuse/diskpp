#pragma once

/***************************************************************
 * Boundary conditions helpers
 */

enum bc {
    none,
    dirichlet,
    neumann,
};

template<typename Mesh>
void
set_boundary(const Mesh& msh, std::vector<bc>& bcs, bc bc_type, size_t bnd)
{
    if (bcs.size() != msh.faces_size()) {
        bcs.resize( msh.faces_size() );
    }

    size_t fcnum = 0;
    for (auto& fc : faces(msh)) {
        auto bi = msh.boundary_info(fc);
        if (bi.is_boundary() and bi.id() == bnd) {
            bcs[fcnum] = bc_type;
        }
        fcnum++;
    }   
}

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
