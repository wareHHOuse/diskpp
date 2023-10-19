#pragma once


#include "diskpp/methods/hho"

template<typename Mesh>
class hho_assembler
{
    using coordinate_type = typename Mesh::coordinate_type;
    using scalar_type = coordinate_type;
    using mesh_type = Mesh;
    using face_type = typename Mesh::face_type;
    using cell_type = typename Mesh::cell_type;

    //Mesh                                msh;
    disk::hho_degree_info               hdi;

    std::vector<Triplet<scalar_type>>   triplets;
    std::vector<bool>                   is_dirichlet;

    std::vector<size_t>                 compress_table;
    std::vector<size_t>                 expand_table;

    const size_t INVALID_OFFSET = (size_t) ~0;

    size_t face_basis_size(void) const {
        return disk::scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    }

    /* Get the offset of a face in the linear system */
    size_t get_system_offset(const Mesh& msh,
        const typename Mesh::face_type& fc)
    {
        auto face_num = offset(msh, fc);
        assert(face_num < msh.faces_size());
        assert(compress_table.size() == msh.faces_size());
        auto cnum = compress_table[face_num];
        assert(cnum != INVALID_OFFSET);
        auto fbs = face_basis_size();
        return cnum*fbs;
    }

    /* Determine if a face should be assembled */
    bool is_in_system(const Mesh& msh,
        const typename Mesh::face_type& fc)
    {
        auto bi = msh.boundary_info(fc);
        auto face_num = offset(msh, fc);
        assert(face_num < msh.faces_size());
        assert(is_dirichlet.size() == msh.faces_size());
        return not (bi.is_boundary() and is_dirichlet[face_num]);
    }

    void make_tables(mesh_type& msh)
    {
        compress_table.resize( msh.faces_size(), INVALID_OFFSET);
        expand_table.resize( sysfcs );

        size_t face_i = 0;
        size_t compressed_ofs = 0;
        for (auto& fc : faces(msh))
        {
            assert(compressed_ofs <= face_i);
            if ( is_in_system(msh, fc) )
            {
                assert(face_i < compress_table.size());
                compress_table[face_i] = compressed_ofs;
                assert(compressed_ofs < expand_table.size());
                expand_table[compressed_ofs] = face_i;
                compressed_ofs++;
            }

            face_i++;
        }

        assert(face_i == msh.faces_size());
    }

public:
    
    SparseMatrix<scalar_type>           LHS;
    Matrix<scalar_type, Dynamic, 1>     RHS;

    size_t                              syssz;
    size_t                              sysfcs;


    hho_assembler()
    {}

    hho_assembler(Mesh& msh, const disk::hho_degree_info& hdi,
                  const std::vector<bool>& is_dirichlet)
    {
        initialize(msh, hdi, is_dirichlet);
    }

    void clear()
    {
        triplets.clear();
        is_dirichlet.clear();
        compress_table.clear();
        expand_table.clear();
        syssz = 0;
        sysfcs = 0;
    }

    void initialize(Mesh& msh, const disk::hho_degree_info& p_hdi,
                    const std::vector<bool>& p_is_dirichlet)
    {
        clear();

        is_dirichlet = p_is_dirichlet;
        hdi = p_hdi;

        auto fbs = face_basis_size();

        auto in_system = [&](const face_type& fc) -> bool {
            auto ofs = offset(msh, fc);
            assert(ofs < is_dirichlet.size());
            return not (msh.is_boundary(fc) and is_dirichlet[ofs]);
        };

        sysfcs = std::count_if(msh.faces_begin(), msh.faces_end(), in_system);
        syssz = fbs*sysfcs;

        make_tables(msh);

        LHS = SparseMatrix<scalar_type>(syssz, syssz);
        RHS = Matrix<scalar_type, Dynamic, 1>::Zero(syssz);

        std::cout << "Assembler initialized: " << sysfcs << " faces in system, ";
        std::cout << syssz << " DoFs" << std::endl;
    }

    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<scalar_type, Dynamic, Dynamic>& lhsc,
             const Matrix<scalar_type, Dynamic, 1>& rhs,
             const Matrix<scalar_type, Dynamic, 1>& dirichlet_data)
    {
        auto fbs = face_basis_size();

        auto fcs = faces(msh, cl);
        for (size_t fi = 0; fi < fcs.size(); fi++)
        {
            if ( not is_in_system(msh, fcs[fi]) )
                continue;

            auto cofsi = get_system_offset(msh, fcs[fi]);
            for (size_t fj = 0; fj < fcs.size(); fj++)
            {
                auto lofsi = fi*fbs;
                auto lofsj = fj*fbs;

                if ( not is_in_system(msh, fcs[fj]) ) 
                {
                    RHS.segment(cofsi, fbs) += -lhsc.block(lofsi, lofsj, fbs, fbs)*dirichlet_data.segment(lofsj, fbs);
                    continue;
                }

                auto cofsj = get_system_offset(msh, fcs[fj]);
                for (size_t i = 0; i < fbs; i++)
                    for(size_t j = 0; j < fbs; j++)
                        triplets.push_back( Triplet<scalar_type>(cofsi+i, cofsj+j, lhsc(lofsi+i, lofsj+j)) );
            }

            RHS.segment(cofsi, fbs) += rhs.segment(fi*fbs, fbs);
        }
    }

    void finalize(void)
    {
        LHS.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();
        std::cout << "Matrix has " << LHS.nonZeros() << " nonzeros." << std::endl; 
    }

    disk::dynamic_vector<scalar_type>
    get_expanded_solution(const Mesh& msh, disk::dynamic_vector<scalar_type>& sol)
    {
        auto fbs = face_basis_size();

        disk::dynamic_vector<scalar_type> ret = 
            disk::dynamic_vector<scalar_type>::Zero( fbs*msh.faces_size() );

        for (size_t i = 0; i < sysfcs; i++)
        {
            auto in_offset = i*fbs;
            auto out_offset = expand_table.at(i)*fbs;
            ret.segment(out_offset, fbs) = sol.segment(in_offset, fbs);
        }

        return ret;
    }

    disk::dynamic_vector<scalar_type>
    get_element_dofs(const Mesh& msh, const typename Mesh::cell& cl, disk::dynamic_vector<scalar_type>& sol)
    {
        auto fbs = face_basis_size();
        auto fcs = faces(msh, cl);
        disk::dynamic_vector<scalar_type> ret = 
            disk::dynamic_vector<scalar_type>::Zero( fbs*fcs.size() );

        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto& fc = fcs[i];
            if ( not is_in_system(msh, fc) )
                continue;

            auto ofs = get_system_offset(msh, fc);
            ret.segment(i*fbs, fbs) = sol.segment(ofs, fbs);
        }

        return ret;
    }
};


