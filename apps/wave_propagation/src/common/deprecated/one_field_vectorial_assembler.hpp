//
//  one_field_vectorial_assembler.hpp
//  acoustics
//
//  Created by Omar Durán on 4/11/20.
//

#pragma once
#ifndef one_field_vectorial_assembler_hpp
#define one_field_vectorial_assembler_hpp

#include "diskpp/bases/bases.hpp"
#include "diskpp/methods/hho"

template<typename Mesh>
class one_field_vectorial_assembler
{
    
    
    typedef disk::BoundaryConditions<Mesh, false>    boundary_type;
    using T = typename Mesh::coordinate_type;

    std::vector<size_t>                 m_compress_indexes;
    std::vector<size_t>                 m_expand_indexes;

    disk::hho_degree_info               m_hho_di;
    boundary_type                       m_bnd;
    std::vector< Triplet<T> >           m_triplets;

    size_t      m_n_edges;
    size_t      m_n_essential_edges;
    bool        m_hho_stabilization_Q;

    class assembly_index
    {
        size_t  idx;
        bool    assem;

    public:
        assembly_index(size_t i, bool as)
            : idx(i), assem(as)
        {}

        operator size_t() const
        {
            if (!assem)
                throw std::logic_error("Invalid assembly_index");

            return idx;
        }

        bool assemble() const
        {
            return assem;
        }

        friend std::ostream& operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };

public:

    SparseMatrix<T>         LHS;
    Matrix<T, Dynamic, 1>   RHS;

    one_field_vectorial_assembler(const Mesh& msh, const disk::hho_degree_info& hho_di, const boundary_type& bnd)
        : m_hho_di(hho_di), m_bnd(bnd), m_hho_stabilization_Q(true)
    {
            
        auto is_dirichlet = [&](const typename Mesh::face& fc) -> bool {

            auto fc_id = msh.lookup(fc);
            return bnd.is_dirichlet_face(fc_id);
        };

        m_n_edges = msh.faces_size();
        m_n_essential_edges = std::count_if(msh.faces_begin(), msh.faces_end(), is_dirichlet);

        m_compress_indexes.resize( m_n_edges );
        m_expand_indexes.resize( m_n_edges - m_n_essential_edges );

        size_t compressed_offset = 0;
        for (size_t i = 0; i < m_n_edges; i++)
        {
            auto fc = *std::next(msh.faces_begin(), i);
            if ( !is_dirichlet(fc) )
            {
                m_compress_indexes.at(i) = compressed_offset;
                m_expand_indexes.at(compressed_offset) = i;
                compressed_offset++;
            }
        }

        size_t n_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);

        size_t system_size = n_cbs * msh.cells_size() + n_fbs * (m_n_edges - m_n_essential_edges);

        LHS = SparseMatrix<T>( system_size, system_size );
        RHS = Matrix<T, Dynamic, 1>::Zero( system_size );
    }

    void scatter_data(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs,
             const Matrix<T, Dynamic, 1>& rhs)
    {
        auto fcs = faces(msh, cl);
        size_t n_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs + n_fbs*fcs.size());

        auto cell_offset        = disk::priv::offset(msh, cl);
        auto cell_LHS_offset    = cell_offset * n_cbs;

        for (size_t i = 0; i < n_cbs; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
        
        Matrix<T, Dynamic, 1> dirichlet_data = Matrix<T, Dynamic, 1>::Zero(n_cbs + fcs.size()*n_fbs);
        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto face_offset = disk::priv::offset(msh, fc);
            auto face_LHS_offset = n_cbs * msh.cells_size() + m_compress_indexes.at(face_offset)*n_fbs;

            auto fc_id = msh.lookup(fc);
            bool dirichlet = m_bnd.is_dirichlet_face(fc_id);

            for (size_t i = 0; i < n_fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
            
            if (dirichlet)
             {
                 auto fb = make_vector_monomial_basis(msh, fc, m_hho_di.face_degree());
                 auto dirichlet_fun  = m_bnd.dirichlet_boundary_func(fc_id);

                 Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb);
                 Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_fun);
                 dirichlet_data.block(n_cbs + face_i*n_fbs, 0, n_fbs, 1) = mass.llt().solve(rhs);
             }
            
        }

        assert( asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if ( asm_map[j].assemble() )
                    m_triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
                else
                    RHS(asm_map[i]) -= lhs(i,j) * dirichlet_data(j);
            }
        }

        for (size_t i = 0; i < rhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;
            RHS(asm_map[i]) += rhs(i);
        }

    } //

    void assemble(const Mesh& msh, std::function<static_vector<double, 2>(const typename Mesh::point_type& )> rhs_fun){
        
        for (auto& cell : msh)
        {
            auto reconstruction_operator   = make_matrix_symmetric_gradrec(msh, cell, m_hho_di);
            auto rec_for_stab   = make_vector_hho_symmetric_laplacian(msh, cell, m_hho_di);
            auto divergence_operator = make_hho_divergence_reconstruction(msh, cell, m_hho_di);
            
            Matrix<T, Dynamic, Dynamic> R_operator = reconstruction_operator.second;
            

            Matrix<T, Dynamic, Dynamic> S_operator;
            if(m_hho_stabilization_Q)
            {
                auto stabilization_operator    = make_vector_hho_stabilization(msh, cell, rec_for_stab.first, m_hho_di);
                S_operator = stabilization_operator;
            }else{
                auto stabilization_operator    = make_vector_hdg_stabilization(msh, cell, m_hho_di);
                S_operator = stabilization_operator;
            }
            
            Matrix<T, Dynamic, Dynamic> vectorial_laplacian_operator_loc = 2.0*(R_operator + S_operator)+divergence_operator.second;
            auto cell_basis   = make_vector_monomial_basis(msh, cell, m_hho_di.cell_degree());
            Matrix<T, Dynamic, 1> f_loc = make_rhs(msh, cell, cell_basis, rhs_fun);
            
            scatter_data(msh, cell, vectorial_laplacian_operator_loc, f_loc);
        }
        finalize();
    }

    void finalize(void)
    {
        LHS.setFromTriplets( m_triplets.begin(), m_triplets.end() );
        m_triplets.clear();
    }

    Matrix<T, Dynamic, 1>
    gather_dof_data(  const Mesh& msh, const typename Mesh::cell_type& cl,
                    const Matrix<T, Dynamic, 1>& x_glob) const
    {
        auto num_faces = howmany_faces(msh, cl);
        auto cell_ofs = disk::priv::offset(msh, cl);
        size_t n_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        
        Matrix<T, Dynamic, 1> x_el(n_cbs + num_faces * n_fbs );
        x_el.block(0, 0, n_cbs, 1) = x_glob.block(cell_ofs * n_cbs, 0, n_cbs, 1);
        auto fcs = faces(msh, cl);
        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto fc = fcs[i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
            const auto face_id                  = eid.second;

            if (m_bnd.is_dirichlet_face( face_id))
            {
                auto fb = make_vector_monomial_basis(msh, fc, m_hho_di.face_degree());
                auto dirichlet_fun  = m_bnd.dirichlet_boundary_func(face_id);
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb);
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_fun);
                x_el.block(n_cbs + i * n_fbs, 0, n_fbs, 1) = mass.llt().solve(rhs);
            }
            else
            {
                auto face_ofs = disk::priv::offset(msh, fc);
                auto global_ofs = n_cbs * msh.cells_size() + m_compress_indexes.at(face_ofs)*n_fbs;
                x_el.block(n_cbs + i*n_fbs, 0, n_fbs, 1) = x_glob.block(global_ofs, 0, n_fbs, 1);
            }
        }
        return x_el;
    }
            
    void set_hdg_stabilization(){
        if(m_hho_di.cell_degree() > m_hho_di.face_degree())
        {
            m_hho_stabilization_Q = false;
            std::cout << "Proceeding with HDG stabilization cell degree is higher than face degree." << std::endl;
            std::cout << "cell degree = " << m_hho_di.cell_degree() << std::endl;
            std::cout << "face degree = " << m_hho_di.face_degree() << std::endl;
        }else{
            std::cout << "Proceeding with HHO stabilization cell and face degree are equal." << std::endl;
            std::cout << "cell degree = " << m_hho_di.cell_degree() << std::endl;
            std::cout << "face degree = " << m_hho_di.face_degree() << std::endl;
        }
    }
    void set_hho_stabilization(){
        m_hho_stabilization_Q = true;
    }
};

#endif /* one_field_vectorial_assembler_hpp */
