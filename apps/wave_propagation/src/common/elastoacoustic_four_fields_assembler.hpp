//
//  elastoacoustic_four_fields_assembler.hpp
//  acoustics
//
//  Created by Omar Durán on 9/7/20.
//  


#pragma once
#ifndef elastoacoustic_four_fields_assembler_hpp
#define elastoacoustic_four_fields_assembler_hpp

#include "diskpp/bases/bases.hpp"
#include "diskpp/methods/hho"
#include "../common/assembly_index.hpp"
#include "../common/acoustic_material_data.hpp"
#include "../common/elastic_material_data.hpp"
#include <map>

#ifdef HAVE_INTEL_TBB
#include <tbb/parallel_for.h>
#endif

template<typename Mesh>
class elastoacoustic_four_fields_assembler
{
    
    typedef disk::BoundaryConditions<Mesh, false>    e_boundary_type;
    typedef disk::BoundaryConditions<Mesh, true>     a_boundary_type;
    using T = typename Mesh::coordinate_type;

    std::vector<size_t>                 m_e_compress_indexes;
    std::vector<size_t>                 m_e_expand_indexes;
    
    std::vector<size_t>                 m_a_compress_indexes;
    std::vector<size_t>                 m_a_expand_indexes;

    disk::hho_degree_info               m_hho_di;
    e_boundary_type                     m_e_bnd;
    a_boundary_type                     m_a_bnd;
    std::vector< Triplet<T> >           m_triplets;
    std::vector< Triplet<T> >           m_triplets_stab;
    std::vector< Triplet<T> >           m_c_triplets;
    std::vector< Triplet<T> >           m_mass_triplets;
    std::map<size_t,elastic_material_data<T>> m_e_material;
    std::map<size_t,acoustic_material_data<T>> m_a_material;
    std::map<size_t,size_t> m_e_cell_index;
    std::map<size_t,size_t> m_a_cell_index;
    std::vector< size_t >               m_e_elements_with_bc_eges;
    std::vector< size_t >               m_a_elements_with_bc_eges;
    std::map<size_t,std::pair<size_t,size_t>>   m_interface_cell_indexes;

    size_t      m_n_elastic_cells;
    size_t      m_n_acoustic_cells;

    size_t      n_e_edges;
    size_t      n_a_edges;

    size_t      m_n_edges;
    size_t      m_n_essential_edges;
    size_t      m_n_elastic_cell_dof;
    size_t      m_n_acoustic_cell_dof;
    size_t      m_n_elastic_face_dof;
    size_t      m_n_acoustic_face_dof;
    bool        m_hho_stabilization_Q;
    bool        m_scaled_stabilization_Q;

public:

    SparseMatrix<T>         LHS;
    SparseMatrix<T>         LHS_STAB;
    Matrix<T, Dynamic, 1>   RHS;
    SparseMatrix<T>         MASS;
    SparseMatrix<T>         COUPLING;

    elastoacoustic_four_fields_assembler(const Mesh& msh, const disk::hho_degree_info& hho_di, const e_boundary_type& e_bnd, const a_boundary_type& a_bnd, std::map<size_t,elastic_material_data<T>> & e_material, std::map<size_t,acoustic_material_data<T>> & a_material)
        : m_hho_di(hho_di), m_e_bnd(e_bnd), m_a_bnd(a_bnd), m_e_material(e_material), m_a_material(a_material), m_hho_stabilization_Q(true), m_scaled_stabilization_Q(false)
    {
            
        auto storage = msh.backend_storage();
        auto is_e_dirichlet = [&](const typename Mesh::face& fc) -> bool {

            auto fc_id = msh.lookup(fc);
            return e_bnd.is_dirichlet_face(fc_id);
        };
        
        auto is_a_dirichlet = [&](const typename Mesh::face& fc) -> bool {

            auto fc_id = msh.lookup(fc);
            return a_bnd.is_dirichlet_face(fc_id);
        };

        size_t n_e_essential_edges = std::count_if(msh.faces_begin(), msh.faces_end(), is_e_dirichlet);
        size_t n_a_essential_edges = std::count_if(msh.faces_begin(), msh.faces_end(), is_a_dirichlet);
        
        std::set<size_t> e_egdes;
        for (auto &chunk : m_e_material) {
            size_t cell_i = chunk.first;
            auto& cell = storage->surfaces[cell_i];
            auto cell_faces = faces(msh,cell);
            for (auto &face : cell_faces) {
                if (!is_e_dirichlet(face)) {
                    auto fc_id = msh.lookup(face);
                    e_egdes.insert(fc_id);
                }
            }
        }
        n_e_edges = e_egdes.size();
        std::set<size_t> a_egdes;
        for (auto &chunk : m_a_material) {
            size_t cell_i = chunk.first;
            auto& cell = storage->surfaces[cell_i];
            auto cell_faces = faces(msh,cell);
            for (auto &face : cell_faces) {
                if (!is_a_dirichlet(face)) {
                    auto fc_id = msh.lookup(face);
                    a_egdes.insert(fc_id);
                }
            }
        }
        n_a_edges = a_egdes.size();
        
        m_n_edges = msh.faces_size();
        m_n_essential_edges = n_e_essential_edges + n_a_essential_edges;

        m_e_compress_indexes.resize( m_n_edges );
        m_e_expand_indexes.resize( m_n_edges - m_n_essential_edges );
        
        m_a_compress_indexes.resize( m_n_edges );
        m_a_expand_indexes.resize( m_n_edges - m_n_essential_edges );

        size_t e_compressed_offset = 0;
        for (auto face_id : e_egdes) {
          m_e_compress_indexes.at(face_id) = e_compressed_offset;
          m_e_expand_indexes.at(e_compressed_offset) = face_id;
          e_compressed_offset++;
        }
        size_t a_compressed_offset = 0;
        for (auto face_id : a_egdes) {
          m_a_compress_indexes.at(face_id) = a_compressed_offset;
          m_a_expand_indexes.at(a_compressed_offset) = face_id;
          a_compressed_offset++;
        }
    
        size_t n_cbs = get_e_cell_basis_data();
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        
        size_t n_s_cbs = get_a_cell_basis_data();
        size_t n_s_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);

        m_n_elastic_cell_dof = (n_cbs * m_e_material.size());
        m_n_acoustic_cell_dof = (n_s_cbs * m_a_material.size());
        
        m_n_elastic_face_dof = (n_fbs * n_e_edges);
        m_n_acoustic_face_dof = (n_s_fbs * n_a_edges);
        size_t system_size = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof + m_n_acoustic_face_dof;

        LHS = SparseMatrix<T>( system_size, system_size );
        LHS_STAB = SparseMatrix<T>( system_size, system_size ); //OPTIONAL
        RHS = Matrix<T, Dynamic, 1>::Zero( system_size );
        MASS = SparseMatrix<T>( system_size, system_size );
        COUPLING = SparseMatrix<T>( system_size, system_size );
        classify_cells(msh);
        build_cells_maps();
    }
    
    void scatter_e_data_stab(size_t e_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs,
             const Matrix<T, Dynamic, 1>& rhs)
    {
        auto fcs = faces(msh, cl);
        size_t n_cbs = get_e_cell_basis_data();
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs + n_fbs*fcs.size());

        auto cell_LHS_offset = e_cell_ind * n_cbs;

        for (size_t i = 0; i < n_cbs; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
            
        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto fc_id = msh.lookup(fc);
            auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_e_compress_indexes.at(fc_id)*n_fbs;
            bool dirichlet = m_e_bnd.is_dirichlet_face(fc_id);

            for (size_t i = 0; i < n_fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
        }
            
        assert( asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if ( asm_map[j].assemble() )
                    m_triplets_stab.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
            }
        }

    }
    

    void scatter_e_data(size_t e_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs,
             const Matrix<T, Dynamic, 1>& rhs)
    {
        auto fcs = faces(msh, cl);
        size_t n_cbs = get_e_cell_basis_data();
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs + n_fbs*fcs.size());

        auto cell_LHS_offset = e_cell_ind * n_cbs;

        for (size_t i = 0; i < n_cbs; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
            
        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto fc_id = msh.lookup(fc);
            auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_e_compress_indexes.at(fc_id)*n_fbs;
            bool dirichlet = m_e_bnd.is_dirichlet_face(fc_id);

            for (size_t i = 0; i < n_fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
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
            }
        }

        for (size_t i = 0; i < rhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;
            RHS(asm_map[i].vidx()) += rhs(i);
        }

    }
    
    void scatter_a_data(size_t a_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl,
    const Matrix<T, Dynamic, Dynamic>& lhs,
    const Matrix<T, Dynamic, 1>& rhs)
    {
        auto fcs = faces(msh, cl);
        size_t n_cbs = get_a_cell_basis_data();
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs + n_fbs*fcs.size());

        auto cell_LHS_offset    = a_cell_ind * n_cbs + m_n_elastic_cell_dof;

        for (size_t i = 0; i < n_cbs; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
        
        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto fc_id = msh.lookup(fc);
            auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof + m_a_compress_indexes.at(fc_id)*n_fbs;

            bool dirichlet = m_a_bnd.is_dirichlet_face(fc_id);

            for (size_t i = 0; i < n_fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
            
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
            }
        }

        for (size_t i = 0; i < rhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;
            RHS(asm_map[i].vidx()) += rhs(i);
        }

    }
    
    void scatter_a_data_stab(size_t a_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl,
    const Matrix<T, Dynamic, Dynamic>& lhs,
    const Matrix<T, Dynamic, 1>& rhs)
    {
        auto fcs = faces(msh, cl);
        size_t n_cbs = get_a_cell_basis_data();
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs + n_fbs*fcs.size());

        auto cell_LHS_offset    = a_cell_ind * n_cbs + m_n_elastic_cell_dof;

        for (size_t i = 0; i < n_cbs; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
        
        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto fc_id = msh.lookup(fc);
            auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof + m_a_compress_indexes.at(fc_id)*n_fbs;

            bool dirichlet = m_a_bnd.is_dirichlet_face(fc_id);

            for (size_t i = 0; i < n_fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
            
        }

        assert( asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if ( asm_map[j].assemble() )
                    m_triplets_stab.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
            }
        }
    }
    
    void scatter_e_mass_data(size_t e_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& mass_matrix)
    {
        
        size_t n_cbs = get_e_cell_basis_data();
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs);
        
        auto cell_LHS_offset    = e_cell_ind * n_cbs;

        for (size_t i = 0; i < n_cbs; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );

        assert( asm_map.size() == mass_matrix.rows() && asm_map.size() == mass_matrix.cols() );

        for (size_t i = 0; i < mass_matrix.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < mass_matrix.cols(); j++)
            {
                if ( asm_map[j].assemble() )
                    m_mass_triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], mass_matrix(i,j)) );
            }
        }

    }
    
    void scatter_a_mass_data(size_t a_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& mass_matrix)
    {
        size_t n_cbs = get_a_cell_basis_data();
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs);

        auto cell_LHS_offset    = a_cell_ind * n_cbs + m_n_elastic_cell_dof;
        
        for (size_t i = 0; i < n_cbs; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );

        assert( asm_map.size() == mass_matrix.rows() && asm_map.size() == mass_matrix.cols() );

        for (size_t i = 0; i < mass_matrix.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < mass_matrix.cols(); j++)
            {
                if ( asm_map[j].assemble() )
                    m_mass_triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], mass_matrix(i,j)) );
            }
        }

    }
    
    void scatter_ea_interface_data(const Mesh& msh, const typename Mesh::face_type& face, const Matrix<T, Dynamic, Dynamic>& interface_matrix) {
        
        auto vfbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        auto sfbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);

        std::vector<assembly_index> asm_map_e, asm_map_a;
        asm_map_e.reserve(vfbs);
        asm_map_a.reserve(sfbs);
        auto fc_id = msh.lookup(face);
        
        auto e_face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_e_compress_indexes.at(fc_id)*vfbs;
        bool e_dirichlet = m_e_bnd.is_dirichlet_face(fc_id);
        for (size_t i = 0; i < vfbs; i++){
            asm_map_e.push_back( assembly_index(e_face_LHS_offset+i, !e_dirichlet) );
        }
            
        auto a_face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof + m_a_compress_indexes.at(fc_id)*sfbs;
        bool a_dirichlet = m_a_bnd.is_dirichlet_face(fc_id);
        for (size_t i = 0; i < sfbs; i++){
            asm_map_a.push_back( assembly_index(a_face_LHS_offset+i, !a_dirichlet) );
        }
        
        assert( asm_map_e.size() == interface_matrix.rows() && asm_map_a.size() == interface_matrix.cols() );

        for (size_t i = 0; i < interface_matrix.rows(); i++) {
            for (size_t j = 0; j < interface_matrix.cols(); j++) {
              m_c_triplets.push_back( Triplet<T>(asm_map_e[i].vidx(), asm_map_a[j].vidx(),   interface_matrix(i,j)) );
              m_c_triplets.push_back( Triplet<T>(asm_map_a[j].vidx(), asm_map_e[i].vidx(), - interface_matrix(i,j)) );
            }
        }

    }
    
    void scatter_e_rhs_data(size_t e_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, 1>& rhs) {
    
        size_t n_cbs = get_e_cell_basis_data();
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs);

        auto cell_LHS_offset    = e_cell_ind * n_cbs;

        for (size_t i = 0; i < n_cbs; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );

        assert( asm_map.size() == rhs.rows() );

        for (size_t i = 0; i < rhs.rows(); i++) {
            if (!asm_map[i].assemble())
                continue;
            RHS(asm_map[i].vidx()) += rhs(i);
        }

    }
    
    void scatter_a_rhs_data(size_t a_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, 1>& rhs)
    {
        size_t n_cbs = get_a_cell_basis_data();
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs);

        auto cell_LHS_offset    = a_cell_ind * n_cbs + m_n_elastic_cell_dof;

        for (size_t i = 0; i < n_cbs; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );

        assert( asm_map.size() == rhs.rows());

        for (size_t i = 0; i < rhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;
            RHS(asm_map[i].vidx()) += rhs(i);
        }

    }
    
    void scatter_e_bc_data(size_t e_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs)
    {
    
        auto fcs = faces(msh, cl);
        size_t n_cbs = get_e_cell_basis_data();
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs + n_fbs*fcs.size());

        auto cell_LHS_offset    = e_cell_ind * n_cbs;

        for (size_t i = 0; i < n_cbs; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
        
        Matrix<T, Dynamic, 1> dirichlet_data = Matrix<T, Dynamic, 1>::Zero(n_cbs + fcs.size()*n_fbs);
        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto fc_id = msh.lookup(fc);
            auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_e_compress_indexes.at(fc_id)*n_fbs;

            bool dirichlet = m_e_bnd.is_dirichlet_face(fc_id);

            for (size_t i = 0; i < n_fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
            
            if (dirichlet)
             {
                 auto fb = make_vector_monomial_basis(msh, fc, m_hho_di.face_degree());
                 auto dirichlet_fun  = m_e_bnd.dirichlet_boundary_func(fc_id);

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
                if ( !asm_map[j].assemble() )
                    RHS(asm_map[i].vidx()) -= lhs(i,j) * dirichlet_data(j);
            }
        }

    }
    
    void scatter_a_bc_data(size_t a_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs)
    {
        auto fcs = faces(msh, cl);
        size_t n_cbs = get_a_cell_basis_data();
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs + n_fbs*fcs.size());
        
        auto cell_LHS_offset    = a_cell_ind * n_cbs + m_n_elastic_cell_dof;

        for (size_t i = 0; i < n_cbs; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
        
        Matrix<T, Dynamic, 1> dirichlet_data = Matrix<T, Dynamic, 1>::Zero(n_cbs + fcs.size()*n_fbs);
        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto fc_id = msh.lookup(fc);
            auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof + m_a_compress_indexes.at(fc_id)*n_fbs;
            
            bool dirichlet = m_a_bnd.is_dirichlet_face(fc_id);

            for (size_t i = 0; i < n_fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
            
            if (dirichlet)
             {
                 auto fb = make_scalar_monomial_basis(msh, fc, m_hho_di.face_degree());
                 auto dirichlet_fun  = m_a_bnd.dirichlet_boundary_func(fc_id);

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
                if ( !asm_map[j].assemble() )
                    RHS(asm_map[i].vidx()) -= lhs(i,j) * dirichlet_data(j);
                    
            }
        }

    }
    
    void assemble(const Mesh& msh, std::function<static_vector<T, 2>(const typename Mesh::point_type& )> e_rhs_fun, std::function<T(const typename Mesh::point_type& )> a_rhs_fun, bool explicit_scheme){
        
        auto storage = msh.backend_storage();
        LHS.setZero();
        RHS.setZero();
        
        // elastic block
        for (auto e_chunk : m_e_material) {
            size_t e_cell_ind = m_e_cell_index[e_chunk.first];
            auto& cell = storage->surfaces[e_chunk.first];
            Matrix<T, Dynamic, Dynamic> mixed_operator_loc = e_mixed_operator(e_chunk.second,msh,cell,explicit_scheme);
            Matrix<T, Dynamic, 1> f_loc = e_mixed_rhs(msh, cell, e_rhs_fun);
            scatter_e_data(e_cell_ind, msh, cell, mixed_operator_loc, f_loc);
        }
        
        // acoustic block
        for (auto a_chunk : m_a_material) {
            size_t a_cell_ind = m_a_cell_index[a_chunk.first];
            auto& cell = storage->surfaces[a_chunk.first];
            Matrix<T, Dynamic, Dynamic> mixed_operator_loc = a_mixed_operator(a_chunk.second, msh, cell,explicit_scheme);
            Matrix<T, Dynamic, 1> f_loc = a_mixed_rhs(msh, cell, a_rhs_fun);
            scatter_a_data(a_cell_ind, msh, cell, mixed_operator_loc, f_loc);
        }
        finalize();
    }
    

    void assemble_coupling_terms(const Mesh& msh){
        
        auto storage = msh.backend_storage();
        COUPLING.setZero();
        // coupling blocks
        for (auto chunk : m_interface_cell_indexes) {
            auto& face = storage->edges[chunk.first];
            auto& e_cell = storage->surfaces[chunk.second.first];
            auto& a_cell = storage->surfaces[chunk.second.second];
            Matrix<T, Dynamic, Dynamic> interface_operator_loc = e_interface_operator(msh, face, e_cell, a_cell);
            scatter_ea_interface_data(msh, face, interface_operator_loc);
        }
        finalize_coupling();
    }

    void assemble_mass(const Mesh& msh, bool add_vector_mass_Q = true){
        
        auto storage = msh.backend_storage();
        MASS.setZero();
        
        // elastic block
        for (auto e_chunk : m_e_material) {
            size_t e_cell_ind = m_e_cell_index[e_chunk.first];
            auto& cell = storage->surfaces[e_chunk.first];
            elastic_material_data<T> material = e_chunk.second;
            Matrix<T, Dynamic, Dynamic> mass_matrix = e_mass_operator(material, msh, cell, add_vector_mass_Q);
            scatter_e_mass_data(e_cell_ind,msh, cell, mass_matrix);
        }
        
        // acoustic block
        for (auto a_chunk : m_a_material) {
            size_t a_cell_ind = m_a_cell_index[a_chunk.first];
            auto& cell = storage->surfaces[a_chunk.first];
            acoustic_material_data<T> material = a_chunk.second;
            Matrix<T, Dynamic, Dynamic> mass_matrix = a_mass_operator(material, msh, cell);
            scatter_a_mass_data(a_cell_ind, msh, cell, mass_matrix);
        }
        
        finalize_mass();
    }
    
    void assemble_rhs(const Mesh& msh, std::function<static_vector<T, 2>(const typename Mesh::point_type& )> e_rhs_fun, std::function<T(const typename Mesh::point_type& )> a_rhs_fun, bool explicit_scheme){
        
        RHS.setZero();
        #ifdef HAVE_INTEL_TBB2
                size_t n_cells = msh.cells_size();
                tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                    [this,&msh,&rhs_fun] (size_t & cell_ind){
                        auto& cell = msh.backend_storage()->surfaces[cell_ind];
                        auto cell_basis   = make_vector_monomial_basis(msh, cell, m_hho_di.cell_degree());
                        Matrix<T, Dynamic, 1> f_loc = make_rhs(msh, cell, cell_basis, rhs_fun);
                        this->scatter_rhs_data(msh, cell, f_loc);
                }
            );
        #else
            auto storage = msh.backend_storage();
             for (auto e_chunk : m_e_material) {
                 size_t e_cell_ind = m_e_cell_index[e_chunk.first];
                 auto& cell = storage->surfaces[e_chunk.first];
                 Matrix<T, Dynamic, 1> f_loc = e_mixed_rhs(msh, cell, e_rhs_fun);
                 scatter_e_rhs_data(e_cell_ind, msh, cell, f_loc);
             }
        
            for (auto a_chunk : m_a_material) {
                size_t a_cell_ind = m_a_cell_index[a_chunk.first];
                auto& cell = storage->surfaces[a_chunk.first];
                Matrix<T, Dynamic, 1> f_loc = a_mixed_rhs(msh, cell, a_rhs_fun);
                scatter_a_rhs_data(a_cell_ind, msh, cell, f_loc);
            }
        
        #endif
        apply_bc(msh, explicit_scheme);
    }
    
    void apply_bc(const Mesh& msh, bool explicit_scheme){
        
        #ifdef HAVE_INTEL_TBB2
                size_t n_cells = m_elements_with_bc_eges.size();
                tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                    [this,&msh] (size_t & i){
                        size_t cell_ind = m_elements_with_bc_eges[i];
                        auto& cell = msh.backend_storage()->surfaces[cell_ind];
                        Matrix<T, Dynamic, Dynamic> laplacian_operator_loc = laplacian_operator(cell_ind, msh, cell);
                        scatter_bc_data(msh, cell, laplacian_operator_loc);
                }
            );
        #else
            auto storage = msh.backend_storage();
            for (auto& cell_ind : m_e_elements_with_bc_eges)
            {
                auto& cell = storage->surfaces[cell_ind];
                size_t e_cell_ind = m_e_cell_index[cell_ind];
                elastic_material_data<T> e_mat = m_e_material.find(cell_ind)->second;
                Matrix<T, Dynamic, Dynamic> mixed_operator_loc = e_mixed_operator(e_mat,msh,cell,explicit_scheme);
                scatter_e_bc_data(e_cell_ind, msh, cell, mixed_operator_loc);
            }
        
            for (auto& cell_ind : m_a_elements_with_bc_eges)
            {
                auto& cell = storage->surfaces[cell_ind];
                size_t a_cell_ind = m_a_cell_index[cell_ind];
                acoustic_material_data<T> a_mat = m_a_material.find(cell_ind)->second;
                Matrix<T, Dynamic, Dynamic> mixed_operator_loc = a_mixed_operator(a_mat, msh, cell, explicit_scheme);
                scatter_a_bc_data(a_cell_ind, msh, cell, mixed_operator_loc);
            }
        
        #endif
        
    }

    Matrix<T, Dynamic, Dynamic> e_mixed_operator(elastic_material_data<T> & material, const Mesh& msh, const typename Mesh::cell_type& cell, bool explicit_scheme){
            
        T rho = material.rho();
        T vp = material.vp();
        T vs = material.vs();
        T mu = rho * vs * vs;
        T lambda = rho * vp * vp - 2*mu;
    
        auto reconstruction_operator = strain_tensor_reconstruction(msh, cell);
        Matrix<T, Dynamic, Dynamic> R_operator = reconstruction_operator.second;
        auto n_rows = R_operator.rows();
        auto n_cols = R_operator.cols();
        Matrix<T, Dynamic, Dynamic> S_operator = Matrix<T, Dynamic, Dynamic>::Zero(n_rows, n_cols);

        if(explicit_scheme) {
            auto stabilization_operator = make_vector_hdg_stabilization(msh, cell, m_hho_di, m_scaled_stabilization_Q);
            auto n_s_rows = stabilization_operator.rows();
            auto n_s_cols = stabilization_operator.cols();
            S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
        }
        else {
            if(m_hho_stabilization_Q) {
                auto rec_for_stab = make_vector_hho_symmetric_laplacian(msh, cell, m_hho_di);
                auto stabilization_operator = make_vector_hho_stabilization(msh, cell, rec_for_stab.first, m_hho_di, m_scaled_stabilization_Q);
                auto n_s_rows = stabilization_operator.rows();
                auto n_s_cols = stabilization_operator.cols();
                S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
            }
            else {
                auto stabilization_operator = make_vector_hdg_stabilization(msh, cell, m_hho_di, m_scaled_stabilization_Q);
                auto n_s_rows = stabilization_operator.rows();
                auto n_s_cols = stabilization_operator.cols();
                S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
            }
        }
        
        // // COEF SPECTRAL RADIUS
        T equal_order_eta = 1.50;
        T mixed_order_eta = 0.95; 

        T coef = 1.0;
        T eta = coef * equal_order_eta;

        return R_operator + eta*(rho*vs)*S_operator;
    }
            
    Matrix<T, Dynamic, Dynamic>
    symmetric_tensor_mass_matrix(const Mesh& msh, const typename Mesh::cell_type& cell) {
            
        size_t dim = Mesh::dimension;
        auto gradeg = m_hho_di.grad_degree();
        auto ten_b = make_sym_matrix_monomial_basis(msh, cell, gradeg);
        auto ten_bs = disk::sym_matrix_basis_size(gradeg, dim, dim);
        Matrix<T, Dynamic, Dynamic> mass_matrix = Matrix<T, Dynamic, Dynamic>::Zero(ten_bs, ten_bs);
        
        auto qps = integrate(msh, cell, 2 * gradeg);

        // number of tensor components
        size_t dec = 0;
         if (dim == 3)
             dec = 6;
         else if (dim == 2)
             dec = 3;
         else
             std::logic_error("Expected 3 >= dim > 1");

         for (auto& qp : qps)
         {
             auto phi = ten_b.eval_functions(qp.point());

             for (size_t j = 0; j < ten_bs; j++)
             {
                 
                auto qp_phi_j = disk::priv::inner_product(qp.weight(), phi[j]);
                for (size_t i = j; i < ten_bs; i += dec){
                         mass_matrix(i, j) += disk::priv::inner_product(phi[i], qp_phi_j);
                }
             }
         }

        for (size_t j = 0; j < ten_bs; j++){
            for (size_t i = 0; i < j; i++){
                 mass_matrix(i, j) = mass_matrix(j, i);
            }
        }
        
        return mass_matrix;
    }
    
    Matrix<T, Dynamic, Dynamic>
    symmetric_tensor_trace_mass_matrix(const Mesh& msh, const typename Mesh::cell_type& cell) {
            
        size_t dim = Mesh::dimension;
        auto gradeg = m_hho_di.grad_degree();
        auto ten_b = make_sym_matrix_monomial_basis(msh, cell, gradeg);
        auto ten_bs = disk::sym_matrix_basis_size(gradeg, dim, dim);
        Matrix<T, Dynamic, Dynamic> mass_matrix = Matrix<T, Dynamic, Dynamic>::Zero(ten_bs, ten_bs);
        
        auto qps = integrate(msh, cell, 2 * gradeg);

        // number of tensor components
        size_t dec = 0;
         if (dim == 3)
             dec = 6;
         else if (dim == 2)
             dec = 3;
         else
             std::logic_error("Expected 3 >= dim > 1");

         for (auto& qp : qps) {

             auto phi = ten_b.eval_functions(qp.point());

             for (size_t j = 0; j < ten_bs; j++) {
                auto identity = phi[j];
                identity.setZero();
                for(size_t d = 0; d < dim; d++)
                    identity(d,d) = 1.0;

                auto trace = phi[j].trace();
                auto trace_phi_j = disk::priv::inner_product(phi[j].trace(), identity);
                auto qp_phi_j = disk::priv::inner_product(qp.weight(), trace_phi_j);
                for (size_t i = 0; i < ten_bs; i ++)
                    mass_matrix(i, j) += disk::priv::inner_product(phi[i], qp_phi_j);
                
             }
         }
        
        return mass_matrix;
    }
    
    Matrix<typename Mesh::coordinate_type, Dynamic, 1>
    e_mixed_rhs(const Mesh& msh, const typename Mesh::cell_type& cell, std::function<static_vector<double, 2>(const typename Mesh::point_type& )> & rhs_fun, size_t di = 0)
    {
        auto recdeg = m_hho_di.grad_degree();
        auto celdeg = m_hho_di.cell_degree();
        auto facdeg = m_hho_di.face_degree();

        auto ten_bs = disk::sym_matrix_basis_size(recdeg, Mesh::dimension, Mesh::dimension);
        auto vec_bs = disk::vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
        size_t n_cbs = ten_bs + vec_bs;
        auto cell_basis   = make_vector_monomial_basis(msh, cell, celdeg);
        using T = typename Mesh::coordinate_type;

        Matrix<T, Dynamic, 1> ret_loc = Matrix<T, Dynamic, 1>::Zero(cell_basis.size());
        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(n_cbs);

        const auto qps = integrate(msh, cell, 2 * (celdeg + di));

        for (auto& qp : qps)
        {
            const auto phi  = cell_basis.eval_functions(qp.point());
            const auto qp_f = disk::priv::inner_product(qp.weight(), rhs_fun(qp.point()));
            ret_loc += disk::priv::outer_product(phi, qp_f);
        }
        ret.block(ten_bs,0,vec_bs,1) = ret_loc;
        return ret;
    }
            
    Matrix<T, Dynamic, Dynamic> e_mass_operator(elastic_material_data<T> & material, const Mesh& msh, const typename Mesh::cell_type& cell, bool add_vector_mass_Q = true){
            
        size_t n_ten_cbs = disk::sym_matrix_basis_size(m_hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t n_cbs = n_ten_cbs + n_vec_cbs;
            
        T rho = material.rho();
        T vp = material.vp();
        T vs = material.vs();
        T mu = rho * vs * vs;
        T lambda = rho * vp * vp - 2*mu;
    
        Matrix<T, Dynamic, Dynamic> mass_matrix = Matrix<T, Dynamic, Dynamic>::Zero(n_cbs, n_cbs);
        
        // Symmetric stress tensor mass block
        
        // Stress tensor
        Matrix<T, Dynamic, Dynamic> mass_matrix_sigma  = symmetric_tensor_mass_matrix(msh, cell);
        
        // Tensor trace
        Matrix<T, Dynamic, Dynamic> mass_matrix_trace_sigma  = symmetric_tensor_trace_mass_matrix(msh, cell);
        
        // Constitutive relationship inverse
        mass_matrix_trace_sigma *= (lambda/(2.0*mu+2.0*lambda));
        mass_matrix_sigma -= mass_matrix_trace_sigma;
        mass_matrix_sigma *= (1.0/(2.0*mu));
        mass_matrix.block(0, 0, n_ten_cbs, n_ten_cbs) = mass_matrix_sigma;
        
        if (add_vector_mass_Q) {
            // vector velocity mass mass block
            auto vec_basis = disk::make_vector_monomial_basis(msh, cell, m_hho_di.cell_degree());
            Matrix<T, Dynamic, Dynamic> mass_matrix_v = disk::make_mass_matrix(msh, cell, vec_basis);
            mass_matrix_v *= rho;
            mass_matrix.block(n_ten_cbs, n_ten_cbs, n_vec_cbs, n_vec_cbs) = mass_matrix_v;
        }

        return mass_matrix;
    }
    
    Matrix<T, Dynamic, Dynamic> a_mass_operator(acoustic_material_data<T> & material, const Mesh& msh, const typename Mesh::cell_type& cell, bool add_scalar_mass_Q = true){
            
        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_cbs = n_scal_cbs + n_vec_cbs;
            
        Matrix<T, Dynamic, Dynamic> mass_matrix = Matrix<T, Dynamic, Dynamic>::Zero(n_cbs, n_cbs);
        
        auto recdeg = m_hho_di.reconstruction_degree();
        auto rec_basis = make_scalar_monomial_basis(msh, cell, recdeg);
        auto rbs = disk::scalar_basis_size(recdeg, Mesh::dimension);

        Matrix<T, Dynamic, Dynamic> mass_matrix_q_full  = make_stiffness_matrix(msh, cell, rec_basis);
        Matrix<T, Dynamic, Dynamic> mass_matrix_q = Matrix<T, Dynamic, Dynamic>::Zero(rbs-1, rbs-1);
        mass_matrix_q = mass_matrix_q_full.block(1, 1, rbs-1, rbs-1);
        mass_matrix_q *= material.rho();
        mass_matrix.block(0, 0, n_vec_cbs, n_vec_cbs) = mass_matrix_q;
        
        if (add_scalar_mass_Q) {
            auto scal_basis = disk::make_scalar_monomial_basis(msh, cell, m_hho_di.cell_degree());
            Matrix<T, Dynamic, Dynamic> mass_matrix_v = disk::make_mass_matrix(msh, cell, scal_basis);
            mass_matrix_v *= (1.0/(material.rho()*material.vp()*material.vp()));
            mass_matrix.block(n_vec_cbs, n_vec_cbs, n_scal_cbs, n_scal_cbs) = mass_matrix_v;
        }

        return mass_matrix;
    }
    
    std::pair<   Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
                 Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>  >
    strain_tensor_reconstruction(const Mesh& msh, const typename Mesh::cell_type& cell) {

        using T        = typename Mesh::coordinate_type;
        typedef Matrix<T, Dynamic, Dynamic> matrix_type;

        const size_t N = Mesh::dimension;

        auto graddeg = m_hho_di.grad_degree();
        auto celdeg  = m_hho_di.cell_degree();
        auto facdeg  = m_hho_di.face_degree();

        auto ten_b = make_sym_matrix_monomial_basis(msh, cell, graddeg);
        auto vec_b = make_vector_monomial_basis(msh, cell, celdeg);

        auto ten_bs = disk::sym_matrix_basis_size(graddeg, Mesh::dimension, Mesh::dimension);
        auto vec_bs = disk::vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
        auto fbs = disk::vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension);

        auto num_faces = howmany_faces(msh, cell);

        matrix_type gr_lhs = matrix_type::Zero(ten_bs, ten_bs);
        matrix_type gr_rhs = matrix_type::Zero(ten_bs, vec_bs + num_faces * fbs);
        
        const auto qps = integrate(msh, cell, 2 * graddeg);

        size_t dec = 0;
         if (N == 3)
             dec = 6;
         else if (N == 2)
             dec = 3;
         else
             std::logic_error("Expected 3 >= dim > 1");

         for (auto& qp : qps) {
             const auto gphi = ten_b.eval_functions(qp.point());
             for (size_t j = 0; j < ten_bs; j++) {
                auto qp_gphi_j = disk::priv::inner_product(qp.weight(), gphi[j]);
                for (size_t i = j; i < ten_bs; i += dec) {
                         gr_lhs(i, j) += disk::priv::inner_product(gphi[i], qp_gphi_j);
                }
             }
         }

         for (size_t j = 0; j < ten_bs; j++)
             for (size_t i = 0; i < j; i++)
                 gr_lhs(i, j) = gr_lhs(j, i);

         if (celdeg > 0) {
             const auto qpc = integrate(msh, cell, graddeg + celdeg - 1);
             for (auto& qp : qpc) {
                 const auto gphi    = ten_b.eval_functions(qp.point());
                 const auto dphi    = vec_b.eval_sgradients(qp.point());
                 const auto qp_dphi = disk::priv::inner_product(qp.weight(), dphi);
                 gr_rhs.block(0, 0, ten_bs, vec_bs) += disk::priv::outer_product(gphi, qp_dphi);
             }
         }

         const auto fcs = faces(msh, cell);
         for (size_t i = 0; i < fcs.size(); i++) {
             const auto fc = fcs[i];
             const auto n  = normal(msh, cell, fc);
             const auto fb = make_vector_monomial_basis(msh, fc, facdeg);
             const auto qps_f = integrate(msh, fc, graddeg + std::max(celdeg, facdeg));
             for (auto& qp : qps_f) {
                 const auto gphi = ten_b.eval_functions(qp.point());
                 const auto cphi = vec_b.eval_functions(qp.point());
                 const auto fphi = fb.eval_functions(qp.point());
                 const auto qp_gphi_n = disk::priv::inner_product(gphi, disk::priv::inner_product(qp.weight(), n));
                 gr_rhs.block(0, vec_bs + i * fbs, ten_bs, fbs) += disk::priv::outer_product(qp_gphi_n, fphi);
                 gr_rhs.block(0, 0, ten_bs, vec_bs) -= disk::priv::outer_product(qp_gphi_n, cphi);
             }
         }
            
        auto n_rows = gr_rhs.cols() + ten_bs;
        auto n_cols = gr_rhs.cols() + ten_bs;

        // Shrinking data
        matrix_type data_mixed = matrix_type::Zero(n_rows,n_cols);
        data_mixed.block(0, (ten_bs), ten_bs, n_cols-(ten_bs)) = -gr_rhs;
        data_mixed.block((ten_bs), 0, n_rows-(ten_bs), ten_bs) = gr_rhs.transpose();

        matrix_type oper = gr_lhs.llt().solve(gr_rhs);
        return std::make_pair(oper, data_mixed);
    }
    
    Matrix<T, Dynamic, Dynamic> a_mixed_operator(acoustic_material_data<T> & material, const Mesh& msh, const typename Mesh::cell_type& cell, bool explicit_scheme){
        
        auto reconstruction_operator = mixed_scalar_reconstruction(msh, cell);
        Matrix<T, Dynamic, Dynamic> R_operator = reconstruction_operator.second;
        auto n_rows = R_operator.rows();
        auto n_cols = R_operator.cols();

        Matrix<T, Dynamic, Dynamic> S_operator = Matrix<T, Dynamic, Dynamic>::Zero(n_rows, n_cols);
        if(explicit_scheme) {
            auto stabilization_operator = make_scalar_hdg_stabilization(msh, cell, m_hho_di, m_scaled_stabilization_Q);
            auto n_s_rows = stabilization_operator.rows();
            auto n_s_cols = stabilization_operator.cols();
            S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
        }
        else {
            if (m_hho_stabilization_Q) {
                auto stabilization_operator = make_scalar_hho_stabilization(msh, cell, reconstruction_operator.first, m_hho_di, m_scaled_stabilization_Q);
                auto n_s_rows = stabilization_operator.rows();
                auto n_s_cols = stabilization_operator.cols();
                S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
            } 
            else {
                auto stabilization_operator = make_scalar_hdg_stabilization(msh, cell, m_hho_di, m_scaled_stabilization_Q);
                auto n_s_rows = stabilization_operator.rows();
                auto n_s_cols = stabilization_operator.cols();
                S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
            }
        }

        // COEF SPECTRAL RADIUS
        T equal_order_eta = 0.8;
        T mixed_order_eta = 0.4; 
        

        T coef = 1.0;
        T eta = coef * equal_order_eta;

        T rho = material.rho();
        T vp  = material.vp();

        return R_operator + (eta/(vp*rho))*S_operator;
    }

    std::pair<   Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
                 Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>  >
    mixed_scalar_reconstruction(const Mesh& msh, const typename Mesh::cell_type& cell)
    {
        using T = typename Mesh::coordinate_type;
        typedef Matrix<T, Dynamic, Dynamic> matrix_type;
        typedef Matrix<T, Dynamic, 1>       vector_type;

        const size_t DIM = Mesh::dimension;

        const auto recdeg = m_hho_di.reconstruction_degree();
        const auto celdeg = m_hho_di.cell_degree();
        const auto facdeg = m_hho_di.face_degree();

        auto cb = make_scalar_monomial_basis(msh, cell, recdeg);

        const auto rbs = disk::scalar_basis_size(recdeg, Mesh::dimension);
        const auto cbs = disk::scalar_basis_size(celdeg, Mesh::dimension);
        const auto fbs = disk::scalar_basis_size(facdeg, Mesh::dimension - 1);

        const auto num_faces = howmany_faces(msh, cell);

        const matrix_type stiff  = make_stiffness_matrix(msh, cell, cb);
        matrix_type gr_lhs = matrix_type::Zero(rbs-1, rbs-1);
        matrix_type gr_rhs = matrix_type::Zero(rbs-1, cbs + num_faces*fbs);

        gr_lhs = stiff.block(1, 1, rbs-1, rbs-1);
        gr_rhs.block(0, 0, rbs-1, cbs) = stiff.block(1, 0, rbs-1, cbs);

        const auto fcs = faces(msh, cell);
        for (size_t i = 0; i < fcs.size(); i++)
        {
            const auto fc = fcs[i];
            const auto n  = normal(msh, cell, fc);
            auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

            auto qps_f = integrate(msh, fc, recdeg - 1 + std::max(facdeg,celdeg));
            for (auto& qp : qps_f)
            {
                vector_type c_phi_tmp = cb.eval_functions(qp.point());
                vector_type c_phi = c_phi_tmp.head(cbs);
                Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic, DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
                vector_type f_phi = fb.eval_functions(qp.point());
                gr_rhs.block(0, cbs+i*fbs, rbs-1, fbs) += qp.weight() * (c_dphi * n) * f_phi.transpose();
                gr_rhs.block(0, 0, rbs-1, cbs) -= qp.weight() * (c_dphi * n) * c_phi.transpose();
            }
        }

        auto vec_cell_size = gr_lhs.cols();
        auto nrows = gr_rhs.cols()+vec_cell_size;
        auto ncols = gr_rhs.cols()+vec_cell_size;
        
        // Shrinking data
        matrix_type data_mixed = matrix_type::Zero(nrows,ncols);
        data_mixed.block(0, vec_cell_size, vec_cell_size, ncols-vec_cell_size) = -gr_rhs;
        data_mixed.block(vec_cell_size, 0, nrows-vec_cell_size, vec_cell_size) = gr_rhs.transpose();
        
        matrix_type oper = gr_lhs.llt().solve(gr_rhs);
        return std::make_pair(oper, data_mixed);
    }
    
    Matrix<typename Mesh::coordinate_type, Dynamic, 1>
    a_mixed_rhs(const Mesh& msh, const typename Mesh::cell_type& cell, std::function<double(const typename Mesh::point_type& )> & rhs_fun, size_t di = 0)
    {
        const auto recdeg = m_hho_di.reconstruction_degree();
        const auto celdeg = m_hho_di.cell_degree();
        const auto rbs = disk::scalar_basis_size(recdeg, Mesh::dimension) - 1;
        const auto cbs = disk::scalar_basis_size(celdeg, Mesh::dimension) + rbs;
        auto cell_basis = make_scalar_monomial_basis(msh, cell, celdeg);
        using T = typename Mesh::coordinate_type;

        Matrix<T, Dynamic, 1> ret_loc = Matrix<T, Dynamic, 1>::Zero(cell_basis.size());
        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs);

        const auto qps = integrate(msh, cell, 2 * (celdeg + di));

        for (auto& qp : qps)
        {
            const auto phi  = cell_basis.eval_functions(qp.point());
            const auto qp_f = disk::priv::inner_product(qp.weight(), rhs_fun(qp.point()));
            ret_loc += disk::priv::outer_product(phi, qp_f);
        }
        ret.block(rbs,0,cell_basis.size(),1) = ret_loc;
        return ret;
    }
    
    Matrix<T, Dynamic, Dynamic> e_interface_operator(const Mesh& msh, const typename Mesh::face_type& face, const typename Mesh::cell_type& e_cell, const typename Mesh::cell_type& a_cell){

        Matrix<T, Dynamic, Dynamic> interface_operator;
        auto facdeg = m_hho_di.face_degree();
        auto vfbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        auto sfbs = disk::scalar_basis_size(facdeg, Mesh::dimension - 1);

        interface_operator = Matrix<T, Dynamic, Dynamic>::Zero(vfbs, sfbs);
        
        auto vfb = make_vector_monomial_basis(msh, face, facdeg);
        auto sfb = make_scalar_monomial_basis(msh, face, facdeg);
        const auto qps = integrate(msh, face, facdeg);
        const auto n = disk::normal(msh, e_cell, face);
        for (auto& qp : qps) {
          const auto v_f_phi = vfb.eval_functions(qp.point());
          const auto s_f_phi = sfb.eval_functions(qp.point());
          assert(v_f_phi.rows() == vfbs);
          assert(s_f_phi.rows() == sfbs);
          const auto n_dot_v_f_phi = disk::priv::inner_product(v_f_phi,disk::priv::inner_product(qp.weight(), n));
          const auto result = -1.0*disk::priv::outer_product(n_dot_v_f_phi, s_f_phi);
          interface_operator += result;
        }
        return interface_operator;
    }
    
    Matrix<T, Dynamic, 1> e_neumman_bc_operator(const Mesh& msh, const typename Mesh::face_type& face, const typename Mesh::cell_type& e_cell, const typename Mesh::cell_type& a_cell, std::function<T(const typename Mesh::point_type& )> a_vel_fun){

        Matrix<T, Dynamic, Dynamic> neumann_operator;
        auto facdeg = m_hho_di.face_degree();
        auto vfbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);

        neumann_operator = Matrix<T, Dynamic, Dynamic>::Zero(vfbs, 1);
        
        auto vfb = make_vector_monomial_basis(msh, face, facdeg);
        auto sfb = make_scalar_monomial_basis(msh, face, facdeg);
        const auto qps = integrate(msh, face, facdeg);
        const auto n = disk::normal(msh, e_cell, face);
        for (auto& qp : qps)
        {
            const auto v_f_phi = vfb.eval_functions(qp.point());
            const auto s_f_phi = sfb.eval_functions(qp.point());

            assert(v_f_phi.rows() == vfbs);
            const auto n_dot_v_f_phi = disk::priv::inner_product(v_f_phi,disk::priv::inner_product(qp.weight(), n));
            const auto result = disk::priv::outer_product(n_dot_v_f_phi, a_vel_fun(qp.point()));
            neumann_operator += result;
            
        }
        return neumann_operator;
        }
    
    Matrix<T, Dynamic, Dynamic> a_neumman_bc_operator(const Mesh& msh, const typename Mesh::face_type& face, const typename Mesh::cell_type& e_cell, const typename Mesh::cell_type& a_cell, std::function<static_vector<T, 2>(const typename Mesh::point_type& )> e_vel_fun){

        Matrix<T, Dynamic, Dynamic> neumann_operator;
        auto facdeg = m_hho_di.face_degree();
        auto sfbs = disk::scalar_basis_size(facdeg, Mesh::dimension - 1);

        neumann_operator = Matrix<T, Dynamic, Dynamic>::Zero(sfbs, 1);
        
        auto vfb = make_vector_monomial_basis(msh, face, facdeg);
        auto sfb = make_scalar_monomial_basis(msh, face, facdeg);
        const auto qps = integrate(msh, face, facdeg);
        const auto n = disk::normal(msh, a_cell, face);
        for (auto& qp : qps)
        {
            const auto v_f_phi = vfb.eval_functions(qp.point());
            const auto s_f_phi = sfb.eval_functions(qp.point());

            assert(s_f_phi.rows() == sfbs);
            const auto n_dot_v_f = disk::priv::inner_product(e_vel_fun(qp.point()),disk::priv::inner_product(qp.weight(), n));
            const auto result = disk::priv::inner_product(n_dot_v_f, s_f_phi);
            neumann_operator += result;
        }
        return neumann_operator;
    }
    
    void classify_cells(const Mesh& msh){

        m_e_elements_with_bc_eges.clear();
        for (auto& cell : msh)
        {
            auto cell_ind = msh.lookup(cell);
            auto face_list = faces(msh, cell);
            for (size_t face_i = 0; face_i < face_list.size(); face_i++)
            {
                auto fc = face_list[face_i];
                auto fc_id = msh.lookup(fc);
                bool is_dirichlet_Q = m_e_bnd.is_dirichlet_face(fc_id);
                if (is_dirichlet_Q)
                {
                    m_e_elements_with_bc_eges.push_back(cell_ind);
                    break;
                }
            }
        }
        
        m_a_elements_with_bc_eges.clear();
        for (auto& cell : msh)
        {
            typename Mesh::point_type bar = barycenter(msh, cell);
            if (bar.x() < 0) {
                continue;
            }
            
            auto cell_ind = msh.lookup(cell);
            auto face_list = faces(msh, cell);
            for (size_t face_i = 0; face_i < face_list.size(); face_i++)
            {
                auto fc = face_list[face_i];
                auto fc_id = msh.lookup(fc);
                bool is_dirichlet_Q = m_a_bnd.is_dirichlet_face(fc_id);
                if (is_dirichlet_Q)
                {
                    m_a_elements_with_bc_eges.push_back(cell_ind);
                    break;
                }
            }
        }
    }
    
    void build_cells_maps(){
        
        // elastic data
        size_t e_cell_ind = 0;
        for (auto chunk : m_e_material) {
            m_e_cell_index.insert(std::make_pair(chunk.first,e_cell_ind));
            e_cell_ind++;
        }
        m_n_elastic_cells = e_cell_ind;

        // acoustic data
        size_t a_cell_ind = 0;
        for (auto chunk : m_a_material) {
            m_a_cell_index.insert(std::make_pair(chunk.first,a_cell_ind));
            a_cell_ind++;
        }
        m_n_acoustic_cells = a_cell_ind;

    }
    
    void project_over_cells(const Mesh& msh, Matrix<T, Dynamic, 1> & x_glob, std::function<static_vector<T, 2>(const typename Mesh::point_type& )> v_fun, std::function<static_matrix<T,2,2>(const typename Mesh::point_type& )> flux_fun, std::function<T(const typename Mesh::point_type& )> v_s_fun, std::function<static_vector<T, 2>(const typename Mesh::point_type& )> flux_s_fun){
        
        auto storage = msh.backend_storage();
        size_t n_dof = MASS.rows();
        x_glob = Matrix<T, Dynamic, 1>::Zero(n_dof);

        // elastic block
    
        size_t n_ten_cbs = disk::sym_matrix_basis_size(m_hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t n_e_cbs = n_ten_cbs + n_vec_cbs;
    
        for (auto e_chunk : m_e_material) {
            size_t e_cell_ind = m_e_cell_index[e_chunk.first];
            auto& cell = storage->surfaces[e_chunk.first];
    
            Matrix<T, Dynamic, 1> x_proj_ten_dof = project_ten_function(msh, cell, flux_fun);
            Matrix<T, Dynamic, 1> x_proj_vec_dof = project_function(msh, cell, m_hho_di.cell_degree(), v_fun);
            
            Matrix<T, Dynamic, 1> x_proj_dof = Matrix<T, Dynamic, 1>::Zero(n_e_cbs);
            x_proj_dof.block(0, 0, n_ten_cbs, 1)          = x_proj_ten_dof;
            x_proj_dof.block(n_ten_cbs, 0, n_vec_cbs, 1)  = x_proj_vec_dof;
            scatter_e_cell_dof_data(e_cell_ind, msh, cell, x_glob, x_proj_dof);
        }
        
        // acoustic block
        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_v_s_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_a_cbs = n_v_s_cbs + n_scal_cbs;
        for (auto a_chunk : m_a_material) {
            size_t a_cell_ind = m_a_cell_index[a_chunk.first];
            auto& cell = storage->surfaces[a_chunk.first];
            
            Matrix<T, Dynamic, 1> x_proj_vec_dof = project_vec_function(msh, cell, flux_s_fun);
            Matrix<T, Dynamic, 1> x_proj_sca_dof = project_function(msh, cell, m_hho_di.cell_degree(), v_s_fun);
            
            Matrix<T, Dynamic, 1> x_proj_dof = Matrix<T, Dynamic, 1>::Zero(n_a_cbs);
            x_proj_dof.block(0, 0, n_v_s_cbs, 1) = x_proj_vec_dof;
            x_proj_dof.block(n_v_s_cbs, 0, n_scal_cbs, 1) = x_proj_sca_dof;
            
            scatter_a_cell_dof_data(a_cell_ind, msh, cell, x_glob, x_proj_dof);
        }
    
    }
    
    Matrix<T, Dynamic, 1> project_vec_function(const Mesh& msh, const typename Mesh::cell_type& cell,
                      std::function<static_vector<double, 2>(const typename Mesh::point_type& )> vec_fun){
    
            auto recdeg = m_hho_di.reconstruction_degree();
            auto rec_basis = make_scalar_monomial_basis(msh, cell, recdeg);
            auto rbs = disk::scalar_basis_size(recdeg, Mesh::dimension);
            Matrix<T, Dynamic, Dynamic> mass_matrix_q_full  = make_stiffness_matrix(msh, cell, rec_basis);
            Matrix<T, Dynamic, Dynamic> mass_matrix_q = Matrix<T, Dynamic, Dynamic>::Zero(rbs-1, rbs-1);
            mass_matrix_q = mass_matrix_q_full.block(1, 1, rbs-1, rbs-1);

            Matrix<T, Dynamic, 1> rhs = Matrix<T, Dynamic, 1>::Zero(rbs-1);
            Matrix<T, 1, 2> f_vec = Matrix<T, Dynamic, Dynamic>::Zero(1, 2);
            const auto qps = integrate(msh, cell, 2*recdeg);
            for (auto& qp : qps) {
                auto dphi = rec_basis.eval_gradients(qp.point());
                f_vec(0,0) = vec_fun(qp.point())[0];
                f_vec(0,1) = vec_fun(qp.point())[1];
                for (size_t i = 0; i < rbs-1; i++){
                    Matrix<T, 2, 1> phi_i = dphi.block(i+1, 0, 1, 2).transpose();
                    rhs(i,0) = rhs(i,0) + (qp.weight() * f_vec*phi_i)(0,0);
                }
            }
            Matrix<T, Dynamic, 1> x_dof = mass_matrix_q.llt().solve(rhs);
            return x_dof;
    }

    Matrix<T, Dynamic, 1> project_ten_function(const Mesh& msh, const typename Mesh::cell_type& cell,
                      std::function<static_matrix<T, 2,2>(const typename Mesh::point_type& )> ten_fun){
    
        Matrix<T, Dynamic, Dynamic> mass_matrix  = symmetric_tensor_mass_matrix(msh, cell);
        size_t dim = Mesh::dimension;
        auto gradeg = m_hho_di.grad_degree();
        auto ten_bs = disk::sym_matrix_basis_size(gradeg, dim, dim);
        auto ten_b = make_sym_matrix_monomial_basis(msh, cell, gradeg);
        Matrix<T, Dynamic, 1> rhs = Matrix<T, Dynamic, 1>::Zero(ten_bs);

        const auto qps = integrate(msh, cell, 2 * gradeg);
        for (auto& qp : qps)
        {
            auto phi = ten_b.eval_functions(qp.point());
            static_matrix<T, 2,2> sigma = ten_fun(qp.point());
            for (size_t i = 0; i < ten_bs; i++){
                auto qp_phi_i = disk::priv::inner_product(qp.weight(), phi[i]);
                rhs(i,0) += disk::priv::inner_product(qp_phi_i,sigma);
            }
        }
        Matrix<T, Dynamic, 1> x_dof = mass_matrix.llt().solve(rhs);
        return x_dof;
    }
    
    void
    scatter_e_cell_dof_data(size_t e_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cell,
                    Matrix<T, Dynamic, 1>& x_glob, Matrix<T, Dynamic, 1> x_proj_dof) const
    {
        size_t n_cbs = get_e_cell_basis_data();
        auto cell_ofs = e_cell_ind * n_cbs;
        x_glob.block(cell_ofs, 0, n_cbs, 1) = x_proj_dof;
    }
    
    void
    scatter_a_cell_dof_data(size_t a_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cell,
                    Matrix<T, Dynamic, 1>& x_glob, Matrix<T, Dynamic, 1> x_proj_dof) const
    {
        size_t n_cbs = get_a_cell_basis_data();
        auto cell_ofs = a_cell_ind * n_cbs + m_n_elastic_cell_dof;
        x_glob.block(cell_ofs, 0, n_cbs, 1) = x_proj_dof;
    }
    
    
    void project_over_faces(const Mesh& msh, Matrix<T, Dynamic, 1> & x_glob, std::function<static_vector<T, 2>(const typename Mesh::point_type& )> vec_fun, std::function<T(const typename Mesh::point_type& )> scal_fun){

        auto storage = msh.backend_storage();

        // elastic block
        for (auto e_chunk : m_e_material) {
            auto& cell = storage->surfaces[e_chunk.first];
            auto fcs = faces(msh, cell);
            for (size_t i = 0; i < fcs.size(); i++)
            {
                auto face = fcs[i];
                auto fc_id = msh.lookup(face);
                bool is_dirichlet_Q = m_e_bnd.is_dirichlet_face(fc_id);
                if (is_dirichlet_Q)
                {
                    continue;
                }
                Matrix<T, Dynamic, 1> x_proj_dof = project_function(msh, face, m_hho_di.face_degree(), vec_fun);
                scatter_e_face_dof_data(msh, face, x_glob, x_proj_dof);
            }
        }
        
        // acoustic block
        for (auto a_chunk : m_a_material) {
            auto& cell = storage->surfaces[a_chunk.first];
            auto fcs = faces(msh, cell);
            for (size_t i = 0; i < fcs.size(); i++)
            {
                auto face = fcs[i];
                auto fc_id = msh.lookup(face);
                bool is_dirichlet_Q = m_a_bnd.is_dirichlet_face(fc_id);
                if (is_dirichlet_Q)
                {
                    continue;
                }
                Matrix<T, Dynamic, 1> x_proj_dof = project_function(msh, face, m_hho_di.face_degree(), scal_fun);
                scatter_a_face_dof_data(msh, face, x_glob, x_proj_dof);
            }
        }
    }
    
    void
    scatter_e_face_dof_data(  const Mesh& msh, const typename Mesh::face_type& face,
                    Matrix<T, Dynamic, 1>& x_glob, Matrix<T, Dynamic, 1> x_proj_dof) const
    {
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        auto fc_id = msh.lookup(face);
        auto glob_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_e_compress_indexes.at(fc_id)*n_fbs;
        x_glob.block(glob_offset, 0, n_fbs, 1) = x_proj_dof;
    }
    
    void
    scatter_a_face_dof_data( const Mesh& msh, const typename Mesh::face_type& face,
                      Matrix<T, Dynamic, 1>& x_glob, Matrix<T, Dynamic, 1> x_proj_dof) const
    {
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        auto fc_id = msh.lookup(face);
        auto glob_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof + m_a_compress_indexes.at(fc_id)*n_fbs;
        x_glob.block(glob_offset, 0, n_fbs, 1) = x_proj_dof;
        
    }
    
    Matrix<T, Dynamic, 1>
    gather_e_dof_data(size_t cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, 1>& x_glob) const {
                        
        auto e_cell_ind = m_e_cell_index.find(cell_ind)->second;
        auto num_faces = howmany_faces(msh, cl);
        size_t n_cbs = get_e_cell_basis_data();
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        
        Matrix<T, Dynamic, 1> x_el(n_cbs + num_faces * n_fbs );
        x_el.block(0, 0, n_cbs, 1) = x_glob.block(e_cell_ind * n_cbs, 0, n_cbs, 1);
        auto fcs = faces(msh, cl);
        for (size_t i = 0; i < fcs.size(); i++) {
            auto fc = fcs[i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
            const auto face_id = eid.second;

            if (m_e_bnd.is_dirichlet_face( face_id)) {
                auto fb = make_vector_monomial_basis(msh, fc, m_hho_di.face_degree());
                auto dirichlet_fun = m_e_bnd.dirichlet_boundary_func(face_id);
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb);
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_fun);
                x_el.block(n_cbs + i * n_fbs, 0, n_fbs, 1) = mass.llt().solve(rhs);
            }
            else {
                auto fc_id = msh.lookup(fc);
                auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_e_compress_indexes.at(fc_id)*n_fbs;
                x_el.block(n_cbs + i*n_fbs, 0, n_fbs, 1) = x_glob.block(face_LHS_offset, 0, n_fbs, 1);
            }
        }
        return x_el;
    }
    
    Matrix<T, Dynamic, 1>
    gather_a_dof_data(size_t cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl,
                    const Matrix<T, Dynamic, 1>& x_glob) const
    {
        auto a_cell_ind = m_a_cell_index.find(cell_ind)->second;
        auto num_faces = howmany_faces(msh, cl);
        size_t n_cbs = get_a_cell_basis_data();
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        
        Matrix<T, Dynamic, 1> x_local(n_cbs + num_faces * n_fbs );
        x_local.block(0, 0, n_cbs, 1) = x_glob.block(a_cell_ind * n_cbs + m_n_elastic_cell_dof, 0, n_cbs, 1);
        auto fcs = faces(msh, cl);
        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto fc = fcs[i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
            const auto face_id                  = eid.second;

            if (m_a_bnd.is_dirichlet_face( face_id))
            {
                auto fb = disk::make_scalar_monomial_basis(msh, fc, m_hho_di.face_degree());
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb, m_hho_di.face_degree());
                auto velocity = m_a_bnd.dirichlet_boundary_func(face_id);
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, velocity, m_hho_di.face_degree());
                x_local.block(n_cbs + i * n_fbs, 0, n_fbs, 1) = mass.llt().solve(rhs);
            }
            else
            {
                auto fc_id = msh.lookup(fc);
                auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof + m_a_compress_indexes.at(fc_id)*n_fbs;
                x_local.block(n_cbs + i*n_fbs, 0, n_fbs, 1) = x_glob.block(face_LHS_offset, 0, n_fbs, 1);
            }
        }
        return x_local;
    }
            
    void finalize()
    {
        LHS.setFromTriplets( m_triplets.begin(), m_triplets.end() );
        m_triplets.clear();
        LHS_STAB.setFromTriplets( m_triplets_stab.begin(), m_triplets_stab.end() );
        m_triplets_stab.clear();
    }
            
    void finalize_mass()
    {
        MASS.setFromTriplets( m_mass_triplets.begin(), m_mass_triplets.end() );
        m_mass_triplets.clear();
    }

    void finalize_coupling()
    {
        COUPLING.setFromTriplets( m_c_triplets.begin(), m_c_triplets.end() );
        m_c_triplets.clear();
    }

    void set_coupling_stabilization() {
        m_hho_stabilization_Q = false;
    }   

    void set_hdg_stabilization() {
        // std::cout << bold << red << "   SUMMARY: " << reset << std::endl;
        if(m_hho_di.cell_degree() > m_hho_di.face_degree()) {
          m_hho_stabilization_Q = false;
      }                                                 
        else{
            m_hho_stabilization_Q = true;
        }
    }
    
    void set_interface_cell_indexes(std::map<size_t,std::pair<size_t,size_t>> & interface_cell_indexes){
        m_interface_cell_indexes = interface_cell_indexes;
    }
            
    void set_hho_stabilization(){
        m_hho_stabilization_Q = true;
    }
      
    void set_scaled_stabilization(){
        m_scaled_stabilization_Q = true;
    }
            
    e_boundary_type & get_e_bc_conditions(){
             return m_e_bnd;
    }
    
    a_boundary_type & get_a_bc_conditions(){
             return m_a_bnd;
    }
    
    std::map<size_t,elastic_material_data<T>> & get_e_material_data(){
             return m_e_material;
    }
    
    std::map<size_t,acoustic_material_data<T>> & get_a_material_data(){
             return m_a_material;
    }
                
    size_t get_a_n_cells_dof() const{
        return m_n_acoustic_cell_dof;
    }
    
    size_t get_e_n_cells_dof() const {
        return m_n_elastic_cell_dof;
    }
    
    size_t get_n_face_dof() const {
        size_t n_face_dof = m_n_elastic_face_dof + m_n_acoustic_face_dof;
        return n_face_dof;
    }

    size_t get_e_face_dof() const {
        return m_n_elastic_face_dof;
    }

    size_t get_a_face_dof() const {
        return m_n_acoustic_face_dof;
    }
    
    size_t get_e_cell_basis_data() const {
        size_t n_ten_cbs = disk::sym_matrix_basis_size(m_hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t n_cbs = n_ten_cbs + n_vec_cbs;
        return n_cbs;
    }
    
    size_t get_a_cell_basis_data() const {
        size_t n_vel_scal_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_cbs = n_vel_scal_cbs + n_scal_cbs;
        return n_cbs;
    }
            
    size_t get_cell_basis_data(){

        // Elastic cell basis
        size_t n_ten_cbs = disk::sym_matrix_basis_size(m_hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t n_cbs1    = n_ten_cbs + n_vec_cbs;
        
        // Acoustic cell basis 
        size_t n_vel_scal_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_scal_cbs     = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_cbs2         = n_vel_scal_cbs + n_scal_cbs;
        
        size_t n_cbs = n_cbs1 + n_cbs2;
        
        return n_cbs;
    }

    size_t get_n_faces() {
        return m_n_edges - m_n_essential_edges;
    }

    size_t get_face_basis_data() {
      size_t n_fbs1 = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
      size_t n_fbs2 = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
      size_t n_fbs = n_fbs1 + n_fbs2;
      return n_fbs;
    }

    size_t get_e_face_basis_data() {
      size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
      return n_fbs;
    }

    size_t get_a_face_basis_data() {
      size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
      return n_fbs;
    }

    size_t get_acoustic_cells() {
        return m_n_acoustic_cells;
    }

    size_t get_elastic_cells() {
        return m_n_elastic_cells;
    }

    size_t get_elastic_faces() {
        return n_e_edges;
    }

    size_t get_acoustic_faces() {
        return n_a_edges;
    }

    std::map<size_t,std::pair<size_t,size_t>>  get_interfaces() {
        return m_interface_cell_indexes;
    }

    std::vector<size_t> get_e_compress() {
        return m_e_compress_indexes;
    }

    std::vector<size_t> get_a_compress() {
        return m_a_compress_indexes;
    }

    std::vector<size_t> get_e_expand() {
        return m_e_expand_indexes;
    }

    std::vector<size_t> get_a_expand() {
        return m_a_expand_indexes;
    }

   std::pair<size_t,size_t> Scc_block_dimension(){
        
        // Elastic block
        auto graddeg = m_hho_di.grad_degree();
        auto ten_bs = disk::sym_matrix_basis_size(graddeg, Mesh::dimension, Mesh::dimension)*m_n_elastic_cells;
        
        // Acoustic block
        auto recdeg = m_hho_di.reconstruction_degree();
        auto rbs = disk::scalar_basis_size(recdeg, Mesh::dimension);
        auto vec_cell_size = (rbs-1)*m_n_acoustic_cells; 

        return std::make_pair(ten_bs, vec_cell_size);
        
    }

};

#endif /* elastoacoustic_four_fields_assembler_hpp */
