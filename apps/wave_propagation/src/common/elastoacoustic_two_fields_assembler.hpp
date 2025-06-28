//
//  elastoacoustic_two_fields_assembler.hpp
//  acoustics
//
//  Created by Omar Durán on 6/3/20.
//

#pragma once
#ifndef elastoacoustic_two_fields_assembler_hpp
#define elastoacoustic_two_fields_assembler_hpp

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
class elastoacoustic_two_fields_assembler
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
    std::vector< Triplet<T> >           m_c_triplets;
    std::vector< Triplet<T> >           m_mass_triplets;
    std::map<size_t,elastic_material_data<T>> m_e_material;
    std::map<size_t,acoustic_material_data<T>> m_a_material;
    std::map<size_t,size_t> m_e_cell_index;
    std::map<size_t,size_t> m_a_cell_index;
    std::vector< size_t >               m_e_elements_with_bc_eges;
    std::vector< size_t >               m_a_elements_with_bc_eges;
    std::map<size_t,std::pair<size_t,size_t>>   m_interface_cell_indexes;

public:

    size_t      m_n_edges;
    size_t      m_n_essential_edges;
    bool        m_hho_stabilization_Q;
    size_t      m_n_elastic_cell_dof;
    size_t      m_n_acoustic_cell_dof;
    size_t      m_n_elastic_face_dof;
    size_t      m_n_acoustic_face_dof;

    SparseMatrix<T>         LHS;
    Matrix<T, Dynamic, 1>   RHS;
    SparseMatrix<T>         MASS;
    SparseMatrix<T>         COUPLING;

    elastoacoustic_two_fields_assembler(const Mesh& msh, const disk::hho_degree_info& hho_di, const e_boundary_type& e_bnd, const a_boundary_type& a_bnd, std::map<size_t,elastic_material_data<T>> & e_material, std::map<size_t,acoustic_material_data<T>> & a_material)
        : m_hho_di(hho_di), m_e_bnd(e_bnd), m_a_bnd(a_bnd), m_e_material(e_material), m_a_material(a_material), m_hho_stabilization_Q(true)
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
        size_t n_e_edges = e_egdes.size();

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
        size_t n_a_edges = a_egdes.size();
        
        m_n_edges = msh.faces_size();
        m_n_essential_edges = n_e_essential_edges + n_a_essential_edges;

        m_e_compress_indexes.resize( m_n_edges );
        m_e_expand_indexes.resize( m_n_edges - m_n_essential_edges );
        
        m_a_compress_indexes.resize( m_n_edges );
        m_a_expand_indexes.resize( m_n_edges - m_n_essential_edges );

        
        
        size_t e_compressed_offset = 0;
        for (auto face_id : e_egdes)
        {
            m_e_compress_indexes.at(face_id) = e_compressed_offset;
            m_e_expand_indexes.at(e_compressed_offset) = face_id;
            e_compressed_offset++;
        }
        
        size_t a_compressed_offset = 0;
        for (auto face_id : a_egdes)
        {
            m_a_compress_indexes.at(face_id) = a_compressed_offset;
            m_a_expand_indexes.at(a_compressed_offset) = face_id;
            a_compressed_offset++;
        }

        size_t n_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        
        size_t n_s_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_s_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);

        m_n_elastic_cell_dof = (n_cbs * m_e_material.size());
        m_n_acoustic_cell_dof = (n_s_cbs * m_a_material.size());
        
        m_n_elastic_face_dof = (n_fbs * n_e_edges);
        m_n_acoustic_face_dof = (n_s_fbs * n_a_edges);
        size_t system_size = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof + m_n_acoustic_face_dof;

        LHS = SparseMatrix<T>( system_size, system_size );
        RHS = Matrix<T, Dynamic, 1>::Zero( system_size );
        MASS = SparseMatrix<T>( system_size, system_size );
        COUPLING = SparseMatrix<T>( system_size, system_size );
            
        classify_cells(msh);
        build_cells_maps();
    }
    
    void scatter_e_data(size_t e_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs,
             const Matrix<T, Dynamic, 1>& rhs)
    {
        auto fcs = faces(msh, cl);
        size_t n_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs + n_fbs*fcs.size());

        auto cell_LHS_offset    = e_cell_ind * n_cbs;

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
            RHS(int(asm_map[i])) += rhs(i);
        }

    }
    
    void scatter_a_data(size_t a_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl,
    const Matrix<T, Dynamic, Dynamic>& lhs,
    const Matrix<T, Dynamic, 1>& rhs)
    {
        auto fcs = faces(msh, cl);
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
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
            RHS(int(asm_map[i])) += rhs(i);
        }

    }
    
    void scatter_e_mass_data(size_t e_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& mass_matrix)
    {
        size_t n_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
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
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
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
    
    void scatter_ea_interface_data(const Mesh& msh, const typename Mesh::face_type& face,
             const Matrix<T, Dynamic, Dynamic>& interface_matrix)
    {
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

        for (size_t i = 0; i < interface_matrix.rows(); i++)
        {
            for (size_t j = 0; j < interface_matrix.cols(); j++)
            {
                    m_c_triplets.push_back( Triplet<T>(asm_map_e[i], asm_map_a[j], interface_matrix(i,j)) );
                    m_c_triplets.push_back( Triplet<T>(asm_map_a[j], asm_map_e[i], - interface_matrix(i,j)) );
            }
        }

    }
    
    void scatter_e_interface_neumann_data(const Mesh& msh, const typename Mesh::face_type& face,
             const Matrix<T, Dynamic, 1>& interface_bc_data)
    {
        auto vfbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);

        std::vector<assembly_index> asm_map_e;
        asm_map_e.reserve(vfbs);
        
        auto fc_id = msh.lookup(face);
        
        auto e_face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_e_compress_indexes.at(fc_id)*vfbs;
        bool e_dirichlet = m_e_bnd.is_dirichlet_face(fc_id);
        for (size_t i = 0; i < vfbs; i++){
            asm_map_e.push_back( assembly_index(e_face_LHS_offset+i, !e_dirichlet) );
        }
                    
        assert( asm_map_e.size() == interface_bc_data.rows());

        for (size_t i = 0; i < interface_bc_data.rows(); i++)
        {
            if (!asm_map_e[i].assemble())
                 continue;
             RHS(int(asm_map_e[i])) += interface_bc_data(i,0);
        }

    }
    
    
    void scatter_a_interface_neumann_data(const Mesh& msh, const typename Mesh::face_type& face,
             const Matrix<T, Dynamic, 1>& interface_bc_data)
    {

        auto sfbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);

        std::vector<assembly_index> asm_map_a;
        asm_map_a.reserve(sfbs);
        
        auto fc_id = msh.lookup(face);
        
        auto a_face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof + m_a_compress_indexes.at(fc_id)*sfbs;
        bool a_dirichlet = m_a_bnd.is_dirichlet_face(fc_id);
        for (size_t i = 0; i < sfbs; i++){
            asm_map_a.push_back( assembly_index(a_face_LHS_offset+i, !a_dirichlet) );
        }
        
        assert( asm_map_a.size() == interface_bc_data.rows() );

        for (size_t i = 0; i < interface_bc_data.rows(); i++)
        {
            if (!asm_map_a[i].assemble())
                 continue;
             RHS(int(asm_map_a[i])) += interface_bc_data(i,0);
        }

    }
    
    void scatter_e_rhs_data(size_t e_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, 1>& rhs)
    {
    
        size_t n_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs);

        auto cell_LHS_offset    = e_cell_ind * n_cbs;

        for (size_t i = 0; i < n_cbs; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );

        assert( asm_map.size() == rhs.rows() );

        for (size_t i = 0; i < rhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;
            RHS(int(asm_map[i])) += rhs(i);
        }

    }
    
    void scatter_a_rhs_data(size_t a_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, 1>& rhs)
    {
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
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
            RHS(int(asm_map[i])) += rhs(i);
        }

    }
    
    void scatter_e_bc_data(size_t e_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs)
    {
    
        auto fcs = faces(msh, cl);
        size_t n_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
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
                    RHS(int(asm_map[i])) -= lhs(i,j) * dirichlet_data(j);
            }
        }

    }
    
    void scatter_a_bc_data(size_t a_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs)
    {
        auto fcs = faces(msh, cl);
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
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
                    RHS(int(asm_map[i])) -= lhs(i,j) * dirichlet_data(j);
                    
            }
        }

    }
    
    void assemble(const Mesh& msh, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> e_rhs_fun, std::function<T(const typename Mesh::point_type& )> a_rhs_fun){
        
        auto storage = msh.backend_storage();
        LHS.setZero();
        RHS.setZero();
        
        // elastic block
        for (auto e_chunk : m_e_material) {
            size_t e_cell_ind = m_e_cell_index[e_chunk.first];
            auto& cell = storage->surfaces[e_chunk.first];
            Matrix<T, Dynamic, Dynamic> vectorial_laplacian_operator_loc = e_laplacian_operator(e_chunk.second,msh,cell);
            auto cell_basis   = make_vector_monomial_basis(msh, cell, m_hho_di.cell_degree());
            Matrix<T, Dynamic, 1> f_loc = make_rhs(msh, cell, cell_basis, e_rhs_fun);
            scatter_e_data(e_cell_ind, msh, cell, vectorial_laplacian_operator_loc, f_loc);
        }
        
        // acoustic block
        for (auto a_chunk : m_a_material) {
            size_t a_cell_ind = m_a_cell_index[a_chunk.first];
            auto& cell = storage->surfaces[a_chunk.first];
            
            Matrix<T, Dynamic, Dynamic> laplacian_operator_loc = a_laplacian_operator(a_chunk.second, msh, cell);
            auto cell_basis   = make_scalar_monomial_basis(msh, cell, m_hho_di.cell_degree());
            Matrix<T, Dynamic, 1> f_loc = make_rhs(msh, cell, cell_basis, a_rhs_fun);
            scatter_a_data(a_cell_ind, msh, cell, laplacian_operator_loc, f_loc);
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
    
    void apply_bc_conditions_on_interface(const Mesh& msh, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> e_vel_fun, std::function<T(const typename Mesh::point_type& )> a_vel_fun){
        
        auto storage = msh.backend_storage();
        // Applying transmission conditions as neumann data for the case of uncoupled problems
        for (auto chunk : m_interface_cell_indexes) {
            auto& face = storage->edges[chunk.first];
            auto& e_cell = storage->surfaces[chunk.second.first];
            auto& a_cell = storage->surfaces[chunk.second.second];
            
            Matrix<T, Dynamic, 1> e_neumann_operator_loc = e_neumman_bc_operator(msh, face, e_cell, a_cell, a_vel_fun);
            scatter_e_interface_neumann_data(msh, face, e_neumann_operator_loc);
            
            Matrix<T, Dynamic, 1> a_neumann_operator_loc = a_neumman_bc_operator(msh, face, e_cell, a_cell, e_vel_fun);
            scatter_a_interface_neumann_data(msh, face, a_neumann_operator_loc);
            
        }
    }

    void assemble_mass(const Mesh& msh){
        
        auto storage = msh.backend_storage();
        MASS.setZero();
        
        // elastic block
        for (auto e_chunk : m_e_material) {
            size_t e_cell_ind = m_e_cell_index[e_chunk.first];
            auto& cell = storage->surfaces[e_chunk.first];
            Matrix<T, Dynamic, Dynamic> mass_matrix = e_mass_operator(e_chunk.second,msh, cell);
            scatter_e_mass_data(e_cell_ind,msh, cell, mass_matrix);
        }
        
        // acoustic block
        for (auto a_chunk : m_a_material) {
            size_t a_cell_ind = m_a_cell_index[a_chunk.first];
            auto& cell = storage->surfaces[a_chunk.first];
            Matrix<T, Dynamic, Dynamic> mass_matrix = a_mass_operator(a_chunk.second,msh, cell);
            scatter_a_mass_data(a_cell_ind,msh, cell, mass_matrix);
        }
        
        finalize_mass();
    }
    
    void assemble_rhs(const Mesh& msh, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> e_rhs_fun, std::function<T(const typename Mesh::point_type& )> a_rhs_fun){
        
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
                 auto cell_basis   = make_vector_monomial_basis(msh, cell, m_hho_di.cell_degree());
                 Matrix<T, Dynamic, 1> f_loc = make_rhs(msh, cell, cell_basis, e_rhs_fun);
                 scatter_e_rhs_data(e_cell_ind, msh, cell, f_loc);
             }
        
            for (auto a_chunk : m_a_material) {
                size_t a_cell_ind = m_a_cell_index[a_chunk.first];
                auto& cell = storage->surfaces[a_chunk.first];
                auto cell_basis   = make_scalar_monomial_basis(msh, cell, m_hho_di.cell_degree());
                Matrix<T, Dynamic, 1> f_loc = make_rhs(msh, cell, cell_basis, a_rhs_fun);
                scatter_a_rhs_data(a_cell_ind, msh, cell, f_loc);
            }
        
        #endif
        apply_bc(msh);
    }
    
    void apply_bc(const Mesh& msh){
        
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
                elastic_material_data<T> e_mat = m_e_material.find(cell_ind)->second;
                Matrix<T, Dynamic, Dynamic> vectorial_laplacian_operator_loc = e_laplacian_operator(e_mat,msh,cell);
                scatter_e_bc_data(cell_ind, msh, cell, vectorial_laplacian_operator_loc);
            }
        
            for (auto& cell_ind : m_a_elements_with_bc_eges)
            {
                auto& cell = storage->surfaces[cell_ind];
                acoustic_material_data<T> a_mat = m_a_material.find(cell_ind)->second;
                Matrix<T, Dynamic, Dynamic> laplacian_operator_loc = a_laplacian_operator(a_mat, msh, cell);
                scatter_a_bc_data(cell_ind, msh, cell, laplacian_operator_loc);
            }
        
        #endif
        
    }
    
    Matrix<T, Dynamic, Dynamic> e_mass_operator(elastic_material_data<T> & material, const Mesh& msh, const typename Mesh::cell_type& cell){
        auto vec_basis = disk::make_vector_monomial_basis(msh, cell, m_hho_di.cell_degree());
        Matrix<T, Dynamic, Dynamic> mass_matrix = disk::make_mass_matrix(msh, cell, vec_basis);
        mass_matrix *= (material.rho());
        return mass_matrix;
    }
    
    Matrix<T, Dynamic, Dynamic> a_mass_operator(acoustic_material_data<T> & material, const Mesh& msh, const typename Mesh::cell_type& cell){

        auto scal_basis = disk::make_scalar_monomial_basis(msh, cell, m_hho_di.cell_degree());
        Matrix<T, Dynamic, Dynamic> mass_matrix = disk::make_mass_matrix(msh, cell, scal_basis);
        mass_matrix *= (1.0/(material.rho()*material.vp()*material.vp()));
        return mass_matrix;
    }
    
    
    Matrix<T, Dynamic, Dynamic> e_laplacian_operator(elastic_material_data<T> & material, const Mesh& msh, const typename Mesh::cell_type& cell){
           
        T mu = material.rho()*material.vs()*material.vs();
        T lambda = material.rho()*material.vp()*material.vp() - 2.0*mu;
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
        return 2.0 * mu * (R_operator + S_operator) + lambda * divergence_operator.second;
    }
    
    Matrix<T, Dynamic, Dynamic> a_laplacian_operator(acoustic_material_data<T> & material, const Mesh& msh, const typename Mesh::cell_type& cell){

        auto reconstruction_operator   = make_scalar_hho_laplacian(msh, cell, m_hho_di);
        Matrix<T, Dynamic, Dynamic> R_operator = reconstruction_operator.second;
        
        Matrix<T, Dynamic, Dynamic> S_operator;
        if(m_hho_stabilization_Q)
        {
            auto stabilization_operator    = make_scalar_hho_stabilization(msh, cell, reconstruction_operator.first, m_hho_di);
            S_operator = stabilization_operator;
        }else{
            auto stabilization_operator    = make_scalar_hdg_stabilization(msh, cell, m_hho_di);
            S_operator = stabilization_operator;
        }
        return (1.0/material.rho())*(R_operator + S_operator);
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
        for (auto& qp : qps)
        {
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
    
    Matrix<T, Dynamic, Dynamic> a_neumman_bc_operator(const Mesh& msh, const typename Mesh::face_type& face, const typename Mesh::cell_type& e_cell, const typename Mesh::cell_type& a_cell, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> e_vel_fun){

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
        
        // acoustic data
        size_t a_cell_ind = 0;
        for (auto chunk : m_a_material) {
            m_a_cell_index.insert(std::make_pair(chunk.first,a_cell_ind));
            a_cell_ind++;
        }
    }
    
    void project_over_cells(const Mesh& msh, Matrix<T, Dynamic, 1> & x_glob, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> vec_fun, std::function<T(const typename Mesh::point_type& )> scal_fun){
        
        auto storage = msh.backend_storage();
        size_t n_dof = MASS.rows();
        x_glob = Matrix<T, Dynamic, 1>::Zero(n_dof);

        // elastic block
        for (auto e_chunk : m_e_material) {
            size_t e_cell_ind = m_e_cell_index[e_chunk.first];
            auto& cell = storage->surfaces[e_chunk.first];
            Matrix<T, Dynamic, 1> x_proj_dof = project_function(msh, cell, m_hho_di.cell_degree(), vec_fun);
            scatter_e_cell_dof_data(e_cell_ind, msh, cell, x_glob, x_proj_dof);
        }
        
        // acoustic block
        for (auto a_chunk : m_a_material) {
            size_t a_cell_ind = m_a_cell_index[a_chunk.first];
            auto& cell = storage->surfaces[a_chunk.first];
            Matrix<T, Dynamic, 1> x_proj_dof = project_function(msh, cell, m_hho_di.cell_degree(), scal_fun);
            scatter_a_cell_dof_data(a_cell_ind, msh, cell, x_glob, x_proj_dof);
        }
    
    }
    
    void
    scatter_e_cell_dof_data(size_t e_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cell,
                    Matrix<T, Dynamic, 1>& x_glob, Matrix<T, Dynamic, 1> x_proj_dof) const
    {
        size_t n_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        auto cell_ofs = e_cell_ind * n_cbs;
        x_glob.block(cell_ofs, 0, n_cbs, 1) = x_proj_dof;
    }
    
    void
    scatter_a_cell_dof_data(size_t a_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cell,
                    Matrix<T, Dynamic, 1>& x_glob, Matrix<T, Dynamic, 1> x_proj_dof) const
    {
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        auto cell_ofs = a_cell_ind * n_cbs + m_n_elastic_cell_dof;
        x_glob.block(cell_ofs, 0, n_cbs, 1) = x_proj_dof;
    }
    
    
    
    void project_over_faces(const Mesh& msh, Matrix<T, Dynamic, 1> & x_glob, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> vec_fun, std::function<T(const typename Mesh::point_type& )> scal_fun){

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
    gather_e_dof_data(size_t cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl,
                    const Matrix<T, Dynamic, 1>& x_glob) const
    {
        auto e_cell_ind = m_e_cell_index.find(cell_ind)->second;
        auto num_faces = howmany_faces(msh, cl);
        size_t n_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        
        Matrix<T, Dynamic, 1> x_el(n_cbs + num_faces * n_fbs );
        x_el.block(0, 0, n_cbs, 1) = x_glob.block(e_cell_ind * n_cbs, 0, n_cbs, 1);
        auto fcs = faces(msh, cl);
        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto fc = fcs[i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
            const auto face_id                  = eid.second;

            if (m_e_bnd.is_dirichlet_face( face_id))
            {
                auto fb = make_vector_monomial_basis(msh, fc, m_hho_di.face_degree());
                auto dirichlet_fun  = m_e_bnd.dirichlet_boundary_func(face_id);
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb);
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_fun);
                x_el.block(n_cbs + i * n_fbs, 0, n_fbs, 1) = mass.llt().solve(rhs);
            }
            else
            {

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
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
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
    
    void set_interface_cell_indexes(std::map<size_t,std::pair<size_t,size_t>> & interface_cell_indexes){
        m_interface_cell_indexes = interface_cell_indexes;
    }
            
    void set_hho_stabilization(){
        m_hho_stabilization_Q = true;
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
                
    size_t get_a_n_cells_dof(){
        return m_n_acoustic_cell_dof;
    }
    
    size_t get_e_n_cells_dof(){
        return m_n_elastic_cell_dof;
    }
    
    size_t get_n_face_dof(){
        size_t n_face_dof = m_n_elastic_face_dof + m_n_acoustic_face_dof;
        return n_face_dof;
    }
    
    size_t get_e_cell_basis_data(){
        size_t n_v_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        return n_v_cbs;
    }
    
    size_t get_a_cell_basis_data(){
        size_t n_s_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        return n_s_cbs;
    }
            
};

#endif /* elastoacoustic_two_fields_assembler_hpp */
