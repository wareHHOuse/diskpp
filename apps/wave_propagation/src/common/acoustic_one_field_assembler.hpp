//
//  acoustic_one_field_assembler.hpp
//  acoustics
//
//  Created by Omar Durán on 4/14/20.
//

#pragma once
#ifndef acoustic_one_field_assembler_hpp
#define acoustic_one_field_assembler_hpp

#include "diskpp/bases/bases.hpp"
#include "diskpp/methods/hho"
#include "../common/assembly_index.hpp"
#include "../common/acoustic_material_data.hpp"

#ifdef HAVE_INTEL_TBB
#include <tbb/parallel_for.h>
#endif

template<typename Mesh>
class acoustic_one_field_assembler
{
    
    
    typedef disk::BoundaryConditions<Mesh, true>    boundary_type;
    using T = typename Mesh::coordinate_type;

    std::vector<size_t>                 m_compress_indexes;
    std::vector<size_t>                 m_expand_indexes;

    disk::hho_degree_info               m_hho_di;
    boundary_type                       m_bnd;
    std::vector< Triplet<T> >           m_triplets;
    std::vector< Triplet<T> >           m_mass_triplets;
    std::vector< acoustic_material_data<T> > m_material;
    std::vector< size_t >               m_elements_with_bc_eges;
    
    std::vector< std::pair<size_t,size_t> >  m_f_l_indexes;

    size_t      m_n_edges;
    size_t      m_n_essential_edges;
    bool        m_hho_stabilization_Q;

public:

    SparseMatrix<T>         LHS;
    Matrix<T, Dynamic, 1>   RHS;
    SparseMatrix<T>         MASS;

    acoustic_one_field_assembler(const Mesh& msh, const disk::hho_degree_info& hho_di, const boundary_type& bnd)
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

        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);

        size_t system_size = n_cbs * msh.cells_size() + n_fbs * (m_n_edges - m_n_essential_edges);

        LHS = SparseMatrix<T>( system_size, system_size );
        RHS = Matrix<T, Dynamic, 1>::Zero( system_size );
        MASS = SparseMatrix<T>( system_size, system_size );
            
        classify_cells(msh);
        expand_triplet(msh);
    }
    
    void scatter_data(const std::pair<size_t,size_t> & dest_indexes, const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs,
             const Matrix<T, Dynamic, 1>& rhs)
    {
        auto fcs = faces(msh, cl);
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs + n_fbs*fcs.size());

        auto cell_offset        = disk::priv::offset(msh, cl);
        auto cell_LHS_offset    = cell_offset * n_cbs;

        for (size_t i = 0; i < n_cbs; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
        
        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto face_offset = disk::priv::offset(msh, fc);
            auto face_LHS_offset = n_cbs * msh.cells_size() + m_compress_indexes.at(face_offset)*n_fbs;

            auto fc_id = msh.lookup(fc);
            bool dirichlet = m_bnd.is_dirichlet_face(fc_id);

            for (size_t i = 0; i < n_fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
            
        }

        assert( asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );

        size_t first_ind = dest_indexes.first;
        size_t last_ind  = dest_indexes.second;
        size_t l = 0;
        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if ( asm_map[j].assemble() ){
                    m_triplets[first_ind+l] = Triplet<T>(asm_map[i], asm_map[j], lhs(i,j));
                    l++;
                }
            }
        }
        
        size_t cell_nnz = (last_ind - first_ind);
        assert( (l-1) == cell_nnz );

        for (size_t i = 0; i < rhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;
            RHS(int(asm_map[i])) += rhs(i);
        }

    }
            
    void scatter_bc_data(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs)
    {
        auto fcs = faces(msh, cl);
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
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
                 auto fb = make_scalar_monomial_basis(msh, fc, m_hho_di.face_degree());
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
                if ( !asm_map[j].assemble() )
		  RHS(int(asm_map[i])) -= lhs(i,j) * dirichlet_data(j);
                    
            }
        }

    }
            
    void scatter_rhs_data(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, 1>& rhs)
    {
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs);

        auto cell_offset        = disk::priv::offset(msh, cl);
        auto cell_LHS_offset    = cell_offset * n_cbs;

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
            
    void scatter_mass_data(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& mass_matrix)
    {
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs);

        auto cell_offset        = disk::priv::offset(msh, cl);
        auto cell_LHS_offset    = cell_offset * n_cbs;

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

    void assemble(const Mesh& msh, std::function<T(const typename Mesh::point_type& )> rhs_fun){
        
        LHS.setZero();
        RHS.setZero();
        #ifdef HAVE_INTEL_TBB
                size_t n_cells = msh.cells_size();
                tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                    [this,&msh,&rhs_fun] (size_t & cell_ind){
                        auto& cell = msh.backend_storage()->surfaces[cell_ind];
                        Matrix<T, Dynamic, Dynamic> laplacian_operator_loc = this->laplacian_operator(cell_ind, msh, cell);
                        auto cell_basis   = make_scalar_monomial_basis(msh, cell, m_hho_di.cell_degree());
                        Matrix<T, Dynamic, 1> f_loc = make_rhs(msh, cell, cell_basis, rhs_fun);
                        this->scatter_data(m_f_l_indexes[cell_ind], msh, cell, laplacian_operator_loc, f_loc);
                }
            );
        #else
            for (size_t cell_ind = 0; cell_ind < msh.cells_size(); cell_ind++)
            {
                auto& cell = msh.backend_storage()->surfaces[cell_ind];
                Matrix<T, Dynamic, Dynamic> laplacian_operator_loc = laplacian_operator(cell_ind, msh, cell);
                auto cell_basis   = make_scalar_monomial_basis(msh, cell, m_hho_di.cell_degree());
                Matrix<T, Dynamic, 1> f_loc = make_rhs(msh, cell, cell_basis, rhs_fun);
                
                scatter_data(m_f_l_indexes[cell_ind], msh, cell, laplacian_operator_loc, f_loc);
            }
        #endif

        finalize();
    }
            
    void apply_bc(const Mesh& msh, T scale = T(1.0)){
        
        #ifdef HAVE_INTEL_TBB
                size_t n_cells = m_elements_with_bc_eges.size();
                tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                    [this,&msh,&scale] (size_t & i){
                        size_t cell_ind = m_elements_with_bc_eges[i];
                        auto& cell = msh.backend_storage()->surfaces[cell_ind];
                        Matrix<T, Dynamic, Dynamic> laplacian_operator_loc = laplacian_operator(cell_ind, msh, cell);
                        laplacian_operator_loc *= scale;
                        scatter_bc_data(msh, cell, laplacian_operator_loc);
                }
            );
        #else
            auto storage = msh.backend_storage();
            for (auto& cell_ind : m_elements_with_bc_eges)
            {
                auto& cell = storage->surfaces[cell_ind];
                Matrix<T, Dynamic, Dynamic> laplacian_operator_loc = laplacian_operator(cell_ind, msh, cell);
                laplacian_operator_loc *= scale;
                scatter_bc_data(msh, cell, laplacian_operator_loc);
            }
        #endif
        
    }
            
    void assemble_rhs(const Mesh& msh, std::function<T(const typename Mesh::point_type& )> rhs_fun){
        
        RHS.setZero();
        #ifdef HAVE_INTEL_TBB
                size_t n_cells = msh.cells_size();
                tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                    [this,&msh,&rhs_fun] (size_t & cell_ind){
                        auto& cell = msh.backend_storage()->surfaces[cell_ind];
                        auto cell_basis   = make_scalar_monomial_basis(msh, cell, m_hho_di.cell_degree());
                        Matrix<T, Dynamic, 1> f_loc = make_rhs(msh, cell, cell_basis, rhs_fun);
                        this->scatter_rhs_data(msh, cell, f_loc);
                }
            );
        #else
            auto contribute = [this,&msh,&rhs_fun] (const typename Mesh::cell_type& cell){
                auto cell_basis   = make_scalar_monomial_basis(msh, cell, m_hho_di.cell_degree());
                Matrix<T, Dynamic, 1> f_loc = make_rhs(msh, cell, cell_basis, rhs_fun);
                this->scatter_rhs_data(msh, cell, f_loc);
            };
            
            for (auto& cell : msh){
                contribute(cell);
            }
        #endif
    }
            
    void assemble_mass(const Mesh& msh){
        
        MASS.setZero();
        size_t cell_ind = 0;
        for (auto& cell : msh)
        {
            Matrix<T, Dynamic, Dynamic> mass_matrix = mass_operator(cell_ind,msh,cell);
            scatter_mass_data(msh, cell, mass_matrix);
            cell_ind++;
        }
        finalize_mass();
    }
            
    Matrix<T, Dynamic, Dynamic> laplacian_operator(size_t & cell_ind, const Mesh& msh, const typename Mesh::cell_type& cell){

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
        acoustic_material_data<T> & material = m_material[cell_ind];
        return (1.0/material.rho())*(R_operator + S_operator);
    }
    
    Matrix<T, Dynamic, Dynamic> mass_operator(size_t & cell_ind, const Mesh& msh, const typename Mesh::cell_type& cell){

        acoustic_material_data<T> & material = m_material[cell_ind];
        auto scal_basis = disk::make_scalar_monomial_basis(msh, cell, m_hho_di.cell_degree());
        Matrix<T, Dynamic, Dynamic> mass_matrix = disk::make_mass_matrix(msh, cell, scal_basis);
        mass_matrix *= (1.0/(material.rho()*material.vp()*material.vp()));
        return mass_matrix;
    }
            
    void classify_cells(const Mesh& msh){

        m_elements_with_bc_eges.clear();
        size_t cell_ind = 0;
        for (auto& cell : msh)
        {
            auto face_list = faces(msh, cell);
            for (size_t face_i = 0; face_i < face_list.size(); face_i++)
            {
                auto fc = face_list[face_i];
                auto fc_id = msh.lookup(fc);
                bool is_dirichlet_Q = m_bnd.is_dirichlet_face(fc_id);
                if (is_dirichlet_Q)
                {
                    m_elements_with_bc_eges.push_back(cell_ind);
                    break;
                }
            }
            cell_ind++;
        }
    }
            
    void project_over_cells(const Mesh& msh, Matrix<T, Dynamic, 1> & x_glob, std::function<T(const typename Mesh::point_type& )> scal_fun){
        
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        size_t n_dof = n_cbs * msh.cells_size() + n_fbs * (m_n_edges - m_n_essential_edges);
        x_glob = Matrix<T, Dynamic, 1>::Zero(n_dof);
        for (auto& cell : msh)
        {
            Matrix<T, Dynamic, 1> x_proj_dof = project_function(msh, cell, m_hho_di.cell_degree(), scal_fun);
            scatter_cell_dof_data(msh, cell, x_glob, x_proj_dof);
        }
    }
            
    void project_over_faces(const Mesh& msh, Matrix<T, Dynamic, 1> & x_glob, std::function<T(const typename Mesh::point_type& )> scal_fun){

        for (auto& cell : msh)
        {
            auto fcs = faces(msh, cell);
            for (size_t i = 0; i < fcs.size(); i++)
            {
                auto face = fcs[i];
                auto fc_id = msh.lookup(face);
                bool is_dirichlet_Q = m_bnd.is_dirichlet_face(fc_id);
                if (is_dirichlet_Q)
                {
                    continue;
                }
                Matrix<T, Dynamic, 1> x_proj_dof = project_function(msh, face, m_hho_di.face_degree(), scal_fun);
                scatter_face_dof_data(msh, face, x_glob, x_proj_dof);
            }
        }
    }
            
    void expand_triplet(const Mesh& msh){
            
        size_t n_cells = msh.cells_size();
        m_f_l_indexes.reserve(n_cells);
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
            
        size_t first_ind = 0;
        size_t last_ind = 0;
        for (auto& cell : msh)
        {
            auto fcs = faces(msh, cell);
            last_ind = 0;
            for (size_t i = 0; i < fcs.size(); i++)
            {
                auto face = fcs[i];
                auto fc_id = msh.lookup(face);
                bool is_dirichlet_Q = m_bnd.is_dirichlet_face(fc_id);
                if (is_dirichlet_Q)
                {
                    continue;
                }
                last_ind += n_fbs;
            }
            last_ind = first_ind + n_cbs*n_cbs + 2.0*n_cbs*last_ind + last_ind*last_ind;
            m_f_l_indexes.push_back(std::make_pair(first_ind,last_ind-1));
            first_ind = last_ind;
        }

        // expanding global triplets
        m_triplets.resize(last_ind);
    }
            
    void clear(void)
    {
        m_elements_with_bc_eges.clear();
        m_f_l_indexes.clear();
        m_triplets.clear();
        m_mass_triplets.clear();
    }

    void finalize(void)
    {
        LHS.setFromTriplets( m_triplets.begin(), m_triplets.end() );
        m_triplets.clear();
    }
            
    void finalize_mass(void)
    {
        MASS.setFromTriplets( m_mass_triplets.begin(), m_mass_triplets.end() );
        m_mass_triplets.clear();
    }
            
    Matrix<T, Dynamic, 1>
    gather_dof_data(  const Mesh& msh, const typename Mesh::cell_type& cl,
                    const Matrix<T, Dynamic, 1>& x_glob) const
    {
        auto num_faces = howmany_faces(msh, cl);
        auto cell_ofs = disk::priv::offset(msh, cl);
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        
        Matrix<T, Dynamic, 1> x_local(n_cbs + num_faces * n_fbs );
        x_local.block(0, 0, n_cbs, 1) = x_glob.block(cell_ofs * n_cbs, 0, n_cbs, 1);
        auto fcs = faces(msh, cl);
        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto fc = fcs[i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
            const auto face_id                  = eid.second;

            if (m_bnd.is_dirichlet_face( face_id))
            {
                auto fb = disk::make_scalar_monomial_basis(msh, fc, m_hho_di.face_degree());
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb, m_hho_di.face_degree());
                auto velocity = m_bnd.dirichlet_boundary_func(face_id);
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, velocity, m_hho_di.face_degree());
                x_local.block(n_cbs + i * n_fbs, 0, n_fbs, 1) = mass.llt().solve(rhs);
            }
            else
            {
                auto face_ofs = disk::priv::offset(msh, fc);
                auto global_ofs = n_cbs * msh.cells_size() + m_compress_indexes.at(face_ofs)*n_fbs;
                x_local.block(n_cbs + i*n_fbs, 0, n_fbs, 1) = x_glob.block(global_ofs, 0, n_fbs, 1);
            }
        }
        return x_local;
    }
            
    void
    scatter_cell_dof_data(  const Mesh& msh, const typename Mesh::cell_type& cell,
                    Matrix<T, Dynamic, 1>& x_glob, Matrix<T, Dynamic, 1> x_proj_dof) const
    {
        auto cell_ofs = disk::priv::offset(msh, cell);
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        x_glob.block(cell_ofs * n_cbs, 0, n_cbs, 1) = x_proj_dof;
    }
    
    void
    scatter_face_dof_data(  const Mesh& msh, const typename Mesh::face_type& face,
                    Matrix<T, Dynamic, 1>& x_glob, Matrix<T, Dynamic, 1> x_proj_dof) const
    {
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(),Mesh::dimension);
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        size_t n_cells = msh.cells_size();
        auto face_offset = disk::priv::offset(msh, face);
        auto glob_offset = n_cbs * n_cells + m_compress_indexes.at(face_offset)*n_fbs;
        x_glob.block(glob_offset, 0, n_fbs, 1) = x_proj_dof;
    }
            
    void load_material_data(const Mesh& msh){
        m_material.clear();
        m_material.reserve(msh.cells_size());
        T rho = 1.0;
        T vp = 1.0;
        acoustic_material_data<T> material(rho,vp);
        for (size_t cell_ind = 0; cell_ind < msh.cells_size(); cell_ind++)
        {
            m_material.push_back(material);
        }
    }
    
    void load_material_data(const Mesh& msh, acoustic_material_data<T> & material){
        m_material.clear();
        m_material.reserve(msh.cells_size());
        for (size_t cell_ind = 0; cell_ind < msh.cells_size(); cell_ind++)
        {
            m_material.push_back(material);
        }
    }
      
    void load_material_data(const Mesh& msh, std::function<std::vector<T>(const typename Mesh::point_type& )> acoustic_mat_fun){
        m_material.clear();
        m_material.reserve(msh.cells_size());
        for (auto& cell : msh)
        {
            auto bar = barycenter(msh, cell);
            std::vector<T> mat_data = acoustic_mat_fun(bar);
            T rho = mat_data[0];
            T vp = mat_data[1];
            acoustic_material_data<T> material(rho,vp);
            m_material.push_back(material);
        }
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
            
    boundary_type & get_bc_conditions(){
             return m_bnd;
    }
            
    std::vector< acoustic_material_data<T> > & get_material_data(){
        return m_material;
    }
    
    size_t get_n_face_dof(){
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        size_t n_face_dof = (m_n_edges - m_n_essential_edges) * n_fbs;
        return n_face_dof;
    }
    
    size_t get_n_cell_basis(){
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(),Mesh::dimension);
        return n_cbs;
    }
    
    size_t get_n_face_basis(){
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        return n_fbs;
    }
    
    std::vector< std::pair<size_t,size_t> > & get_f_l_indexes(){
        return m_f_l_indexes;
    }
    
    std::vector< Triplet<T> > & get_triplets(){
        return m_triplets;
    }
    
    size_t get_cell_basis_data(){
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(),Mesh::dimension);
        return n_cbs;
    }
    
};

#endif /* acoustic_one_field_assembler_hpp */
