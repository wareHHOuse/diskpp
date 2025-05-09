//
//  acoustic_two_fields_assembler.hpp
//  acoustics
//
//  Created by Omar Durán on 4/21/20.
//

#pragma once
#ifndef acoustic_two_fields_assembler_hpp
#define acoustic_two_fields_assembler_hpp

#include "diskpp/bases/bases.hpp"
#include "diskpp/methods/hho"
#include "../common/assembly_index.hpp"
#include "../common/acoustic_material_data.hpp"

#ifdef HAVE_INTEL_TBB
#include <tbb/parallel_for.h>
#endif

template<typename Mesh>
class acoustic_two_fields_assembler
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
    
    size_t      m_n_edges;
    size_t      m_n_essential_edges;
    bool        m_hho_stabilization_Q;
    bool        m_scaled_stabilization_Q;

public:

    SparseMatrix<T>         LHS;
    Matrix<T, Dynamic, 1>   RHS;
    SparseMatrix<T>         MASS;
            
    acoustic_two_fields_assembler(const Mesh& msh, const disk::hho_degree_info& hho_di, const boundary_type& bnd)
        : m_hho_di(hho_di), m_bnd(bnd), m_hho_stabilization_Q(true), m_scaled_stabilization_Q(false)
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

        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_cbs  = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_cbs = n_scal_cbs + n_vec_cbs;
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);

        size_t system_size = n_cbs * msh.cells_size() + n_fbs * (m_n_edges - m_n_essential_edges);

        LHS = SparseMatrix<T>( system_size, system_size );
        RHS = Matrix<T, Dynamic, 1>::Zero( system_size );
        MASS = SparseMatrix<T>( system_size, system_size );
        
        classify_cells(msh);
    }

    void scatter_data(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs,
             const Matrix<T, Dynamic, 1>& rhs)
    {
    
        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_cbs = n_scal_cbs + n_vec_cbs;
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        auto fcs = faces(msh, cl);
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
            RHS(size_t(asm_map[i])) += rhs(i);
        }

    }
            
    void scatter_bc_data(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs)
    {
    
        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_cbs = n_scal_cbs + n_vec_cbs;
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        auto fcs = faces(msh, cl);
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
    
        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_cbs = n_scal_cbs + n_vec_cbs;
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
        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_cbs = n_scal_cbs + n_vec_cbs;
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

    void assemble(const Mesh& msh, std::function<double(const typename Mesh::point_type& )> rhs_fun){
        
        LHS.setZero();
        RHS.setZero();
        size_t cell_ind = 0;
        for (auto& cell : msh)
        {
            Matrix<T, Dynamic, 1> f_loc = mixed_rhs(msh, cell, rhs_fun);
            Matrix<T, Dynamic, Dynamic> mixed_operator_loc = mixed_operator(cell_ind, msh, cell);
            scatter_data(msh, cell, mixed_operator_loc, f_loc);
            cell_ind++;
        }
        finalize();
    }
    
    Matrix<T, Dynamic, Dynamic> mixed_operator(size_t & cell_ind, const Mesh& msh, const typename Mesh::cell_type& cell){
        acoustic_material_data<T> & material = m_material[cell_ind];
        auto reconstruction_operator   = mixed_scalar_reconstruction(msh, cell);
        Matrix<T, Dynamic, Dynamic> R_operator = reconstruction_operator.second;
        auto n_rows = R_operator.rows();
        auto n_cols = R_operator.cols();

        Matrix<T, Dynamic, Dynamic> S_operator = Matrix<T, Dynamic, Dynamic>::Zero(n_rows, n_cols);
        if(m_hho_stabilization_Q)
        {
            auto stabilization_operator    = make_scalar_hho_stabilization(msh, cell, reconstruction_operator.first, m_hho_di, m_scaled_stabilization_Q);
            auto n_s_rows = stabilization_operator.rows();
            auto n_s_cols = stabilization_operator.cols();
            S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
        }else{
            auto stabilization_operator    = make_scalar_hdg_stabilization(msh, cell, m_hho_di, m_scaled_stabilization_Q);
            auto n_s_rows = stabilization_operator.rows();
            auto n_s_cols = stabilization_operator.cols();
            S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
        }

        T rho = material.rho();
        T vp = material.vp();
        return R_operator + ((1.0)/(vp*rho))*S_operator;
    }
            
    void apply_bc(const Mesh& msh){
        
        #ifdef HAVE_INTEL_TBB
                size_t n_cells = m_elements_with_bc_eges.size();
                tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                    [this,&msh] (size_t & i){
                        size_t cell_ind = m_elements_with_bc_eges[i];
                        auto& cell = msh.backend_storage()->surfaces[cell_ind];
                        Matrix<T, Dynamic, Dynamic> mixed_operator_loc = mixed_operator(cell_ind, msh, cell);
                        scatter_bc_data(msh, cell, mixed_operator_loc);
                }
            );
        #else
            auto storage = msh.backend_storage();
            for (auto& cell_ind : m_elements_with_bc_eges)
            {
                auto& cell = storage->surfaces[cell_ind];
                Matrix<T, Dynamic, Dynamic> mixed_operator_loc = mixed_operator(cell_ind, msh, cell);
                scatter_bc_data(msh, cell, mixed_operator_loc);
            }
        #endif
        
    }
            
    void assemble_rhs(const Mesh& msh, std::function<double(const typename Mesh::point_type& )> & rhs_fun, size_t di = 0){
        
        RHS.setZero();
         
    #ifdef HAVE_INTEL_TBB
            size_t n_cells = msh.cells_size();
            tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                [this,&msh,&rhs_fun,&di] (size_t & cell_ind){
                    auto& cell = msh.backend_storage()->surfaces[cell_ind];
                    Matrix<T, Dynamic, 1> f_loc = this->mixed_rhs(msh, cell, rhs_fun,di);
                    this->scatter_rhs_data(msh, cell, f_loc);
            }
        );
    #else
        auto contribute = [this,&msh,&rhs_fun,&di] (const typename Mesh::cell_type& cell){
            Matrix<T, Dynamic, 1> f_loc = this->mixed_rhs(msh, cell, rhs_fun,di);
            this->scatter_rhs_data(msh, cell, f_loc);
        };
        
        for (auto& cell : msh){
            contribute(cell);
        }
    #endif
    }
            
    void assemble_mass(const Mesh& msh, bool add_scalar_mass_Q = true){
        
        MASS.setZero();
        for (size_t cell_ind = 0; cell_ind < msh.cells_size(); cell_ind++)
        {
            auto& cell = msh.backend_storage()->surfaces[cell_ind];
            Matrix<T, Dynamic, Dynamic> mass_matrix = mass_operator(cell_ind, msh, cell, add_scalar_mass_Q);
            scatter_mass_data(msh, cell, mass_matrix);
        }
        finalize_mass();
    }
            
    Matrix<T, Dynamic, Dynamic> mass_operator(size_t & cell_ind, const Mesh& msh, const typename Mesh::cell_type& cell, bool add_scalar_mass_Q = true){
            
        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_cbs = n_scal_cbs + n_vec_cbs;
            
        acoustic_material_data<T> & material = m_material[cell_ind];
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
            
    void classify_cells(const Mesh& msh){

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

    void project_over_cells(const Mesh& msh, Matrix<T, Dynamic, 1> & x_glob, std::function<double(const typename Mesh::point_type& )> scal_fun, std::function<std::vector<double>(const typename Mesh::point_type& )> vec_fun){
        size_t n_dof = MASS.rows();
        x_glob = Matrix<T, Dynamic, 1>::Zero(n_dof);
        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_cbs = n_scal_cbs + n_vec_cbs;
        for (auto& cell : msh)
        {
    
            Matrix<T, Dynamic, 1> x_proj_vec_dof = project_vec_function(msh, cell, vec_fun);
            Matrix<T, Dynamic, 1> x_proj_sca_dof = project_function(msh, cell, m_hho_di.cell_degree(), scal_fun);
            
            Matrix<T, Dynamic, 1> x_proj_dof = Matrix<T, Dynamic, 1>::Zero(n_cbs);
            x_proj_dof.block(0, 0, n_vec_cbs, 1) = x_proj_vec_dof;
            x_proj_dof.block(n_vec_cbs, 0, n_scal_cbs, 1) = x_proj_sca_dof;
            scatter_cell_dof_data(msh, cell, x_glob, x_proj_dof);
        }
    }
            
    Matrix<T, Dynamic, 1> project_vec_function(const Mesh& msh, const typename Mesh::cell_type& cell,
                      std::function<std::vector<double>(const typename Mesh::point_type& )> vec_fun){
    
            auto recdeg = m_hho_di.reconstruction_degree();
            auto rec_basis = make_scalar_monomial_basis(msh, cell, recdeg);
            auto rbs = disk::scalar_basis_size(recdeg, Mesh::dimension);
            Matrix<T, Dynamic, Dynamic> mass_matrix_q_full  = make_stiffness_matrix(msh, cell, rec_basis);
            Matrix<T, Dynamic, Dynamic> mass_matrix_q = Matrix<T, Dynamic, Dynamic>::Zero(rbs-1, rbs-1);
            mass_matrix_q = mass_matrix_q_full.block(1, 1, rbs-1, rbs-1);

            Matrix<T, Dynamic, 1> rhs = Matrix<T, Dynamic, 1>::Zero(rbs-1);
            Matrix<T, 1, 2> f_vec = Matrix<T, Dynamic, Dynamic>::Zero(1, 2);
            const auto qps = integrate(msh, cell, 2*recdeg);
            for (auto& qp : qps)
            {
              auto dphi = rec_basis.eval_gradients(qp.point());
              std::vector<double> flux = vec_fun(qp.point());
              f_vec(0,0) = flux[0];
              f_vec(0,1) = flux[1];
              for (size_t i = 0; i < rbs-1; i++){
              Matrix<T, 2, 1> phi_i = dphi.block(i+1, 0, 1, 2).transpose();
                  rhs(i,0) = rhs(i,0) + (qp.weight() * f_vec*phi_i)(0,0);
              }
            }
            Matrix<T, Dynamic, 1> x_dof = mass_matrix_q.llt().solve(rhs);
            return x_dof;
    }
            
    void project_over_faces(const Mesh& msh, Matrix<T, Dynamic, 1> & x_glob, std::function<double(const typename Mesh::point_type& )> scal_fun){

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
        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_cbs = n_scal_cbs + n_vec_cbs;
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        
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
                auto fb = disk::make_scalar_monomial_basis(msh, fc, m_hho_di.face_degree());
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb, m_hho_di.face_degree());
                auto velocity = m_bnd.dirichlet_boundary_func(face_id);
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, velocity, m_hho_di.face_degree());
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
       
    void scatter_cell_dof_data(  const Mesh& msh, const typename Mesh::cell_type& cell,
                    Matrix<T, Dynamic, 1>& x_glob, Matrix<T, Dynamic, 1> &x_proj_dof) const
    {
        auto cell_ofs = disk::priv::offset(msh, cell);
        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_cbs = n_scal_cbs + n_vec_cbs;
        x_glob.block(cell_ofs * n_cbs, 0, n_cbs, 1) = x_proj_dof;
    }
    
    void scatter_face_dof_data(  const Mesh& msh, const typename Mesh::face_type& face,
                    Matrix<T, Dynamic, 1>& x_glob, Matrix<T, Dynamic, 1> &x_proj_dof) const
    {
        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_cbs = n_scal_cbs + n_vec_cbs;
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        size_t n_cells = msh.cells_size();
        auto face_offset = disk::priv::offset(msh, face);
        auto glob_offset = n_cbs * n_cells + m_compress_indexes.at(face_offset)*n_fbs;
        x_glob.block(glob_offset, 0, n_fbs, 1) = x_proj_dof;
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
    mixed_rhs(const Mesh& msh, const typename Mesh::cell_type& cell, std::function<double(const typename Mesh::point_type& )> & rhs_fun, size_t di = 0)
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
    
    void load_material_data(const Mesh& msh, acoustic_material_data<T> material){
        m_material.clear();
        m_material.reserve(msh.cells_size());
        for (size_t cell_ind = 0; cell_ind < msh.cells_size(); cell_ind++)
        {
            m_material.push_back(material);
        }
    }
      
    void load_material_data(const Mesh& msh, std::function<std::vector<double>(const typename Mesh::point_type& )> acoustic_mat_fun){
        m_material.clear();
        m_material.reserve(msh.cells_size());
        for (auto& cell : msh)
        {
            auto bar = barycenter(msh, cell);
            std::vector<double> mat_data = acoustic_mat_fun(bar);
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
      
    void set_scaled_stabilization(){
        m_scaled_stabilization_Q = true;
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
    
    size_t get_n_faces(){
        return m_n_edges - m_n_essential_edges;
    }
    
    size_t get_face_basis_data(){
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        return n_fbs;
    }
    
    size_t get_cell_basis_data(){
        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_cbs = n_scal_cbs + n_vec_cbs;
        return n_cbs;
    }
    
};

#endif /* acoustic_two_fields_assembler_hpp */
