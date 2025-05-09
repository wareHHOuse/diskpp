//
//  three_fields_vectorial_assembler.hpp
//  acoustics
//
//  Created by Omar Durán on 4/12/20.
//

#pragma once
#ifndef three_fields_vectorial_assembler_hpp
#define three_fields_vectorial_assembler_hpp

#include "diskpp/bases/bases.hpp"
#include "diskpp/methods/hho"

template<typename Mesh>
class three_fields_vectorial_assembler
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

    three_fields_vectorial_assembler(const Mesh& msh, const disk::hho_degree_info& hho_di, const boundary_type& bnd)
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

        size_t n_ten_cbs = disk::sym_matrix_basis_size(m_hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_sca_cbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t n_cbs = n_ten_cbs + n_sca_cbs + n_vec_cbs;
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
        size_t n_ten_cbs = disk::sym_matrix_basis_size(m_hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_sca_cbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t n_cbs = n_ten_cbs + n_sca_cbs + n_vec_cbs;
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

    }

    void assemble(const Mesh& msh, std::function<static_vector<double, 2>(const typename Mesh::point_type& )> rhs_fun){
        
        size_t n_ten_cbs = disk::sym_matrix_basis_size(m_hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_sca_cbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension);
        for (auto& cell : msh)
        {
            auto reconstruction_operator   = strain_tensor_reconstruction(msh, cell);
            Matrix<T, Dynamic, Dynamic> R_operator = reconstruction_operator.second;
            auto n_rows = R_operator.rows();
            auto n_cols = R_operator.cols();
            
            // Weak hydrostatic component
            auto divergence_operator = div_vector_reconstruction(msh, cell);
            Matrix<T, Dynamic, Dynamic> D_operator = divergence_operator.second;
            
            Matrix<T, Dynamic, Dynamic> S_operator = Matrix<T, Dynamic, Dynamic>::Zero(n_rows, n_cols);
            if(m_hho_stabilization_Q)
            {
                auto rec_for_stab   = make_vector_hho_symmetric_laplacian(msh, cell, m_hho_di);
                auto stabilization_operator    = make_vector_hho_stabilization(msh, cell, rec_for_stab.first, m_hho_di);
                auto n_s_rows = stabilization_operator.rows();
                auto n_s_cols = stabilization_operator.cols();
                S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
            }else{
                auto stabilization_operator    = make_vector_hdg_stabilization(msh, cell, m_hho_di);
                auto n_s_rows = stabilization_operator.rows();
                auto n_s_cols = stabilization_operator.cols();
                S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
            }
            
            // For the sake of clarity ...
            T lambda, mu;
            lambda = 1.0;
            mu = 1.0;
            
            R_operator.block(0, 0, n_ten_cbs, n_ten_cbs) *= (1.0/(2.0*mu));
            D_operator.block(0, 0, n_sca_cbs, n_sca_cbs) *= (1.0/(lambda));
            
            Matrix<T, Dynamic, Dynamic> mixed_elastic_hho_operator = R_operator + D_operator + 2.0*mu*S_operator;
            Matrix<T, Dynamic, 1> f_loc = mixed_rhs(msh, cell, rhs_fun);
            
            scatter_data(msh, cell, mixed_elastic_hho_operator, f_loc);
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
        size_t n_ten_cbs = disk::sym_matrix_basis_size(m_hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_sca_cbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t n_cbs = n_ten_cbs + n_sca_cbs + n_vec_cbs;
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        
        Matrix<T, Dynamic, 1> x_el(n_vec_cbs + num_faces * n_fbs );
        x_el.block(0, 0, n_vec_cbs, 1) = x_glob.block(cell_ofs * n_cbs + n_ten_cbs + n_sca_cbs, 0, n_vec_cbs, 1);
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
                x_el.block(n_vec_cbs + i * n_fbs, 0, n_fbs, 1) = mass.llt().solve(rhs);
            }
            else
            {
                auto face_ofs = disk::priv::offset(msh, fc);
                auto global_ofs = n_cbs * msh.cells_size() + m_compress_indexes.at(face_ofs)*n_fbs;
                x_el.block(n_vec_cbs + i*n_fbs, 0, n_fbs, 1) = x_glob.block(global_ofs, 0, n_fbs, 1);
            }
        }
        return x_el;
    }
            
    std::pair<   Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
                 Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>  >
    strain_tensor_reconstruction(const Mesh& msh, const typename Mesh::cell_type& cell)
    {

        using T        = typename Mesh::coordinate_type;
        typedef Matrix<T, Dynamic, Dynamic> matrix_type;

        const size_t N = Mesh::dimension;

        auto graddeg = m_hho_di.grad_degree();
        auto celdeg  = m_hho_di.cell_degree();
        auto facdeg  = m_hho_di.face_degree();

        auto ten_b = make_sym_matrix_monomial_basis(msh, cell, graddeg);
        auto vec_b = make_vector_monomial_basis(msh, cell, celdeg);

        auto ten_bs = disk::sym_matrix_basis_size(graddeg, Mesh::dimension, Mesh::dimension);
        auto sca_bs = disk::scalar_basis_size(facdeg, Mesh::dimension);
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

         for (auto& qp : qps)
         {
             const auto gphi = ten_b.eval_functions(qp.point());

             for (size_t j = 0; j < ten_bs; j++)
             {
                 
                auto qp_gphi_j = disk::priv::inner_product(qp.weight(), gphi[j]);
                for (size_t i = j; i < ten_bs; i += dec){
                         gr_lhs(i, j) += disk::priv::inner_product(gphi[i], qp_gphi_j);
                }
             }
         }

        // upper part
         for (size_t j = 0; j < ten_bs; j++)
             for (size_t i = 0; i < j; i++)
                 gr_lhs(i, j) = gr_lhs(j, i);

         // compute rhs
         if (celdeg > 0)
         {
             const auto qpc = integrate(msh, cell, graddeg + celdeg - 1);
             for (auto& qp : qpc)
             {
                 const auto gphi    = ten_b.eval_functions(qp.point());
                 const auto dphi    = vec_b.eval_sgradients(qp.point());
                 const auto qp_dphi = disk::priv::inner_product(qp.weight(), dphi);

                 gr_rhs.block(0, 0, ten_bs, vec_bs) += disk::priv::outer_product(gphi, qp_dphi);

             } // end qp
         }

         const auto fcs = faces(msh, cell);
         for (size_t i = 0; i < fcs.size(); i++)
         {
             const auto fc = fcs[i];
             const auto n  = normal(msh, cell, fc);
             const auto fb = make_vector_monomial_basis(msh, fc, facdeg);

             const auto qps_f = integrate(msh, fc, graddeg + std::max(celdeg, facdeg));
             for (auto& qp : qps_f)
             {
                 const auto gphi = ten_b.eval_functions(qp.point());
                 const auto cphi = vec_b.eval_functions(qp.point());
                 const auto fphi = fb.eval_functions(qp.point());

                 const auto qp_gphi_n = disk::priv::inner_product(gphi, disk::priv::inner_product(qp.weight(), n));
                 gr_rhs.block(0, vec_bs + i * fbs, ten_bs, fbs) += disk::priv::outer_product(qp_gphi_n, fphi);
                 gr_rhs.block(0, 0, ten_bs, vec_bs) -= disk::priv::outer_product(qp_gphi_n, cphi);
             }
         }
            
        auto n_rows = gr_rhs.cols() + ten_bs + sca_bs;
        auto n_cols = gr_rhs.cols() + ten_bs + sca_bs;

        // Shrinking data
        matrix_type data_mixed = matrix_type::Zero(n_rows,n_cols);
        data_mixed.block(0, 0, ten_bs, ten_bs) = gr_lhs;
        data_mixed.block(0, (ten_bs + sca_bs), ten_bs, n_cols-(ten_bs + sca_bs)) = -gr_rhs;
        data_mixed.block((ten_bs + sca_bs), 0, n_rows-(ten_bs + sca_bs), ten_bs) = gr_rhs.transpose();

        matrix_type oper = gr_lhs.llt().solve(gr_rhs);
        return std::make_pair(oper, data_mixed);
    }
    
    std::pair<   Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
                 Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>  >
    div_vector_reconstruction(const Mesh& msh, const typename Mesh::cell_type& cell)
    {
        using T = typename Mesh::coordinate_type;
        typedef Matrix<T, Dynamic, Dynamic> matrix_type;
            
        auto graddeg = m_hho_di.grad_degree();
        auto facedeg  = m_hho_di.face_degree();
        auto ten_bs = disk::sym_matrix_basis_size(graddeg, Mesh::dimension, Mesh::dimension);
        auto sca_bs = disk::scalar_basis_size(facedeg, Mesh::dimension);
            
        auto cbas_s = make_scalar_monomial_basis(msh, cell, m_hho_di.face_degree());
        auto dr_lhs = make_mass_matrix(msh, cell, cbas_s);
        auto dr_rhs = make_hho_divergence_reconstruction_rhs(msh, cell, m_hho_di);

        assert( dr_lhs.rows() == sca_bs && dr_lhs.cols() == sca_bs );
            
        // Shrinking data
        auto n_rows = dr_rhs.cols() + ten_bs + sca_bs;
        auto n_cols = dr_rhs.cols() + ten_bs + sca_bs;
        matrix_type data_mixed = matrix_type::Zero(n_rows,n_cols);
        data_mixed.block(ten_bs, ten_bs, sca_bs, sca_bs) = dr_lhs;
        data_mixed.block(ten_bs, (ten_bs + sca_bs), sca_bs, n_cols-(ten_bs + sca_bs)) = -dr_rhs;
        data_mixed.block((ten_bs + sca_bs), ten_bs, n_rows-(ten_bs + sca_bs), sca_bs) = dr_rhs.transpose();

        matrix_type oper = dr_lhs.ldlt().solve(dr_rhs);
        return std::make_pair(oper, data_mixed);
    }
            
    Matrix<typename Mesh::coordinate_type, Dynamic, 1>
    mixed_rhs(const Mesh& msh, const typename Mesh::cell_type& cell, std::function<static_vector<double, 2>(const typename Mesh::point_type& )> rhs_fun, size_t di = 0)
    {
        auto recdeg = m_hho_di.grad_degree();
        auto celdeg = m_hho_di.cell_degree();
        auto facdeg = m_hho_di.face_degree();

        auto ten_bs = disk::sym_matrix_basis_size(recdeg, Mesh::dimension, Mesh::dimension);
        auto sca_bs = disk::scalar_basis_size(facdeg, Mesh::dimension);
        auto vec_bs = disk::vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
        size_t n_cbs = ten_bs + sca_bs + vec_bs;
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
        ret.block(ten_bs + sca_bs,0,vec_bs,1) = ret_loc;
        return ret;
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

#endif /* three_fields_vectorial_assembler_hpp */
