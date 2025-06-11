//
//  two_fields_assembler.hpp
//  acoustics
//
//  Created by Omar Durán on 4/10/20.
//

#ifndef two_fields_assembler_hpp
#define two_fields_assembler_hpp

#include "diskpp/bases/bases.hpp"
#include "diskpp/methods/hho"

template<typename Mesh>
class two_fields_assembler
{
    
    
    typedef disk::BoundaryConditions<Mesh, true>    boundary_type;
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

    two_fields_assembler(const Mesh& msh, const disk::hho_degree_info& hho_di, const boundary_type& bnd)
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

        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_cbs = n_scal_cbs + n_vec_cbs;
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);

        size_t system_size = n_cbs * msh.cells_size() + n_fbs * (m_n_edges - m_n_essential_edges);

        LHS = SparseMatrix<T>( system_size, system_size );
        RHS = Matrix<T, Dynamic, 1>::Zero( system_size );
    }

    void scatter_data(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs,
             const Matrix<T, Dynamic, 1>& rhs)
    {
        auto fcs = faces(msh, cl);
        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_cbs = n_scal_cbs + n_vec_cbs;
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs + n_fbs*fcs.size());

        auto cell_offset        = disk::offset(msh, cl);
        auto cell_LHS_offset    = cell_offset * n_cbs;

        for (size_t i = 0; i < n_cbs; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
        
        Matrix<T, Dynamic, 1> dirichlet_data = Matrix<T, Dynamic, 1>::Zero(n_cbs + fcs.size()*n_fbs);
        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto face_offset = disk::offset(msh, fc);
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

    void assemble(const Mesh& msh, std::function<double(const typename Mesh::point_type& )> rhs_fun){
        
        for (auto& cell : msh)
        {
            auto reconstruction_operator   = mixed_scalar_reconstruction(msh, cell);
            Matrix<T, Dynamic, Dynamic> R_operator = reconstruction_operator.second;
            auto n_rows = R_operator.rows();
            auto n_cols = R_operator.cols();

            Matrix<T, Dynamic, Dynamic> S_operator = Matrix<T, Dynamic, Dynamic>::Zero(n_rows, n_cols);
            if(m_hho_stabilization_Q)
            {
                auto stabilization_operator    = make_scalar_hho_stabilization(msh, cell, reconstruction_operator.first, m_hho_di);
                auto n_s_rows = stabilization_operator.rows();
                auto n_s_cols = stabilization_operator.cols();
                S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
            }else{
                auto stabilization_operator    = make_scalar_hdg_stabilization(msh, cell, m_hho_di);
                auto n_s_rows = stabilization_operator.rows();
                auto n_s_cols = stabilization_operator.cols();
                S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
            }
    
            Matrix<T, Dynamic, Dynamic> mixed_operator_loc = R_operator + S_operator;
            Matrix<T, Dynamic, 1> f_loc = mixed_rhs(msh, cell, rhs_fun);
            
            scatter_data(msh, cell, mixed_operator_loc, f_loc);
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
        auto cell_ofs = disk::offset(msh, cl);
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
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
                auto face_ofs = disk::offset(msh, fc);
                auto global_ofs = n_cbs * msh.cells_size() + m_compress_indexes.at(face_ofs)*n_fbs;
                x_el.block(n_cbs + i*n_fbs, 0, n_fbs, 1) = x_glob.block(global_ofs, 0, n_fbs, 1);
            }
        }
        return x_el;
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
        data_mixed.block(0, 0, vec_cell_size, vec_cell_size) = gr_lhs;
        data_mixed.block(0, vec_cell_size, vec_cell_size, ncols-vec_cell_size) = -gr_rhs;
        data_mixed.block(vec_cell_size, 0, nrows-vec_cell_size, vec_cell_size) = gr_rhs.transpose();
        
        matrix_type oper = gr_lhs.llt().solve(gr_rhs);
        return std::make_pair(oper, data_mixed);
    }
            
    Matrix<typename Mesh::coordinate_type, Dynamic, 1>
    mixed_rhs(const Mesh& msh, const typename Mesh::cell_type& cell, std::function<double(const typename Mesh::point_type& )> rhs_fun, size_t di = 0)
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

#endif /* two_fields_assembler_hpp */
