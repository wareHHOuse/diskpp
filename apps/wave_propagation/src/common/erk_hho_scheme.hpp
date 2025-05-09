//
//  erk_hho_scheme.hpp
//  acoustics
//
//  Created by Omar Durán on 5/8/20.
//

#ifndef erk_hho_scheme_hpp
#define erk_hho_scheme_hpp

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

template<typename T>
class erk_hho_scheme
{
    private:

    SparseMatrix<T> m_Mc;
    SparseMatrix<T> m_Kcc;
    SparseMatrix<T> m_Kcf;
    SparseMatrix<T> m_Kfc;
    SparseMatrix<T> m_Sff;
    SparseMatrix<T> m_Mc_inv;
    SparseMatrix<T> m_Sff_inv;
    Matrix<T, Dynamic, 1> m_Fc;

    #ifdef HAVE_INTEL_MKL
        PardisoLDLT<Eigen::SparseMatrix<T>>  m_analysis_f;
    #else
        SimplicialLDLT<SparseMatrix<T>> m_analysis_f;
    #endif
    
    ConjugateGradient<SparseMatrix<T>> m_analysis_cg;
    
    size_t m_n_c_dof;
    size_t m_n_f_dof;
    bool m_sff_is_block_diagonal_Q;
    bool m_iterative_solver_Q;
    
    public:
    
    erk_hho_scheme(SparseMatrix<T> & Kg, Matrix<T, Dynamic, 1> & Fg, SparseMatrix<T> & Mg, size_t n_f_dof){
        
        
        m_n_c_dof = Kg.rows() - n_f_dof;
        m_n_f_dof = n_f_dof;
        
        // Building blocks
        m_Mc = Mg.block(0, 0, m_n_c_dof, m_n_c_dof);
        m_Kcc = Kg.block(0, 0, m_n_c_dof, m_n_c_dof);
        m_Kcf = Kg.block(0, m_n_c_dof, m_n_c_dof, n_f_dof);
        m_Kfc = Kg.block(m_n_c_dof,0, n_f_dof, m_n_c_dof);
        m_Sff = Kg.block(m_n_c_dof,m_n_c_dof, n_f_dof, n_f_dof);
        m_Fc = Fg.block(0, 0, m_n_c_dof, 1);
        m_sff_is_block_diagonal_Q   = false;
        m_iterative_solver_Q        = false;
    }
    
    void setIterativeSolver(T tolerance = 1.0e-11){
        m_iterative_solver_Q = true;
        m_analysis_cg.setTolerance(tolerance);
    }
        
    void DecomposeFaceTerm(){
        
        if (m_iterative_solver_Q) {
            m_analysis_cg.compute(m_Sff);
            m_analysis_cg.setMaxIterations(m_Sff.rows());
        }else{
            m_analysis_f.analyzePattern(m_Sff);
            m_analysis_f.factorize(m_Sff);
        }
        m_sff_is_block_diagonal_Q = false;
    }
    

    #ifdef HAVE_INTEL_MKL
        PardisoLDLT<Eigen::SparseMatrix<T>> & FacesAnalysis(){
            return m_analysis_f;
        }
    #else
        SimplicialLDLT<SparseMatrix<T>> & FacesAnalysis(){
            return m_analysis_f;
        }
    #endif
    
    SparseMatrix<T> & Mc(){
        return m_Mc;
    }

    SparseMatrix<T> & Kcc(){
        return m_Kcc;
    }
    
    SparseMatrix<T> & Kcf(){
        return m_Kcf;
    }
    
    SparseMatrix<T> & Kfc(){
        return m_Kfc;
    }
    
    SparseMatrix<T> & Sff(){
        return m_Sff;
    }
    
    Matrix<T, Dynamic, 1> & Fc(){
        return m_Fc;
    }
    
    void SetFg(Matrix<T, Dynamic, 1> & Fg){
        m_Fc = Fg.block(0, 0, m_n_c_dof, 1);
    }
    
    void Kcc_inverse(std::pair<size_t,size_t> cell_basis_data){
                
        size_t n_cells = cell_basis_data.first;
        size_t n_cbs   = cell_basis_data.second;
        size_t nnz_cc = n_cbs*n_cbs*n_cells;
        std::vector< Triplet<T> > triplets_cc;
        triplets_cc.resize(nnz_cc);
        m_Mc_inv = SparseMatrix<T>( m_n_c_dof, m_n_c_dof );
        #ifdef HAVE_INTEL_TBB
                tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                    [this,&triplets_cc,&n_cbs] (size_t & cell_ind){
                    
                    size_t stride_eq = cell_ind * n_cbs;
                    size_t stride_l = cell_ind * n_cbs * n_cbs;

                    SparseMatrix<T> m_Mc_loc = m_Mc.block(stride_eq, stride_eq, n_cbs, n_cbs);
                    SparseLU<SparseMatrix<T>> analysis_cc;
                    analysis_cc.analyzePattern(m_Mc_loc);
                    analysis_cc.factorize(m_Mc_loc);
                    Matrix<T, Dynamic, Dynamic> m_Mc_inv_loc = analysis_cc.solve(Matrix<T, Dynamic, Dynamic>::Identity(n_cbs, n_cbs));
            
                    size_t l = 0;
                    for (size_t i = 0; i < m_Mc_inv_loc.rows(); i++)
                    {
                        for (size_t j = 0; j < m_Mc_inv_loc.cols(); j++)
                        {
                            triplets_cc[stride_l+l] = Triplet<T>(stride_eq+i, stride_eq+j, m_Mc_inv_loc(i,j));
                            l++;
                        }
                    }
                }
            );
        #else

            for (size_t cell_ind = 0; cell_ind < n_cells; cell_ind++)
            {
                size_t stride_eq = cell_ind * n_cbs;
                size_t stride_l = cell_ind * n_cbs * n_cbs;
                
                SparseMatrix<T> m_Mc_loc = m_Mc.block(stride_eq, stride_eq, n_cbs, n_cbs);
                SparseLU<SparseMatrix<T>> analysis_cc;
                analysis_cc.analyzePattern(m_Mc_loc);
                analysis_cc.factorize(m_Mc_loc);
                Matrix<T, Dynamic, Dynamic> m_Mc_inv_loc = analysis_cc.solve(Matrix<T, Dynamic, Dynamic>::Identity(n_cbs, n_cbs));
        
                size_t l = 0;
                for (size_t i = 0; i < m_Mc_inv_loc.rows(); i++)
                {
                    for (size_t j = 0; j < m_Mc_inv_loc.cols(); j++)
                    {
                        triplets_cc[stride_l+l] = Triplet<T>(stride_eq+i, stride_eq+j, m_Mc_inv_loc(i,j));
                        l++;
                    }
                }

            }
        #endif
        m_Mc_inv.setFromTriplets( triplets_cc.begin(), triplets_cc.end() );
        triplets_cc.clear();
        return;

    }
    
    void Sff_inverse(std::pair<size_t,size_t> face_basis_data){
                
        size_t n_faces = face_basis_data.first;
        size_t n_fbs   = face_basis_data.second;
        size_t nnz_ff = n_fbs*n_fbs*n_faces;
        std::vector< Triplet<T> > triplets_ff;
        triplets_ff.resize(nnz_ff);
        m_Sff_inv = SparseMatrix<T>( m_n_f_dof, m_n_f_dof );
        #ifdef HAVE_INTEL_TBB
                tbb::parallel_for(size_t(0), size_t(n_faces), size_t(1),
                    [this,&triplets_ff,&n_fbs] (size_t & face_ind){
                    
                    size_t stride_eq = face_ind * n_fbs;
                    size_t stride_l = face_ind * n_fbs * n_fbs;

                    SparseMatrix<T> S_ff_loc = m_Sff.block(stride_eq, stride_eq, n_fbs, n_fbs);
                    SparseLU<SparseMatrix<T>> analysis_ff;
                    analysis_ff.analyzePattern(S_ff_loc);
                    analysis_ff.factorize(S_ff_loc);
                    Matrix<T, Dynamic, Dynamic> S_ff_inv_loc = analysis_ff.solve(Matrix<T, Dynamic, Dynamic>::Identity(n_fbs, n_fbs));
            
                    size_t l = 0;
                    for (size_t i = 0; i < S_ff_inv_loc.rows(); i++)
                    {
                        for (size_t j = 0; j < S_ff_inv_loc.cols(); j++)
                        {
                            triplets_ff[stride_l+l] = Triplet<T>(stride_eq+i, stride_eq+j, S_ff_inv_loc(i,j));
                            l++;
                        }
                    }
                }
            );
        #else

            for (size_t face_ind = 0; face_ind < n_faces; face_ind++)
            {
                size_t stride_eq = face_ind * n_fbs;
                size_t stride_l = face_ind * n_fbs * n_fbs;

                SparseMatrix<T> S_ff_loc = m_Sff.block(stride_eq, stride_eq, n_fbs, n_fbs);
                SparseLU<SparseMatrix<T>> analysis_ff;
                analysis_ff.analyzePattern(S_ff_loc);
                analysis_ff.factorize(S_ff_loc);
                Matrix<T, Dynamic, Dynamic> S_ff_inv_loc = analysis_ff.solve(Matrix<T, Dynamic, Dynamic>::Identity(n_fbs, n_fbs));
                size_t l = 0;
                for (size_t i = 0; i < S_ff_inv_loc.rows(); i++)
                {
                    for (size_t j = 0; j < S_ff_inv_loc.cols(); j++)
                    {
                        triplets_ff[stride_l+l] = Triplet<T>(stride_eq+i, stride_eq+j, S_ff_inv_loc(i,j));
                        l++;
                    }
                }

            }
        #endif
        m_Sff_inv.setFromTriplets( triplets_ff.begin(), triplets_ff.end() );
        triplets_ff.clear();
        m_sff_is_block_diagonal_Q = true;
        return;

    }
    
    void refresh_faces_unknowns(Matrix<T, Dynamic, 1> & x){
        
        Matrix<T, Dynamic, 1> x_c_dof = x.block(0, 0, m_n_c_dof, 1);
    
        // Faces update from cells data
        Matrix<T, Dynamic, 1> RHSf = Kfc()*x_c_dof;
        if (m_sff_is_block_diagonal_Q) {
            x.block(m_n_c_dof, 0, m_n_f_dof, 1) = - m_Sff_inv * RHSf;
        }else{
            if (m_iterative_solver_Q) {
                x.block(m_n_c_dof, 0, m_n_f_dof, 1) = -m_analysis_cg.solve(RHSf); // new state
            }else{
                x.block(m_n_c_dof, 0, m_n_f_dof, 1) = -FacesAnalysis().solve(RHSf); // new state
            }
        }
    
    }
    
    void erk_weight(Matrix<T, Dynamic, 1> & y, Matrix<T, Dynamic, 1> & k){
        
        k=y;
        Matrix<T, Dynamic, 1> y_c_dof = y.block(0, 0, m_n_c_dof, 1);
        Matrix<T, Dynamic, 1> y_f_dof = y.block(m_n_c_dof, 0, m_n_f_dof, 1);
        
        // Cells update
        Matrix<T, Dynamic, 1> RHSc = Fc() - Kcc()*y_c_dof - Kcf()*y_f_dof;
        Matrix<T, Dynamic, 1> k_c_dof = m_Mc_inv * RHSc;
        k.block(0, 0, m_n_c_dof, 1) = k_c_dof;
        
        // Faces update
        Matrix<T, Dynamic, 1> RHSf = Kfc()*k_c_dof;
        if (m_sff_is_block_diagonal_Q) {
            k.block(m_n_c_dof, 0, m_n_f_dof, 1) = - m_Sff_inv * RHSf;
        }else{
            if (m_iterative_solver_Q) {
                k.block(m_n_c_dof, 0, m_n_f_dof, 1) = -m_analysis_cg.solve(RHSf); // new state
                std::cout << "Number of iterations (CG): " << m_analysis_cg.iterations() << std::endl;
                std::cout << "Estimated error: " << m_analysis_cg.error() << std::endl;
            }else{
                k.block(m_n_c_dof, 0, m_n_f_dof, 1) = -FacesAnalysis().solve(RHSf); // new state
            }
        }
    
    }
    
};

#endif /* erk_hho_scheme_hpp */
