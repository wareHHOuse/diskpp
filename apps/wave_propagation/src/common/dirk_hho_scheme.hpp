//
//  dirk_hho_scheme.hpp
//  acoustics
//
//  Created by Omar Durán on 4/21/20.
//

#ifndef dirk_hho_scheme_hpp
#define dirk_hho_scheme_hpp

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include "../common/linear_solver.hpp"

template<typename T>
class dirk_hho_scheme {

    private:

    T m_scale;
    SparseMatrix<T> m_Mg;
    SparseMatrix<T> m_Kg;
    Matrix<T, Dynamic, 1> m_Fg;
    linear_solver<T>  m_analysis;
    std::vector<std::pair<size_t,size_t>> m_cell_basis_data;
    size_t m_n_f_dof;
    bool m_global_sc_Q;
    bool m_iteraive_solver_Q;
    
    public:
    
    dirk_hho_scheme(SparseMatrix<T> & Kg, Matrix<T, Dynamic, 1> & Fg, SparseMatrix<T> & Mg){
        
        m_Mg = Mg;
        m_Kg = Kg;
        m_Fg = Fg;
        m_scale = 0.0;
        m_global_sc_Q = false;
        m_iteraive_solver_Q = false;
    }
    
    dirk_hho_scheme(SparseMatrix<T> & Kg, Matrix<T, Dynamic, 1> & Fg, SparseMatrix<T> & Mg, T scale){
        m_Mg = Mg;
        m_Kg = Kg;
        m_Fg = Fg;
        m_scale = scale;
        m_global_sc_Q = false;
        m_iteraive_solver_Q = false;
    }
    
    void set_static_condensation_data(std::pair<size_t,size_t> cell_basis_data, size_t n_f_dof){
        std::vector<std::pair<size_t,size_t>> vec_cell_basis_data(1);
        vec_cell_basis_data[0] = cell_basis_data;
        m_cell_basis_data = vec_cell_basis_data;
        m_n_f_dof = n_f_dof;
        m_global_sc_Q = true;
    }
    
    void set_static_condensation_data(std::vector<std::pair<size_t,size_t>> cell_basis_data, size_t n_f_dof){
        m_cell_basis_data = cell_basis_data;
        m_n_f_dof = n_f_dof;
        m_global_sc_Q = true;
    }
    
    void condense_equations() {
        SparseMatrix<T> K = m_Mg + m_scale * m_Kg;
        m_analysis.set_Kg(K,m_n_f_dof);
        m_analysis.condense_equations(m_cell_basis_data);
    }
    
    void ComposeMatrix() {
        if(m_global_sc_Q) 
            condense_equations();
        else{
            SparseMatrix<T> K = m_Mg + m_scale * m_Kg;
            m_analysis.set_Kg(K);
        }
    }
    
    void setIterativeSolver() {
        m_iteraive_solver_Q = true;
        m_analysis.set_iterative_solver(false, 1.0e-1);
    }
        
    void DecomposeMatrix() {
        m_analysis.factorize();
    }

    linear_solver<T> & DirkAnalysis() {
        return m_analysis;
    }
    
    SparseMatrix<T> & Mg(){
        return m_Mg;
    }
    
    SparseMatrix<T> & Kg(){
        return m_Kg;
    }
    
    Matrix<T, Dynamic, 1> & Fg(){
        return m_Fg;
    }
    
    void SetScale(T & scale){
        m_scale = scale;
    }
    
    void SetFg(Matrix<T, Dynamic, 1> & Fg) {
        m_Fg = Fg;
    }
    
    void irk_weight(Matrix<T, Dynamic, 1> & y, Matrix<T, Dynamic, 1> & k, T dt, T a, bool is_sdirk_Q){
    
        m_Fg -= Kg()*y;
        if (is_sdirk_Q) 
            k = DirkAnalysis().solve(m_Fg);
        else {
            T scale = a * dt;
            SetScale(scale);
            ComposeMatrix();
            if (m_iteraive_solver_Q) 
                DirkAnalysis().set_iterative_solver();
            DecomposeMatrix();
            k = DirkAnalysis().solve(m_Fg);
        }
    }
    
};

#endif /* dirk_hho_scheme_hpp */
