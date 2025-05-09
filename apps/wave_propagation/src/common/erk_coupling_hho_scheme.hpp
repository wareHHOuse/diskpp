#ifndef erk_coupling_hho_scheme_hpp
#define erk_coupling_hho_scheme_hpp

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include "../common/assembly_index.hpp"
#include <iomanip>  // Pour std::setw


template<typename T>
class erk_coupling_hho_scheme
{
    private:

    SparseMatrix<T> m_Mc;
    SparseMatrix<T> m_Kcc;
    SparseMatrix<T> m_Kcf;
    SparseMatrix<T> m_Kfc;
    SparseMatrix<T> m_Sff;
    SparseMatrix<T> m_Scc;
    SparseMatrix<T> m_Mc_inv;
    SparseMatrix<T> m_Sff_inv;
    SparseMatrix<T> m_inv_Sff;
    SparseMatrix<T> m_Cg;

    Matrix<T, Dynamic, 1> m_Fc;

    #ifdef HAVE_INTEL_MKL
        PardisoLDLT<Eigen::SparseMatrix<T>>  m_analysis_f;
    #else
        SimplicialLDLT<SparseMatrix<T>> m_analysis_f;
    #endif
    
    ConjugateGradient<SparseMatrix<T>> m_analysis_cg;
    
    size_t m_n_ec_dof;
    size_t m_n_ac_dof;
    size_t m_n_c_dof;
    size_t m_n_ef_dof;
    size_t m_n_af_dof;
    size_t m_n_f_dof;

    bool m_sff_is_block_diagonal_Q;
    bool m_iterative_solver_Q;
    
    public:
    
    erk_coupling_hho_scheme(SparseMatrix<T> & Kg, Matrix<T, Dynamic, 1> & Fg, SparseMatrix<T> & Mg, SparseMatrix<T> & Cg, size_t elastic_cell_dofs, size_t acoustic_cell_dofs, size_t e_face_dofs, size_t a_face_dofs){
        
        m_n_ec_dof = elastic_cell_dofs;
        m_n_ac_dof = acoustic_cell_dofs;
        // m_n_c_dof  = m_n_ec_dof + m_n_ac_dof;
        m_n_ef_dof = e_face_dofs;
        m_n_af_dof = a_face_dofs;
        m_n_f_dof  = e_face_dofs + a_face_dofs;
        m_n_c_dof  = Kg.rows() - m_n_f_dof;

        m_Mc  = Mg.block(        0,         0, m_n_c_dof, m_n_c_dof);
        m_Kcc = Kg.block(        0,         0, m_n_c_dof, m_n_c_dof);
        m_Kcf = Kg.block(        0, m_n_c_dof, m_n_c_dof, m_n_f_dof);
        m_Kfc = Kg.block(m_n_c_dof,         0, m_n_f_dof, m_n_c_dof);
        m_Sff = Kg.block(m_n_c_dof, m_n_c_dof, m_n_f_dof, m_n_f_dof);

        m_Cg  = Cg.block(m_n_c_dof, m_n_c_dof, m_n_f_dof, m_n_f_dof);
        
        m_Fc  = Fg.block(0, 0, m_n_c_dof, 1);

        m_sff_is_block_diagonal_Q   = true;
        m_iterative_solver_Q        = false;

    }
    
    void Mcc_inverse(size_t e_cells, size_t a_cells, size_t e_cbs, size_t a_cbs) {
                    
        size_t nnz_cc  = e_cbs*e_cbs*e_cells + a_cbs*a_cbs*a_cells;
        std::vector< Triplet<T> > triplets_cc;
        triplets_cc.resize(nnz_cc);
        m_Mc_inv = SparseMatrix<T>(m_n_c_dof, m_n_c_dof);
            
        #ifdef HAVE_INTEL_TBB
            
            tbb::parallel_for(size_t(0), size_t(e_cells), size_t(1),
                              [this,&triplets_cc,&e_cbs] (size_t & cell_ind) {
                    
                size_t stride_eq = cell_ind * e_cbs;
                size_t stride_l  = cell_ind * e_cbs * e_cbs;

                SparseMatrix<T> m_Mc_loc = m_Mc.block(stride_eq, stride_eq, e_cbs, e_cbs);
                SparseLU<SparseMatrix<T>> analysis_cc;
                analysis_cc.analyzePattern(m_Mc_loc);
                analysis_cc.factorize(m_Mc_loc);
                Matrix<T, Dynamic, Dynamic> m_Mc_inv_loc = analysis_cc.solve(Matrix<T, Dynamic, Dynamic>::Identity(e_cbs, e_cbs));
                    
                size_t l = 0;
                for (size_t i = 0; i < m_Mc_inv_loc.rows(); i++) {
                  for (size_t j = 0; j < m_Mc_inv_loc.cols(); j++) {
                    triplets_cc[stride_l+l] = Triplet<T>(stride_eq+i, stride_eq+j, m_Mc_inv_loc(i,j));
                    l++;
                  }
                }
            });
            tbb::parallel_for(size_t(0), size_t(a_cells), size_t(1), 
                              [this,&triplets_cc,&a_cbs,&e_cbs,&e_cells] (size_t & cell_ind) {
                    
                size_t stride_eq = cell_ind*a_cbs       + e_cells*e_cbs;
                size_t stride_l  = cell_ind*a_cbs*a_cbs + e_cells*e_cbs*e_cbs;
                        
                SparseMatrix<T> m_Mc_loc = m_Mc.block(stride_eq, stride_eq, a_cbs, a_cbs);
                SparseLU<SparseMatrix<T>> analysis_cc;
                analysis_cc.analyzePattern(m_Mc_loc);
                analysis_cc.factorize(m_Mc_loc);
                Matrix<T, Dynamic, Dynamic> m_Mc_inv_loc = analysis_cc.solve(Matrix<T, Dynamic, Dynamic>::Identity(a_cbs, a_cbs));
                
                size_t l = 0;
                for (size_t i = 0; i < m_Mc_inv_loc.rows(); i++) {
                  for (size_t j = 0; j < m_Mc_inv_loc.cols(); j++) {
                    triplets_cc[stride_l+l] = Triplet<T>(stride_eq+i, stride_eq+j, m_Mc_inv_loc(i,j));
                    l++;
                  }
                }
            });

        #else

            for (size_t cell_ind = 0; cell_ind < e_cells; cell_ind++) {
                size_t stride_eq = cell_ind * e_cbs;
                size_t stride_l  = cell_ind * e_cbs * e_cbs;
                        
                SparseMatrix<T> m_Mc_loc = m_Mc.block(stride_eq, stride_eq, e_cbs, e_cbs);
                SparseLU<SparseMatrix<T>> analysis_cc;
                analysis_cc.analyzePattern(m_Mc_loc);
                analysis_cc.factorize(m_Mc_loc);
                Matrix<T, Dynamic, Dynamic> m_Mc_inv_loc = analysis_cc.solve(Matrix<T, Dynamic, Dynamic>::Identity(e_cbs, e_cbs));
                
                size_t l = 0;
                for (size_t i = 0; i < m_Mc_inv_loc.rows(); i++) {
                  for (size_t j = 0; j < m_Mc_inv_loc.cols(); j++) {
                    triplets_cc[stride_l+l] = Triplet<T>(stride_eq+i, stride_eq+j, m_Mc_inv_loc(i,j));
                    l++;
                  }
                }
            }
            for (size_t cell_ind = 0; cell_ind < a_cells; cell_ind++) {
                size_t stride_eq = cell_ind*a_cbs       + e_cells*e_cbs;
                size_t stride_l  = cell_ind*a_cbs*a_cbs + e_cells*e_cbs*e_cbs;
                        
                SparseMatrix<T> m_Mc_loc = m_Mc.block(stride_eq, stride_eq, a_cbs, a_cbs);
                SparseLU<SparseMatrix<T>> analysis_cc;
                analysis_cc.analyzePattern(m_Mc_loc);
                analysis_cc.factorize(m_Mc_loc);
                Matrix<T, Dynamic, Dynamic> m_Mc_inv_loc = analysis_cc.solve(Matrix<T, Dynamic, Dynamic>::Identity(a_cbs, a_cbs));
                
                size_t l = 0;
                for (size_t i = 0; i < m_Mc_inv_loc.rows(); i++) {
                  for (size_t j = 0; j < m_Mc_inv_loc.cols(); j++) {
                    triplets_cc[stride_l+l] = Triplet<T>(stride_eq+i, stride_eq+j, m_Mc_inv_loc(i,j));
                    l++;
                  }
                }
            }

            #endif
            
            m_Mc_inv.setFromTriplets(triplets_cc.begin(), triplets_cc.end());
            triplets_cc.clear();
            return;

        }


    void Sff_inverse(size_t e_faces, size_t a_faces, size_t e_fbs, size_t a_fbs, std::vector<size_t> e_compress, std::vector<size_t> a_compress, std::set<size_t> elastic_internal_faces, std::set<size_t> acoustic_internal_faces, std::set<size_t> interfaces_index) {

        size_t n_interfaces = interfaces_index.size();                                          // Number of interfaces
        size_t nnz_ff = e_fbs*e_fbs*e_faces + a_fbs*a_fbs*a_faces + 2*e_fbs*a_fbs*n_interfaces; // Number of nonzeros
        std::vector< Triplet<T> > triplets_ff;
        triplets_ff.resize(nnz_ff);
        m_Sff_inv = SparseMatrix<T>(m_n_f_dof, m_n_f_dof);                                      // size: number of faces x number of faces

        // Inversion of elastic stabilization 
        for (size_t face_ind = 0; face_ind < e_faces; face_ind++) {
          // std::cout << "Elastic face: " << face_ind << std::endl << std::endl;
          size_t stride_eq = face_ind * e_fbs;
          size_t stride_l  = face_ind * e_fbs * e_fbs;
          SparseMatrix<T> S_ff_loc = m_Sff.block(stride_eq, stride_eq, e_fbs, e_fbs);
          SparseLU<SparseMatrix<T>> analysis_ff;
          analysis_ff.analyzePattern(S_ff_loc);
          analysis_ff.factorize(S_ff_loc);
          Matrix<T, Dynamic, Dynamic> S_ff_inv_loc = analysis_ff.solve(Matrix<T, Dynamic, Dynamic>::Identity(e_fbs, e_fbs));
          size_t l = 0;
          for (size_t i = 0; i < S_ff_inv_loc.rows(); i++) {
            for (size_t j = 0; j < S_ff_inv_loc.cols(); j++) {
              triplets_ff[stride_l+l] = Triplet<T>(stride_eq+i, stride_eq+j, S_ff_inv_loc(i,j));
              l++;
            }
          }
        }

        // Inversion of acoustic stabilization 
        for (size_t face_ind = 0; face_ind < a_faces; face_ind++) {
          // std::cout << "Acoutic face: " << e_faces + face_ind  << std::endl << std::endl; 
          size_t stride_eq = e_faces*e_fbs       + face_ind*a_fbs ;
          size_t stride_l  = e_faces*e_fbs*e_fbs + face_ind*a_fbs*a_fbs;   
          SparseMatrix<T> S_ff_loc = m_Sff.block(stride_eq, stride_eq, a_fbs, a_fbs);
          SparseLU<SparseMatrix<T>> analysis_ff;
          analysis_ff.analyzePattern(S_ff_loc);
          analysis_ff.factorize(S_ff_loc);
          Matrix<T, Dynamic, Dynamic> S_ff_inv_loc = analysis_ff.solve(Matrix<T, Dynamic, Dynamic>::Identity(a_fbs, a_fbs));  
          size_t l = 0;
          for (size_t i = 0; i < S_ff_inv_loc.rows(); i++) {
            for (size_t j = 0; j < S_ff_inv_loc.cols(); j++) {
              triplets_ff[stride_l+l] = Triplet<T>(stride_eq+i, stride_eq+j, S_ff_inv_loc(i,j));
              l++;
            }
          }
        }
          
        // Inversion of coupling terms 
        size_t cpt = 0;
        for (auto face : interfaces_index) {                                     // Parcours des interfaces
          size_t e_face_LHS_offset = e_compress.at(face)*e_fbs;                  // Indice de la face elastique
          size_t a_face_LHS_offset = e_faces*e_fbs + a_compress.at(face)*a_fbs;  // Indice de la face acoustique
          // std::cout << "Interface: " << face  << std::endl;
          // std::cout << "Elastic interface: "  << e_compress.at(face) << std::endl;
          // std::cout << "Acoustic interface: " << e_faces + a_compress.at(face) << std::endl << std::endl; 
          size_t fbs = e_fbs + a_fbs;
          size_t e_stride_l = e_compress.at(face)*e_fbs*e_fbs;
          size_t a_stride_l = e_faces*e_fbs*e_fbs + a_compress.at(face)*a_fbs*a_fbs;
          size_t i_stride_l = e_faces*e_fbs*e_fbs + a_faces*a_fbs*a_fbs + 2*cpt*a_fbs*e_fbs;
          
          // Extraction du bloc stabilisation local
          Matrix<T, Dynamic, Dynamic> dense_SC_ff(fbs, fbs);
          SparseMatrix<T> elastic_stab  = m_Sff.block(e_face_LHS_offset, e_face_LHS_offset, e_fbs, e_fbs);
          SparseMatrix<T> acoustic_stab = m_Sff.block(a_face_LHS_offset, a_face_LHS_offset, a_fbs, a_fbs);
          dense_SC_ff.block(0, 0, e_fbs, e_fbs)         = elastic_stab;
          dense_SC_ff.block(e_fbs, e_fbs, a_fbs, a_fbs) = acoustic_stab;

          // Extraction du bloc coupling
          SparseMatrix<T> coupling_ela  = m_Cg.block(e_face_LHS_offset, a_face_LHS_offset, e_fbs, a_fbs);
          SparseMatrix<T> coupling_acou = m_Cg.block(a_face_LHS_offset, e_face_LHS_offset, a_fbs, e_fbs);
          dense_SC_ff.block(0, e_fbs, e_fbs, a_fbs) = coupling_ela;
          dense_SC_ff.block(e_fbs, 0, a_fbs, e_fbs) = coupling_acou;

          // Inversion
          SparseMatrix<T> SC_ff_loc = dense_SC_ff.sparseView();
          SparseLU<SparseMatrix<T>> analysis_ff;
          analysis_ff.analyzePattern(SC_ff_loc);
          analysis_ff.factorize(SC_ff_loc);
          Matrix<T, Dynamic, Dynamic> SC_ff_inv_loc = analysis_ff.solve(Matrix<T, Dynamic, Dynamic>::Identity(fbs, fbs));

          size_t l = 0;
          for (size_t i = 0; i < e_fbs; i++) {
            for (size_t j = 0; j < e_fbs; j++) {
              triplets_ff[e_stride_l+l] = Triplet<T>(e_face_LHS_offset+i, e_face_LHS_offset+j, SC_ff_inv_loc(i,j));
              l++;
            }
          } 
          l = 0;
          for (size_t i = 0; i < a_fbs; i++) {
            for (size_t j = 0; j < a_fbs; j++) {
              triplets_ff[a_stride_l+l] = Triplet<T>(a_face_LHS_offset+i, a_face_LHS_offset+j, SC_ff_inv_loc(e_fbs+i,e_fbs+j));
              l++;
            }
          }
          l = 0;
          // Upper right bloc
          for (size_t i = 0; i < e_fbs; i++) {
            for (size_t j = 0; j < a_fbs; j++) {
              triplets_ff[i_stride_l+l] = Triplet<T>(e_face_LHS_offset+i, a_face_LHS_offset+j, SC_ff_inv_loc(i,e_fbs+j));
              l++;
            }
          }
          // Lower left bloc
          for (size_t i = 0; i < a_fbs; i++) {
            for (size_t j = 0; j < e_fbs; j++) {
              triplets_ff[i_stride_l+l] = Triplet<T>(a_face_LHS_offset+i, e_face_LHS_offset+j, SC_ff_inv_loc(e_fbs+i,j));
              l++;
            }
          }
          cpt++;
        }        

        m_Sff_inv.setFromTriplets(triplets_ff.begin(), triplets_ff.end());
        triplets_ff.clear();
        return;

    }

    
    void inverse_Sff() {

        // Bi CG
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solverBiCG;
        solverBiCG.compute(m_Sff);
        if (solverBiCG.info() != Eigen::Success) 
            std::cout << "Error: Matrix decomposition failed, the matrix may not be invertible";
        Eigen::SparseMatrix<double> identity(m_Sff.rows(), m_Sff.cols());
        identity.setIdentity();
        m_inv_Sff = solverBiCG.solve(identity);
        if (solverBiCG.info() != Eigen::Success) 
            std::cout << "Error: Solving the system failed, the matrix may not be invertible";

      return;

    }
    
    void setIterativeSolver(T tolerance = 1.0e-11){
        m_iterative_solver_Q = true;
        m_analysis_cg.setTolerance(tolerance);
    }
        
    void DecomposeFaceTerm(){
        
        if (m_iterative_solver_Q) {
            m_analysis_cg.compute(m_Sff);
            m_analysis_cg.setMaxIterations(m_Sff.rows());
        }

        else {
            m_analysis_f.analyzePattern(m_Sff);
            m_analysis_f.factorize(m_Sff);
        }
        m_sff_is_block_diagonal_Q = false;
    }
    
    void refresh_faces_unknowns(Matrix<T, Dynamic, 1> & x) {
      
      Matrix<T, Dynamic, 1> x_c_dof = x.block(0, 0, m_n_c_dof, 1);
    
      // Faces update from cells data
      Matrix<T, Dynamic, 1> RHSf = Kfc()*x_c_dof;
      if (m_sff_is_block_diagonal_Q) {
        x.block(m_n_c_dof, 0, m_n_f_dof, 1) = - m_Sff_inv * RHSf;
      }
      else { 
        inverse_Sff();
        x.block(m_n_c_dof, 0, m_n_f_dof, 1) = - m_inv_Sff * RHSf;
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
      Matrix<T, Dynamic, 1> RHSf = Kfc()*k_c_dof ;
      // To comment to test eta = 0
      if (m_sff_is_block_diagonal_Q) {
        k.block(m_n_c_dof, 0, m_n_f_dof, 1) = - m_Sff_inv * RHSf; 
      }
      else {
        k.block(m_n_c_dof, 0, m_n_f_dof, 1) = - m_inv_Sff * RHSf; 
      }
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

    void compute_eigenvalues(std::ostream & simulation_log = std::cout){

        SparseMatrix<T> A_SCHUR = m_Kcc - m_Kcf*m_Sff_inv*m_Kfc;
        SparseMatrix<T> A = A_SCHUR.transpose()*m_Mc_inv*A_SCHUR;
        Spectra::SparseSymMatProd<double> opA(A);
        Spectra::SparseCholesky<double>   opB(m_Mc);
        Spectra::SymGEigsSolver<Spectra::SparseSymMatProd<double>, Spectra::SparseCholesky<double>, Spectra::GEigsMode::Cholesky> eigs(opA, opB, 1, 10);
        // Initialize and compute
        eigs.init();
        eigs.compute(Spectra::SortRule::LargestMagn);
        std::cout << std::endl << bold << red << "   Computation of the eigenvalues: " << reset;
        std::cout << std::endl << bold << cyan << "      State of the computation: " << reset;
        bool debug = true;
        if (debug) {
            if(eigs.info() == Spectra::CompInfo::Successful)
                std::cout << cyan << bold << "Successful\n";
            if(eigs.info() == Spectra::CompInfo::NotComputed)
                std::cout << cyan << bold << "NotComputed\n";
            if(eigs.info() == Spectra::CompInfo::NotConverging)
                std::cout << cyan << bold << "NotConverging\n";
            if(eigs.info() == Spectra::CompInfo::NumericalIssue)
                std::cout << cyan << bold << "NumericalIssue\n";
        }
        eigs.eigenvalues();
        std::cout << bold << cyan << "      Eigenvalue found: " << reset << cyan << eigs.eigenvalues() << std::endl << std::endl; 
        simulation_log << "Eigenvalue found: " << eigs.eigenvalues() << std::endl;
        
    }

    void compute_eigenvalues_bis(SparseMatrix<T> LHS_STAB, std::pair<size_t,size_t> block_dimension, std::ostream & simulation_log = std::cout){
        
        auto ten_bs = block_dimension.first;         // Elastic block
        auto vec_cell_size = block_dimension.second; // Acoustic block
        SparseMatrix<T> m_Scc = LHS_STAB.block(0,0, m_n_c_dof, m_n_c_dof); // Stabilisation

        SparseMatrix<T> Delta = m_Kcf*m_Sff_inv*m_Kfc;
        SparseMatrix<T> A_SCHUR = m_Kcc - Delta;
        SparseMatrix<T> S_SCHUR = m_Scc - Delta;
        SparseMatrix<T> A = A_SCHUR.transpose()*m_Mc_inv*A_SCHUR - S_SCHUR.transpose()*m_Mc_inv*S_SCHUR;
        Spectra::SparseSymMatProd<double> opA(A);
        Spectra::SparseCholesky<double>   opB(m_Mc);
        Spectra::SymGEigsSolver<Spectra::SparseSymMatProd<double>, Spectra::SparseCholesky<double>, Spectra::GEigsMode::Cholesky> eigs(opA, opB, 1, 10);
        // Initialize and compute
        eigs.init();
        eigs.compute(Spectra::SortRule::LargestMagn);
        bool debug = true;
        if (debug) {
            if(eigs.info() == Spectra::CompInfo::Successful)
                std::cout << "Successful\n";
            if(eigs.info() == Spectra::CompInfo::NotComputed)
                std::cout << "NotComputed\n";
            if(eigs.info() == Spectra::CompInfo::NotConverging)
                std::cout << "NotConverging\n";
            if(eigs.info() == Spectra::CompInfo::NumericalIssue)
                std::cout << "NumericalIssue\n";
        }
        eigs.eigenvalues();
        std::cout << std::endl;
        std::cout << bold << red << "   Eigenvalue found: " << reset << eigs.eigenvalues();
        simulation_log << "Eigenvalue found: " << eigs.eigenvalues() << std::endl;
        
    }

// Delta = A_{TF} A_{FF}^{-1} A_{FT}
// A' = A_{TT}-Delta
// S'=S_{TT}-Delta

// Pb spectral pour A'MA' - S'MS'

};




#endif /* erk_hho_scheme_hpp */
