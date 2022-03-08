/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *  
 * Matteo Cicuttin (C) 2020,2021
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group  
 */

#include <iostream>
#include <iomanip>
#include <regex>

#include <sstream>
#include <string>

#include <unistd.h>

#include "mumps_new.hpp"

#include "bases/bases.hpp"
#include "quadratures/quadratures.hpp"
#include "methods/hho"
#include "methods/implementation_hho/curl.hpp"
#include "core/loaders/loader.hpp"
#include "output/silo.hpp"
#include "solvers/solver.hpp"


template<typename Mesh, typename ScalT = typename Mesh::coordinate_type>
class maxwell_assembler_condensed
{
    using CoordT = typename Mesh::coordinate_type;

    Mesh                    msh;
    disk::hho_degree_info   chdi, ghdi;

    std::vector<Triplet<ScalT>> triplets;
    std::vector<Triplet<ScalT>> triplets_cond1;
    std::vector<Triplet<ScalT>> triplets_cond2;
    std::vector<Triplet<ScalT>> triplets_decond1;
    std::vector<Triplet<ScalT>> triplets_decond2;
    std::vector<Triplet<ScalT>> triplets_decond3;

    std::vector<bool>       is_dirichlet;

    std::vector<size_t>     compress_table;
    std::vector<size_t>     expand_table;

    const size_t INVALID_OFFSET = (size_t) ~0;

    size_t cell_basis_size() const
    {
        return disk::vector_basis_size(chdi.cell_degree(), 3, 3);
    }

    size_t face_basis_size() const
    {
        if (Mesh::dimension == 3)
            return disk::vector_basis_size(chdi.face_degree(), 2, 2);
        if (Mesh::dimension == 2)
            return disk::scalar_basis_size(chdi.face_degree(), 1);
        return 0;
    }

    size_t get_system_offset(const Mesh& msh,
        const typename Mesh::face_type& fc)
    {
        auto face_num = offset(msh, fc);
        assert(face_num < msh.faces_size());
        assert(compress_table.size() == msh.faces_size());
        auto cnum = compress_table[face_num];
        assert(cnum != INVALID_OFFSET);
        auto fbs = face_basis_size();
        return cnum*fbs;
    }

    bool is_in_system(const Mesh& msh,
        const typename Mesh::face_type& fc)
    {
        auto bi = msh.boundary_info(fc);
        auto face_num = offset(msh, fc);
        assert(face_num < msh.faces_size());
        assert(is_dirichlet.size() == msh.faces_size());
        return not (bi.is_boundary() and is_dirichlet[face_num]);
    }

public:
    
    SparseMatrix<ScalT>         LHS;
    SparseMatrix<ScalT>         Cond1;
    SparseMatrix<ScalT>         Cond2;
    SparseMatrix<ScalT>         Decond1;
    SparseMatrix<ScalT>         Decond2;
    SparseMatrix<ScalT>         Decond3;
    Matrix<ScalT, Dynamic, 1>   RHS;
    Matrix<ScalT, Dynamic, 1>   X_prev;
    Matrix<ScalT, Dynamic, 1>   X_curr;

    size_t                  syssz;
    size_t                  sysfcs;
    size_t                  fullsz;

    maxwell_assembler_condensed(const Mesh& pmsh,
                      const disk::hho_degree_info& pchdi,
                      const std::vector<bool> p_is_dirichlet)
        : msh(pmsh), chdi(pchdi), is_dirichlet(p_is_dirichlet)
    {
        auto cbs = cell_basis_size();
        auto fbs = face_basis_size();

        using face_type = typename Mesh::face_type;
        auto in_system = [&](const face_type& fc) -> bool {
            auto ofs = offset(msh, fc);
            return not (msh.is_boundary(fc) and is_dirichlet.at(ofs));
        };

        sysfcs = std::count_if(msh.faces_begin(), msh.faces_end(), in_system);
        syssz = fbs*sysfcs;
        fullsz = cbs*msh.cells_size() + syssz; 

        LHS = SparseMatrix<ScalT>(syssz, syssz);
        Cond1 = SparseMatrix<ScalT>(syssz, fullsz);
        Cond2 = SparseMatrix<ScalT>(syssz, fullsz);
        Decond1 = SparseMatrix<ScalT>(fullsz, fullsz);
        Decond2 = SparseMatrix<ScalT>(fullsz, fullsz);
        Decond3 = SparseMatrix<ScalT>(fullsz, syssz);
        RHS = Matrix<ScalT, Dynamic, 1>::Zero(syssz);
        X_prev = Matrix<ScalT, Dynamic, 1>::Zero(fullsz);
        X_curr = Matrix<ScalT, Dynamic, 1>::Zero(fullsz);

        compress_table.resize( msh.faces_size(), INVALID_OFFSET);
        expand_table.resize( sysfcs );

        size_t face_i = 0;
        size_t compressed_ofs = 0;
        for (auto& fc : faces(msh))
        {
            assert(compressed_ofs <= face_i);
            if ( is_in_system(msh, fc) )
            {
                assert(face_i < compress_table.size());
                compress_table[face_i] = compressed_ofs;
                assert(compressed_ofs < expand_table.size());
                expand_table[compressed_ofs] = face_i;
                compressed_ofs++;
            }

            face_i++;
        }

        assert(face_i == msh.faces_size());
    }

    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<ScalT, Dynamic, 1> local,
             const Matrix<ScalT, Dynamic, 1> global)
    {
        auto cbs = cell_basis_size();
        auto fbs = face_basis_size();
    }

    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<ScalT, Dynamic, Dynamic>& lhsc,
             const Matrix<ScalT, Dynamic, Dynamic>& cond1,
             const Matrix<ScalT, Dynamic, Dynamic>& cond2,
             const Matrix<ScalT, Dynamic, Dynamic>& decond1,
             const Matrix<ScalT, Dynamic, Dynamic>& decond2,
             const Matrix<ScalT, Dynamic, Dynamic>& decond3,
             const Matrix<ScalT, Dynamic, 1>& x_prev,
             const Matrix<ScalT, Dynamic, 1>& x_curr)
    {
        auto cbs = cell_basis_size();
        auto cell_ofs = cbs*offset(msh, cl);
        auto faces_start = cbs*msh.cells_size();
        auto fbs = face_basis_size();

        auto fcs = faces(msh, cl);
        for (size_t fi = 0; fi < fcs.size(); fi++)
        {
            if ( not is_in_system(msh, fcs[fi]) )
                continue;

            auto cofsi = get_system_offset(msh, fcs[fi]);
            for (size_t fj = 0; fj < fcs.size(); fj++)
            {
                auto lofsi = fi*fbs;
                auto lofsj = fj*fbs;

                if ( not is_in_system(msh, fcs[fj]) ) 
                {
                    //RHS.segment(cofsi, fbs) += -lhsc.block(lofsi, lofsj, fbs, fbs)*dirichlet_data.segment(lofsj, fbs);
                    continue;
                }

                auto cofsj = get_system_offset(msh, fcs[fj]);
                for (size_t i = 0; i < fbs; i++)
                {
                    for(size_t j = 0; j < fbs; j++)
                    {
                        triplets.push_back( Triplet<ScalT>(cofsi+i, cofsj+j, lhsc(lofsi+i, lofsj+j)) );
                        triplets_cond1.push_back( Triplet<ScalT>(cofsi+i, faces_start+cofsj+j, cond1(lofsi+i, cbs+lofsj+j)) );
                        triplets_cond2.push_back( Triplet<ScalT>(cofsi+i, faces_start+cofsj+j, cond2(lofsi+i, cbs+lofsj+j)) );
                    }
                }
            }

            for (size_t i = 0; i < fbs; i++)
            {
                for (size_t j = 0; j < cbs; j++)
                {
                    auto lofsi = fi*fbs;
                    triplets_cond1.push_back( Triplet<ScalT>(cofsi+i, cell_ofs+j, cond1(lofsi+i, j)) );
                    triplets_cond2.push_back( Triplet<ScalT>(cofsi+i, cell_ofs+j, cond2(lofsi+i, j)) );
                }
            }

            X_prev.segment(faces_start+cofsi, fbs) = x_prev.segment(cbs+fi*fbs, fbs);
            X_curr.segment(faces_start+cofsi, fbs) = x_curr.segment(cbs+fi*fbs, fbs);
        }

        /* Decondensation stuff */
        for (size_t i = 0; i < cbs; i++)
        {
            for (size_t j = 0; j < cbs; j++)
            {
                triplets_decond1.push_back( Triplet<ScalT>(cell_ofs+i, cell_ofs+j, decond1(i, j)) );
                triplets_decond2.push_back( Triplet<ScalT>(cell_ofs+i, cell_ofs+j, decond2(i, j)) );
            }

            X_prev.segment(cell_ofs, cbs) = x_prev.head(cbs);
            X_curr.segment(cell_ofs, cbs) = x_curr.head(cbs);
        }

        for (size_t fj = 0; fj < fcs.size(); fj++)
        {
            if ( not is_in_system(msh, fcs[fj]) ) 
                continue;
            
            auto cofsj = get_system_offset(msh, fcs[fj]);
            
            for (size_t i = 0; i < cbs; i++)
            {
                for (size_t j = 0; j < fbs; j++)
                {
                    auto lofsj = fj*fbs;
                    triplets_decond1.push_back( Triplet<ScalT>(cell_ofs+i, faces_start+cofsj+j, decond1(i, cbs+lofsj+j)) );
                    triplets_decond2.push_back( Triplet<ScalT>(cell_ofs+i, faces_start+cofsj+j, decond2(i, cbs+lofsj+j)) );
                    triplets_decond3.push_back( Triplet<ScalT>(cell_ofs+i, cofsj+j, decond3(i, lofsj+j)) );
                }
            }
        }
    }

    void finalize(void)
    {
        LHS.setFromTriplets(triplets.begin(), triplets.end());
        Cond1.setFromTriplets(triplets_cond1.begin(), triplets_cond1.end());
        Cond2.setFromTriplets(triplets_cond2.begin(), triplets_cond2.end());
        Decond1.setFromTriplets(triplets_decond1.begin(), triplets_decond1.end());
        Decond2.setFromTriplets(triplets_decond2.begin(), triplets_decond2.end());
        Decond3.setFromTriplets(triplets_decond3.begin(), triplets_decond3.end());
        triplets.clear();
        triplets_cond1.clear();
        triplets_cond2.clear();
        triplets_decond1.clear();
        triplets_decond2.clear();
        triplets_decond3.clear();
        std::cout << "Solving for " << syssz << " DOFs. ";
        std::cout << "Matrix has " << LHS.nonZeros() << " nonzeros." << std::endl; 
    }

    disk::dynamic_vector<ScalT>
    get_expanded_solution(const Mesh& msh, disk::dynamic_vector<ScalT>& sol)
    {
        auto fbs = face_basis_size();

        disk::dynamic_vector<ScalT> ret = disk::dynamic_vector<ScalT>::Zero( fbs*msh.faces_size() );

        for (size_t i = 0; i < sysfcs; i++)
        {
            auto in_offset = i*fbs;
            auto out_offset = expand_table.at(i)*fbs;
            ret.segment(out_offset, fbs) = sol.segment(in_offset, fbs);
        }

        return ret;
    }

    disk::dynamic_vector<ScalT>
    get_element_dofs(const Mesh& msh, const typename Mesh::cell& cl, disk::dynamic_vector<ScalT>& sol)
    {
        auto fbs = face_basis_size();
        auto fcs = faces(msh, cl);
        disk::dynamic_vector<ScalT> ret = disk::dynamic_vector<ScalT>::Zero( fbs*fcs.size() );

        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto& fc = fcs[i];
            if ( not is_in_system(msh, fc) )
                continue;

            auto ofs = get_system_offset(msh, fc);
            ret.segment(i*fbs, fbs) = sol.segment(ofs, fbs);
        }

        return ret;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
vector_wave_solver(Mesh<T,3,Storage>& msh, size_t order)
{
    typedef Mesh<T,3,Storage>                   mesh_type;
    typedef typename mesh_type::cell_type       cell_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;
    
    disk::hho_degree_info chdi( { .rd = (size_t) order,
                                  .cd = (size_t) order,
                                  .fd = (size_t) order } );

    T alpha = 10.;
    T delta_t = 1e-12;
    T epsilon = 8.85e-12;
    T mu = 4e-7*M_PI;
    T sigma = 0;

    T beta = 0.25;
    T gamma = 0.5;

    T dt = delta_t;
    T dt2 = delta_t*delta_t;

    auto sol_fun = [&](const point_type& pt) -> Matrix<T, 3, 1> {
        Matrix<T, 3, 1> ret;
        ret(0) = 0.0;
        ret(1) = 0.0;
        ret(2) = std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
        return ret;
    };

    std::vector<bool> is_dirichlet(msh.faces_size(), true);
    maxwell_assembler_condensed<mesh_type> assm(msh, chdi, is_dirichlet);

    std::cout << "Assembling to triplets" << std::endl;
    for (auto& cl : msh)
    {
        auto cb = make_vector_monomial_basis(msh, cl, chdi.cell_degree());
        auto cbs = cb.size();
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();
        auto fbs = disk::vector_basis_size(chdi.face_degree(), 2, 2);
        auto nfd = num_faces * fbs; 
        auto CR = disk::curl_reconstruction_pk(msh, cl, chdi);
        Matrix<T, Dynamic, Dynamic> ST = disk::curl_hdg_stabilization(msh, cl, chdi);

        Matrix<T, Dynamic, Dynamic> mK = (1./mu)*(CR.second + alpha*ST);

        Matrix<T, Dynamic, Dynamic> mass = make_vector_mass_oper(msh, cl, chdi);
        Matrix<T, Dynamic, Dynamic> mM = epsilon*mass;
        Matrix<T, Dynamic, Dynamic> mC = sigma*mass;

        Matrix<T, Dynamic, Dynamic> mA = beta*dt2*mK + gamma*dt*mC + mM;
        Matrix<T, Dynamic, Dynamic> mP = 2.0*mM + (2.0*gamma - 1.0)*dt*mC - (0.5 - 2.0*beta + gamma)*dt2*mK;
        Matrix<T, Dynamic, Dynamic> mQ = -mM - (gamma - 1.0)*dt*mC - (0.5 + beta - gamma)*dt2*mK; 

        Matrix<T, Dynamic, Dynamic> mA_TT = mA.block(  0,   0, cbs, cbs);
        Matrix<T, Dynamic, Dynamic> mA_TF = mA.block(  0, cbs, cbs, nfd);
        Matrix<T, Dynamic, Dynamic> mA_FT = mA.block(cbs,   0, nfd, cbs);
        Matrix<T, Dynamic, Dynamic> mA_FF = mA.block(cbs, cbs, nfd, nfd);

        Matrix<T, Dynamic, Dynamic> mP_TT = mP.block(  0,   0, cbs, cbs);
        Matrix<T, Dynamic, Dynamic> mP_TF = mP.block(  0, cbs, cbs, nfd);
        Matrix<T, Dynamic, Dynamic> mP_FT = mP.block(cbs,   0, nfd, cbs);
        Matrix<T, Dynamic, Dynamic> mP_FF = mP.block(cbs, cbs, nfd, nfd);

        Matrix<T, Dynamic, Dynamic> mQ_TT = mQ.block(  0,   0, cbs, cbs);
        Matrix<T, Dynamic, Dynamic> mQ_TF = mQ.block(  0, cbs, cbs, nfd);
        Matrix<T, Dynamic, Dynamic> mQ_FT = mQ.block(cbs,   0, nfd, cbs);
        Matrix<T, Dynamic, Dynamic> mQ_FF = mQ.block(cbs, cbs, nfd, nfd);

        LDLT<Matrix<T, Dynamic, Dynamic>> mA_TT_ldlt(mA_TT);

        /* Condensation for x^{n} */
        Matrix<T, Dynamic, Dynamic> cond1 = Matrix<T, Dynamic, Dynamic>::Zero(nfd, cbs+nfd);
        cond1.block(0,   0, nfd, cbs) += mP_FT;
        cond1.block(0, cbs, nfd, nfd) += mP_FF;

        Matrix<T, Dynamic, Dynamic> tmp = Matrix<T, Dynamic, Dynamic>::Zero(cbs, cbs+nfd);
        tmp.block(0,   0, cbs, cbs) += mP_TT;
        tmp.block(0, cbs, cbs, nfd) += mP_TF;
        cond1 -= mA_FT*mA_TT_ldlt.solve(tmp);

        /* Decondensation for x^{n} */
        Matrix<T, Dynamic, Dynamic> decond1 = mA_TT_ldlt.solve(tmp);

        /* Condensation for x^{n-1} */
        Matrix<T, Dynamic, Dynamic> cond2 = Matrix<T, Dynamic, Dynamic>::Zero(nfd, cbs+nfd);
        cond2.block(0,   0, nfd, cbs) += mQ_FT;
        cond2.block(0, cbs, nfd, nfd) += mQ_FF;

        tmp = Matrix<T, Dynamic, Dynamic>::Zero(cbs, cbs+nfd);
        tmp.block(0,   0, cbs, cbs) += mQ_TT;
        tmp.block(0, cbs, cbs, nfd) += mQ_TF;
        cond2 -= mA_FT*mA_TT_ldlt.solve(tmp);

        /* Decondensation for x^{n-1} */
        Matrix<T, Dynamic, Dynamic> decond2 = mA_TT_ldlt.solve(tmp);

        /* Decondensation for x^{n+1} */
        Matrix<T, Dynamic, Dynamic> decond3 = -mA_TT_ldlt.solve(mA_TF);

        /* LHS */
        Matrix<T, Dynamic, Dynamic> mA_cond = mA_FF - mA_FT*mA_TT_ldlt.solve(mA_TF);

        Matrix<T, Dynamic, 1> x_curr = project_tangent(msh, cl, chdi, sol_fun);
        auto freq = 211.98528005809e6; //M_SQRT1_2 * M_PI * M_PI;
        Matrix<T, Dynamic, 1> x_prev = x_curr*std::cos(-2.0*M_PI*freq*dt);

        assm.assemble(msh, cl, mA_cond, cond1, cond2, decond1, decond2, decond3, x_prev, x_curr);
    }

    assm.finalize();

    mumps_solver<double> solver;

    //disk::dump_sparse_matrix(assm.LHS, "lhs_hho.txt");

    auto fact_start = std::chrono::steady_clock::now();
    solver.factorize(assm.LHS);
    auto fact_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = fact_end - fact_start;
    std::cout << "Factorization time: " << elapsed_seconds.count() << " s\n";


    disk::dynamic_vector<T> RHS = disk::dynamic_vector<T>::Zero(assm.RHS.rows());
    disk::dynamic_vector<T> X_prev = assm.X_prev;
    disk::dynamic_vector<T> X_curr = assm.X_curr;
    disk::dynamic_vector<T> X_next = disk::dynamic_vector<T>::Zero(assm.RHS.rows());
    disk::dynamic_vector<T> B_prev = disk::dynamic_vector<T>::Zero(assm.RHS.rows());
    disk::dynamic_vector<T> B_curr = disk::dynamic_vector<T>::Zero(assm.RHS.rows());
    disk::dynamic_vector<T> B_next = disk::dynamic_vector<T>::Zero(assm.RHS.rows());

    disk::dynamic_vector<T> X_next_decond;

    size_t num_ts = 5000;
    for (size_t i = 0; i < num_ts; i++)
    {
        auto start = std::chrono::steady_clock::now();
        std::cout << "Iteration " << i << std::endl;
        RHS = assm.Cond1 * X_curr
            + assm.Cond2 * X_prev;
            //+ dt2*beta*B_next
            //+ dt2*(0.5 + gamma - 2.0*beta)*B_curr
            //+ dt2*(0.5 - gamma + beta)*B_prev;

        X_next = solver.solve(RHS);

        X_next_decond = assm.Decond1*X_curr 
                      + assm.Decond2*X_prev
                      + assm.Decond3*X_next;


        if (i%100 == 0)
        {
            std::vector<scalar_type> data_ex( msh.cells_size() );
            std::vector<scalar_type> data_ey( msh.cells_size() );
            std::vector<scalar_type> data_ez( msh.cells_size() );
            size_t cell_i = 0;
            for (auto& cl : msh)
            {
                auto cb = make_vector_monomial_basis(msh, cl, chdi.cell_degree());
                auto cbs = cb.size();
                auto bar = barycenter(msh, cl);
                auto esolseg = X_next_decond.segment(cell_i*cbs, cbs);
                auto phi = cb.eval_functions(bar);
                auto ls = phi.transpose()*esolseg;
                data_ex[cell_i] = ls(0);
                data_ey[cell_i] = ls(1);
                data_ez[cell_i] = ls(2);
                cell_i++;
            }

            disk::silo_database silo_db;

            std::stringstream ss;
            ss << "maxwell_td_" << i << ".silo";
            silo_db.create(ss.str());
            silo_db.add_mesh(msh, "mesh");

            disk::silo_zonal_variable<scalar_type> ex("ex", data_ex);
            silo_db.add_variable("mesh", ex);

            disk::silo_zonal_variable<scalar_type> ey("ey", data_ey);
            silo_db.add_variable("mesh", ey);

            disk::silo_zonal_variable<scalar_type> ez("ez", data_ez);
            silo_db.add_variable("mesh", ez);
        }

        X_prev = X_curr;
        X_curr = X_next_decond;
        auto end = std::chrono::steady_clock::now();

        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "elapsed time: " << elapsed_seconds.count() << " s, ";
        std::cout << std::scientific << double(assm.LHS.rows())/elapsed_seconds.count() << " DoFs/s \n";
    }
}

int main(int argc, char **argv)
{
    rusage_monitor rm;

    using coord_T   = double;
    using scalar_T  = std::complex<double>;

    using T = double;

    T           stab_param = 1.0;
    T           omega = M_PI;
    bool        solve_eigvals = false;
    size_t      degree = 1;
    char *      mesh_filename = nullptr;
    char *      param_filename = nullptr;
    bool        use_ho_stab = false;

    int ch;
    while ( (ch = getopt(argc, argv, "a:ek:m:w:")) != -1 )
    {
        switch(ch)
        {
            case 'a':
                stab_param = std::stod(optarg);
                break;

            case 'e':
                solve_eigvals = true;
                break;

            case 'k':
                degree = std::stoi(optarg);
                break;

            case 'm':
                mesh_filename = optarg;
                break;

            case 'w':
                omega = M_PI*std::stod(optarg);
                break;

            case '?':
            default:
                std::cout << "Invalid option" << std::endl;
                return 1;
        }
    }

    if (mesh_filename == nullptr)
    {
        std::cout << "forgot -m?" << std::endl;
        return 1;
    }


    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        auto msh = disk::load_netgen_3d_mesh<coord_T>(mesh_filename);

        //if (solve_eigvals)
        //    maxwell_eigenvalue_solver(msh, degree, stab_param);
        //else
            //vector_wave_solver(msh, degree, stab_param, omega, false);
            vector_wave_solver(msh, degree);

        return 0;
    }

    /* GMSH */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo$") ))
    {
        std::cout << "Guessed mesh format: GMSH" << std::endl;
        disk::simplicial_mesh<T,3> msh;
        disk::gmsh_geometry_loader< disk::simplicial_mesh<T,3> > loader;
        //disk::generic_mesh<T,3> msh;
        //disk::gmsh_geometry_loader< disk::generic_mesh<T,3> > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        vector_wave_solver(msh, degree);

        return 0;
    }

    /* FVCA6 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
        disk::generic_mesh<coord_T,3> msh;
        
        disk::load_mesh_fvca6_3d<coord_T>(mesh_filename, msh);
        
        vector_wave_solver(msh, degree);
        
        return 0;
    }

}