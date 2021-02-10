/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *  
 * Matteo Cicuttin (C) 2020
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

#include "bases/bases.hpp"
#include "quadratures/quadratures.hpp"
#include "methods/hho"
#include "methods/implementation_hho/curl.hpp"
#include "core/loaders/loader.hpp"
#include "output/silo.hpp"
#include "solvers/solver.hpp"
#include "solvers/mumps.hpp"

#include "compinfo.h"

#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

template<typename Mesh>
class maxwell_assembler
{
    using T = typename Mesh::coordinate_type;

    Mesh                    msh;
    disk::hho_degree_info   chdi, ghdi;

    std::vector<Triplet<T>> triplets;

    std::vector<size_t>     compress_table;
    std::vector<size_t>     expand_table;

public:
    
    SparseMatrix<T>         LHS;
    Matrix<T, Dynamic, 1>   RHS;

    size_t                  syssz;
    size_t                  s_begin;

    maxwell_assembler(const Mesh& pmsh,
                      const disk::hho_degree_info& pchdi)
        : msh(pmsh), chdi(pchdi)
    {
        auto vcbs = disk::vector_basis_size(chdi.cell_degree(), Mesh::dimension, Mesh::dimension);
        auto vfbs = disk::vector_basis_size(chdi.face_degree(), Mesh::dimension-1, Mesh::dimension-1);
        auto v_end = vcbs*msh.cells_size() + vfbs*( msh.internal_faces_size() );
        auto scbs = disk::scalar_basis_size(ghdi.cell_degree(), Mesh::dimension);
        auto sfbs = disk::scalar_basis_size(ghdi.face_degree(), Mesh::dimension-1);
        auto s_end = v_end + scbs*msh.cells_size() + sfbs*( msh.internal_faces_size() );

        std::cout << vcbs << " " << vfbs << " " << scbs << " " << sfbs << std::endl;

        syssz = v_end;
        s_begin = v_end;

        LHS = SparseMatrix<T>(syssz, syssz);
        RHS = Matrix<T, Dynamic, 1>::Zero(syssz);

        compress_table.resize( msh.faces_size() );
        expand_table.resize( msh.internal_faces_size() );

        size_t face_i = 0;
        size_t compressed_ofs = 0;
        for (auto& fc : faces(msh))
        {
            if ( not msh.is_boundary(fc) )
            {
                compress_table[face_i] = compressed_ofs;
                expand_table[compressed_ofs] = face_i;
                compressed_ofs++;
            }

            face_i++;
        }
    }

    bool is_homo_dirichlet(size_t)
    {
        return true;
    }

    template<typename Function>
    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs,
             const Matrix<T, Dynamic, 1>& rhs,
             const Function& gD)
    {
        auto num_internal_faces = msh.internal_faces_size();
        auto vcbs = disk::vector_basis_size(chdi.cell_degree(), Mesh::dimension, Mesh::dimension);
        auto vfbs = disk::vector_basis_size(chdi.face_degree(), Mesh::dimension-1, Mesh::dimension-1);
        auto v_end = vcbs*msh.cells_size() + vfbs*num_internal_faces;

        auto v_cell_ofs = offset(msh, cl) * vcbs;

        auto fcs = faces(msh, cl);

        /* Curl-curl term, cell part */
        for (size_t i = 0; i < vcbs; i++)
        {
            auto gi = v_cell_ofs + i;
            for (size_t j = 0; j < vcbs; j++)
            {
                auto gj = v_cell_ofs + j;
                triplets.push_back( Triplet<T>(gi, gj, lhs(i,j)) );
            }

            RHS(gi) = rhs(i);
        }

        /* Offdiag parts */
        for (size_t i = 0; i < vcbs; i++)
        {
            auto gi = v_cell_ofs + i;
            for (size_t fj = 0; fj < fcs.size(); fj++)
            {
                if ( msh.is_boundary(fcs[fj]) )
                    continue;

                for (size_t j = 0; j < vfbs; j++)
                {
                    auto gj = vcbs * msh.cells_size() + vfbs * compress_table.at( offset(msh, fcs[fj]) ) + j;

                    auto li = i;
                    auto lj = vcbs + vfbs*fj + j;
                    triplets.push_back( Triplet<T>(gi, gj, lhs(li,lj)) );
                }
            }
        }

        for (size_t fi = 0; fi < fcs.size(); fi++)
        {
            if ( msh.is_boundary(fcs[fi]) )
                continue;

            for (size_t i = 0; i < vfbs; i++)
            {
                auto gi = vcbs * msh.cells_size() + vfbs * compress_table.at( offset(msh, fcs[fi]) ) + i;
                for (size_t j = 0; j < vcbs; j++)
                {
                    auto gj = v_cell_ofs + j;
                    auto li = vcbs + vfbs*fi + i;
                    auto lj = j;
                    triplets.push_back( Triplet<T>(gi, gj, lhs(li,lj)) );
                }
            }

            for (size_t fj = 0; fj < fcs.size(); fj++)
            {
                if ( msh.is_boundary(fcs[fj]) )
                {
                    auto bid = msh.boundary_id(fcs[fj]);
                    if ( is_homo_dirichlet(bid) )
                        continue;

                    using namespace std::placeholders;
                    auto gD_fun = std::bind(gD, bid, _1);
                    auto fb = make_vector_monomial_tangential_basis(msh, fcs[fj], chdi.face_degree());
                    disk::dynamic_vector<T> dirichlet_data = project_function(msh, fcs[fj], fb, gD_fun);
                    auto lofsi = vcbs+fi*vfbs;
                    auto lofsj = vcbs+fj*vfbs;
                    auto cofsj = vcbs*msh.cells_size() + vfbs * compress_table.at( offset(msh, fcs[fj]) );
                    RHS.segment(cofsj, vfbs) -= lhs.block(lofsi, lofsj, vfbs, vfbs)*dirichlet_data;
                }
                else
                {
                    auto cofsi = vcbs*msh.cells_size() + vfbs * compress_table.at( offset(msh, fcs[fi]) );
                    auto cofsj = vcbs*msh.cells_size() + vfbs * compress_table.at( offset(msh, fcs[fj]) );

                    for (size_t i = 0; i < vfbs; i++)
                    {
                        for(size_t j = 0; j < vfbs; j++)
                        {
                            auto lofsi = vcbs+fi*vfbs;
                            auto lofsj = vcbs+fj*vfbs;
                            triplets.push_back( Triplet<T>(cofsi+i, cofsj+j, lhs(lofsi+i, lofsj+j)) );
                        }
                    }
                }
            }
        }
    }

    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhsc,
             const Matrix<T, Dynamic, 1>& rhs)
    {
        auto gD = [&](size_t, const typename Mesh::point_type&) -> Matrix<T,3,1> {
            return  Matrix<T,3,1>::Zero();
        };
        assemble(msh, cl, lhsc, rhs, gD);
    }

    void finalize(void)
    {
        LHS.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();
        std::cout << "Solving for " << syssz << " DOFs. ";
        std::cout << "Matrix has " << LHS.nonZeros() << " nonzeros." << std::endl; 
    }
};

template<typename Mesh, typename ScalT = typename Mesh::coordinate_type>
class maxwell_assembler_condensed
{
    using CoordT = typename Mesh::coordinate_type;

    Mesh                    msh;
    disk::hho_degree_info   chdi, ghdi;

    std::vector<Triplet<ScalT>> triplets;
    std::vector<bool>       is_dirichlet;

    std::vector<size_t>     compress_table;
    std::vector<size_t>     expand_table;

    const size_t INVALID_OFFSET = (size_t) ~0;

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
    Matrix<ScalT, Dynamic, 1>   RHS;

    size_t                  syssz;
    size_t                  sysfcs;

    maxwell_assembler_condensed(const Mesh& pmsh,
                      const disk::hho_degree_info& pchdi,
                      const std::vector<bool> p_is_dirichlet)
        : msh(pmsh), chdi(pchdi), is_dirichlet(p_is_dirichlet)
    {
        auto fbs = face_basis_size();

        using face_type = typename Mesh::face_type;
        auto in_system = [&](const face_type& fc) -> bool {
            auto ofs = offset(msh, fc);
            return not (msh.is_boundary(fc) and is_dirichlet.at(ofs));
        };

        sysfcs = std::count_if(msh.faces_begin(), msh.faces_end(), in_system);
        syssz = fbs*sysfcs;

        LHS = SparseMatrix<ScalT>(syssz, syssz);
        RHS = Matrix<ScalT, Dynamic, 1>::Zero(syssz);

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
             const Matrix<ScalT, Dynamic, Dynamic>& lhsc,
             const Matrix<ScalT, Dynamic, 1>& rhs)
    {
        auto fbs = face_basis_size();

        auto fcs = faces(msh, cl);
        for (size_t fi = 0; fi < fcs.size(); fi++)
        {
            if ( not is_in_system(msh, fcs[fi]) )
                continue;

            auto cofsi = get_system_offset(msh, fcs[fi]);
            for (size_t fj = 0; fj < fcs.size(); fj++)
            {
                if ( not is_in_system(msh, fcs[fj]) ) 
                    continue;

                auto cofsj = get_system_offset(msh, fcs[fj]);

                for (size_t i = 0; i < fbs; i++)
                {
                    for(size_t j = 0; j < fbs; j++)
                    {
                        auto lofsi = fi*fbs;
                        auto lofsj = fj*fbs;
                        triplets.push_back( Triplet<ScalT>(cofsi+i, cofsj+j, lhsc(lofsi+i, lofsj+j)) );
                    }
                }
            }

            RHS.segment(cofsi, fbs) += rhs.segment(fi*fbs, fbs);
        }
    }

    void finalize(void)
    {
        LHS.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();
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



template<typename Mesh>
class maxwell_eigenvalue_assembler
{
    using T = typename Mesh::coordinate_type;

    Mesh                    msh;
    disk::hho_degree_info   hdi;

    std::vector<Triplet<T>> trpM, trpK;

    std::vector<size_t>     compress_table;
    std::vector<size_t>     expand_table;

public:
    SparseMatrix<T>         gM, gK;

    size_t                  syssz;

    maxwell_eigenvalue_assembler(const Mesh& pmsh, const disk::hho_degree_info& phdi)
        : msh(pmsh), hdi(phdi)
    {
        auto cbs = disk::vector_basis_size(hdi.cell_degree(), Mesh::dimension, Mesh::dimension);
        auto fbs = disk::vector_basis_size(hdi.face_degree(), Mesh::dimension-1, Mesh::dimension-1);
        syssz = cbs*msh.cells_size() + fbs*( msh.internal_faces_size() );

        gK = SparseMatrix<T>(syssz, syssz);
        gM = SparseMatrix<T>(syssz, syssz);

        compress_table.resize( msh.faces_size() );
        expand_table.resize( msh.internal_faces_size() );

        size_t face_i = 0;
        size_t compressed_ofs = 0;
        for (auto& fc : faces(msh))
        {
            if ( not msh.is_boundary(fc) )
            {
                compress_table[face_i] = compressed_ofs;
                expand_table[compressed_ofs] = face_i;
                compressed_ofs++;
            }

            face_i++;
        }
    }

    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& K,
             const Matrix<T, Dynamic, Dynamic>& M)
    {
        auto cbs = disk::vector_basis_size(hdi.cell_degree(), Mesh::dimension, Mesh::dimension);
        auto fbs = disk::vector_basis_size(hdi.face_degree(), Mesh::dimension-1, Mesh::dimension-1);

        auto cell_ofs = offset(msh, cl);
        auto cell_dofs_ofs = offset(msh, cl) * cbs;

        /* Cell-cell part */
        for (size_t i = 0; i < cbs; i++)
        {
            for (size_t j = 0; j < cbs; j++)
            {
                auto gi = cell_dofs_ofs + i;
                auto gj = cell_dofs_ofs + j;
                trpK.push_back( Triplet<T>(gi, gj, K(i,j)) );
                trpM.push_back( Triplet<T>(gi, gj, M(i,j)) );
            }
        }

        auto fcs = faces(msh, cl);

        /* Offdiag parts: Cell-face */
        for (size_t i = 0; i < cbs; i++)
        {
            auto gi = cell_dofs_ofs + i;
            for (size_t fj = 0; fj < fcs.size(); fj++)
            {
                if ( msh.is_boundary(fcs[fj]) )
                    continue;

                for (size_t j = 0; j < fbs; j++)
                {
                    auto gj = cbs * msh.cells_size() + fbs * compress_table.at( offset(msh, fcs[fj]) ) + j;

                    auto li = i;
                    auto lj = cbs + fbs*fj + j;
                    trpK.push_back( Triplet<T>(gi, gj, K(li,lj)) );
                }
            }
        }

        /* Offdiag parts: Face-cell */
        for (size_t fi = 0; fi < fcs.size(); fi++)
        {
            if ( msh.is_boundary(fcs[fi]) )
                continue;

            for (size_t i = 0; i < fbs; i++)
            {
                auto gi = cbs * msh.cells_size() + fbs * compress_table.at( offset(msh, fcs[fi]) ) + i;
                for (size_t j = 0; j < cbs; j++)
                {
                    auto gj = cell_dofs_ofs + j;
                    auto li = cbs + fbs*fi + i;
                    auto lj = j;
                    trpK.push_back( Triplet<T>(gi, gj, K(li,lj)) );
                }
            }
        }

        /* Face-face part */
        for (size_t fi = 0; fi < fcs.size(); fi++)
        {
            if ( msh.is_boundary(fcs[fi]) )
                continue;

            for (size_t fj = 0; fj < fcs.size(); fj++)
            {
                if ( msh.is_boundary(fcs[fj]) )
                    continue;
            
                auto dofsi = cbs*msh.cells_size() + fbs * compress_table.at( offset(msh, fcs[fi]) );
                auto dofsj = cbs*msh.cells_size() + fbs * compress_table.at( offset(msh, fcs[fj]) );

                for (size_t i = 0; i < fbs; i++)
                    for(size_t j = 0; j < fbs; j++)
                        trpK.push_back( Triplet<T>(dofsi+i, dofsj+j, K(cbs+fi*fbs+i, cbs+fj*fbs+j)) );

            }
        }  
    }

    void finalize(void)
    {
        gK.setFromTriplets(trpK.begin(), trpK.end());
        trpK.clear();
        gM.setFromTriplets(trpM.begin(), trpM.end());
        trpM.clear();

        std::cout << "Solving for " << syssz << " DOFs. ";
        std::cout << "K has " << gK.nonZeros() << " nonzeros, M has ";
        std::cout << gM.nonZeros() << " nonzeros." << std::endl; 
    }
};



#define USE_STATIC_CONDENSATION

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
computation_info<typename Mesh<T,2,Storage>::coordinate_type>
vector_wave_solver(Mesh<T,2,Storage>& msh, size_t order,
                   const typename Mesh<T,2,Storage>::coordinate_type& alpha,
                   const typename Mesh<T,2,Storage>::coordinate_type& omega)
{
    typedef Mesh<T,2,Storage>                   mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;
    
    disk::hho_degree_info chdi( { .rd = (size_t) order+1,
                                  .cd = (size_t) order+1,
                                  .fd = (size_t) order } );

    auto rhs_fun = [&](const point_type& pt) -> T {
        return omega * omega * std::sin(omega*pt.x())*std::sin(omega*pt.y());
    };

    auto curl_sol_fun = [&](const point_type& pt) -> Matrix<T,2,1> {
        Matrix<T,2,1> ret;
        ret(0) =   omega * std::sin(omega*pt.x())*std::cos(omega*pt.y());
        ret(1) = - omega * std::cos(omega*pt.x())*std::sin(omega*pt.y());
        return ret;
    };

    auto sol_fun = [&](const point_type& pt) -> T {
        return std::sin(omega*pt.x())*std::sin(omega*pt.y());
    };

#ifdef USE_STATIC_CONDENSATION
    std::vector<bool> is_dirichlet(msh.faces_size(), true);
    maxwell_assembler_condensed<mesh_type> assm(msh, chdi, is_dirichlet);
#else
    maxwell_assembler<mesh_type> assm(msh, chdi);
#endif

    std::cout << "Assembling to triplets" << std::endl;
    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, chdi.cell_degree());

        auto CR = disk::curl_reconstruction(msh, cl, chdi);
        auto ST = disk::curl_hdg_stabilization(msh, cl, chdi);

        Matrix<T, Dynamic, Dynamic> MM = Matrix<T, Dynamic, Dynamic>::Zero(ST.rows(), ST.cols());
        MM.block(0,0,cb.size(),cb.size()) = make_mass_matrix(msh, cl, cb);

        Matrix<T, Dynamic, Dynamic> lhs = CR.second + alpha*ST - (omega*omega)*MM;
        Matrix<T, Dynamic, 1> rhs = Matrix<T, Dynamic, 1>::Zero(lhs.rows());

        rhs.segment(0, cb.size()) = make_rhs(msh, cl, cb, rhs_fun, 1);


#ifdef USE_STATIC_CONDENSATION
        auto [LC, bC] = disk::static_condensation(lhs, rhs, cb.size());
        assm.assemble(msh, cl, LC, bC);
#else
        assm.assemble(msh, cl, lhs, rhs);
#endif
    }

    std::cout << "Triplets to matrix" << std::endl;
    assm.finalize();

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.syssz);

    std::cout << "Running pardiso" << std::endl;
    disk::solvers::pardiso_params<T> pparams;
    pparams.report_factorization_Mflops = true;
    mkl_pardiso(pparams, assm.LHS, assm.RHS, sol);


    std::vector<T> data_uz;

    T l2_err_e = 0.0;
    T l2_err_h = 0.0;
    T l2_err_Rh = 0.0;
    size_t cell_i = 0;
    T energy_err = 0.0;
    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, chdi.cell_degree());
        auto rb = make_scalar_monomial_basis(msh, cl, chdi.reconstruction_degree());

#ifdef USE_STATIC_CONDENSATION
        auto CR = disk::curl_reconstruction(msh, cl, chdi);
        auto ST = disk::curl_hdg_stabilization(msh, cl, chdi);
        
        Matrix<T, Dynamic, Dynamic> MM = Matrix<T, Dynamic, Dynamic>::Zero(ST.rows(), ST.cols());
        MM.block(0,0,cb.size(),cb.size()) = make_mass_matrix(msh, cl, cb);

        Matrix<T, Dynamic, Dynamic> lhs = CR.second + alpha*ST - (omega*omega)*MM;
        Matrix<T, Dynamic, 1> rhs = Matrix<T, Dynamic, 1>::Zero(lhs.rows());
        rhs.segment(0, cb.size()) = make_rhs(msh, cl, cb, rhs_fun, 1);
        auto edofs = assm.get_element_dofs(msh, cl, sol);
        auto esol = disk::static_decondensation(lhs, rhs, edofs);
        Matrix<T, Dynamic, 1> lsol = esol.segment(0, cb.size());

        Matrix<T, Dynamic, 1> asol = disk::project_function(msh, cl, cb, sol_fun);

        Matrix<T, Dynamic, 1> prj = project_tangent(msh, cl, chdi, sol_fun);

        Matrix<T, Dynamic, 1> hsol = Matrix<T, Dynamic, 1>::Zero(rb.size());
        
        hsol.segment(1, rb.size()-1) = CR.first * esol;

        energy_err += (prj-esol).dot(lhs*(prj-esol));
#else
        Matrix<T, Dynamic, Dynamic> MMe = disk::make_mass_matrix(msh, cl, cb);
        Matrix<T, Dynamic, 1> arhs = disk::make_rhs(msh, cl, cb, sol_fun);
        Matrix<T, Dynamic, 1> asol = MMe.llt().solve(arhs);
        Matrix<T, Dynamic, 1> lsol = sol.segment(cell_i*cb.size(), cb.size());
#endif


        auto qps = integrate(msh, cl, 2*chdi.reconstruction_degree());
        for (auto& qp : qps)
        {
            Matrix<T, Dynamic, 2> hphi  = rb.eval_curls2(qp.point());
            Matrix<T, Dynamic, 2> hphi2 = hphi.block(0,0,cb.size(),2); 
            Matrix<T, Dynamic, 1> ephi  = cb.eval_functions(qp.point());
            //Matrix<T, Dynamic, 3> rphi = rb.eval_functions(qp.point());

            Matrix<T,Dynamic,1> esolseg = esol.segment(0, cb.size());
            auto hdiff = disk::eval(esolseg, hphi2) - curl_sol_fun(qp.point());
            auto Rhdiff = disk::eval(hsol, hphi) - curl_sol_fun(qp.point());
            auto ediff = disk::eval(esolseg, ephi) - sol_fun(qp.point());

            l2_err_h += qp.weight() * hdiff.dot(hdiff);
            l2_err_Rh += qp.weight() * Rhdiff.dot(Rhdiff);
            l2_err_e += qp.weight() * ediff*ediff;
        }

        cell_i++;
    }

    std::cout << "h = " << disk::average_diameter(msh) << " ";
    std::cout << "√(Σ‖e - eₜ‖²) = " << std::sqrt(l2_err_e) << ", ";
    std::cout << "√(Σ‖∇×(e - eₜ)‖²) = " << std::sqrt(l2_err_h) << ", ";
    std::cout << "‖∇×(u - C(uₕ))‖ = " << std::sqrt(l2_err_Rh) << ", ";
    std::cout << "aₕ(I(u)-uₕ,I(u)-uₕ) = " << std::sqrt(energy_err) << std::endl;

    /*
    disk::silo_database silo_db;
    silo_db.create("maxwell.silo");
    silo_db.add_mesh(msh, "mesh");

    disk::silo_zonal_variable<T> uz("uz", data_uz);
    silo_db.add_variable("mesh", uz);

    silo_db.close();
 */

    return computation_info<T>({
            .l2_error_e = std::sqrt(l2_err_e),
            .l2_error_h = std::sqrt(l2_err_h),
            .nrg_error = 0,
            .mflops = pparams.mflops,
            .dofs = assm.syssz,
            .nonzeros = assm.LHS.nonZeros()
        });
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
computation_info<typename Mesh<T,3,Storage>::coordinate_type>
vector_wave_solver(Mesh<T,3,Storage>& msh, size_t order,
                   const typename Mesh<T,3,Storage>::coordinate_type& alpha,
                   const typename Mesh<T,3,Storage>::coordinate_type& omega,
                   bool ho_stab)
{
    typedef Mesh<T,3,Storage>                   mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;
    
    disk::hho_degree_info chdi( { .rd = (size_t) order+1,
                                  .cd = (size_t) order,
                                  .fd = (size_t) order } );

    auto rhs_fun = [&](const point_type& pt) -> Matrix<T, 3, 1> {
        Matrix<T, 3, 1> ret;
        ret(0) = 0.0;
        ret(1) = 0.0;
        ret(2) = omega*omega * std::sin(omega*pt.x())*std::sin(omega*pt.y());

        //ret(0) = M_PI * M_PI * std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y())*std::sin(M_PI*pt.z());
        //ret(1) = M_PI * M_PI * std::cos(M_PI*pt.x())*std::cos(M_PI*pt.y())*std::sin(M_PI*pt.z());
        //ret(2) = M_PI * M_PI * std::cos(M_PI*pt.x())*std::sin(M_PI*pt.y())*std::cos(M_PI*pt.z());

        return ret;
    };

    auto sol_fun = [&](const point_type& pt) -> Matrix<T, 3, 1> {
        Matrix<T, 3, 1> ret;
        ret(0) = 0.0;
        ret(1) = 0.0;
        ret(2) = std::sin(omega*pt.x())*std::sin(omega*pt.y());

        //ret(0) = std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y())*std::sin(M_PI*pt.z());
        //ret(1) = 0.0;
        //ret(2) = 0.0;

        return ret;
    };

    auto curl_sol_fun = [&](const point_type& pt) -> Matrix<T, 3, 1> {
        Matrix<T, 3, 1> ret;
        ret(0) =  omega*std::sin(omega*pt.x())*std::cos(omega*pt.y());
        ret(1) = -omega*std::cos(omega*pt.x())*std::sin(omega*pt.y());
        ret(2) = 0.0;

        //ret(0) =   0.0;
        //ret(1) =   M_PI * std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y())*std::cos(M_PI*pt.z());
        //ret(2) = - M_PI * std::sin(M_PI*pt.x())*std::cos(M_PI*pt.y())*std::sin(M_PI*pt.z());

        return ret;
    };

#ifdef USE_STATIC_CONDENSATION
    std::vector<bool> is_dirichlet(msh.faces_size(), true);
    maxwell_assembler_condensed<mesh_type> assm(msh, chdi, is_dirichlet);
#else
    maxwell_assembler<mesh_type> assm(msh, chdi);
#endif

    std::cout << "Assembling to triplets" << std::endl;
    for (auto& cl : msh)
    {
        auto CR = disk::curl_reconstruction(msh, cl, chdi);

        Matrix<T, Dynamic, Dynamic> ST;
        //if (ho_stab)
        //    ST = disk::curl_hho_stabilization(msh, cl, CR.first, chdi);
        //else
            ST = disk::curl_hdg_stabilization(msh, cl, chdi);

        auto MM = make_vector_mass_oper(msh, cl, chdi);

        Matrix<T, Dynamic, Dynamic> lhs = CR.second + alpha*ST - (omega*omega)*MM;
        Matrix<T, Dynamic, 1> rhs = Matrix<T, Dynamic, 1>::Zero(lhs.rows());

        auto cb = make_vector_monomial_basis(msh, cl, chdi.cell_degree());
        rhs.segment(0, cb.size()) = make_rhs(msh, cl, cb, rhs_fun, 1);

#ifdef USE_STATIC_CONDENSATION
        auto [LC, bC] = disk::static_condensation(lhs, rhs, cb.size());
        assm.assemble(msh, cl, LC, bC);
#else
        assm.assemble(msh, cl, lhs, rhs);
#endif
    }

    std::cout << "Triplets to matrix" << std::endl;
    assm.finalize();

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.syssz);

    ///*
    std::cout << "Running pardiso" << std::endl;
    disk::solvers::pardiso_params<T> pparams;
    pparams.report_factorization_Mflops = true;
    pparams.out_of_core = PARDISO_OUT_OF_CORE_IF_NEEDED;
    mkl_pardiso(pparams, assm.LHS, assm.RHS, sol);
    //*/

    /*
    std::cout << "Running MUMPS" << std::endl;
    mumps_solver<T> mumps;
    sol = mumps.solve(assm.LHS, assm.RHS);
    */

    std::vector<T> data_ex, data_ey, data_ez;
    std::vector<T> data_hx, data_hy, data_hz;

    T l2_err_h = 0.0;
    T l2_err_Rh = 0.0;
    T l2_err_e = 0.0;
    T energy_err = 0.0;
    size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto cb = make_vector_monomial_basis(msh, cl, chdi.cell_degree());
        auto rb = make_vector_monomial_basis(msh, cl, chdi.reconstruction_degree());

#ifdef USE_STATIC_CONDENSATION
        auto CR = disk::curl_reconstruction(msh, cl, chdi);

        Matrix<T, Dynamic, Dynamic> ST;
        //if (ho_stab)
        //    ST = disk::curl_hho_stabilization(msh, cl, CR.first, chdi);
        //else
            ST = disk::curl_hdg_stabilization(msh, cl, chdi);
        
        auto MM = make_vector_mass_oper(msh, cl, chdi);
        Matrix<T, Dynamic, Dynamic> lhs = CR.second + alpha*ST - (omega*omega)*MM;
        Matrix<T, Dynamic, 1> rhs = Matrix<T, Dynamic, 1>::Zero(lhs.rows());
        rhs.segment(0, cb.size()) = make_rhs(msh, cl, cb, rhs_fun);
        auto edofs = assm.get_element_dofs(msh, cl, sol);

        Matrix<T, Dynamic, 1> esol = disk::static_decondensation(lhs, rhs, edofs);
        Matrix<T, Dynamic, 1> asol = project_tangent(msh, cl, chdi, sol_fun, 1);
        Matrix<T, Dynamic, 1> hsol = Matrix<T, Dynamic, 1>::Zero(rb.size());

        hsol.segment(3, rb.size()-3) = CR.first * esol;

        auto qps = integrate(msh, cl, 2*chdi.reconstruction_degree());
        for (auto& qp : qps)
        {
            Matrix<T, Dynamic, 3> hphi = rb.eval_curls2(qp.point());
            Matrix<T, Dynamic, 3> hphi2 = cb.eval_curls2(qp.point());
            Matrix<T, Dynamic, 3> ephi = cb.eval_functions(qp.point());
            Matrix<T, Dynamic, 3> rphi = rb.eval_functions(qp.point());

            Matrix<T,Dynamic,1> esolseg = esol.segment(0, cb.size());
            auto hdiff = disk::eval(esolseg, hphi2) - curl_sol_fun(qp.point());
            auto Rhdiff = disk::eval(hsol, hphi) - curl_sol_fun(qp.point());
            auto ediff = disk::eval(esolseg, ephi) - sol_fun(qp.point());

            l2_err_h += qp.weight() * hdiff.dot(hdiff);
            l2_err_Rh += qp.weight() * Rhdiff.dot(Rhdiff);
            l2_err_e += qp.weight() * ediff.dot(ediff);
        }

        Matrix<T, Dynamic, Dynamic> lhs2 = CR.second + alpha*ST - (omega*omega)*MM;
        energy_err += (esol-asol).dot(lhs2*(esol-asol));

#else
        Matrix<T, Dynamic, Dynamic> MMe = disk::make_mass_matrix(msh, cl, cb);
        Matrix<T, Dynamic, 1> arhs = disk::make_rhs(msh, cl, cb, sol_fun);
        Matrix<T, Dynamic, 1> asol = MMe.llt().solve(arhs);
        Matrix<T, Dynamic, 1> lsol = sol.segment(cell_i*cb.size(), cb.size());
#endif

        auto bar = barycenter(msh, cl);
        auto phi = cb.eval_functions(bar);
        auto cphi = cb.eval_curls2(bar);
        Matrix<T,Dynamic,1> esolseg = esol.segment(0, cb.size());
        auto locsol_e = disk::eval( esolseg, phi );
        auto locsol_h = disk::eval( esolseg, cphi );


        data_ex.push_back( locsol_e(0) );
        data_ey.push_back( locsol_e(1) );
        data_ez.push_back( locsol_e(2) );
        data_hx.push_back( locsol_h(0) );
        data_hy.push_back( locsol_h(1) );
        data_hz.push_back( locsol_h(2) );

        cell_i++;
    }

    std::cout << "h = " << disk::average_diameter(msh) << " ";
    std::cout << "√(Σ‖e - eₜ‖²) = " << std::sqrt(l2_err_e) << ", ";
    std::cout << "√(Σ‖∇×(e - eₜ)‖²) = " << std::sqrt(l2_err_h) << ", ";
    std::cout << "‖∇×(u - C(uₕ))‖ = " << std::sqrt(l2_err_Rh) << ", ";
    std::cout << "aₕ(I(u)-uₕ,I(u)-uₕ) = " << std::sqrt(energy_err) << std::endl;
    disk::silo_database silo_db;
    silo_db.create("maxwell.silo");
    silo_db.add_mesh(msh, "mesh");

    disk::silo_zonal_variable<T> ex("ex", data_ex);
    silo_db.add_variable("mesh", ex);

    disk::silo_zonal_variable<T> ey("ey", data_ey);
    silo_db.add_variable("mesh", ey);

    disk::silo_zonal_variable<T> ez("ez", data_ez);
    silo_db.add_variable("mesh", ez);

    silo_db.add_expression("e", "{ex, ey, ez}", DB_VARTYPE_VECTOR);

    disk::silo_zonal_variable<T> hx("hx", data_hx);
    silo_db.add_variable("mesh", hx);

    disk::silo_zonal_variable<T> hy("hy", data_hy);
    silo_db.add_variable("mesh", hy);

    disk::silo_zonal_variable<T> hz("hz", data_hz);
    silo_db.add_variable("mesh", hz);

    silo_db.add_expression("h", "{hx, hy, hz}", DB_VARTYPE_VECTOR);

    silo_db.close();

    return computation_info<T>({
            .l2_error_e = std::sqrt(l2_err_e),
            .l2_error_h = std::sqrt(l2_err_h),
            .nrg_error = std::sqrt(energy_err),
            .mflops = pparams.mflops,
            //.mflops = (int)mumps.get_Mflops(),
            .dofs = assm.syssz,
            .nonzeros = assm.LHS.nonZeros()
        });
}


template<template<typename, size_t, typename> class Mesh, typename CoordT, typename Storage>
void
vector_wave_solver_complex(Mesh<CoordT,3,Storage>& msh, size_t order)
{
    typedef Mesh<CoordT,3,Storage>              mesh_type;
    typedef std::complex<double>                scalar_type;
    //typedef double                              scalar_type;
    typedef typename mesh_type::point_type      point_type;
    typedef typename mesh_type::cell_type       cell_type;
    typedef typename mesh_type::face_type       face_type;

    double omega = 2*M_PI*1.55e10;
    double alpha = 5;
    
    disk::hho_degree_info chdi( { .rd = (size_t) order,
                                  .cd = (size_t) order,
                                  .fd = (size_t) order } );

    auto rhs_fun = [&](const point_type& pt) -> Matrix<scalar_type, 3, 1> {
        Matrix<scalar_type, 3, 1> ret;
        ret(0) = 0.0;
        ret(1) = 0.0;
        ret(2) = 0.0;//omega*omega * std::sin(omega*pt.x())*std::sin(omega*pt.y());
        return ret;
    };

    size_t air_tag = 13;

    size_t port_tag = 1188;
    size_t port_id;

    size_t imp_tag = 1193;
    size_t imp_id;

    std::vector<bool> is_dirichlet(msh.faces_size(), false);

    for (auto& fc : faces(msh))
    {
        auto bi = msh.boundary_info(fc);
        if (not bi.is_boundary())
            continue;

        if (bi.is_internal())
            continue;

        if (bi.tag() != port_tag and bi.tag() != imp_tag)
            is_dirichlet.at( offset(msh, fc) ) = true;
        else
        {
            if (bi.tag() == port_tag)
                port_id = bi.id();
            else
                imp_id = bi.id();
        }
    }

    maxwell_assembler_condensed<mesh_type, scalar_type> assm(msh, chdi, is_dirichlet);

    auto eps0 = 8.8541878128e-12;
    auto mu0 = 4*M_PI*1e-7;

    auto compute_local_contribution = [&](const cell_type& cl) -> auto {
        auto di = msh.domain_info(cl);
        double eps_r;
        if (di.tag() == air_tag)
            eps_r = 1;
        else
            eps_r = 24;

        auto eps = eps_r * eps0;
        auto mu = mu0;
        auto Z = std::sqrt(mu/eps);
        auto CR = disk::curl_reconstruction_pk(msh, cl, chdi);
        auto ST = disk::curl_hdg_stabilization(msh, cl, chdi);
        auto MM = disk::make_vector_mass_oper(msh, cl, chdi);
        auto [ZZ, zr] = disk::make_impedance_term<scalar_type>(msh, cl, chdi, port_id, 1./Z);
        auto [ZB, zb] = disk::make_impedance_term<scalar_type>(msh, cl, chdi, imp_id, 1./Z);

        Matrix<scalar_type, Dynamic, Dynamic> lhs =
            (1./mu) * CR.second +
            alpha * ST -
            (eps*omega)*omega*MM -
            std::complex<double>(0,omega)*(ZZ+ZB);

        Matrix<scalar_type, Dynamic, 1> rhs = Matrix<scalar_type, Dynamic, 1>::Zero(lhs.rows());
        rhs += std::complex<double>(0,2*omega)*zr;
        return std::make_pair(lhs, rhs);
    };

    std::cout << "Assembling to triplets" << std::endl;
    for (auto& cl : msh)
    {
        auto cbs = disk::vector_basis_size(chdi.cell_degree(), 3, 3);
        auto [lhs, rhs] = compute_local_contribution(cl);
        auto [LC, bC] = disk::static_condensation(lhs, rhs, cbs);
        assm.assemble(msh, cl, LC, bC);
    }

    std::cout << "Triplets to matrix" << std::endl;
    assm.finalize();

    disk::dynamic_vector<scalar_type> sol = disk::dynamic_vector<scalar_type>::Zero(assm.syssz);

    ///*
    std::cout << "Running pardiso" << std::endl;
    disk::solvers::pardiso_params<scalar_type> pparams;
    pparams.report_factorization_Mflops = true;
    pparams.out_of_core = PARDISO_OUT_OF_CORE_IF_NEEDED;
    mkl_pardiso(pparams, assm.LHS, assm.RHS, sol);
    //*/

    /*
    std::cout << "Running MUMPS" << std::endl;
    mumps_solver<T> mumps;
    sol = mumps.solve(assm.LHS, assm.RHS);
    */

    std::vector<scalar_type> data_ex, data_ey, data_ez;
    std::vector<std::pair<scalar_type, int>> nodal_ez;
    nodal_ez.resize(msh.points_size(), std::make_pair(0.0, 0));

    std::vector<std::pair<scalar_type, int>> nodal_ez_fc;
    nodal_ez_fc.resize(msh.points_size(), std::make_pair(0.0, 0));

    scalar_type err = 0.0;
    size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto [lhs, rhs] = compute_local_contribution(cl);

        auto cb = disk::make_vector_monomial_basis<mesh_type, cell_type, scalar_type>(msh, cl, chdi.cell_degree());

        auto edofs = assm.get_element_dofs(msh, cl, sol);

        Matrix<scalar_type, Dynamic, 1> esol = disk::static_decondensation(lhs, rhs, edofs);

        auto bar = barycenter(msh, cl);
        Matrix<scalar_type, Dynamic, 3> phi = cb.eval_functions(bar);
        Matrix<scalar_type,Dynamic,1> esolseg = esol.segment(0, cb.size());
        auto locsol_e = disk::eval( esolseg, phi );

        data_ex.push_back( locsol_e(0) );
        data_ey.push_back( locsol_e(1) );
        data_ez.push_back( locsol_e(2) );


        auto fcs = faces(msh, cl);
        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto& fc = fcs[i];
            if (not msh.is_boundary(fc))
                continue;

            auto fb = make_vector_monomial_tangential_basis(msh, fc, chdi.face_degree());
            auto fbs = fb.size();
            auto face_dofs = edofs.segment(fbs*i, fbs);
            auto pts = points(msh, fc);
            auto ptids = fc.point_ids();

            for (size_t i = 0; i < pts.size(); i++)
            {
                auto fphi = fb.eval_functions(pts[i]);
                auto fls = fphi.transpose()*face_dofs;
                auto cphi = cb.eval_functions(pts[i]);
                auto cls = cphi.transpose()*esolseg;
                auto ptid = ptids.at(i);
                nodal_ez_fc.at(ptid).first += std::abs(fls.norm() - cls.norm());
                nodal_ez_fc.at(ptid).second++;
            }
        }

        cell_i++;
    }

    std::vector<scalar_type> nez;
    for (auto& ez : nodal_ez_fc)
        if (ez.second != 0)
            nez.push_back( ez.first/double(ez.second) );
        else
            nez.push_back(0);


    /*
    std::ofstream ofs("facesol.dat");
    size_t face_i = 0;
    for (auto& fc : faces(msh))
    {
        if ( msh.is_boundary(fc) and is_dirichlet.at(face_i) )
            continue;

        auto fb = disk::make_vector_monomial_tangential_basis<mesh_type, face_type, scalar_type>(msh, fc, chdi.face_degree());
        auto fbs = fb.size();
        auto bar = barycenter(msh, fc);

        Matrix<scalar_type, Dynamic, 3> phi = fb.eval_functions(bar);
        Matrix<scalar_type,Dynamic,1> solseg = sol.segment(face_i*fbs, fbs);
        auto locsol = disk::eval(solseg, phi);

        auto xmag = std::abs(locsol(0));
        auto ymag = std::abs(locsol(1));
        auto zmag = std::abs(locsol(2));
        auto vmag = std::sqrt( xmag*xmag + ymag*ymag + zmag*zmag );

        ofs << bar.x() << " " << bar.y() << " " << bar.z() << " " << vmag << std::endl;

        face_i++;
    }
    */

    disk::silo_database silo_db;
    silo_db.create("maxwell.silo");
    silo_db.add_mesh(msh, "mesh");

    disk::silo_zonal_variable<scalar_type> ex("ex", data_ex);
    silo_db.add_variable("mesh", ex);

    disk::silo_zonal_variable<scalar_type> ey("ey", data_ey);
    silo_db.add_variable("mesh", ey);

    disk::silo_zonal_variable<scalar_type> ez("ez", data_ez);
    silo_db.add_variable("mesh", ez);

    disk::silo_nodal_variable<scalar_type> ezn("ez_nodal", nez);
    silo_db.add_variable("mesh", ezn);

    silo_db.add_expression("e", "{ex, ey, ez}", DB_VARTYPE_VECTOR);

    silo_db.close();
}



template<typename Mesh>
void maxwell_eigenvalue_solver(Mesh& msh, size_t order, const typename Mesh::coordinate_type& alpha)
{
    using T = typename Mesh::coordinate_type;

    disk::hho_degree_info chdi( { .rd = (size_t) order+1,
                                  .cd = (size_t) order,
                                  .fd = (size_t) order } );


    maxwell_eigenvalue_assembler<Mesh> assm(msh, chdi);

    std::cout << "Assembling to triplets" << std::endl;
    for (auto& cl : msh)
    {
        auto CR = disk::curl_reconstruction<double>(msh, cl, chdi);
        auto ST = disk::curl_hdg_stabilization(msh, cl, chdi);
        auto MM = make_vector_mass_oper(msh, cl, chdi);
        assm.assemble(msh, cl, CR.second+alpha*ST, MM);
    }

    assm.finalize();

    disk::feast_eigensolver_params<T> fep;

    fep.verbose = true;
    fep.tolerance = 8;
    fep.min_eigval = 1;
    fep.max_eigval = 100;
    fep.subspace_size = 50;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    hho_eigvecs;
    Eigen::Matrix<T, Eigen::Dynamic, 1>                 hho_eigvals;

    //disk::dump_sparse_matrix(assm.gK, "gK.txt");
    //disk::dump_sparse_matrix(assm.gM, "gM.txt");

    generalized_eigenvalue_solver(fep, assm.gK, assm.gM, hho_eigvecs, hho_eigvals);


    for (size_t i = 0; i < fep.eigvals_found; i++)
    {
        std::vector<T> data_ux, data_uy, data_uz;

        std::stringstream fn;
        fn << "eigenfun_" << i << ".silo";

        disk::silo_database silo_db;
        silo_db.create(fn.str().c_str());
        silo_db.add_mesh(msh, "mesh");

        size_t cell_i = 0;
        for (auto& cl : msh)
        {
            auto cb = make_vector_monomial_basis(msh, cl, chdi.cell_degree());
            Matrix<T, Dynamic, 1> lsol = hho_eigvecs.block(cell_i*cb.size(), i, cb.size(), 1);

            Matrix<T,1,3> acc = Matrix<T,1,3>::Zero();
            auto bar = barycenter(msh, cl);
            for (size_t i = 0; i < cb.size(); i++)
            {
                auto phi = cb.eval_functions(bar);
                acc += lsol(i)*phi.row(i);
            }

            data_ux.push_back( acc(0) );
            data_uy.push_back( acc(1) );
            data_uz.push_back( acc(2) );

            cell_i++;
        }

        disk::silo_zonal_variable<T> ux("ux", data_ux);
        silo_db.add_variable("mesh", ux);

        disk::silo_zonal_variable<T> uy("uy", data_uy);
        silo_db.add_variable("mesh", uy);

        disk::silo_zonal_variable<T> uz("uz", data_uz);
        silo_db.add_variable("mesh", uz);
    
        silo_db.add_expression("u", "{ux, uy, uz}", DB_VARTYPE_VECTOR);
        silo_db.add_expression("mag_u", "magnitude(u)", DB_VARTYPE_SCALAR);
    }

    for (size_t i = 0; i < fep.eigvals_found; i++)
    {
        std::cout << std::setprecision(8) << hho_eigvals(i) << " -> ";
        std::cout << std::sqrt( hho_eigvals(i) ) << std::endl;
    }
    
}


void autotest_alpha(size_t order)
{
    using T = double;
    
    using Mesh = disk::simplicial_mesh<T,3>;

    std::stringstream ssl2;
    ssl2 << "l2log_k" << order << ".txt";

    std::stringstream ssnrg;
    ssnrg << "nrglog_k" << order << ".txt";

    std::ofstream ofs_l2( ssl2.str() );
    std::ofstream ofs_nrg( ssnrg.str() );

    for (T sp = 1.0; sp < 16.0; sp += 1.0)
    {
        ofs_l2  << sp << " ";
        ofs_nrg << sp << " ";
        
        for (size_t i = 1; i < 5; i++)
        {
            Mesh msh;
            std::stringstream ss;
            ss << "../../../meshes/3D_tetras/netgen/cube" << i << ".mesh";

            disk::load_mesh_netgen(ss.str().c_str(), msh);

            auto diam = disk::average_diameter(msh);

            auto ci = vector_wave_solver(msh, order, sp, M_PI, false);  
            ofs_l2  << ci.l2_error_e << " " << ci.l2_error_h << " ";
            ofs_nrg << ci.nrg_error << " ";
        }

        ofs_l2 << std::endl;
        ofs_nrg << std::endl;
    }
}

void autotest_convergence(size_t order_min, size_t order_max)
{
    using T = double;
    
    using Mesh = disk::simplicial_mesh<T,3>;

    std::ofstream ofs( "hho_convergence_convt.txt" );

    for (size_t i = 1; i < 5; i++)
    {
        Mesh msh;
        std::stringstream ss;
        ss << "../../../meshes/3D_tetras/netgen/convt0" << i << ".mesh";

        disk::load_mesh_netgen(ss.str().c_str(), msh);
        auto diam = disk::average_diameter(msh);

        ofs << diam << " ";

        for (size_t order = order_min; order <= order_max; order++)
        {
            if (order >= 2 and i > 4)
                break;

            auto ci = vector_wave_solver(msh, order, 1.0, M_PI, false);
            ofs << ci.l2_error_e << " " << ci.l2_error_h << " " << ci.nrg_error << " " << ci.mflops << " ";
            ofs << ci.dofs << " " << ci.nonzeros << " ";
        }

        ofs << std::endl;
    }
}

void autotest(size_t order)
{
    //for (size_t i = 0; i <= order; i++)
    //    autotest_alpha(i);

    autotest_convergence(0, order);
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
    bool        use_ho_stab = false;

    int ch;
    while ( (ch = getopt(argc, argv, "Aa:ek:m:w:H")) != -1 )
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

            case 'A':
                autotest(degree);
                return 0;

            case 'w':
                omega = M_PI*std::stod(optarg);
                break;

            case 'H':
                use_ho_stab = true;
                break;

            case '?':
            default:
                std::cout << "Invalid option" << std::endl;
                return 1;
        }
    }

#if 0
    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<coord_T>(mesh_filename);

        //if (solve_eigvals)
        //    maxwell_eigenvalue_solver(msh, degree, stab_param);
        //else
            vector_wave_solver(msh, degree, stab_param, omega);

        return 0;
    }

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        auto msh = disk::load_fvca5_2d_mesh<coord_T>(mesh_filename);

        //if (solve_eigvals)
        //    maxwell_eigenvalue_solver(msh, degree, stab_param);
        //else
            vector_wave_solver(msh, degree, stab_param, omega);

        return 0;
    }
#endif

    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        auto msh = disk::load_netgen_3d_mesh<coord_T>(mesh_filename);

        //if (solve_eigvals)
        //    maxwell_eigenvalue_solver(msh, degree, stab_param);
        //else
            //vector_wave_solver(msh, degree, stab_param, omega, false);
            vector_wave_solver_complex(msh, degree);

        return 0;
    }

    /* GMSH */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo$") ))
    {
        std::cout << "Guessed mesh format: GMSH" << std::endl;
        disk::generic_mesh<T,3> msh;
        disk::gmsh_geometry_loader< disk::generic_mesh<T,3> > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        vector_wave_solver_complex(msh, degree);

        return 0;
    }

#if 0
    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        auto msh = disk::load_cartesian_3d_mesh<coord_T>(mesh_filename);
        
        if (solve_eigvals)
            maxwell_eigenvalue_solver(msh, degree, stab_param);
        else
            //vector_wave_solver(msh, degree, stab_param, omega);
            vector_wave_solver_complex(msh, degree);
        
        return 0;
    }

    /* FVCA6 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
        disk::generic_mesh<coord_T,3> msh;
        
        disk::load_mesh_fvca6_3d<coord_T>(mesh_filename, msh);
        
        if (solve_eigvals)
            maxwell_eigenvalue_solver(msh, degree, stab_param);
        else
            //vector_wave_solver(msh, degree, stab_param, omega);
            vector_wave_solver_complex(msh, degree);
        
        return 0;
    }
#endif
}

