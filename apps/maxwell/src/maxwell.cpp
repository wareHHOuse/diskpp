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

#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/hho"
#include "diskpp/methods/implementation_hho/curl.hpp"
#include "diskpp/output/silo.hpp"

#include "mumps.hpp"

#include "compinfo.h"

//#include <Eigen/IterativeLinearSolvers>
//#include <unsupported/Eigen/IterativeSolvers>

#include "paramloader.hpp"

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
             const Matrix<ScalT, Dynamic, 1>& rhs,
             const Matrix<ScalT, Dynamic, 1>& dirichlet_data)
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
                auto lofsi = fi*fbs;
                auto lofsj = fj*fbs;

                if ( not is_in_system(msh, fcs[fj]) )
                {
                    RHS.segment(cofsi, fbs) += -lhsc.block(lofsi, lofsj, fbs, fbs)*dirichlet_data.segment(lofsj, fbs);
                    continue;
                }

                auto cofsj = get_system_offset(msh, fcs[fj]);
                for (size_t i = 0; i < fbs; i++)
                    for(size_t j = 0; j < fbs; j++)
                        triplets.push_back( Triplet<ScalT>(cofsi+i, cofsj+j, lhsc(lofsi+i, lofsj+j)) );
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

#if 0
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
computation_info<typename Mesh<T,2,Storage>::coordinate_type>
vector_wave_solver(Mesh<T,2,Storage>& msh, size_t order,
                   const typename Mesh<T,2,Storage>::coordinate_type& alpha,
                   const typename Mesh<T,2,Storage>::coordinate_type& omega)
{
    typedef Mesh<T,2,Storage>                   mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;

    disk::hho_degree_info chdi( { .rd = (size_t) order,
                                  .cd = (size_t) order,
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

    std::cout << "Running MUMPS" << std::endl;
    sol = mumps_lu(assm.LHS, assm.RHS);


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
            Matrix<T, Dynamic, 2> hphi  = rb.eval_curls(qp.point());
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
vector_wave_solver(Mesh<T,3,Storage>& msh, size_t order)
{
    typedef Mesh<T,3,Storage>                   mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;

    disk::hho_degree_info chdi( { .rd = (size_t) order,
                                  .cd = (size_t) order+1,
                                  .fd = (size_t) order } );

    T omega = std::sqrt(5);

    auto rhs_fun = [&](const point_type& pt) -> Matrix<T, 3, 1> {
        Matrix<T, 3, 1> ret;
        ret(0) = 0.0;
        ret(1) = 0.0;
        //ret(2) = (2*M_PI*M_PI - omega*omega) * std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
        ret(2) = - pt.y()*std::sin(M_PI*pt.y()) * (2*M_PI*std::cos(M_PI*pt.x()) - pt.x()*M_PI*M_PI*std::sin(M_PI*pt.x()) )
                 - pt.x()*std::sin(M_PI*pt.x()) * (2*M_PI*std::cos(M_PI*pt.y()) - pt.y()*M_PI*M_PI*std::sin(M_PI*pt.y()) )
                 - omega*omega*pt.x() * pt.y() * std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
        return ret;
    };

    auto sol_fun = [&](const point_type& pt) -> Matrix<T, 3, 1> {
        Matrix<T, 3, 1> ret;
        ret(0) = 0.0;
        ret(1) = 0.0;
        ret(2) = pt.x() * pt.y() * std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
        return ret;
    };

    auto curl_sol_fun = [&](const point_type& pt) -> Matrix<T, 3, 1> {
        Matrix<T, 3, 1> ret;
        ret(0) =  pt.x()*std::sin(M_PI*pt.x())*( std::sin(M_PI*pt.y()) + M_PI*pt.y()*std::cos(M_PI*pt.y()) );
        ret(1) = -pt.y()*std::sin(M_PI*pt.y())*( std::sin(M_PI*pt.x()) + M_PI*pt.x()*std::cos(M_PI*pt.x()) );
        ret(2) =  0.0;
        return ret;
    };

    std::vector<bool> is_dirichlet(msh.faces_size(), true);
    maxwell_assembler_condensed<mesh_type> assm(msh, chdi, is_dirichlet);

    std::cout << "Assembling to triplets" << std::endl;
    for (auto& cl : msh)
    {
        auto CR = disk::curl_reconstruction_pk(msh, cl, chdi);

        Matrix<T, Dynamic, Dynamic> ST = disk::curl_hdg_stabilization(msh, cl, chdi);
        auto MM = make_vector_mass_oper(msh, cl, chdi);

        Matrix<T, Dynamic, Dynamic> lhs = CR.second + omega*(ST - omega*MM);
        Matrix<T, Dynamic, 1> rhs = Matrix<T, Dynamic, 1>::Zero(lhs.rows());

        auto cb = make_vector_monomial_basis(msh, cl, chdi.cell_degree());
        rhs.segment(0, cb.size()) = make_rhs(msh, cl, cb, rhs_fun, 1);

        auto [LC, bC] = disk::static_condensation(lhs, rhs, cb.size());

        disk::dynamic_vector<T> dd = disk::dynamic_vector<T>::Zero(lhs.rows());

        assm.assemble(msh, cl, LC, bC, dd);
    }

    std::cout << "Triplets to matrix" << std::endl;
    assm.finalize();

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.syssz);

    std::cout << "Running MUMPS" << std::endl;
    sol = mumps_lu(assm.LHS, assm.RHS);

    std::vector<T> data_ex, data_ey, data_ez;
    std::vector<T> data_hx, data_hy, data_hz;

    T l2_err_e = 0.0;
    T l2_err_R = 0.0;
    T mass_err = 0.0;
    T curl_err = 0.0;
    size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto cb = make_vector_monomial_basis(msh, cl, chdi.cell_degree());
        auto rb = make_vector_monomial_basis(msh, cl, chdi.reconstruction_degree());

        auto CR = disk::curl_reconstruction_pk(msh, cl, chdi);

        Matrix<T, Dynamic, Dynamic> ST = disk::curl_hdg_stabilization(msh, cl, chdi);

        auto MM = make_vector_mass_oper(msh, cl, chdi);
        Matrix<T, Dynamic, Dynamic> lhs = CR.second + omega*(ST - omega*MM);
        Matrix<T, Dynamic, 1> rhs = Matrix<T, Dynamic, 1>::Zero(lhs.rows());
        rhs.segment(0, cb.size()) = make_rhs(msh, cl, cb, rhs_fun, 1);
        auto edofs = assm.get_element_dofs(msh, cl, sol);

        Matrix<T, Dynamic, 1> esol = disk::static_decondensation(lhs, rhs, edofs);
        Matrix<T,Dynamic,1> esolseg = esol.segment(0, cb.size());
        Matrix<T,Dynamic,1> rsol = CR.first*esol;

        Matrix<T,Dynamic,1> asol = disk::project_tangent(msh, cl, chdi, sol_fun, 1);

        auto qps = integrate(msh, cl, 2*chdi.cell_degree()+2);
        for (auto& qp : qps)
        {
            Matrix<T, Dynamic, 3> ephi = cb.eval_functions(qp.point());
            disk::dynamic_vector<T> ediff = disk::eval(esolseg, ephi) - sol_fun(qp.point());
            l2_err_e += qp.weight() * ediff.dot(ediff);

            Matrix<T, Dynamic, 3> rphi = rb.eval_functions(qp.point());
            disk::dynamic_vector<T> rdiff = disk::eval(rsol, rphi) - curl_sol_fun(qp.point());
            l2_err_R += qp.weight() * rdiff.dot(rdiff);
        }

        Matrix<T, Dynamic, Dynamic> CC = CR.second + omega*ST;
        Matrix<T,Dynamic,1> diff = esol - asol;
        mass_err += diff.dot(MM*diff);
        curl_err += diff.dot(CC*diff);

        cell_i++;
    }

    std::cout << "h = " << disk::average_diameter(msh) << " ";
    std::cout << "√(Σ‖e - eₜ‖²) = " << std::sqrt(l2_err_e) << ", ";
    //std::cout << "√(Σ‖∇×(e - eₜ)‖²) = " << std::sqrt(l2_err_h) << ", ";
    std::cout << "‖∇×(u - C(uₕ))‖ = " << std::sqrt(l2_err_R) << std::endl;// << ", ";
    //std::cout << "aₕ(I(u)-uₕ,I(u)-uₕ) = " << std::sqrt(energy_err) << std::endl;
    std::cout << "mass_err = " << std::sqrt(mass_err) << " curl_err = " << std::sqrt(curl_err) << std::endl;

    /*
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

    */

    return computation_info<T>({
            .l2_error_e = std::sqrt(l2_err_e),
            .l2_error_h = std::sqrt(l2_err_R),
            .nrg_error = 0.0,//std::sqrt(energy_err),
            .mflops = pparams.mflops,
            //.mflops = (int)mumps.get_Mflops(),
            .dofs = assm.syssz,
            .nonzeros = assm.LHS.nonZeros()
        });
}
#endif

template<template<typename, size_t, typename> class Mesh, typename CoordT, typename Storage>
void
vector_wave_solver_complex(Mesh<CoordT,3,Storage>& msh, parameter_loader<Mesh<CoordT,3,Storage>, double> &pl)
{
    typedef Mesh<CoordT,3,Storage>              mesh_type;
    typedef std::complex<double>                scalar_type;
    //typedef double                              scalar_type;
    typedef typename mesh_type::point_type      point_type;
    typedef typename mesh_type::cell_type       cell_type;
    typedef typename mesh_type::face_type       face_type;

    size_t order = pl.order();
    double omega = 2*M_PI*pl.frequency();
    double alpha = 1;
    auto eps0 = parameter_loader<mesh_type, double>::eps0;
    auto mu0 = parameter_loader<mesh_type, double>::mu0;
    auto jw = scalar_type(0,omega);
    auto jwmu0 = scalar_type(0,omega*mu0);
    auto ksq = (eps0*omega)*(mu0*omega);

    auto& lua = pl.state();

    double wgz = lua["wgz"];

    disk::hho_degree_info chdi( disk::priv::hdi_named_args{
                                  .rd = (size_t) order,
                                  .cd = (size_t) order,
                                  .fd = (size_t) order
                                } );

    const auto fbs = disk::vector_basis_size(chdi.face_degree(),2, 2);
    const auto cbs = disk::vector_basis_size(chdi.cell_degree(), 3, 3);

    auto src_y = [&](const point_type& pt) -> Matrix<std::complex<double>,3,1> {
        Matrix<std::complex<double>,3,1> ret;
        ret(0) = 0;
        ret(1) = 1;
        ret(2) = 0;
        return ret;
    };

    auto tfsf_E = [&](const point_type& pt) -> Matrix<std::complex<double>,3,1> {
        Matrix<std::complex<double>,3,1> ret;
        ret(0) = 0;
        ret(1) = 0;
        ret(2) = std::sin(M_PI*pt.y()/0.0229);
        return ret;
    };

    auto tfsf_H = [&](const point_type& pt) -> Matrix<std::complex<double>,3,1> {
        Matrix<std::complex<double>,3,1> ret;
        ret(0) = 0;
        ret(1) = std::sin(M_PI*pt.y()/0.0229);
        ret(2) = 0;
        return ret;
    };

    std::vector<bool> is_dirichlet(msh.faces_size(), false);

    for (size_t i = 0; i < msh.faces_size(); i++)
    {
        auto fc = *(msh.faces_begin() + i);

        auto bi = msh.boundary_info(fc);
        if (not bi.is_boundary())
            continue;

        if (bi.is_internal())
            continue;

        auto face_tag = bi.tag();
        if ( pl.is_impedance_like(face_tag) )
            continue;

        if ( pl.is_magnetic_like(face_tag) )
            continue;

        is_dirichlet[i] = true;
    }

    maxwell_assembler_condensed<mesh_type, scalar_type> assm(msh, chdi, is_dirichlet);

    auto compute_local_contribution = [&](const cell_type& cl) -> auto {
        auto di = msh.domain_info(cl);

        auto epsr = pl.epsilon( di.tag() );
        auto mur = pl.mu( di.tag() );
        auto sigma = pl.sigma( di.tag() );
        auto Z = std::sqrt( (jw*mur*mu0) / (sigma + jw*epsr*eps0) );
        auto CR = disk::curl_reconstruction_pk(msh, cl, chdi);
        auto ST = disk::curl_hdg_stabilization(msh, cl, chdi);
        auto MM = disk::make_vector_mass_oper(msh, cl, chdi);

        auto stabparam = omega*mu0*std::sqrt(real(epsr)/real(mur));// omega*(mu0/Z);

        Matrix<scalar_type, Dynamic, Dynamic> lhs =
            (1./mur) * CR.second - (ksq*epsr - jwmu0*sigma)*MM + stabparam*ST;

        Matrix<scalar_type, Dynamic, Dynamic> lhst =
            (1./mur) * CR.second + stabparam*ST;

        Matrix<scalar_type, Dynamic, 1> rhs =
            Matrix<scalar_type, Dynamic, 1>::Zero(lhs.rows());

        const auto cb = make_vector_monomial_basis(msh, cl, chdi.reconstruction_degree());
        auto fcs = faces(msh, cl);

        Matrix<scalar_type, Dynamic, 1> dirichlet_data =
            Matrix<scalar_type, Dynamic, 1>::Zero(fcs.size() * fbs);
        /*
        if (di.tag() == 14)
        {
            disk::dynamic_vector<scalar_type> tfsf =
                        disk::dynamic_vector<scalar_type>::Zero(lhs.rows());

            auto qps = disk::integrate(msh, cl, 2*chdi.cell_degree());
            for (const auto& qp : qps)
            {
                auto phi = cb.eval_functions(qp.point());
                tfsf.head(cbs) += qp.weight() * phi * tfsf_E(qp.point())*std::exp( std::complex<double>(0.0, -1.0)*std::sqrt(ksq)*qp.point().x() );
            }

            rhs.head(cbs) += lhst*tfsf;
        }
        */

        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto& fc = fcs[i];
            auto bi = msh.boundary_info(fc);
            if (not bi.is_boundary())
                continue;

            auto tag = bi.tag();

            if (bi.is_internal())
            {
                auto n = normal(msh, cl, fc);
                if (di.tag() == 14 and bi.tag() == 1188)
                //if (di.tag() == 5)
                {
                    const auto fb = disk::make_vector_monomial_tangential_basis<mesh_type, face_type, scalar_type>(msh, fc, chdi.face_degree());

                    disk::dynamic_matrix<scalar_type> tfsf_mass =
                        disk::dynamic_matrix<scalar_type>::Zero(fbs, fbs);

                    disk::dynamic_vector<scalar_type> tfsf_rhs =
                        disk::dynamic_vector<scalar_type>::Zero(fbs);

                    disk::dynamic_vector<scalar_type> tfsf_rhsH =
                        disk::dynamic_vector<scalar_type>::Zero(fbs);

                    auto qps = disk::integrate(msh, fc, 2*chdi.face_degree());
                    for (const auto& qp : qps)
                    {
                        auto f_phi = fb.eval_functions(qp.point());
                        tfsf_mass += qp.weight() * f_phi * f_phi.transpose();
                        tfsf_rhs += qp.weight() * f_phi * tfsf_E(qp.point());
                        tfsf_rhsH += qp.weight() * f_phi * tfsf_H(qp.point());
                    }
                    //if (bi.tag() == 7)
                    {
                        disk::dynamic_vector<scalar_type> s = disk::dynamic_vector<scalar_type>::Zero(lhs.rows());
                        s.segment(cbs+i*fbs, fbs) = tfsf_mass.ldlt().solve(tfsf_rhs);
                        rhs += lhst*s;
                        //auto Y = disk::make_impedance_term(msh, fc, chdi.face_degree());
                        rhs.segment(cbs+i*fbs, fbs) -= (jwmu0/wgz)*tfsf_rhs; //Y*tfsf_mass.ldlt().solve(tfsf_rhs);
                        //rhs.segment(cbs+i*fbs, fbs) += (omega*mu0/wgz)*tfsf_rhs;
                    }

                }
            }

            if ( pl.is_impedance(tag) )
            {
                auto bnd_Z = Z;
                auto [specified, new_Z] = pl.impedance(tag);
                if (specified)
                    bnd_Z = new_Z;

                auto Y = disk::make_impedance_term(msh, fc, chdi.face_degree());
                lhs.block(cbs+i*fbs, cbs+i*fbs, fbs, fbs) += (jwmu0/bnd_Z)*Y;
            }

            if ( pl.is_plane_wave(tag) )
            {
                auto f = [&](const point_type& pt) -> Matrix<std::complex<double>,3,1> {
                    return pl.plane_wave_source(tag, pt);
                };
                auto jw = scalar_type(0,omega);
                auto [Y, y] = disk::make_plane_wave_term<std::complex<double>>(msh, fc, chdi.face_degree(), f);
                lhs.block(cbs+i*fbs, cbs+i*fbs, fbs, fbs) += (jwmu0/Z)*Y;
                rhs.segment(cbs+i*fbs, fbs) += 2.0*(jwmu0/Z)*y;
            }

            if ( pl.is_dirichlet(tag) )
            {
                auto f = [&](const point_type& pt) -> Matrix<std::complex<double>,3,1> {
                    return pl.dirichlet_data(tag, pt);
                };

                const auto fb = disk::make_vector_monomial_tangential_basis<mesh_type, face_type, scalar_type>(msh, fc, chdi.face_degree());
                const auto fbs = fb.size();
                Matrix<scalar_type, Dynamic, Dynamic> mass = Matrix<scalar_type, Dynamic, Dynamic>::Zero(fbs, fbs);
                Matrix<scalar_type, Dynamic, 1> pf = Matrix<scalar_type, Dynamic, 1>::Zero(fbs);

                const auto qps_f = integrate(msh, fc, 2*chdi.face_degree());
                for (auto& qp : qps_f)
                {
                    Matrix<scalar_type, Dynamic, 3> f_phi = fb.eval_functions(qp.point());
                    mass += qp.weight() * f_phi * f_phi.transpose();
                    pf += qp.weight() * f_phi * f(qp.point());
                }

                dirichlet_data.segment(i*fbs, fbs) += mass.ldlt().solve(pf);
            }
        }

        return std::make_tuple(lhs, rhs, dirichlet_data);
    };

    std::cout << "Assembling to triplets" << std::endl;
    for (auto& cl : msh)
    {
        auto cbs = disk::vector_basis_size(chdi.cell_degree(), 3, 3);
        auto [lhs, rhs, dirichlet_data] = compute_local_contribution(cl);
        auto [LC, bC] = disk::static_condensation(lhs, rhs, cbs);
        assm.assemble(msh, cl, LC, bC, dirichlet_data);
    }

    std::cout << "Triplets to matrix" << std::endl;
    assm.finalize();

    disk::dynamic_vector<scalar_type> sol = disk::dynamic_vector<scalar_type>::Zero(assm.syssz);

    std::cout << "Running MUMPS" << std::endl;
    sol = mumps_lu(assm.LHS, assm.RHS);

    std::vector<std::pair<scalar_type, int>> data_ex, data_ey, data_ez, data_diff, data_z;
    data_ex.resize(msh.points_size(), std::make_pair(0.0, 0));
    data_ey.resize(msh.points_size(), std::make_pair(0.0, 0));
    data_ez.resize(msh.points_size(), std::make_pair(0.0, 0));
    data_diff.resize(msh.points_size(), std::make_pair(0.0, 0));
    data_z.resize(msh.points_size(), std::make_pair(0.0, 0));

    scalar_type P_incident = 0.0;
    scalar_type P_reflected = 0.0;

    scalar_type err = 0.0;
    size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto di = msh.domain_info(cl);
        auto epsr = pl.epsilon( di.tag() );
        auto mur = pl.mu( di.tag() );
        auto sigma = pl.sigma( di.tag() );
        auto Z = std::sqrt( (jw*mur*mu0) / (sigma + jw*epsr*eps0) );

        auto [lhs, rhs, dd] = compute_local_contribution(cl);

        auto cb = disk::make_vector_monomial_basis<mesh_type, cell_type, scalar_type>(msh, cl, chdi.cell_degree());

        auto edofs = assm.get_element_dofs(msh, cl, sol);
        edofs += dd;

        Matrix<scalar_type, Dynamic, 1> esol = disk::static_decondensation(lhs, rhs, edofs);

        auto bar = barycenter(msh, cl);
        Matrix<scalar_type,Dynamic,1> esolseg = esol.segment(0, cb.size());

        auto pts = points(msh, cl);
        auto ptids = cl.point_ids();

        for (size_t i = 0; i < pts.size(); i++)
        {
            auto phi = cb.eval_functions(pts[i]);
            auto cphi = cb.eval_curls(pts[i]);
            auto ls = phi.transpose()*esolseg;
            auto ptid = ptids.at(i);
            data_ex[ptid].first += ls(0);
            data_ex[ptid].second++;
            data_ey[ptid].first += ls(1);
            data_ey[ptid].second++;
            data_ez[ptid].first += ls(2);
            data_ez[ptid].second++;

            auto cls = (1./(-jw*mur*mu0))*cphi.transpose()*esolseg;
            data_z[ptid].first += ls.norm()/cls.norm();
            data_z[ptid].second += 1;
        }





        auto fcs = faces(msh, cl);
        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto& fc = fcs[i];
            const auto fb = disk::make_vector_monomial_tangential_basis<mesh_type, face_type, scalar_type>(msh, fc, chdi.face_degree());
            const auto n  = normal(msh, cl, fc);

            auto fbs = fb.size();
            auto cbs = cb.size();
            const auto num_faces_dofs = fcs.size() * fbs;
            auto offset = cbs + i*fbs;

            Matrix<scalar_type, Dynamic, Dynamic> mass = Matrix<scalar_type, Dynamic, Dynamic>::Zero(fbs, fbs);
            Matrix<scalar_type, Dynamic, Dynamic> trace = Matrix<scalar_type, Dynamic, Dynamic>::Zero(fbs, cbs);
            Matrix<scalar_type, Dynamic, Dynamic> rhs1 = Matrix<scalar_type, Dynamic, Dynamic>::Zero(fbs, cbs + num_faces_dofs);
            Matrix<scalar_type, Dynamic, Dynamic> rhs2 = Matrix<scalar_type, Dynamic, Dynamic>::Zero(fbs, cbs + num_faces_dofs);
            rhs1.block(0, offset, fbs, fbs) = (jwmu0/Z)*Matrix<scalar_type, Dynamic, Dynamic>::Identity(fbs, fbs);

            const auto qps_f = integrate(msh, fc, 2*std::max(chdi.cell_degree(),chdi.face_degree()));
            for (auto& qp : qps_f)
            {
                Matrix<scalar_type, Dynamic, 3> f_phi = fb.eval_functions(qp.point());
                Matrix<scalar_type, Dynamic, 3> c_cphi_tmp = cb.eval_curls(qp.point());
                Matrix<scalar_type, Dynamic, 3> c_cphi = disk::vcross(c_cphi_tmp, n);

                mass += qp.weight() * f_phi * f_phi.transpose();
                trace += qp.weight() * f_phi * c_cphi.transpose();
            }
            rhs2.block(0,0,fbs,cbs) = (1./mur)*mass.ldlt().solve(trace);

            auto fc_pts = points(msh, fc);
            auto fc_ptids = fc.point_ids();

            for (size_t j = 0; j < fc_pts.size(); j++)
            {
                auto fphi = fb.eval_functions(fc_pts[j]);
                auto ls1 = fphi.transpose()*(rhs1*esol);
                auto ls2 = fphi.transpose()*(rhs2*esol);
                auto ls = fphi.transpose()*( (rhs1+rhs2)*esol );

                Matrix<scalar_type, Dynamic, 3> cphi_tmp = cb.eval_functions(fc_pts[j]);
                Matrix<scalar_type, Dynamic, 3> n_x_cphi_x_n = disk::vcross(n, disk::vcross(cphi_tmp, n));
                Matrix<scalar_type, Dynamic, 3> ccphi_tmp = cb.eval_curls(fc_pts[j]);
                Matrix<scalar_type, Dynamic, 3> ccphi_x_n = disk::vcross(ccphi_tmp, n);

                Matrix<scalar_type, Dynamic, 3> imp = (1./mur)*ccphi_x_n + (jwmu0/Z)*n_x_cphi_x_n;

                auto ptid = fc_ptids.at(j);
                data_diff.at(ptid).first += ls.norm();//std::complex<double>(ls1.norm(), ls2.norm());
                data_diff.at(ptid).second++;
            }

            auto bi = msh.boundary_info(fc);
            auto tag = bi.tag();

            //if (bi.is_internal())
            {
                if (bi.tag() == 1188)
                {
                    Matrix<scalar_type, Dynamic, 1> fdofs = esol.segment(offset, fbs);
                    for (auto& qp : qps_f)
                    {
                        Matrix<scalar_type, Dynamic, 3> fphi = fb.eval_functions(qp.point());
                        Matrix<scalar_type, 3, 1> fval = Matrix<scalar_type, 3, 1>::Zero();

                        for (size_t kk = 0; kk < fdofs.size(); kk++)
                        {
                            fval += fdofs(kk)*fphi.row(kk);
                        }

                        auto Finc = tfsf_E(qp.point());
                        P_reflected += qp.weight() * (fval).dot(Finc);
                        P_incident += qp.weight() * (Finc).dot(Finc);
                    }
                }
            }
        }

        cell_i++;
    }

    std::cout << P_reflected << std::endl;
    std::cout << P_incident << std::endl;

    std::cout << "S11 = " << 10.0*log10(  abs(P_reflected/P_incident) ) << std::endl;

    lua["s11"] = 10.0*log10( abs(P_reflected/P_incident) );

    auto tran = [](const std::pair<scalar_type, int>& nd) -> auto {
        if (nd.second == 0)
            return scalar_type(0.0);

        return nd.first/scalar_type(nd.second);
    };

    std::vector<scalar_type> nodal_ex( data_ex.size() );
    std::transform(data_ex.begin(), data_ex.end(), nodal_ex.begin(), tran);

    std::vector<scalar_type> nodal_ey( data_ey.size() );
    std::transform(data_ey.begin(), data_ey.end(), nodal_ey.begin(), tran);

    std::vector<scalar_type> nodal_ez( data_ez.size() );
    std::transform(data_ez.begin(), data_ez.end(), nodal_ez.begin(), tran);

    std::vector<double> nodal_mag( data_ex.size() );
    for (size_t i = 0; i < nodal_mag.size(); i++)
    {
        auto ex = nodal_ex[i];
        auto ey = nodal_ey[i];
        auto ez = nodal_ez[i];
        nodal_mag[i] = real(std::sqrt( ex*conj(ex) + ey*conj(ey) + ez*conj(ez) ));
    }

    std::vector<scalar_type> nodal_diff( data_diff.size() );
    std::transform(data_diff.begin(), data_diff.end(), nodal_diff.begin(), tran);

    std::vector<scalar_type> nodal_z( data_z.size() );
    std::transform(data_z.begin(), data_z.end(), nodal_z.begin(), tran);

    std::string silo_fn = lua["silo_fn"];

    disk::silo_database silo_db;
    silo_db.create(silo_fn);
    silo_db.add_mesh(msh, "mesh");

    disk::silo_nodal_variable<scalar_type> ex("ex", nodal_ex);
    silo_db.add_variable("mesh", ex);

    disk::silo_nodal_variable<scalar_type> ey("ey", nodal_ey);
    silo_db.add_variable("mesh", ey);

    disk::silo_nodal_variable<scalar_type> ez("ez", nodal_ez);
    silo_db.add_variable("mesh", ez);

    disk::silo_nodal_variable<double> emag("e_mag", nodal_mag);
    silo_db.add_variable("mesh", emag);

    disk::silo_nodal_variable<scalar_type> diff("diff", nodal_diff);
    silo_db.add_variable("mesh", diff);

    disk::silo_nodal_variable<scalar_type> z("Z", nodal_z);
    silo_db.add_variable("mesh", z);

    silo_db.close();
}


template<template<typename, size_t, typename> class Mesh, typename CoordT, typename Storage>
void
vector_wave_solver_complex_driver(Mesh<CoordT,3,Storage>& msh, const std::string& cfg_fn)
{
    using mesh_type = Mesh<CoordT,3,Storage>;
    parameter_loader<mesh_type, double> pl;

    auto run_solver = [&](){ vector_wave_solver_complex(msh, pl); };

    sol::state& lua = pl.state();
    lua["run_solver"] = run_solver;

    bool ok = pl.load(msh, cfg_fn);
    if (!ok)
        return;
}

#if 0
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
#endif

#if 0
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
#endif

#if 0
void autotest_convergence(size_t order_min, size_t order_max)
{
    using T = double;

    using Mesh = disk::simplicial_mesh<T,3>;

    std::ofstream ofs( "hho_convergence_convt.txt" );

    for (size_t i = 1; i < 4; i++)
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

            auto ci = vector_wave_solver(msh, order);
            ofs << ci.l2_error_e << " " << ci.l2_error_h << " ";// << ci.nrg_error << " " << ci.mflops << " ";
            //ofs << ci.dofs << " " << ci.nonzeros << " ";
        }

        ofs << std::endl;
    }
}

void autotest(size_t order)
{
    //for (size_t i = 0; i <= order; i++)
    //    autotest_alpha(i);

    autotest_convergence(1, order);
}

int main2(void)
{
    //_MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    autotest(3);
}
#endif

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
    while ( (ch = getopt(argc, argv, "Aa:ek:m:w:HP:")) != -1 )
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

            //case 'A':
            //    autotest(degree);
            //    return 0;

            case 'w':
                omega = M_PI*std::stod(optarg);
                break;

            case 'H':
                use_ho_stab = true;
                break;

            case 'P':
                param_filename = optarg;
                break;

            case '?':
            default:
                std::cout << "Invalid option" << std::endl;
                return 1;
        }
    }

    if (mesh_filename == nullptr or param_filename == nullptr)
    {
        std::cout << "forgot -m or -P?" << std::endl;
        return 1;
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
        disk::simplicial_mesh<coord_T, 3> msh;
        disk::load_mesh_netgen<coord_T>(mesh_filename, msh);

        //if (solve_eigvals)
        //    maxwell_eigenvalue_solver(msh, degree, stab_param);
        //else
            //vector_wave_solver(msh, degree, stab_param, omega, false);
            vector_wave_solver_complex_driver(msh, param_filename);

        return 0;
    }

    /* GMSH */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo$") ))
    {
        std::cout << "Guessed mesh format: GMSH" << std::endl;
        //disk::simplicial_mesh<T,3> msh;
        //disk::gmsh_geometry_loader< disk::simplicial_mesh<T,3> > loader;
        disk::generic_mesh<T,3> msh;
        disk::gmsh_geometry_loader< disk::generic_mesh<T,3> > loader;

        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        vector_wave_solver_complex_driver(msh, param_filename);

        return 0;
    }

    /* GMSH */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3s$") ))
    {
        std::cout << "Guessed mesh format: GMSH" << std::endl;
        disk::simplicial_mesh<T,3> msh;
        disk::gmsh_geometry_loader< disk::simplicial_mesh<T,3> > loader;
        //disk::generic_mesh<T,3> msh;
        //disk::gmsh_geometry_loader< disk::generic_mesh<T,3> > loader;

        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        vector_wave_solver_complex_driver(msh, param_filename);

        return 0;
    }

#if 0
    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        disk::cartesian_mesh<coord_T, 3> msh;
        disk::load_mesh_diskpp_cartesian<coord_T>(mesh_filename, msh);

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

    return 0;
}
