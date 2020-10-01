/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
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

template<typename Mesh>
class laplacian_H0t_assembler
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

    laplacian_H0t_assembler(const Mesh& pmsh,
                            const disk::hho_degree_info& pchdi)
        : msh(pmsh), chdi(pchdi)
    {
        auto cbs   = disk::vector_basis_size(chdi.cell_degree(), Mesh::dimension, Mesh::dimension);
        auto fbs_b = disk::scalar_basis_size(chdi.face_degree(), Mesh::dimension-1);
        auto fbs_i = disk::vector_basis_size(chdi.face_degree(), Mesh::dimension-1, Mesh::dimension);

        syssz = cbs * msh.cells_size();
        for (auto& fc : faces(msh))
        {
            if (msh.is_boundary(fc))
                syssz += fbs_b;
            else
                syssz += fbs_i;
        }

        LHS = SparseMatrix<T>(syssz, syssz);
        RHS = Matrix<T, Dynamic, 1>::Zero(syssz);

        compress_table.resize( msh.faces_size() );
        expand_table.resize( msh.internal_faces_size() );

        size_t original_ofs = 0;
        size_t compressed_ofs = 0;
        size_t face_i = 0;
        for (auto& fc : faces(msh))
        {
            if ( not msh.is_boundary(fc) )
            {
                compress_table[face_i] = compressed_ofs;
                //expand_table[compressed_ofs] = original_ofs;
                compressed_ofs += fbs_i;
                original_ofs += fbs_i;
            }
            else
            {
                original_ofs += fbs_b;
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
        auto cbs = disk::vector_basis_size(chdi.cell_degree(), Mesh::dimension, Mesh::dimension);
        auto fbs_b = disk::scalar_basis_size(chdi.face_degree(), Mesh::dimension-1);
        auto fbs_i = disk::vector_basis_size(chdi.face_degree(), Mesh::dimension-1, Mesh::dimension);

        auto cell_ofs = offset(msh, cl) * cbs;

        auto fcs = faces(msh, cl);
        std::vector<size_t> fbss;
        std::vector<size_t> fbos;
        size_t fbo = cbs;
        for (auto& fc : fcs)
        {
            if (msh.is_boundary(fc))
            {
                fbss.push_back(fbs_b);
                fbos.push_back(fbo);
                fbo += fbs_b;
            }
            else
            {
                fbss.push_back(fbs_i);
                fbos.push_back(fbo);
                fbo += fbs_i;
            }
        }


        /* Curl-curl term, cell part */
        for (size_t i = 0; i < cbs; i++)
        {
            auto gi = cell_ofs + i;
            for (size_t j = 0; j < cbs; j++)
            {
                auto gj = cell_ofs + j;
                triplets.push_back( Triplet<T>(gi, gj, lhs(i,j)) );
            }

            RHS(gi) = rhs(i);
        }

        /* Offdiag parts */
        for (size_t i = 0; i < cbs; i++)
        {
            auto gi = cell_ofs + i;
            for (size_t fj = 0; fj < fcs.size(); fj++)
            {
                if ( false and msh.is_boundary(fcs[fj]) )
                    continue;

                auto fbs_j = fbss[fj];
                for (size_t j = 0; j < fbs_j; j++)
                {
                    auto gj = cbs * msh.cells_size() + compress_table.at( offset(msh, fcs[fj]) ) + j;

                    auto li = i;
                    auto lj = fbos[fj] + j;
                    triplets.push_back( Triplet<T>(gi, gj, lhs(li,lj)) );
                }
            }
        }

        for (size_t fi = 0; fi < fcs.size(); fi++)
        {
            if ( msh.is_boundary(fcs[fi]) )
                continue;

            auto fbs_i = fbss[fi];

            for (size_t i = 0; i < fbs_i; i++)
            {
                auto gi = cbs * msh.cells_size() + compress_table.at( offset(msh, fcs[fi]) ) + i;
                for (size_t j = 0; j < cbs; j++)
                {
                    auto gj = cell_ofs + j;
                    auto li = fbos[fi] + i;
                    auto lj = j;
                    triplets.push_back( Triplet<T>(gi, gj, lhs(li,lj)) );
                }
            }

            for (size_t fj = 0; fj < fcs.size(); fj++)
            {
                auto fbs_j = fbss[fj];
                if ( false and msh.is_boundary(fcs[fj]) )
                {
                    auto bid = msh.boundary_id(fcs[fj]);
                    if ( is_homo_dirichlet(bid) )
                        continue;

                    using namespace std::placeholders;
                    auto gD_fun = std::bind(gD, bid, _1);
                    auto fb = make_vector_monomial_tangential_basis(msh, fcs[fj], chdi.face_degree());
                    disk::dynamic_vector<T> dirichlet_data = project_function(msh, fcs[fj], fb, gD_fun);
                    auto lofsi = fbos[fi];
                    auto lofsj = fbos[fj];
                    auto cofsj = cbs*msh.cells_size() + compress_table.at( offset(msh, fcs[fj]) );
                    RHS.segment(cofsj, fbs_j) -= lhs.block(lofsi, lofsj, fbs_i, fbs_j)*dirichlet_data;
                }
                else
                {
                    auto cofsi = cbs*msh.cells_size() + compress_table.at( offset(msh, fcs[fi]) );
                    auto cofsj = cbs*msh.cells_size() + compress_table.at( offset(msh, fcs[fj]) );

                    for (size_t i = 0; i < fbs_i; i++)
                    {
                        for(size_t j = 0; j < fbs_j; j++)
                        {
                            auto lofsi = fbos[fi];
                            auto lofsj = fbos[fj];
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

template<typename Mesh>
class maxwell_assembler_condensed
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

    maxwell_assembler_condensed(const Mesh& pmsh,
                      const disk::hho_degree_info& pchdi)
        : msh(pmsh), chdi(pchdi)
    {
        auto fbs = disk::vector_basis_size(chdi.face_degree(), Mesh::dimension-1, Mesh::dimension-1);
        syssz = fbs*msh.internal_faces_size();

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

    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhsc,
             const Matrix<T, Dynamic, 1>& rhs)
    {
        auto num_internal_faces = msh.internal_faces_size();
        auto fbs = disk::vector_basis_size(chdi.face_degree(), Mesh::dimension-1, Mesh::dimension-1);

        auto fcs = faces(msh, cl);
        for (size_t fi = 0; fi < fcs.size(); fi++)
        {
            if ( msh.is_boundary(fcs[fi]) )
                continue;

            for (size_t fj = 0; fj < fcs.size(); fj++)
            {
                if ( msh.is_boundary(fcs[fj]) )
                    continue;

                auto cofsi = fbs * compress_table.at( offset(msh, fcs[fi]) );
                auto cofsj = fbs * compress_table.at( offset(msh, fcs[fj]) );

                for (size_t i = 0; i < fbs; i++)
                {
                    for(size_t j = 0; j < fbs; j++)
                    {
                        auto lofsi = fi*fbs;
                        auto lofsj = fj*fbs;
                        triplets.push_back( Triplet<T>(cofsi+i, cofsj+j, lhsc(lofsi+i, lofsj+j)) );
                    }
                }
            }
        }
    }

    void finalize(void)
    {
        LHS.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();
        std::cout << "Solving for " << syssz << " DOFs. ";
        std::cout << "Matrix has " << LHS.nonZeros() << " nonzeros." << std::endl; 
    }

    disk::dynamic_vector<T>
    get_expanded_solution(const Mesh& msh, disk::dynamic_vector<T>& sol)
    {
        auto num_faces = msh.internal_faces_size();
        auto fbs = disk::vector_basis_size(chdi.face_degree(), Mesh::dimension-1, Mesh::dimension-1);

        disk::dynamic_vector<T> ret = disk::dynamic_vector<T>::Zero( fbs*msh.faces_size() );

        for (size_t i = 0; i < num_faces; i++)
        {
            auto in_offset = i*fbs;
            auto out_offset = expand_table.at(i)*fbs;
            ret.segment(out_offset, fbs) = sol.segment(in_offset, fbs);
        }

        return ret;
    }

    disk::dynamic_vector<T>
    get_element_dofs(const Mesh& msh, const typename Mesh::cell& cl, disk::dynamic_vector<T>& sol)
    {
        auto fbs = disk::vector_basis_size(chdi.face_degree(), Mesh::dimension-1, Mesh::dimension-1);
        auto fcs = faces(msh, cl);
        disk::dynamic_vector<T> ret = disk::dynamic_vector<T>::Zero( fbs*fcs.size() );

        size_t i = 0;
        for (auto& fc : fcs)
        {
            if (not msh.is_boundary(fc))
            {
                auto ofs = fbs*compress_table.at( offset(msh, fc) );
                ret.segment(i*fbs, fbs) = sol.segment(ofs, fbs);
            }
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




//#define USE_STATIC_CONDENSATION

template<typename Mesh>
void laplacian_H0t_solver(Mesh& msh, size_t order)
{
    using T = typename Mesh::coordinate_type;

    disk::hho_degree_info chdi( { .rd = (size_t) order+1,
                                  .cd = (size_t) order,
                                  .fd = (size_t) order } );

    auto rhs_fun = [&](const typename Mesh::point_type& pt) -> Matrix<T, 3, 1> {
        Matrix<T, 3, 1> ret;
        ret(0) = 0.0;
        ret(1) = 0.0;
        ret(2) = -2*M_PI * M_PI * std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
        return ret;
    };

    auto sol_fun = [&](const typename Mesh::point_type& pt) -> Matrix<T, 3, 1> {
        Matrix<T, 3, 1> ret;
        ret(0) = 0.0;
        ret(1) = 0.0;
        ret(2) = std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
        return ret;
    };

    laplacian_H0t_assembler<Mesh> assm(msh, chdi);

    T omega = M_PI;

    std::cout << "Assembling to triplets" << std::endl;
    for (auto& cl : msh)
    {

        auto LP = make_vector_lapl_H0t_oper(msh, cl, chdi);
        auto LS = make_vector_lapl_H0t_stab(msh, cl, chdi);

        auto cb = make_vector_monomial_basis(msh, cl, chdi.cell_degree());
        Matrix<T, Dynamic, Dynamic> lhs = -LP.second + LS;
        lhs.block(0,0,cb.size(),cb.size()) += omega*omega*make_mass_matrix(msh, cl, cb);
        Matrix<T, Dynamic, 1> rhs(lhs.rows());

        rhs.segment(0, cb.size()) = make_rhs(msh, cl, cb, rhs_fun, 1);

        assm.assemble(msh, cl, lhs, rhs);
    }

    std::cout << "Triplets to matrix" << std::endl;
    assm.finalize();

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.syssz);

    std::cout << "Running pardiso" << std::endl;
    disk::solvers::pardiso_params<T> pparams;
    pparams.report_factorization_Mflops = true;
    mkl_pardiso(pparams, assm.LHS, assm.RHS, sol);

    /*
    disk::solvers::conjugated_gradient_params<T> cgp;
    cgp.max_iter = assm.LHS.rows();
    cgp.verbose = true;
    conjugated_gradient(cgp, assm.LHS, assm.RHS, sol);
    */

    std::vector<T> data_ux, data_uy, data_uz;

    T err = 0.0; size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto cb = make_vector_monomial_basis(msh, cl, chdi.cell_degree());

        Matrix<T, Dynamic, Dynamic> MMe = disk::make_mass_matrix(msh, cl, cb, 1);
        Matrix<T, Dynamic, 1> arhs = disk::make_rhs(msh, cl, cb, sol_fun, 1);
        Matrix<T, Dynamic, 1> asol = MMe.llt().solve(arhs);
        Matrix<T, Dynamic, 1> lsol = sol.segment(cell_i*cb.size(), cb.size());

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

        Matrix<T, Dynamic, 1> diff = lsol - asol;

        err += diff.dot(MMe*diff);

        cell_i++;
    }

    std::cout << "h = " << disk::average_diameter(msh) << " ";
    std::cout << "err = " << std::sqrt(err) << std::endl;

    disk::silo_database silo_db;
    silo_db.create("maxwell.silo");
    silo_db.add_mesh(msh, "mesh");

    disk::silo_zonal_variable<T> ux("ux", data_ux);
    silo_db.add_variable("mesh", ux);

    disk::silo_zonal_variable<T> uy("uy", data_uy);
    silo_db.add_variable("mesh", uy);

    disk::silo_zonal_variable<T> uz("uz", data_uz);
    silo_db.add_variable("mesh", uz);

    silo_db.add_expression("u", "{ux, uy, uz}", DB_VARTYPE_VECTOR);
    silo_db.add_expression("mag_u", "magnitude(u)", DB_VARTYPE_SCALAR);

    silo_db.close();
}

template<typename Mesh>
void vector_wave_solver(Mesh& msh, size_t order)
{
    using T = typename Mesh::coordinate_type;

    disk::hho_degree_info chdi( { .rd = (size_t) order+1,
                                  .cd = (size_t) order,
                                  .fd = (size_t) order } );

    auto rhs_fun = [&](const typename Mesh::point_type& pt) -> Matrix<T, 3, 1> {
        Matrix<T, 3, 1> ret;
        ret(0) = 0.0;
        ret(1) = 0.0;
        ret(2) = M_PI * M_PI * std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
        return ret;
    };

    auto sol_fun = [&](const typename Mesh::point_type& pt) -> Matrix<T, 3, 1> {
        Matrix<T, 3, 1> ret;
        ret(0) = 0.0;//M_PI*std::sin(M_PI*pt.x());
        ret(1) = 0.0;
        ret(2) = std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
        return ret;
    };

#ifdef USE_STATIC_CONDENSATION
    maxwell_assembler_condensed<Mesh> assm(msh, chdi);
#else
    maxwell_assembler<Mesh> assm(msh, chdi);
#endif

    T omega = M_PI;

    T norm_C  = 0.0;
    T norm_Ct = 0.0;
    T norm_S  = 0.0;
    T norm_M  = 0.0;

    std::cout << "Assembling to triplets" << std::endl;
    for (auto& cl : msh)
    {

        //auto LP = make_vector_hho_laplacian(msh, cl, chdi);
        //auto LS = make_vector_hdg_stabilization(msh, cl, chdi);

        auto CR = disk::make_vector_hho_curl_impl(msh, cl, chdi);
        auto ST = make_vector_hho_curl_stab(msh, cl, chdi);
        auto MM = make_vector_mass_oper(msh, cl, chdi);

        Matrix<T, Dynamic, Dynamic> lhs = CR.second + ST - (omega*omega)*MM;
        Matrix<T, Dynamic, 1> rhs(lhs.rows());

        auto cb = make_vector_monomial_basis(msh, cl, chdi.cell_degree());
        rhs.segment(0, cb.size()) = make_rhs(msh, cl, cb, rhs_fun, 1);

        //VectorXcd eivals = lhs.block(0,0,cb.size(), cb.size()).eigenvalues();
        //std::cout << eivals.transpose() << std::endl;

#ifdef USE_STATIC_CONDENSATION
        auto [LC, bC] = disk::static_condensation(lhs, rhs, cb.size());
        assm.assemble(msh, cl, LC, bC);
#else
        assm.assemble(msh, cl, lhs, rhs);
#endif

        Matrix<T, Dynamic, 1> prj = project_tangent(msh, cl, chdi, sol_fun);

        norm_C += prj.dot(CR.second*prj);
        norm_S += prj.dot(ST*prj);
        norm_M += prj.dot(MM*prj);


        Matrix<T, Dynamic, Dynamic> CC = make_curl_curl_matrix(msh, cl, cb);
        Matrix<T, Dynamic, 1> pc = project_function(msh, cl, chdi.cell_degree(), sol_fun);
        norm_Ct += pc.dot(CC*pc);
    }

    std::cout << std::sqrt(norm_C) << " " << std::sqrt(norm_S) << " " << std::sqrt(norm_M) << std::endl;
    std::cout << std::sqrt(norm_Ct) << std::endl;

    std::cout << "Triplets to matrix" << std::endl;
    assm.finalize();

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.syssz);

    
    std::cout << "Running pardiso" << std::endl;
    disk::solvers::pardiso_params<T> pparams;
    pparams.report_factorization_Mflops = true;
    mkl_pardiso(pparams, assm.LHS, assm.RHS, sol);
    
    /*
    disk::solvers::conjugated_gradient_params<T> cgp;
    cgp.max_iter = assm.LHS.rows();
    cgp.verbose = true;
    conjugated_gradient(cgp, assm.LHS, assm.RHS, sol);
    */

    std::vector<T> data_ux, data_uy, data_uz;

    T err = 0.0; size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto cb = make_vector_monomial_basis(msh, cl, chdi.cell_degree());

#ifdef USE_STATIC_CONDENSATION
        auto CR = disk::make_vector_hho_curl_impl(msh, cl, chdi);
        auto ST = make_vector_hho_curl_stab(msh, cl, chdi);
        auto MM = make_vector_mass_oper(msh, cl, chdi);
        Matrix<T, Dynamic, Dynamic> lhs = CR.second + ST - (omega*omega)*MM;
        Matrix<T, Dynamic, 1> rhs(lhs.rows());
        rhs.segment(0, cb.size()) = make_rhs(msh, cl, cb, rhs_fun, 1);
        auto edofs = assm.get_element_dofs(msh, cl, sol);
        auto esol = disk::static_decondensation(lhs, rhs, edofs);
        Matrix<T, Dynamic, 1> lsol = esol.segment(0, cb.size());
#else
        Matrix<T, Dynamic, Dynamic> MMe = disk::make_mass_matrix(msh, cl, cb, 1);
        Matrix<T, Dynamic, 1> arhs = disk::make_rhs(msh, cl, cb, sol_fun, 1);
        Matrix<T, Dynamic, 1> asol = MMe.llt().solve(arhs);
        Matrix<T, Dynamic, 1> lsol = sol.segment(cell_i*cb.size(), cb.size());
#endif

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

        Matrix<T, Dynamic, 1> diff = lsol - asol;

        err += diff.dot(MMe*diff);

        cell_i++;
    }

    std::cout << "h = " << disk::average_diameter(msh) << " ";
    std::cout << "err = " << std::sqrt(err) << std::endl;

    disk::silo_database silo_db;
    silo_db.create("maxwell.silo");
    silo_db.add_mesh(msh, "mesh");

    disk::silo_zonal_variable<T> ux("ux", data_ux);
    silo_db.add_variable("mesh", ux);

    disk::silo_zonal_variable<T> uy("uy", data_uy);
    silo_db.add_variable("mesh", uy);

    disk::silo_zonal_variable<T> uz("uz", data_uz);
    silo_db.add_variable("mesh", uz);

    silo_db.add_expression("u", "{ux, uy, uz}", DB_VARTYPE_VECTOR);
    silo_db.add_expression("mag_u", "magnitude(u)", DB_VARTYPE_SCALAR);

    silo_db.close();
}

template<typename Mesh>
void maxwell_eigenvalue_solver(Mesh& msh, size_t order)
{
    using T = typename Mesh::coordinate_type;

    disk::hho_degree_info chdi( { .rd = (size_t) order+1,
                                  .cd = (size_t) order,
                                  .fd = (size_t) order } );


    maxwell_eigenvalue_assembler<Mesh> assm(msh, chdi);

    std::cout << "Assembling to triplets" << std::endl;
    for (auto& cl : msh)
    {
        auto CR = disk::make_vector_hho_curl_impl(msh, cl, chdi);
        auto ST = make_vector_hho_curl_stab(msh, cl, chdi);
        auto MM = make_vector_mass_oper(msh, cl, chdi);
        assm.assemble(msh, cl, CR.second+ST, MM);
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
        std::cout << hho_eigvals(i) << std::endl;
    
}


int main(int argc, char **argv)
{
    using T = double;

    if (argc != 3)
    {
        std::cout << "params!" << std::endl;
        return 1;
    }

    int order = std::stoi(argv[1]);
    if (order < 0)
    {
        std::cout << "Order must be >= 0" << std::endl;
        return 1;
    }

    int mesh = std::stoi(argv[2]);

    using Mesh = disk::simplicial_mesh<T,3>;
    //using Mesh = disk::cartesian_mesh<T, 3>;
    //using Mesh = disk::generic_mesh<T, 3>;

    Mesh msh;

    

    std::stringstream ss;
    ss << "../../../meshes/3D_tetras/netgen/cube" << mesh << ".mesh";
    disk::load_mesh_netgen(ss.str().c_str(), msh);

    //disk::load_mesh_diskpp_cartesian("../../../meshes/3D_hexa/diskpp/testmesh-16-16-16.hex", msh);


    //disk::load_mesh_fvca6_3d("../../../meshes/3D_general/fvca6/dbls_10.msh", msh);
    //disk::load_mesh_fvca6_3d("../../../meshes/3D_hexa/fvca6/hexa_2x2x2.msh", msh);
    //disk::load_mesh_fvca6_3d("../../../meshes/3D_tetras/fvca6/tet.0.msh", msh);



    std::cout << "Cells: " << msh.cells_size() << std::endl;

    

    

    //laplacian_H0t_solver(msh, order);
    vector_wave_solver(msh, order);
    //maxwell_eigenvalue_solver(msh, order);


}

