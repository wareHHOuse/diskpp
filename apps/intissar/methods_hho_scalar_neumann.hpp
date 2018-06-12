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
 * Karol Cascavita (C) 2018                     klcascavitam@unal.edu.co
 * Nicolas Pignet  (C) 2018                     nicolas.pignet@enpc.fr
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

#pragma once

#include <iterator>

#include "core/common/eigen.hpp"
#include "core/revolution/bases"
#include "core/revolution/quadratures"
#include "core/revolution/implementation_bits/methods_hho.hpp"
#include "BoundaryConditionsRobin.hpp"

using namespace Eigen;

namespace revolution
{


/********************************************************************/
template<typename Mesh, typename T>
class diffusion_condensed_assembler2
{
    typedef disk::mechanics::BoundaryConditionsScalar2<Mesh> boundary_type;
    boundary_type                       m_bnd;
    std::vector<size_t>                 compress_table;
    std::vector<size_t>                 expand_table;

    hho_degree_info                     di;

    std::vector< Triplet<T> >           triplets;

    size_t      num_all_faces, num_dirichlet_faces, num_neumann_faces, num_robin_faces, num_other_faces;

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

    diffusion_condensed_assembler2(const Mesh& msh, hho_degree_info hdi, const boundary_type& bnd)
        : di(hdi), m_bnd(bnd)
    {
        auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
            const auto face_id=eid.second;
            return m_bnd.is_dirichlet_face(face_id);
        };
        auto is_neumann = [&](const typename Mesh::face_type& fc) -> bool {
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
            const auto face_id=eid.second;
            return m_bnd.is_neumann_face(face_id);
        };

         auto is_robin = [&](const typename Mesh::face_type& fc) -> bool {
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
            const auto face_id=eid.second;
            return m_bnd.is_robin_face(face_id);
        };


        num_all_faces = msh.faces_size();
        num_dirichlet_faces = std::count_if(msh.faces_begin(), msh.faces_end(), is_dirichlet);
        num_neumann_faces = std::count_if(msh.faces_begin(), msh.faces_end(), is_neumann);
        num_robin_faces = std::count_if(msh.faces_begin(), msh.faces_end(), is_robin);
        num_other_faces = num_all_faces - num_dirichlet_faces ;

        compress_table.resize( num_all_faces );
        expand_table.resize( num_other_faces );

        size_t compressed_offset = 0;
        for (size_t i = 0; i < num_all_faces; i++)
        {
            auto fc = *std::next(msh.faces_begin(), i);
            if ((!is_dirichlet(fc)) )
            {
                compress_table.at(i) = compressed_offset;
                expand_table.at(compressed_offset) = i;
                compressed_offset++;
            }
        }

        auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);

        auto system_size = fbs * num_other_faces;

        LHS = SparseMatrix<T>( system_size, system_size );
        RHS = Matrix<T, Dynamic, 1>::Zero( system_size );
    }

    #if 0
    void dump_tables() const
    {
        std::cout << "Compress table: " << std::endl;
        for (size_t i = 0; i < compress_table.size(); i++)
            std::cout << i << " -> " << compress_table.at(i) << std::endl;
    }
    #endif

    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs,
             const Matrix<T, Dynamic, 1>& rhs)
    {
        auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
            const auto face_id=eid.second;
            return m_bnd.is_dirichlet_face(face_id);
        };

        auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);
        auto fcs = faces(msh, cl);

        std::vector<assembly_index> asm_map;
        asm_map.reserve(fcs.size()*fbs);

        Matrix<T, Dynamic, 1> data = Matrix<T, Dynamic, 1>::Zero(fcs.size()*fbs);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto face_offset = priv::offset(msh, fc);
            auto face_LHS_offset = compress_table.at(face_offset)*fbs;

            bool dirichlet = is_dirichlet(fc);


            for (size_t i = 0; i < fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !(dirichlet)) );


            if (dirichlet)
            {
                auto fb = make_scalar_monomial_basis(msh, fc, di.face_degree());
                auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
                if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
                const auto face_id=eid.second;
                auto dirichlet_bf = m_bnd.dirichlet_boundary_func(face_id);
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_bf, di.face_degree());

                data.block(face_i*fbs, 0, fbs, 1) = mass.llt().solve(rhs);
            }

        }

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if ( asm_map[j].assemble() )
                    triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
                else
                    RHS(asm_map[i]) -= lhs(i,j) * data(j);
            }

            RHS(asm_map[i]) += rhs(i);
        }
    } // assemble()

    void neumann_boundary_function (const Mesh& msh, const boundary_type& bnd)
    {

        if (bnd.nb_faces_neumann() > 0)
        {
            for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
            {
                const auto bfc = *itor;
                const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
                if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
                const auto face_id = eid.second;
                ///on a la face_id

                if (bnd.is_neumann_face(face_id))
                {
                    if (bnd.is_dirichlet_face(face_id))
                    {
                        switch (bnd.dirichlet_boundary_type(face_id))
                        {
                            case disk::mechanics::DIRICHLET:
                                throw std::invalid_argument("You tried to impose both Dirichlet and Neumann conditions on the same face");
                                break;
                            default:
                                throw std::logic_error("Unknown Dirichlet Conditions");
                            break;
                        }
                    }
                    else
                    {
                        // Add Warning for Robin condition

                        const size_t face_degree   = di.face_degree();
                        const size_t num_face_dofs = scalar_basis_size(face_degree, Mesh::dimension - 1);
                        std::vector<assembly_index> asm_map;
                        asm_map.reserve(num_face_dofs);

                        auto face_offset = face_id;
                        auto face_LHS_offset = face_offset*num_face_dofs;

                        bool isneumann = bnd.is_neumann_face(face_id);

                        for (size_t i = 0; i < num_face_dofs; i++)
                            asm_map.push_back( assembly_index(face_LHS_offset+i, !(isneumann)) );

                        auto fb = make_scalar_monomial_basis(msh, bfc, face_degree);
                        Matrix<T, Dynamic, 1> neumann = make_rhs(msh, bfc, fb, bnd.neumann_boundary_func(face_id));

                        assert(neumann.size() == num_face_dofs);

                        for (size_t i = 0; i < num_face_dofs; i++)
                        {
                            if (!asm_map[i].assemble())
                                continue;

                            RHS(asm_map[i]) += neumann[i];
                        }
                    }
                }
            }
        }
    }

    void robin_boundary_function (const Mesh& msh, const boundary_type& bnd)
    {

        if (bnd.nb_faces_robin() > 0)
        {
            for (auto itor = msh.boundary_faces_begin(); itor!= msh.boundary_faces_end(); itor++)
            {
                const auto bfc =*itor ;
                const auto eid = find_element_id (msh.faces_begin(), msh.faces_end(), bfc);
                if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
                const auto face_id = eid.second;

                if (bnd.is_robin_face(face_id))
                {
                    if (bnd.is_neumann_face(face_id))
                    {
                        switch (bnd.neumann_boundary_type(face_id))
                        {
                            case disk::mechanics::NEUMANN:
                                throw std::invalid_argument("You tried to impose both Neumann and Robin conditions on the same face");
                                break;
                            default:
                                throw std::logic_error("Unknown Neumann Conditions");
                            break;

                        }
                    }
                    else if (bnd.is_dirichlet_face(face_id))
                    {
                        switch (bnd.dirichlet_boundary_type(face_id))
                        {
                            case disk::mechanics::DIRICHLET:
                                throw std::invalid_argument("You tried to impose both Dirichlet and Robin conditions on the same face");
                                break;
                            default:
                                throw std::logic_error("Unknown Dirichlet Conditions");
                            break;
                        }
                    }
                    else
                    {
                        const size_t face_degree = di.face_degree();
                        const size_t num_face_dofs = scalar_basis_size(face_degree, Mesh::dimension -1);
                        std::vector<assembly_index> asm_map;
                        asm_map.reserve(num_face_dofs);

                        auto face_offset = face_id;
                        auto face_LHS_offset = face_offset*num_face_dofs;

                        for (size_t i = 0; i < num_face_dofs; i++)
                            asm_map.push_back( assembly_index(face_LHS_offset+i, true) );

                        auto fb = make_scalar_monomial_basis(msh, bfc, face_degree);
                        Matrix<T, Dynamic, 1> robin = make_rhs(msh,bfc,fb,bnd.robin_boundary_func(face_id));
                        assert (robin.size() == num_face_dofs);

                        Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, bfc, fb, face_degree);

                        for (size_t i = 0; i < num_face_dofs; i++)
                        {
                            if (asm_map[i].assemble())
                                continue;

                            RHS(asm_map[i]) += robin[i];

                            for (size_t j = 0; j < num_face_dofs; j++)
                            {
                                if ( asm_map[j].assemble() )
                                    triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], mass(i,j)) );
                                //else
                                //    RHS(asm_map[i]) -= lhs(i,j) * data(j); //dirichlet
                            }
                        }
                    }
                }
            }
        }
    }



    Matrix<T, Dynamic,1>
    take_local_data(const Mesh& msh, const typename Mesh::cell_type& cl,
    const Matrix<T, Dynamic, 1>& solution)
    {
        auto facdeg = di.face_degree();
        auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);
        auto fcs = faces(msh, cl);

        auto num_faces = fcs.size();

        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(num_faces*fbs);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];

            auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
                auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
                if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
                const auto face_id=eid.second;
                return m_bnd.is_dirichlet_face(face_id);
            };

            bool dirichlet = is_dirichlet(fc);


             if (dirichlet)
            {
                auto fb = make_scalar_monomial_basis(msh, fc, di.face_degree());
                auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
                if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
                const auto face_id=eid.second;
                auto dirichlet_bf=m_bnd.dirichlet_boundary_func(face_id);
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_bf, di.face_degree());
                ret.block(face_i*fbs, 0, fbs, 1) = mass.llt().solve(rhs);
            }
            else
            {
                auto face_offset = priv::offset(msh, fc);
                auto face_SOL_offset = compress_table.at(face_offset)*fbs;
                ret.block(face_i*fbs, 0, fbs, 1) = solution.block(face_SOL_offset, 0, fbs, 1);
            }
        }

        return ret;
    }

    void finalize(void)
    {
        LHS.setFromTriplets( triplets.begin(), triplets.end() );
        triplets.clear();

        dump_sparse_matrix(LHS, "diff.dat");
    }

    size_t num_assembled_faces() const
    {
        return num_other_faces;
    }
};


/********************************************************************/
template<typename Mesh, typename BoundaryType>
auto make_diffusion_assembler2(const Mesh& msh, hho_degree_info hdi,
         const BoundaryType& m_bnd)
{
    return diffusion_condensed_assembler2<Mesh, typename Mesh::scalar_type> (msh, hdi, m_bnd);
}

}// revolution
