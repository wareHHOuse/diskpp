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
 * Karol Cascavita (C) 2018                     karol.cascavita@enpc.fr
 * Nicolas Pignet  (C) 2018, 2019               nicolas.pignet@enpc.fr
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

#include "bases/bases.hpp"
#include "boundary_conditions/boundary_conditions.hpp"
#include "common/eigen.hpp"
#include "quadratures/quadratures.hpp"

using namespace Eigen;

namespace disk
{

namespace priv
{

template<typename Mesh>
size_t
offset(const Mesh& msh, const typename Mesh::cell_type& cl)
{
    auto itor = std::lower_bound(msh.cells_begin(), msh.cells_end(), cl);
    if (itor == msh.cells_end())
        throw std::logic_error("Cell not found: this is likely a bug.");

    return std::distance(msh.cells_begin(), itor);
}

template<typename Mesh>
size_t
offset(const Mesh& msh, const typename Mesh::face_type& fc)
{
    auto itor = std::lower_bound(msh.faces_begin(), msh.faces_end(), fc);
    if (itor == msh.faces_end())
        throw std::logic_error("Face not found: this is likely a bug.");

    return std::distance(msh.faces_begin(), itor);
}

} // priv

// assembler for scalar primal problem with HHO like diffusion problem
template<typename Mesh>
class diffusion_condensed_assembler
{
    using T = typename Mesh::coordinate_type;
    typedef disk::scalar_boundary_conditions<Mesh> boundary_type;

    std::vector<size_t>     compress_table;
    std::vector<size_t>     expand_table;
    hho_degree_info         di;
    std::vector<Triplet<T>> triplets;
    bool                    use_bnd;

    size_t num_all_faces, num_dirichlet_faces, num_other_faces, system_size;

    class assembly_index
    {
        size_t idx;
        bool   assem;

      public:
        assembly_index(size_t i, bool as) : idx(i), assem(as) {}

        operator size_t() const
        {
            if (!assem)
                throw std::logic_error("Invalid assembly_index");

            return idx;
        }

        bool
        assemble() const
        {
            return assem;
        }

        friend std::ostream&
        operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };

  public:
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    SparseMatrix<T> LHS;
    vector_type     RHS;

    diffusion_condensed_assembler(const Mesh& msh, hho_degree_info hdi) : di(hdi), use_bnd(false)
    {
        auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool { return msh.is_boundary(fc); };

        num_all_faces       = msh.faces_size();
        num_dirichlet_faces = std::count_if(msh.faces_begin(), msh.faces_end(), is_dirichlet);
        num_other_faces     = num_all_faces - num_dirichlet_faces;

        compress_table.resize(num_all_faces);
        expand_table.resize(num_other_faces);

        size_t compressed_offset = 0;
        for (size_t i = 0; i < num_all_faces; i++)
        {
            const auto fc = *std::next(msh.faces_begin(), i);
            if (!is_dirichlet(fc))
            {
                compress_table.at(i)               = compressed_offset;
                expand_table.at(compressed_offset) = i;
                compressed_offset++;
            }
        }

        const auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension - 1);
        system_size    = fbs * num_other_faces;

        LHS = SparseMatrix<T>(system_size, system_size);
        RHS = vector_type::Zero(system_size);
    }

    diffusion_condensed_assembler(const Mesh& msh, hho_degree_info hdi, const boundary_type& bnd) :
      di(hdi), use_bnd(true)
    {
        auto is_dirichlet = [&](const typename Mesh::face& fc) -> bool {
            const auto fc_id = msh.lookup(fc);
            return bnd.is_dirichlet_face(fc_id);
        };

        num_all_faces       = msh.faces_size();
        num_dirichlet_faces = std::count_if(msh.faces_begin(), msh.faces_end(), is_dirichlet);
        num_other_faces     = num_all_faces - num_dirichlet_faces;

        compress_table.resize(num_all_faces);
        expand_table.resize(num_other_faces);

        size_t compressed_offset = 0;
        for (size_t i = 0; i < num_all_faces; i++)
        {
            const auto fc = *std::next(msh.faces_begin(), i);
            if (!is_dirichlet(fc))
            {
                compress_table.at(i)               = compressed_offset;
                expand_table.at(compressed_offset) = i;
                compressed_offset++;
            }
        }

        const auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension - 1);
        system_size    = fbs * num_other_faces;

        LHS = SparseMatrix<T>(system_size, system_size);
        RHS = vector_type::Zero(system_size);
    }

    template<typename Function>
    void
    assemble(const Mesh&                     msh,
             const typename Mesh::cell_type& cl,
             const matrix_type&              lhs,
             const vector_type&              rhs,
             const Function&                 dirichlet_bf)
    {
        if (use_bnd)
            throw std::invalid_argument("diffusion_assembler: you have to use boundary type");

        auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool { return msh.is_boundary(fc); };

        const auto fbs    = scalar_basis_size(di.face_degree(), Mesh::dimension - 1);
        const auto fcs    = faces(msh, cl);
        const auto fcs_id = faces_id(msh, cl);

        std::vector<assembly_index> asm_map;
        asm_map.reserve(fcs.size() * fbs);

        vector_type dirichlet_data = vector_type::Zero(fcs.size() * fbs);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            const auto fc              = fcs[face_i];
            const auto face_offset     = fcs_id[face_i]; // priv::offset(msh, fc);
            const auto face_LHS_offset = compress_table.at(face_offset) * fbs;

            const bool dirichlet = is_dirichlet(fc);

            for (size_t i = 0; i < fbs; i++)
                asm_map.push_back(assembly_index(face_LHS_offset + i, !dirichlet));

            if (dirichlet)
            {
                dirichlet_data.block(face_i * fbs, 0, fbs, 1) =
                  project_function(msh, fc, di.face_degree(), dirichlet_bf, di.face_degree());
            }
        }

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if (asm_map[j].assemble())
                    triplets.push_back(Triplet<T>(asm_map[i], asm_map[j], lhs(i, j)));
                else
                    RHS(asm_map[i]) -= lhs(i, j) * dirichlet_data(j);
            }

            RHS(asm_map[i]) += rhs(i);
        }
    } // assemble()

    void
    assemble(const Mesh&                     msh,
             const typename Mesh::cell_type& cl,
             const boundary_type&            bnd,
             const matrix_type&              lhs,
             const vector_type&              rhs)
    {
        if (!use_bnd)
            throw std::invalid_argument("diffusion_assembler: you have to use boundary type in the constructor");

        const auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension - 1);
        const auto fcs = faces(msh, cl);

        std::vector<assembly_index> asm_map;
        asm_map.reserve(fcs.size() * fbs);

        vector_type dirichlet_data = vector_type::Zero(fcs.size() * fbs);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            const auto fc              = fcs[face_i];
            const auto face_offset     = priv::offset(msh, fc);
            const auto face_LHS_offset = compress_table.at(face_offset) * fbs;

            const auto face_id   = msh.lookup(fc);
            const bool dirichlet = bnd.is_dirichlet_face(face_id);

            for (size_t i = 0; i < fbs; i++)
                asm_map.push_back(assembly_index(face_LHS_offset + i, !dirichlet));

            if (dirichlet)
            {
                auto dirichlet_fun = bnd.dirichlet_boundary_func(face_id);

                dirichlet_data.block(face_i * fbs, 0, fbs, 1) =
                  project_function(msh, fc, di.face_degree(), dirichlet_fun, di.face_degree());
            }
        }

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if (asm_map[j].assemble())
                    triplets.push_back(Triplet<T>(asm_map[i], asm_map[j], lhs(i, j)));
                else
                    RHS(asm_map[i]) -= lhs(i, j) * dirichlet_data(j);
            }

            RHS(asm_map[i]) += rhs(i);
        }
    } // assemble()

    void
    impose_neumann_boundary_conditions(const Mesh& msh, const boundary_type& bnd)
    {
        if (!use_bnd)
            throw std::invalid_argument("diffusion_assembler: you have to use boundary type in the constructor");

        if (bnd.nb_faces_neumann() > 0)
        {
            for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
            {
                const auto bfc     = *itor;
                const auto face_id = msh.lookup(bfc);

                if (bnd.is_neumann_face(face_id))
                {
                    if (bnd.is_dirichlet_face(face_id))
                    {
                        throw std::invalid_argument("You tried to impose"
                                                    "both Dirichlet and Neumann conditions on the same face");
                    }
                    else if (bnd.is_robin_face(face_id))
                    {
                        throw std::invalid_argument("You tried to impose"
                                                    "both Robin and Neumann conditions on the same face");
                    }
                    else
                    {
                        const size_t                face_degree   = di.face_degree();
                        const size_t                num_face_dofs = scalar_basis_size(face_degree, Mesh::dimension - 1);
                        std::vector<assembly_index> asm_map;
                        asm_map.reserve(num_face_dofs);

                        auto face_offset     = face_id;
                        auto face_LHS_offset = compress_table.at(face_offset) * num_face_dofs;

                        for (size_t i = 0; i < num_face_dofs; i++)
                        {
                            asm_map.push_back(assembly_index(face_LHS_offset + i, true));
                        }

                        auto        fb      = make_scalar_monomial_basis(msh, bfc, face_degree);
                        vector_type neumann = make_rhs(msh, bfc, fb, bnd.neumann_boundary_func(face_id), face_degree);

                        assert(neumann.size() == num_face_dofs);
                        for (size_t i = 0; i < neumann.size(); i++)
                        {
                            RHS(asm_map[i]) += neumann[i];
                        }
                    }
                }
            }
        }
        else
            throw std::invalid_argument("There are no Neumann faces");
    }

    void
    impose_robin_boundary_conditions(const Mesh& msh, const boundary_type& bnd)
    {
        if (!use_bnd)
            throw std::invalid_argument("diffusion_assembler: you have to use boundary type in the constructor");
        if (bnd.nb_faces_robin() > 0)
        {
            for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
            {
                const auto bfc = *itor;
                const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
                if (!eid.first)
                    throw std::invalid_argument("This is a bug: face not found");
                const auto face_id = eid.second;

                if (bnd.is_robin_face(face_id))
                {
                    if (bnd.is_neumann_face(face_id))
                    {
                        switch (bnd.neumann_boundary_type(face_id))
                        {
                            case disk::NEUMANN:
                                throw std::invalid_argument(
                                  "You tried to impose both Neumann and Robin conditions on the same face");
                                break;
                            default: throw std::logic_error("Unknown Neumann Conditions"); break;
                        }
                    }
                    else if (bnd.is_dirichlet_face(face_id))
                    {
                        switch (bnd.dirichlet_boundary_type(face_id))
                        {
                            case disk::DIRICHLET:
                                throw std::invalid_argument(
                                  "You tried to impose both Dirichlet and Robin conditions on the same face");
                                break;
                            default: throw std::logic_error("Unknown Dirichlet Conditions"); break;
                        }
                    }
                    else
                    {
                        const size_t                face_degree   = di.face_degree();
                        const size_t                num_face_dofs = scalar_basis_size(face_degree, Mesh::dimension - 1);
                        std::vector<assembly_index> asm_map;
                        asm_map.reserve(num_face_dofs);

                        auto face_offset     = face_id;
                        auto face_LHS_offset = compress_table.at(face_offset) * num_face_dofs;

                        for (size_t i = 0; i < num_face_dofs; i++)
                            asm_map.push_back(assembly_index(face_LHS_offset + i, true));

                        auto        fb    = make_scalar_monomial_basis(msh, bfc, face_degree);
                        vector_type robin = make_rhs(msh, bfc, fb, bnd.robin_boundary_func(face_id), face_degree);
                        assert(robin.size() == num_face_dofs);

                        matrix_type mass = make_mass_matrix(msh, bfc, fb);

                        for (size_t i = 0; i < num_face_dofs; i++)
                        {
                            RHS(asm_map[i]) += robin[i];

                            for (size_t j = 0; j < num_face_dofs; j++)
                            {
                                if (asm_map[j].assemble())
                                    triplets.push_back(Triplet<T>(asm_map[i], asm_map[j], mass(i, j)));
                            }
                        }
                    }
                }
            }
        }
    }

    template<typename Function>
    vector_type
    take_local_data(const Mesh&                     msh,
                    const typename Mesh::cell_type& cl,
                    const vector_type&              solution,
                    const Function&                 dirichlet_bf)
    {
        const auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension - 1);
        const auto fcs = faces(msh, cl);

        const auto num_faces = fcs.size();

        vector_type ret = vector_type::Zero(num_faces * fbs);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            const auto fc = fcs[face_i];

            auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool { return msh.is_boundary(fc); };

            const bool dirichlet = is_dirichlet(fc);

            if (dirichlet)
            {
                ret.block(face_i * fbs, 0, fbs, 1) =
                  project_function(msh, fc, di.face_degree(), dirichlet_bf, di.face_degree());
            }
            else
            {
                const auto face_offset     = priv::offset(msh, fc);
                const auto face_SOL_offset = compress_table.at(face_offset) * fbs;

                ret.block(face_i * fbs, 0, fbs, 1) = solution.block(face_SOL_offset, 0, fbs, 1);
            }
        }

        return ret;
    }

    vector_type
    take_local_data(const Mesh&                     msh,
                    const typename Mesh::cell_type& cl,
                    const boundary_type&            bnd,
                    const vector_type&              solution)
    {
        const auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension - 1);
        const auto fcs = faces(msh, cl);

        const auto num_faces = fcs.size();

        vector_type ret = vector_type::Zero(num_faces * fbs);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            const auto fc = fcs[face_i];

            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");
            const auto face_id = eid.second;

            const bool dirichlet = bnd.is_dirichlet_face(face_id);

            if (dirichlet)
            {
                const auto dirichlet_bf = bnd.dirichlet_boundary_func(face_id);

                ret.block(face_i * fbs, 0, fbs, 1) =
                  project_function(msh, fc, di.face_degree(), dirichlet_bf, di.face_degree());
            }
            else
            {
                const auto face_offset             = priv::offset(msh, fc);
                const auto face_SOL_offset         = compress_table.at(face_offset) * fbs;
                ret.block(face_i * fbs, 0, fbs, 1) = solution.block(face_SOL_offset, 0, fbs, 1);
            }
        }

        return ret;
    }

    void
    finalize(void)
    {
        LHS.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();

        // dump_sparse_matrix(LHS, "diff.dat");
    }

    size_t
    num_assembled_faces() const
    {
        return num_other_faces;
    }
};

template<typename Mesh>
auto
make_diffusion_assembler(const Mesh& msh, const hho_degree_info& hdi)
{
    return diffusion_condensed_assembler<Mesh>(msh, hdi);
}

template<typename Mesh>
auto
make_diffusion_assembler(const Mesh& msh, const hho_degree_info& hdi, const scalar_boundary_conditions<Mesh>& bnd)
{
    return diffusion_condensed_assembler<Mesh>(msh, hdi, bnd);
}

template<typename Mesh>
class stokes_assembler
{
    using T = typename Mesh::coordinate_type;
    typedef disk::vector_boundary_conditions<Mesh> boundary_type;

    std::vector<size_t> compress_table;
    std::vector<size_t> expand_table;

    boundary_type           m_bnd;
    hho_degree_info         di;
    std::vector<Triplet<T>> triplets;

    size_t num_all_faces, num_dirichlet_faces, num_other_faces;
    size_t cbs_A, cbs_B, fbs_A;
    size_t system_size;

    class assembly_index
    {
        size_t idx;
        bool   assem;

      public:
        assembly_index(size_t i, bool as) : idx(i), assem(as) {}

        operator size_t() const
        {
            if (!assem)
                throw std::logic_error("Invalid assembly_index");

            return idx;
        }

        bool
        assemble() const
        {
            return assem;
        }

        friend std::ostream&
        operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };

  public:
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    SparseMatrix<T> LHS;
    vector_type     RHS;

    stokes_assembler(const Mesh& msh, const hho_degree_info& hdi, const boundary_type& bnd) : di(hdi), m_bnd(bnd)
    {
        auto is_dirichlet = [&](const typename Mesh::face& fc) -> bool {
            auto fc_id = msh.lookup(fc);
            return bnd.is_dirichlet_face(fc_id);
        };

        num_all_faces       = msh.faces_size();
        num_dirichlet_faces = std::count_if(msh.faces_begin(), msh.faces_end(), is_dirichlet);
        num_other_faces     = num_all_faces - num_dirichlet_faces;

        compress_table.resize(num_all_faces);
        expand_table.resize(num_other_faces);

        size_t compressed_offset = 0;
        for (size_t i = 0; i < num_all_faces; i++)
        {
            auto fc = *std::next(msh.faces_begin(), i);
            if (!is_dirichlet(fc))
            {
                compress_table.at(i)               = compressed_offset;
                expand_table.at(compressed_offset) = i;
                compressed_offset++;
            }
        }

        cbs_A = vector_basis_size(di.cell_degree(), Mesh::dimension, Mesh::dimension);
        fbs_A = vector_basis_size(di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        cbs_B = scalar_basis_size(di.face_degree(), Mesh::dimension);

        system_size = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * msh.cells_size() + 1;

        LHS = SparseMatrix<T>(system_size, system_size);
        RHS = vector_type::Zero(system_size);

        // testing Boundary module of Nicolas
        auto   num_face_dofs = vector_basis_size(di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        size_t face_dofs     = 0;
        for (size_t face_id = 0; face_id < msh.faces_size(); face_id++)
            face_dofs += num_face_dofs - m_bnd.dirichlet_imposed_dofs(face_id, di.face_degree());

        assert(face_dofs == fbs_A * num_other_faces);
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
    initialize()
    {
        LHS = SparseMatrix<T>(system_size, system_size);
        RHS = vector_type::Zero(system_size);
        return;
    }

    void
    assemble(const Mesh&                     msh,
             const typename Mesh::cell_type& cl,
             const matrix_type&              lhs_A,
             const matrix_type&              lhs_B,
             const vector_type&              rhs)
    {
        auto fcs = faces(msh, cl);

        std::vector<assembly_index> asm_map;
        asm_map.reserve(cbs_A + fcs.size() * fbs_A);

        auto cell_offset     = priv::offset(msh, cl);
        auto cell_LHS_offset = cell_offset * cbs_A;

        auto B_offset = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * cell_offset;

        for (size_t i = 0; i < cbs_A; i++)
            asm_map.push_back(assembly_index(cell_LHS_offset + i, true));

        vector_type dirichlet_data = vector_type::Zero(cbs_A + fcs.size() * fbs_A);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc              = fcs[face_i];
            auto face_offset     = priv::offset(msh, fc);
            auto face_LHS_offset = cbs_A * msh.cells_size() + compress_table.at(face_offset) * fbs_A;

            auto fc_id     = msh.lookup(fc);
            bool dirichlet = m_bnd.is_dirichlet_face(fc_id);

            for (size_t i = 0; i < fbs_A; i++)
                asm_map.push_back(assembly_index(face_LHS_offset + i, !dirichlet));

            if (dirichlet)
            {
                const auto face_id = msh.lookup(fc);

                auto dirichlet_fun = m_bnd.dirichlet_boundary_func(face_id);

                dirichlet_data.block(cbs_A + face_i * fbs_A, 0, fbs_A, 1) =
                  project_function(msh, fc, di.face_degree(), dirichlet_fun, di.face_degree());
            }
        }

        assert(asm_map.size() == lhs_A.rows() && asm_map.size() == lhs_A.cols());

        for (size_t i = 0; i < lhs_A.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs_A.cols(); j++)
            {
                if (asm_map[j].assemble())
                    triplets.push_back(Triplet<T>(asm_map[i], asm_map[j], lhs_A(i, j)));
                else
                    RHS(asm_map[i]) -= lhs_A(i, j) * dirichlet_data(j);
            }
        }

        for (size_t i = 0; i < lhs_B.rows(); i++)
        {
            for (size_t j = 0; j < lhs_B.cols(); j++)
            {
                auto global_i = B_offset + i;
                auto global_j = asm_map[j];
                if (asm_map[j].assemble())
                {
                    triplets.push_back(Triplet<T>(global_i, global_j, lhs_B(i, j)));
                    triplets.push_back(Triplet<T>(global_j, global_i, lhs_B(i, j)));
                }
                else
                    RHS(global_i) -= lhs_B(i, j) * dirichlet_data(j);
            }
        }

        auto        scalar_cell_basis = make_scalar_monomial_basis(msh, cl, di.face_degree());
        auto        qps               = integrate(msh, cl, di.face_degree());
        vector_type mult              = vector_type::Zero(scalar_cell_basis.size());
        for (auto& qp : qps)
        {
            auto phi = scalar_cell_basis.eval_functions(qp.point());
            mult += qp.weight() * phi;
        }
        auto mult_offset = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * msh.cells_size();

        for (size_t i = 0; i < mult.rows(); i++)
        {
            triplets.push_back(Triplet<T>(B_offset + i, mult_offset, mult(i)));
            triplets.push_back(Triplet<T>(mult_offset, B_offset + i, mult(i)));
        }

        RHS.block(cell_LHS_offset, 0, cbs_A, 1) += rhs.block(0, 0, cbs_A, 1);

    } // assemble()

    void
    assemble_alg(const Mesh&                     msh,
                 const typename Mesh::cell_type& cl,
                 const matrix_type&              lhs_A,
                 const matrix_type&              lhs_B,
                 const vector_type&              rhs)
    {
        auto fcs = faces(msh, cl);

        std::vector<assembly_index> asm_map;
        asm_map.reserve(cbs_A + fcs.size() * fbs_A);

        auto cell_offset     = priv::offset(msh, cl);
        auto cell_LHS_offset = cell_offset * cbs_A;

        auto B_offset = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * cell_offset;

        for (size_t i = 0; i < cbs_A; i++)
            asm_map.push_back(assembly_index(cell_LHS_offset + i, true));

        vector_type dirichlet_data = vector_type::Zero(cbs_A + fcs.size() * fbs_A);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc              = fcs[face_i];
            auto face_offset     = priv::offset(msh, fc);
            auto face_LHS_offset = cbs_A * msh.cells_size() + compress_table.at(face_offset) * fbs_A;

            auto fc_id     = msh.lookup(fc);
            bool dirichlet = m_bnd.is_dirichlet_face(fc_id);

            for (size_t i = 0; i < fbs_A; i++)
                asm_map.push_back(assembly_index(face_LHS_offset + i, !dirichlet));

            if (dirichlet)
            {
                const auto face_id = msh.lookup(fc);

                auto dirichlet_fun = m_bnd.dirichlet_boundary_func(face_id);

                dirichlet_data.block(cbs_A + face_i * fbs_A, 0, fbs_A, 1) =
                  project_function(msh, fc, di.face_degree(), dirichlet_fun, di.face_degree());
            }
        }

        assert(asm_map.size() == lhs_A.rows() && asm_map.size() == lhs_A.cols());
        assert(asm_map.size() == rhs.rows());

        for (size_t i = 0; i < lhs_A.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            RHS(asm_map[i]) += rhs(i);

            for (size_t j = 0; j < lhs_A.cols(); j++)
            {
                if (asm_map[j].assemble())
                    triplets.push_back(Triplet<T>(asm_map[i], asm_map[j], lhs_A(i, j)));
                else
                    RHS(asm_map[i]) -= lhs_A(i, j) * dirichlet_data(j);
            }
        }

        for (size_t i = 0; i < lhs_B.rows(); i++)
        {
            for (size_t j = 0; j < lhs_B.cols(); j++)
            {
                auto global_i = B_offset + i;
                auto global_j = asm_map[j];
                if (asm_map[j].assemble())
                {
                    triplets.push_back(Triplet<T>(global_i, global_j, lhs_B(i, j)));
                    triplets.push_back(Triplet<T>(global_j, global_i, lhs_B(i, j)));
                }
                else
                    RHS(global_i) -= lhs_B(i, j) * dirichlet_data(j);
            }
        }

        auto        scalar_cell_basis = make_scalar_monomial_basis(msh, cl, di.face_degree());
        auto        qps               = integrate(msh, cl, di.face_degree());
        vector_type mult              = vector_type::Zero(scalar_cell_basis.size());
        for (auto& qp : qps)
        {
            auto phi = scalar_cell_basis.eval_functions(qp.point());
            mult += qp.weight() * phi;
        }
        auto mult_offset = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * msh.cells_size();

        for (size_t i = 0; i < mult.rows(); i++)
        {
            triplets.push_back(Triplet<T>(B_offset + i, mult_offset, mult(i)));
            triplets.push_back(Triplet<T>(mult_offset, B_offset + i, mult(i)));
        }

        // RHS.block(cell_LHS_offset, 0, cbs_A, 1) += rhs.block(0, 0, cbs_A, 1);

    } // assemble_alg()

    void
    finalize(void)
    {
        LHS.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();
    }
    size_t
    num_assembled_faces() const
    {
        return num_other_faces;
    }

    size_t
    global_system_size() const
    {
        return system_size;
    }

    vector_type
    take_velocity(const Mesh& msh, const typename Mesh::cell_type& cl, const vector_type& sol) const
    {
        auto num_faces = howmany_faces(msh, cl);
        auto dim       = Mesh::dimension;
        auto cell_ofs  = priv::offset(msh, cl);

        vector_type svel(cbs_A + num_faces * fbs_A);
        svel.block(0, 0, cbs_A, 1) = sol.block(cell_ofs * cbs_A, 0, cbs_A, 1);
        auto fcs                   = faces(msh, cl);
        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto       fc      = fcs[i];
            const auto face_id = msh.lookup(fc);

            if (m_bnd.is_dirichlet_face(face_id))
            {
                auto velocity = m_bnd.dirichlet_boundary_func(face_id);

                svel.block(cbs_A + i * fbs_A, 0, fbs_A, 1) =
                  project_function(msh, fc, di.face_degree(), velocity, di.face_degree());
            }
            else
            {
                auto face_ofs   = priv::offset(msh, fc);
                auto global_ofs = cbs_A * msh.cells_size() + compress_table.at(face_ofs) * fbs_A;
                svel.block(cbs_A + i * fbs_A, 0, fbs_A, 1) = sol.block(global_ofs, 0, fbs_A, 1);
            }
        }
        return svel;
    }

    vector_type
    take_pressure(const Mesh& msh, const typename Mesh::cell_type& cl, const vector_type& sol) const
    {
        auto cell_ofs = priv::offset(msh, cl);
        auto pres_ofs = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * cell_ofs;

        vector_type spres = sol.block(pres_ofs, 0, cbs_B, 1);
        return spres;
    }

    auto
    global_face_offset(const Mesh& msh, const typename Mesh::face_type& fc)
    {
        auto cbs_A = vector_basis_size(di.cell_degree(), Mesh::dimension, Mesh::dimension);
        auto fbs_A = vector_basis_size(di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        auto cbs_B = scalar_basis_size(di.face_degree(), Mesh::dimension);

        auto face_offset = priv::offset(msh, fc);
        return cbs_A * msh.cells_size() + compress_table.at(face_offset) * fbs_A;
    }

    auto
    global_face_offset(const Mesh& msh, const typename Mesh::face_type& fc) const
    {
        auto cbs_A = vector_basis_size(di.cell_degree(), Mesh::dimension, Mesh::dimension);
        auto fbs_A = vector_basis_size(di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        auto cbs_B = scalar_basis_size(di.face_degree(), Mesh::dimension);

        auto face_offset = priv::offset(msh, fc);
        return cbs_A * msh.cells_size() + compress_table.at(face_offset) * fbs_A;
    }
};

template<typename Mesh, typename BoundaryType>
auto
make_stokes_assembler(const Mesh& msh, hho_degree_info hdi, const BoundaryType& bnd)
{
    return stokes_assembler<Mesh>(msh, hdi, bnd);
}

template<typename Mesh>
class stokes_assembler_alg
{
    using T = typename Mesh::coordinate_type;
    typedef disk::vector_boundary_conditions<Mesh> boundary_type;

    std::vector<size_t> compress_table;
    std::vector<size_t> expand_table;

    boundary_type           m_bnd;
    hho_degree_info         di;
    std::vector<Triplet<T>> triplets;
    Matrix<T, Dynamic, 1>   RHS_DIRICHLET;

    size_t num_all_faces, num_dirichlet_faces, num_other_faces;
    size_t cbs_A, cbs_B, fbs_A;
    size_t system_size;

    class assembly_index
    {
        size_t idx;
        bool   assem;

      public:
        assembly_index(size_t i, bool as) : idx(i), assem(as) {}

        operator size_t() const
        {
            if (!assem)
                throw std::logic_error("Invalid assembly_index");

            return idx;
        }

        bool
        assemble() const
        {
            return assem;
        }

        friend std::ostream&
        operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };

  public:
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    SparseMatrix<T> LHS;
    vector_type     RHS;

    stokes_assembler_alg(const Mesh& msh, const hho_degree_info& hdi, const boundary_type& bnd) : di(hdi), m_bnd(bnd)
    {
        auto is_dirichlet = [&](const typename Mesh::face& fc) -> bool {
            auto fc_id = msh.lookup(fc);
            return bnd.is_dirichlet_face(fc_id);
        };

        num_all_faces       = msh.faces_size();
        num_dirichlet_faces = std::count_if(msh.faces_begin(), msh.faces_end(), is_dirichlet);
        num_other_faces     = num_all_faces - num_dirichlet_faces;

        compress_table.resize(num_all_faces);
        expand_table.resize(num_other_faces);

        size_t compressed_offset = 0;
        for (size_t i = 0; i < num_all_faces; i++)
        {
            auto fc = *std::next(msh.faces_begin(), i);
            if (!is_dirichlet(fc))
            {
                compress_table.at(i)               = compressed_offset;
                expand_table.at(compressed_offset) = i;
                compressed_offset++;
            }
        }

        cbs_A = vector_basis_size(di.cell_degree(), Mesh::dimension, Mesh::dimension);
        fbs_A = vector_basis_size(di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        cbs_B = scalar_basis_size(di.face_degree(), Mesh::dimension);

        system_size = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * msh.cells_size() + 1;

        LHS           = SparseMatrix<T>(system_size, system_size);
        RHS           = vector_type::Zero(system_size);
        RHS_DIRICHLET = vector_type::Zero(system_size);

        // testing Boundary module of Nicolas
        auto   num_face_dofs = vector_basis_size(di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        size_t face_dofs     = 0;
        for (size_t face_id = 0; face_id < msh.faces_size(); face_id++)
            face_dofs += num_face_dofs - m_bnd.dirichlet_imposed_dofs(face_id, di.face_degree());

        assert(face_dofs == fbs_A * num_other_faces);
    }

    void
    initialize_lhs()
    {
        LHS           = SparseMatrix<T>(system_size, system_size);
        RHS_DIRICHLET = vector_type::Zero(system_size);
        return;
    }

    void
    initialize_rhs()
    {
        RHS = vector_type::Zero(system_size);
        return;
    }

    void
    assemble_lhs(const Mesh&                     msh,
                 const typename Mesh::cell_type& cl,
                 const matrix_type&              lhs_A,
                 const matrix_type&              lhs_B)
    {
        auto                        fcs = faces(msh, cl);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(cbs_A + fcs.size() * fbs_A);

        auto cell_offset     = priv::offset(msh, cl);
        auto cell_LHS_offset = cell_offset * cbs_A;

        auto B_offset = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * cell_offset;

        for (size_t i = 0; i < cbs_A; i++)
            asm_map.push_back(assembly_index(cell_LHS_offset + i, true));

        vector_type dirichlet_data = vector_type::Zero(cbs_A + fcs.size() * fbs_A);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc              = fcs[face_i];
            auto face_offset     = priv::offset(msh, fc);
            auto face_LHS_offset = cbs_A * msh.cells_size() + compress_table.at(face_offset) * fbs_A;

            auto fc_id     = msh.lookup(fc);
            bool dirichlet = m_bnd.is_dirichlet_face(fc_id);

            for (size_t i = 0; i < fbs_A; i++)
                asm_map.push_back(assembly_index(face_LHS_offset + i, !dirichlet));

            if (dirichlet)
            {
                auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
                if (!eid.first)
                    throw std::invalid_argument("This is a bug: face not found");

                const auto face_id = eid.second;

                auto dirichlet_fun = m_bnd.dirichlet_boundary_func(face_id);

                dirichlet_data.block(cbs_A + face_i * fbs_A, 0, fbs_A, 1) =
                  project_function(msh, fc, di.face_degree(), dirichlet_fun, di.face_degree());
            }
        }

        assert(asm_map.size() == lhs_A.rows() && asm_map.size() == lhs_A.cols());

        for (size_t i = 0; i < lhs_A.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs_A.cols(); j++)
            {
                if (asm_map[j].assemble())
                    triplets.push_back(Triplet<T>(asm_map[i], asm_map[j], lhs_A(i, j)));
                else
                    RHS_DIRICHLET(asm_map[i]) -= lhs_A(i, j) * dirichlet_data(j);
            }
        }

        for (size_t i = 0; i < lhs_B.rows(); i++)
        {
            for (size_t j = 0; j < lhs_B.cols(); j++)
            {
                auto global_i = B_offset + i;
                auto global_j = asm_map[j];
                if (asm_map[j].assemble())
                {
                    triplets.push_back(Triplet<T>(global_i, global_j, lhs_B(i, j)));
                    triplets.push_back(Triplet<T>(global_j, global_i, lhs_B(i, j)));
                }
                else
                    RHS_DIRICHLET(global_i) -= lhs_B(i, j) * dirichlet_data(j);
            }
        }

        auto        scalar_cell_basis = make_scalar_monomial_basis(msh, cl, di.face_degree());
        auto        qps               = integrate(msh, cl, di.face_degree());
        vector_type mult              = vector_type::Zero(scalar_cell_basis.size());
        for (auto& qp : qps)
        {
            auto phi = scalar_cell_basis.eval_functions(qp.point());
            mult += qp.weight() * phi;
        }
        auto mult_offset = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * msh.cells_size();

        for (size_t i = 0; i < mult.rows(); i++)
        {
            triplets.push_back(Triplet<T>(B_offset + i, mult_offset, mult(i)));
            triplets.push_back(Triplet<T>(mult_offset, B_offset + i, mult(i)));
        }

        // RHS.block(cell_LHS_offset, 0, cbs_A, 1) += rhs.block(0, 0, cbs_A, 1);

    } // assemble_alg()
    void
    assemble_rhs(const Mesh& msh, const typename Mesh::cell_type& cl, const vector_type& rhs)
    {
        auto                        fcs = faces(msh, cl);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(cbs_A + fcs.size() * fbs_A);

        auto cell_offset     = priv::offset(msh, cl);
        auto cell_LHS_offset = cell_offset * cbs_A;

        auto B_offset = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * cell_offset;

        for (size_t i = 0; i < cbs_A; i++)
            asm_map.push_back(assembly_index(cell_LHS_offset + i, true));

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc              = fcs[face_i];
            auto face_offset     = priv::offset(msh, fc);
            auto face_LHS_offset = cbs_A * msh.cells_size() + compress_table.at(face_offset) * fbs_A;

            auto fc_id     = msh.lookup(fc);
            bool dirichlet = m_bnd.is_dirichlet_face(fc_id);

            for (size_t i = 0; i < fbs_A; i++)
                asm_map.push_back(assembly_index(face_LHS_offset + i, !dirichlet));
        }

        assert(asm_map.size() == rhs.rows());

        for (size_t i = 0; i < rhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            RHS(asm_map[i]) += rhs(i);
        }
    } // assemble_alg()

    void
    finalize_lhs(void)
    {
        LHS.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();
    }

    void
    finalize_rhs(void)
    {
        RHS += RHS_DIRICHLET;
    }

    size_t
    global_system_size() const
    {
        return system_size;
    }

    vector_type
    take_velocity(const Mesh& msh, const typename Mesh::cell_type& cl, const vector_type& sol) const
    {
        auto num_faces = howmany_faces(msh, cl);
        auto dim       = Mesh::dimension;
        auto cell_ofs  = priv::offset(msh, cl);

        vector_type svel(cbs_A + num_faces * fbs_A);
        svel.block(0, 0, cbs_A, 1) = sol.block(cell_ofs * cbs_A, 0, cbs_A, 1);
        auto fcs                   = faces(msh, cl);
        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto       fc      = fcs[i];
            const auto face_id = msh.lookup(fc);

            if (m_bnd.is_dirichlet_face(face_id))
            {
                auto velocity = m_bnd.dirichlet_boundary_func(face_id);

                svel.block(cbs_A + i * fbs_A, 0, fbs_A, 1) =
                  project_function(msh, fc, di.face_degree(), velocity, di.face_degree());
            }
            else
            {
                auto face_ofs   = priv::offset(msh, fc);
                auto global_ofs = cbs_A * msh.cells_size() + compress_table.at(face_ofs) * fbs_A;
                svel.block(cbs_A + i * fbs_A, 0, fbs_A, 1) = sol.block(global_ofs, 0, fbs_A, 1);
            }
        }
        return svel;
    }

    vector_type
    take_pressure(const Mesh& msh, const typename Mesh::cell_type& cl, const vector_type& sol) const
    {
        auto cell_ofs = priv::offset(msh, cl);
        auto pres_ofs = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * cell_ofs;

        vector_type spres = sol.block(pres_ofs, 0, cbs_B, 1);
        return spres;
    }
};

template<typename Mesh, typename BoundaryType>
auto
make_stokes_assembler_alg(const Mesh& msh, hho_degree_info hdi, const BoundaryType& bnd)
{
    return stokes_assembler_alg<Mesh>(msh, hdi, bnd);
}

// assembler for vector primal problem with HHO like vector_laplacian problem and mechanics
template<typename Mesh>
class assembler_mechanics
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    typedef disk::vector_boundary_conditions<Mesh> bnd_type;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;
    typedef sparse_matrix<scalar_type>  sparse_type;
    typedef triplet<scalar_type>        triplet_type;

    const static size_t dimension = mesh_type::dimension;

    std::vector<triplet_type> m_triplets;
    size_t                    m_num_unknowns;
    std::vector<size_t>       face_compress_map;
    hho_degree_info           m_hdi;

  public:
    sparse_type LHS;
    vector_type RHS;

    assembler_mechanics() {}

    assembler_mechanics(const mesh_type& msh, const hho_degree_info& hdi, const bnd_type& bnd) : m_hdi(hdi)
    {
        const auto num_face_dofs = vector_basis_size(m_hdi.face_degree(), dimension - 1, dimension);

        face_compress_map.resize(msh.faces_size());

        size_t total_dofs = 0;
        for (size_t face_id = 0; face_id < msh.faces_size(); face_id++)
        {

            face_compress_map.at(face_id) = total_dofs;
            const auto free_dofs          = num_face_dofs - bnd.dirichlet_imposed_dofs(face_id, m_hdi.face_degree());
            total_dofs += free_dofs;
        }
        m_num_unknowns = total_dofs;
        LHS            = sparse_type(m_num_unknowns, m_num_unknowns);
        RHS            = vector_type::Zero(m_num_unknowns);
    }

    // don't forget to reset RHS at each Newton iteration
    void
    setZeroRhs()
    {
        RHS.setZero();
    }

    template<typename LocalContrib>
    void
    assemble(const mesh_type& msh, const cell_type& cl, const bnd_type& bnd, const LocalContrib& lc, int di = 0)
    {
        const size_t      face_degree   = m_hdi.face_degree();
        const auto        num_face_dofs = vector_basis_size(face_degree, dimension - 1, dimension);
        const scalar_type zero          = 0;

        const auto          fcs = faces(msh, cl);
        std::vector<size_t> l2g(fcs.size() * num_face_dofs);
        std::vector<bool>   l2l(fcs.size() * num_face_dofs, true);
        vector_type         rhs_bc = vector_type::Zero(fcs.size() * num_face_dofs);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            const auto fc                       = fcs[face_i];
            const auto face_id                  = msh.lookup(fc);
            const bool fc_is_dirichlet_boundary = bnd.is_dirichlet_face(face_id);
            const auto face_offset              = face_compress_map.at(face_id);
            const auto pos                      = face_i * num_face_dofs;

            if (!fc_is_dirichlet_boundary)
            {
                for (size_t i = 0; i < num_face_dofs; i++)
                {
                    l2g.at(pos + i) = face_offset + i;
                }
            }
            else
            {
                size_t ind_sol = 0;

                vector_type proj_bcf = project_function(msh, fc, face_degree, bnd.dirichlet_boundary_func(face_id), di);

                bool ind_ok = false;
                for (size_t face_j = 0; face_j < fcs.size(); face_j++)
                {
                    const auto fcj  = fcs[face_j];
                    auto       eidj = find_element_id(msh.faces_begin(), msh.faces_end(), fcj);
                    if (!eidj.first)
                        throw std::invalid_argument("This is a bug: face not found");

                    matrix_type mat_Fj = lc.first.block(face_j * num_face_dofs, pos, num_face_dofs, num_face_dofs);

                    switch (bnd.dirichlet_boundary_type(face_id))
                    {
                        case disk::DIRICHLET:
                        {
                            if (!ind_ok)
                            {
                                for (size_t i = 0; i < num_face_dofs; i++)
                                {
                                    l2g.at(pos + i) = 0xDEADBEEF;
                                    l2l.at(pos + i) = false;
                                }
                                ind_ok = true;
                            }
                            break;
                        }
                        case disk::CLAMPED:
                        {
                            proj_bcf.setZero();
                            mat_Fj.setZero();
                            if (!ind_ok)
                            {
                                for (size_t i = 0; i < num_face_dofs; i++)
                                {
                                    l2g.at(pos + i) = 0xDEADBEEF;
                                    l2l.at(pos + i) = false;
                                }
                                ind_ok = true;
                            }
                            break;
                        }
                        case disk::DX:
                        {
                            for (size_t i = 0; i < num_face_dofs; i += dimension)
                            {
                                mat_Fj.col(i + 1).setZero();
                                proj_bcf(i + 1) = zero;
                                if (dimension == 3)
                                {
                                    mat_Fj.col(i + 2).setZero();
                                    proj_bcf(i + 2) = zero;
                                }
                                if (!ind_ok)
                                {
                                    l2g.at(pos + i)     = 0xDEADBEEF;
                                    l2l.at(pos + i)     = false;
                                    l2g.at(pos + i + 1) = face_offset + ind_sol++;
                                    if (dimension == 3)
                                    {
                                        l2g.at(pos + i + 2) = face_offset + ind_sol++;
                                    }
                                }
                            }
                            ind_ok = true;
                            break;
                        }
                        case disk::DY:
                        {
                            for (size_t i = 0; i < num_face_dofs; i += dimension)
                            {
                                mat_Fj.col(i).setZero();
                                proj_bcf(i) = zero;
                                if (dimension == 3)
                                {
                                    mat_Fj.col(i + 2).setZero();
                                    proj_bcf(i + 2) = zero;
                                }
                                if (!ind_ok)
                                {
                                    l2g.at(pos + i)     = face_offset + ind_sol++;
                                    l2g.at(pos + i + 1) = 0xDEADBEEF;
                                    l2l.at(pos + i + 1) = false;
                                    if (dimension == 3)
                                    {
                                        l2g.at(pos + i + 2) = face_offset + ind_sol++;
                                    }
                                }
                            }
                            ind_ok = true;
                            break;
                        }
                        case disk::DZ:
                        {
                            if (dimension != 3)
                                throw std::invalid_argument("You are not in 3D");
                            for (size_t i = 0; i < num_face_dofs; i += dimension)
                            {
                                mat_Fj.col(i).setZero();
                                proj_bcf(i) = zero;
                                mat_Fj.col(i + 1).setZero();
                                proj_bcf(i + 1) = zero;
                                if (!ind_ok)
                                {
                                    l2g.at(pos + i)     = face_offset + ind_sol++;
                                    l2g.at(pos + i + 1) = face_offset + ind_sol++;
                                    l2g.at(pos + i + 2) = 0xDEADBEEF;
                                    l2l.at(pos + i + 2) = false;
                                }
                            }
                            ind_ok = true;
                            break;
                        }
                        case disk::DXDY:
                        {
                            for (size_t i = 0; i < num_face_dofs; i += dimension)
                            {
                                if (dimension == 3)
                                {
                                    mat_Fj.col(i + 2).setZero();
                                    proj_bcf(i + 2) = zero;
                                }
                                if (!ind_ok)
                                {
                                    l2g.at(pos + i)     = 0xDEADBEEF;
                                    l2g.at(pos + i + 1) = 0xDEADBEEF;
                                    l2l.at(pos + i)     = false;
                                    l2l.at(pos + i + 1) = false;
                                    if (dimension == 3)
                                    {
                                        l2g.at(pos + i + 2) = face_offset + ind_sol++;
                                    }
                                }
                            }
                            ind_ok = true;
                            break;
                        }
                        case disk::DXDZ:
                        {
                            if (dimension != 3)
                                throw std::invalid_argument("You are not in 3D");
                            for (size_t i = 0; i < num_face_dofs; i += dimension)
                            {
                                mat_Fj.col(i + 1).setZero();
                                proj_bcf(i + 1) = zero;
                                if (!ind_ok)
                                {
                                    l2g.at(pos + i)     = 0xDEADBEEF;
                                    l2l.at(pos + i)     = false;
                                    l2g.at(pos + i + 1) = face_offset + ind_sol++;
                                    l2g.at(pos + i + 2) = 0xDEADBEEF;
                                    l2l.at(pos + i + 2) = false;
                                }
                            }
                            ind_ok = true;
                            break;
                        }
                        case disk::DYDZ:
                        {
                            if (dimension != 3)
                                throw std::invalid_argument("You are not in 3D");
                            for (size_t i = 0; i < num_face_dofs; i += dimension)
                            {
                                mat_Fj.col(i).setZero();
                                proj_bcf(i) = zero;
                                if (!ind_ok)
                                {
                                    l2g.at(pos + i)     = face_offset + ind_sol++;
                                    l2g.at(pos + i + 1) = 0xDEADBEEF;
                                    l2g.at(pos + i + 2) = 0xDEADBEEF;
                                    l2l.at(pos + i + 1) = false;
                                    l2l.at(pos + i + 2) = false;
                                }
                            }
                            ind_ok = true;
                            break;
                        }
                        default:
                        {
                            throw std::logic_error("Unknown Dirichlet Conditions (assembler)");
                            break;
                        }
                    }

                    rhs_bc.segment(face_j * num_face_dofs, num_face_dofs) += mat_Fj * proj_bcf;
                }
            }
        }
        assert(lc.first.rows() == lc.first.cols());
        assert(lc.first.rows() == lc.second.size());
        assert(lc.second.size() == l2g.size());
        assert(lc.second.size() == rhs_bc.size());

#ifdef FILL_COLMAJOR
        for (size_t j = 0; j < lc.first.cols(); j++)
        {
            if (l2l[j])
            {
                for (size_t i = 0; i < lc.first.rows(); i++)
                {
                    if (l2l[i])
                    {
                        m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));
                    }
                }
                RHS(l2g.at(j)) += lc.second(j) - rhs_bc(j);
            }
        }
#else
        for (size_t i = 0; i < lc.first.rows(); i++)
        {
            if (l2l[i])
            {
                for (size_t j = 0; j < lc.first.cols(); j++)
                {
                    if (l2l[j])
                    {
                        m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));
                    }
                }
                RHS(l2g.at(i)) += lc.second(i) - rhs_bc(i);
            }
        }
#endif
    }

    vector_type
    take_local_data(const Mesh&                     msh,
                    const typename Mesh::cell_type& cl,
                    const bnd_type&                 bnd,
                    const vector_type&              solution,
                    size_t                          di = 0)
    {
        const auto face_degree   = m_hdi.face_degree();
        const auto num_face_dofs = vector_basis_size(face_degree, dimension - 1, dimension);

        const auto fcs       = faces(msh, cl);
        const auto num_faces = fcs.size();

        vector_type ret = vector_type::Zero(num_face_dofs * num_faces);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            const auto fc              = fcs[face_i];
            const auto face_id         = msh.lookup(fc);
            const auto face_offset     = face_i * num_face_dofs;
            const auto compress_offset = face_compress_map.at(face_id);

            if (bnd.is_dirichlet_face(face_id))
            {
                size_t sol_ind = 0;

                const vector_type proj_bcf =
                  project_function(msh, fc, face_degree, bnd.dirichlet_boundary_func(face_id), di);

                assert(proj_bcf.size() == num_face_dofs);

                switch (bnd.dirichlet_boundary_type(face_id))
                {
                    case disk::DIRICHLET:
                    {
                        ret.segment(face_offset, num_face_dofs) = proj_bcf;
                        break;
                    }
                    case disk::CLAMPED:
                    {
                        ret.segment(face_offset, num_face_dofs).setZero();
                        break;
                    }
                    case disk::DX:
                    {

                        for (size_t i = 0; i < num_face_dofs; i += dimension)
                        {
                            ret(face_offset + i)     = proj_bcf(i);
                            ret(face_offset + i + 1) = solution(compress_offset + sol_ind++);
                            if (dimension == 3)
                            {
                                ret(face_offset + i + 2) = solution(compress_offset + sol_ind++);
                            }
                        }
                        break;
                    }
                    case disk::DY:
                    {
                        for (size_t i = 0; i < num_face_dofs; i += dimension)
                        {
                            ret(face_offset + i)     = solution(compress_offset + sol_ind++);
                            ret(face_offset + i + 1) = proj_bcf(i + 1);
                            if (dimension == 3)
                            {
                                ret(face_offset + i + 2) = solution(compress_offset + sol_ind++);
                            }
                        }
                        break;
                    }
                    case disk::DZ:
                    {
                        if (dimension != 3)
                            throw std::invalid_argument("You are not in 3D");
                        for (size_t i = 0; i < num_face_dofs; i += dimension)
                        {
                            ret(face_offset + i)     = solution(compress_offset + sol_ind++);
                            ret(face_offset + i + 1) = solution(compress_offset + sol_ind++);
                            ret(face_offset + i + 2) = proj_bcf(i + 2);
                        }
                        break;
                    }
                    case disk::DXDY:
                    {
                        for (size_t i = 0; i < num_face_dofs; i += dimension)
                        {
                            ret(face_offset + i)     = proj_bcf(i);
                            ret(face_offset + i + 1) = proj_bcf(i + 1);
                            if (dimension == 3)
                            {
                                ret(face_offset + i + 2) = solution(compress_offset + sol_ind++);
                            }
                        }
                        break;
                    }
                    case disk::DXDZ:
                    {
                        if (dimension != 3)
                            throw std::invalid_argument("You are not in 3D");
                        for (size_t i = 0; i < num_face_dofs; i += dimension)
                        {
                            ret(face_offset + i)     = proj_bcf(i);
                            ret(face_offset + i + 1) = solution(compress_offset + sol_ind++);
                            ret(face_offset + i + 2) = proj_bcf(i + 2);
                        }
                        break;
                    }
                    case disk::DYDZ:
                    {
                        if (dimension != 3)
                            throw std::invalid_argument("You are not in 3D");
                        for (size_t i = 0; i < num_face_dofs; i += dimension)
                        {
                            ret(face_offset + i)     = solution(compress_offset + sol_ind++);
                            ret(face_offset + i + 1) = proj_bcf(i + 1);
                            ret(face_offset + i + 2) = proj_bcf(i + 2);
                        }
                        break;
                    }
                    default:
                    {
                        throw std::logic_error("Unknown Dirichlet Conditions (assembler)");
                        break;
                    }
                }
            }
            else
            {
                ret.segment(face_offset, num_face_dofs) = solution.segment(compress_offset, num_face_dofs);
            }
        }

        return ret;
    }

    vector_type
    expand_solution(const mesh_type& msh, const bnd_type& bnd, const vector_type& solution, int di = 0)
    {
        assert(solution.size() == m_num_unknowns);
        const auto face_degree   = m_hdi.face_degree();
        const auto num_face_dofs = vector_basis_size(face_degree, dimension - 1, dimension);

        vector_type ret = vector_type::Zero(num_face_dofs * msh.faces_size());

        for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
        {
            const auto bfc             = *itor;
            const auto face_id         = msh.lookup(bfc);
            const auto face_offset     = face_id * num_face_dofs;
            const auto compress_offset = face_compress_map.at(face_id);

            if (bnd.is_dirichlet_face(face_id))
            {
                size_t sol_ind = 0;

                const vector_type proj_bcf =
                  project_function(msh, bfc, face_degree, bnd.dirichlet_boundary_func(face_id), di);

                assert(proj_bcf.size() == num_face_dofs);

                switch (bnd.dirichlet_boundary_type(face_id))
                {
                    case disk::DIRICHLET:
                    {
                        ret.segment(face_offset, num_face_dofs) = proj_bcf;
                        break;
                    }
                    case disk::CLAMPED:
                    {
                        ret.segment(face_offset, num_face_dofs).setZero();
                        break;
                    }
                    case disk::DX:
                    {

                        for (size_t i = 0; i < num_face_dofs; i += dimension)
                        {
                            ret(face_offset + i)     = proj_bcf(i);
                            ret(face_offset + i + 1) = solution(compress_offset + sol_ind++);
                            if (dimension == 3)
                            {
                                ret(face_offset + i + 2) = solution(compress_offset + sol_ind++);
                            }
                        }
                        break;
                    }
                    case disk::DY:
                    {
                        for (size_t i = 0; i < num_face_dofs; i += dimension)
                        {
                            ret(face_offset + i)     = solution(compress_offset + sol_ind++);
                            ret(face_offset + i + 1) = proj_bcf(i + 1);
                            if (dimension == 3)
                            {
                                ret(face_offset + i + 2) = solution(compress_offset + sol_ind++);
                            }
                        }
                        break;
                    }
                    case disk::DZ:
                    {
                        if (dimension != 3)
                            throw std::invalid_argument("You are not in 3D");
                        for (size_t i = 0; i < num_face_dofs; i += dimension)
                        {
                            ret(face_offset + i)     = solution(compress_offset + sol_ind++);
                            ret(face_offset + i + 1) = solution(compress_offset + sol_ind++);
                            ret(face_offset + i + 2) = proj_bcf(i + 2);
                        }
                        break;
                    }
                    case disk::DXDY:
                    {
                        for (size_t i = 0; i < num_face_dofs; i += dimension)
                        {
                            ret(face_offset + i)     = proj_bcf(i);
                            ret(face_offset + i + 1) = proj_bcf(i + 1);
                            if (dimension == 3)
                            {
                                ret(face_offset + i + 2) = solution(compress_offset + sol_ind++);
                            }
                        }
                        break;
                    }
                    case disk::DXDZ:
                    {
                        if (dimension != 3)
                            throw std::invalid_argument("You are not in 3D");
                        for (size_t i = 0; i < num_face_dofs; i += dimension)
                        {
                            ret(face_offset + i)     = proj_bcf(i);
                            ret(face_offset + i + 1) = solution(compress_offset + sol_ind++);
                            ret(face_offset + i + 2) = proj_bcf(i + 2);
                        }
                        break;
                    }
                    case disk::DYDZ:
                    {
                        if (dimension != 3)
                            throw std::invalid_argument("You are not in 3D");
                        for (size_t i = 0; i < num_face_dofs; i += dimension)
                        {
                            ret(face_offset + i)     = solution(compress_offset + sol_ind++);
                            ret(face_offset + i + 1) = proj_bcf(i + 1);
                            ret(face_offset + i + 2) = proj_bcf(i + 2);
                        }
                        break;
                    }
                    default:
                    {
                        throw std::logic_error("Unknown Dirichlet Conditions (assembler)");
                        break;
                    }
                }
            }
            else
            {
                ret.segment(face_offset, num_face_dofs) = solution.segment(compress_offset, num_face_dofs);
            }
        }

        return ret;
    }

    template<typename LocalContrib>
    void
    assemble_nl(const mesh_type&                msh,
                const cell_type&                cl,
                const bnd_type&                 bnd,
                const LocalContrib&             lc,
                const std::vector<vector_type>& sol_F,
                int                             di = 0)
    {
        assert(sol_F.size() == msh.faces_size());
        const size_t      face_degree   = m_hdi.face_degree();
        const auto        num_face_dofs = vector_basis_size(face_degree, dimension - 1, dimension);
        const scalar_type zero          = 0;

        const auto          fcs = faces(msh, cl);
        std::vector<size_t> l2g(fcs.size() * num_face_dofs);
        std::vector<bool>   l2l(fcs.size() * num_face_dofs, true);
        vector_type         rhs_bc = vector_type::Zero(fcs.size() * num_face_dofs);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            const auto fc                       = fcs[face_i];
            const auto face_id                  = msh.lookup(fc);
            const bool fc_is_dirichlet_boundary = bnd.is_dirichlet_face(face_id);
            const auto face_offset              = face_compress_map.at(face_id);
            const auto pos                      = face_i * num_face_dofs;

            if (!fc_is_dirichlet_boundary)
            {
                for (size_t i = 0; i < num_face_dofs; i++)
                {
                    l2g.at(pos + i) = face_offset + i;
                }
            }
            else
            {
                size_t ind_sol = 0;

                const vector_type proj_bcf =
                  project_function(msh, fc, face_degree, bnd.dirichlet_boundary_func(face_id), di);
                assert(proj_bcf.size() == sol_F[face_id].size());

                vector_type incr   = proj_bcf - sol_F[face_id];
                bool        ind_ok = false;
                for (size_t face_j = 0; face_j < fcs.size(); face_j++)
                {
                    const auto fcj  = fcs[face_j];
                    auto       eidj = find_element_id(msh.faces_begin(), msh.faces_end(), fcj);
                    if (!eidj.first)
                        throw std::invalid_argument("This is a bug: face not found");

                    matrix_type mat_Fj = lc.first.block(face_j * num_face_dofs, pos, num_face_dofs, num_face_dofs);

                    switch (bnd.dirichlet_boundary_type(face_id))
                    {
                        case disk::DIRICHLET:
                        {
                            if (!ind_ok)
                            {
                                for (size_t i = 0; i < num_face_dofs; i++)
                                {
                                    l2g.at(pos + i) = 0xDEADBEEF;
                                    l2l.at(pos + i) = false;
                                }
                                ind_ok = true;
                            }
                            break;
                        }
                        case disk::CLAMPED:
                        {
                            incr = -sol_F[face_id];
                            if (!ind_ok)
                            {
                                for (size_t i = 0; i < num_face_dofs; i++)
                                {
                                    l2g.at(pos + i) = 0xDEADBEEF;
                                    l2l.at(pos + i) = false;
                                }
                                ind_ok = true;
                            }
                            break;
                        }
                        case disk::DX:
                        {
                            for (size_t i = 0; i < num_face_dofs; i += dimension)
                            {
                                mat_Fj.col(i + 1).setZero();
                                incr(i + 1) = zero;
                                if (dimension == 3)
                                {
                                    mat_Fj.col(i + 2).setZero();
                                    incr(i + 2) = zero;
                                }
                                if (!ind_ok)
                                {
                                    l2g.at(pos + i)     = 0xDEADBEEF;
                                    l2l.at(pos + i)     = false;
                                    l2g.at(pos + i + 1) = face_offset + ind_sol++;
                                    if (dimension == 3)
                                    {
                                        l2g.at(pos + i + 2) = face_offset + ind_sol++;
                                    }
                                }
                            }
                            ind_ok = true;
                            break;
                        }
                        case disk::DY:
                        {
                            for (size_t i = 0; i < num_face_dofs; i += dimension)
                            {
                                mat_Fj.col(i).setZero();
                                incr(i) = zero;
                                if (dimension == 3)
                                {
                                    mat_Fj.col(i + 2).setZero();
                                    incr(i + 2) = zero;
                                }
                                if (!ind_ok)
                                {
                                    l2g.at(pos + i)     = face_offset + ind_sol++;
                                    l2g.at(pos + i + 1) = 0xDEADBEEF;
                                    l2l.at(pos + i + 1) = false;
                                    if (dimension == 3)
                                    {
                                        l2g.at(pos + i + 2) = face_offset + ind_sol++;
                                    }
                                }
                            }
                            ind_ok = true;
                            break;
                        }
                        case disk::DZ:
                        {
                            if (dimension != 3)
                                throw std::invalid_argument("You are not in 3D");
                            for (size_t i = 0; i < num_face_dofs; i += dimension)
                            {
                                mat_Fj.col(i).setZero();
                                incr(i) = zero;
                                mat_Fj.col(i + 1).setZero();
                                incr(i + 1) = zero;
                                if (!ind_ok)
                                {
                                    l2g.at(pos + i)     = face_offset + ind_sol++;
                                    l2g.at(pos + i + 1) = face_offset + ind_sol++;
                                    l2g.at(pos + i + 2) = 0xDEADBEEF;
                                    l2l.at(pos + i + 2) = false;
                                }
                            }
                            ind_ok = true;
                            break;
                        }
                        case disk::DXDY:
                        {
                            for (size_t i = 0; i < num_face_dofs; i += dimension)
                            {
                                if (dimension == 3)
                                {
                                    mat_Fj.col(i + 2).setZero();
                                    incr(i + 2) = zero;
                                }
                                if (!ind_ok)
                                {
                                    l2g.at(pos + i)     = 0xDEADBEEF;
                                    l2g.at(pos + i + 1) = 0xDEADBEEF;
                                    l2l.at(pos + i)     = false;
                                    l2l.at(pos + i + 1) = false;
                                    if (dimension == 3)
                                    {
                                        l2g.at(pos + i + 2) = face_offset + ind_sol++;
                                    }
                                }
                            }
                            ind_ok = true;
                            break;
                        }
                        case disk::DXDZ:
                        {
                            if (dimension != 3)
                                throw std::invalid_argument("You are not in 3D");
                            for (size_t i = 0; i < num_face_dofs; i += dimension)
                            {
                                mat_Fj.col(i + 1).setZero();
                                incr(i + 1) = zero;
                                if (!ind_ok)
                                {
                                    l2g.at(pos + i)     = 0xDEADBEEF;
                                    l2l.at(pos + i)     = false;
                                    l2g.at(pos + i + 1) = face_offset + ind_sol++;
                                    l2g.at(pos + i + 2) = 0xDEADBEEF;
                                    l2l.at(pos + i + 2) = false;
                                }
                            }
                            ind_ok = true;
                            break;
                        }
                        case disk::DYDZ:
                        {
                            if (dimension != 3)
                                throw std::invalid_argument("You are not in 3D");
                            for (size_t i = 0; i < num_face_dofs; i += dimension)
                            {
                                mat_Fj.col(i).setZero();
                                incr(i) = zero;
                                if (!ind_ok)
                                {
                                    l2g.at(pos + i)     = face_offset + ind_sol++;
                                    l2g.at(pos + i + 1) = 0xDEADBEEF;
                                    l2g.at(pos + i + 2) = 0xDEADBEEF;
                                    l2l.at(pos + i + 1) = false;
                                    l2l.at(pos + i + 2) = false;
                                }
                            }
                            ind_ok = true;
                            break;
                        }
                        default:
                        {
                            throw std::logic_error("Unknown Dirichlet Conditions");
                            break;
                        }
                    }

                    rhs_bc.segment(face_j * num_face_dofs, num_face_dofs) += mat_Fj * incr;
                }
            }
        }
        assert(lc.first.rows() == lc.first.cols());
        assert(lc.first.rows() == lc.second.size());
        assert(lc.second.size() == l2g.size());
        assert(lc.second.size() == rhs_bc.size());

#ifdef FILL_COLMAJOR
        for (size_t j = 0; j < lc.first.cols(); j++)
        {
            if (l2l[j])
            {
                for (size_t i = 0; i < lc.first.rows(); i++)
                {
                    if (l2l[i])
                    {
                        m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));
                    }
                }
                RHS(l2g.at(j)) += lc.second(j) - rhs_bc(j);
            }
        }
#else
        for (size_t i = 0; i < lc.first.rows(); i++)
        {
            if (l2l[i])
            {
                for (size_t j = 0; j < lc.first.cols(); j++)
                {
                    if (l2l[j])
                    {
                        m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));
                    }
                }
                RHS(l2g.at(i)) += lc.second(i) - rhs_bc(i);
            }
        }
#endif
    }

    vector_type
    expand_solution_nl(const mesh_type&                msh,
                       const bnd_type&                 bnd,
                       const vector_type&              solution,
                       const std::vector<vector_type>& sol_F,
                       int                             di = 0)
    {
        assert(solution.size() == m_num_unknowns);
        assert(sol_F.size() == msh.faces_size());
        const auto face_degree   = m_hdi.face_degree();
        const auto num_face_dofs = vector_basis_size(face_degree, dimension - 1, dimension);

        vector_type ret = vector_type::Zero(num_face_dofs * msh.faces_size());

        for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
        {
            const auto bfc             = *itor;
            const auto face_id         = msh.lookup(bfc);
            const auto face_offset     = face_id * num_face_dofs;
            const auto compress_offset = face_compress_map.at(face_id);

            if (bnd.is_dirichlet_face(face_id))
            {
                size_t sol_ind = 0;

                const vector_type proj_bcf =
                  project_function(msh, bfc, face_degree, bnd.dirichlet_boundary_func(face_id), di);
                vector_type incr = proj_bcf - sol_F[face_id];
                assert(proj_bcf.size() == num_face_dofs);

                switch (bnd.dirichlet_boundary_type(face_id))
                {
                    case disk::DIRICHLET:
                    {
                        ret.segment(face_offset, num_face_dofs) = incr;
                        break;
                    }
                    case disk::CLAMPED:
                    {
                        incr                                    = -sol_F[face_id];
                        ret.segment(face_offset, num_face_dofs) = incr;
                        break;
                    }
                    case disk::DX:
                    {

                        for (size_t i = 0; i < num_face_dofs; i += dimension)
                        {
                            ret(face_offset + i)     = incr(i);
                            ret(face_offset + i + 1) = solution(compress_offset + sol_ind++);
                            if (dimension == 3)
                            {
                                ret(face_offset + i + 2) = solution(compress_offset + sol_ind++);
                            }
                        }
                        break;
                    }
                    case disk::DY:
                    {
                        for (size_t i = 0; i < num_face_dofs; i += dimension)
                        {
                            ret(face_offset + i)     = solution(compress_offset + sol_ind++);
                            ret(face_offset + i + 1) = incr(i + 1);
                            if (dimension == 3)
                            {
                                ret(face_offset + i + 2) = solution(compress_offset + sol_ind++);
                            }
                        }
                        break;
                    }
                    case disk::DZ:
                    {
                        if (dimension != 3)
                            throw std::invalid_argument("You are not in 3D");
                        for (size_t i = 0; i < num_face_dofs; i += dimension)
                        {
                            ret(face_offset + i)     = solution(compress_offset + sol_ind++);
                            ret(face_offset + i + 1) = solution(compress_offset + sol_ind++);
                            ret(face_offset + i + 2) = incr(i + 2);
                        }
                        break;
                    }
                    case disk::DXDY:
                    {
                        for (size_t i = 0; i < num_face_dofs; i += dimension)
                        {
                            ret(face_offset + i)     = incr(i);
                            ret(face_offset + i + 1) = incr(i + 1);
                            if (dimension == 3)
                            {
                                ret(face_offset + i + 2) = solution(compress_offset + sol_ind++);
                            }
                        }
                        break;
                    }
                    case disk::DXDZ:
                    {
                        if (dimension != 3)
                            throw std::invalid_argument("You are not in 3D");
                        for (size_t i = 0; i < num_face_dofs; i += dimension)
                        {
                            ret(face_offset + i)     = incr(i);
                            ret(face_offset + i + 1) = solution(compress_offset + sol_ind++);
                            ret(face_offset + i + 2) = incr(i + 2);
                        }
                        break;
                    }
                    case disk::DYDZ:
                    {
                        if (dimension != 3)
                            throw std::invalid_argument("You are not in 3D");
                        for (size_t i = 0; i < num_face_dofs; i += dimension)
                        {
                            ret(face_offset + i)     = solution(compress_offset + sol_ind++);
                            ret(face_offset + i + 1) = incr(i + 1);
                            ret(face_offset + i + 2) = incr(i + 2);
                        }
                        break;
                    }
                    default:
                    {
                        throw std::logic_error("Unknown Dirichlet Conditions");
                        break;
                    }
                }
            }
            else
            {
                ret.segment(face_offset, num_face_dofs) = solution.segment(compress_offset, num_face_dofs);
            }
        }

        return ret;
    }

    void
    impose_neumann_boundary_conditions(const mesh_type& msh, const bnd_type& bnd)
    {
        const auto face_degree   = m_hdi.face_degree();
        const auto num_face_dofs = vector_basis_size(face_degree, dimension - 1, dimension);

        if (bnd.nb_faces_neumann() > 0)
        {
            for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
            {
                const auto bfc     = *itor;
                const auto face_id = msh.lookup(bfc);

                if (bnd.is_neumann_face(face_id))
                {
                    const size_t      face_offset = face_compress_map.at(face_id);
                    auto              fb          = make_vector_monomial_basis(msh, bfc, face_degree);
                    const vector_type neumann = make_rhs(msh, bfc, fb, bnd.neumann_boundary_func(face_id), face_degree);

                    assert(neumann.size() == num_face_dofs);

                    if (bnd.is_dirichlet_face(face_id))
                    {
                        switch (bnd.dirichlet_boundary_type(face_id))
                        {
                            case disk::DIRICHLET:
                            {
                                throw std::invalid_argument("You tried to impose both Dirichlet and "
                                                            "Neumann conditions on the same face");
                                break;
                            }
                            case disk::CLAMPED:
                            {
                                throw std::invalid_argument("You tried to impose both Dirichlet and "
                                                            "Neumann conditions on the same face");
                                break;
                            }
                            case disk::DX:
                            {
                                for (size_t i = 0; i < num_face_dofs; i += dimension)
                                {
                                    RHS(face_offset + i + 1) += neumann(i + 1);
                                    if (dimension == 3)
                                    {
                                        RHS(face_offset + i + 2) += neumann(i + 2);
                                    }
                                }
                                break;
                            }
                            case disk::DY:
                            {
                                for (size_t i = 0; i < num_face_dofs; i += dimension)
                                {
                                    RHS(face_offset + i) = neumann(i);
                                    if (dimension == 3)
                                    {
                                        RHS(face_offset + i + 2) += neumann(i + 2);
                                    }
                                }

                                break;
                            }
                            case disk::DZ:
                            {
                                if (dimension != 3)
                                    throw std::invalid_argument("You are not in 3D");
                                for (size_t i = 0; i < num_face_dofs; i += dimension)
                                {
                                    RHS(face_offset + i) += neumann(i);
                                    RHS(face_offset + i + 1) += neumann(i + 1);
                                }
                                break;
                            }
                            case disk::DXDY:
                            {
                                for (size_t i = 0; i < num_face_dofs; i += dimension)
                                {
                                    if (dimension == 3)
                                    {
                                        RHS(face_offset + i + 2) += neumann(i + 2);
                                    }
                                }
                                break;
                            }
                            case disk::DXDZ:
                            {
                                if (dimension != 3)
                                    throw std::invalid_argument("You are not in 3D");
                                for (size_t i = 0; i < num_face_dofs; i += dimension)
                                {
                                    RHS(face_offset + i + 1) += neumann(i + 1);
                                }
                                break;
                            }
                            case disk::DYDZ:
                            {
                                if (dimension != 3)
                                    throw std::invalid_argument("You are not in 3D");
                                for (size_t i = 0; i < num_face_dofs; i += dimension)
                                {
                                    RHS(face_offset + i) += neumann(i);
                                }
                                break;
                            }
                            default:
                            {
                                throw std::logic_error("Unknown Dirichlet Conditions");
                                break;
                            }
                        }
                    }
                    else
                    {
                        RHS.segment(face_offset, num_face_dofs) += neumann;
                    }
                }
            }
        }
    }

    void
    finalize()
    {
        LHS.setFromTriplets(m_triplets.begin(), m_triplets.end());
        m_triplets.clear();
    }
};

template<typename Mesh, typename BoundaryType>
auto
make_mechanics_assembler(const Mesh& msh, const hho_degree_info hdi, const BoundaryType& bnd)
{
    return assembler_mechanics<Mesh>(msh, hdi, bnd);
}

/**
 * @brief Assembler for HHO methods where the discrete problem is scalar and formulated only in terms of primal unknowns \f$(u_T, u_{\partial T}) \f$
 *
 * @tparam Mesh
 */
template<typename Mesh>
class scalar_primal_hho_assembler
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::face::id_type   face_id_type;

    typedef disk::scalar_boundary_conditions<Mesh> boundary_type;

    std::vector<ident_raw_t> compress_table;

    hho_degree_info                   di;
    std::vector<Triplet<scalar_type>> triplets;

    size_t num_all_faces, num_dirichlet_faces, num_other_faces;
    size_t fbs;
    size_t system_size;

    class assembly_index
    {
        size_t idx;
        bool   assem;

      public:
        assembly_index(size_t i, bool as) : idx(i), assem(as) {}

        operator size_t() const
        {
            if (!assem)
                throw std::logic_error("Invalid assembly_index");

            return idx;
        }

        bool
        assemble() const
        {
            return assem;
        }

        friend std::ostream&
        operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };

  public:
    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;


    SparseMatrix<scalar_type> LHS;
    vector_type               RHS;


    scalar_primal_hho_assembler(const Mesh& msh, const hho_degree_info& hdi, const boundary_type& bnd) : di(hdi)
    {
        num_all_faces       = msh.faces_size();
        num_dirichlet_faces = bnd.nb_faces_dirichlet();
        num_other_faces     = num_all_faces - num_dirichlet_faces;

        compress_table.resize(num_all_faces);

        size_t compressed_offset = 0;
        for (size_t i = 0; i < num_all_faces; i++)
        {
            compress_table.at(i) = compressed_offset;
            if (!bnd.is_dirichlet_face(i))
            {
                compressed_offset++;
            }
        }
        assert(compressed_offset == num_other_faces);

        fbs = scalar_basis_size(di.face_degree(), Mesh::dimension - 1);

        system_size = fbs * num_other_faces;

        this->initialize();

        // preallocate memory
        triplets.reserve(2 * (di.face_degree() + 1) * system_size);
    }

    void
    initialize()
    {
        initialize_lhs();
        initialize_rhs();
    }

    void
    initialize_lhs()
    {
        LHS = SparseMatrix<scalar_type>(system_size, system_size);
        return;
    }

    void
    initialize_rhs()
    {
        RHS = vector_type::Zero(system_size);
        return;
    }

    void
    assemble(const Mesh&          msh,
             const cell_type&     cl,
             const boundary_type& bnd,
             const matrix_type&   lhs,
             const vector_type&   rhs)
    {
        const auto                  fcs_id = faces_id(msh, cl);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(fcs_id.size() * fbs);

        vector_type dirichlet_data = vector_type::Zero(fcs_id.size() * fbs);

        for (size_t face_i = 0; face_i < fcs_id.size(); face_i++)
        {
            const auto face_id         = fcs_id[face_i];
            const auto face_LHS_offset = compress_table.at(face_id) * fbs;

            const bool dirichlet = bnd.is_dirichlet_face(face_id);

            for (size_t i = 0; i < fbs; i++)
            {
                asm_map.push_back(assembly_index(face_LHS_offset + i, !dirichlet));
            }

            if (dirichlet)
            {
                auto dirichlet_fun = bnd.dirichlet_boundary_func(face_id);

                const auto fc = *std::next(msh.faces_begin(), face_id);

                dirichlet_data.segment(face_i * fbs, fbs) =
                  project_function(msh, fc, di.face_degree(), dirichlet_fun, di.face_degree());
            }
        }

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if (asm_map[j].assemble())
                    triplets.push_back(Triplet<scalar_type>(asm_map[i], asm_map[j], lhs(i, j)));
                else
                    RHS(asm_map[i]) -= lhs(i, j) * dirichlet_data(j);
            }

            RHS(asm_map[i]) += rhs(i);
        }
    }

    void
    finalize(void)
    {
        LHS.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();
    }

    size_t
    global_system_size() const
    {
        return system_size;
    }

    vector_type
    take_local_solution(const Mesh& msh, const cell_type& cl, const boundary_type& bnd, const vector_type& sol) const
    {
        const auto fcs_id    = faces_id(msh, cl);
        const auto num_faces = fcs_id.size();

        vector_type ret = vector_type::Zero(num_faces * fbs);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            const auto face_id = fcs_id[face_i];

            if (bnd.is_dirichlet_face(face_id))
            {
                const auto dirichlet_bf = bnd.dirichlet_boundary_func(face_id);
                const auto fc           = *std::next(msh.faces_begin(), face_id);
                ret.segment(face_i * fbs, fbs) =
                  project_function(msh, fc, di.face_degree(), dirichlet_bf, di.face_degree());
            }
            else
            {
                const auto face_SOL_offset     = compress_table.at(face_id) * fbs;
                ret.segment(face_i * fbs, fbs) = sol.segment(face_SOL_offset, fbs);
            }
        }
        return ret;
    }

    void
    impose_neumann_boundary_conditions(const Mesh& msh, const boundary_type& bnd)
    {
        if (bnd.nb_faces_neumann() > 0)
        {
            for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
            {
                const auto bfc     = *itor;
                const auto face_id = msh.lookup(bfc);

                if (bnd.is_neumann_face(face_id))
                {
                    if (bnd.is_dirichlet_face(face_id))
                    {
                        throw std::invalid_argument("You tried to impose"
                                                    "both Dirichlet and Neumann conditions on the same face");
                    }
                    else if (bnd.is_robin_face(face_id))
                    {
                        throw std::invalid_argument("You tried to impose"
                                                    "both Robin and Neumann conditions on the same face");
                    }
                    else
                    {
                        const size_t                face_degree = di.face_degree();
                        std::vector<assembly_index> asm_map;
                        asm_map.reserve(fbs);

                        auto face_LHS_offset = compress_table.at(face_id) * fbs;

                        for (size_t i = 0; i < fbs; i++)
                        {
                            asm_map.push_back(assembly_index(face_LHS_offset + i, true));
                        }

                        const auto        fb = make_scalar_monomial_basis(msh, bfc, face_degree);
                        const vector_type neumann =
                          make_rhs(msh, bfc, fb, bnd.neumann_boundary_func(face_id), face_degree);

                        assert(neumann.size() == fbs);
                        for (size_t i = 0; i < neumann.size(); i++)
                        {
                            RHS(asm_map[i]) += neumann(i);
                        }
                    }
                }
            }
        }
    }

    void
    impose_robin_boundary_conditions(const Mesh& msh, const boundary_type& bnd)
    {
        if (bnd.nb_faces_robin() > 0)
        {
            for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
            {
                const auto bfc     = *itor;
                const auto face_id = msh.lookup(bfc);

                if (bnd.is_robin_face(face_id))
                {
                    if (bnd.is_neumann_face(face_id))
                    {
                        switch (bnd.neumann_boundary_type(face_id))
                        {
                            case disk::NEUMANN:
                                throw std::invalid_argument(
                                  "You tried to impose both Neumann and Robin conditions on the same face");
                                break;
                            default: throw std::logic_error("Unknown Neumann Conditions"); break;
                        }
                    }
                    else if (bnd.is_dirichlet_face(face_id))
                    {
                        switch (bnd.dirichlet_boundary_type(face_id))
                        {
                            case disk::DIRICHLET:
                                throw std::invalid_argument(
                                  "You tried to impose both Dirichlet and Robin conditions on the same face");
                                break;
                            default: throw std::logic_error("Unknown Dirichlet Conditions"); break;
                        }
                    }
                    else
                    {
                        const size_t                face_degree = di.face_degree();
                        std::vector<assembly_index> asm_map;
                        asm_map.reserve(fbs);

                        const auto face_LHS_offset = compress_table.at(face_id) * fbs;

                        for (size_t i = 0; i < fbs; i++)
                            asm_map.push_back(assembly_index(face_LHS_offset + i, true));

                        const auto        fb    = make_scalar_monomial_basis(msh, bfc, face_degree);
                        const vector_type robin = make_rhs(msh, bfc, fb, bnd.robin_boundary_func(face_id), face_degree);
                        assert(robin.size() == fbs);

                        const matrix_type mass = make_mass_matrix(msh, bfc, fb);

                        for (size_t i = 0; i < fbs; i++)
                        {
                            RHS(asm_map[i]) += robin(i);

                            for (size_t j = 0; j < fbs; j++)
                            {
                                triplets.push_back(Triplet<scalar_type>(asm_map[i], asm_map[j], mass(i, j)));
                            }
                        }
                    }
                }
            }
        }
    }
};

template<typename Mesh>
auto
make_scalar_primal_hho_assembler(const Mesh&                             msh,
                                 const hho_degree_info&                  hdi,
                                 const scalar_boundary_conditions<Mesh>& bnd)
{
    return scalar_primal_hho_assembler<Mesh>(msh, hdi, bnd);
}

} // end disk
