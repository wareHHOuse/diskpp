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
#include <regex>
#include <sstream>
#include <iomanip>

#include <unistd.h>

#include "loaders/loader.hpp"
#include "methods/hho"
#include "solvers/solver.hpp"
#include "output/silo.hpp"
#include "timecounter.h"
#include "colormanip.h"

using namespace disk;
using namespace Eigen;


namespace obs_priv {
/* This way of estimating mesh h was used long time ago in disk++ but later
 * it was removed because it is incorrect. I put it here only to be able to get
 * the the same numbers I had in the first revision of the paper (and correct
 * it), but otherwise this function must not be used.
 *
 * See geometry_all.hpp at commit e278738.
 */
template<typename Mesh>
typename Mesh::coordinate_type
mesh_h(const Mesh& msh)
{
    typename Mesh::coordinate_type h{};
    for (auto itor = msh.cells_begin(); itor != msh.cells_end(); itor++)
    {
        auto cell = *itor;
        auto cell_measure = measure(msh, cell);

        auto fcs = faces(msh, cell);
        typename Mesh::coordinate_type face_sum{};
        for (auto& f : fcs)
        {
            auto m = measure(msh, f);
            face_sum += m;
        }
        h = std::max(h, cell_measure/face_sum);
    }

    return h;
}
} //namespace obs_priv;


/* Right-hand side definition */
template<typename Mesh>
struct rhs_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T,2,Storage> >
{
    typedef Mesh<T,2,Storage>                       mesh_type;
    typedef typename mesh_type::coordinate_type     scalar_type;
    typedef typename mesh_type::point_type          point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto r0 = 0.7;
        auto r = sqrt( pt.x()*pt.x() + pt.y()*pt.y() );

        if ( r > r0 )
            return -16 *r*r + 8*r0*r0;
        else
            return -8.0*( r0*r0*(r0*r0 + 1) ) + 8*r0*r0*r*r;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T,3,Storage> >
{
    typedef Mesh<T,3,Storage>                       mesh_type;
    typedef typename mesh_type::coordinate_type     scalar_type;
    typedef typename mesh_type::point_type          point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto r0 = 0.7;
        auto r = sqrt( pt.x()*pt.x() + pt.y()*pt.y() + pt.z()*pt.z() );

        if ( r > r0 )
            return -4 * (2*r*r + 3*(r*r - r0*r0));
        else
            return -8.0*r0*r0*(1 - r*r + r0*r0);
    }
};

template<typename Mesh>
auto make_rhs_function(const Mesh& msh)
{
    return rhs_functor<Mesh>();
}

/* Solution definition */
template<typename Mesh>
struct sol_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct sol_functor< Mesh<T,2,Storage> >
{
    typedef Mesh<T,2,Storage>                       mesh_type;
    typedef typename mesh_type::coordinate_type     scalar_type;
    typedef typename mesh_type::point_type          point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto r0 = 0.7;
        auto r = sqrt(pt.x()*pt.x() + pt.y()*pt.y());
        auto s = r*r - r0*r0;
        auto t = std::max(s, 0.0);
        return t*t;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct sol_functor< Mesh<T,3,Storage> >
{
    typedef Mesh<T,3,Storage>                       mesh_type;
    typedef typename mesh_type::coordinate_type     scalar_type;
    typedef typename mesh_type::point_type          point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto r0 = 0.7;
        auto r = sqrt( pt.x()*pt.x() + pt.y()*pt.y() + pt.z()*pt.z() );
        auto s = r*r - r0*r0;
        auto t = std::max(s, 0.0);
        return t*t;
    }
};

template<typename Mesh>
auto make_solution_function(const Mesh& msh)
{
    return sol_functor<Mesh>();
}

/* Boundary conditions definition */
template<typename Mesh>
auto make_boundary_function(const Mesh& msh)
{
    return sol_functor<Mesh>();
}

/* Obstacle definition */
template<typename Mesh>
struct obstacle_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct obstacle_functor< Mesh<T,2,Storage> >
{
    typedef Mesh<T,2,Storage>                       mesh_type;
    typedef typename mesh_type::coordinate_type     scalar_type;
    typedef typename mesh_type::point_type          point_type;

    scalar_type operator()(const point_type& pt) const
    {
        return 0.0;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct obstacle_functor< Mesh<T,3,Storage> >
{
    typedef Mesh<T,3,Storage>                       mesh_type;
    typedef typename mesh_type::coordinate_type     scalar_type;
    typedef typename mesh_type::point_type          point_type;

    scalar_type operator()(const point_type& pt) const
    {
        return 0.0;
    }
};

template<typename Mesh>
auto make_obstacle_function(const Mesh& msh)
{
    return obstacle_functor<Mesh>();
}

/* Assembler for the obstacle problem (see "Bubbles enriched quadratic finite
 * element method for the 3D-elliptic obstacle problem - S. Gaddam, T. Gudi",
 * eqn. 5.1 onwards) */
template<typename Mesh>
class obstacle_assembler
{
    using T = typename Mesh::coordinate_type;
    std::vector<size_t>                 face_ct, A_ct, B_ct; //compress tables
    std::vector<size_t>                 face_et, A_et, B_et; //expand tables
    std::vector<bool>                   is_in_set_A;

    std::vector< Triplet<T> >           triplets;

    hho_degree_info                     hdi;

    size_t      num_all_faces, num_dirichlet_faces, num_other_faces;
    size_t      num_all_cells, num_A_cells, num_I_cells;

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

        bool assemble() const { return assem; }

        friend std::ostream& operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };


public:

    SparseMatrix<T>         LHS;
    Matrix<T, Dynamic, 1>   RHS;

    obstacle_assembler(const Mesh& msh,
                       const std::vector<bool>& in_A,
                       const hho_degree_info& p_hdi)
        : hdi(p_hdi), is_in_set_A(in_A)
    {
        /* Boundary faces are dirichlet by default */
        auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
            return msh.is_boundary(fc);
        };

        num_all_faces = msh.faces_size();
        num_dirichlet_faces = std::count_if(msh.faces_begin(), msh.faces_end(), is_dirichlet);
        num_other_faces = num_all_faces - num_dirichlet_faces;

        assert( is_in_set_A.size() == msh.cells_size() );

        num_all_cells = msh.cells_size();
        num_A_cells = std::count(is_in_set_A.begin(), is_in_set_A.end(), true);
        num_I_cells = std::count(is_in_set_A.begin(), is_in_set_A.end(), false);

        /* Make A tables: keep the unknowns of cells in set I */
        A_ct.resize( num_all_cells );
        A_et.resize( num_I_cells );
        for (size_t i = 0, co = 0; i < num_all_cells; i++)
        {
            auto cl = *std::next(msh.cells_begin(), i);
            if ( !is_in_set_A.at(i) )
            {
                A_ct.at(i) = co;
                A_et.at(co) = i;
                co++;
            }
        }

        /* Make face tables */
        face_ct.resize( num_all_faces );
        face_et.resize( num_other_faces );
        for (size_t i = 0, co = 0; i < num_all_faces; i++)
        {
            auto fc = *std::next(msh.faces_begin(), i);
            if ( !is_dirichlet(fc) )
            {
                face_ct.at(i) = co;
                face_et.at(co) = i;
                co++;
            }
        }

        /* Make B tables: keep the unknowns of cells in set A */
        B_ct.resize( num_all_cells );
        B_et.resize( num_A_cells );
        for (size_t i = 0, co = 0; i < num_all_cells; i++)
        {
            auto cl = *std::next(msh.cells_begin(), i);
            if ( is_in_set_A.at(i) )
            {
                B_ct.at(i) = co;
                B_et.at(co) = i;
                co++;
            }
        }

        auto celdeg = hdi.cell_degree();
        auto facdeg = hdi.face_degree();

        auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
        auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);

        auto system_size = cbs * (num_I_cells + num_A_cells) + fbs * num_other_faces;

        assert( system_size == cbs * msh.cells_size() + fbs * num_other_faces );

        LHS = SparseMatrix<T>( system_size, system_size );
        RHS = Matrix<T, Dynamic, 1>::Zero( system_size );
    }

    void dump_tables() const
    {
        std::cout << "Compress table for A: " << std::endl;
        for (size_t i = 0; i < A_ct.size(); i++)
            std::cout << i << " -> " << A_ct.at(i) << std::endl;

        std::cout << "Compress table for faces: " << std::endl;
        for (size_t i = 0; i < face_ct.size(); i++)
            std::cout << i << " -> " << face_ct.at(i) << std::endl;

        std::cout << "Compress table for B: " << std::endl;
        for (size_t i = 0; i < B_ct.size(); i++)
            std::cout << i << " -> " << B_ct.at(i) << std::endl;
    }

    template<typename Function>
    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs, const Matrix<T, Dynamic, 1>& rhs,
             const Matrix<T, Dynamic, 1>& gamma, const Function& dirichlet_bf)
    {
        auto celdeg = hdi.cell_degree();
        auto facdeg = hdi.face_degree();

        auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
        auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);

        std::vector<assembly_index> asm_map_row, asm_map_col;
        //asm_map.reserve(cbs + 4*fbs);

        /* Cell dofs local to global */
        auto cell_offset        = offset(msh, cl);
        auto cell_LHS_offset    = A_ct.at(cell_offset) * cbs;
        bool cell_needs_asm_A   = !is_in_set_A.at(cell_offset);
        bool cell_needs_asm_B   = is_in_set_A.at(cell_offset);

        for (size_t i = 0; i < cbs; i++)
        {
            asm_map_row.push_back( assembly_index(cell_offset+i, true) );
            asm_map_col.push_back( assembly_index(cell_LHS_offset+i, cell_needs_asm_A) );
        }

        /* Face dofs local to global */
        auto fcs = faces(msh, cl);
        auto numfcs = fcs.size();
        Matrix<T, Dynamic, 1> dirichlet_data = Matrix<T, Dynamic, 1>::Zero(cbs + numfcs*fbs);
        for (size_t face_i = 0; face_i < numfcs; face_i++)
        {
            auto fc = fcs[face_i];
            auto face_offset = offset(msh, fc);
            auto face_LHS_offset_rows = cbs * num_all_cells + face_ct.at(face_offset)*fbs;
            auto face_LHS_offset_cols = cbs * num_I_cells + face_ct.at(face_offset)*fbs;

            bool dirichlet = msh.is_boundary(fc);

            for (size_t i = 0; i < fbs; i++)
            {
                asm_map_row.push_back( assembly_index(face_LHS_offset_rows+i, !dirichlet) );
                asm_map_col.push_back( assembly_index(face_LHS_offset_cols+i, !dirichlet) );
            }

            if (dirichlet)
            {
                auto fb = make_scalar_monomial_basis(msh, fc, hdi.face_degree());
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb);
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_bf);
                dirichlet_data.block(cbs+face_i*fbs, 0, fbs, 1) = mass.llt().solve(rhs);
            }
        }

        //assert( asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if ( !asm_map_row[i].assemble() )
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if ( asm_map_col[j].assemble() )
                    triplets.push_back( Triplet<T>(asm_map_row[i], asm_map_col[j], lhs(i,j)) );
                else
                {
                    if (j < cbs)
                        RHS(asm_map_row[i]) -= lhs(i,j)*gamma(cell_offset);
                    else
                        RHS(asm_map_row[i]) -= lhs(i,j)*dirichlet_data(j);
                }
            }
        }


        /* Needed both in case A and I */
        RHS.block(cell_offset, 0, cbs, 1) += rhs.block(0, 0, cbs, 1);

        if (cell_needs_asm_B)
        {
            auto offset_row = cell_offset * cbs;
            auto offset_col = num_I_cells * cbs + num_other_faces * fbs + B_ct.at(cell_offset);
            triplets.push_back( Triplet<T>(offset_row, offset_col, 1.0) );
        }

    } // assemble_A()

    template<typename Function>
    void
    expand_solution(const Mesh& msh,
    const Matrix<T, Dynamic, 1>& solution, const Function& dirichlet_bf,
    const Matrix<T, Dynamic, 1>& gamma,
    Matrix<T, Dynamic, 1>& alpha, Matrix<T, Dynamic, 1>& beta)
    {
        auto celdeg = hdi.cell_degree();
        auto facdeg = hdi.face_degree();

        auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
        auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);

        alpha.block(0, 0, num_all_cells*cbs, 1) = gamma;
        for (size_t i = 0; i < num_I_cells; i++)
        {
            auto exp_ofs = A_et.at(i);
            alpha.block(exp_ofs*cbs, 0, cbs, 1) = solution.block(i*cbs, 0, cbs, 1);
        }

        beta = Matrix<T, Dynamic, 1>::Zero(msh.cells_size());
        for (size_t i = 0; i < num_A_cells; i++)
        {
            auto exp_ofs = B_et.at(i);
            beta.block(exp_ofs*cbs, 0, cbs, 1) = solution.block(num_I_cells*cbs + num_other_faces*fbs + i*cbs, 0, cbs, 1);
        }

        for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
        {
            auto fc = *itor;
            auto face_ofs = offset(msh, fc);

            bool dirichlet = msh.is_boundary(fc);

            if (dirichlet)
            {
                auto fb = make_scalar_monomial_basis(msh, fc, hdi.face_degree());
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb);
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_bf);
                alpha.block(num_all_cells*cbs + face_ofs*fbs, 0, fbs, 1) = mass.llt().solve(rhs);
            }
            else
            {
                auto face_offset = offset(msh, fc);
                auto face_SOL_offset = cbs * num_I_cells + face_ct.at(face_offset)*fbs;
                alpha.block(num_all_cells*cbs + face_ofs*fbs, 0, fbs, 1) = solution.block(face_SOL_offset, 0, fbs, 1);
            }
        }
    }

    void finalize(void)
    {
        LHS.setFromTriplets( triplets.begin(), triplets.end() );
        triplets.clear();
    }

};

/* Assembler for the obstacle problem (see "Bubbles enriched quadratic finite
 * element method for the 3D-elliptic obstacle problem - S. Gaddam, T. Gudi",
 * eqn. 5.1 onwards) */
template<typename Mesh>
class obstacle_assembler_nitsche
{
    using T = typename Mesh::coordinate_type;
    std::vector<size_t>                 A_ct, B_ct; //compress tables
    std::vector<size_t>                 A_et, B_et; //expand tables
    std::vector<bool>                   is_in_set_A;

    std::vector< Triplet<T> >           triplets;

    hho_degree_info                     hdi;

    size_t      num_faces;
    size_t      num_all_cells, num_A_cells, num_I_cells;

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

        bool assemble() const { return assem; }

        friend std::ostream& operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };


public:

    SparseMatrix<T>         LHS;
    Matrix<T, Dynamic, 1>   RHS;

    obstacle_assembler_nitsche(const Mesh& msh,
                               const std::vector<bool>& in_A,
                               const hho_degree_info& p_hdi)
        : hdi(p_hdi), is_in_set_A(in_A)
    {
        num_faces = msh.faces_size();

        assert( is_in_set_A.size() == msh.cells_size() );

        num_all_cells = msh.cells_size();
        num_A_cells = std::count(is_in_set_A.begin(), is_in_set_A.end(), true);
        num_I_cells = std::count(is_in_set_A.begin(), is_in_set_A.end(), false);

        /* Make A tables: keep the unknowns of cells in set I */
        A_ct.resize( num_all_cells );
        A_et.resize( num_I_cells );
        for (size_t i = 0, co = 0; i < num_all_cells; i++)
        {
            auto cl = *std::next(msh.cells_begin(), i);
            if ( !is_in_set_A.at(i) )
            {
                A_ct.at(i) = co;
                A_et.at(co) = i;
                co++;
            }
        }

        /* Make B tables: keep the unknowns of cells in set A */
        B_ct.resize( num_all_cells );
        B_et.resize( num_A_cells );
        for (size_t i = 0, co = 0; i < num_all_cells; i++)
        {
            auto cl = *std::next(msh.cells_begin(), i);
            if ( is_in_set_A.at(i) )
            {
                B_ct.at(i) = co;
                B_et.at(co) = i;
                co++;
            }
        }

        auto celdeg = hdi.cell_degree();
        auto facdeg = hdi.face_degree();

        auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
        auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);

        auto system_size = cbs * (num_I_cells + num_A_cells) + fbs * num_faces;

        assert( system_size == cbs * msh.cells_size() + fbs * num_faces );

        LHS = SparseMatrix<T>( system_size, system_size );
        RHS = Matrix<T, Dynamic, 1>::Zero( system_size );
    }

    void dump_tables() const
    {
        std::cout << "Compress table for A: " << std::endl;
        for (size_t i = 0; i < A_ct.size(); i++)
            std::cout << i << " -> " << A_ct.at(i) << std::endl;

        std::cout << "Compress table for B: " << std::endl;
        for (size_t i = 0; i < B_ct.size(); i++)
            std::cout << i << " -> " << B_ct.at(i) << std::endl;
    }

    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs, const Matrix<T, Dynamic, 1>& rhs,
             const Matrix<T, Dynamic, 1>& gamma)
    {
        auto celdeg = hdi.cell_degree();
        auto facdeg = hdi.face_degree();

        auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
        auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);

        std::vector<assembly_index> asm_map_row, asm_map_col;
        //asm_map.reserve(cbs + 4*fbs);

        /* Cell dofs local to global */
        auto cell_offset        = offset(msh, cl);
        auto cell_LHS_offset    = A_ct.at(cell_offset) * cbs;
        bool cell_needs_asm_A   = !is_in_set_A.at(cell_offset);
        bool cell_needs_asm_B   = is_in_set_A.at(cell_offset);

        for (size_t i = 0; i < cbs; i++)
        {
            asm_map_row.push_back( assembly_index(cell_offset+i, true) );
            asm_map_col.push_back( assembly_index(cell_LHS_offset+i, cell_needs_asm_A) );
        }

        /* Face dofs local to global */
        auto fcs = faces(msh, cl);
        auto numfcs = fcs.size();
        Matrix<T, Dynamic, 1> dirichlet_data = Matrix<T, Dynamic, 1>::Zero(cbs + numfcs*fbs);
        for (size_t face_i = 0; face_i < numfcs; face_i++)
        {
            auto fc = fcs[face_i];
            auto face_offset = offset(msh, fc);
            auto face_LHS_offset_rows = cbs * num_all_cells + face_offset*fbs;
            auto face_LHS_offset_cols = cbs * num_I_cells + face_offset*fbs;

            for (size_t i = 0; i < fbs; i++)
            {
                asm_map_row.push_back( assembly_index(face_LHS_offset_rows+i, true) );
                asm_map_col.push_back( assembly_index(face_LHS_offset_cols+i, true) );
            }
        }

        //assert( asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if ( !asm_map_row[i].assemble() )
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if ( asm_map_col[j].assemble() )
                    triplets.push_back( Triplet<T>(asm_map_row[i], asm_map_col[j], lhs(i,j)) );
                else
                    if (j < cbs)
                        RHS(asm_map_row[i]) -= lhs(i,j)*gamma(cell_offset);
            }

            RHS(asm_map_row[i]) += rhs(i);
        }

        if (cell_needs_asm_B)
        {
            auto offset_row = cell_offset * cbs;
            auto offset_col = num_I_cells * cbs + num_faces * fbs + B_ct.at(cell_offset);
            triplets.push_back( Triplet<T>(offset_row, offset_col, 1.0) );
        }

    } // assemble_A()

    void
    expand_solution(const Mesh& msh,
    const Matrix<T, Dynamic, 1>& solution, const Matrix<T, Dynamic, 1>& gamma,
    Matrix<T, Dynamic, 1>& alpha, Matrix<T, Dynamic, 1>& beta)
    {
        auto celdeg = hdi.cell_degree();
        auto facdeg = hdi.face_degree();

        auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
        auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);

        alpha.block(0, 0, num_all_cells*cbs, 1) = gamma;
        for (size_t i = 0; i < num_I_cells; i++)
        {
            auto exp_ofs = A_et.at(i);
            alpha.block(exp_ofs*cbs, 0, cbs, 1) = solution.block(i*cbs, 0, cbs, 1);
        }

        beta = Matrix<T, Dynamic, 1>::Zero(msh.cells_size());
        for (size_t i = 0; i < num_A_cells; i++)
        {
            auto exp_ofs = B_et.at(i);
            beta.block(exp_ofs*cbs, 0, cbs, 1) = solution.block(num_I_cells*cbs + num_faces*fbs + i*cbs, 0, cbs, 1);
        }

        for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
        {
            auto fc = *itor;
            auto face_ofs = offset(msh, fc);
            auto face_SOL_offset = cbs * num_I_cells + face_ofs*fbs;
            alpha.block(num_all_cells*cbs + face_ofs*fbs, 0, fbs, 1) = solution.block(face_SOL_offset, 0, fbs, 1);
        }
    }

    void finalize(void)
    {
        LHS.setFromTriplets( triplets.begin(), triplets.end() );
        triplets.clear();
    }

};

template<typename T, typename Mesh>
Matrix<T, Dynamic, 1>
take_local_data(const Mesh& msh, const typename Mesh::cell_type& cl,
                const hho_degree_info& hdi,
                const Matrix<T, Dynamic, 1>& expanded_solution)
{
    auto celdeg = hdi.cell_degree();
    auto facdeg = hdi.face_degree();

    auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);

    auto cell_offset        = offset(msh, cl);
    auto cell_SOL_offset    = cell_offset * cbs;

    auto fcs = faces(msh, cl);
    auto numfcs = fcs.size();

    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs + numfcs*fbs);
    ret.block(0, 0, cbs, 1) = expanded_solution.block(cell_SOL_offset, 0, cbs, 1);

    for (size_t face_i = 0; face_i < numfcs; face_i++)
    {
        auto fc = fcs[face_i];
        auto face_offset = offset(msh, fc);
        auto face_SOL_offset = cbs * msh.cells_size() + face_offset*fbs;
        ret.block(cbs+face_i*fbs, 0, fbs, 1) = expanded_solution.block(face_SOL_offset, 0, fbs, 1);
    }

    return ret;
}

template<typename Mesh>
auto make_obstacle_assembler_strong(const Mesh& msh, const std::vector<bool>& in_A, const hho_degree_info& hdi)
{
    return obstacle_assembler<Mesh>(msh, in_A, hdi);
}

template<typename Mesh>
auto make_obstacle_assembler_nitsche(const Mesh& msh, const std::vector<bool>& in_A, const hho_degree_info& hdi)
{
    return obstacle_assembler_nitsche<Mesh>(msh, in_A, hdi);
}


template<typename Mesh>
void
export_to_silo(const Mesh&,
               const Matrix<typename Mesh::coordinate_type, Dynamic, 1>&,
               int cycle = -1)
{}

/*
template<typename T>
export_to_silo(const tetrahedral_mesh<T>& msh,
               const Matrix<T, Dynamic, 1>& data)
{
    disk::silo_database silo;
    silo.create("obstacle.silo");
    silo.add_mesh(msh, "mesh");

    silo.add_variable("mesh", "alpha", data, disk::zonal_variable_t );
    silo.close();
}
*/

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
void
export_to_silo(const Mesh<T, 2, Storage>& msh,
               const Matrix<T, Dynamic, 1>& data, int cycle = -1)
{
    disk::silo_database silo;

    if (cycle == -1)
        silo.create("obstacle.silo");
    else
    {
        std::stringstream ss;
        ss << "out_" << cycle << ".silo";
        silo.create(ss.str());
    }

    silo.add_mesh(msh, "mesh");

    silo.add_variable("mesh", "alpha", data, disk::zonal_variable_t );
    silo.close();
}

/**
 * Apply Dirichlet boundary conditions via Nitsche trick, version for diffusion.
 * rec_op: the reconstruction operator
 * lhs: the element lhs contribution
 * rhs: the element rhs contribution
 *
 * lhs and rhs get updated in order to apply the boundary conditions.
 */

template<typename T>
using DM = Matrix<T, Dynamic, Dynamic>;

template<typename T>
using DV = Matrix<T, Dynamic, 1>;

template<typename Mesh, typename Function>
void
apply_dirichlet_via_nitsche(const Mesh& msh,
                            const typename Mesh::cell_type& cl,
                            DM<typename Mesh::coordinate_type>& rec_op,
                            DM<typename Mesh::coordinate_type>& lhs,
                            DV<typename Mesh::coordinate_type>& rhs,
                            const hho_degree_info& hdi,
                            const Function& bcs_fun,
                            typename Mesh::coordinate_type penalization = 1.0)
{
    const size_t DIM = Mesh::dimension;
    typedef typename Mesh::coordinate_type              T;
    typedef Matrix<T, Dynamic, Dynamic>       matrix_type;
    typedef Matrix<T, Dynamic, 1>             vector_type;
    typedef Matrix<T, Dynamic, DIM>           grad_type;

    auto rb = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
    auto cbs = scalar_basis_size(hdi.cell_degree(), DIM);
    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    for (size_t face_i = 0; face_i < fcs.size(); face_i++)
    {
        auto fc = fcs[face_i];

        if ( !msh.is_boundary(fc) )
            continue;

        auto fb = make_scalar_monomial_basis(msh, fc, hdi.face_degree());

        auto rbs = rb.size();
        auto fbs = fb.size();

        matrix_type mass = matrix_type::Zero(fbs, fbs);
        matrix_type trace = matrix_type::Zero(fbs, rbs-1);
        vector_type vec1 = vector_type::Zero(rbs-1);
        auto hf = measure(msh, fc);
        auto blkofs = cbs+fbs*face_i;
        auto qps = integrate(msh, fc, 2*hdi.reconstruction_degree());
        for (auto& qp : qps)
        {
            auto n = normal(msh, cl, fc);
            grad_type   r_dphi = rb.eval_gradients( qp.point() );
            vector_type r_dphi_n = r_dphi.block(1,0,rbs-1,DIM)*n;
            vector_type f_phi = fb.eval_functions( qp.point() );

            trace += qp.weight() * f_phi * r_dphi_n.transpose();
            mass += qp.weight() * f_phi * f_phi.transpose();
            vec1 += qp.weight() * r_dphi_n * bcs_fun(qp.point());

            rhs.block(blkofs, 0, fbs, 1) += qp.weight() * f_phi * (penalization/hf) * bcs_fun(qp.point());
        }

        matrix_type l = trace * rec_op;

        assert(l.cols() == lhs.cols());
        lhs.block(blkofs, 0, l.rows(), l.cols()) -= l;
        lhs.block(0, blkofs, l.cols(), l.rows()) -= l.transpose();
        lhs.block(blkofs, blkofs, mass.rows(), mass.cols()) += (penalization/hf)*mass;

        rhs -= rec_op.transpose() * vec1;
    }
}

enum class bc_mode
{
    BC_STRONG,
    BC_NITSCHE
};

template<typename T>
struct obstacle_solver_config
{
    obstacle_solver_config()
        : degree(0), maxiter(100), tolerance(1e-7), penalization(1.0),
          symmetry(-1.0), mode(bc_mode::BC_STRONG)
    {}

    size_t          degree;         /* Method degree */
    size_t          maxiter;        /* Max number of iterations */
    T               tolerance;      /* Iteration tolerance */
    T               penalization;   /* Nitsche penalization value*/
    T               symmetry;       /* Nitsche symmetry term sign */
    bc_mode         mode;           /* Boundary condition mode */
};

template<typename T>
struct run_info
{
    T   mesh_h;
    T   error;
};

template<typename Mesh>
using CT = typename Mesh::coordinate_type;

template<typename Mesh>
bool
obstacle_solver_strong(const Mesh& msh,
                       const obstacle_solver_config< CT<Mesh> >& config,
                       run_info< CT<Mesh> >& ri)
{
    typedef typename Mesh::coordinate_type              scalar_type;
    typedef Matrix<scalar_type, Dynamic, Dynamic>       matrix_type;
    typedef Matrix<scalar_type, Dynamic, 1>             vector_type;

    if ( (config.degree != 0) and (config.degree != 1) )
    {
        std::cout << "Method degree must be 0 or 1." << std::endl;
        return false;
    }

    hho_degree_info hdi(0, config.degree, 1);

    auto num_cells = msh.cells_size();
    auto num_faces = msh.faces_size();
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto num_alpha_dofs = num_cells + fbs*num_faces;

    vector_type     alpha = vector_type::Zero( num_alpha_dofs );
    vector_type     beta  = vector_type::Ones( num_cells );
    vector_type     gamma = vector_type::Zero( num_cells );

    scalar_type     c = 1.0;

    auto obstacle_fun = make_obstacle_function(msh);

    size_t i = 0;
    for (auto& cl : msh)
    {
        auto bar = barycenter(msh, cl);
        gamma(i++) = obstacle_fun(bar);
    }

    auto rhs_fun = make_rhs_function(msh);
    auto bcs_fun = make_boundary_function(msh);

    std::cout << green << "Running obstacle solver, k = ";
    std::cout << config.degree << ", " << "boundary conditions enforced strongly";
    std::cout << nocolor << std::endl;
    for (size_t iter = 0; iter < config.maxiter; iter++)
    {
        std::cout << yellow << "  Iteration " << iter+1 << nocolor << std::endl;

        /* Compute the beta quantity (see "Bubbles enriched quadratic finite element method for
        *  the 3D-elliptic obstacle problem - S. Gaddam, T. Gudi", eqn. 5.1 onwards) */
        vector_type diff = beta + c * ( alpha.head(num_cells) - gamma );
        std::vector<bool> in_A;
        in_A.resize(diff.size());

        for (size_t i = 0; i < diff.size(); i++)
            in_A.at(i) = (diff(i) < 0);

        timecounter tc;

        auto assembler  = make_obstacle_assembler_strong(msh, in_A, hdi);
        
        tc.tic();
        for (auto& cl : msh)
        {
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
            auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);
            auto rhs    = make_rhs(msh, cl, cb, rhs_fun, 1);

            matrix_type lhs = gr.second + stab;

            assembler.assemble(msh, cl, lhs, rhs, gamma, bcs_fun);
        }

        assembler.finalize();

        tc.toc();

        std::cout << "    Assembly time: " << tc << ", ";
        std::cout << assembler.LHS.rows() << " DoFs" << std::endl;

        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        vector_type sol = vector_type::Zero(systsz);

        tc.tic();
        disk::solvers::pardiso_params<scalar_type> pparams;
        pparams.out_of_core = PARDISO_OUT_OF_CORE_IF_NEEDED;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);
        tc.toc();
        std::cout << "    Solution time: " << tc << std::endl;

        auto sol_fun = make_solution_function(msh);

        vector_type alpha_prev = alpha;
        assembler.expand_solution(msh, sol, bcs_fun, gamma, alpha, beta);

        auto eps = (alpha_prev - alpha).norm();
        std::cout << "    Current eps: " << eps << std::endl;

        if ( eps < config.tolerance )
            break;
    }

    scalar_type error = 0.0;
    for (auto& cl : msh)
    {
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
        auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);
        matrix_type lc = gr.second + stab;

        auto sol_fun = make_solution_function(msh);
        vector_type local = take_local_data(msh, cl, hdi, alpha);
        vector_type proj  = project_function(msh, cl, hdi, sol_fun, 1);

        vector_type diff = local - proj;
        error += diff.dot(lc*diff);
    }

    vector_type data = alpha.head(msh.cells_size());
    export_to_silo(msh, data);

    auto correct_h = average_diameter(msh);

    std::cout << "wrong h = " << obs_priv::mesh_h(msh) << " ";
    std::cout << "correct h = " << correct_h << " ";
    std::cout << "error(energy) = " << std::sqrt(error) << std::endl;

    ri.mesh_h = correct_h;
    ri.error = std::sqrt(error);

    return true;
}

template<typename Mesh>
bool
obstacle_solver_nitsche(const Mesh& msh,
                        const obstacle_solver_config< CT<Mesh> >& config,
                        run_info< CT<Mesh> >& ri)
{
    typedef typename Mesh::coordinate_type              scalar_type;
    typedef Matrix<scalar_type, Dynamic, Dynamic>       matrix_type;
    typedef Matrix<scalar_type, Dynamic, 1>             vector_type;

    if ( (config.degree != 0) and (config.degree != 1) )
    {
        std::cout << "Method degree must be 0 or 1." << std::endl;
        return false;
    }

    hho_degree_info hdi(0, config.degree, 1);


    auto num_cells = msh.cells_size();
    auto num_faces = msh.faces_size();
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto num_alpha_dofs = num_cells + fbs*num_faces;

    vector_type     alpha = vector_type::Zero( num_alpha_dofs );
    vector_type     beta  = vector_type::Ones( num_cells );
    vector_type     gamma = vector_type::Zero( num_cells );

    scalar_type c = 1.0;

    auto obstacle_fun = make_obstacle_function(msh);

    size_t i = 0;
    for (auto& cl : msh)
    {
        auto bar = barycenter(msh, cl);
        gamma(i++) = obstacle_fun(bar);
    }

    auto rhs_fun = make_rhs_function(msh);
    auto bcs_fun = make_boundary_function(msh);

    std::cout << green << "Running obstacle solver, k = ";
    std::cout << config.degree << ", " << "boundary conditions enforced via Nitsche";
    std::cout << nocolor << std::endl;

    for (size_t iter = 0; iter < config.maxiter; iter++)
    {
        std::cout << yellow << "  Iteration " << iter+1 << nocolor << std::endl;

        /* Compute the beta quantity (see "Bubbles enriched quadratic finite element method for
        *  the 3D-elliptic obstacle problem - S. Gaddam, T. Gudi", eqn. 5.1 onwards) */
        vector_type diff = beta + c * ( alpha.head(num_cells) - gamma );
        std::vector<bool> in_A;
        in_A.resize(diff.size());

        for (size_t i = 0; i < diff.size(); i++)
            in_A.at(i) = (diff(i) < 0);

        timecounter tc;

        auto assembler  = make_obstacle_assembler_nitsche(msh, in_A, hdi);
        
        tc.tic();
        for (auto& cl : msh)
        {
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
            auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);

            matrix_type lhs = gr.second + stab;

            vector_type rhs = vector_type::Zero(lhs.cols());
            rhs.head(cbs) = make_rhs(msh, cl, cb, rhs_fun, 1);

            apply_dirichlet_via_nitsche(msh, cl, gr.first, lhs, rhs, hdi, bcs_fun, config.penalization);

            assembler.assemble(msh, cl, lhs, rhs, gamma);
        }

        assembler.finalize();

        tc.toc();

        std::cout << "    Assembly time: " << tc << ", ";
        std::cout << assembler.LHS.rows() << " DoFs" << std::endl;

        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        vector_type sol = vector_type::Zero(systsz);

        tc.tic();
        disk::solvers::pardiso_params<scalar_type> pparams;
        pparams.out_of_core = PARDISO_OUT_OF_CORE_IF_NEEDED;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);
        tc.toc();
        std::cout << "    Solution time: " << tc << std::endl;

        auto sol_fun = make_solution_function(msh);

        vector_type alpha_prev = alpha;
        assembler.expand_solution(msh, sol, gamma, alpha, beta);

        auto eps = (alpha_prev - alpha).norm();
        std::cout << "    Current eps: " << eps << std::endl;

        if ( eps < config.tolerance )
            break;
    }

    scalar_type error = 0.0;
    for (auto& cl : msh)
    {
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
        auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);
        matrix_type lc = gr.second + stab;

        vector_type rhs = vector_type::Zero(lc.cols());
        //apply_dirichlet_via_nitsche(msh, cl, gr.first, lc, rhs, hdi, bcs_fun, config.penalization);

        auto sol_fun = make_solution_function(msh);
        vector_type local = take_local_data(msh, cl, hdi, alpha);
        vector_type proj  = project_function(msh, cl, hdi, sol_fun, 1);

        vector_type diff = local - proj;

        error += diff.dot(lc*diff);
    }

    vector_type data = alpha.head(msh.cells_size());
    export_to_silo(msh, data);

    auto correct_h = average_diameter(msh);

    std::cout << "wrong h = " << obs_priv::mesh_h(msh) << " ";
    std::cout << "correct h = " << correct_h << " ";
    std::cout << "error(energy) = " << std::sqrt(error) << std::endl;

    ri.mesh_h = correct_h;
    ri.error = std::sqrt(error);

    return true;
}



template<typename Mesh>
bool
obstacle_solver(const Mesh& msh,
                const obstacle_solver_config< CT<Mesh> >& config,
                run_info< CT<Mesh> >& ri)
{
    switch (config.mode)
    {
        case bc_mode::BC_STRONG:
            return obstacle_solver_strong(msh, config, ri);

        case bc_mode::BC_NITSCHE:
            return obstacle_solver_nitsche(msh, config, ri);
    }

    return false;
}


template<typename Mesh>
bool
hho_solver_xx(const Mesh& msh, size_t degree, size_t maxiter, typename Mesh::coordinate_type penalization)
{
    typedef typename Mesh::coordinate_type  scalar_type;
    typedef typename Mesh::point_type       point_type;

    hho_degree_info hdi(degree, degree, degree+1);


    auto num_cells = msh.cells_size();
    auto num_faces = msh.faces_size();


    SparseMatrix<scalar_type>           LHS;
    Matrix<scalar_type, Dynamic, 1>     RHS;
    std::vector<Triplet<scalar_type>>   triplets;

    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto system_size = cbs*num_cells + fbs*num_faces;

    LHS = SparseMatrix<scalar_type>(system_size, system_size);
    RHS = Matrix<scalar_type, Dynamic, 1>::Zero(system_size);

    timecounter tc;

    auto rhs_fun = [](const point_type& pt) -> auto { return 0.0; };
    auto bcs_fun = [](const point_type& pt) -> auto { return std::exp(pt.x()); };

    scalar_type dt = 0.01;

    tc.tic();
    size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto fcs    = faces(msh, cl);
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
        auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);
        Matrix<scalar_type, Dynamic, Dynamic> lhs = gr.second + stab;

        Matrix<scalar_type, Dynamic, 1> rhs = Matrix<scalar_type, Dynamic, 1>::Zero(lhs.cols());
        rhs.head(cbs) = make_rhs(msh, cl, cb, rhs_fun, 1);

        apply_dirichlet_via_nitsche(msh, cl, gr.first, lhs, rhs, hdi, bcs_fun, penalization);

        Matrix<scalar_type, Dynamic, Dynamic> cell_mass = make_mass_matrix(msh, cl, cb);
        lhs = lhs*dt;
        lhs.block(0,0,cbs,cbs) += cell_mass;

        /* Make local-to-global mapping */
        std::vector<size_t> l2g(cbs + fcs.size() * fbs);
        for (size_t i = 0; i < cbs; i++)
            l2g[i] = cell_i*cbs + i;

        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto f_ofs = offset(msh, fcs[i]);
            for (size_t j = 0; j < fbs; j++)
                l2g[cbs + i*fbs+j] = cbs*num_cells + f_ofs*fbs+j;
        }

        /* Assemble cell contributions */
        for (size_t i = 0; i < lhs.rows(); i++)
        {
            for (size_t j = 0; j < lhs.cols(); j++)
                triplets.push_back( Triplet<scalar_type>(l2g[i], l2g[j], lhs(i,j)) );

            RHS(l2g[i]) += rhs(i);
        }

        cell_i++;
    }

    LHS.setFromTriplets(triplets.begin(), triplets.end());

    tc.toc();
    std::cout << " Assembly time: " << tc << std::endl;

    Matrix<scalar_type, Dynamic, 1> u_prev = Matrix<scalar_type, Dynamic, 1>::Zero(system_size);
    Matrix<scalar_type, Dynamic, 1> u;

    scalar_type t = 0.0;

    for (size_t i = 0; t < 2.0; i++, t += dt)
    {
        std::cout << "Step " << i << std::endl;
        disk::solvers::pardiso_params<scalar_type> pparams;
        Matrix<scalar_type, Dynamic, 1> upp = u_prev;
        upp.tail(msh.faces_size() * fbs) = Matrix<scalar_type, Dynamic, 1>::Zero(msh.faces_size() * fbs);
        upp = upp + dt*RHS;
        mkl_pardiso(pparams, LHS, upp, u);

        Matrix<scalar_type, Dynamic, 1> sol_silo = Matrix<scalar_type, Dynamic, 1>::Zero(msh.cells_size());

        for (size_t i = 0; i < msh.cells_size(); i++)
            sol_silo(i) = u(i*cbs);

        export_to_silo( msh, sol_silo, i );

        u_prev = Matrix<scalar_type, Dynamic, 1>::Zero(system_size);

        cell_i = 0;
        for (auto& cl : msh)
        {
            auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            Matrix<scalar_type, Dynamic, Dynamic> cell_mass = make_mass_matrix(msh, cl, cb);
            u_prev.block(cbs*cell_i, 0, cbs, 1) = cell_mass*u.block(cbs*cell_i, 0, cbs, 1);
            cell_i++;
        }
   } 


    /*
    size_t systsz = LHS.rows();
    size_t nnz = LHS.nonZeros();


    Matrix<scalar_type, Dynamic, 1> sol = Matrix<scalar_type, Dynamic, 1>::Zero(systsz);

    tc.tic();
    disk::solvers::pardiso_params<scalar_type> pparams;
    mkl_pardiso(pparams, LHS, RHS, sol);
    tc.toc();
    std::cout << " Solution time: " << tc << std::endl;

    Matrix<scalar_type, Dynamic, 1> sol_silo = Matrix<scalar_type, Dynamic, 1>::Zero(msh.cells_size());

    for (size_t i = 0; i < msh.cells_size(); i++)
        sol_silo(i) = sol(i*cbs);

    export_to_silo( msh, sol_silo );
    */

    return true;
}


template<typename T>
bool
solver_driver(const char *filename,
              const obstacle_solver_config<T>& config,
              run_info<T>& ri)
{
    if (std::regex_match(filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;

        typedef disk::generic_mesh<T, 2>  mesh_type;

        mesh_type msh;
        disk::fvca5_mesh_loader<T, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        auto tr = [](const typename mesh_type::point_type& pt) -> auto {

            auto px = -1 * ( 1-pt.x() ) + 1 * pt.x();
            auto py = -1 * ( 1-pt.y() ) + 1 * pt.y();
            return typename mesh_type::point_type({px, py});
        };

        msh.transform(tr);

        return obstacle_solver(msh, config, ri);
    }

    if (std::regex_match(filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;

        typedef disk::simplicial_mesh<T, 2>  mesh_type;

        mesh_type msh;
        disk::netgen_mesh_loader<T, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        auto tr = [](const typename mesh_type::point_type& pt) -> auto {

            auto px = -1 * ( 1-pt.x() ) + 1 * pt.x();
            auto py = -1 * ( 1-pt.y() ) + 1 * pt.y();
            return typename mesh_type::point_type({px, py});
        };

        msh.transform(tr);

        return obstacle_solver(msh, config, ri);
    }

    if (std::regex_match(filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: Cartesian 2D" << std::endl;

        typedef disk::cartesian_mesh<T, 2>   mesh_type;

        mesh_type msh;
        disk::cartesian_mesh_loader<T, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        auto tr = [](const typename mesh_type::point_type& pt) -> auto {

            auto px = -1 * ( 1-pt.x() ) + 1 * pt.x();
            auto py = -1 * ( 1-pt.y() ) + 1 * pt.y();
            return typename mesh_type::point_type({px, py});
        };

        msh.transform(tr);

        return obstacle_solver(msh, config, ri);
    }

    if (std::regex_match(filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;

        typedef disk::generic_mesh<T, 3>   mesh_type;

        mesh_type msh;
        disk::fvca6_mesh_loader<T, 3> loader;


        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }

        loader.populate_mesh(msh);

        return obstacle_solver(msh, config, ri);
    }

    if (std::regex_match(filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;

        typedef disk::simplicial_mesh<T, 3>   mesh_type;

        mesh_type msh;
        disk::netgen_mesh_loader<T, 3> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        return obstacle_solver(msh, config, ri);
    }

    if (std::regex_match(filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: Cartesian 3D" << std::endl;

        typedef disk::cartesian_mesh<T, 3>   mesh_type;

        mesh_type msh;
        disk::cartesian_mesh_loader<T, 3> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        return obstacle_solver(msh, config, ri);
    }

    return false;
}

template<typename T>
void
do_autotest(const std::vector<std::string>& meshes, const char *outfile, size_t degree)
{
    std::vector<run_info<T>>    ri_strong;
    std::vector<run_info<T>>    ri_nitsche;

    obstacle_solver_config<T> config;
    config.degree = degree;
    config.penalization = 100;

    run_info<T> ri;

    for (size_t i = 0; i < meshes.size(); i++)
    {
        auto mfn = meshes[i].c_str();

        config.mode = bc_mode::BC_STRONG;
        solver_driver(mfn, config, ri);
        ri_strong.push_back(ri);

        config.mode = bc_mode::BC_NITSCHE;
        solver_driver(mfn, config, ri);
        ri_nitsche.push_back(ri);
    }

    auto compute_rate = [](const run_info<T>& prev, const run_info<T>& cur) -> T {
        return std::log(prev.error/cur.error)/std::log(prev.mesh_h/cur.mesh_h);
    };

    std::ofstream ofs(outfile);

    auto fmtfl = [](std::ostream& os) -> std::ostream& {
        os << std::setw(10) << std::setprecision(3) << std::scientific;
        return os;
    };

    auto fmtrate = [](std::ostream& os) -> std::ostream& {
        os << std::setw(10) << std::setprecision(3) << std::defaultfloat;
        return os;
    };

    for (size_t i = 0; i < meshes.size(); i++)
    {
        ofs << fmtfl << ri_strong[i].mesh_h << "  ";
        ofs << fmtfl << ri_strong[i].error << "  ";

        auto strong_rate = compute_rate(ri_strong[i-1], ri_strong[i]);

        if (i > 0)
            ofs << fmtrate << strong_rate << "  ";
        else
            ofs << std::setw(10) << "   ---    ";

       
        ofs << fmtfl << ri_nitsche[i].error << "  ";

        auto nitsche_rate = compute_rate(ri_nitsche[i-1], ri_nitsche[i]);

        if (i > 0)
            ofs << fmtrate << nitsche_rate;
        else
            ofs << std::setw(10) << "   ---   ";

        ofs << std::endl;
    }

}

template<typename T>
void
launch_autotests(void)
{
    std::vector<std::string> meshfiles;
/*
    meshfiles.push_back("../../../meshes/2D_triangles/netgen/tri01.mesh2d");
    meshfiles.push_back("../../../meshes/2D_triangles/netgen/tri02.mesh2d");
    meshfiles.push_back("../../../meshes/2D_triangles/netgen/tri03.mesh2d");
    meshfiles.push_back("../../../meshes/2D_triangles/netgen/tri04.mesh2d");
    meshfiles.push_back("../../../meshes/2D_triangles/netgen/tri05.mesh2d");

    do_autotest<T>(meshfiles, "triangles_k0.txt", 0);
    do_autotest<T>(meshfiles, "triangles_k1.txt", 1);

    meshfiles.clear();
    meshfiles.push_back("../../../meshes/2D_quads/fvca5/mesh2_2.typ1");
    meshfiles.push_back("../../../meshes/2D_quads/fvca5/mesh2_3.typ1");
    meshfiles.push_back("../../../meshes/2D_quads/fvca5/mesh2_4.typ1");
    meshfiles.push_back("../../../meshes/2D_quads/fvca5/mesh2_5.typ1");
    meshfiles.push_back("../../../meshes/2D_quads/fvca5/mesh2_6.typ1");

    do_autotest<T>(meshfiles, "quads_k0.txt", 0);
    do_autotest<T>(meshfiles, "quads_k1.txt", 1);

    meshfiles.clear();
    meshfiles.push_back("../../../meshes/2D_hex/fvca5/hexagonal_1.typ1");
    meshfiles.push_back("../../../meshes/2D_hex/fvca5/hexagonal_2.typ1");
    meshfiles.push_back("../../../meshes/2D_hex/fvca5/hexagonal_3.typ1");
    meshfiles.push_back("../../../meshes/2D_hex/fvca5/hexagonal_4.typ1");
    meshfiles.push_back("../../../meshes/2D_hex/fvca5/hexagonal_5.typ1");

    do_autotest<T>(meshfiles, "hexagons_k0.txt", 0);
    do_autotest<T>(meshfiles, "hexagons_k1.txt", 1);

    meshfiles.clear();
    meshfiles.push_back("../../../meshes/3D_tetras/fvca6/tet.1.msh");
    meshfiles.push_back("../../../meshes/3D_tetras/fvca6/tet.2.msh");
    meshfiles.push_back("../../../meshes/3D_tetras/fvca6/tet.3.msh");
    meshfiles.push_back("../../../meshes/3D_tetras/fvca6/tet.4.msh");
    meshfiles.push_back("../../../meshes/3D_tetras/fvca6/tet.5.msh");

    do_autotest<T>(meshfiles, "tetras_k0.txt", 0);
    do_autotest<T>(meshfiles, "tetras_k1.txt", 1);

    meshfiles.clear();
    meshfiles.push_back("../../../meshes/3D_hexa/diskpp/testmesh-4-4-4.hex");
    meshfiles.push_back("../../../meshes/3D_hexa/diskpp/testmesh-8-8-8.hex");
    meshfiles.push_back("../../../meshes/3D_hexa/diskpp/testmesh-16-16-16.hex");
    meshfiles.push_back("../../../meshes/3D_hexa/diskpp/testmesh-32-32-32.hex");
    //meshfiles.push_back("../../../meshes/3D_hexa/diskpp/testmesh-64-64-64.hex");

    do_autotest<T>(meshfiles, "hexas_k0.txt", 0);
    do_autotest<T>(meshfiles, "hexas_k1.txt", 1);
*/
    meshfiles.clear();
    meshfiles.push_back("../../../meshes/3D_general/fvca6/dbls_10.msh");
    meshfiles.push_back("../../../meshes/3D_general/fvca6/dbls_20.msh");
    meshfiles.push_back("../../../meshes/3D_general/fvca6/dbls_30.msh");
    meshfiles.push_back("../../../meshes/3D_general/fvca6/dbls_40.msh");

    //do_autotest<T>(meshfiles, "prisms_k0.txt", 0);
    do_autotest<T>(meshfiles, "prisms_k1.txt", 1);
}

int main(int argc, char **argv)
{
    using T = double;

    char    *filename       = nullptr;
    int ch;

    bool autotest_mode = false;

    obstacle_solver_config<T> config;
    run_info<T> ri;

    while ( (ch = getopt(argc, argv, "k:i:p:s:NSa")) != -1 )
    {
        switch(ch)
        {
            case 'k':
                config.degree = atoi(optarg);
                if (config.degree != 0 && config.degree != 1)
                {
                    std::cout << "Degree can be 0 or 1." << std::endl;
                    return 1;
                }
                break;

            case 'i':
                config.maxiter = atoi(optarg);
                break;

            case 'p':
                config.penalization = atof(optarg);
                break;

            case 's':
                config.symmetry = atof(optarg);
                break;

            case 'N':
                config.mode = bc_mode::BC_NITSCHE;
                break;

            case 'S':
                config.mode = bc_mode::BC_STRONG;
                break;

            case 'a':
                autotest_mode = true;
                break;

            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    if (autotest_mode)
    {
        launch_autotests<T>();
        return;
    }

    if (argc != 1)
    {
        std::cout << "Please specify a 2D or 3D mesh" << std::endl;
        return 0;
    }

    filename = argv[0];

    solver_driver(filename, config, ri);


    return 0;
}
