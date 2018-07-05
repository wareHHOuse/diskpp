/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Karol Cascavita (C) 2018                     karol.cascavita@enpc.fr
 * Intissar Addali (C) 2018                     intissar.addali@inria.fr
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
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include <map>

#include "colormanip.h"

#include "geometry/geometry.hpp"
#include "loaders/loader.hpp"
#include "revolution/methods/hho"
#include "solvers/solver.hpp"


/***************************************************************************/
/* RHS definition */
template<typename Mesh>
struct rhs_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        return   - 2.* M_PI *  std::sin(2. * M_PI * pt.x());
    }
};

template<typename Mesh>
auto make_rhs_function(const Mesh& msh)
{
    return rhs_functor<Mesh>();
}

/* Expected solution definition */
template<typename Mesh>
struct dirichlet_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct dirichlet_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        return 0.;
    }
};


template<typename Mesh>
auto make_dirichlet_function(const Mesh& msh)
{
    return dirichlet_functor<Mesh>();
}

/* NEUMANN */
template<typename Mesh>
struct neumann_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct neumann_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        return 0.;
    }
};

template<typename Mesh>
 auto make_neumann_function(const Mesh& msh)
{
    return neumann_functor<Mesh>();
}

using namespace revolution;

template<typename Mesh>
std::vector<size_t>
full_offset(const Mesh& msh, const hho_degree_info& hdi)
{
    size_t  i = 1, dofs = 0;
    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    std::vector<size_t> vec(msh.cells_size());
    for (auto itor = msh.cells_begin(); itor != msh.cells_end() -1; itor++)
    {
        auto cl = *itor;
        dofs += howmany_faces(msh, cl) * fbs + cbs;
        vec.at(i++) = dofs;
    }
    return vec;
}

template<typename Mesh>
bool
fix_point_solver(const Mesh& msh, const hho_degree_info& hdi,
                const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
                const typename Mesh::coordinate_type & gamma_0,
                const std::vector<size_t>& is_contact_vector)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    auto max_iter = 1000;
    auto tol = 1.e-6;
    auto rhs_fun = make_rhs_function(msh);

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                        - fbs*msh.boundary_faces_size() ;
    vector_type       full_sol_old = dynamic_vector<T>::Zero(num_full_dofs);
    vector_type       full_sol     = dynamic_vector<T>::Zero(num_full_dofs);

    auto offset_vector = full_offset(msh, hdi);

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        auto cl_count = 0;
        full_sol_old = full_sol;
        auto assembler = make_diffusion_assembler2(msh, hdi, bnd);
        for (auto& cl : msh)
        {
            auto num_dofs = cbs + howmany_faces(msh,cl) * fbs;
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
            auto rhs_f  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            matrix_type A = gr.second + stab;
            vector_type rhs = vector_type::Zero(num_dofs);
            rhs.block(0, 0, cbs, 1) = rhs_f;

            if (is_contact_vector.at(cl_count))
            {
                auto cell_ofs = offset_vector.at(cl_count);
                vector_type u_full = full_sol_old.block(cell_ofs, 0, num_dofs, 1);

            	A -= make_contact_unnamed(msh, cl, hdi, gr.first, gamma_0, bnd);
                rhs -= make_contact_negative(msh, cl, hdi, gr.first, gamma_0, bnd, u_full);
            }

            auto sc = diffusion_static_condensation_compute_full(msh, cl, hdi, A, rhs);
            assembler.assemble(msh, cl, sc.first, sc.second);
            cl_count++;
        }

        assembler.impose_neumann_boundary_conditions(msh, bnd);
        assembler.finalize();

        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
        std::cout << "Mesh faces: " << msh.faces_size() << std::endl;
        std::cout << "Dofs: " << systsz << std::endl;

        dynamic_vector<T> sol = dynamic_vector<T>::Zero(systsz);

        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = true;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

        T error = 0.0 ;
        std::ofstream ofs("sol.dat");

        if(!ofs.is_open())
        {
            std::cout<< "Error opening file"<<std::endl;
            return false;
        }

        cl_count = 0;
        for (auto& cl : msh)
        {
            auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
            vector_type rhs = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            matrix_type A   = gr.second + stab;

            if (is_contact_vector.at(cl_count))
            {
                auto cell_ofs = offset_vector.at(cl_count);
                vector_type u_full = full_sol_old.block(cell_ofs, 0, num_total_dofs, 1);

            	A -= make_contact_unnamed(msh, cl, hdi, gr.first, gamma_0, bnd );
                vector_type rhs_contact =
                    make_contact_negative(msh, cl, hdi, gr.first, gamma_0, bnd, u_full);

                rhs -=  rhs_contact.block(0, 0, cbs, 1);
            }

            vector_type  u_faces = assembler.take_local_data(msh, cl, sol);
            vector_type  u_full  = diffusion_static_condensation_recover(msh, cl, hdi, A, rhs, u_faces);

            auto cell_ofs = offset_vector.at(cl_count);
            vector_type  u_full_old = full_sol_old.block(cell_ofs, 0, num_total_dofs ,1);
            full_sol.block(cell_ofs, 0, num_total_dofs ,1) = u_full;

            vector_type  diff = u_full - u_full_old;
            error += diff.dot(A * diff);

            auto bar = barycenter(msh, cl);
            ofs << bar.x() << " " << bar.y() <<" "<< u_full(0) << std::endl;
            cl_count++;
        }
        ofs.close();

        std::cout << "iter : "<< iter << "   ; error = "<< std::sqrt(error)<< std::endl;
        if( std::sqrt(error)  < tol)
            return true;
    }

    std::cout << "No convergence of the fix point solver" << std::endl;
    return false;
}

template<typename Mesh>
void
run_signorini(const Mesh& msh, const typename Mesh::coordinate_type& gamma_0)
{
    using T = typename Mesh::coordinate_type;
    const size_t degree = 0 ;
    hho_degree_info hdi(degree);

    auto g_N = make_neumann_function(msh);
    auto u_D = make_dirichlet_function(msh);

    typedef disk::mechanics::BoundaryConditionsScalar<Mesh> boundary_type;
    boundary_type  m_bnd(msh);

    m_bnd.addDirichletBC(disk::mechanics::DIRICHLET,1,u_D); //TOP
    m_bnd.addNeumannBC(disk::mechanics::NEUMANN,3,g_N); //
    m_bnd.addNeumannBC(disk::mechanics::NEUMANN,4,g_N);
    m_bnd.addContactBC(disk::mechanics::SIGNORINI,2);

    //cells with contact faces
    auto num_cells = msh.cells_size();
    std::vector<size_t> is_contact_vector = std::vector<size_t>(num_cells);
    size_t i =0;
    for (auto& cl : msh)
    {
        auto fcs=faces(msh, cl);
	    for (auto& fc:fcs)
	    {
	    	auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
            const auto face_id=eid.second;

		    if (m_bnd.is_contact_face(face_id))
		    {
			    is_contact_vector.at(i) = 1;
			    continue;
		    }
	    }
	    i++;
    }

    fix_point_solver( msh, hdi, m_bnd, gamma_0, is_contact_vector);
}

int main(int argc, char **argv)
{
    char    *filename       = nullptr;
    using T = double;

    int ch;
    T gamma_0 = 0.1;

    while ( (ch = getopt(argc, argv, "g:")) != -1 )
    {
        switch(ch)
        {
            case 'g':
                gamma_0 = atof(optarg);
                if (gamma_0 <= 0)
                {
                    std::cout << "gamma_0 must be >0. Falling back to 0.1" << std::endl;
                    gamma_0 = 0.1;
                }
                break;

            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    filename = argv[0];

    /* Medit 2d*/
    if (std::regex_match(filename, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        typedef disk::generic_mesh<T, 2>  mesh_type;
        mesh_type msh = disk::load_medit_2d_mesh<T>(filename);
        run_signorini(msh, gamma_0);
    }

    return 0;
}
