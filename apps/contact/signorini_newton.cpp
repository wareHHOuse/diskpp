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
        return - 2.* M_PI *  std::sin(2. * M_PI * pt.x());
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
auto
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
auto
newton_solver(const Mesh& msh, const hho_degree_info& hdi,
                const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
                const typename Mesh::coordinate_type & gamma_0,
                const typename Mesh::coordinate_type & theta,
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
    dynamic_vector<T> full_sol_old = dynamic_vector<T>::Zero(num_full_dofs);
    dynamic_vector<T> full_sol = dynamic_vector<T>::Zero(num_full_dofs);

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
            rhs.block(0, 0, cbs, 1) -= rhs_f;

            auto cell_ofs = offset_vector.at(cl_count);
            vector_type u_full = full_sol_old.block(cell_ofs, 0, num_dofs, 1);

            if (is_contact_vector.at(cl_count))
            {
                //matrix An(du, v)
            	A -= theta * make_contact_unnamed(msh, cl, hdi, gr.first, gamma_0, bnd);

                //rhs: int [Pn(u_old,v)]_Pn[v]
                rhs += = make_contact_negative_trace(msh, cl, hdi, gr.first, gamma_0, theta, bnd, u_full);
            }
            //rhs: An(u_old,v)
            rhs += A * u_full;

            //matrix Heaviside term
            if(in_set_A)
               A -=  make_contact_heaveside(msh, cl, hdi, gr.first, gamma_0, theta, bnd);

            auto sc = diffusion_static_condensation_compute_full(msh, cl, hdi, A, rhs);
            assembler.assemble(msh, cl, sc.first, sc.second);
            cl_count++;
        }

        assembler.impose_neumann_boundary_conditions(msh, bnd);
        assembler.finalize();

        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        //std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
        //std::cout << "Mesh faces: " << msh.faces_size() << std::endl;
        //std::cout << "Dofs: " << systsz << std::endl;

        dynamic_vector<T> sol = dynamic_vector<T>::Zero(systsz);

        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = true;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

        T error = 0.0 ;
        std::ofstream ofs("sol.dat");

        if(!ofs.is_open())
            std::cout<< "Error opening file"<<std::endl;

        cl_count = 0;
        for (auto& cl : msh)
        {
            auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
            auto rhs_force = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            matrix_type A   = gr.second + stab;

            vector_type rhs = vector_type::Zero(num_dofs);
            rhs.block(0, 0, cbs, 1) -= rhs_force;

            auto cell_ofs = offset_vector.at(cl_count);
            vector_type u_full = full_sol_old.block(cell_ofs, 0, num_total_dofs, 1);

            if (is_contact_vector.at(cl_count))
            {
                //matrix An(du, v)
            	A -= theta * make_contact_unnamed(msh, cl, hdi, gr.first, gamma_0, bnd);

                //rhs: int [Pn(u_old,v)]_Pn[v]
                rhs += = make_contact_negative_trace(msh, cl, hdi, gr.first, gamma_0, theta, bnd, u_full);
            }
            //rhs: An(u_old,v)
            rhs += A * u_full;

            //matrix Heaviside term ::: Dont use this before rhs += A * u_full;
            if(in_set_A)
               A -=  make_contact_heaveside(msh, cl, hdi, gr.first, gamma_0, theta, bnd);

            vector_type cell_rhs = rhs.block(0, 0, cbs, 1);

            Matrix<T, Dynamic, 1> u_faces =
                                        assembler.take_local_data(msh, cl, sol);
            Matrix<T, Dynamic, 1> u_full =
                diffusion_static_condensation_recover(msh, cl, hdi, A, cell_rhs, u_faces);

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

        std::cout << "  "<< iter << "  "<< std::sqrt(error)<< std::endl;
        if( std::sqrt(error)  < tol)
            return 0;
    }
}

template<typename Mesh>
void
run_signorini(const Mesh& msh,
              const typename Mesh::coordinate_type& gamma_0,
              const typename Mesh::coordinate_type& theta)
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
    m_bnd.addNeumannBC(disk::mechanics::NEUMANN,4,g_N); //
    m_bnd.addContactBC(disk::mechanics::SIGNORINI,2); //BOTTOM

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

    newton_solver( msh, hdi, m_bnd, gamma_0, theta, is_contact_vector);

    return;
}

int main(int argc, char **argv)
{
    char    *filename       = nullptr;
    using T = double;

    int ch;
    size_t degree = 1;
    T gamma_0 = 0.1;
    T theta = 1.;

    while ( (ch = getopt(argc, argv, "k:g:npz")) != -1 )
    {
        switch(ch)
        {
            case 'k':
                degree = atoi(optarg);
                if (degree != 0 && degree != 1)
                {
                    std::cout << "Degree can be 0 or 1. Falling back to 1" << std::endl;
                    degree = 1;
                }
                break;
                case 'n':
                    theta = -1.;
                    break;
                case 'p':
                    theta =  1.;
                    break;
                case 'z':
                    theta = 0.;
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
        run_signorini(msh, gamma_0, theta);
    }

    return 0;
}
