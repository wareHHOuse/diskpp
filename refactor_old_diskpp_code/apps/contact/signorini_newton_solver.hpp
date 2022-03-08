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


#include "methods/hho"
#include "output/silo.hpp"
#include "common.hpp"
#include "solvers/solver.hpp"
#include "contrib/colormanip.h"

template<typename Mesh, typename Function, typename Analytical>
std::pair<typename Mesh::coordinate_type, typename Mesh::coordinate_type>
solve_faces(const Mesh&  msh, const Function& rhs_fun, const Analytical& sol_fun,
    const algorithm_parameters<typename Mesh::coordinate_type>& ap,
    const disk::scalar_boundary_conditions<Mesh>& bnd)
{
    std::cout << "INSIDE FACE-BASED TRACE" << std::endl;

    using T = typename Mesh::coordinate_type;

    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    matrix_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>                 vector_type;

    auto is_contact_vector = make_is_contact_vector(msh, bnd);


    hho_degree_info      hdi(ap.degree, ap.degree);

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                    - fbs*msh.boundary_faces_size() ;

    auto offset_vector = full_offset(msh, hdi);

    disk::dynamic_vector<T>  full_sol = disk::dynamic_vector<T>::Zero(num_full_dofs);


    auto max_iter = 1000;
    auto tol = 1.e-9;

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        auto cl_count = 0;
        auto assembler =  make_contact_face_assembler_new(msh, hdi, bnd);
        assembler.imposed_dirichlet_boundary_conditions(msh, bnd, full_sol);

        for (auto& cl : msh)
        {
            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs  = cbs + howmany_faces(msh, cl) * fbs;
            vector_type  u_full  = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
            auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);

            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            matrix_type Ah  = gr.second + stab;

            vector_type Bnegative  = vector_type::Zero(num_total_dofs);
            matrix_type Anitsche   = matrix_type::Zero(num_total_dofs, num_total_dofs);
            matrix_type Aheaviside = matrix_type::Zero(num_total_dofs, num_total_dofs);

            if (is_contact_vector.at(cl_count) == 1)
            {
                Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                Bnegative  = make_hho_negative_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
                Aheaviside = make_hho_heaviside_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
            }

            matrix_type A =   Ah - Anitsche + Aheaviside;
            vector_type b = -(Ah - Anitsche) * u_full - Bnegative;
            b.block(0, 0, cbs, 1) += Lh;

            auto sc = make_scalar_static_condensation(msh, cl, hdi, A, b);
            assembler.assemble(msh, cl, sc.first, sc.second);
            cl_count++;
        }

        assembler.impose_neumann_boundary_conditions(msh, bnd);
        assembler.finalize();

        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        disk::dynamic_vector<T> dsol = disk::dynamic_vector<T>::Zero(systsz);

        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = false;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, dsol);

        T H1_increment = 0.0 ;
        T L2_increment = 0.0 ;

        disk::dynamic_vector<T> diff_sol = disk::dynamic_vector<T>::Zero(num_full_dofs);

        cl_count = 0;

        for (auto& cl : msh)
        {
            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;
            vector_type  u_full = full_sol.block(cell_ofs, 0, num_total_dofs, 1);


            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
            auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);

            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            matrix_type Ah  = gr.second + stab;

            vector_type Bnegative  = vector_type::Zero(num_total_dofs);
            matrix_type Anitsche   = matrix_type::Zero(num_total_dofs, num_total_dofs);
            matrix_type Aheaviside = matrix_type::Zero(num_total_dofs, num_total_dofs);

            if (is_contact_vector.at(cl_count))
            {
                Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                Bnegative  = make_hho_negative_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
                Aheaviside = make_hho_heaviside_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
            }

            matrix_type A =   Ah - Anitsche + Aheaviside;
            vector_type b = -(Ah - Anitsche) * u_full - Bnegative;
            b.block(0, 0, cbs, 1) += Lh;

            vector_type cell_rhs = b.block(0, 0, cbs, 1);
            vector_type du_faces = assembler.take_local_data_increment(msh, cl, dsol);
            vector_type du_full  =
                make_scalar_static_decondensation(msh, cl, hdi, A, cell_rhs, du_faces);

            diff_sol.block(cell_ofs, 0, num_total_dofs ,1) = du_full;

            H1_increment += du_full.dot(A * du_full);

            matrix_type mass  = make_mass_matrix(msh, cl, cb);//, hdi.cell_degree());

            vector_type u_diff = du_full.block(0, 0, cbs, 1);
            L2_increment += u_diff.dot(mass * u_diff);


            cl_count++;
        }

        full_sol += diff_sol;

        std::cout << "  "<< iter << "  "<< std::scientific<< std::setprecision(3);
        std::cout << std::sqrt(H1_increment)<< "   "<< std::sqrt(L2_increment)<<std::endl;

        if( std::sqrt(H1_increment)  < tol)
        {
            std::ofstream efs("solution_whho_faces_i" + tostr(iter) + ".dat");

            if(!efs.is_open())
                std::cout<< "Error opening file"<<std::endl;


            T H1_error  = 0.0 ;
            T L2_error  = 0.0 ;


            auto cl_count = 0;
            for(auto& cl : msh)
            {
                const auto cell_ofs = offset_vector.at(cl_count);
                const auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;
                vector_type  u_full = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

                auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
                auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);

                matrix_type Ah  = gr.second + stab;
                matrix_type Anitsche   = matrix_type::Zero(num_total_dofs, num_total_dofs);
                matrix_type Aheaviside = matrix_type::Zero(num_total_dofs, num_total_dofs);

                if (is_contact_vector.at(cl_count))
                {
                    Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                    Aheaviside = make_hho_heaviside_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
                }

                matrix_type A =   Ah - Anitsche + Aheaviside;

                vector_type realsol = project_function(msh, cl, hdi, sol_fun, 2);

                vector_type diff = realsol - u_full;
                H1_error += diff.dot(Ah*diff);

                auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
                matrix_type mass  = make_mass_matrix(msh, cl, cb, hdi.cell_degree());
                vector_type u_diff = diff.block(0, 0, cbs, 1);
                L2_error += u_diff.dot(mass * u_diff);

                //plot
                auto bar = barycenter(msh, cl);
                efs << bar.x() << " " << bar.y() <<" "<< u_full(0) << std::endl;

                auto fcs = faces(msh, cl);
                auto face_count = 0;
                for(auto fc: fcs)
                {
                    auto offset_face = cbs + face_count * fbs;
                    auto barf = barycenter(msh, fc);
                    efs << barf.x() << " " << barf.y() <<" "<< u_full(offset_face)<< std::endl;
                    //efs << " "<< du_full(offset_face) << std::endl;
                    face_count++;
                }
                cl_count++;

            }
            efs.close();

            //std::cout << "(H1, L2) : "<< std::sqrt(H1_error) <<" - " <<std::sqrt(L2_error) << std::endl;
            std::cout << "  "<< ap.gamma_0 << "    " << iter << "   "<< std::sqrt(H1_error) << "       "<< std::scientific<< std::setprecision(3);
            std::cout << std::sqrt(H1_increment)<< std::endl;

            return std::make_pair(std::sqrt(H1_error), std::sqrt(L2_error));
        }
    }
    return std::make_pair(0., 0.);
}

template<typename Mesh, typename Function, typename Analytical>
dynamic_vector<typename Mesh::coordinate_type>
solve_faces_hier(const Mesh&  msh, const Function& rhs_fun, const Analytical& sol_fun,
    const algorithm_parameters<typename Mesh::coordinate_type>& ap,
    const disk::scalar_boundary_conditions<Mesh>& bnd,
    const hho_degree_info& hdi)
{
    std::cout << "INSIDE FACE-BASED TRACE" << std::endl;
    std::cout << ap << std::endl;

    using T = typename Mesh::coordinate_type;

    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    matrix_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>                 vector_type;

    auto is_contact_vector = make_is_contact_vector(msh, bnd);

    bool warning = (hdi.cell_degree() - 1 != hdi.face_degree())? true :false;
    if (warning)    std::cout << magenta;

    std::cout << " * cell degree :"<< hdi.cell_degree() << std::endl;
    std::cout << " * face degree :"<< hdi.face_degree() << std::endl;

    if (warning)
        std::cout << reset;

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                    - fbs*msh.boundary_faces_size() ;

    auto offset_vector = full_offset(msh, hdi);

    disk::dynamic_vector<T>  full_sol = disk::dynamic_vector<T>::Zero(num_full_dofs);


    //-------------------------------------------------------------------------
    //
    // Initialization with projection of the anylitcal solution
    //
    //-------------------------------------------------------------------------
    #if 0
    auto cl_count_0 = 0;
    std::ofstream efs("solution_whho_fnew_init.dat");

    if(!efs.is_open())
        std::cout<< "Error opening file"<<std::endl;

    for(auto& cl : msh)
    {
        auto cell_ofs = offset_vector.at(cl_count_0);
        auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;

        vector_type  proj_fun = project_function(msh, cl, hdi, sol_fun);
        full_sol.block(cell_ofs, 0, num_total_dofs, 1)= proj_fun;
        cl_count_0++;

        auto bar = barycenter(msh, cl);
        auto fcs = faces(msh, cl);
        auto face_count = 0;
        vector_type u_full = full_sol.block(cell_ofs, 0, num_total_dofs, 1);
        for(auto fc: fcs)
        {
            auto offset_face = cbs + face_count * fbs;
            auto barf = barycenter(msh, fc);
            efs << barf.x() << " " << barf.y() <<" "<< u_full(offset_face)<<std::endl;
            face_count++;
        }
    }
    efs.close();
    #endif
    //-------------------------------------------------------------------------



    auto max_iter = 1000;
    auto tol = 1.e-10;

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        auto cl_count = 0;
        auto assembler =  make_contact_face_assembler_new(msh, hdi, bnd);
        assembler.imposed_dirichlet_boundary_conditions(msh, bnd, full_sol);

        for (auto& cl : msh)
        {
            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs  = cbs + howmany_faces(msh, cl) * fbs;
            vector_type  u_full  = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
            auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);

            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            matrix_type Ah  = gr.second + stab;

            vector_type Bnegative  = vector_type::Zero(num_total_dofs);
            matrix_type Anitsche   = matrix_type::Zero(num_total_dofs, num_total_dofs);
            matrix_type Aheaviside = matrix_type::Zero(num_total_dofs, num_total_dofs);

            if (is_contact_vector.at(cl_count) == 1)
            {
                Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                Bnegative  = make_hho_negative_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
                Aheaviside = make_hho_heaviside_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
            }

            matrix_type A =   Ah - Anitsche + Aheaviside;
            vector_type b = -(Ah - Anitsche) * u_full - Bnegative;
            b.block(0, 0, cbs, 1) += Lh;

            auto sc = make_scalar_static_condensation(msh, cl, hdi, A, b);
            assembler.assemble(msh, cl, sc.first, sc.second);
            cl_count++;
        }

        //assembler.impose_neumann_boundary_conditions(msh, bnd);
        assembler.finalize();

        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        disk::dynamic_vector<T> dsol = disk::dynamic_vector<T>::Zero(systsz);

        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = false;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, dsol);

        T H1_increment = 0.0 ;
        T L2_increment = 0.0 ;

        disk::dynamic_vector<T> diff_sol = disk::dynamic_vector<T>::Zero(num_full_dofs);

        cl_count = 0;

        for (auto& cl : msh)
        {
            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;
            vector_type  u_full = full_sol.block(cell_ofs, 0, num_total_dofs, 1);


            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
            auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);

            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            matrix_type Ah  = gr.second + stab;

            vector_type Bnegative  = vector_type::Zero(num_total_dofs);
            matrix_type Anitsche   = matrix_type::Zero(num_total_dofs, num_total_dofs);
            matrix_type Aheaviside = matrix_type::Zero(num_total_dofs, num_total_dofs);

            if (is_contact_vector.at(cl_count))
            {
                Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                Bnegative  = make_hho_negative_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
                Aheaviside = make_hho_heaviside_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
            }

            matrix_type A =   Ah - Anitsche + Aheaviside;
            vector_type b = -(Ah - Anitsche) * u_full - Bnegative;
            b.block(0, 0, cbs, 1) += Lh;

            vector_type cell_rhs = b.block(0, 0, cbs, 1);
            vector_type du_faces = assembler.take_local_data_increment(msh, cl, dsol);
            vector_type du_full  =
                make_scalar_static_decondensation(msh, cl, hdi, A, cell_rhs, du_faces);

            diff_sol.block(cell_ofs, 0, num_total_dofs ,1) = du_full;

            H1_increment += du_full.dot(A * du_full);

            matrix_type mass  = make_mass_matrix(msh, cl, cb);//, hdi.cell_degree());

            vector_type u_diff = du_full.block(0, 0, cbs, 1);
            L2_increment += u_diff.dot(mass * u_diff);


            cl_count++;
        }

        full_sol += diff_sol;

        std::cout << "  "<< iter << "  "<< std::scientific<< std::setprecision(3);
        std::cout << std::sqrt(H1_increment)<< "   "<< std::sqrt(L2_increment)<<  "  "; //std::endl;

        if( std::sqrt(H1_increment)  < tol)
        {
            std::ofstream efs("solution_whho_face_i" + tostr(iter) + ".dat");

            if(!efs.is_open())
                std::cout<< "Error opening file"<<std::endl;


            T H1_error  = 0.0 ;
            T L2_error  = 0.0 ;


            auto cl_count = 0;
            for(auto& cl : msh)
            {
                const auto cell_ofs = offset_vector.at(cl_count);
                const auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;
                vector_type  u_full = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

                auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
                auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);

                matrix_type Ah  = gr.second + stab;
                matrix_type Anitsche   = matrix_type::Zero(num_total_dofs, num_total_dofs);
                matrix_type Aheaviside = matrix_type::Zero(num_total_dofs, num_total_dofs);

                if (is_contact_vector.at(cl_count))
                {
                    Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                    Aheaviside = make_hho_heaviside_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
                }

                matrix_type A =   Ah - Anitsche + Aheaviside;

                vector_type realsol = project_function(msh, cl, hdi, sol_fun, 2);

                vector_type diff = realsol - u_full;
                H1_error += diff.dot(Ah*diff);

                auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
                matrix_type mass  = make_mass_matrix(msh, cl, cb, hdi.cell_degree());
                vector_type u_diff = diff.block(0, 0, cbs, 1);
                L2_error += u_diff.dot(mass * u_diff);

                //plot
                auto bar = barycenter(msh, cl);
                efs << bar.x() << " " << bar.y() <<" "<< u_full(0) << std::endl;

                auto fcs = faces(msh, cl);
                auto face_count = 0;
                for(auto fc: fcs)
                {
                    auto offset_face = cbs + face_count * fbs;
                    auto barf = barycenter(msh, fc);
                    efs << barf.x() << " " << barf.y() <<" "<< u_full(offset_face)<< std::endl;
                    face_count++;
                }
                cl_count++;

            }
            efs.close();

            std::cout << "(H1, L2) : "<< std::sqrt(H1_error) <<" - " <<std::sqrt(L2_error) << std::endl;

            return full_sol;
        }
    }
    return full_sol;
}

template<typename Mesh, typename Function, typename Analytical>
std::pair<typename Mesh::coordinate_type, typename Mesh::coordinate_type>
solve_cells_full(const Mesh&  msh, const Function& rhs_fun, const Analytical& sol_fun,
    const algorithm_parameters<typename Mesh::coordinate_type>& ap,
    const disk::scalar_boundary_conditions<Mesh>& bnd)
{
    std::cout << "INSIDE CELL-BASED TRACE" << std::endl;
    std::cout << ap << std::endl;
    using T = typename Mesh::coordinate_type;

    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    matrix_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>                 vector_type;

    auto is_contact_vector = make_is_contact_vector(msh, bnd);

    hho_degree_info      hdi(ap.degree+1, ap.degree); //Not allow (degree, degree)

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);


    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                    - fbs*msh.boundary_faces_size() ;

    auto offset_vector = full_offset(msh, hdi);

    disk::dynamic_vector<T>  full_sol = disk::dynamic_vector<T>::Zero(num_full_dofs);


    auto max_iter = 500;
    auto tol = 1.e-9;

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        auto assembler = make_contact_full_assembler_new(msh, hdi, bnd);

        assembler.imposed_dirichlet_boundary_conditions(msh, bnd, full_sol);

        auto cl_count = 0;

        auto res = 0.;

        for (auto& cl : msh)
        {
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());

            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;

            vector_type  u_full = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            vector_type b  = vector_type::Zero(num_total_dofs);
            matrix_type A  = matrix_type::Zero(num_total_dofs, num_total_dofs);

            if (is_contact_vector.at(cl_count) == 1)
            {
                auto gr   = make_hho_contact_scalar_laplacian(msh, cl, hdi, bnd);
                auto stab = make_scalar_hdg_stabilization(msh, cl, hdi);

                matrix_type Ah  = gr.second + stab;
                vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun);//, hdi.cell_degree());

                matrix_type  Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                vector_type  Bnegative  = make_hho_negative(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
                matrix_type  Aheaviside = make_hho_heaviside(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);


                //Original
                A =   Ah - Anitsche + Aheaviside;
                b = -(Ah - Anitsche) * u_full - Bnegative;
                b.block(0, 0, cbs, 1) += Lh;
                res += b.dot(vector_type::Ones(num_total_dofs));
            }
            else
            {
                auto gr   = make_scalar_hho_laplacian(msh, cl, hdi);
                auto stab = make_scalar_hdg_stabilization(msh, cl, hdi);

                matrix_type Ah = gr.second + stab;
                vector_type Lh = make_rhs(msh, cl, cb, rhs_fun);//, hdi.cell_degree());

                A = Ah;
                b = -Ah * u_full;
                b.block(0, 0, cbs, 1) += Lh;

                res += b.dot(vector_type::Ones(num_total_dofs));
            }

            assembler.assemble(msh, cl, A, b);

            cl_count++;
        }

        assembler.impose_neumann_boundary_conditions(msh, bnd);
        assembler.finalize();

        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        disk::dynamic_vector<T> dsol = disk::dynamic_vector<T>::Zero(systsz);

        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = false;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, dsol);


        T H1_increment  = 0.0 ;
        T L2_increment  = 0.0 ;

        cl_count = 0;
        disk::dynamic_vector<T> diff_sol = disk::dynamic_vector<T>::Zero(num_full_dofs);

        for (auto& cl : msh)
        {
            const auto cell_ofs = offset_vector.at(cl_count);
            const auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;

            vector_type  u_full = full_sol.block(cell_ofs, 0, num_total_dofs, 1);
            matrix_type  Ah  = matrix_type::Zero(num_total_dofs, num_total_dofs);
            vector_type  du_full  = vector_type::Zero(num_total_dofs);

            auto stab = make_scalar_hdg_stabilization(msh, cl, hdi);

            if (is_contact_vector.at(cl_count)==1)
            {
                auto gr  = make_hho_contact_scalar_laplacian(msh, cl, hdi, bnd);
                Ah  = gr.second + stab;
            }
            else
            {
                auto gr  = make_scalar_hho_laplacian(msh, cl, hdi);
                Ah  = gr.second + stab;
            }

            //Erase this function, since it coulb be used to take data from full_sol
            //dsol has zero dirichlet conditions so nothing is computed in this faces/
            du_full = assembler.take_local_data_increment(msh, cl, dsol);
            diff_sol.block(cell_ofs, 0, num_total_dofs ,1) = du_full;

            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            matrix_type mass  = make_mass_matrix(msh, cl, cb);//, hdi.cell_degree());
            vector_type u_diff = du_full.block(0, 0, cbs, 1);

            H1_increment += du_full.dot(Ah * du_full);
            L2_increment += u_diff.dot(mass * u_diff);

            cl_count++;
        }

        full_sol += diff_sol;

        std::cout << "  "<< iter << "  "<< std::scientific<< std::setprecision(3);
        std::cout << std::sqrt(H1_increment)<< "   "<< std::sqrt(L2_increment) << std::endl;

        if( std::sqrt(H1_increment)  < tol)
        {
            std::ofstream efs("solution_whho_cfull_i" + tostr(iter) + ".dat");

            if(!efs.is_open())
                std::cout<< "Error opening file"<<std::endl;

            T H1_error  = 0.0 ;
            T H1_error_ref  = 0.0 ;
            T L2_error  = 0.0 ;

            auto cl_count = 0;
            for(auto& cl : msh)
            {
                const auto cell_ofs = offset_vector.at(cl_count);
                const auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;

                vector_type  u_full   = full_sol.block(cell_ofs, 0, num_total_dofs, 1);
                matrix_type  Ah  = matrix_type::Zero(num_total_dofs, num_total_dofs);

                auto stab = make_scalar_hdg_stabilization(msh, cl, hdi);


                if (is_contact_vector.at(cl_count)==1)
                {
                    auto gr  = make_hho_contact_scalar_laplacian(msh, cl, hdi, bnd);
                    Ah  = gr.second + stab;
                }
                else
                {
                    auto gr  = make_scalar_hho_laplacian(msh, cl, hdi);
                    Ah  = gr.second + stab;
                }

                vector_type realsol = project_function(msh, cl, hdi, sol_fun, 2);

                vector_type diff = realsol - u_full;
                H1_error += diff.dot(Ah*diff);

                #if 0
                //H1_error_ref +=  realsol.dot(Ah * realsol);
                {

                    auto cb = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
                    auto qps = integrate(msh, cl, 4 * hdi.reconstruction_degree());
                    for (auto& qp : qps)
                    {

                        //1.2. Grad reference
                        auto ref_dphi = cb.eval_gradients(qp.point());

                        Eigen::Matrix<T,1,2> ref_grad = Eigen::Matrix<T,1,2>::Zero();
                        for (size_t i = 1; i < cb.size(); i++)
                            ref_grad += realsol(i) * ref_dphi.block(i,0,1,2);

                        //1.3. H1-error
                        Eigen::Matrix<T,1,2> diff =  ref_grad;

                        H1_error_ref += qp.weight() * diff.dot(diff);
                    }
                }
                #endif

                auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
                matrix_type mass  = make_mass_matrix(msh, cl, cb, hdi.cell_degree());

                vector_type u_diff = diff.block(0, 0, cbs, 1);
                L2_error += u_diff.dot(mass * u_diff);

                vector_type  du_full = assembler.take_local_data_increment(msh, cl, dsol);

                //plot
                auto bar = barycenter(msh, cl);
                efs << bar.x() << " " << bar.y() <<" "<< u_full(0) << " ";
                efs<< du_full(0)<< std::endl;

                auto fcs = faces(msh, cl);
                auto face_count = 0;
                for(auto fc: fcs)
                {
                    auto offset_face = cbs + face_count * fbs;
                    auto barf = barycenter(msh, fc);
                    efs << barf.x() << " " << barf.y() <<" "<< u_full(offset_face);
                    efs << " "<< du_full(offset_face) << std::endl;
                    face_count++;
                }
                cl_count++;

            }
            efs.close();

            std::cout << "(H1, L2) : "<< std::sqrt(H1_error) <<" - " <<std::sqrt(L2_error) << std::endl;
            //std::cout << "  " << ap.gamma_0 << "    " << iter << "   "<< std::sqrt(H1_error);
            //std::cout << "  " << std::sqrt(H1_error_ref)<< "       " << std::scientific<< std::setprecision(3);
            //std::cout << std::sqrt(H1_increment)<< std::endl;

            return std::make_pair(std::sqrt(H1_error), std::sqrt(L2_error));
        }

    }
    return std::make_pair(0., 0.);
}

template<typename Mesh, typename Function, typename Analytical>
//dynamic_vector<typename Mesh::coordinate_type> //This is for hierarchical (solve conflicts!!)
dynamic_vector<typename Mesh::coordinate_type>
solve_cells_full_hier(const Mesh&  msh, const Function& rhs_fun, const Analytical& sol_fun,
    const algorithm_parameters<typename Mesh::coordinate_type>& ap,
    const disk::scalar_boundary_conditions<Mesh>& bnd,
    const hho_degree_info& hdi)
{
    std::cout << "INSIDE CELL-BASED TRACE HIERARCHICAL" << std::endl;
    std::cout << ap << std::endl;
    using T = typename Mesh::coordinate_type;

    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    matrix_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>                 vector_type;

    auto is_contact_vector = make_is_contact_vector(msh, bnd);

    bool warning = (hdi.cell_degree() - 1 != hdi.face_degree())? true :false;
    if (warning)    std::cout << red;

    std::cout << " * cell degree :"<< hdi.cell_degree() << std::endl;
    std::cout << " * face degree :"<< hdi.face_degree() << std::endl;

    if (warning)
        std::cout << reset;

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);


    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                    - fbs*msh.boundary_faces_size() ;

    auto offset_vector = full_offset(msh, hdi);

    disk::dynamic_vector<T>  full_sol = disk::dynamic_vector<T>::Zero(num_full_dofs);


    auto max_iter = 100;
    auto tol = 1.e-10;

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        auto assembler = make_contact_full_assembler_new(msh, hdi, bnd);

        assembler.imposed_dirichlet_boundary_conditions(msh, bnd, full_sol);

        auto cl_count = 0;

        auto res = 0.;

        for (auto& cl : msh)
        {
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());

            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;

            vector_type  u_full = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            vector_type b  = vector_type::Zero(num_total_dofs);
            matrix_type A  = matrix_type::Zero(num_total_dofs, num_total_dofs);

            if (is_contact_vector.at(cl_count) == 1)
            {
                auto gr   = make_hho_contact_scalar_laplacian(msh, cl, hdi, bnd);
                auto stab = make_scalar_hdg_stabilization(msh, cl, hdi);

                matrix_type Ah  = gr.second + stab;
                vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun);//, hdi.cell_degree());

                matrix_type  Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                vector_type  Bnegative  = make_hho_negative(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
                matrix_type  Aheaviside = make_hho_heaviside(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);

                //Original
                A =   Ah - Anitsche + Aheaviside;
                b = -(Ah - Anitsche) * u_full - Bnegative;
                b.block(0, 0, cbs, 1) += Lh;
                res += b.dot(vector_type::Ones(num_total_dofs));
            }
            else
            {
                auto gr   = make_scalar_hho_laplacian(msh, cl, hdi);
                auto stab = make_scalar_hdg_stabilization(msh, cl, hdi);

                matrix_type Ah = gr.second + stab;
                vector_type Lh = make_rhs(msh, cl, cb, rhs_fun);//, hdi.cell_degree());

                A = Ah;
                b = -Ah * u_full;
                b.block(0, 0, cbs, 1) += Lh;

                res += b.dot(vector_type::Ones(num_total_dofs));
            }

            assembler.assemble(msh, cl, A, b);

            cl_count++;
        }

        assembler.impose_neumann_boundary_conditions(msh, bnd);
        assembler.finalize();

        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        disk::dynamic_vector<T> dsol = disk::dynamic_vector<T>::Zero(systsz);

        std::cout << "here" << std::endl;
        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = false;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, dsol);

        //disk::solvers::conjugated_gradient_params<T> cgparams;
        //conjugated_gradient(cgparams, assembler.LHS, assembler.RHS, dsol);

        T H1_increment  = 0.0 ;
        T L2_increment  = 0.0 ;

        cl_count = 0;
        disk::dynamic_vector<T> diff_sol = disk::dynamic_vector<T>::Zero(num_full_dofs);

        for (auto& cl : msh)
        {
            const auto cell_ofs = offset_vector.at(cl_count);
            const auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;

            vector_type  u_full = full_sol.block(cell_ofs, 0, num_total_dofs, 1);
            matrix_type  Ah  = matrix_type::Zero(num_total_dofs, num_total_dofs);
            vector_type  du_full  = vector_type::Zero(num_total_dofs);

            auto stab = make_scalar_hdg_stabilization(msh, cl, hdi);

            if (is_contact_vector.at(cl_count)==1)
            {
                auto gr  = make_hho_contact_scalar_laplacian(msh, cl, hdi, bnd);
                Ah  = gr.second + stab;
            }
            else
            {
                auto gr  = make_scalar_hho_laplacian(msh, cl, hdi);
                Ah  = gr.second + stab;
            }

            //Erase this function, since it coulb be used to take data from full_sol
            //dsol has zero dirichlet conditions so nothing is computed in this faces/
            du_full = assembler.take_local_data_increment(msh, cl, dsol);
            diff_sol.block(cell_ofs, 0, num_total_dofs ,1) = du_full;

            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            matrix_type mass  = make_mass_matrix(msh, cl, cb);//, hdi.cell_degree());
            vector_type u_diff = du_full.block(0, 0, cbs, 1);

            H1_increment += du_full.dot(Ah * du_full);
            L2_increment += u_diff.dot(mass * u_diff);

            cl_count++;
        }

        full_sol += diff_sol;

        std::cout << "  "<< iter << "  "<< std::scientific<< std::setprecision(3);
        std::cout << std::sqrt(H1_increment)<< "   "<< std::sqrt(L2_increment)<< "  "<< std::endl;

        if( std::sqrt(H1_increment)  < tol)
        {
            std::ofstream efs("solution_whho_cnew_i" + tostr(iter) + ".dat");

            if(!efs.is_open())
                std::cout<< "Error opening file"<<std::endl;


            T H1_error  = 0.0 ;
            T L2_error  = 0.0 ;


            auto cl_count = 0;
            for(auto& cl : msh)
            {
                const auto cell_ofs = offset_vector.at(cl_count);
                const auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;

                vector_type  u_full   = full_sol.block(cell_ofs, 0, num_total_dofs, 1);
                matrix_type  Ah  = matrix_type::Zero(num_total_dofs, num_total_dofs);

                auto stab = make_scalar_hdg_stabilization(msh, cl, hdi);


                if (is_contact_vector.at(cl_count)==1)
                {
                    auto gr  = make_hho_contact_scalar_laplacian(msh, cl, hdi, bnd);
                    Ah  = gr.second + stab;
                }
                else
                {
                    auto gr  = make_scalar_hho_laplacian(msh, cl, hdi);
                    Ah  = gr.second + stab;
                }

                vector_type realsol = project_function(msh, cl, hdi, sol_fun, 2);

                vector_type diff = realsol - u_full;
                H1_error += diff.dot(Ah*diff);

                auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
                matrix_type mass  = make_mass_matrix(msh, cl, cb, hdi.cell_degree());

                vector_type u_diff = diff.block(0, 0, cbs, 1);
                L2_error += u_diff.dot(mass * u_diff);

                vector_type  du_full = assembler.take_local_data_increment(msh, cl, dsol);

                //plot
                auto bar = barycenter(msh, cl);
                efs << bar.x() << " " << bar.y() <<" "<< u_full(0) << " ";
                efs<< du_full(0)<< std::endl;

                auto fcs = faces(msh, cl);
                auto face_count = 0;
                for(auto fc: fcs)
                {
                    auto offset_face = cbs + face_count * fbs;
                    auto barf = barycenter(msh, fc);
                    efs << barf.x() << " " << barf.y() <<" "<< u_full(offset_face);
                    efs << " "<< du_full(offset_face) << std::endl;
                    face_count++;
                }
                cl_count++;

            }
            efs.close();

            std::cout << "  "<< iter << "  "<< std::sqrt(H1_error) << "   -    "<< std::scientific<< std::setprecision(3);
            std::cout << std::sqrt(H1_increment)<< "   "<< std::sqrt(L2_increment)<< "  "<< std::endl;

            //std::cout << "(H1, L2) : "<< std::sqrt(H1_error) <<" - " <<std::sqrt(L2_error) << std::endl;
            //if( std::sqrt(H1_increment)  < tol)
                return full_sol;
        }

    }
    return full_sol;
}

template<typename Mesh>
std::pair<typename Mesh::coordinate_type, typename Mesh::coordinate_type>
run_signorini_unknown( Mesh& msh,
     const algorithm_parameters<typename Mesh::coordinate_type>& ap)
{
    typedef typename Mesh::point_type  point_type;
    using T =  typename Mesh::coordinate_type;

    renumber_boundaries(msh);
    dump_to_matlab(msh,"mesh.m");

    std::cout << "mesh size:" <<average_diameter(msh) << std::endl;

    auto force = [](const point_type& p) -> T {
        return - 2.* M_PI *  std::sin(2. * M_PI * p.x());
    };

    auto zero_fun = [](const point_type& p) -> T {
        return 0.;
    };

    typedef disk::scalar_boundary_conditions<Mesh> boundary_type;
    boundary_type  bnd(msh);

    /*--------------------------------------------------------------------
    *        __1__
    *   4   |     | 2
    *       |_____|
    *          3
    *-------------------------------------------------------------------*/

    bnd.addDirichletBC(disk::DIRICHLET,1, zero_fun);
    bnd.addNeumannBC(disk::NEUMANN,  2, zero_fun);
    bnd.addNeumannBC(disk::NEUMANN,  4, zero_fun);
    bnd.addContactBC(disk::SIGNORINI,3);

    switch (ap.solver)
    {
        case EVAL_ON_FACES:
            return solve_faces(msh, force, zero_fun, ap, bnd);
            break;
        case EVAL_IN_CELLS_FULL:
            return solve_cells_full(msh, force, zero_fun, ap, bnd);
            break;
        default:
            throw std::invalid_argument("Invalid solver");
    }
}

template<typename Mesh>
std::pair<typename Mesh::coordinate_type, typename Mesh::coordinate_type>
run_signorini_analytical(Mesh& msh,
     const algorithm_parameters<typename Mesh::coordinate_type>& ap)
{
    typedef typename Mesh::point_type  point_type;
    using T =  typename Mesh::coordinate_type;

    renumber_boundaries(msh, 0., 1., -1., -1.);

    dump_to_matlab(msh,"mesh.m");

    auto force = [](const point_type& p) -> T {
        return 0.;
    };

    auto zero_fun = [](const point_type& p) -> T {
        return 0.;
    };

    auto left = [](const point_type& p) -> T {
        T radio = std::sqrt(p.x()*p.x() + p.y()*p.y());
        T theta = std::atan2(p.y(), p.x());
        T sintcos = std::sin(theta) *std::cos(4.5 *theta);
        T sincost = std::sin(4.5 *theta) * std::cos(theta);
         return  -4.5 * std::pow(radio, 3.5) *( sincost - sintcos);
    };

    auto right = [](const point_type& p) -> T {
        T radio = std::sqrt(p.x()*p.x() + p.y()*p.y());
        T theta = std::atan2(p.y(), p.x());
        T sintcos = std::sin(theta) *std::cos(4.5 *theta);
        T sincost = std::sin(4.5 *theta) * std::cos(theta);
         return  4.5 * std::pow(radio, 3.5) *( sincost - sintcos);
    };

    auto myatan = [](const  point_type& p) -> T {

        T v = 0;

        if ( std::abs(p.y()) < 1.e-8 & p.x() < 0)
            v = -M_PI;
        else
            v = std::atan2( p.y(), p.x());

        return v;
    };

    auto realmod = [](const T& a, const T& b) -> T {

        T result = fmod(a, b);
        return result >= 0 ? result : result + b;
    };


    auto fun = [&](const point_type& p) -> T {
        T radio = std::sqrt(p.x()*p.x() + p.y()*p.y());
        T theta = myatan(p);

        return  -std::pow(radio, 5.5) * std::sin(5.5 *theta);
    };

    typedef disk::scalar_boundary_conditions<Mesh> boundary_type;
    boundary_type  bnd(msh);


    bnd.addDirichletBC(disk::DIRICHLET,4, fun);
    bnd.addDirichletBC(disk::DIRICHLET,2, fun);
    bnd.addDirichletBC(disk::DIRICHLET,3, fun);
    bnd.addContactBC(disk::SIGNORINI,1);

    switch (ap.solver)
    {
        case EVAL_ON_FACES:
            return solve_faces(msh, force, fun, ap, bnd);
            break;
        case EVAL_IN_CELLS_FULL:
            return solve_cells_full(msh, force, fun, ap, bnd);
            break;

        default:
            throw std::invalid_argument("Invalid trace-type. Choose for face-based (f) or cell-based (l) trace versions.");
    }
}

template<typename Mesh>
std::pair<typename Mesh::coordinate_type, typename Mesh::coordinate_type>
run_signorini(  Mesh& msh, const algorithm_parameters<typename Mesh::coordinate_type>& ap,
                const bool& run_exact)
{
    if(run_exact)
        return run_signorini_analytical(msh, ap);
    else
        return run_signorini_unknown(msh, ap);
}
