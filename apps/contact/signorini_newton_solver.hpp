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

#include "revolution/bases"
#include "revolution/quadratures"
#include "revolution/methods/hho"


#include "output/silo.hpp"
#include "common.hpp"
#include "solvers/solver.hpp"

template<typename Mesh, typename Function>
auto
solve_faces(const Mesh&  msh, const Function& rhs_fun,
    const algorithm_parameters<typename Mesh::scalar_type>& ap,
    const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd)
{
    using T = typename Mesh::scalar_type;

    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    matrix_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>                 vector_type;

    auto is_contact_vector = make_is_contact_vector(msh, bnd);

    hho_degree_info      hdi(ap.degree + 1, ap.degree); //Also allow (degree + 1, degree)
    std::cout << " * cell degree :"<< hdi.cell_degree() << std::endl;
    std::cout << " * face degree :"<< hdi.face_degree() << std::endl;

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                    - fbs*msh.boundary_faces_size() ;

    auto offset_vector = full_offset(msh, hdi);

    dynamic_vector<T>  full_sol = dynamic_vector<T>::Zero(num_full_dofs);
    auto max_iter = 1000;
    auto tol = 1.e-6;

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        auto cl_count = 0;
        auto assembler = make_diffusion_assembler2(msh, hdi, bnd);

        for (auto& cl : msh)
        {
            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs  = cbs + howmany_faces(msh, cl) * fbs;
            vector_type  u_full  = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);

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

            auto sc = diffusion_static_condensation_compute_full(msh, cl, hdi, A, b);
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

        dynamic_vector<T> dsol = dynamic_vector<T>::Zero(systsz);

        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = true;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, dsol);

        T H1_error = 0.0 ;
        T L2_error = 0.0 ;

        dynamic_vector<T> diff_sol = dynamic_vector<T>::Zero(num_full_dofs);

        cl_count = 0;

        for (auto& cl : msh)
        {
            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;
            vector_type  u_full = full_sol.block(cell_ofs, 0, num_total_dofs, 1);


            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);

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
            vector_type du_faces = assembler.take_local_data(msh, cl, dsol);
            vector_type du_full  =
                diffusion_static_condensation_recover(msh, cl, hdi, A, cell_rhs, du_faces);

            diff_sol.block(cell_ofs, 0, num_total_dofs ,1) = du_full;

            H1_error += du_full.dot(A * du_full);

            matrix_type mass  = make_mass_matrix(msh, cl, cb);//, hdi.cell_degree());

            vector_type u_diff = du_full.block(0, 0, cbs, 1);
            L2_error += u_diff.dot(mass * u_diff);


            cl_count++;
        }

        full_sol += diff_sol;

        std::cout << "  "<< iter << "  "<< std::sqrt(H1_error) << "  "<< std::sqrt(L2_error) << std::endl;
        if( std::sqrt(H1_error)  < tol)
        {
            std::ofstream efs("solution_whho_faces.dat");

            if(!efs.is_open())
                std::cout<< "Error opening file"<<std::endl;

            auto cl_count = 0;
            for(auto& cl : msh)
            {
                auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
                auto cell_ofs = offset_vector.at(cl_count++);
                vector_type u_bar = full_sol.block(cell_ofs, 0, 1, 1);
                auto bar = barycenter(msh, cl);
                efs << bar.x() << " " << bar.y() <<" "<< u_bar(0) << std::endl;
            }

            efs.close();
            return 0;
        }

        save_data(full_sol, "faces_output.dat");
    }
    return 1;
}



template<typename Mesh, typename Function>
auto
solve_faces_borrar(const Mesh&  msh, const Function& rhs_fun,
    const algorithm_parameters<typename Mesh::scalar_type>& ap,
    const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd)
{
    using T = typename Mesh::scalar_type;

    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    matrix_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>                 vector_type;

    auto is_contact_vector = make_is_contact_vector(msh, bnd);

    hho_degree_info      hdi(ap.degree + 1, ap.degree); //Also allow (degree + 1, degree)
    std::cout << " * cell degree :"<< hdi.cell_degree() << std::endl;
    std::cout << " * face degree :"<< hdi.face_degree() << std::endl;

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                    - fbs*msh.boundary_faces_size() ;

    auto offset_vector = full_offset(msh, hdi);

    dynamic_vector<T>  full_sol = dynamic_vector<T>::Zero(num_full_dofs);

    //Starting newton with solution using faces
    solve_faces(msh, rhs_fun, ap, bnd);
    full_sol = read_data<T>("faces_output.dat");

    auto max_iter = 1000;
    auto tol = 1.e-6;

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        auto cl_count = 0;
        auto assembler = make_diffusion_assembler2(msh, hdi, bnd);

        for (auto& cl : msh)
        {
            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs  = cbs + howmany_faces(msh, cl) * fbs;
            vector_type  u_full  = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);

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

            auto sc = diffusion_static_condensation_compute_full(msh, cl, hdi, A, b);
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

        dynamic_vector<T> dsol = dynamic_vector<T>::Zero(systsz);

        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = true;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, dsol);

        T H1_error = 0.0 ;
        T L2_error = 0.0 ;

        dynamic_vector<T> diff_sol = dynamic_vector<T>::Zero(num_full_dofs);
        std::vector<std::tuple<T, T, T>> residue_vector(msh.cells_size());
        std::vector<T> L2_vector(msh.cells_size());

        cl_count = 0;

        for (auto& cl : msh)
        {
            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;
            vector_type  u_full = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);

            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            matrix_type Ah  = gr.second + stab;

            vector_type Bnegative  = vector_type::Zero(num_total_dofs);
            matrix_type Anitsche   = matrix_type::Zero(num_total_dofs, num_total_dofs);
            matrix_type Aheaviside = matrix_type::Zero(num_total_dofs, num_total_dofs);

            if (is_contact_vector.at(cl_count))
            {
                std::cout << "u_full ("<< cl_count << ") = [";
                for(size_t i = 0; i < u_full.size(); i++)
                    std::cout << u_full(i)<< " ";
                std::cout << "]" << std::endl;

                Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                Bnegative  = make_hho_negative_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
                Aheaviside = make_hho_heaviside_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
            }

            matrix_type A =   Ah - Anitsche + Aheaviside;
            vector_type b = -(Ah - Anitsche) * u_full - Bnegative;
            b.block(0, 0, cbs, 1) += Lh;


            vector_type cell_rhs = b.block(0, 0, cbs, 1);
            vector_type du_faces = assembler.take_local_data(msh, cl, dsol);
            vector_type du_full  =
                diffusion_static_condensation_recover(msh, cl, hdi, A, cell_rhs, du_faces);

            diff_sol.block(cell_ofs, 0, num_total_dofs ,1) = du_full;

            H1_error += du_full.dot(A * du_full);

            matrix_type mass  = make_mass_matrix(msh, cl, cb);//, hdi.cell_degree());

            vector_type u_diff = du_full.block(0, 0, cbs, 1);
            L2_error += u_diff.dot(mass * u_diff);
            L2_vector.at(cl_count) = u_diff.dot(mass * u_diff);

            vector_type bb_face = -Ah * u_full + Anitsche * u_full - Bnegative; // + bb_cell_all.block(cbs, 0, fcs.size() * fbs, 1);
            vector_type bb_cell = Lh;

            residue_vector.at(cl_count) = std::make_tuple( bb_cell.norm(), bb_face.norm(), b.norm());

            //#if 0
            if (is_contact_vector.at(cl_count))
            {
                std::cout << "du_full ("<< cl_count << ") = [";
                for(size_t i = 0; i < du_full.size(); i++)
                    std::cout << du_full(i)<< " ";
                std::cout << "]" << std::endl;
            }
            cl_count++;
        }

        full_sol += diff_sol;

        std::cout << "  "<< iter << "  "<< std::sqrt(H1_error) << "  "<< std::sqrt(L2_error) << std::endl;
        if( std::sqrt(H1_error)  < tol)
        {
            std::ofstream efs("solution_whho_faces_borrar.dat");

            if(!efs.is_open())
                std::cout<< "Error opening file"<<std::endl;

            auto cl_count = 0;
            for(auto& cl : msh)
            {
                auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
                auto cell_ofs = offset_vector.at(cl_count);
                vector_type u_bar = full_sol.block(cell_ofs, 0, 1, 1);
                auto bar = barycenter(msh, cl);
                auto res = residue_vector.at(cl_count);
                efs << bar.x() << " " << bar.y() <<" "<< u_bar(0) << " ";
                efs <<  std::setprecision(16)<< L2_vector.at(cl_count)<< " ";
                efs << std::get<0>(res)<< " "<< std::get<1>(res)<< " "<< std::get<2>(res)<< " "<< std::endl;
                cl_count++;
            }

            efs.close();
            return 0;
        }

        save_data(full_sol, "faces_output.dat");
    }
    return 1;
}


template<typename Mesh, typename Function>
auto
solve_cells(const Mesh&  msh, const Function& rhs_fun,
    const algorithm_parameters<typename Mesh::scalar_type>& ap,
    const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd)
{
    using T = typename Mesh::scalar_type;

    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    matrix_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>                 vector_type;

    auto is_contact_vector = make_is_contact_vector(msh, bnd);
    size_t vcount = 0;
    std::cout << "Contact cells: " << std::endl;
    for(auto v : is_contact_vector)
    {
        if(v)
            std::cout << vcount <<" ";
        vcount++;
    }
    std::cout << std::endl;


    hho_degree_info      hdi(ap.degree +1, ap.degree); //Not allow (degree, degree)

    std::cout << " * cell degree :"<< hdi.cell_degree() << std::endl;
    std::cout << " * face degree :"<< hdi.face_degree() << std::endl;

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                    - fbs*msh.boundary_faces_size() ;

    auto offset_vector = full_offset(msh, hdi);

    dynamic_vector<T>  full_sol = dynamic_vector<T>::Zero(num_full_dofs);
    auto max_iter = 15000;
    auto tol = 1.e-8;

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        auto cl_count = 0;
        //ERROR**** This doesnt work for signorini with cells, since vF dofs have to be taken out from the matrix
        auto assembler = make_diffusion_assembler2(msh, hdi, bnd);

        for (auto& cl : msh)
        {
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());

            if (is_contact_vector.at(cl_count) == 1)
            {
                auto cell_ofs = offset_vector.at(cl_count);
                auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;
                vector_type  u_full = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

                auto gr   = make_hho_contact_scalar_laplacian(msh, cl, hdi, bnd);
                auto stab = make_hdg_contact_stabilization(msh, cl, hdi, bnd);

                matrix_type Ah  = gr.second + stab;
                vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());

                matrix_type  Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                vector_type  Bnegative  = make_hho_negative(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
                matrix_type  Aheaviside = make_hho_heaviside(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);

                matrix_type A =   Ah - Anitsche + Aheaviside;
                vector_type b = -(Ah - Anitsche) * u_full - Bnegative;
                b.block(0, 0, cbs, 1) += Lh;

                auto sc = diffusion_static_condensation_compute_full(msh, cl, hdi, A, b);
                assembler.assemble(msh, cl, sc.first, sc.second);

            }
            else
            {
                auto gr   = make_hho_scalar_laplacian(msh, cl, hdi);
                auto stab = make_hdg_scalar_stabilization(msh, cl, hdi);

                vector_type Lh = make_rhs(msh, cl, cb, rhs_fun);//, hdi.cell_degree());
                matrix_type Ah = gr.second + stab;

                auto sc = diffusion_static_condensation_compute(msh, cl, hdi, Ah, Lh);
                assembler.assemble(msh, cl, sc.first, sc.second);
            }
            cl_count++;
        }

        assembler.impose_neumann_boundary_conditions(msh, bnd);
        assembler.finalize();

        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();


        dynamic_vector<T> dsol = dynamic_vector<T>::Zero(systsz);

        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = true;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, dsol);

        T error  = 0.0 ;
        cl_count = 0;
        dynamic_vector<T> diff_sol = dynamic_vector<T>::Zero(num_full_dofs);

        for (auto& cl : msh)
        {
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            const auto cell_ofs = offset_vector.at(cl_count);
            const auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;

            vector_type  u_full = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            if (is_contact_vector.at(cl_count)==1)
            {
                auto gr  = make_hho_contact_scalar_laplacian(msh, cl, hdi, bnd);
                auto stab = make_hdg_contact_stabilization(msh, cl, hdi, bnd);

                matrix_type Ah  = gr.second + stab;
                vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());

                matrix_type  Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                vector_type  Bnegative  = make_hho_negative(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
                matrix_type  Aheaviside = make_hho_heaviside(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);

                matrix_type A =   Ah - Anitsche + Aheaviside;
                vector_type b = -(Ah - Anitsche) * u_full - Bnegative;
                b.block(0, 0, cbs, 1) += Lh;

                vector_type cell_rhs = b.block(0, 0, cbs, 1);
                vector_type du_faces = assembler.take_local_data(msh, cl, dsol);
                vector_type du_full  =
                    diffusion_static_condensation_recover(msh, cl, hdi, A, cell_rhs, du_faces);

                diff_sol.block(cell_ofs, 0, num_total_dofs ,1) = du_full;
                error += du_full.dot(A * du_full);

                //#if 0
                std::cout << "du_full ("<< cl_count << ") = [";
                for(size_t i = 0; i < du_full.size(); i++)
                    std::cout << du_full(i)<< " ";
                std::cout << "]" << std::endl;
            }
            else
            {
                auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
                auto stab   = make_hdg_scalar_stabilization(msh, cl, hdi);


                vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun);//, hdi.cell_degree());
                matrix_type Ah  = gr.second + stab;

                vector_type du_faces = assembler.take_local_data(msh, cl, dsol);
                vector_type du_full  =
                    diffusion_static_condensation_recover(msh, cl, hdi, Ah, Lh, du_faces);

                diff_sol.block(cell_ofs, 0, num_total_dofs ,1) = du_full;
                error += du_full.dot(Ah * du_full);
            }
            cl_count++;
        }

        full_sol += diff_sol;

        std::cout << "  "<< iter << "  "<< std::sqrt(error)<< std::endl;
        if( std::sqrt(error)  < tol)
        {
            std::ofstream efs("solution_whho_cell.dat");

            if(!efs.is_open())
                std::cout<< "Error opening file"<<std::endl;

            auto cl_count = 0;
            for(auto& cl : msh)
            {
                auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
                auto cell_ofs = offset_vector.at(cl_count++);
                vector_type u_bar = full_sol.block(cell_ofs, 0, 1, 1);
                auto bar = barycenter(msh, cl);
                efs << bar.x() << " " << bar.y() <<" "<< u_bar(0) << std::endl;
            }

            efs.close();
            return 0;
        }
    }
    return 1;
}

template<typename Mesh, typename Function>
auto
solve_cells_full(const Mesh&  msh, const Function& rhs_fun,
    const algorithm_parameters<typename Mesh::scalar_type>& ap,
    const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd)
{
    std::cout << "Im in CELLS FULL" << std::endl;
    using T = typename Mesh::scalar_type;

    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    matrix_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>                 vector_type;

    auto is_contact_vector = make_is_contact_vector(msh, bnd);

    size_t vcount = 0;
    std::cout << "Contact cells: " << std::endl;
    for(auto v : is_contact_vector)
    {
        if(v)
            std::cout << vcount <<" ";
        vcount++;
    }
    std::cout << std::endl;


    hho_degree_info      hdi(ap.degree +1, ap.degree); //Not allow (degree, degree)
    std::cout << " * cell degree :"<< hdi.cell_degree() << std::endl;
    std::cout << " * face degree :"<< hdi.face_degree() << std::endl;

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                    - fbs*msh.boundary_faces_size() ;

    auto offset_vector = full_offset(msh, hdi);

    dynamic_vector<T>  full_sol = dynamic_vector<T>::Zero(num_full_dofs);

    //__________________________________________________________________________
    //Starting newton with solution using faces
    solve_faces(msh, rhs_fun, ap, bnd);
    full_sol = read_data<T>("faces_output.dat");


    {
        std::ofstream efs("solution_whho_c_input.dat");

        if(!efs.is_open())
            std::cout<< "Error opening file"<<std::endl;

        auto cl_count = 0;
        for(auto& cl : msh)
        {
            auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
            auto cell_ofs = offset_vector.at(cl_count);
            vector_type u_bar = full_sol.block(cell_ofs, 0, 1, 1);
            auto bar = barycenter(msh, cl);


            efs << bar.x() << " " << bar.y() <<" "<< u_bar(0) <<" "<< std::endl;
            cl_count++;
        }


        efs.close();
    }
    //__________________________________________________________________________

    auto max_iter = 3;
    auto tol = 1.e-6;

    for(size_t iter = 0; iter < max_iter; iter++)
    {

        auto assembler = make_contact_full_assembler(msh, hdi, bnd);

        auto cl_count = 0;
        std::vector<std::tuple<T, T, T>> residue_vector(msh.cells_size());
        auto dofs = 0;
        for (auto& cl : msh)
        {
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());

            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;

            vector_type  u_full = full_sol.block(dofs, 0, num_total_dofs, 1);//cell_ofs, 0, num_total_dofs, 1);
            dofs += num_total_dofs;

            if (is_contact_vector.at(cl_count) == 1)
            {

                auto fcs = faces(msh, cl);
                for(size_t fi = 0; fi < fcs.size(); fi++)
                {
                    auto face_id = msh.lookup(fcs[fi]);
                    if(bnd.is_contact_face(face_id))
                        u_full.block(cbs + fbs * fi, 0, fbs, 1) = vector_type::Zero(fbs);
                }

                #if 0
                std::cout << "u_full ("<< cl_count << ") = [";
                for(size_t i = 0; i < u_full.size(); i++)
                    std::cout << u_full(i)<< " ";
                std::cout << "]" << std::endl;
                #endif

                auto gr   = make_hho_contact_scalar_laplacian(msh, cl, hdi, bnd);
                auto stab = make_hdg_contact_stabilization(msh, cl, hdi, bnd);

                matrix_type Ah  = gr.second + stab;
                vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun);//, hdi.cell_degree());

                matrix_type  Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                //vector_type  Bnegative  = make_hho_negative_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full, 1.);
                matrix_type  Aheaviside = make_hho_heaviside_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full, 1.);

                vector_type  Bnegative  = make_hho_negative(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
                //matrix_type  Aheaviside = make_hho_heaviside(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);

                for(size_t fi = 0; fi < fcs.size(); fi++)
                {
                    auto face_id = msh.lookup(fcs[fi]);
                    if(bnd.is_contact_face(face_id))
                    {
                        Ah.block(cbs + fbs * fi, 0, fbs, Ah.cols()) = matrix_type::Zero(fbs, Ah.cols());
                        Ah.block(0, cbs + fbs * fi, Ah.rows(), fbs) = matrix_type::Zero(Ah.rows(), fbs);

                        Anitsche.block(cbs + fbs * fi, 0, fbs, Anitsche.cols()) = matrix_type::Zero(fbs, Anitsche.cols());
                        Anitsche.block(0, cbs + fbs * fi, Anitsche.rows(), fbs) = matrix_type::Zero(Anitsche.rows(), fbs);

                        Aheaviside.block(cbs + fbs * fi, 0, fbs, Aheaviside.cols()) = matrix_type::Zero(fbs, Aheaviside.cols());
                        Aheaviside.block(0, cbs + fbs * fi, Aheaviside.rows(), fbs) = matrix_type::Zero(Aheaviside.rows(), fbs);

                        Bnegative.block(cbs + fbs * fi, 0, fbs, 1) = vector_type::Zero(fbs);
                    }
                }
                //Original
                matrix_type A =   Ah - Anitsche + Aheaviside;
                vector_type b = -(Ah - Anitsche) * u_full - Bnegative;
                b.block(0, 0, cbs, 1) += Lh;
                //assert(b.rows() == num_total_dofs);
                //assembler.assemble(msh, cl, A, b);

                //Test
                //matrix_type A =   Ah - Anitsche + Aheaviside;
                //vector_type b  = vector_type::Zero(num_total_dofs); // - (Ah * u_full).block(0, 0, cbs, 1) + Lh;
                //b.block(0, 0, cbs, 1) += Lh;

                for(size_t fi = 0; fi < fcs.size(); fi++)
                {
                    auto face_id = msh.lookup(fcs[fi]);
                    if(bnd.is_contact_face(face_id))
                    {
                        A.block(cbs + fbs * fi, 0, fbs, A.cols()) = matrix_type::Zero(fbs, A.cols());
                        A.block(0, cbs + fbs * fi, A.rows(), fbs) = matrix_type::Zero(A.rows(), fbs);

                        b.block(cbs + fbs * fi, 0, fbs, 1) = vector_type::Zero(fbs);
                    }
                }

                vector_type bb_face = -Ah * u_full + Anitsche * u_full - Bnegative; // + bb_cell_all.block(cbs, 0, fcs.size() * fbs, 1);
                vector_type bb_cell = Lh;

                residue_vector.at(cl_count) = std::make_tuple( bb_cell.norm(), bb_face.norm(), b.norm());

                assembler.assemble(msh, cl, A, b);
                //#endif
                //***************************************************************


            }
            else
            {
                auto gr   = make_hho_scalar_laplacian(msh, cl, hdi);
                auto stab = make_hdg_scalar_stabilization(msh, cl, hdi);

                matrix_type Ah = gr.second + stab;
                vector_type Lh = make_rhs(msh, cl, cb, rhs_fun);//, hdi.cell_degree());

                //original
                vector_type b = -Ah * u_full;
                b.block(0, 0, cbs, 1) += Lh;

                //vector_type b  = vector_type::Zero(num_total_dofs);
                assembler.assemble(msh, cl, Ah, b);

                vector_type bb_face = -Ah * u_full;
                vector_type bb_cell = Lh;

                residue_vector.at(cl_count) = std::make_tuple( bb_cell.norm(), bb_face.norm(), b.norm());
            }
            cl_count++;
        }

        //assembler.impose_neumann_boundary_conditions(msh, bnd);
        assembler.finalize();

        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        //std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
        //std::cout << "Mesh faces: " << msh.faces_size() << std::endl;
        //std::cout << "Dofs: " << systsz << std::endl;

        dynamic_vector<T> dsol = dynamic_vector<T>::Zero(systsz);

        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = true;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, dsol);


        dump_sparse_matrix(assembler.LHS, "Amatrix_i" + tostr(iter) + ".dat");
        T H1_error  = 0.0 ;
        T L2_error  = 0.0 ;

        cl_count = 0;
        dynamic_vector<T> diff_sol = dynamic_vector<T>::Zero(num_full_dofs);
        std::vector<T> L2_vector(msh.cells_size());

        for (auto& cl : msh)
        {
            const auto cell_ofs = offset_vector.at(cl_count);
            const auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;

            vector_type  u_full = full_sol.block(cell_ofs, 0, num_total_dofs, 1);
            matrix_type  Ah  = matrix_type::Zero(num_total_dofs, num_total_dofs);
            vector_type  du_full  = vector_type::Zero(num_total_dofs);

            if (is_contact_vector.at(cl_count)==1)
            {
                auto gr   = make_hho_contact_scalar_laplacian(msh, cl, hdi, bnd);
                auto stab = make_hdg_contact_stabilization(msh, cl, hdi, bnd);

                Ah  = gr.second + stab;

                auto fcs = faces(msh, cl);

                for(size_t fi = 0; fi < fcs.size(); fi++)
                {
                    auto face_id = msh.lookup(fcs[fi]);
                    if(bnd.is_contact_face(face_id))
                    {
                        Ah.block(cbs + fbs * fi, 0, fbs, Ah.cols()) = matrix_type::Zero(fbs, Ah.cols());
                        Ah.block(0, cbs + fbs * fi, Ah.rows(), fbs) = matrix_type::Zero(Ah.rows(), fbs);

                        u_full.block(cbs + fbs * fi, 0, fbs, 1) = vector_type::Zero(fbs);
                    }
                }

                du_full = assembler.take_local_data(msh, cl, dsol);

                #if 0
                std::cout << "du_full ("<< cl_count << ") = [";
                for(size_t i = 0; i < du_full.size(); i++)
                    std::cout << du_full(i)<< " ";
                std::cout << "]" << std::endl;
                #endif

            }
            else
            {
                auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
                auto stab   = make_hdg_scalar_stabilization(msh, cl, hdi);

                Ah  = gr.second + stab;
                du_full = assembler.take_local_data(msh, cl, dsol);
            }


            diff_sol.block(cell_ofs, 0, num_total_dofs ,1) = du_full;
            H1_error += du_full.dot(Ah * du_full);

            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            matrix_type mass  = make_mass_matrix(msh, cl, cb);//, hdi.cell_degree());

            vector_type u_diff = du_full.block(0, 0, cbs, 1);
            L2_error += u_diff.dot(mass * u_diff);

            L2_vector.at(cl_count) = u_diff.dot(mass * u_diff);

            cl_count++;
        }


        full_sol += diff_sol;

        std::cout << "  "<< iter << "  "<< std::sqrt(H1_error)<< "   "<< std::sqrt(L2_error)<< std::endl;
        //if( std::sqrt(error)  < tol)
        {
            std::ofstream efs("solution_whho_cfull_i" + tostr(iter) + ".dat");

            if(!efs.is_open())
                std::cout<< "Error opening file"<<std::endl;

            auto cl_count = 0;
            for(auto& cl : msh)
            {
                auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
                auto cell_ofs = offset_vector.at(cl_count);
                vector_type u_bar = full_sol.block(cell_ofs, 0, 1, 1);
                auto bar = barycenter(msh, cl);


                auto res =  residue_vector.at(cl_count);
                efs << bar.x() << " " << bar.y() <<" "<< u_bar(0) <<" ";
                efs << std::setprecision(16)<< L2_vector.at(cl_count)<< " ";
                efs << std::get<0>(res)<< " "<< std::get<1>(res)<< " "<< std::get<2>(res)<< " "<< std::endl;
                cl_count++;
            }


            efs.close();

            if( std::sqrt(H1_error)  < tol)
                return 0;
        }
    }
    return 1;
}


template<typename Mesh, typename Function>
dynamic_vector<typename Mesh::scalar_type>
diffusion(const Mesh&  msh, const Function& rhs_fun,
        const algorithm_parameters<typename Mesh::scalar_type>& ap,
        const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd)
{
    using T = typename Mesh::scalar_type;

    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    matrix_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>                 vector_type;

    auto is_contact_vector = make_is_contact_vector(msh, bnd);
    hho_degree_info     hdi(ap.degree +1, ap.degree); //Also allow (degree + 1, degree)
    std::cout << " * cell degree :"<< hdi.cell_degree() << std::endl;
    std::cout << " * face degree :"<< hdi.face_degree() << std::endl;

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                    - fbs*msh.boundary_faces_size() ;

    auto offset_vector = full_offset(msh, hdi);

    dynamic_vector<T>  diff_sol = dynamic_vector<T>::Zero(num_full_dofs);

    auto cl_count = 0;
    auto assembler = make_diffusion_assembler2(msh, hdi, bnd);

    for (auto& cl : msh)
    {
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
        auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);

        vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
        matrix_type Ah  = gr.second + stab;

        auto sc = diffusion_static_condensation_compute(msh, cl, hdi, Ah, Lh);
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

    dynamic_vector<T> diffusion_sol = dynamic_vector<T>::Zero(num_full_dofs);

    cl_count = 0;
    for (auto& cl : msh)
    {
        auto cell_ofs = offset_vector.at(cl_count);
        auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;

        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
        auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);

        vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
        matrix_type Ah  = gr.second + stab;

        vector_type u_faces = assembler.take_local_data(msh, cl, sol);
        vector_type u_full  =
        diffusion_static_condensation_recover(msh, cl, hdi, Ah, Lh, u_faces);

        diffusion_sol.block(cell_ofs, 0, num_total_dofs ,1) = u_full;
        cl_count++;
    }
    return diffusion_sol;
}

template<typename Mesh, typename Function>
auto
solve_cell_as_faces(const Mesh&  msh, const Function& rhs_fun,
    const algorithm_parameters<typename Mesh::scalar_type>& ap,
    const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd)
{
    using T = typename Mesh::scalar_type;

    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    matrix_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>                 vector_type;

    auto is_contact_vector = make_is_contact_vector(msh, bnd);
    hho_degree_info     hdi(ap.degree +1, ap.degree); //Also allow (degree + 1, degree)
    std::cout << " * cell degree :"<< hdi.cell_degree() << std::endl;
    std::cout << " * face degree :"<< hdi.face_degree() << std::endl;

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                        - fbs*msh.boundary_faces_size() ;

    auto offset_vector = full_offset(msh, hdi);

    auto max_iter = 13;
    auto tol = 1.e-6;

    dynamic_vector<T>  full_sol = dynamic_vector<T>::Zero(num_full_dofs);

    full_sol = diffusion(msh, rhs_fun, ap, bnd);

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        auto cl_count = 0;
        auto assembler = make_diffusion_assembler2(msh, hdi, bnd);

        for (auto& cl : msh)
        {
            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs  = cbs + howmany_faces(msh, cl) * fbs;
            vector_type  u_full  = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);

            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            matrix_type Ah  = gr.second + stab;

            vector_type Bnegative  = vector_type::Zero(num_total_dofs);
            matrix_type Anitsche   = matrix_type::Zero(num_total_dofs, num_total_dofs);
            matrix_type Aheaviside = matrix_type::Zero(num_total_dofs, num_total_dofs);

            if (is_contact_vector.at(cl_count) == 1)
            {
                Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                Bnegative  = make_hho_negative_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full, 1.);
                Aheaviside = make_hho_heaviside_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full, 1.);
            }

            matrix_type A =   Ah - Anitsche + Aheaviside;
            vector_type b = -(Ah - Anitsche) * u_full - Bnegative;
            b.block(0, 0, cbs, 1) += Lh;

            auto sc = diffusion_static_condensation_compute_full(msh, cl, hdi, A, b);
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

        dynamic_vector<T> dsol = dynamic_vector<T>::Zero(systsz);

        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = true;
        mkl_pardiso_ldlt(pparams, assembler.LHS, assembler.RHS, dsol);

        T error = 0.0 ;

        dynamic_vector<T> diff_sol = dynamic_vector<T>::Zero(num_full_dofs);

        cl_count = 0;

        for (auto& cl : msh)
        {
            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;
            vector_type  u_full = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);

            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            matrix_type Ah  = gr.second + stab;

            vector_type Bnegative  = vector_type::Zero(num_total_dofs);
            matrix_type Anitsche   = matrix_type::Zero(num_total_dofs, num_total_dofs);
            matrix_type Aheaviside = matrix_type::Zero(num_total_dofs, num_total_dofs);

            if (is_contact_vector.at(cl_count))
            {
                Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                Bnegative  = make_hho_negative_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full, 1.);
                Aheaviside = make_hho_heaviside_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full, 1.);
            }

            matrix_type A =   Ah - Anitsche + Aheaviside;
            vector_type b = -(Ah - Anitsche) * u_full - Bnegative;
            b.block(0, 0, cbs, 1) += Lh;

            vector_type cell_rhs = b.block(0, 0, cbs, 1);
            vector_type du_faces = assembler.take_local_data(msh, cl, dsol);
            vector_type du_full  =
            diffusion_static_condensation_recover(msh, cl, hdi, A, cell_rhs, du_faces);

            diff_sol.block(cell_ofs, 0, num_total_dofs ,1) = du_full;

            error += du_full.dot(A * du_full);
            cl_count++;
        }

        full_sol += diff_sol;

        std::cout << "  "<< iter << "  "<< std::sqrt(error)<< std::endl;
        if( std::sqrt(error)  < tol)
        {
            std::ofstream efs("solution_whho_cell_as_faces" + tostr(iter) + ".dat");
            if(!efs.is_open())
                std::cout<< "Error opening file"<<std::endl;

            auto cl_count = 0;
            for(auto& cl : msh)
            {
                auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
                auto cell_ofs = offset_vector.at(cl_count++);
                vector_type u_bar = full_sol.block(cell_ofs, 0, 1, 1);
                auto bar = barycenter(msh, cl);
                efs << bar.x() << " " << bar.y() <<" "<< u_bar(0) << std::endl;
            }

            efs.close();
            return 0;
        }

    }
    return 1;
}

template<typename Mesh, typename Function>
dynamic_vector<typename Mesh::scalar_type>
solve_cell_as_faces_test(const Mesh&  msh, const Function& rhs_fun,
    const algorithm_parameters<typename Mesh::scalar_type>& ap,
    const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd)
{
    std::cout << "INSIDE CELL AS FACES TEST WITH H = "<< average_diameter(msh) << std::endl;
    std::cout << "GAMMA = "<< ap.gamma_0 << std::endl;
    std::cout << "THETA = "<< ap.theta << std::endl;
    dump_to_matlab(msh,"mesh.m");
    using T = typename Mesh::scalar_type;

    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    matrix_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>                 vector_type;

    auto is_contact_vector = make_is_contact_vector(msh, bnd);
    hho_degree_info     hdi(ap.degree + 1, ap.degree); //Also allow (degree + 1, degree)
    //std::cout << " * cell degree :"<< hdi.cell_degree() << std::endl;
    //std::cout << " * face degree :"<< hdi.face_degree() << std::endl;

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                        - fbs*msh.boundary_faces_size() ;

    auto offset_vector = full_offset(msh, hdi);

    save_data(offset_vector, "ofs_int.data");

    auto max_iter = 13;
    auto tol = 1.e-6;

    dynamic_vector<T>  full_sol = dynamic_vector<T>::Zero(num_full_dofs);

    full_sol = diffusion(msh, rhs_fun, ap, bnd);

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        auto cl_count = 0;
        auto assembler = make_diffusion_assembler2(msh, hdi, bnd);

        for (auto& cl : msh)
        {
            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs  = cbs + howmany_faces(msh, cl) * fbs;
            vector_type  u_full  = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());

            vector_type Bnegative  = vector_type::Zero(num_total_dofs);
            matrix_type Ah   = matrix_type::Zero(num_total_dofs, num_total_dofs);
            matrix_type Anitsche   = matrix_type::Zero(num_total_dofs, num_total_dofs);
            matrix_type Aheaviside = matrix_type::Zero(num_total_dofs, num_total_dofs);

            if (is_contact_vector.at(cl_count) == 1)
            {
                auto gr     = make_hho_contact_scalar_laplacian(msh, cl, hdi, bnd);
                auto stab   = make_hdg_scalar_stabilization(msh, cl, hdi);

                Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                Bnegative  = make_hho_negative_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full, 1.);
                Aheaviside = make_hho_heaviside_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full, 1.);

                Ah  = gr.second + stab;
            }
            else
            {
                auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
                auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);

                Ah  = gr.second + stab;
            }

            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());

            matrix_type A =   Ah - Anitsche + Aheaviside;
            vector_type b = -(Ah - Anitsche) * u_full - Bnegative;
            b.block(0, 0, cbs, 1) += Lh;

            auto sc = diffusion_static_condensation_compute_full(msh, cl, hdi, A, b);
            assembler.assemble(msh, cl, sc.first, sc.second);
            cl_count++;
        }

        assembler.impose_neumann_boundary_conditions(msh, bnd);
        assembler.finalize();
    dump_to_matlab(msh,"mesh.m");
        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        //std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
        //std::cout << "Mesh faces: " << msh.faces_size() << std::endl;
        //std::cout << "Dofs: " << systsz << std::endl;

        dynamic_vector<T> dsol = dynamic_vector<T>::Zero(systsz);

        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = true;
        mkl_pardiso_ldlt(pparams, assembler.LHS, assembler.RHS, dsol);

        T error = 0.0 ;

        dynamic_vector<T> diff_sol = dynamic_vector<T>::Zero(num_full_dofs);

        cl_count = 0;

        for (auto& cl : msh)
        {
            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;
            vector_type  u_full = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);

            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            matrix_type Ah  = matrix_type::Zero(num_total_dofs, num_total_dofs);

            vector_type Bnegative  = vector_type::Zero(num_total_dofs);
            matrix_type Anitsche   = matrix_type::Zero(num_total_dofs, num_total_dofs);
            matrix_type Aheaviside = matrix_type::Zero(num_total_dofs, num_total_dofs);

            if (is_contact_vector.at(cl_count) == 1)
            {
                auto gr     = make_hho_contact_scalar_laplacian(msh, cl, hdi, bnd);
                auto stab   = make_hdg_scalar_stabilization(msh, cl, hdi);

                Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                //Bnegative  = make_hho_negative_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full, 1.);
                //Aheaviside = make_hho_heaviside_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full, 1.);

                Bnegative  = make_hho_negative(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
                Aheaviside = make_hho_heaviside(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);

                Ah  = gr.second + stab;
            }
            else
            {
                auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
                auto stab   = make_hdg_scalar_stabilization(msh, cl, hdi);

                Ah  = gr.second + stab;
            }

            matrix_type A =   Ah - Anitsche + Aheaviside;
            vector_type b = -(Ah - Anitsche) * u_full - Bnegative;
            b.block(0, 0, cbs, 1) += Lh;

            vector_type cell_rhs = b.block(0, 0, cbs, 1);
            vector_type du_faces = assembler.take_local_data(msh, cl, dsol);
            vector_type du_full  =
            diffusion_static_condensation_recover(msh, cl, hdi, A, cell_rhs, du_faces);

            diff_sol.block(cell_ofs, 0, num_total_dofs ,1) = du_full;

            error += du_full.dot(Ah * du_full);
            cl_count++;
        }

        full_sol += diff_sol;

        std::cout << "  "<< iter << "  "<< std::sqrt(error)<< std::endl;
        if( std::sqrt(error)  < tol)
        {
            std::ofstream efs("solution_whho_cell_as_faces_t_" + tostr(iter) + ".dat");
            if(!efs.is_open())
                std::cout<< "Error opening file"<<std::endl;

            auto cl_count = 0;
            for(auto& cl : msh)
            {
                auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
                auto cell_ofs = offset_vector.at(cl_count++);
                vector_type u_bar = full_sol.block(cell_ofs, 0, 1, 1);
                auto bar = barycenter(msh, cl);
                efs << bar.x() << " " << bar.y() <<" "<< u_bar(0) << std::endl;
            }

            efs.close();

            return full_sol;
        }

    }

    assert(false);
    return full_sol;
}


template<typename Mesh, typename Function>
auto
solve_with_parameter(const Mesh&  msh, const Function& rhs_fun,
    const algorithm_parameters<typename Mesh::scalar_type>& ap,
    const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
    const typename Mesh::coordinate_type& parameter)
{
    using T = typename Mesh::scalar_type;

    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    matrix_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>                 vector_type;

    std::cout << "INSIDE SOLVE_ WITH PARAMETER" << std::endl;
    hho_degree_info      hdi(ap.degree + 1, ap.degree); //Also allow (degree + 1, degree)
    std::cout << " * cell degree : "<< hdi.cell_degree() << std::endl;
    std::cout << " * face degree : "<< hdi.face_degree() << std::endl;
    std::cout << " * parameter   : "<< parameter  << std::endl;

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                    - fbs*msh.boundary_faces_size() ;

    auto offset_vector = full_offset(msh, hdi);

    dynamic_vector<T>  full_sol = dynamic_vector<T>::Zero(num_full_dofs);
    auto max_iter = 1000;
    auto tol = 1.e-6;

    auto is_contact_vector = make_is_contact_vector(msh, bnd);

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        auto cl_count = 0;
        auto assembler = make_diffusion_assembler2(msh, hdi, bnd);

        for (auto& cl : msh)
        {
            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs  = cbs + howmany_faces(msh, cl) * fbs;
            vector_type  u_full  = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
            //auto stab   = make_hdg_scalar_stabilization(msh, cl, hdi);

            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun);//, hdi.cell_degree());
            matrix_type Ah  = gr.second + stab;

            vector_type Bnegative  = vector_type::Zero(num_total_dofs);
            matrix_type Anitsche   = matrix_type::Zero(num_total_dofs, num_total_dofs);
            matrix_type Aheaviside = matrix_type::Zero(num_total_dofs, num_total_dofs);

            if (is_contact_vector.at(cl_count) == 1)
            {
                Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                Bnegative  = make_hho_negative_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full, parameter);
                Aheaviside = make_hho_heaviside_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full, parameter);
            }

            matrix_type A =   Ah - Anitsche + Aheaviside;
            vector_type b = -(Ah - Anitsche) * u_full - Bnegative;
            b.block(0, 0, cbs, 1) += Lh;

            auto sc = diffusion_static_condensation_compute_full(msh, cl, hdi, A, b);
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

        dynamic_vector<T> dsol = dynamic_vector<T>::Zero(systsz);

        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = true;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, dsol);

        T error = 0.0 ;

        dynamic_vector<T> diff_sol = dynamic_vector<T>::Zero(num_full_dofs);

        cl_count = 0;

        for (auto& cl : msh)
        {
            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;
            vector_type  u_full = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
            //auto stab   = make_hdg_scalar_stabilization(msh, cl, hdi);

            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun);//, hdi.cell_degree());
            matrix_type Ah  = gr.second + stab;

            vector_type Bnegative  = vector_type::Zero(num_total_dofs);
            matrix_type Anitsche   = matrix_type::Zero(num_total_dofs, num_total_dofs);
            matrix_type Aheaviside = matrix_type::Zero(num_total_dofs, num_total_dofs);

            if (is_contact_vector.at(cl_count))
            {
                Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                Bnegative  = make_hho_negative_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full, parameter);
                Aheaviside = make_hho_heaviside_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full, parameter);
            }

            matrix_type A =   Ah - Anitsche + Aheaviside;
            vector_type b = -(Ah - Anitsche) * u_full - Bnegative;
            b.block(0, 0, cbs, 1) += Lh;

            vector_type cell_rhs = b.block(0, 0, cbs, 1);
            vector_type du_faces = assembler.take_local_data(msh, cl, dsol);
            vector_type du_full  =
                diffusion_static_condensation_recover(msh, cl, hdi, A, cell_rhs, du_faces);

            diff_sol.block(cell_ofs, 0, num_total_dofs ,1) = du_full;

            error += du_full.dot(A * du_full);
            cl_count++;
        }

        full_sol += diff_sol;

        std::cout << "  "<< iter << "  "<< std::sqrt(error)<< std::endl;
        if( std::sqrt(error)  < tol)
        {
            std::ofstream efs("solution_whho_faces.dat");

            if(!efs.is_open())
                std::cout<< "Error opening file"<<std::endl;

            auto cl_count = 0;
            for(auto& cl : msh)
            {
                auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
                auto cell_ofs = offset_vector.at(cl_count++);
                vector_type u_bar = full_sol.block(cell_ofs, 0, 1, 1);
                auto bar = barycenter(msh, cl);
                efs << bar.x() << " " << bar.y() <<" "<< u_bar(0) << std::endl;
            }

            efs.close();
            return 0;
        }

    }
    return 1;
}


template<typename Mesh, typename Function>
auto
solve_with_parameter_full(const Mesh&  msh, const Function& rhs_fun,
    const algorithm_parameters<typename Mesh::scalar_type>& ap,
    const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
    const typename Mesh::coordinate_type& parameter)
{
    using T = typename Mesh::scalar_type;

    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    matrix_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>                 vector_type;

    std::cout << "INSIDE SOLVE_ WITH PARAMETER" << std::endl;
    hho_degree_info      hdi(ap.degree + 1, ap.degree); //Also allow (degree + 1, degree)
    std::cout << " * cell degree : "<< hdi.cell_degree() << std::endl;
    std::cout << " * face degree : "<< hdi.face_degree() << std::endl;
    std::cout << " * parameter   : "<< parameter  << std::endl;

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                    - fbs*msh.boundary_faces_size() ;

    auto offset_vector = full_offset(msh, hdi);

    dynamic_vector<T>  full_sol = dynamic_vector<T>::Zero(num_full_dofs);
    auto max_iter = 15;
    auto tol = 1.e-6;

    auto is_contact_vector = make_is_contact_vector(msh, bnd);


    for(size_t iter = 0; iter < max_iter; iter++)
    {
        auto cl_count = 0;
        auto assembler = make_diffusion_assembler2(msh, hdi, bnd);

        for (auto& cl : msh)
        {
            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs  = cbs + howmany_faces(msh, cl) * fbs;
            vector_type  u_full  = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hdg_scalar_stabilization(msh, cl, hdi);

            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun);//, hdi.cell_degree());
            matrix_type Ah  = gr.second + stab;

            vector_type Bnegative  = vector_type::Zero(num_total_dofs);
            matrix_type Anitsche   = matrix_type::Zero(num_total_dofs, num_total_dofs);
            matrix_type Aheaviside = matrix_type::Zero(num_total_dofs, num_total_dofs);

            if (is_contact_vector.at(cl_count) == 1)
            {
                Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                Bnegative  = make_hho_negative_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full, parameter);
                Aheaviside = make_hho_heaviside_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full, parameter);
            }

            matrix_type A =   Ah - Anitsche + Aheaviside;
            vector_type b = -(Ah - Anitsche) * u_full - Bnegative;
            b.block(0, 0, cbs, 1) += Lh;

            auto sc = diffusion_static_condensation_compute_full(msh, cl, hdi, A, b);
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

        dynamic_vector<T> dsol = dynamic_vector<T>::Zero(systsz);

        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = true;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, dsol);

        T error = 0.0 ;

        dynamic_vector<T> diff_sol = dynamic_vector<T>::Zero(num_full_dofs);

        cl_count = 0;

        for (auto& cl : msh)
        {
            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;
            vector_type  u_full = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hdg_scalar_stabilization(msh, cl, hdi);

            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun);//, hdi.cell_degree());
            matrix_type Ah  = gr.second + stab;

            vector_type Bnegative  = vector_type::Zero(num_total_dofs);
            matrix_type Anitsche   = matrix_type::Zero(num_total_dofs, num_total_dofs);
            matrix_type Aheaviside = matrix_type::Zero(num_total_dofs, num_total_dofs);

            if (is_contact_vector.at(cl_count))
            {
                Anitsche   = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                Bnegative  = make_hho_negative_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full, parameter);
                Aheaviside = make_hho_heaviside_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full, parameter);
            }

            matrix_type A =   Ah - Anitsche + Aheaviside;
            vector_type b = -(Ah - Anitsche) * u_full - Bnegative;
            b.block(0, 0, cbs, 1) += Lh;

            vector_type cell_rhs = b.block(0, 0, cbs, 1);
            vector_type du_faces = assembler.take_local_data(msh, cl, dsol);
            vector_type du_full  =
                diffusion_static_condensation_recover(msh, cl, hdi, A, cell_rhs, du_faces);

            diff_sol.block(cell_ofs, 0, num_total_dofs ,1) = du_full;

            error += du_full.dot(A * du_full);
            cl_count++;
        }

        full_sol += diff_sol;

        std::cout << "  "<< iter << "  "<< std::sqrt(error)<< std::endl;
        if( std::sqrt(error)  < tol)
        {
            std::ofstream efs("solution_whho_cpf_i" + tostr(iter) + ".dat");

            if(!efs.is_open())
                std::cout<< "Error opening file"<<std::endl;

            auto cl_count = 0;
            for(auto& cl : msh)
            {
                auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
                auto cell_ofs = offset_vector.at(cl_count++);
                vector_type u_bar = full_sol.block(cell_ofs, 0, 1, 1);
                auto bar = barycenter(msh, cl);
                efs << bar.x() << " " << bar.y() <<" "<< u_bar(0) << std::endl;
            }

            efs.close();
            return 0;
        }

    }
    return 1;
}


template<typename Mesh>
void
run_signorini(  const Mesh& msh, const algorithm_parameters<typename Mesh::scalar_type>& ap,
                const typename Mesh::scalar_type& parameter)
{
    typedef typename Mesh::point_type  point_type;
    using T =  typename Mesh::scalar_type;


    dump_to_matlab(msh,"mesh.m");

    auto force = [](const point_type& p) -> T {
        return - 2.* M_PI *  std::sin(2. * M_PI * p.x());
    };

    auto zero_fun = [](const point_type& p) -> T {
        return 0.;
    };

    typedef disk::mechanics::BoundaryConditionsScalar<Mesh> boundary_type;
    boundary_type  bnd(msh);

    /*--------------------------------------------------------------------------
    *  Check boundary labels for the unitary square domain
    *          Netgen     _____          Medit     _____
    *                4   |     | 2                |     |
    *                    |_____|                  |_____|
    *                       3                        2
    *-------------------------------------------------------------------------*/

    bnd.addDirichletBC(disk::mechanics::DIRICHLET,1, zero_fun); //TOP
    bnd.addNeumannBC(disk::mechanics::NEUMANN, 2, zero_fun); //
    bnd.addNeumannBC(disk::mechanics::NEUMANN, 4, zero_fun); //
    //bnd.addNeumannBC(disk::mechanics::NEUMANN, 3, zero_fun); //TOP
    bnd.addContactBC(disk::mechanics::SIGNORINI,3); //BOTTOM

    switch (ap.solver)
    {
        case EVAL_IN_CELLS:
            solve_cells(msh, force, ap, bnd);
            break;
        case EVAL_IN_CELLS_FULL:
            solve_cells_full(msh, force, ap, bnd);
            break;
        case EVAL_ON_FACES:
            solve_faces_borrar(msh, force, ap, bnd);
            break;
        case EVAL_IN_CELLS_AS_FACES:
            solve_cell_as_faces_test(msh, force, ap, bnd);
            break;
        case EVAL_WITH_PARAMETER:
            solve_with_parameter_full(msh, force, ap, bnd, parameter);
            break;
        default:
            throw std::invalid_argument("Invalid solver");
    }

    return;
}
