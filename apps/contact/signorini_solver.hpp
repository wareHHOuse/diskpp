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

#include <unistd.h>

#include "revolution/bases"
#include "revolution/quadratures"
#include "revolution/methods/hho"
#include "contact_hho.hpp"

#include "core/loaders/loader.hpp"

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"
#include "output/postMesh.hpp"

#include "output/silo.hpp"
#include "solvers/solver.hpp"

using namespace revolution;

template<typename T>
struct algorithm_parameters
{
    algorithm_parameters():gamma_0(0.1), theta(1.), solver("fix_point")
    {}

    T gamma_0;
    T theta;
    std::string solver;
};

template< typename T>
std::string
tostr(const T & var)
{
    std::ostringstream  ostr;
    ostr << var;
    return ostr.str();
}

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
make_is_contact_vector(const Mesh& msh,
                const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd)
{
    //cells with contact faces
    auto num_cells = msh.cells_size();
    std::vector<size_t> ret = std::vector<size_t>(num_cells);
    size_t i =0;
    for (auto& cl : msh)
    {
        auto fcs = faces(msh, cl);
        for (auto& fc:fcs)
        {
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
            const auto face_id=eid.second;

            if (bnd.is_contact_face(face_id))
            {
                ret.at(i) = 1;
                continue;
            }
        }
        i++;
    }
    return ret;
}

template<typename Mesh, typename T, typename Function>
auto
fix_point_solver_faces(const Mesh& msh, const hho_degree_info& hdi,
                const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
                const algorithm_parameters<T>& ap,
                const Function& rhs_fun)
{
    std::cout << "Solver : Fix point" << std::endl;
    std::cout << "Theta  : "<< ap.theta << std::endl;
    std::cout << "Gamma0 : "<< ap.gamma_0 << std::endl;

    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    auto is_contact_vector = make_is_contact_vector(msh, bnd);
    auto max_iter = 1000;
    auto tol = 1.e-6;

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                        - fbs*msh.boundary_faces_size() ;
    dynamic_vector<T> full_sol_old = dynamic_vector<T>::Zero(num_full_dofs);

    auto offset_vector = full_offset(msh, hdi);

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        int cl_count = 0;
        auto assembler = make_diffusion_assembler2(msh, hdi, bnd);
        for (auto& cl : msh)
        {
            auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
            auto rhs_f  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            matrix_type A = gr.second + stab;
            vector_type rhs = vector_type::Zero(num_total_dofs);
            rhs.block(0, 0, cbs, 1) = rhs_f;

            if (is_contact_vector.at(cl_count) == 1)
            {
                auto cell_ofs = offset_vector.at(cl_count);
                vector_type u_full = full_sol_old.block(cell_ofs, 0, num_total_dofs, 1);

            	A   -= make_contact_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                rhs -= make_contact_negative_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
            }

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
        std::ofstream ofs("sol_"+ tostr(iter) +".dat");

        if(!ofs.is_open())
            std::cout<< "Error opening file"<<std::endl;

        cl_count = 0;
        dynamic_vector<T> full_sol = dynamic_vector<T>::Zero(num_full_dofs);

        for (auto& cl : msh)
        {
            auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
            vector_type rhs = vector_type::Zero(num_total_dofs);
            rhs.block(0, 0, cbs, 1) = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            matrix_type A   = gr.second + stab;

            if (is_contact_vector.at(cl_count))
            {
                auto cell_ofs = offset_vector.at(cl_count);
                vector_type u_full_old = full_sol_old.block(cell_ofs, 0, num_total_dofs, 1);

            	A -= make_contact_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                rhs -=
                    make_contact_negative_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full_old);
            }
            vector_type cell_rhs = rhs.block(0, 0, cbs, 1);
            vector_type u_faces  = assembler.take_local_data(msh, cl, sol);
            vector_type u_full   =
                diffusion_static_condensation_recover(msh, cl, hdi, A, cell_rhs, u_faces);

            auto cell_ofs = offset_vector.at(cl_count);
            full_sol.block(cell_ofs, 0, num_total_dofs, 1) = u_full;

            vector_type u_full_old = full_sol_old.block(cell_ofs, 0, num_total_dofs ,1);
            vector_type diff = u_full - u_full_old;

            error += diff.dot(A * diff);

            auto bar = barycenter(msh, cl);
            ofs << bar.x() << " " << bar.y() <<" "<< u_full(0) << std::endl;
            cl_count++;
        }
        ofs.close();

        std::cout << "  "<< iter << "  "<< std::sqrt(error)<< std::endl;
        if( std::sqrt(error)  < tol)
            return full_sol;
        full_sol_old = full_sol;
    }
}

template<typename Mesh, typename T, typename Function>
auto
fix_point_solver(const Mesh& msh, const hho_degree_info& hdi,
                const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
                const algorithm_parameters<T>& ap,
                const Function& rhs_fun)
{
    if(hdi.cell_degree() != hdi.face_degree() + 1)
        throw std::invalid_argument("Change solver.");

    std::cout << "Solver : Fix point" << std::endl;
    std::cout << "Theta  : "<< ap.theta << std::endl;
    std::cout << "Gamma0 : "<< ap.gamma_0 << std::endl;

    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    auto is_contact_vector = make_is_contact_vector(msh, bnd);
    auto max_iter = 1000;
    auto tol = 1.e-6;

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                        - fbs*msh.boundary_faces_size() ;
    dynamic_vector<T> full_sol_old = dynamic_vector<T>::Zero(num_full_dofs);

    auto offset_vector = full_offset(msh, hdi);

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        int cl_count = 0;
        auto assembler = make_diffusion_assembler2(msh, hdi, bnd);
        for (auto& cl : msh)
        {
            auto cb  = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());

            if (is_contact_vector.at(cl_count) == 1)
            {
                auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
                auto cell_ofs = offset_vector.at(cl_count);
                vector_type u_full = full_sol_old.block(cell_ofs, 0, num_total_dofs, 1);

                vector_type rhs_force = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
                vector_type rhs      = vector_type::Zero(num_total_dofs);
                rhs.block(0, 0, cbs, 1) = rhs_force;

                auto gr   = make_contact_hho_scalar_laplacian(msh, cl, hdi, bnd);
                auto stab = make_hho_scalar_stabilization_simple(msh, cl, hdi);
                matrix_type A  = gr.second + stab;

            	A -= make_contact_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                rhs -= make_contact_negative(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);

                auto sc = diffusion_static_condensation_compute_full(msh, cl, hdi, A, rhs);
                assembler.assemble(msh, cl, sc.first, sc.second);
            }
            else
            {
                auto gr   = make_hho_scalar_laplacian(msh, cl, hdi);
                auto stab = make_hho_scalar_stabilization_simple(msh, cl, hdi);
                vector_type rhs = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
                matrix_type A   = gr.second + stab;
                auto sc = diffusion_static_condensation_compute(msh, cl, hdi, A, rhs);
                assembler.assemble(msh, cl, sc.first, sc.second);
            }
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
        std::ofstream ofs("sol_"+ tostr(iter) +".dat");

        if(!ofs.is_open())
            std::cout<< "Error opening file"<<std::endl;

        cl_count = 0;
        dynamic_vector<T> full_sol = dynamic_vector<T>::Zero(num_full_dofs);

        for (auto& cl : msh)
        {
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;

            auto cell_ofs = offset_vector.at(cl_count);
            vector_type u_full_old = full_sol_old.block(cell_ofs, 0, num_total_dofs, 1);
            vector_type u_full = vector_type::Zero(num_total_dofs, 1);
            matrix_type A = matrix_type::Zero(num_total_dofs, num_total_dofs);

            if (is_contact_vector.at(cl_count))
            {
                vector_type rhs_force = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
                vector_type rhs      = vector_type::Zero(num_total_dofs);
                rhs.block(0, 0, cbs, 1) = rhs_force;

                auto gr   = make_contact_hho_scalar_laplacian(msh, cl, hdi, bnd);
                auto stab = make_hho_scalar_stabilization_simple(msh, cl, hdi);
                A  = gr.second + stab;
                A -= make_contact_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                rhs -= make_contact_negative(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full_old);

                vector_type cell_rhs = rhs.block(0, 0, cbs, 1);
                vector_type u_faces  = assembler.take_local_data(msh, cl, sol);
                u_full = diffusion_static_condensation_recover(msh, cl, hdi, A, cell_rhs, u_faces);
            }
            else
            {
                auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
                auto stab   = make_hho_scalar_stabilization_simple(msh, cl, hdi);
                A   = gr.second + stab;
                vector_type rhs = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
                vector_type u_faces  = assembler.take_local_data(msh, cl, sol);
                u_full = diffusion_static_condensation_recover(msh, cl, hdi, A, rhs, u_faces);
            }

            vector_type diff = u_full - u_full_old;
            error += diff.dot(A * diff);

            full_sol.block(cell_ofs, 0, num_total_dofs, 1) = u_full;

            auto bar = barycenter(msh, cl);
            ofs << bar.x() << " " << bar.y() <<" "<< u_full(0) << std::endl;
            cl_count++;
        }
        ofs.close();

        std::cout << "  "<< iter << "  "<< std::sqrt(error)<< std::endl;
        if( std::sqrt(error)  < tol)
            return full_sol;
        full_sol_old = full_sol;
    }
}

template<typename Mesh, typename T, typename Function>
auto
newton_solver(const Mesh& msh, const hho_degree_info& hdi,
                const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
                const algorithm_parameters<T>& ap,
                const Function& rhs_fun,
                dynamic_vector<T>& full_sol)
{
    std::cout << "Solver : Generalized Newton" << std::endl;
    std::cout << "Theta  : "<< ap.theta << std::endl;
    std::cout << "Gamma0 : "<< ap.gamma_0 << std::endl;

    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    auto max_iter = 1000;
    auto tol = 1.e-6;
    auto is_contact_vector = make_is_contact_vector(msh, bnd);

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                        - fbs*msh.boundary_faces_size() ;

    auto offset_vector = full_offset(msh, hdi);

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        auto cl_count = 0;

        auto assembler = make_diffusion_assembler2(msh, hdi, bnd);
        for (auto& cl : msh)
        {
            auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
            auto rhs_force  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            matrix_type A   = gr.second + stab;
            vector_type rhs = vector_type::Zero(num_total_dofs);
            rhs.block(0, 0, cbs, 1) = rhs_force;

            auto cell_ofs = offset_vector.at(cl_count);
            vector_type u_full_old = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            if (is_contact_vector.at(cl_count))
            {
                //matrix An(du, v)
            	A -= make_contact_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );

                //rhs: int [Pn(u_old,v)]_Pn[v]
                rhs -= make_contact_negative(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full_old);
            }
            //rhs: An(u_old,v)
            rhs -= A * u_full_old;

            if (is_contact_vector.at(cl_count))
               A -=  make_contact_heaviside_other(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full_old);

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

        dynamic_vector<T> dsol = dynamic_vector<T>::Zero(systsz);

        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = true;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, dsol);

        T error = 0.0 ;
        std::ofstream ofs("sol_"+ tostr(iter) +".dat");

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

            vector_type rhs = vector_type::Zero(num_total_dofs);
            rhs.block(0, 0, cbs, 1) = rhs_force;

            auto cell_ofs = offset_vector.at(cl_count);
            vector_type u_full_old = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            if (is_contact_vector.at(cl_count))
            {
                //matrix An(du, v)
            	A -= make_contact_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );

                //rhs: int [Pn(u_old,v)]_Pn[v]
                rhs -= make_contact_negative(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full_old);
            }
            //rhs: An(u_old,v)
            rhs -= A * u_full_old;

            //matrix Heaviside term: Dont use this before rhs += A * u_full;
            if (is_contact_vector.at(cl_count))
               A -=  make_contact_heaviside_other(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full_old);

            vector_type cell_rhs = rhs.block(0, 0, cbs, 1);
            vector_type du_faces = assembler.take_local_data(msh, cl, dsol);
            vector_type du_full  =
                diffusion_static_condensation_recover(msh, cl, hdi, A, cell_rhs, du_faces);

            vector_type  diff = du_full;
            error += diff.dot(A * diff);

            full_sol.block(cell_ofs, 0, num_total_dofs ,1) += du_full;

            auto bar = barycenter(msh, cl);
            ofs << bar.x() << " " << bar.y() <<" "<< du_full(0) << std::endl;
            cl_count++;
        }
        ofs.close();

        std::cout << "  "<< iter << "  "<< std::sqrt(error)<< std::endl;
        if( std::sqrt(error)  < tol)
            return 0;
    }
}
