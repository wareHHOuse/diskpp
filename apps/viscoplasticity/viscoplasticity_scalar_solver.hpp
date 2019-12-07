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

#include "bases/bases.hpp"
#include "quadratures/quadratures.hpp"
#include "methods/hho"

#include "core/loaders/loader.hpp"

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"
#include "output/postMesh.hpp"

#include "output/silo.hpp"
#include "solvers/solver.hpp"
#include "viscoplasticity_utils.hpp"

using namespace disk;


template<typename Mesh>
auto
make_bnd(const Mesh& msh, const scalar_problem_type& problem)
{
    using T = typename Mesh::coordinate_type;


    auto zero_fun  = [](const typename Mesh::point_type& p) -> T {
        return 0.;
    };

    disk::scalar_boundary_conditions<Mesh> bnd(msh);

    switch(problem)
    {
        case ANNULUS:
            //bnd.addDirichletEverywhere(zero_fun); //TOP
            //#if 0
            bnd.addDirichletBC(disk::DIRICHLET,1, zero_fun); //TOP
            bnd.addNeumannBC(disk::NEUMANN, 2, zero_fun); //
            //#endif
            break;
        case CIRCULAR:
            bnd.addDirichletBC(disk::DIRICHLET,2, zero_fun); //TOP
            bnd.addNeumannBC(disk::NEUMANN, 1, zero_fun); //
            break;
        case SQUARE:
            bnd.addDirichletEverywhere(zero_fun); //TOP
            break;
        default:
            throw std::invalid_argument("No  defined problem");
    }
    return bnd;
}
template<typename Mesh>
auto
make_rhs_functor(const Mesh& msh, const scalar_problem_type& problem)
{
    using T = typename Mesh::coordinate_type;
    auto rhs_fun  = [](const typename Mesh::point_type& p) -> T {
        return 1.;
    };
    return rhs_fun;
}


template<typename Mesh>
class augmented_lagrangian_viscoplasticity
{
    typedef Mesh mesh_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;
    typedef typename mesh_type::coordinate_type T;
    typedef typename mesh_type::point_type  point_type;
    typedef disk::scalar_boundary_conditions<mesh_type> boundary_type;

    typedef Matrix<T, Mesh::dimension, 1>      tensor_type;
    typedef Matrix<T, Dynamic, Dynamic>        matrix_type;
    typedef Matrix<T, Dynamic, 1>              vector_type;

    typedef std::function<T   (const point_type &)>             scalar_funtion_type;

    hho_degree_info di;
    T             yield;
    size_t        cbs, fbs, sbs;

    tensors_at_quad_pts_utils<mesh_type>    tsr_utils;
    std::vector<std::pair<size_t, size_t>>  tsr_offsets_vector;
    std::vector<size_t>                     sol_offset_map;
    viscoplasticity_data<T>                 vp;

    disk::dynamic_matrix<T>     multiplier, auxiliar, auxiliar_old;

public:
    disk::dynamic_vector<T>   full_sol, full_sol_old;
    std::pair<T, T>     convergence;
    //boundary_type           bnd;

    augmented_lagrangian_viscoplasticity(const Mesh& msh,
                            const hho_degree_info & hdi,
                            const viscoplasticity_data<T>& vpst_parameters):
                            di(hdi), vp(vpst_parameters)
    {
        yield   =  vp.Bn * vp.f * vp.Lref / 2.; // * vp.mu;// * omegaExt; //* f * Lref;

        const auto dim =  Mesh::dimension;

        cbs = scalar_basis_size(di.cell_degree(), dim);
        fbs = scalar_basis_size(di.face_degree(), dim - 1);
        sbs = vector_basis_size(di.face_degree(), dim, dim);

        size_t quad_degree = 2. * di.face_degree();
        tsr_utils = tensors_at_quad_pts_utils<mesh_type>(msh, quad_degree);
        tsr_offsets_vector = tsr_utils.offsets_vector();
        auto pair_offset_map = make_scalar_solution_offset(msh, di);
        auto total_dofs = pair_offset_map.first;
        sol_offset_map = pair_offset_map.second;

        full_sol = disk::dynamic_vector<T>::Zero(total_dofs);
        full_sol_old = disk::dynamic_vector<T>::Zero(total_dofs);
    };

    matrix_type
    compute_auxiliar(   const mesh_type& msh,
                        const cell_type& cl,
                        const disk::dynamic_vector<T>& full_solution)
    {
        auto cl_id = msh.lookup( cl);
        auto num_total_dofs = cbs + fbs * howmany_faces(msh, cl);
        auto sol_offset     = sol_offset_map.at(cl_id);
        vector_type u_TF    = full_solution.block(sol_offset, 0, num_total_dofs, 1);

        auto value = 1./(vp.mu + vp.alpha);
        auto G = make_hho_gradrec_vector(msh, cl, di);
        vector_type   Gu = G.first * u_TF;

        auto qps = integrate(msh, cl, tsr_utils.quad_degree());

        auto offset = tsr_offsets_vector.at(cl_id).first;
        assert(tsr_offsets_vector.at(cl_id).second == qps.size()); //Take out this after

        auto sb = make_vector_monomial_basis(msh, cl, di.face_degree());

        matrix_type stress = multiplier.block(0, offset, sbs, qps.size());
        matrix_type gamma  = matrix_type::Zero(sbs, qps.size());

        size_t qp_count    = 0;
        for(auto& qp: qps)
        {
            auto s_phi  = sb.eval_functions(qp.point());
            vector_type stress_qp = stress.block( 0, qp_count, sbs, 1);
            vector_type theta     = stress_qp  +  vp.alpha * Gu;

            tensor_type theta_eval = eval(theta, s_phi);

            T theta_norm  = theta_eval.norm();
            T tol = 1.e-8;

            // A. Liquid
            if( (theta_norm + tol) >   yield)
                gamma.block( 0, qp_count, sbs, 1) = value * theta * (1. - (yield/theta_norm));

            qp_count++;
        }
        assert(qps.size() == qp_count);

        return gamma;
    }

    auto
    update_multiplier(const mesh_type& msh)
    {
        T conv_stress = 0.;
        T conv_gamma = 0.;

        const auto dim = Mesh::dimension;

        size_t cl_id = 0;
        for(auto cl: msh)
        {
            auto sb = make_vector_monomial_basis(msh, cl, di.face_degree());

            auto num_total_dofs = cbs + fbs * howmany_faces(msh, cl);
            auto sol_offset     = sol_offset_map.at(cl_id++);
            vector_type u_TF    = full_sol.block(sol_offset, 0, num_total_dofs, 1);
            auto G = make_hho_gradrec_vector(msh, cl, di);
            vector_type Gu = G.first * u_TF;

            auto cl_id = msh.lookup(cl);
            auto offset =  tsr_offsets_vector.at(cl_id).first;
            auto qps = integrate(msh, cl, tsr_utils.quad_degree());

            matrix_type     gamma_now = compute_auxiliar(msh, cl, full_sol_old); //gamm^n computed with (u,sigma)^n-1
            matrix_type     gamma_old = auxiliar.block(0, offset, sbs, qps.size());

            //gamma updating for next iteration
            auxiliar.block(0, offset, sbs, qps.size()) = gamma_now;

            matrix_type     diff_gamma  =   vp.alpha * (gamma_now - gamma_old);
            matrix_type     diff_stress = - vp.alpha *  gamma_now;

            //Update multiplier
            for(size_t i = 0; i < qps.size(); i++)
                diff_stress.block(0, i, sbs, 1) += vp.alpha *  Gu;
            multiplier.block(0, offset, sbs, qps.size()) += diff_stress;

            //Error computations
            size_t qp_count = 0;
            for(auto& qp: qps)
            {
                auto s_phi  = sb.eval_functions(qp.point());
                matrix_type mm =  priv::outer_product(s_phi, s_phi);

                vector_type diff_stress_qp  = diff_stress.block(0, qp_count, sbs, 1);
                vector_type diff_gamma_qp   = diff_gamma.block(0, qp_count, sbs, 1);

                conv_stress += diff_stress_qp.dot(mm * diff_stress_qp);
                conv_gamma  += diff_gamma_qp.dot(mm * diff_gamma_qp);

                qp_count++;
            }
        }
        convergence = std::make_pair(conv_stress, conv_gamma);

        auxiliar_old = auxiliar;
        return;
    }

    template<typename Assembler>
    Matrix<T, Dynamic, 1>
    make_rhs_alg(   const mesh_type& msh,
                    const cell_type& cl,
                    const Assembler& assembler)
    {
        auto G  = make_hho_gradrec_vector(msh, cl, di);
        auto cb = make_scalar_monomial_basis(msh, cl, di.cell_degree());
        auto sb = make_vector_monomial_basis(msh, cl, di.face_degree());

        auto cell_ofs =  priv::offset(msh, cl);
        auto num_faces = howmany_faces(msh, cl);

        //(stress - alpha * gamma, Gv)
        auto cl_id  =  msh.lookup(cl);
        auto offset =  tsr_offsets_vector.at(cl_id).first;
        auto qps    =  integrate(msh, cl, tsr_utils.quad_degree());
        auto tsr_size = tsr_offsets_vector.at(cl_id).second;
        assert( tsr_size == qps.size()); //Take out this after

        matrix_type stress = multiplier.block(0, offset, sbs, qps.size());
        matrix_type gamma  = compute_auxiliar(msh, cl, full_sol_old);

        matrix_type str_agam = stress -  vp.alpha * gamma;
        vector_type rhs = vector_type::Zero(cbs + num_faces * fbs);

        size_t qp_count = 0;

        //#if 0
        for(auto& qp: qps)
        {
            auto s_phi  = sb.eval_functions(qp.point());
            matrix_type mm =  priv::outer_product(s_phi, s_phi);

            vector_type str_agam_qp  = str_agam.block(0, qp_count, sbs, 1);

            rhs -= qp.weight() * G.first.transpose() * mm * str_agam_qp;
            qp_count++;
        }
        //#endif

        //(f, v_T)
        auto rhs_fun = make_rhs_functor(msh, vp.problem);
        rhs.block( 0, 0, cbs, 1) += make_rhs(msh, cl, cb, rhs_fun);

        return rhs;
    }

    template<typename Assembler>
    void
    make_global_rhs(const mesh_type& msh, Assembler& assembler)
    {
        //std::cout << "START GLOBAL RHS" << std::endl;

        assembler.initialize_rhs();

        for (auto cl : msh)
        {
            auto G  = make_hho_gradrec_vector(msh, cl, di);
            auto gr = make_hho_scalar_laplacian(msh, cl, di);
            matrix_type stab = make_hho_scalar_stabilization(msh, cl, gr.first, di);
            matrix_type A   = vp.alpha * G.second + vp.mu * stab;
            vector_type rhs = make_rhs_alg(msh, cl, assembler);

            auto sc  = diffusion_static_condensation_compute_full(msh, cl, di, A, rhs);
            assembler.assemble_rhs(msh, cl, sc.second);
        }
        assembler.finalize_rhs();

        //std::cout << "END GLOBAL RHS" << std::endl;

        return;
    }

    template<typename Assembler>
    void
    make_global_matrix(const mesh_type& msh, Assembler& assembler)
    {
        //std::cout << "START GLOBAL MATTRIX" << std::endl;

        assembler.initialize_lhs();

        for (auto cl : msh)
        {
            auto G  = make_hho_gradrec_vector(msh, cl, di);
            auto gr = make_hho_scalar_laplacian(msh, cl, di);
            matrix_type stab   = make_hho_scalar_stabilization(msh, cl, gr.first, di);
            matrix_type A = vp.alpha * G.second + vp.mu * stab;

            vector_type rhs_zero = vector_type::Zero(A.rows());
            auto sc  = diffusion_static_condensation_compute(msh, cl, di, A, rhs_zero);
            assembler.assemble_lhs(msh, cl, sc.first);
        }

        assembler.finalize_lhs();

        //std::cout << "END GLOBAL MATTRIX" << std::endl;

        return;
    }

    template<typename Assembler>
    auto
    recover_solution(const mesh_type& msh, const Assembler& assembler, const disk::dynamic_vector<T>& sol)
    {
        //std::cout << "START RECOVER SOLUTION" << std::endl;

        size_t cl_id  = 0;
        for (auto& cl : msh)
        {
            auto G  = make_hho_gradrec_vector(msh, cl, di);
            auto gr = make_hho_scalar_laplacian(msh, cl, di);
            matrix_type stab = make_hho_scalar_stabilization(msh, cl, gr.first, di);
            matrix_type A   = vp.alpha * G.second + vp.mu * stab;
            vector_type rhs = make_rhs_alg(msh, cl, assembler);

            vector_type cell_rhs = rhs.block( 0, 0, cbs, 1);
            vector_type uloc   = assembler.take_local_data(msh, cl, sol);
            vector_type ufull  =
                diffusion_static_condensation_recover(msh, cl, di,  A, cell_rhs, uloc);

            auto sol_offset     = sol_offset_map.at(cl_id++);
            auto num_total_dofs = cbs + fbs * howmany_faces(msh, cl);
            full_sol.block(sol_offset, 0, num_total_dofs, 1) = ufull;
        }
        //std::cout << "END RECOVER SOLUTION" << std::endl;

        return;
    }

    template<typename Assembler>
    void
    post_processing(const mesh_type& msh, const Assembler& assembler,
                    const std::string & info)
    {
        //std::cout << "START POST-PROCESSING" << std::endl;
        auto dim = Mesh::dimension;
        auto rbs = scalar_basis_size(di.reconstruction_degree(), dim);

        disk::dynamic_vector<T> cell_sol(cbs * msh.cells_size());
        disk::dynamic_vector<T> cell_rec_sol(rbs * msh.cells_size());

        std::ofstream ofs("data_" + info + ".data");
        if (!ofs.is_open())
            std::cout << "Error opening file"<<std::endl;

        size_t cl_id = 0;
        for(auto cl : msh)
        {
            auto gr  = make_hho_scalar_laplacian(msh, cl, di);
            auto num_total_dofs = cbs + fbs * howmany_faces(msh, cl);

            auto sol_offset    = sol_offset_map.at(cl_id++);
            vector_type ufull  = full_sol.block(sol_offset, 0, num_total_dofs, 1);
            assert((gr.first * ufull).rows() == rbs - dim);

            auto cell_ofs =  priv::offset(msh, cl);
            cell_rec_sol.block(cell_ofs * rbs + dim, 0, rbs - dim, 1) = gr.first * ufull;
            cell_rec_sol.block(cell_ofs * rbs, 0, dim, 1) = ufull.block(0,0, dim, 1);
            cell_sol.block(cell_ofs * cbs, 0, cbs, 1)     = ufull.block(0,0, cbs, 1);

            //plot using only the value in the barycenter
            auto bar = barycenter(msh, cl);

            //Velocity
            vector_type cell_vel = ufull.block(0,0, cbs, 1);
            auto cb  = make_scalar_monomial_basis(msh, cl, di.cell_degree());
            auto phi = cb.eval_functions(bar);
            T ueval = eval(cell_vel, phi);

            ofs << bar.x()<< " "<< bar.y()<< " " << ueval << std::endl;
        }
        ofs.close();

        //std::cout << "END POST-PROCESSING" << std::endl;
        return;
    }

    bool
    run_alg(const mesh_type& msh)
    {
        std::string info = "_k" + tostr(di.face_degree()) + "_a" + tostr(vp.alpha);

        boundary_type bnd = make_bnd(msh, vp.problem);

        auto assembler = make_diffusion_assembler_alg(msh, di, bnd);

        auto systsz = assembler.global_system_size();
        disk::dynamic_vector<T> sol =  disk::dynamic_vector<T>::Zero(systsz);

        auto num_total_quads = tsr_utils.num_total_quad_points();

        multiplier   = disk::dynamic_matrix<T>::Zero(sbs, num_total_quads);
        auxiliar     = disk::dynamic_matrix<T>::Zero(sbs, num_total_quads);
        auxiliar_old = disk::dynamic_matrix<T>::Zero(sbs, num_total_quads);

        auto max_iters = 50000;
        auto Ninf = 1.e+4;
        auto tolerance = 1.e-8;

        for(size_t iter = 0; iter < max_iters; iter++)
        {
            //WARNINGS: This one must go after make_global_matrix!!!!!!
            //WARNINGS: This one must go after make_global_matrix!!!!!!
            //WARNINGS: This one must go after make_global_matrix!!!!!!
            //WARNINGS: This one must go after make_global_matrix!!!!!!
            make_global_rhs(msh, assembler);
            if(iter == 0)
                make_global_matrix(msh, assembler);

            size_t systsz = assembler.LHS.rows();
            size_t nnz = assembler.LHS.nonZeros();

            disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(systsz);
            disk::solvers::pardiso_params<T> pparams;
            pparams.report_factorization_Mflops = true;

            mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

            recover_solution(msh, assembler, sol);
            update_multiplier(msh);

            //---------------------------------------------------------------------
            T cvg_total = std::sqrt(convergence.first + convergence.second);

            if(iter % 100 == 0)
                std::cout << "  i : "<< iter <<"  - " << std::sqrt(cvg_total)<<std::endl;

            assert(cvg_total < Ninf);
            if( cvg_total < tolerance)
            {
                std::cout << "  i : "<< iter <<"  - " << cvg_total <<std::endl;
                post_processing( msh, assembler, info);
                return true;
            }
            //sol_old = sol;

        }
        return false;
    }
};
