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

#include "methods/hho"

#include "core/loaders/loader.hpp"

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"
#include "output/postMesh.hpp"

#include "output/silo.hpp"
#include "solvers/solver.hpp"
#include "viscoplasticity_utils.hpp"

using namespace disk;

enum problem_type
{
    DRIVEN
};

template<typename Mesh>
class augmented_lagrangian_viscoplasticity
{
    typedef Mesh mesh_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;
    typedef typename mesh_type::coordinate_type T;

    typedef disk::mechanics::BoundaryConditions<mesh_type> boundary_type;

    using point_type = typename mesh_type::point_type;

    typedef Matrix<T, Mesh::dimension, Mesh::dimension>         tensor_type;
    typedef Matrix<T, Mesh::dimension, 1>                       vector2d_type;

    typedef Matrix<T, Dynamic, Dynamic>                 matrix_type;
    typedef Matrix<T, Dynamic, 1>                       vector_type;

    typedef std::function<vector2d_type (const point_type &)>   vector_funtion_type;
    typedef std::function<T   (const point_type &)>             scalar_funtion_type;

    hho_degree_info di;
    T             viscosity;
    T             alpha;
    T             yield;
    size_t        cbs, fbs, pbs, sbs;
    vector_funtion_type                     rhs_fun;
    tensors_at_quad_pts_utils<mesh_type>    tsr_utils;
    std::vector<std::pair<size_t, size_t>>  tsr_offsets_vector;

    dynamic_matrix<T>     multiplier, auxiliar, auxiliar_old;

public:
    dynamic_vector<T>   sol_old;
    std::pair<T, T>     convergence;
    //boundary_type           bnd;

    augmented_lagrangian_viscoplasticity(const Mesh& msh,
                            const hho_degree_info & hdi,
                            const T& alpha_ext):
                            di(hdi), alpha(alpha_ext)
    {
        viscosity = 1.;
        T f     = 1;
        T Lref  = 1.;
        T Bn    =  2* std::sqrt(2);
        yield   =  Bn * f * Lref; // * viscosity;// * omegaExt; //* f * Lref;

        const auto dim =  Mesh::dimension;

        cbs = vector_basis_size(di.cell_degree(), dim, dim);
        fbs = vector_basis_size(di.face_degree(), dim - 1, dim);
        pbs = scalar_basis_size(di.face_degree(), dim);
        sbs = sym_matrix_basis_size(di.face_degree(), dim, dim);

        size_t quad_degree = 2. * di.face_degree();
        tsr_utils = tensors_at_quad_pts_utils<mesh_type>(msh, quad_degree);
        tsr_offsets_vector = tsr_utils.offsets_vector();
    };

    #if 0
    auto
    define_problem(const mesh_type& msh, const problem_type& problem )
    {   }
    #endif

    template<typename Assembler>
    matrix_type
    compute_auxiliar(   const mesh_type& msh,
                        const cell_type& cl,
                        const Assembler& assembler,
                        const vector_type& velocity_dofs)
    {
        vector_type u_TF  = assembler.take_velocity(msh, cl, velocity_dofs);
        auto value = 1./(2. * (viscosity + alpha));
        auto G = make_hlow_stokes(msh, cl, di, true);
        vector_type   Gu = G.first * u_TF;

        auto qps = integrate(msh, cl, tsr_utils.quad_degree());

        auto cl_id = msh.lookup(cl);
        auto offset = tsr_offsets_vector.at(cl_id).first;
        assert( tsr_offsets_vector.at(cl_id).second == qps.size()); //Take out this after

        auto sb = make_sym_matrix_monomial_basis(msh, cl, di.face_degree());

        matrix_type stress = multiplier.block(0, offset, sbs, qps.size());
        matrix_type gamma  = matrix_type::Zero(sbs, qps.size());
        size_t qp_count    = 0;

        for(auto& qp: qps)
        {
            auto s_phi  = sb.eval_functions(qp.point());
            vector_type stress_qp = stress.block( 0, qp_count, sbs, 1);
            vector_type theta     = stress_qp  +  2. * alpha * Gu;

            tensor_type theta_eval = eval(theta, s_phi);

            T theta_norm  = theta_eval.norm();
            T tol = 1.e-8;

            // A. Liquid
            if( (theta_norm + tol) >  std::sqrt(2) * yield)
                gamma.block( 0, qp_count, sbs, 1) = value * theta *
                                            (1. - std::sqrt(2)*(yield/theta_norm));

            qp_count++;
        }
        assert(qps.size() == qp_count);

        return gamma;

    }

    template<typename Assembler>
    auto
    update_multiplier(const mesh_type& msh, const Assembler& assembler, const dynamic_vector<T>& sol)
    {
        T conv_stress = 0.;
        T conv_gamma = 0.;

        const auto dim = Mesh::dimension;

        for(auto cl: msh)
        {
            auto sb = make_sym_matrix_monomial_basis(msh, cl, di.face_degree());

            vector_type u_TF = assembler.take_velocity(msh, cl, sol);
            auto G = make_hlow_stokes(msh, cl, di, true);
            vector_type Gu = G.first * u_TF;

            auto cl_id = msh.lookup(cl);
            auto offset =  tsr_offsets_vector.at(cl_id).first;

            auto qps = integrate(msh, cl, tsr_utils.quad_degree());

            matrix_type     gamma_old = auxiliar_old.block(0, offset, sbs, qps.size());
            matrix_type     gamma     = auxiliar.block(0, offset, sbs, qps.size());
            matrix_type     diff_gamma  = alpha * (gamma - gamma_old);
            matrix_type     diff_stress = - 2. * alpha *  gamma;

            //Update multiplier
            for(size_t i = 0; i < qps.size(); i++)
                diff_stress.block(0, i++, sbs, 1) += 2. * alpha *  Gu;
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
        auto G = make_hlow_stokes(msh, cl, di, true);
        auto cb = make_vector_monomial_basis(msh, cl, di.cell_degree());
        auto sb = make_sym_matrix_monomial_basis(msh, cl, di.face_degree());

        auto cell_ofs =  priv::offset(msh, cl);
        auto num_faces = howmany_faces(msh, cl);

        //(stress - alpha * gamma, Gv)
        auto cl_id = msh.lookup(cl);
        auto offset =  tsr_offsets_vector.at(cl_id).first;
        auto qps = integrate(msh, cl, tsr_utils.quad_degree());
        auto tsr_size = tsr_offsets_vector.at(cl_id).second;
        assert( tsr_size == qps.size()); //Take out this after

        matrix_type stress = multiplier.block(0, offset, sbs, qps.size());
        matrix_type gamma  = compute_auxiliar( msh, cl, assembler, sol_old);

        auxiliar.block(0, offset, sbs, qps.size()) = gamma;

        matrix_type str_agam = stress - 2. * alpha * gamma;

        vector_type rhs = vector_type::Zero(cbs + num_faces * fbs);
        size_t qp_count = 0;

        for(auto& qp: qps)
        {
            auto s_phi  = sb.eval_functions(qp.point());
            matrix_type mm =  priv::outer_product(s_phi, s_phi);

            vector_type str_agam_qp  = str_agam.block(0, qp_count, sbs, 1);

            rhs -= qp.weight() * G.first.transpose() * mm * str_agam_qp;
            qp_count++;
        }


        //(f, v_T)
        rhs.block( 0, 0, cbs, 1) += make_rhs(msh, cl, cb, rhs_fun);

        return rhs;
    }


    template<typename Assembler>
    void
    make_global_rhs(const mesh_type& msh, Assembler& assembler)
    {
        assembler.initialize_rhs();

        for (auto cl : msh)
        {
            vector_type local_rhs = make_rhs_alg(msh, cl, assembler);
            assembler.assemble_rhs(msh, cl, local_rhs);
        }
        assembler.finalize_rhs();

        return;
    }

    template<typename Assembler>
    void
    make_global_matrix(const mesh_type& msh, Assembler& assembler)
    {
        assembler.initialize_lhs();

        for (auto cl : msh)
        {
            auto G  = make_hlow_stokes(msh, cl, di, true);
            auto gr = make_hho_stokes(msh, cl, di, true);
            matrix_type stab = make_hho_vector_stabilization(msh, cl, gr.first, di);
            auto dr = make_hho_divergence_reconstruction_stokes_rhs(msh, cl, di);

            matrix_type A = 2.*(alpha * G.second + viscosity * stab);

            assembler.assemble_lhs(msh, cl, A, -dr);
        }

        assembler.finalize_lhs();

        return;
    }

    template<typename Assembler>
    void
    post_processing(const mesh_type& msh, const Assembler& assembler,
                    const std::string & info,
                    const problem_type& problem)
    {
        auto dim = Mesh::dimension;
        auto rbs = vector_basis_size(di.reconstruction_degree(), dim, dim);

        dynamic_vector<T> cell_sol(cbs * msh.cells_size());
        dynamic_vector<T> cell_rec_sol(rbs * msh.cells_size());
        dynamic_vector<T> press_vec(pbs * msh.cells_size());

        std::ofstream ofs("data_" + info + ".data");
        if (!ofs.is_open())
            std::cout << "Error opening file"<<std::endl;

        for(auto cl : msh)
        {
            auto gr  = make_hho_stokes(msh, cl, di, true);
            auto cell_ofs = priv::offset(msh, cl);
            vector_type svel =  assembler.take_velocity(msh, cl, sol_old);
            assert((gr.first * svel).rows() == rbs - dim);
            cell_rec_sol.block(cell_ofs * rbs + dim, 0, rbs - dim, 1) = gr.first * svel;
            cell_rec_sol.block(cell_ofs * rbs, 0, dim, 1) = svel.block(0,0, dim, 1);
            cell_sol.block(cell_ofs * cbs, 0, cbs, 1)     = svel.block(0,0, cbs, 1);

            //this is only for k = 0, since there is only one velocity;
            auto bar = barycenter(msh, cl);

            //Velocity
            vector_type cell_vel = svel.block(0,0, cbs, 1);
            auto cb  = make_vector_monomial_basis(msh, cl, di.cell_degree());
            auto phi = cb.eval_functions(bar);
            vector_type ueval = eval(cell_vel, phi);

            //Pressure
            vector_type spress =  assembler.take_pressure(msh, cl, sol_old);
            auto pb  = make_scalar_monomial_basis(msh, cl, di.face_degree());
            auto p_phi = pb.eval_functions(bar);
            T peval =  p_phi.dot(spress);

            /***********************************************************************/
            #if 0
            //Stress
            auto sb = make_sym_matrix_monomial_basis(msh, cl, di.face_degree());
            auto s_phi  = sb.eval_functions(bar);
            auto G = make_hlow_stokes(msh, cl, di, use_sym_grad);
            vector_type Gu = G.first * svel;

            //*************************************************************
            vector_type stress = multiplier.block(cell_ofs * sbs, 0, sbs, 1);
            vector_type theta  = stress  +  factor * alpha * Gu;
            tensor_type theta_eval = eval(theta, s_phi);
            tensor_type sigma_eval = eval(stress, s_phi);
            //*************************************************************

            //div
            tensor_type grad_eval = eval(Gu, s_phi);
            T divu = grad_eval(0,0) + grad_eval(1,1);

            // tr_stress
            T tr_stress = sigma_eval(0,0) + sigma_eval(1,1);
            #endif

            ofs << ueval(0)   << " " << ueval(1) << " " << peval<< " ";
            //ofs << theta_eval.norm() << " " << sigma_eval.norm()   << " ";
            //ofs << divu << " "<< tr_stress<<std::endl;
            ofs <<std::endl;
        }
        ofs.close();

        std::pair<point_type, point_type> p_x, p_y;
        auto eps = 1.e-4;
        p_y = std::make_pair(point_type({0.5 + eps, 0.0 + eps}), point_type({0.5 + eps, 1.0 + eps}));

        //plot_over_line(msh, p_x, cell_rec_sol, di.reconstruction_degree(), "plot_over_x_" + info + ".data");
        //plot_over_line(msh, p_y, cell_rec_sol, di.reconstruction_degree(), "plot_over_y_" + info + ".data");
        plot_over_line(msh, p_y, cell_sol, di.cell_degree(), "plot_over_y_" + info + ".data");
        compute_discontinuous_velocity( msh, cell_sol, di, "velocity_" + info +".msh");
        save_coords(msh, "Coords_"+ info + ".data");
        quiver( msh, sol_old, assembler, di, "quiver_"+ info + ".data");
        return;
    }

    bool
    run_alg(const mesh_type& msh, const std::string& info, const problem_type& problem )
    {
        boundary_type bnd(msh);

        auto wall = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
            return Matrix<T, Mesh::dimension, 1>::Zero();
        };
        auto movingWall = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
            return Matrix<T, Mesh::dimension, 1>{1,0};
        };

        switch (problem)
        {
            case DRIVEN:

                rhs_fun  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    return Matrix<T, Mesh::dimension, 1>::Zero();
                };

                bnd.addDirichletBC(0, 1, movingWall);
                bnd.addDirichletBC(0, 2, wall);
                bnd.addDirichletBC(0, 3, wall);
                bnd.addDirichletBC(0, 4, wall);
               break;
            default:
                throw std::invalid_argument("Invalid problem");
        }

        auto assembler = make_stokes_assembler_alg(msh, di, bnd);
        auto systsz = assembler.global_system_size();
        sol_old = vector_type::Zero(systsz);
        dynamic_vector<T> sol =  dynamic_vector<T>::Zero(systsz);

        auto num_total_quads = tsr_utils.num_total_quad_points();

        multiplier   = dynamic_matrix<T>::Zero(sbs, num_total_quads);
        auxiliar     = dynamic_matrix<T>::Zero(sbs, num_total_quads);
        auxiliar_old = dynamic_matrix<T>::Zero(sbs, num_total_quads);

        auto max_iters = 50000;
        auto Ninf = 1.e+4;
        auto tolerance = 1.e-8;

        for(size_t iter = 0; iter < max_iters; iter++)
        {
            make_global_rhs(msh,assembler);
            if(iter == 0)
                make_global_matrix(msh, assembler);

            //dump_sparse_matrix(assembler.LHS, "stokes.txt");
            size_t systsz = assembler.LHS.rows();
            size_t nnz = assembler.LHS.nonZeros();

            dynamic_vector<T> sol = dynamic_vector<T>::Zero(systsz);
            disk::solvers::pardiso_params<T> pparams;
            pparams.report_factorization_Mflops = true;
            mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

            update_multiplier(msh, assembler, sol);
            //---------------------------------------------------------------------
            T cvg_total = std::sqrt(convergence.first + convergence.second);

            if(iter % 500 == 0)
                std::cout << "  i : "<< iter <<"  - " << std::sqrt(cvg_total)<<std::endl;

            assert(cvg_total < Ninf);
            if( cvg_total < tolerance)
            {
                std::cout << "  i : "<< iter <<"  - " << cvg_total <<std::endl;
                post_processing( msh, assembler, info, problem);
                return true;
            }
            sol_old = sol;

        }
        return false;
    }

};
