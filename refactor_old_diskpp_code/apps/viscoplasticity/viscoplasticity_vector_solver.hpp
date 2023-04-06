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

template<typename Mesh>
class ADMM
{
    typedef Mesh mesh_type;
    typedef typename mesh_type::coordinate_type     T;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;
    typedef typename mesh_type::point_type  point_type;

    typedef disk::vector_boundary_conditions<mesh_type> boundary_type;
    typedef Matrix<T, Mesh::dimension, Mesh::dimension>         tensor_type;
    typedef Matrix<T, Mesh::dimension, 1>               vector2d_type;
    typedef Matrix<T, Dynamic, Dynamic>                 matrix_type;
    typedef Matrix<T, Dynamic, 1>                       vector_type;

    hho_degree_info     di;
    size_t              cbs, fbs, pbs, sbs;
    matrix_type         multiplier;
    matrix_type         auxiliar;
    matrix_type         auxiliar_old;
    bingham_data<T, vector_problem>     vp;

    tensors_at_quad_pts_utils<mesh_type>    tsr_utils;
    std::vector<std::pair<size_t, size_t>>  tsr_offsets_vector;


public:
    disk::dynamic_vector<T>   sol_old;
    std::pair<T, T>     convergence;

    ADMM(const  Mesh& msh,
        const   hho_degree_info & hdi,
        const   bingham_data< T , vector_problem> & vp_ext):
        di(hdi), vp(vp_ext)
    {
        const auto dim =  Mesh::dimension;

        cbs = vector_basis_size(di.cell_degree(), dim, dim);
        fbs = vector_basis_size(di.face_degree(), dim - 1, dim);
        pbs = scalar_basis_size(di.face_degree(), dim);
        sbs = sym_matrix_basis_size(di.face_degree(), dim, dim);

        size_t quad_degree = 2. * di.face_degree();
        tsr_utils = tensors_at_quad_pts_utils<mesh_type>(msh, quad_degree);
        tsr_offsets_vector = tsr_utils.offsets_vector();
    };

    template<typename Assembler>
    matrix_type
    compute_auxiliar(   const mesh_type & msh,
                        const cell_type & cl,
                        const Assembler & assembler,
                        const vector_type& velocity_dofs)
    {
        vector_type u_TF  = assembler.take_velocity(msh, cl, velocity_dofs);
        auto value = 1./(2. * (vp.mu + vp.alpha));
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
            vector_type theta     = stress_qp  +  2. * vp.alpha * Gu;

            tensor_type theta_eval = eval(theta, s_phi);

            T theta_norm  = theta_eval.norm();
            T tol = 1.e-8;

            // A. Liquid
            if( (theta_norm + tol) >  std::sqrt(2) * vp.yield)
            {
                gamma.block( 0, qp_count, sbs, 1) = value * theta *
                                        (1. - std::sqrt(2) * (vp.yield/theta_norm));
            }
            qp_count++;
        }
        assert(qps.size() == qp_count);

        return gamma;

    }

    template<typename Assembler>
    auto
    update_multiplier(const mesh_type& msh, const Assembler& assembler,
                      const size_t iter,
                      const disk::dynamic_vector<T>& velocity_dofs,
                      const std::string& name)
    {
        T conv_stress = 0.;
        T conv_gamma = 0.;
        T conv_total = 0;
        const auto dim = Mesh::dimension;

        for(auto cl: msh)
        {
            auto sb = make_sym_matrix_monomial_basis(msh, cl, di.face_degree());

            vector_type u_TF = assembler.take_velocity(msh, cl, velocity_dofs);
            auto G = make_hlow_stokes(msh, cl, di, true);
            vector_type Gu = G.first * u_TF;

            auto cl_id = msh.lookup(cl);
            auto offset =  tsr_offsets_vector.at(cl_id).first;

            auto qps = integrate(msh, cl, tsr_utils.quad_degree());

            matrix_type     gamma_old = auxiliar_old.block(0, offset, sbs, qps.size());
            matrix_type     gamma     = auxiliar.block(0, offset, sbs, qps.size());
            matrix_type     diff_gamma  = vp.alpha * (gamma - gamma_old);
            matrix_type     diff_stress = - 2. * vp.alpha *  gamma;

            //Update multiplier
            for(size_t i = 0; i < qps.size(); i++)
                diff_stress.block(0, i, sbs, 1) += 2. * vp.alpha *  Gu;
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
                conv_total  += conv_stress + conv_gamma;

                qp_count++;
            }
        }

        std::ofstream ifs(name, std::ofstream::app);
        if(!ifs.is_open())
            std::cout << "Error opening   " << std::endl;

        ifs << tostr(iter) <<"  "<< std::sqrt(conv_total) << "  ";
        ifs << std::sqrt(conv_stress) <<"  "<< std::sqrt(conv_gamma)<< std::endl;
        ifs.close();

        convergence = std::make_pair(std::sqrt(conv_total), std::sqrt(conv_stress));

        auxiliar_old = auxiliar;
        return;
    }

    template<typename Assembler>
    Matrix<T, Dynamic, 1>
    make_rhs_admm(   const mesh_type& msh,
                    const cell_type& cl,
                    const Assembler& assembler)
    {
        auto G = make_hlow_stokes(msh, cl, di, true);
        auto cb = make_vector_monomial_basis(msh, cl, di.cell_degree());
        auto sb = make_sym_matrix_monomial_basis(msh, cl, di.face_degree());

        auto cell_ofs =  offset(msh, cl);
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

        matrix_type str_agam = stress - 2. * vp.alpha * gamma;

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

        auto rhs_fun_test  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
            return Matrix<T, Mesh::dimension, 1>::Zero();
        };

        //(f, v_T)
        rhs.block( 0, 0, cbs, 1) += make_rhs(msh, cl, cb, rhs_fun_test);

        return rhs;
    }

    template<typename Assembler>
    void
    make_global_rhs(const mesh_type& msh, Assembler& assembler)
    {
        assembler.initialize_rhs();

        for (auto cl : msh)
        {
            vector_type local_rhs = make_rhs_admm(msh, cl, assembler);
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

            matrix_type A = 2.* vp.alpha * ( G.second + stab);

            assembler.assemble_lhs(msh, cl, A, -dr);
        }

        assembler.finalize_lhs();

        return;
    }

    template<typename Assembler>
    void
    post_processing(const mesh_type& msh, const Assembler& assembler, const size_t iter)
    {
        auto dim = Mesh::dimension;
        auto rbs = vector_basis_size(di.reconstruction_degree(), dim, dim);

        disk::dynamic_vector<T> press_vec = disk::dynamic_vector<T>::Zero(pbs * msh.cells_size());
        disk::dynamic_vector<T> cell_sol  = disk::dynamic_vector<T>::Zero(cbs * msh.cells_size());
        disk::dynamic_vector<T> cell_rec_sol = disk::dynamic_vector<T>::Zero(rbs * msh.cells_size());

        std::string  ext =  vp.info + "_i"+ tostr(iter) + ".data";
        std::ofstream ofs("data_" + ext );
        if (!ofs.is_open())
            std::cout << "Error opening file"<<std::endl;

        std::vector<T> ux;        ux.reserve(msh.cells_size());
        std::vector<T> uy;        uy.reserve(msh.cells_size());
        std::vector<T> e_press;   e_press.reserve(msh.cells_size());
        std::vector<T> e_theta;   e_theta.reserve(msh.cells_size());
        std::vector<T> e_sigma;   e_sigma.reserve(msh.cells_size());

        auto cell_id = 0;
        for(auto cl : msh)
        {
            auto gr  = disk::make_hho_stokes(msh, cl, di, true);

            vector_type svel =  assembler.take_velocity(msh, cl, sol_old);
            assert((gr.first * svel).rows() == rbs - dim);

            cell_rec_sol.block(cell_id * rbs + dim, 0, rbs - dim, 1) = gr.first * svel;
            cell_rec_sol.block(cell_id * rbs, 0, dim, 1) = svel.block(0,0, dim, 1);
            cell_sol.block(cell_id * cbs, 0, cbs, 1)     = svel.block(0,0, cbs, 1);

            auto G  = disk::make_hlow_stokes(msh, cl, di, true);
            auto cb = disk::make_vector_monomial_basis(msh, cl, di.cell_degree());
            auto pb = disk::make_scalar_monomial_basis(msh, cl, di.face_degree());
            auto sb = disk::make_sym_matrix_monomial_basis(msh, cl, di.face_degree());


            auto offset = tsr_offsets_vector.at(cell_id).first;
            auto qps    = integrate(msh, cl, tsr_utils.quad_degree());

            vector_type Gu  = G.first * svel;
            vector_type cell_vel = svel.block(0,0, cbs, 1);
            vector_type spress   = assembler.take_pressure(msh, cl, sol_old);
            matrix_type stress   = multiplier.block(0, offset, sbs, qps.size());

            auto qps_count = 0;

            for(auto& qp : qps)
            {
                auto v_phi  = cb.eval_functions(qp.point());
                auto p_phi  = pb.eval_functions(qp.point());
                auto s_phi  = sb.eval_functions(qp.point());

                vector_type vel_eval = disk::eval(cell_vel, v_phi);
                T press_eval =  p_phi.dot(spress);

                //Stress
                vector_type stress_qp  = stress.block(0, qps_count, sbs, 1);
                vector_type theta_qp   = stress_qp  +  2. * vp.alpha * Gu;
                tensor_type theta_eval = disk::eval(theta_qp,  s_phi);
                tensor_type sigma_eval = disk::eval(stress_qp, s_phi);
                tensor_type grad_eval  = disk::eval(Gu, s_phi);

                T divergence   = grad_eval(0,0)  + grad_eval(1,1);
                T trace_stress = sigma_eval(0,0) + sigma_eval(1,1);

                ofs << qp.point().x() << " " << qp.point().y() << " ";
                ofs << vel_eval(0)    << " " << vel_eval(1) << " " << press_eval<< " ";
                ofs << theta_eval.norm() << " " << sigma_eval.norm()   << " ";
                ofs << divergence << " " << trace_stress << std::endl;

                if( di.cell_degree() == 0 && di.face_degree() == 0 )
                {
                    e_theta.push_back( theta_eval.norm() );
                    e_sigma.push_back( sigma_eval.norm() );
                }

                qps_count++;
            }


            auto bar = barycenter(msh, cl);
            auto v_phi  = cb.eval_functions(bar);
            auto p_phi  = pb.eval_functions(bar);

            auto u_bar_x = cell_vel(0) * v_phi(0,0);
            auto u_bar_y = cell_vel(1) * v_phi(1,1);
            auto p_bar   = spress(0) * p_phi(0);

            ux.push_back( u_bar_x);
            uy.push_back( u_bar_y);
            e_press.push_back( p_bar );

        }
        ofs.close();

        std::string silo_db_name = "bingham_" + tostr(iter) + ".silo";
        disk::silo_database silo;
        silo.create(silo_db_name);
        silo.add_mesh(msh, "mesh");

        disk::silo_zonal_variable<T> silo_ux("ux", ux);
        disk::silo_zonal_variable<T> silo_uy("uy", uy);
        disk::silo_zonal_variable<T> silo_press("pressure", e_press);
        disk::silo_zonal_variable<T> silo_theta("theta", e_theta);
        disk::silo_zonal_variable<T> silo_sigma("sigma", e_sigma);
        //disk::silo_zonal_variable<T> silo_div("pressure", div);

        silo.add_variable("mesh", silo_ux);
        silo.add_variable("mesh", silo_uy);
        silo.add_variable("mesh", silo_press);
        silo.add_variable("mesh", silo_theta);
        silo.add_variable("mesh", silo_sigma);

        silo.close();

        std::pair<point_type, point_type> p_x, p_y;
        auto eps = 1.e-4;
        p_y = std::make_pair(point_type({0.5 + eps, 0.0 + eps}), point_type({0.5 + eps, 1.0 + eps}));

        //plot_over_line(msh, p_x, cell_rec_sol, di.reconstruction_degree(), "plot_over_x_" + info + ".data");
        //plot_over_line(msh, p_y, cell_rec_sol, di.reconstruction_degree(), "plot_over_y_" + info + ".data");
        plot_over_line(msh, p_y, cell_sol, di.cell_degree(), "plot_over_y_" + ext);

        //compute_discontinuous_velocity( msh, cell_sol, di, "velocity_" + vp.info + "_i"+ tostr(iter) +".msh");
        save_coords(msh, "Coords_"+ ext);
        //quiver( msh, sol_old, assembler, di, "quiver_"+ vp.info + "_i"+ tostr(iter) + ".data");

        paraview<mesh_type> pw(di.cell_degree());
        //pp.paraview(msh, "Velocity_Magnitud", Vmagnitud, "scalar");
        pw.make_file(msh, "Velocity_"+ vp.info + "_i"+ tostr(iter), cell_sol, "vector");
        std::cout << "Velocity_"+ vp.info + "_i"+ tostr(iter) << std::endl;

        return;
    }

    bool
    run(const mesh_type& msh, const boundary_type& bnd)
    {
        auto assembler = disk::make_stokes_assembler_alg(msh, di, bnd);
        auto systsz = assembler.global_system_size();

        sol_old = vector_type::Zero(systsz);
        vector_type sol =  vector_type::Zero(systsz);

        size_t num_total_quads = tsr_utils.num_total_quad_points();

        multiplier   = matrix_type::Zero(sbs, num_total_quads);
        auxiliar     = matrix_type::Zero(sbs, num_total_quads);
        auxiliar_old = matrix_type::Zero(sbs, num_total_quads);

        auto Ninf = 1.e+10;
        auto max_iters = 500000;
        auto tolerance = 1.e-8;

        std::string data_filename = "convg" + vp.info + ".data";
        std::ofstream ofs(data_filename);
        if( !ofs.is_open())
            std::cout << "Error opening data-file" << std::endl;
        ofs.close();

        Eigen::PardisoLU<Eigen::SparseMatrix<T>>  solver;

        for(size_t iter = 0; iter < max_iters; iter++)
        {
            //------------------------------------------------------------------
            // Stokes-like system

            if(iter == 0)
            {
                make_global_matrix(msh, assembler);
                solver.analyzePattern(assembler.LHS);
                solver.factorize(assembler.LHS);
            }
            //WARNINGS: This one must go after make_global_matrix!!!!!!
            make_global_rhs(msh,assembler);

            vector_type sol = solver.solve(assembler.RHS);
            //------------------------------------------------------------------
            update_multiplier(msh, assembler, iter, sol,  data_filename);
            //------------------------------------------------------------------
            //Erase cout
            std::cout << "  i : "<< iter <<"  - " << std::sqrt(convergence.first) << " ";
                std::cout << std::sqrt(convergence.second) << std::endl;

            sol_old = sol;

            assert(convergence.first < Ninf);
            bool stop = convergence.first < tolerance;
            if(! (iter % 100 == 0 ||  stop))
                continue;

            std::cout << "  i : "<< iter <<"  - " << convergence.first <<std::endl;
            post_processing( msh, assembler, iter);

            if(stop)
                return true;
            //------------------------------------------------------------------
        }
        return false;
    }

};
