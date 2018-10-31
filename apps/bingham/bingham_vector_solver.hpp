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

#include "core/loaders/loader.hpp"

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"
#include "output/postMesh.hpp"

#include "output/silo.hpp"
#include "solvers/solver.hpp"
#include "viscoplasticity_utils.hpp"


template<typename Mesh>
class ADDM
{
    typedef Mesh mesh_type;
    typedef typename mesh_type::coordinate_type     T;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;
    typedef typename mesh_type::point_type      point_type;

    typedef disk::mechanics::BoundaryConditions<mesh_type> boundary_type;
    typedef Matrix<T, Dynamic, Dynamic>         matrix_type;
    typedef Matrix<T, Dynamic, 1>               vector_type;

    vector_type             multiplier, auxiliar;
    bingham_data<T, vector_problem>     vp;
    hho_degree_info         di;
    size_t        cbs, fbs, pbs, sbs, dim;

public:
    vector_type             sol, sol_old;
    std::pair<T, T>         convergence;

    ADDM(const Mesh  & msh,
         const hho_degree_info  & hdi,
         const bingham_data< typename Mesh::coordinate_type , vector_problem> & vp_ext):
         di(hdi), vp(vp_ext)
    {
        dim =  Mesh::dimension;

        cbs = revolution::vector_basis_size(di.cell_degree(), dim, dim);
        fbs = revolution::vector_basis_size(di.face_degree(), dim - 1, dim);
        pbs = revolution::scalar_basis_size(di.face_degree(), dim);
        sbs = revolution::sym_matrix_basis_size(di.face_degree(), dim, dim);
    };

    template<typename Assembler>
    vector_type
    compute_auxiliar(   const mesh_type& msh,
                        const cell_type& cl,
                        const Assembler& assembler,
                        const vector_type& velocity)
    {
        vector_type u_TF  = assembler.take_velocity(msh, cl, velocity);
        auto value = 1./(2. * (vp.mu + vp.alpha));
        auto G = revolution::make_hlow_stokes(msh, cl, di, true);

        auto cell_ofs    = revolution::priv::offset(msh, cl);
        vector_type   Gu = G.first * u_TF;
        vector_type   stress = multiplier.block(cell_ofs * sbs, 0, sbs, 1);

        //Theta
        auto sb = revolution::make_sym_matrix_monomial_basis(msh, cl, di.face_degree());
        //barycenter only for k = 0;
        auto bar   = barycenter(msh, cl);
        auto s_phi = sb.eval_functions(bar);

        vector_type  theta  = stress  +  2. * vp.alpha * Gu;
        matrix_type  theta_eval = revolution::eval(theta, s_phi);
        T theta_norm  = std::sqrt((theta_eval.cwiseProduct(theta_eval)).sum());
        T theta_eigen = theta_eval.norm();

        assert(abs(theta_norm  - theta_eigen) < 1.e-7);
        assert(abs(theta_norm  - theta.norm())< 1.e-7);

        vector_type gamma = vector_type::Zero(sbs);
        T epsilon = 1.e-6 *std::sqrt(2) * vp.yield ;
        if( theta_norm > std::sqrt(2) * vp.yield + epsilon )
            gamma =  value * theta * (1. - std::sqrt(2) * (vp.yield/theta_norm));

        return gamma;
    }

    template<typename Assembler>
    vector_type
    make_rhs_alg(   const mesh_type& msh,
                    const cell_type& cl,
                    const Assembler& assembler)
    {
        auto G  = revolution::make_hlow_stokes(msh, cl, di, true);
        auto cb = revolution::make_vector_monomial_basis(msh, cl, di.cell_degree());
        auto sb = revolution::make_sym_matrix_monomial_basis(msh, cl, di.face_degree());

        auto cell_ofs =  revolution::priv::offset(msh, cl);
        auto num_faces = howmany_faces(msh, cl);

        vector_type rhs = vector_type::Zero(cbs + fbs * num_faces);

        //1. (stress - 2 alpha * gamma, Gv)
        vector_type   stress = multiplier.block( sbs * cell_ofs,  0, sbs, 1);
        vector_type   gamma  = compute_auxiliar( msh,  cl, assembler, sol_old);
        auxiliar.block(cell_ofs * sbs, 0, sbs, 1) = gamma;

        vector_type   str_agam = stress - 2. * vp.alpha * gamma;
        matrix_type   mm = revolution::make_mass_matrix(msh, cl, sb);

        rhs = -G.first.transpose() * mm * str_agam;

        //2. (f, v_T)
        auto rhs_fun_test  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
            return Matrix<T, Mesh::dimension, 1>::Zero();
        };
        rhs.block( 0, 0, cbs, 1) += make_rhs(msh, cl, cb, rhs_fun_test);

        return rhs;
    }

    template<typename Assembler>
    void
    make_global_rhs(const mesh_type& msh, Assembler& assembler)
    {
        assembler.initialize_rhs();

        for (auto& cl : msh)
        {
            vector_type local_rhs = make_rhs_alg(msh, cl, assembler);
            assembler.assemble_rhs(msh, cl, local_rhs);
        }
        assembler.finalize_rhs();
    }

    template<typename Assembler>
    void
    make_global_matrix(const mesh_type& msh, Assembler& assembler)
    {
        assembler.initialize_lhs();

        for (auto& cl : msh)
        {
            auto gr = revolution::make_hho_stokes(msh, cl, di, true);
            auto G  = revolution::make_hlow_stokes(msh, cl, di, true);
            matrix_type dr = make_hho_divergence_reconstruction_stokes_rhs(msh, cl, di);

            matrix_type stab = make_hho_vector_stabilization(msh, cl, gr.first, di);
            matrix_type A    = 2. * (vp.alpha * G.second + vp.mu * stab);

            assembler.assemble_lhs(msh, cl, A, -dr);
        }

        assembler.finalize_lhs();
    }

    template<typename Assembler>
    void
    run_stokes_like(const mesh_type& msh, Assembler& assembler,
                    const size_t iter)
    {

        make_global_rhs(msh, assembler);
        if(iter == 0)
            make_global_matrix(msh, assembler);

        //dump_sparse_matrix(assembler.LHS, "stokes.txt");
        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        vector_type solution = vector_type::Zero(systsz);
        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = false;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, solution);

        sol = vector_type::Zero(systsz);
        sol = solution;
    }

    template<typename Assembler>
    void
    update_multiplier(const mesh_type& msh, const Assembler& assembler,
                        const size_t iter,
                        const std::string& name)
    {
        T conv_residue = 0.;
        T conv_stress  = 0.;
        T conv_grad    = 0.;

        auto cell_id = 0;
        const auto dim = Mesh::dimension;

        for(auto& cl : msh)
        {
            vector_type u_now = assembler.take_velocity(msh, cl, sol);
            vector_type u_old = assembler.take_velocity(msh, cl, sol_old);

            auto G = revolution::make_hlow_stokes(msh, cl, di, true);

            vector_type  Guold = G.first * u_old;
            vector_type  Gunow = G.first * u_now;

            vector_type  gamma = auxiliar.block(cell_id * sbs, 0, sbs, 1);
            vector_type  diff_stress = 2. * vp.alpha * (Gunow - gamma);

            //Updating multiplier
            multiplier.block(cell_id * sbs, 0, sbs, 1) += diff_stress;

            //convergence:  | sigma^{n+1} -  sigma^n|^2
            auto sb = revolution::make_sym_matrix_monomial_basis(msh, cl, di.face_degree());
            matrix_type mass = make_mass_matrix(msh, cl, sb);

            vector_type  diff_grad =  vp.alpha * (Gunow - Guold);

            conv_stress  += diff_stress.dot(mass * diff_stress);
            conv_grad    += diff_grad.dot(mass * diff_grad);
            conv_residue += conv_stress +  conv_grad;
            cell_id++;
        }

        std::ofstream ifs(name, std::ofstream::app);
        if(!ifs.is_open())
            std::cout << "Error opening   " << std::endl;

        ifs << tostr(iter) <<"  "<< std::sqrt(conv_residue) << "  ";
        ifs << std::sqrt(conv_stress) <<"  "<< std::sqrt(conv_grad)<< std::endl;
        ifs.close();

        convergence = std::make_pair(std::sqrt(conv_residue), std::sqrt(conv_stress));
    }


    template<typename Assembler>
    void
    post_processing(const mesh_type& msh, Assembler& assembler, const size_t iter)
    {
        auto dim = Mesh::dimension;
        auto rbs = vector_basis_size(di.reconstruction_degree(), dim, dim);

        vector_type sol_cell = vector_type::Zero(cbs * msh.cells_size());
        vector_type sol_rec  = vector_type::Zero(rbs * msh.cells_size());
        vector_type pressure = vector_type::Zero(pbs * msh.cells_size());

        std::string  ext =  vp.info + "_i"+ tostr(iter) + ".data";
        std::ofstream ofs("data_" + ext );
        if (!ofs.is_open())
            std::cout << "Error opening file"<<std::endl;

        std::vector<T> ux;        ux.reserve(msh.cells_size());
        std::vector<T> uy;        uy.reserve(msh.cells_size());
        std::vector<T> e_press;   e_press.reserve(msh.cells_size());
        std::vector<T> e_theta;   e_theta.reserve(msh.cells_size());
        std::vector<T> e_sigma;   e_sigma.reserve(msh.cells_size());
        //std::vector<T> div;     div.reserve(msh.cells_size());

        auto cell_id = 0;
        for(auto cl : msh)
        {
            auto gr  = make_hho_stokes(msh, cl, di, true);
            vector_type svel =  assembler.take_velocity(msh, cl, sol);
            assert((gr.first * svel).rows() == rbs - dim);

            sol_rec.block(cell_id * rbs + dim, 0, rbs - dim, 1) = gr.first * svel;
            sol_rec.block(cell_id * rbs, 0, dim, 1) = svel.block(0,0, dim, 1);
            sol_cell.block(cell_id * cbs, 0, cbs, 1) = svel.block(0,0, cbs, 1);

            //this is only for k = 0, since there is only one velocity;
            auto bar = barycenter(msh, cl);

            //Velocity
            vector_type cell_vel = svel.block(0,0, cbs, 1);
            auto cb  = revolution::make_vector_monomial_basis(msh, cl, di.cell_degree());
            auto phi = cb.eval_functions(bar);
            vector_type ueval = revolution::eval(cell_vel, phi);

            //Pressure
            vector_type spress =  assembler.take_pressure(msh, cl, sol);
            auto pb  = revolution::make_scalar_monomial_basis(msh, cl, di.face_degree());
            auto p_phi = pb.eval_functions(bar);
            T peval =  p_phi.dot(spress);

            //Stress
            auto sb = revolution::make_sym_matrix_monomial_basis(msh, cl, di.face_degree());
            auto s_phi  = sb.eval_functions(bar);
            auto G = revolution::make_hlow_stokes(msh, cl, di, true);
            vector_type   Gu = G.first * svel;
            vector_type   stress = multiplier.block(cell_id * sbs, 0, sbs, 1);
            vector_type   theta  = stress  +  2. * vp.alpha * Gu;
            matrix_type   theta_eval = revolution::eval(theta, s_phi);
            matrix_type   sigma_eval = revolution::eval(stress, s_phi);
            matrix_type   Gu_eval = revolution::eval(Gu, s_phi);

            assert(abs(theta_eval.norm() - theta.norm() )<1.e-7);
            assert(abs(sigma_eval.norm() - stress.norm())<1.e-7);
            assert(abs(Gu_eval.norm() - Gu.norm())<1.e-7);

            T divu = Gu_eval(0,0) + Gu_eval(1,1);
            T tr_stress = sigma_eval(0,0) + sigma_eval(1,1);

            ux.push_back( ueval(0) );
            uy.push_back( ueval(1) );
            e_press.push_back( peval );
            e_theta.push_back( theta_eval.norm() );
            e_sigma.push_back( sigma_eval.norm() );
            //div.push_back( divu );

            ofs << ueval(0)   << " " << ueval(1) << " " << peval<< " ";
            ofs << theta_eval.norm() << " " << sigma_eval.norm()   << " ";
            ofs << divu << " "<< tr_stress<<std::endl;

            cell_id++;
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

        p_x = std::make_pair(point_type({0.0 + eps, 0.5 + eps}), point_type({1.0 + eps, 0.5 + eps}));
        p_y = std::make_pair(point_type({0.5 + eps, 0.0 + eps}), point_type({0.5 + eps, 1.0 + eps}));

        //plot_over_line(msh, p_x, cell_rec_sol, di.reconstruction_degree(), "plot_over_x_" + ext);
        //plot_over_line(msh, p_y, cell_rec_sol, di.reconstruction_degree(), "plot_over_y_" + ext);

        if (vp.problem == DRIVEN)
        {
            std::vector<point_type> pts;
            std::array<T, 21> py = {1., 0.975, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7,
                                    0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35,
                                    0.3, 0.25, 0.2, 0.15, 0.1, 0.05};

            for(size_t i = 0; i< 21; i++)
                pts.push_back(point_type({0.5, py[i]}));

            find_values_at_points(msh, pts, sol_rec, di.reconstruction_degree(), "error_rec" + ext);
            find_values_at_points(msh, pts, sol_cell, di.cell_degree(), "error_sol" + ext);
            //plot_over_line(msh, p_y, sol_cell, di.cell_degree(), "plot_over_y_" + ext);
        }
        //compute_discontinuous_velocity( msh, sol_cell, di, "velocity_" + vp.info +".msh");
        //save_coords(msh, "Coords_"+ ext);
        //quiver( msh, sol, assembler, di, "quiver_"+ ext);

        return;
    }


    bool
    run(const mesh_type& msh, const boundary_type& bnd)
    {
        auto assembler = revolution::make_stokes_assembler_alg(msh, di, bnd);
        auto systsz    = assembler.global_system_size();

        sol         =  vector_type::Zero(systsz);
        sol_old     =  vector_type::Zero(systsz);
        multiplier  =  vector_type::Zero(sbs * msh.cells_size());
        auxiliar    =  vector_type::Zero(sbs * msh.cells_size());

        size_t max_iters = 50000;
        T Ninf      = 100000;
        T tolerance = 1.e-8;

        std::string data_filename = "convg" + vp.info + ".data";
        std::ofstream ofs(data_filename);
        if( !ofs.is_open())
            std::cout << "Error opening data-file" << std::endl;
        ofs.close();

        Eigen::PardisoLU<Eigen::SparseMatrix<T>>  solver;

        for(size_t iter = 0; iter < max_iters; iter++)
        {
            sol_old   = sol;
            auxiliar  = vector_type::Zero(sbs * msh.cells_size());
            //------------------------------------------------------------------
            //run_stokes_like(msh, assembler, iter);
            make_global_rhs(msh, assembler);
            if (iter == 0)
            {
                make_global_matrix(msh, assembler);
                solver.analyzePattern(assembler.LHS);
                solver.factorize(assembler.LHS);
            }
            sol = solver.solve(assembler.RHS);
            //------------------------------------------------------------------
            update_multiplier(msh, assembler, iter, data_filename);
            //------------------------------------------------------------------
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
