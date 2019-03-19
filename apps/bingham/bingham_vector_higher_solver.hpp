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
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;
    typedef typename mesh_type::point_type      point_type;

    typedef disk::BoundaryConditions<mesh_type, false>              boundary_type;
    typedef Matrix<T, Mesh::dimension, Mesh::dimension>         tensor_type;
    typedef Matrix<T, Dynamic, Dynamic>                         matrix_type;
    typedef Matrix<T, Dynamic, 1>                               vector_type;

    hho_degree_info         di;
    bingham_data<T, vector_problem>     vp;
    eigen_compatible_stdvector<matrix_type>         multiplier;
    eigen_compatible_stdvector<matrix_type>         auxiliar;
    eigen_compatible_stdvector<matrix_type>         auxiliar_old;

public:
    vector_type             sol_old, sol;
    std::tuple<T, T, T, T, T>         convergence;

    ADMM(const Mesh  & msh,
         const hho_degree_info  & hdi,
         const bingham_data< typename Mesh::coordinate_type , vector_problem> & vp_ext):
         di(hdi), vp(vp_ext)
    {}

    template<typename Assembler>
    void
    compute_auxiliar(   const mesh_type& msh,
                        const Assembler& assembler,
                        const size_t iter)
    {
        const auto dim =  Mesh::dimension;
        auto sbs = sym_matrix_basis_size(di.face_degree(), dim, dim);

        for(auto& cl : msh)
        {
            vector_type u_TF  = assembler.take_velocity(msh, cl, sol_old);
            auto value = 1./(2. * (vp.mu + vp.alpha));
            //auto G = make_hlow_stokes(msh, cl, di, true);
            auto G  = make_matrix_symmetric_gradrec(msh, cl, di);
            auto cl_id = msh.lookup(cl);

            auto sb  = make_sym_matrix_monomial_basis(msh, cl, di.face_degree());
            auto qps = integrate(msh, cl, 2 * di.face_degree());

            vector_type   Gu = G.first * u_TF;
            matrix_type   stress = multiplier.at(cl_id);
            matrix_type   gamma  = matrix_type::Zero(sbs, qps.size());

            //Theta
            matrix_type  theta = stress;
            for(size_t i = 0; i < qps.size(); i++)
                theta.block(0, i, sbs, 1) += 2. * vp.alpha *  Gu;

            size_t qp_count    = 0;
            for(auto& qp: qps)
            {
                vector_type theta_qp = theta.block( 0, qp_count, sbs, 1);

                auto s_phi  = sb.eval_functions(qp.point());
                tensor_type theta_eval = disk::eval(theta_qp, s_phi);
                T theta_norm  = theta_eval.norm();

                T epsilon = 1.e-8 ;

                bool running_stokes = (vp.yield < 1.e-12 || iter == 0)? true : false;

                if( theta_norm > std::sqrt(2) * vp.yield * (1 + epsilon) && !running_stokes)
                {
                    gamma.block( 0, qp_count, sbs, 1) =  value * theta_qp * (1. - std::sqrt(2) * (vp.yield/theta_norm));
                }

                qp_count++;
            }
            auxiliar.at(cl_id)  = gamma;
        }

        return;
    }

    template<typename Assembler>
    void
    initialize_auxiliar(const mesh_type& msh, const Assembler& assembler)
    {
        auto dim =  Mesh::dimension;
        auto sbs = disk::sym_matrix_basis_size(di.face_degree(), dim, dim);
        for (auto& cl : msh)
        {
            auto cl_id = msh.lookup(cl);
            auto G   = disk::make_hlow_stokes(msh, cl, di, true);
            //auto G   = make_hho_sym_gradrec_matrix(msh, cl, di);
            auto qps = integrate(msh, cl, 2 * di.face_degree());
            vector_type unow = assembler.take_velocity(msh, cl, sol);
            matrix_type gamma = matrix_type::Zero(sbs, qps.size());

            auto qp_count = 0;
            for(auto& qp: qps)
            {
                gamma.block( 0, qp_count, sbs, 1) = G.first * unow;
                auxiliar.at(cl_id)  = gamma;
                qp_count++;
            }
        }
        return;
    }

    template<typename Assembler>
    auto
    update_multiplier(const mesh_type& msh, const Assembler& assembler,
                      const size_t iter,
                      const std::string& name)
    {
        T residue_sig_gam(0), residue_sig_grad(0);
        T residue_stress(0),  residue_gamma(0),  residue_grad(0);

        const auto dim = Mesh::dimension;
        auto cbs = vector_basis_size(di.cell_degree(), dim, dim);
        auto fbs = vector_basis_size(di.face_degree(), dim - 1, dim);
        auto sbs = sym_matrix_basis_size(di.face_degree(), dim, dim);

        auto coef_admm = (iter == 0)? vp.mu : vp.alpha;

        auto cell_id = 0;
        for(auto cl: msh)
        {
            auto qps = integrate(msh, cl, 2 * di.face_degree());

            vector_type u_now = assembler.take_velocity(msh, cl, sol);
            vector_type u_old = assembler.take_velocity(msh, cl, sol_old);

            //auto G = disk::make_hlow_stokes(msh, cl, di, true);
            auto G = make_matrix_symmetric_gradrec(msh, cl, di);
            vector_type  Guold = G.first * u_old;
            vector_type  Gunow = G.first * u_now;
            vector_type  diff_grad =  coef_admm * (Gunow - Guold);

            matrix_type  gamma_old = auxiliar_old.at(cell_id);
            matrix_type  gamma     = auxiliar.at(cell_id);
            matrix_type  diff_gamma  = coef_admm * (gamma - gamma_old);

            matrix_type  diff_stress = - 2. * coef_admm *  gamma;
            for(size_t i = 0; i < qps.size(); i++)
                diff_stress.block(0, i, sbs, 1) += 2. * coef_admm *  Gunow;

            //Updating multiplier
            multiplier.at(cell_id) += diff_stress;

            //Convergence
            auto sb = make_sym_matrix_monomial_basis(msh, cl, di.face_degree());

            size_t qp_count = 0;
            for(auto& qp: qps)
            {
                auto s_phi  = sb.eval_functions(qp.point());
                matrix_type mm =  priv::outer_product(s_phi, s_phi);
                vector_type diff_stress_qp  = diff_stress.block(0, qp_count, sbs, 1);
                vector_type diff_gamma_qp   = diff_gamma.block(0, qp_count, sbs, 1);

                //1. | sigma^{n+1} -  sigma^n|^2
                residue_stress += qp.weight() * diff_stress_qp.dot(mm * diff_stress_qp);

                //2. | Gamma^{n+1} -  Gamma^n|^2
                residue_gamma  += qp.weight() * diff_gamma_qp.dot(mm * diff_gamma_qp);

                //3. | Gu^{n+1} -  Gu^n|^2
                residue_grad   += qp.weight() * diff_grad.dot(mm * diff_grad);

                qp_count++;
            }

            cell_id++;
        }

        //4. | sigma^{n+1} -  sigma^n|^2 + coef^2 | Gamma^{n+1} -  Gamma^n|^2
        residue_sig_gam =  residue_stress + residue_gamma;

        //5. | sigma^{n+1} -  sigma^n|^2 +  coef^2 | Gu^{n+1} -  Gu^n|^2
        residue_sig_grad =  residue_stress + residue_grad;

        std::ofstream ifs(name, std::ofstream::app);
        if(!ifs.is_open())
            std::cout << "Error opening   " << std::endl;


        auto sqrt_error = [] (const T& a ) -> T {
            if(a > -1.e-14)
                return std::sqrt(std::abs(a));
            else
                throw std::runtime_error("Negative sign for errors");
        };

        auto error_sig_gam  = sqrt_error(residue_sig_gam);
        auto error_sig_grad = sqrt_error(residue_sig_grad);
        auto error_stress   = sqrt_error(residue_stress);
        auto error_gamma    = sqrt_error(residue_gamma);
        auto error_grad     = sqrt_error(residue_grad);

        ifs << tostr(iter) <<"  "<< error_sig_gam << "  " << error_sig_grad;
        ifs << "  " << error_stress << "  "<< error_gamma << "  "<< error_grad;
        ifs << std::endl;
        ifs.close();

        convergence = std::make_tuple(error_sig_grad, error_sig_gam, error_stress,
            error_gamma, error_grad);
        return;
    }

    template<typename Assembler>
    vector_type
    make_rhs_admm(  const mesh_type& msh,
                    const cell_type& cl,
                    const Assembler& assembler,
                    const size_t iter)
    {
        auto dim =  Mesh::dimension;
        auto cbs = vector_basis_size(di.cell_degree(), dim, dim);
        auto fbs = vector_basis_size(di.face_degree(), dim - 1, dim);
        auto sbs = sym_matrix_basis_size(di.face_degree(), dim, dim);

        auto G  = make_hlow_stokes(msh, cl, di, true);
        //auto G  = make_hho_sym_gradrec_matrix(msh, cl, di);
        auto cb = make_vector_monomial_basis(msh, cl, di.cell_degree());
        auto sb = make_sym_matrix_monomial_basis(msh, cl, di.face_degree());

        auto cl_id = msh.lookup(cl);
        auto num_faces = howmany_faces(msh, cl);

        vector_type rhs = vector_type::Zero(cbs + fbs * num_faces);

        //1. (f, v_T)
        auto rhs_fun_test  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
            return Matrix<T, Mesh::dimension, 1>::Zero();
        };
        rhs.block( 0, 0, cbs, 1) = make_rhs(msh, cl, cb, rhs_fun_test);

        //2. (stress - 2 alpha * gamma, Gv)
        matrix_type   stress = multiplier.at(cl_id);
        matrix_type   gamma  = auxiliar.at(cl_id);
        assert(stress.cols() == gamma.cols() && stress.rows() == gamma.rows());
        matrix_type str_agam = stress - 2. * vp.alpha * gamma;

        size_t qp_count = 0;
        auto qps = integrate(msh, cl, 2 * di.face_degree());
        for(auto& qp: qps)
        {
            auto s_phi  = sb.eval_functions(qp.point());
            matrix_type mm = priv::outer_product(s_phi, s_phi);
            vector_type str_agam_qp  = str_agam.block(0, qp_count, sbs, 1);
            rhs -= qp.weight() * G.first.transpose() * mm * str_agam_qp;
            qp_count++;
        }

        return rhs;
    }

    template<typename Assembler>
    void
    make_global_rhs(const mesh_type& msh, Assembler& assembler, const size_t iter)
    {
        assembler.initialize_rhs();

        for (auto& cl : msh)
        {
            vector_type local_rhs = make_rhs_admm(msh, cl, assembler, iter);
            assembler.assemble_rhs(msh, cl, local_rhs);
        }
        assembler.finalize_rhs();
    }

    template<typename Assembler>
    void
    make_global_matrix(const mesh_type& msh, Assembler& assembler, const size_t iter)
    {
        assembler.initialize_lhs();

        //iter = 0 => Stokes, iter > 0 => Bingham
        auto coef_admm = (iter == 0)? vp.mu : vp.alpha;

        for (auto& cl : msh)
        {
            auto gr = disk::make_hho_stokes(msh, cl, di, true);
            //auto G  = disk::make_hlow_stokes(msh, cl, di, true);
            auto G  = make_matrix_symmetric_gradrec(msh, cl, di);
            matrix_type dr   = make_hho_divergence_reconstruction_rhs(msh, cl, di);
            matrix_type stab = make_vector_hho_stabilization(msh, cl, gr.first, di);
            matrix_type A    = 2. * ( coef_admm * G.second +  vp.mu * stab);

            assembler.assemble_lhs(msh, cl, A, -dr);
        }

        assembler.finalize_lhs();
    }

    template<typename Assembler>
    void
    post_processing(const mesh_type& msh, Assembler& assembler, const size_t iter,
                    const bool stop)
    {

        auto dim = Mesh::dimension;
        auto rbs = vector_basis_size(di.reconstruction_degree(), dim, dim);
        auto cbs = vector_basis_size(di.cell_degree(), dim, dim);
        auto fbs = vector_basis_size(di.face_degree(), dim - 1, dim);
        auto pbs = scalar_basis_size(di.face_degree(), dim);
        auto sbs = sym_matrix_basis_size(di.face_degree(), dim, dim);

        vector_type cell_sol = vector_type::Zero(cbs * msh.cells_size());
        vector_type cell_rec_sol  = vector_type::Zero(rbs * msh.cells_size());
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

            cell_rec_sol.block(cell_id * rbs + dim, 0, rbs - dim, 1) = gr.first * svel;
            cell_rec_sol.block(cell_id * rbs, 0, dim, 1)  = svel.block(0,0, dim, 1);
            cell_sol.block(cell_id * cbs, 0, cbs, 1) = svel.block(0,0, cbs, 1);

            //auto G  = make_hlow_stokes(msh, cl, di, true);
            auto G  = make_matrix_symmetric_gradrec(msh, cl, di);
            auto cb = make_vector_monomial_basis(msh, cl, di.cell_degree());
            auto pb = make_scalar_monomial_basis(msh, cl, di.face_degree());
            auto sb = make_sym_matrix_monomial_basis(msh, cl, di.face_degree());

            auto qps    = integrate(msh, cl, 2 * di.face_degree());

            vector_type Gu  = G.first * svel;
            vector_type cell_vel = svel.block(0,0, cbs, 1);
            vector_type spress   = assembler.take_pressure(msh, cl, sol);
            matrix_type stress   = multiplier.at(cell_id);

            auto qps_count = 0;

            for(auto& qp : qps)
            {
                auto v_phi  = cb.eval_functions(qp.point());
                auto p_phi  = pb.eval_functions(qp.point());
                auto s_phi  = sb.eval_functions(qp.point());

                vector_type ueval = eval(cell_vel, v_phi);
                T peval =  p_phi.dot(spress);

                //Stress
                vector_type stress_qp  = stress.block(0, qps_count, sbs, 1);
                vector_type theta_qp   = stress_qp  +  2. * vp.alpha * Gu;
                tensor_type theta_eval = eval(theta_qp,  s_phi);
                tensor_type sigma_eval = eval(stress_qp, s_phi);
                tensor_type grad_eval  = eval(Gu, s_phi);

                T divergence   = grad_eval(0,0)  + grad_eval(1,1);
                T trace_stress = sigma_eval(0,0) + sigma_eval(1,1);

                ofs << qp.point().x() << " " << qp.point().y() << " ";
                ofs << ueval(0)    << " " << ueval(1) << " " << peval<< " ";
                ofs << theta_eval.norm() << " " << sigma_eval.norm()   << " ";
                ofs << divergence << " " << trace_stress << std::endl;

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

            if( di.face_degree() == 0)
            {
                //Stress
                auto s_phi  = sb.eval_functions(bar);
                vector_type stress_bar  = stress.block(0, 0, sbs, 1);
                vector_type theta_bar  = stress_bar  +  2. * vp.alpha * Gu;
                tensor_type theta_eval = eval(theta_bar,  s_phi);
                tensor_type sigma_eval = eval(stress_bar, s_phi);

                e_theta.push_back( theta_eval.norm() );
                e_sigma.push_back( sigma_eval.norm() );
            }

            cell_id++;
        }

        ofs.close();

        std::string silo_db_name = "bingham_higher_" + tostr(iter) + ".silo";
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

        if (vp.problem == DRIVEN && (iter % 1000 == 0 ||  stop ))
        {
            std::vector<point_type> pts;
            std::array<T, 21> py = {1., 0.975, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7,
                                    0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35,
                                    0.3, 0.25, 0.2, 0.15, 0.1, 0.05};

            for(size_t i = 0; i< 21; i++)
                pts.push_back(point_type({0.5, py[i]}));

            find_values_at_points(msh, pts, cell_rec_sol, di.reconstruction_degree(), "error_rec" + ext);
            find_values_at_points(msh, pts, cell_sol, di.cell_degree(), "error_sol" + ext);
            //plot_over_line(msh, p_y, cell_rec_sol, di.reconstruction_degree(), "plot_over_y_" + ext);
            plot_over_line(msh, p_y, cell_sol, di.cell_degree(), "plot_over_y_" + ext);

            paraview<mesh_type> pw(di.cell_degree());
            //pp.paraview(msh, "Velocity_Magnitud", Vmagnitud, "scalar");
            pw.make_file(msh, "Velocity_"+ vp.info + "_i"+ tostr(iter), cell_rec_sol, "vector");
        }

        //compute_discontinuous_velocity( msh, cell_sol, di, "depl2d.msh");

        return;
    }


    bool
    run(const mesh_type& msh, const boundary_type& bnd)
    {
        auto assembler = disk::make_stokes_assembler_alg(msh, di, bnd);
        auto systsz    = assembler.global_system_size();

        sol_old = vector_type::Zero(systsz);
        sol =  vector_type::Zero(systsz);

        auto tensor_degree = di.face_degree();
        auto quad_degree = 2 * di.face_degree();
        multiplier   = tensor_initialize(msh, quad_degree, tensor_degree);
        auxiliar     = tensor_initialize(msh, quad_degree, tensor_degree);
        auxiliar_old = tensor_initialize(msh, quad_degree, tensor_degree);

        auto Ninf = 1.e+10;
        auto max_iters = 50000;
        auto tolerance = 1.e-8;

        std::string data_filename = "convg" + vp.info + ".data";
        std::ofstream ofs(data_filename);
        if( !ofs.is_open())
            std::cout << "Error opening data-file" << std::endl;
        ofs.close();

        Eigen::PardisoLDLT<Eigen::SparseMatrix<T>>  solver;

        for(size_t iter = 0; iter < max_iters; iter++)
        {
            sol_old   = sol;
            auxiliar_old = auxiliar;
            auxiliar  = tensor_initialize(msh, quad_degree, tensor_degree);
            //------------------------------------------------------------------
            compute_auxiliar(msh, assembler, iter);
            //------------------------------------------------------------------
            //solve stokes-like system
            if (iter == 0 || iter == 1)
            {
                make_global_matrix(msh, assembler, iter);
                solver.analyzePattern(assembler.LHS);
                solver.factorize(assembler.LHS);
            }
            //WARNINGS: This one must go after make_global_matrix!!!!!!
            make_global_rhs(msh, assembler, iter);
            sol = solver.solve(assembler.RHS);
            //------------------------------------------------------------------
            update_multiplier(msh, assembler, iter,  data_filename);
            //------------------------------------------------------------------
            if(iter == 0 )
                initialize_auxiliar(msh, assembler); // Just to decrease gamma error when iter = 1 (this doesnt change the convergence or behaviour afterwards)

            assert(std::get<0>(convergence) < Ninf && std::get<1>(convergence) < Ninf);

            bool stop = std::get<0>(convergence)< tolerance;
            if(! (iter % 100 == 0 ||  stop || iter == max_iters -1))
                continue;

            std::cout << "  i : "<< iter <<"  - " << std::get<0>(convergence)<<" ";
            std::cout << std::get<1>(convergence) <<"  "<< std::get<2>(convergence);
            std::cout << std::endl;
            post_processing( msh, assembler, iter, stop);

            if(stop)
                return true;
            //------------------------------------------------------------------
        }
        return false;
    }
};


template<typename Mesh>
class STOKES
{
    typedef Mesh mesh_type;
    typedef typename mesh_type::coordinate_type     T;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;
    typedef typename mesh_type::point_type      point_type;

    typedef disk::BoundaryConditions<mesh_type>      boundary_type;
    typedef Matrix<T, Mesh::dimension, Mesh::dimension>         tensor_type;
    typedef Matrix<T, Dynamic, Dynamic>                         matrix_type;
    typedef Matrix<T, Dynamic, 1>                               vector_type;

    hho_degree_info         di;
    bingham_data<T, vector_problem>     vp;

public:
    vector_type             sol_old, sol;
    std::tuple<T, T, T, T, T>         convergence;

    STOKES(const Mesh  & msh,
         const hho_degree_info  & hdi,
         const bingham_data< typename Mesh::coordinate_type , vector_problem> & vp_ext):
         di(hdi), vp(vp_ext)
         {}

        template<typename Assembler>
        void
        make_global_rhs(const mesh_type& msh, Assembler& assembler, const size_t iter)
        {
            assembler.initialize_rhs();
            assembler.finalize_rhs();
        }

        template<typename Assembler>
        void
        make_global_matrix(const mesh_type& msh, Assembler& assembler, const size_t iter)
        {
            assembler.initialize_lhs();

            //iter = 0 => Stokes, iter > 0 => Bingham
            auto coef_admm = (iter == 0)? vp.mu : vp.alpha;

            for (auto& cl : msh)
            {
                auto gr = disk::make_hho_stokes(msh, cl, di, true);
                auto G  = disk::make_hlow_stokes(msh, cl, di, true);
                //auto G  = make_matrix_symmetric_gradrec(msh, cl, di);
                matrix_type dr   = make_hho_divergence_reconstruction_rhs(msh, cl, di);
                matrix_type stab = make_vector_hho_stabilization(msh, cl, gr.first, di);
                matrix_type A    = 2. *  vp.mu * ( G.second +   stab);

                assembler.assemble_lhs(msh, cl, A, -dr);
            }

            assembler.finalize_lhs();
        }

        template<typename Assembler>
        void
        post_processing(const mesh_type& msh, Assembler& assembler, const size_t iter,
                        const bool stop)
        {
            auto dim = Mesh::dimension;
            auto rbs = vector_basis_size(di.reconstruction_degree(), dim, dim);
            auto cbs = vector_basis_size(di.cell_degree(), dim, dim);

            vector_type cell_sol = vector_type::Zero(cbs * msh.cells_size());
            vector_type cell_rec_sol  = vector_type::Zero(rbs * msh.cells_size());

            auto cell_id = 0;
            for(auto cl : msh)
            {
                auto gr  = make_hho_stokes(msh, cl, di, true);
                vector_type svel =  assembler.take_velocity(msh, cl, sol);
                assert((gr.first * svel).rows() == rbs - dim);

                cell_rec_sol.block(cell_id * rbs + dim, 0, rbs - dim, 1) = gr.first * svel;
                cell_rec_sol.block(cell_id * rbs, 0, dim, 1)  = svel.block(0,0, dim, 1);
                cell_sol.block(cell_id * cbs, 0, cbs, 1) = svel.block(0,0, cbs, 1);
                cell_id++;
            }

            compute_discontinuous_velocity( msh, cell_sol, di, "depl2d.msh");

            return;
        }

        bool
        run(const mesh_type& msh, const boundary_type& bnd)
        {
            auto assembler = disk::make_stokes_assembler_alg(msh, di, bnd);
            auto systsz    = assembler.global_system_size();

            sol_old = vector_type::Zero(systsz);
            sol =  vector_type::Zero(systsz);

            Eigen::PardisoLU<Eigen::SparseMatrix<T>>  solver;

            size_t iter = 0;

            make_global_matrix(msh, assembler, iter);
            //WARNINGS: This one must go after make_global_matrix!!!!!!
            make_global_rhs(msh, assembler, iter);

            disk::solvers::pardiso_params<T> pparams;
            mkl_pardiso_ldlt(pparams, assembler.LHS, assembler.RHS, sol);

            post_processing( msh, assembler, iter, true);

            return true;
        }
    };
