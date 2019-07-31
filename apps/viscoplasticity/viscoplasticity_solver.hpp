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

enum problem_type
{
    DRIVEN,
    COUETTE,
    POISEUILLE,
    VANE
};

template<typename Mesh>
class augmented_lagrangian_viscoplasticity
{
    typedef Mesh mesh_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;
    typedef typename mesh_type::coordinate_type T;

    using point_type = typename mesh_type::point_type;

    typedef Matrix<T, Mesh::dimension, Mesh::dimension>         tensor_type;
    typedef Matrix<T, Mesh::dimension, 1>                       vector2d_type;
    typedef Matrix<T, Dynamic, 1>                               vector_type ;
    typedef Matrix<T, Dynamic, Dynamic>                         matrix_type;

    typedef std::function<vector2d_type (const point_type &)>   vector_funtion_type;
    typedef std::function<T   (const point_type &)>             scalar_funtion_type;

    typedef disk::vector_boundary_conditions<mesh_type>          boundary_type;

    vector_funtion_type     rhs_fun, velocity;
    dynamic_vector<T>       multiplier, auxiliar, auxiliar_old;

    typename disk::hho_degree_info di;
    T             factor;
    T             viscosity;
    T             alpha;
    T             yield;
    size_t        cbs, fbs, pbs, sbs, dim;

public:
    dynamic_vector<T>       sol, sol_old;
    std::tuple<T, T, T>     convergence;
    bool                    use_sym_grad;

    augmented_lagrangian_viscoplasticity(const mesh_type& msh,
                            const typename disk::hho_degree_info & hdi,
                            const T& alpha_ext):
                            di(hdi), alpha(alpha_ext)
    {
        use_sym_grad = true;
        factor = (use_sym_grad)? 2. : 1.;
        T omegaExt = 2.;
        T f = 1;
        //T Bn  =  0.; //std::sqrt(2) * 10.;
        //Driven
        //yield =  Bn * f * Lref;//Bn * viscosity;

        //Couette
        //yield =  Bn * omegaExt;

        //VANE
        T Lref = 1.; //R
        viscosity = 1.;
        T omega = 1;
        T Vref = omega * Lref;
        yield = Bn * (viscosity * Vref) / Lref; // Bn/std::sqrt(2); // ;

        dim =  Mesh::dimension;

        cbs = disk::vector_basis_size(di.cell_degree(), dim, dim);
        fbs = disk::vector_basis_size(di.face_degree(), dim - 1, dim);
        pbs = disk::scalar_basis_size(di.face_degree(), dim);
        sbs = disk::sym_matrix_basis_size(di.face_degree(), dim, dim);
    };

    auto
    define_problem(const mesh_type& msh, const problem_type& problem )
    {
        boundary_type bnd(msh);

        auto wall = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
            return Matrix<T, Mesh::dimension, 1>::Zero();
        };
        auto movingWall = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
            return Matrix<T, Mesh::dimension, 1>{1,0};
        };
        auto symmetryPlane = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
            return Matrix<T, Mesh::dimension, 1>{0,0};
        };

        switch (problem)
		{
            case DRIVEN:
                velocity  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    if( std::abs(p.y() - 1.) < 1.e-8 )
                        return Matrix<T, Mesh::dimension, 1>{1,0};
                    else
                        return Matrix<T, Mesh::dimension, 1>{0,0};
                };
                rhs_fun  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    return Matrix<T, Mesh::dimension, 1>::Zero();
                };

                #if 0
                bnd.addDirichletBC(0, 1, movingWall);
                bnd.addDirichletBC(0, 2, wall);
                bnd.addDirichletBC(0, 3, wall);
                bnd.addDirichletBC(0, 4, wall);
                #endif
                bnd.addDirichletEverywhere(velocity);

               break;


            case COUETTE:

                rhs_fun  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    return Matrix<T, Mesh::dimension, 1>::Zero();
                };

                velocity  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    Matrix<T, Mesh::dimension, 1> ret = Matrix<T, Mesh::dimension, 1>::Zero();

                    auto theta  = std::atan2(p.y() , p.x());

                    T Rint(0.5), Rext(1.), omegaInt(2.),omegaExt(2.);
                    T r = std::sqrt(p.x()*p.x() + p.y()*p.y());

                    //1.All Solid
                    //auto u_theta_solid = (omegaExt/ omegaInt) * (r / Rint);
                    #if 0
                    //2.All Liquid
                    auto eta = Rint/Rext;
                    auto sigma_i = (2./(eta * eta - 1.)) * ( (1. - (omegaExt*eta / omegaInt))*((1 - eta)/eta)
                        + sgn((omegaInt*eta / omegaExt) - 1.) * Bn * std::log(eta) );
                    auto term_1 = u_theta_solid;
                    auto term_2 = 0.5 *(Rint * Rint) * sigma_i * r * ( 1./(Rext* Rext) - 1/(r *r)) ;
                    auto term_3 = Bn * r * std::log(Rext/r) * sgn(sigma_i);
                    u_theta = term_1 + term_2 + term_3;
                    //3. with plug
                    #endif

                    //ret(0) =  -std::sin(theta) * 2.* r;
                    //ret(1) =  std::cos(theta) * 2.* r;
                    return ret;
                };

                auto v4  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    Matrix<T, Mesh::dimension, 1> ret = Matrix<T, Mesh::dimension, 1>::Zero();

                    auto theta  = std::atan2(p.y() , p.x());
                    //T r = std::sqrt(p.x()*p.x() + p.y()*p.y());
                    //w1 = 1; r = 0.5
                    ret(0) =  -std::sin(theta) * 0.5 ;
                    ret(1) =  std::cos(theta) * 0.5;

                    return ret;
                };
                auto v3  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    Matrix<T, Mesh::dimension, 1> ret = Matrix<T, Mesh::dimension, 1>::Zero();

                    auto theta  = std::atan2(p.y() , p.x());
                    //T r = std::sqrt(p.x()*p.x() + p.y()*p.y());
                    //w1 = 2; r = 1
                    ret(0) =  -std::sin(theta);// * 2. ;
                    ret(1) =  std::cos(theta);// * 2.;

                    //std::cout << "(x,y) = ("<< p.x() << ";"<< p.y() <<");  theta = "  << theta;
                    //std::cout << "  vel = ("<< ret(0) << "; "<< ret(1) << ")"  << std::endl;

                    return ret;
                };

                bnd.addNeumannBC(10, 1, symmetryPlane);
                bnd.addNeumannBC(10, 2, symmetryPlane);

                bnd.addDirichletBC(0, 3, v3); //velocity);
                bnd.addDirichletBC(0, 4, v4);//velocity);

                break;

            case VANE:
                std::cout << " I'm in VANE" << std::endl;
                rhs_fun  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    return Matrix<T, Mesh::dimension, 1>::Zero();
                };

                auto rot  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    T omega = 1;
                    return Matrix<T, Mesh::dimension, 1>{omega * p.y(), -omega * p.x()};
                };

                bnd.addDirichletBC( 0, 2, wall);
                bnd.addDirichletBC( 0, 1, rot);
                bnd.addNeumannBC(10, 3, symmetryPlane);

            default:
                throw std::invalid_argument("No problem defined");
        }

        auto assembler = disk::make_stokes_assembler_alg(msh, di, bnd);

        return assembler;
    }

    template<typename Assembler>
    auto
    initialize(const mesh_type& msh, const Assembler& assembler)
    {
        auto systsz = assembler.global_system_size();
        sol = dynamic_vector<T>::Zero(systsz);
        sol_old = dynamic_vector<T>::Zero(systsz);

        multiplier = dynamic_vector<T>::Zero(msh.cells_size() * sbs);
        auxiliar   = dynamic_vector<T>::Zero(msh.cells_size() * sbs);
        auxiliar_old = dynamic_vector<T>::Zero(msh.cells_size() * sbs);

        return;
    }

    template<typename Assembler>
    auto
    compute_errors( const mesh_type& msh,
                    const Assembler& assembler,
                    const bool& alg_finished)
    {
        T error_vel  = 0.;
        T error_temp = 0.;

        std::ofstream ofs;

        if(alg_finished)
        {
            ofs = std::ofstream("Gu_norm.data");
            if (!ofs.is_open())
                std::cout << "Error opening errors "<<std::endl;
        }

        for (auto& cl : msh)
        {
        	auto bar = barycenter(msh, cl);

            //energy error
            vector_type svel = assembler.take_velocity(msh, cl, sol);
            vector_type pvel = disk::project_function(msh, cl, di, velocity);
            //vector_type  svel_old =  assembler.take_velocity(msh, cl, sol_old);
            vector_type diff_vel = svel - pvel;
            auto gr = disk::make_hho_stokes(msh, cl, di, use_sym_grad);
            matrix_type stab;
            stab = make_vector_hho_stabilization(msh, cl, gr.first, di);
            auto G = disk::make_hlow_stokes(msh, cl, di, use_sym_grad);

            matrix_type B = factor * (viscosity*G.second + viscosity*stab);

            error_vel += diff_vel.dot(B * diff_vel);
            auto Gu_norm = svel.dot(B * svel);
            error_temp += Gu_norm;

            if(alg_finished)
                ofs << bar.x()<< " " << bar.y() << " " << Gu_norm<< std::endl;
        }
        if(alg_finished)
            ofs.close();
        return std::make_pair(std::sqrt(error_vel), std::sqrt(error_temp));
    }

    template<typename Assembler>
    Matrix<T, Dynamic, 1>
    compute_auxiliar(   const mesh_type& msh,
                        const cell_type& cl,
                        const Assembler& assembler,
                        const dynamic_vector<T>& velocity_dofs)
    {
        vector_type u_TF  = assembler.take_velocity(msh, cl, velocity_dofs);
        auto value = 1./(factor * (viscosity + alpha));
        auto G = disk::make_hlow_stokes(msh, cl, di, use_sym_grad);

        auto cell_ofs    = disk::priv::offset(msh, cl);
        vector_type Gu = G.first * u_TF;
        vector_type stress = multiplier.block(cell_ofs * sbs, 0, sbs, 1);

        //Theta
        auto sb = disk::make_sym_matrix_monomial_basis(msh, cl, di.face_degree());
        //barycenter only for k = 0; fix this for higher orders
        auto bar = barycenter(msh, cl);
        auto s_phi  = sb.eval_functions(bar);

        vector_type  theta  = stress  +  factor * alpha * Gu;
        Matrix<T, Mesh::dimension, Mesh::dimension> theta_eval = disk::eval(theta, s_phi);
        T theta_norm  = std::sqrt((theta_eval.cwiseProduct(theta_eval)).sum());

        T theta_eigen = theta_eval.norm();
        assert(theta_norm == theta_eigen);

        //#if 0
        //Gamma Driven
        // A. Solid
        if(theta_norm <=  std::sqrt(2) * yield ||  std::abs(theta_norm - std::sqrt(2) * yield) < 1.e-8)
            return  Matrix<T, Dynamic, 1>::Zero(sbs);
        else  // B. liquid
            return  value * theta * (1. - std::sqrt(2) *(yield/theta_norm));
        //#endif
        #if 0
        // A. Solid
        if(theta_norm <=   yield ||  std::abs(theta_norm -  yield) < 1.e-8)
            return  Matrix<T, Dynamic, 1>::Zero(sbs);
        else  // B. liquid
            return  value * theta * (1. - (yield/theta_norm));
        #endif
    }

    template<typename Assembler>
    auto
    update_multiplier(const mesh_type& msh, const Assembler& assembler)
    {
        T conv_stress = 0.;
        T conv_gamma = 0.;

        const auto dim = Mesh::dimension;

        for(auto cl: msh)
        {
            vector_type u_TF = assembler.take_velocity(msh, cl, sol);
            auto G = disk::make_hlow_stokes(msh, cl, di, use_sym_grad);
            auto cell_ofs = disk::priv::offset(msh, cl);

            vector_type Gu = G.first * u_TF;
            vector_type gamma_old = auxiliar_old.block(cell_ofs *sbs, 0, sbs, 1);
            vector_type gamma = auxiliar.block(cell_ofs *sbs, 0, sbs, 1);

            vector_type diff_stress = factor * alpha * (Gu - gamma);
            vector_type diff_gamma  = alpha * (gamma - gamma_old);

            multiplier.block(cell_ofs * sbs, 0, sbs, 1) += diff_stress;

            auto sb = disk::make_sym_matrix_monomial_basis(msh, cl, di.face_degree());
            matrix_type mass = make_mass_matrix(msh, cl, sb);

            conv_stress += diff_stress.dot(mass * diff_stress);
            conv_gamma  += diff_gamma.dot(mass * diff_gamma);
        }

        convergence = std::make_tuple(conv_stress + conv_gamma, conv_stress, conv_gamma);

        auxiliar_old = auxiliar;
        return;
    }

    template<typename Assembler>
    Matrix<T, Dynamic, 1>
    make_rhs_alg(   const mesh_type& msh,
                    const cell_type& cl,
                    const Assembler& assembler)
    {
        auto G = disk::make_hlow_stokes(msh, cl, di, use_sym_grad);
        auto cb = disk::make_vector_monomial_basis(msh, cl, di.cell_degree());
        auto sb = disk::make_sym_matrix_monomial_basis(msh, cl, di.face_degree());

        auto cell_ofs =  disk::priv::offset(msh, cl);
        auto num_faces = howmany_faces(msh, cl);

        vector_type rhs = vector_type::Zero(cbs + fbs * num_faces);

        //(f, v_T)
        rhs.block( 0, 0, cbs, 1) = make_rhs(msh, cl, cb, rhs_fun);

        //(stress - alpha * gamma, Gv)
        vector_type stress = multiplier.block( sbs * cell_ofs,  0, sbs, 1);
        vector_type gamma  = compute_auxiliar( msh,  cl, assembler, sol_old); //or sol at this point it's the same
        auxiliar.block(cell_ofs * sbs, 0, sbs, 1) = gamma;
        vector_type str_agam = stress - factor * alpha * gamma;
        matrix_type mm = disk::make_mass_matrix(msh, cl, sb);

        rhs -=  G.first.transpose() * mm * str_agam;

        return rhs;
    }

    template<typename Assembler>
    auto
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
    auto
    make_global_matrix(const mesh_type& msh, Assembler& assembler)
    {
        assembler.initialize_lhs();

        for (auto cl : msh)
        {
            auto G  = disk::make_hlow_stokes(msh, cl, di, use_sym_grad);
            auto gr = disk::make_hho_stokes(msh, cl, di, use_sym_grad);
            matrix_type stab = make_vector_hho_stabilization(msh, cl, gr.first, di);
            auto dr = make_hho_divergence_reconstruction_rhs(msh, cl, di);

            matrix_type A = factor *(alpha * G.second + viscosity * stab);

            assembler.assemble_lhs(msh, cl, A, -dr);
        }

        assembler.finalize_lhs();

        return;
    }

    template<typename Assembler>
    auto
    run_stokes_like(const mesh_type& msh, Assembler& assembler, const size_t iter)
    {
        sol_old = sol;

        make_global_rhs(msh,assembler);

        if(iter == 0)
            make_global_matrix(msh, assembler);

        //dump_sparse_matrix(assembler.LHS, "stokes.txt");
        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        sol = dynamic_vector<T>::Zero(systsz);
        disk::solvers::pardiso_params<T> pparams;
        mkl_pardiso_ldlt(pparams, assembler.LHS, assembler.RHS, sol);

        return;
    }

    template<typename Assembler>
    void
    post_processing(const mesh_type& msh, Assembler& assembler,
                    const std::string & info,
                    const problem_type& problem)
    {
        auto dim = Mesh::dimension;
        auto rbs = disk::vector_basis_size(di.reconstruction_degree(), dim, dim);

        dynamic_vector<T> cell_sol(cbs * msh.cells_size());
        dynamic_vector<T> cell_rec_sol(rbs * msh.cells_size());
        dynamic_vector<T> press_vec(pbs * msh.cells_size());

        std::ofstream ofs("data_" + info + ".data");
        if (!ofs.is_open())
            std::cout << "Error opening file"<<std::endl;

        for(auto cl : msh)
        {
            auto gr  = disk::make_hho_stokes(msh, cl, di, use_sym_grad);
            auto cell_ofs = disk::priv::offset(msh, cl);
            vector_type svel =  assembler.take_velocity(msh, cl, sol);
            assert((gr.first * svel).rows() == rbs - dim);
            cell_rec_sol.block(cell_ofs * rbs + dim, 0, rbs - dim, 1) = gr.first * svel;
            cell_rec_sol.block(cell_ofs * rbs, 0, dim, 1) = svel.block(0,0, dim, 1);
            cell_sol.block(cell_ofs * cbs, 0, cbs, 1) = svel.block(0,0, cbs, 1);

            //this is only for k = 0, since there is only one velocity;
            auto bar = barycenter(msh, cl);

            //Velocity
            vector_type cell_vel = svel.block(0,0, cbs, 1);
            auto cb  = disk::make_vector_monomial_basis(msh, cl, di.cell_degree());
            auto phi = cb.eval_functions(bar);
            Matrix<T, Mesh::dimension, 1> ueval = disk::eval(cell_vel, phi);

            //Pressure
            vector_type spress =  assembler.take_pressure(msh, cl, sol);
            auto pb  = disk::make_scalar_monomial_basis(msh, cl, di.face_degree());
            auto p_phi = pb.eval_functions(bar);
            T peval =  p_phi.dot(spress);

            //Stress
            auto sb = disk::make_sym_matrix_monomial_basis(msh, cl, di.face_degree());
            auto s_phi  = sb.eval_functions(bar);
            auto G = disk::make_hlow_stokes(msh, cl, di, use_sym_grad);
            vector_type Gu = G.first * svel;
            vector_type stress = multiplier.block(cell_ofs * sbs, 0, sbs, 1);
            vector_type theta  = stress  +  factor * alpha * Gu;
            tensor_type theta_eval = disk::eval(theta, s_phi);
            tensor_type sigma_eval = disk::eval(stress, s_phi);
            tensor_type grad_eval = disk::eval(Gu, s_phi);

            T divu = grad_eval(0,0) + grad_eval(1,1);
            T tr_stress = sigma_eval(0,0) + sigma_eval(1,1);
            ofs << ueval(0)   << " " << ueval(1) << " " << peval<< " ";
            ofs << theta_eval.norm() << " " << sigma_eval.norm()   << " ";
            ofs << divu << " "<< tr_stress<<std::endl;
        }
        ofs.close();

        typedef point<T,2>             point_type;
        std::pair<point_type, point_type> p_x, p_y;
        auto eps = 1.e-4;
        switch(problem)
        {
            case COUETTE:
                p_x = std::make_pair(point_type({0., 0.}), point_type({0.866, 0.5}));
                p_y = std::make_pair(point_type({0.0, 0.0}), point_type({0.0, 0.0}));
                break;
            default:
                p_x = std::make_pair(point_type({0.0 + eps, 0.5 + eps}), point_type({1.0 + eps, 0.5 + eps}));
                p_y = std::make_pair(point_type({0.5 + eps, 0.0 + eps}), point_type({0.5 + eps, 1.0 + eps}));
                break;
        }
        plot_over_line(msh, p_x, cell_rec_sol, di.reconstruction_degree(), "plot_over_x_" + info + ".data");
        //plot_over_line(msh, p_y, cell_rec_sol, di.reconstruction_degree(), "plot_over_y_" + info + ".data");
        plot_over_line(msh, p_y, cell_sol, di.cell_degree(), "plot_over_y_" + info + ".data");
        compute_discontinuous_velocity( msh, cell_sol, di, "velocity_" + info +".msh");
        save_coords(msh, "Coords_"+ info + ".data");
        quiver( msh, sol, assembler, di, "quiver_"+ info + ".data");
        return;
    }

};
