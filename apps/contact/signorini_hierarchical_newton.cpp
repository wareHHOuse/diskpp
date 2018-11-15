/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
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

#include <iostream>
#include <sstream>
#include <iomanip>
#include <regex>
#include <type_traits>
#include <cassert>

#include <Eigen/Eigenvalues>

#include "loaders/loader.hpp"
#include "cfem/cfem.hpp"
#include "methods/hho"
#include "mesh/mesh_hierarchy.hpp"

//#include "output/silo.hpp"

#include "contrib/sol2/sol.hpp"
#include "contrib/timecounter.h"
#include "contrib/colormanip.h"

#include "core/output/hdf5_io.hpp"
 #include "signorini_newton_solver.hpp"

template<typename T>
class hierarchical_contact_solver
{
    typedef disk::simplicial_mesh<T, 2>     mesh_type;
    typedef Eigen::SparseMatrix<T>          sparse_matrix_type;
    typedef Eigen::Triplet<T>               triplet_type;
    typedef disk::mechanics::BoundaryConditionsScalar<mesh_type> boundary_type;

    typedef static_vector<T, 3>             fem_vector;
    typedef static_matrix<T, 3, 3>          fem_matrix;

    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>                 hho_vector;
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    hho_matrix;

    std::vector<size_t> fem_compress_map, fem_expand_map;

    T hho_gamma_0, hho_theta, hho_tolerance;
    T fem_gamma_0, fem_theta, fem_tolerance;
    eval_solver_type hho_trace_type;

    hho_degree_info hdi, ref_hdi;
    size_t  hierarchy_levels;
    size_t  sol_level_min, sol_level_max;
    bool compute_ref;

    mesh_type					initial_mesh;
	disk::mesh_hierarchy<T> 	mesh_hier;

    disk::silo_database silo_db;

    dynamic_vector<T> full_sol;
    dynamic_vector<T> e_gx;
    std::vector<dynamic_vector<T>> rec_refs;

    std::vector<std::pair<T,std::pair<T,T>>> errdiams;

    bool
    init(const std::string& config_fn)
    {
        std::cout << "Inside init" << std::endl;
        sol::state lua;
    	lua.open_libraries(sol::lib::math, sol::lib::base);
    	lua["config"] = lua.create_table();
        lua["hs"]  = lua.create_table();
        lua["fem"] = lua.create_table();

    	disk::solvers::init_lua(lua);
    	lua["solver"]["feast"] = lua.create_table();

    	auto r = lua.do_file(config_fn);

    	if ( !r.valid() )
    	{
        	std::cout << "Problems opening configuration file" << std::endl;
        	return false;
    	}

    	std::string     input_mesh  = lua["config"]["input_mesh"];

        std::cout << "input_mesh: " << input_mesh << std::endl;

        size_t hho_degree_cell   = lua["config"]["degree_cell"].get_or(0);
        size_t hho_degree_face   = lua["config"]["degree_face"].get_or(0);

        size_t ref_degree_cell   = lua["config"]["ref_degree_cell"].get_or(0);
        size_t ref_degree_face   = lua["config"]["ref_degree_face"].get_or(0);

        compute_ref       = lua["config"]["compute_reference"].get_or(true);

        hdi = hho_degree_info(hho_degree_cell, hho_degree_face);
        ref_hdi = hho_degree_info(ref_degree_cell, ref_degree_face);

        hierarchy_levels    = lua["hs"]["levels"].get_or(4);
        std::cout << "hierarchy_levels : "<< hierarchy_levels << std::endl;
        sol_level_min       = lua["hs"]["sol_level_min"].get_or(0);
        sol_level_max       = lua["hs"]["sol_level_max"].get_or(0);

        //fem----------------------------------------------------------------------------------
        std::string sfem_theta   = lua["hs"]["fem_theta"];

        if( sfem_theta == "p")     { fem_theta     =  1; }
        else if(sfem_theta == "n") { fem_theta     = -1; }
        else if(sfem_theta == "z") { fem_theta     =  0; }
        else
        {
            std::cout << "Theta must be in {-1,0,1}. Choose only n,z,p (respectively). Falling back to p" << std::endl;
            fem_theta     = 1;
        }

        T EXP_GAMMA    = lua["hs"]["fem_gamma"].get_or(0);
        if (EXP_GAMMA >= 0 && EXP_GAMMA < 10)
        {
            fem_gamma_0 = std::pow(10, EXP_GAMMA);
        }
        else
        {
            std::cout << "EXP_GAMMA : "<< EXP_GAMMA << std::endl;
            std::cout << "Invalid gamma exponent. Falling back to 100" << std::endl;
            fem_gamma_0 = 100;
        }

        std::cout << "fem gamma = "<< fem_gamma_0 << std::endl;
        std::cout << "fem_theta = "<< fem_theta<< std::endl;
        //hho----------------------------------------------------------------------------------
        std::string shho_theta   = lua["hs"]["hho_theta"];

        if( shho_theta == "p")     { hho_theta     =  1; }
        else if(shho_theta == "n") { hho_theta     = -1; }
        else if(shho_theta == "z") { hho_theta     =  0; }
        else
        {
            std::cout << "Theta must be in {-1,0,1}. Choose only n,z,p (respectively). Falling back to p" << std::endl;
            hho_theta     = 1;
        }

        std::string shho_trace   = lua["hs"]["hho_trace"];

        if( shho_trace == "f")     { hho_trace_type =  EVAL_ON_FACES; }
        else if(shho_trace == "l") { hho_trace_type =  EVAL_IN_CELLS_FULL; }
        else
            throw std::invalid_argument("Invalid trace-type. Choose for face-based (f) or cell-based (l) trace versions.");

        auto HEXP_GAMMA    = lua["hs"]["hho_gamma"].get_or(0);
        if (HEXP_GAMMA >= 0 && HEXP_GAMMA < 10)
            hho_gamma_0 = std::pow(10, HEXP_GAMMA);
        else
        {
            std::cout << "HEXP_GAMMA : "<< HEXP_GAMMA << std::endl;
            std::cout << "Invalid gamma exponent. Falling back to 100" << std::endl;
            hho_gamma_0 = 100;
        }

        std::cout << "hho gamma = "<< hho_gamma_0 << std::endl;
        std::cout << "hho_theta = "<< hho_theta<< std::endl;

        //hho----------------------------------------------------------------------------------

        bool verbose = lua["solver"]["feast"]["verbose"].get_or(false);
        auto tol     = lua["solver"]["fem_tolerance"].get_or(0);
        if (tol > 0 && tol < 16)
            fem_tolerance = std::exp(-tol);
        else
        {
            std::cout << "Invalid tolerance. Falling back to tol = 1.e-8" << std::endl;
            fem_tolerance = 1.e-10;
        }


	    if (std::regex_match(input_mesh, std::regex(".*\\.mesh2d$") ))
	    {
	        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;

	        mesh_type msh;
	        disk::netgen_mesh_loader<T, 2> loader;
	        loader.verbose(verbose);
	        if (!loader.read_mesh(input_mesh))
	        {
	            std::cout << "Problem loading mesh." << std::endl;
	            return false;
	        }

	        loader.populate_mesh(msh);

	        mesh_hier.build_hierarchy(msh, hierarchy_levels);

	        std::cout << "Mesh avg. diameter: " << average_diameter(msh) << std::endl;
	    }


	    return true;
	}

    auto
    set_test1(const mesh_type& msh)
    {
        typedef typename mesh_type::point_type  point_type;

        auto f = [](const point_type& pt) -> T {
            return - 2. * M_PI * std::sin( 2. * M_PI * pt.x());
        };

        auto zero_fun = [](const point_type& p) -> T {
            return 0;
        };

        disk::mechanics::BoundaryConditionsScalar<mesh_type> bnd(msh);

        /*--------------------------------------------------------------------------
        *  Check boundary labels for the unitary square domain
        *          Netgen     _____          Medit     _____
        *                4   |     | 2                |     |
        *                    |_____|                  |_____|
        *                       3                        2
        *-------------------------------------------------------------------------*/

        bnd.addDirichletBC(disk::mechanics::DIRICHLET,1, zero_fun); //TOP
        bnd.addNeumannBC(disk::mechanics::NEUMANN, 2,zero_fun); //
        bnd.addNeumannBC(disk::mechanics::NEUMANN, 4,zero_fun); //
        bnd.addContactBC(disk::mechanics::SIGNORINI,3); //BOTTOM


        return std::make_pair(f, bnd);
    }

    auto
    cfem_newton_solver(const mesh_type& msh)
    {
        auto pair = set_test1(msh);
        auto bnd = pair.second;
        auto f = pair.first;

        auto dirichlet_nodes = std::vector<bool>( msh.points_size() );
        for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
        {
            auto fc = *itor;
            auto fc_id = msh.lookup(fc);
            if(bnd.is_dirichlet_face(fc_id))
            {
                auto ptids = fc.point_ids();
                dirichlet_nodes.at(ptids[0]) = true;
                dirichlet_nodes.at(ptids[1]) = true;
            }
        }


        fem_compress_map.resize( msh.points_size() );
        auto system_size = std::count_if(dirichlet_nodes.begin(), dirichlet_nodes.end(), [](bool d) -> bool {return !d;});
        fem_expand_map.resize( system_size );

        auto nnum = 0;
        for (size_t i = 0; i < msh.points_size(); i++)
        {
            if ( dirichlet_nodes.at(i) )
                continue;

            fem_expand_map.at(nnum) = i;
            fem_compress_map.at(i) = nnum++;
        }


        e_gx = dynamic_vector<T>(msh.points_size());

        T maxiter   = 10000;

        for(size_t iter= 0; iter < maxiter ; iter++)
        {
            sparse_matrix_type              gA(system_size, system_size);
            dynamic_vector<T>               gb(system_size), gx(system_size);
            gb = dynamic_vector<T>::Zero(system_size);
            std::vector<triplet_type>       triplets;

            auto is_contact_vector = make_is_contact_vector(msh, bnd);
            for (auto& cl : msh)
            {
                auto ptids = cl.point_ids();
                fem_vector uloc;
                for (size_t i = 0; i < 3; i++)
                    uloc(i) = e_gx(ptids[i]);

                fem_matrix Ah = disk::cfem::stiffness_matrix(msh, cl);
                fem_vector Lh = disk::cfem::make_rhs(msh, cl, f);

                fem_vector Bnegative  = fem_vector::Zero();
                fem_matrix Anitsche   = fem_matrix::Zero();
                fem_matrix Aheaviside = fem_matrix::Zero();

                auto cl_id = msh.lookup(cl);
                if(is_contact_vector.at(cl_id))
                {
                    Anitsche   = make_fem_nitsche(msh, cl, bnd, fem_gamma_0, fem_theta);
                    Bnegative  = make_fem_negative(msh, cl, bnd, fem_gamma_0, fem_theta, uloc);
                    Aheaviside = make_fem_heaviside(msh, cl, bnd, fem_gamma_0, fem_theta, uloc);
                }

                fem_matrix A = Ah - Anitsche + Aheaviside;
                fem_vector b = Lh - (Ah - Anitsche) * uloc - Bnegative;

                for (size_t i = 0; i < A.rows(); i++)
                {
                    if ( dirichlet_nodes.at(ptids[i]) )
                        continue;

                    for (size_t j = 0; j < A.cols(); j++)
                    {
                        if ( dirichlet_nodes.at(ptids[j]) )
                            continue;

                        triplets.push_back( triplet_type(fem_compress_map.at(ptids[i]),
                                                         fem_compress_map.at(ptids[j]),
                                                         A(i,j)) );
                    }

                    gb(fem_compress_map.at(ptids[i])) += b(i);
                }
            }

            gA.setFromTriplets(triplets.begin(), triplets.end());

        #ifdef HAVE_INTEL_MKL
            Eigen::PardisoLU<Eigen::SparseMatrix<T>>  solver;
            //solver.pardisoParameterArray()[59] = 0; //out-of-core
        #else
            Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
        #endif

            size_t systsz = gA.rows();
            size_t nnz = gA.nonZeros();

            //std::cout << "Starting linear solver..." << std::endl;
            //std::cout << " * Solving for " << systsz << " unknowns." << std::endl;
            //std::cout << " * Matrix fill: " << 100.0*double(nnz)/(systsz*systsz) << "%" << std::endl;

            solver.analyzePattern(gA);
            solver.factorize(gA);
            gx = solver.solve(gb);

            dynamic_vector<T> diff_gx = dynamic_vector<T>::Zero(msh.points_size());

            for (size_t i = 0; i < gx.size(); i++)
                diff_gx( fem_expand_map.at(i) ) = gx(i);

            T erroru (0);
            T errord (0);

            for (auto& cl : msh)
            {
                fem_matrix     mass_matrix = disk::cfem::mass_matrix(msh, cl);
                fem_vector     diff_loc, uloc;

                auto ptids = cl.point_ids();
                for (size_t i = 0; i < 3; i++)
                {
                    uloc(i)     = e_gx(ptids[i]);
                    diff_loc(i) = diff_gx(ptids[i]);
                }

                erroru = std::sqrt( uloc.transpose() * mass_matrix * uloc);
                errord = std::sqrt( diff_loc.transpose() * mass_matrix * diff_loc);
            }

            e_gx = e_gx + diff_gx;

            //if(errord/erroru < tolerance)
            if(errord < fem_tolerance)
            {
                std::cout <<  iter << "  " << errord << std::endl;
                break;
            }
        }

        save_data(e_gx, "hier_fem_ref_solution.dat");

        std::ofstream ofs("hier_fem_ref.dat");

        size_t cellnum = 0;
        for (auto& cl : msh)
        {
            auto bar = barycenter(msh, cl);
            auto phi = disk::cfem::eval_basis(msh, cl, bar);
            auto ptids = cl.point_ids();

            double val = 0.0;
            for (size_t i = 0; i < 3; i++)
                val += e_gx(ptids[i]) * phi(i);

            ofs << bar.x() << " " << bar.y() << " " << val << std::endl;
            cellnum++;
        }

        ofs.close();

        return 0;
    }

    auto
    hho_newton_solver(const mesh_type& ref_msh)
    {
        using namespace disk;

        std::cout << green << "*  HHO-REFERENCE  *" << reset << std::endl;
        std::cout << "* Mesh size: " << average_diameter(ref_msh) << std::endl;

        algorithm_parameters<T> ap;
        ap.theta   = hho_theta;
        ap.gamma_0 = hho_gamma_0;
        ap.solver  = hho_trace_type;

        auto pair = set_test1(ref_msh);
        auto bnd = pair.second;
        auto f = pair.first;

        auto zero_fun = [](const typename mesh_type::point_type& p) -> T {
            return 0.;
        };
        dynamic_vector<T> ref_sol;

        std::cout <<  green << "* RUN SIGNORINI" << std::endl;
        size_t num_cell_dofs = scalar_basis_size(ref_hdi.cell_degree(), 2);
        size_t num_face_dofs = scalar_basis_size(ref_hdi.face_degree(), 1);

        auto offset_vector = full_offset(ref_msh, ref_hdi);
        auto is_contact_vector = make_is_contact_vector(ref_msh, bnd);

        if(compute_ref)
        {
            //Solve newton
            switch (ap.solver)
            {
                case EVAL_ON_FACES:
                    ref_sol = solve_faces_hier(ref_msh, f, zero_fun, ap, bnd, ref_hdi);
                    break;
                case EVAL_IN_CELLS_FULL:
                    ref_sol = solve_cells_full_hier(ref_msh, f, zero_fun, ap, bnd, ref_hdi);
                    break;
                default:
                    throw std::invalid_argument("Invalid trace-type. Choose for face-based (f) or cell-based (l) trace versions.");
            }
            std::cout << green <<"* SOLVED" << std::endl;

            save_data(ref_sol, "hier_hho_ref_solution.dat");
        }
        else
            ref_sol = read_data<T>("hier_hho_ref_solution.dat");

        std::cout << green <<"* COMPUTE RECONST" << std::endl;

        size_t cl_count = 0;
        for (auto& cl : ref_msh)
        {
            auto cell_ofs  = offset_vector.at(cl_count);
            auto num_faces = howmany_faces(ref_msh, cl);
            auto num_total_dofs   = num_cell_dofs + num_faces * num_face_dofs;
            hho_vector sol_ufull  = ref_sol.block(cell_ofs, 0, num_total_dofs, 1);

            if (is_contact_vector.at(cl_count) == 1 && ap.solver == EVAL_IN_CELLS_FULL)
            {
                auto gr = make_hho_contact_scalar_laplacian(ref_msh, cl, ref_hdi, bnd);
                rec_refs.push_back( gr.first * sol_ufull );
            }
            else
            {
                auto gr = make_hho_scalar_laplacian( ref_msh, cl, ref_hdi);
                rec_refs.push_back( gr.first * sol_ufull );
            }

            cl_count++;
        }
        std::cout << green << "* DONE" << std::endl;
    }


public:

    hierarchical_contact_solver(const std::string& config_fn)
	{
		init(config_fn);
	}

    void run_fem_hho(size_t level,
                 const disk::simplicial_mesh<T,2>& ref_msh,
                 const disk::simplicial_mesh<T,2>& sol_msh)
    {
        algorithm_parameters<T> ap;
        ap.theta   = hho_theta;
        ap.gamma_0 = hho_gamma_0;

        auto pair = set_test1(sol_msh);
        auto bnd = pair.second;
        auto f = pair.first;

        auto zero_fun = [](const typename disk::simplicial_mesh<T,2>::point_type& p) -> T {
            return 0.;
        };

        //Solve newton
        switch (ap.solver)
        {
            case EVAL_ON_FACES:
                full_sol = solve_faces_hier(sol_msh, f, zero_fun, ap, bnd, hdi);
                break;
            case EVAL_IN_CELLS_FULL:
                full_sol = solve_cells_full_hier(sol_msh, f, zero_fun, ap, bnd, hdi);
                break;
            default:
                throw std::invalid_argument("Invalid trace-type. Choose for face-based (f) or cell-based (l) trace versions.");
        }


        std::vector<T> fem_grad_x, fem_grad_y, hho_grad_x, hho_grad_y;
        fem_grad_x.reserve( ref_msh.cells_size() );
        fem_grad_y.reserve( ref_msh.cells_size() );
        hho_grad_x.reserve( ref_msh.cells_size() );
        hho_grad_y.reserve( ref_msh.cells_size() );

        size_t num_cell_dofs = scalar_basis_size(hdi.cell_degree(), 2);
        size_t num_face_dofs = scalar_basis_size(hdi.face_degree(), 1);

        std::vector<dynamic_vector<T>> rec_sols;

        auto offset_vector = full_offset(sol_msh, hdi);

        save_data(offset_vector, "ofs_ext.data");

        auto is_contact_vector = make_is_contact_vector(sol_msh, bnd);
        size_t sol_cl_count = 0;

        for (auto& sol_cl : sol_msh)
        {
            auto cell_ofs = offset_vector.at(sol_cl_count);
            auto num_total_dofs = num_cell_dofs + 3 * num_face_dofs;
            hho_vector  sol_ufull = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            if (is_contact_vector.at(sol_cl_count) == 1  && ap.solver == EVAL_IN_CELLS_FULL)
            {
                auto gr = make_hho_contact_scalar_laplacian(sol_msh, sol_cl, hdi, bnd);
                rec_sols.push_back( gr.first * sol_ufull );
            }
            else
            {
                auto gr = make_hho_scalar_laplacian(sol_msh, sol_cl, hdi);
                rec_sols.push_back( gr.first * sol_ufull );
            }

            sol_cl_count++;
        }

        T H1_error = 0.0; T L2_error = 0.0; T Linf_error = 0.0;
        size_t cell_i = 0;
        for (auto& ref_cl : ref_msh)
        {
            auto cfem_dphi = disk::cfem::eval_basis_grad(ref_msh, ref_cl);
            auto ptids = ref_cl.point_ids();

            Eigen::Matrix<T,1,2> fem_grad = Eigen::Matrix<T,1,2>::Zero();
            for (size_t i = 0; i < 3; i++)
                fem_grad += e_gx(ptids[i]) * cfem_dphi.block(i,0,1,2);


            fem_grad_x.push_back( fem_grad(0) );
            fem_grad_y.push_back( fem_grad(1) );

            /* find parent cell */
            auto bar = barycenter(ref_msh, ref_cl);
            size_t sol_cl_ofs = mesh_hier.locate_point(bar, level);
            auto sol_cl = *std::next(sol_msh.cells_begin(), sol_cl_ofs);

            /*hho*/

            auto sol_cl_id = sol_msh.lookup(sol_cl);
            auto cell_ofs = offset_vector.at(sol_cl_id);
            auto num_total_dofs = num_cell_dofs + 3 * num_face_dofs;
            hho_vector  sol_ufull = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            dynamic_vector<T> rec_sol = rec_sols.at(sol_cl_ofs);

            auto sol_cb = make_scalar_monomial_basis(sol_msh, sol_cl, hdi.reconstruction_degree());
            auto qps = integrate(ref_msh, ref_cl, 2 * (std::max(hdi.face_degree(), hdi.cell_degree())));
            for (auto& qp : qps)
            {
                auto hho_dphi = sol_cb.eval_gradients(qp.point());

                Eigen::Matrix<T,1,2> hho_grad = Eigen::Matrix<T,1,2>::Zero();
                for (size_t i = 1; i < sol_cb.size(); i++)
                    hho_grad += rec_sol(i-1) * hho_dphi.block(i,0,1,2);

                auto hho_phi = sol_cb.eval_functions(qp.point());

                T hho_ucell = sol_ufull(0);
                for (size_t i = 1; i < sol_cb.size(); i++)
                    hho_ucell += rec_sol(i-1) * hho_phi(i,0);

                Eigen::Matrix<T,1,2> diff = hho_grad - fem_grad;


                T fem_cell = 0.;
                Eigen::Matrix<T,Eigen::Dynamic,1> cfem_phi = disk::cfem::eval_basis(ref_msh, ref_cl, qp.point());
                for (size_t i = 0; i < 3; i++)
                    fem_cell += e_gx(ptids[i]) * cfem_phi(i,0);

                T ucell = hho_ucell - fem_cell;

                H1_error += qp.weight() * diff.dot(diff);
                L2_error += qp.weight() * ucell * ucell;
            }

            {
                auto hho_dphi = sol_cb.eval_gradients(bar);

                Eigen::Matrix<T,1,2> hho_grad = Eigen::Matrix<T,1,2>::Zero();
                for (size_t i = 1; i < sol_cb.size(); i++)
                    hho_grad += rec_sol(i-1) * hho_dphi.block(i,0,1,2);

                hho_grad_x.push_back( hho_grad(0) );
                hho_grad_y.push_back( hho_grad(1) );
            }

            cell_i++;
        }

        disk::silo_zonal_variable<T> rgx("ref_grad_x", fem_grad_x);
        silo_db.add_variable("refmesh", rgx);
        disk::silo_zonal_variable<T> rgy("ref_grad_y", fem_grad_y);
        silo_db.add_variable("refmesh", rgy);

        disk::silo_zonal_variable<T> sgx("sol_grad_x", hho_grad_x);
        silo_db.add_variable("refmesh", sgx);
        disk::silo_zonal_variable<T> sgy("sol_grad_y", hho_grad_y);
        silo_db.add_variable("refmesh", sgy);

        std::cout << "(H1, L2): " << std::sqrt(H1_error);
        std::cout << "  -  " << std::sqrt(L2_error) << std::endl;

        std::ofstream efs("hier_hho__solution_H"+ tostr(level) + ".dat");
        if(!efs.is_open())
            std::cout<< "Error opening file"<<std::endl;

        auto cl_count = 0;
        for(auto& cl : sol_msh)
        {
            auto num_total_dofs = num_cell_dofs + howmany_faces(sol_msh, cl) * num_face_dofs;
            auto cell_ofs = offset_vector.at(cl_count++);
            hho_vector u_bar = full_sol.block(cell_ofs, 0, 1, 1);
            auto bar = barycenter(sol_msh, cl);

            efs << bar.x() << " " << bar.y() <<" "<< u_bar(0) << std::endl;
        }


        efs.close();

    }

    void run_hho_hho(size_t level,
                 const disk::simplicial_mesh<T,2>& ref_msh,
                 const disk::simplicial_mesh<T,2>& sol_msh)
    {
        using namespace disk;

        auto diam = average_diameter(sol_msh);
        std::cout << green << "Computed solution, HHO" << reset << std::endl;
        std::cout << "Mesh size: " << diam << std::endl;

        algorithm_parameters<T> ap;
        ap.theta   = hho_theta;
        ap.gamma_0 = hho_gamma_0;
        ap.solver  = hho_trace_type;

        auto pair = set_test1(sol_msh);
        auto bnd = pair.second;
        auto f = pair.first;

        auto zero_fun = [](const typename disk::simplicial_mesh<T,2>::point_type& p) -> T {
            return 0.;
        };
        //Solve newton
        switch (ap.solver)
        {
            case EVAL_ON_FACES:
                full_sol = solve_faces_hier(sol_msh, f, zero_fun, ap, bnd, hdi);
                break;
            case EVAL_IN_CELLS_FULL:
                full_sol = solve_cells_full_hier(sol_msh, f, zero_fun, ap, bnd, hdi);
                break;
            default:
                throw std::invalid_argument("Invalid trace-type. Choose for face-based (f) or cell-based (l) trace versions.");
        }

        std::vector<T> ref_grad_x, ref_grad_y, hho_grad_x, hho_grad_y;
        ref_grad_x.reserve( ref_msh.cells_size() );
        ref_grad_y.reserve( ref_msh.cells_size() );
        hho_grad_x.reserve( ref_msh.cells_size() );
        hho_grad_y.reserve( ref_msh.cells_size() );

        size_t num_cell_dofs = scalar_basis_size(hdi.cell_degree(), 2);
        size_t num_face_dofs = scalar_basis_size(hdi.face_degree(), 1);

        std::vector<dynamic_vector<T>> rec_sols;

        auto offset_vector = full_offset(sol_msh, hdi);

        auto is_contact_vector = make_is_contact_vector(sol_msh, bnd);

        size_t sol_cl_count = 0;
        for (auto& sol_cl : sol_msh)
        {
            auto cell_ofs  = offset_vector.at(sol_cl_count);
            auto num_faces = howmany_faces(sol_msh, sol_cl);
            auto num_total_dofs   = num_cell_dofs + num_faces * num_face_dofs;
            hho_vector sol_ufull  = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            if (is_contact_vector.at(sol_cl_count) == 1 && ap.solver == EVAL_IN_CELLS_FULL)
            {
                auto gr = make_hho_contact_scalar_laplacian(sol_msh, sol_cl, hdi, bnd);
                rec_sols.push_back( gr.first * sol_ufull );
            }
            else
            {
                auto gr = make_hho_scalar_laplacian(sol_msh, sol_cl, hdi);
                rec_sols.push_back( gr.first * sol_ufull );
            }

            sol_cl_count++;
        }

        T H1_error = 0.0; T L2_error = 0.0; T Linf_error = 0.0;
        size_t cell_i = 0;
        for (auto& ref_cl : ref_msh)
        {
            /* find parent cell */
            auto bar = barycenter(ref_msh, ref_cl);
            size_t ref_cl_ofs = cell_i;
            size_t sol_cl_ofs = mesh_hier.locate_point(bar, level);
            auto sol_cl = *std::next(sol_msh.cells_begin(), sol_cl_ofs);

            dynamic_vector<T> rec_ref = rec_refs.at(ref_cl_ofs);
            dynamic_vector<T> rec_sol = rec_sols.at(sol_cl_ofs);

            auto sol_cb = make_scalar_monomial_basis(sol_msh, sol_cl, hdi.reconstruction_degree());
            auto ref_cb = make_scalar_monomial_basis(ref_msh, ref_cl, ref_hdi.reconstruction_degree());
            auto qps = integrate(ref_msh, ref_cl, 2 * ref_hdi.reconstruction_degree());
            for (auto& qp : qps)
            {
                //1.1. Grad current solution
                auto hho_dphi = sol_cb.eval_gradients(qp.point());

                Eigen::Matrix<T,1,2> hho_grad = Eigen::Matrix<T,1,2>::Zero();
                for (size_t i = 1; i < sol_cb.size(); i++)
                    hho_grad += rec_sol(i-1) * hho_dphi.block(i,0,1,2);

                //1.2. Grad reference
                auto ref_dphi = ref_cb.eval_gradients(qp.point());

                Eigen::Matrix<T,1,2> ref_grad = Eigen::Matrix<T,1,2>::Zero();
                for (size_t i = 1; i < ref_cb.size(); i++)
                    ref_grad += rec_ref(i-1) * ref_dphi.block(i,0,1,2);

                //1.3. H1-error
                Eigen::Matrix<T,1,2> diff = hho_grad - ref_grad;

                //2.1. Current solution
                auto hho_phi = sol_cb.eval_functions(qp.point());
                T hho_ucell = 0;
                for (size_t i = 0; i < sol_cb.size(); i++)
                    hho_ucell += rec_sol(i) * hho_phi(i,0);

                //2.2. Reference solution
                auto ref_phi = ref_cb.eval_functions(qp.point());
                T ref_ucell = 0;
                for (size_t i = 0; i < ref_cb.size(); i++)
                    ref_ucell += rec_ref(i) * ref_phi(i,0);

                //1.2. L2-error
                T ucell = hho_ucell - ref_ucell;

                H1_error += qp.weight() * diff.dot(diff);
                L2_error += qp.weight() * ucell * ucell;
            }

            {
                auto hho_dphi = sol_cb.eval_gradients(bar);

                Eigen::Matrix<T,1,2> hho_grad = Eigen::Matrix<T,1,2>::Zero();
                for (size_t i = 1; i < sol_cb.size(); i++)
                    hho_grad += rec_sol(i-1) * hho_dphi.block(i,0,1,2);

                auto ref_dphi = ref_cb.eval_gradients(bar);

                Eigen::Matrix<T,1,2> ref_grad = Eigen::Matrix<T,1,2>::Zero();
                for (size_t i = 1; i < ref_cb.size(); i++)
                    ref_grad += rec_ref(i-1) * ref_dphi.block(i,0,1,2);

                hho_grad_x.push_back( hho_grad(0) );
                hho_grad_y.push_back( hho_grad(1) );
                ref_grad_x.push_back( ref_grad(0) );
                ref_grad_y.push_back( ref_grad(1) );
            }

            cell_i++;
        }

        disk::silo_zonal_variable<T> rgx("ref_grad_x", ref_grad_x);
        silo_db.add_variable("refmesh", rgx);
        disk::silo_zonal_variable<T> rgy("ref_grad_y", ref_grad_y);
        silo_db.add_variable("refmesh", rgy);

        disk::silo_zonal_variable<T> sgx("sol_grad_x", hho_grad_x);
        silo_db.add_variable("refmesh", sgx);
        disk::silo_zonal_variable<T> sgy("sol_grad_y", hho_grad_y);
        silo_db.add_variable("refmesh", sgy);

        std::cout << std::endl;
        std::cout << "(H1 error, L2 error) : " << std::sqrt(H1_error);
        std::cout << " - " << std::sqrt(L2_error) << std::endl;

        auto error_full = std::make_pair(std::sqrt(H1_error), std::sqrt(L2_error));

        errdiams.push_back( std::make_pair(diam, error_full) );

        std::ofstream efs("hier_hho__solution_H"+ tostr(level) + ".dat");
        if(!efs.is_open())
            std::cout<< "Error opening file"<<std::endl;

        auto cl_count = 0;
        for(auto& cl : sol_msh)
        {
            auto num_faces = howmany_faces(sol_msh, cl);
            auto num_total_dofs   = num_cell_dofs + num_faces * num_face_dofs;
            auto cell_ofs  = offset_vector.at(cl_count++);
            hho_vector u_bar = full_sol.block(cell_ofs, 0, 1, 1);
            auto bar = barycenter(sol_msh, cl);

            efs << bar.x() << " " << bar.y() <<" "<< u_bar(0) << std::endl;
        }


        efs.close();

    }

    void
    verify_convergence()
    {
        bool pass       = true;
        bool warning    = false;
        bool high, low, ok;

        {
            auto error_full_i = errdiams[0].second;

            auto H1_error = error_full_i.first;
            auto L2_error = error_full_i.second;

            std::cout << std::fixed << std::setprecision(3) << errdiams[0].first <<" :    ";
            std::cout << " " <<  std::scientific<< std::setprecision(3)<< H1_error;
            std::cout << "  "<<  std::fixed << std::setprecision(3)<<"        -   ";
            std::cout << " " <<  std::scientific<< std::setprecision(3)<< L2_error <<std::endl;
        }

        T expected_rate = hdi.face_degree() + 1;

        for (size_t i = 1; i < errdiams.size(); i++)
        {

            auto error_full_i = errdiams[i].second;
            auto error_full_i_1 = errdiams[i-1].second;

            auto H1_e = log2(error_full_i_1.first/error_full_i.first);
            auto L2_e = log2(error_full_i_1.second/error_full_i.second);
            auto d = log2(errdiams[i-1].first/errdiams[i].first);

            auto H1_rate    = H1_e/d;
            auto L2_rate    = L2_e/d;

            ok   = (std::abs(expected_rate - H1_rate) < 0.4); /* Test passed */
            low  = ((expected_rate - H1_rate) > 0.2); /* a bit too low, warn */
            high = ((H1_rate - expected_rate) > 0.2); /* a bit too high, warn */

            if (low)    std::cout << magenta;
            if (high)   std::cout << cyan;


            auto H1_error = error_full_i.first;
            auto L2_error = error_full_i.second;

            std::cout << std::fixed << std::setprecision(3) << errdiams[i].first <<" :    ";
            std::cout << " " <<  std::scientific<< std::setprecision(3)<< H1_error;
            std::cout << "  "<<  std::fixed << std::setprecision(3)<< H1_rate << "   -   ";
            std::cout << " " <<  std::scientific<< std::setprecision(3)<< L2_error;
            std::cout << "  "<<  std::fixed << std::setprecision(3)<< L2_rate <<std::endl;
            if (low or high)
            {
                std::cout << reset;
                warning = true;
            }
        }

        std::string             okfail = "[\x1b[31;1mFAIL\x1b[0m]";
        if (ok && not warning)  okfail = "[\x1b[32;1m OK \x1b[0m]";
        if (ok && warning)      okfail = "[\x1b[33;1mWARN\x1b[0m]";

        std::cout << okfail << std::endl;

    }

    void
    run_computations()
    {
        auto ref_msh = mesh_hier.finest_mesh();

        std::string visit_output_filename = "hierarchical.silo";
        silo_db.create(visit_output_filename);
        std::cout << green << "Reference solution, HHO" << reset << std::endl;
        std::cout << "Mesh size: " << average_diameter(ref_msh) << std::endl;

        //cfem_newton_solver(ref_msh);
        //e_gx = read_data<T>("hier_fem_ref_solution.dat");

        hho_newton_solver(ref_msh);

        for (size_t level = sol_level_min; level <= sol_level_max; level++)
        {
            std::cout << "* Level : "<< level << std::endl;
            auto sol_msh = *std::next(mesh_hier.meshes_begin(), level);
            //run_fem_hho(level, ref_msh, sol_msh);
            run_hho_hho(level, ref_msh, sol_msh);
        }
        silo_db.close();

        verify_convergence();
    }
};



int main(int argc, char **argv)
{
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    using scalar_type = double;

    if (argc != 2)
    {
        std::cout << "Please specify configuration file" << std::endl;
        return 1;
    }

    hierarchical_contact_solver<scalar_type> hcs(argv[1]);

    hcs.run_computations();

}
