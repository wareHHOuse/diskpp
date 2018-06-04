/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
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
#include "hho/hho.hpp"
#include "revolution/methods/hho"
#include "mesh/mesh_hierarchy.hpp"

#include "output/silo.hpp"
#include "solvers/solver.hpp"

#include "contrib/sol2/sol.hpp"
#include "contrib/timecounter.h"
#include "contrib/colormanip.h"

template<typename T>
class hierarchical_eigval_solver
{
	typedef disk::simplicial_mesh<T, 2>  mesh_type;

	mesh_type					initial_mesh;
	disk::mesh_hierarchy<T> 	mesh_hier;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    fem_eigvecs;
    Eigen::Matrix<T, Eigen::Dynamic, 1>                 fem_eigvals;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    hho_eigvecs;
    Eigen::Matrix<T, Eigen::Dynamic, 1>                 hho_eigvals;

    feast_eigensolver_params<T> fep;

    disk::silo_database silo_db;

    std::vector<size_t> fem_compress_map, fem_expand_map;
    std::vector<size_t> hho_compress_map, hho_expand_map;

    T       stab_weight;
    size_t  hho_degree;
    size_t  hierarchy_levels;
    size_t  which_eigvec;
    size_t  which_level;

    bool
    cfem_eigenvalue_solver(const disk::simplicial_mesh<T, 2>& msh)
    {
        typedef Eigen::SparseMatrix<T>  sparse_matrix_type;
        typedef Eigen::Triplet<T>       triplet_type;

        std::vector<bool> dirichlet_nodes( msh.points_size() );
        for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
        {
            auto fc = *itor;
            auto ptids = fc.point_ids();
            dirichlet_nodes.at(ptids[0]) = true;
            dirichlet_nodes.at(ptids[1]) = true;
        }


        fem_compress_map.resize( msh.points_size() );
        size_t system_size = std::count_if(dirichlet_nodes.begin(), dirichlet_nodes.end(), [](bool d) -> bool {return !d;});
        fem_expand_map.resize( system_size );

        auto nnum = 0;
        for (size_t i = 0; i < msh.points_size(); i++)
        {
            if ( dirichlet_nodes.at(i) )
                continue;

            fem_expand_map.at(nnum) = i;
            fem_compress_map.at(i) = nnum++;
        }

        sparse_matrix_type  gK(system_size, system_size);
        sparse_matrix_type  gM(system_size, system_size);

        std::vector<triplet_type> gtK, gtM;

        timecounter tc;
        tc.tic();
        for (auto& cl : msh)
        {
            auto bar = barycenter(msh, cl);


            auto K = disk::cfem::stiffness_matrix(msh, cl);
            auto M = disk::cfem::mass_matrix(msh, cl);


            auto ptids = cl.point_ids();

            for (size_t i = 0; i < K.rows(); i++)
            {
                if ( dirichlet_nodes.at(ptids[i]) )
                    continue;

                for (size_t j = 0; j < K.cols(); j++)
                {
                    if ( dirichlet_nodes.at(ptids[j]) )
                        continue;

                    gtK.push_back( triplet_type(fem_compress_map.at(ptids[i]),
                                                fem_compress_map.at(ptids[j]),
                                                K(i,j)) );

                    gtM.push_back( triplet_type(fem_compress_map.at(ptids[i]),
                                                fem_compress_map.at(ptids[j]),
                                                M(i,j)) );
                }
            }
        }

        gK.setFromTriplets(gtK.begin(), gtK.end());
        gM.setFromTriplets(gtM.begin(), gtM.end());

        tc.toc();
        std::cout << "Assembly time: " << tc << " seconds (" << msh.cells_size() << " elems)" << std::endl;

        tc.tic();
        bool solver_ok = generalized_eigenvalue_solver(fep, gK, gM, fem_eigvecs, fem_eigvals);
        tc.toc();
        if (!solver_ok)
        {
            std::cout << "Problem during solution phase" << std::endl;
            return false;
        }
        std::cout << "Eigensolver time: " << tc << " seconds (" << gK.rows() << " DoFs)" << std::endl;

        /* Output solutions via silo */
        silo_db.add_mesh(msh, "refmesh");
        for (size_t i = 0; i < fep.eigvals_found; i++)
        {
            std::vector<T> solution_vals;
            solution_vals.resize(msh.points_size());

            dynamic_vector<T> gx = fem_eigvecs.block(0,i,fem_eigvecs.rows(),1);
            for (size_t i = 0; i < gx.size(); i++)
                solution_vals.at( fem_expand_map.at(i) ) = gx(i);

            std::stringstream ss;
            ss << "u_ref" << i;

            disk::silo_nodal_variable<T> u(ss.str(), solution_vals);
            silo_db.add_variable("refmesh", u);
        }

        /* Dump actual eigenvalues */
        std::string fname = "ref_eigenvalues.txt";
        std::ofstream ofs(fname);
        for (size_t i = 0; i < fep.eigvals_found; i++)
        {
            ofs << std::setprecision(10) << fem_eigvals(i) << " ";

            dynamic_vector<T> gx = fem_eigvecs.block(0,i,fem_eigvecs.rows(),1);
            ofs << gx.dot(gM*gx) << std::endl;
        }
        ofs.close();

    }

    template<typename Mesh>
    bool
    hho_eigenvalue_solver(const Mesh& msh, size_t degree)
    {
        using namespace revolution;

        typedef Eigen::Triplet<T> triplet_type;

        std::cout << "Assembly" << std::endl;

        size_t num_cell_dofs = scalar_basis_size(degree, Mesh::dimension);
        size_t num_face_dofs = scalar_basis_size(degree, Mesh::dimension-1);
        size_t num_global_cell_dofs = msh.cells_size() * num_cell_dofs;
        size_t num_global_face_dofs = msh.internal_faces_size() * num_face_dofs;
        size_t num_global_dofs = num_global_cell_dofs + num_global_face_dofs;

        /* Strong imposition of dirichlet bc */
        hho_compress_map.resize( msh.faces_size() );
        size_t num_non_dirichlet_faces = msh.internal_faces_size();
        hho_expand_map.resize( num_non_dirichlet_faces );

        size_t fnum = 0;
        size_t face_i = 0;
        for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
        {
            if ( msh.is_boundary(*itor) )
            {
                face_i++;
                continue;
            }

            hho_expand_map.at(fnum) = face_i;
            hho_compress_map.at(face_i) = fnum++;
            face_i++;
        }

        sparse_matrix<T> gK(num_global_dofs, num_global_dofs);
        sparse_matrix<T> gM(num_global_dofs, num_global_dofs);

        std::vector<triplet_type> gtK, gtM;

        hho_degree_info hdi(degree);

        size_t elem_i = 0;
        for (auto& cl : msh)
        {
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            dynamic_matrix<T> stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
            dynamic_matrix<T> mass   = make_hho_eigval_mass_matrix(msh, cl, hdi.cell_degree());

            dynamic_matrix<T> K = gr.second + stab_weight * stab;
            dynamic_matrix<T> M = mass;

            auto fcs = faces(msh, cl);
            auto num_faces = fcs.size();

            std::vector<size_t> offset_table;
            offset_table.resize(num_cell_dofs + num_face_dofs*num_faces);

            for (size_t i = 0; i < num_cell_dofs; i++)
                offset_table.at(i) = elem_i*num_cell_dofs + i;

            for (size_t i = 0; i < num_faces; i++)
            {
                auto fc = fcs[i];

                if ( msh.is_boundary(fc) )
                {
                    for (size_t j = 0; j < num_face_dofs; j++)
                        offset_table.at(num_cell_dofs + i*num_face_dofs + j) = 0xdeadbeef;
                }
                else
                {
                    auto face_id = msh.lookup(fc);
                    auto face_pos = hho_compress_map.at(face_id);
                    //std::cout << face_id << " -> " << face_pos << std::endl;
                    for (size_t j = 0; j < num_face_dofs; j++)
                        offset_table.at(num_cell_dofs + i*num_face_dofs + j) = num_cell_dofs*msh.cells_size() + num_face_dofs*face_pos+j;
                }
            }

            /* K has cell and face parts */
            for (size_t i = 0; i < K.rows(); i++)
            {
                size_t i_ofs = offset_table.at(i);
                if (i_ofs == 0xdeadbeef)
                    continue;

                for (size_t j = 0; j < K.cols(); j++)
                {
                    size_t j_ofs = offset_table.at(j);
                    if (j_ofs == 0xdeadbeef)
                        continue;

                    assert(i_ofs < num_global_dofs);
                    assert(j_ofs < num_global_dofs);
                    gtK.push_back( triplet_type(i_ofs, j_ofs, K(i,j)) );
                }
            }

            /* M has only cell part */
            for (size_t i = 0; i < M.rows(); i++)
            {
                for (size_t j = 0; j < M.cols(); j++)
                {
                    size_t i_ofs = offset_table.at(i);
                    size_t j_ofs = offset_table.at(j);
                    assert(i_ofs < num_global_dofs);
                    assert(j_ofs < num_global_dofs);
                    gtM.push_back( triplet_type(i_ofs, j_ofs, M(i,j)) );
                }
            }

            elem_i++;
        }

        gK.setFromTriplets(gtK.begin(), gtK.end());
        gM.setFromTriplets(gtM.begin(), gtM.end());

        generalized_eigenvalue_solver(fep, gK, gM, hho_eigvecs, hho_eigvals);

        silo_db.add_mesh(msh, "solmesh");
        for (size_t i = 0; i < fep.eigvals_found; i++)
        {
            std::vector<T> solution_vals;
            solution_vals.resize(msh.cells_size());

            dynamic_vector<T> gx = hho_eigvecs.block(0,i,hho_eigvecs.rows(),1);
            for (size_t i = 0; i < msh.cells_size(); i++)
                solution_vals.at(i) = gx(i*num_cell_dofs);

            std::stringstream ss;
            ss << "u" << i;

            disk::silo_zonal_variable<T> u(ss.str(), solution_vals);
            silo_db.add_variable("solmesh", u);
        }

        std::string fname = "sol_eigenvalues.txt";
        std::ofstream ofs(fname);
        for (size_t i = 0; i < fep.eigvals_found; i++)
        {
            ofs << std::setprecision(10) << hho_eigvals(i) << " ";

            dynamic_vector<T> gx = hho_eigvecs.block(0,i,hho_eigvecs.rows(),1);
            ofs << gx.dot(gM*gx) << std::endl;
        }
        ofs.close();
    }


	bool init(const std::string& config_fn)
	{
		sol::state lua;
    	lua.open_libraries(sol::lib::math, sol::lib::base);
    	lua["config"] = lua.create_table();
        lua["hs"] = lua.create_table();
    	disk::solvers::init_lua(lua);
    	lua["solver"]["feast"] = lua.create_table();

    	auto r = lua.do_file(config_fn);

    	if ( !r.valid() )
    	{
        	std::cout << "Problems opening configuration file" << std::endl;
        	return false;
    	}

    	std::string     input_mesh  = lua["config"]["input_mesh"];
    	bool            verbose     = lua["config"]["verbose"].get_or(false);

        hho_degree          = lua["config"]["degree"].get_or(0);

        hierarchy_levels    = lua["hs"]["levels"].get_or(4);
        which_eigvec        = lua["hs"]["eigvec"].get_or(0);
        which_level         = lua["hs"]["sol_level"].get_or(0);
        stab_weight         = lua["hs"]["stab_weight"].get_or(1);



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

        if ( !setup_feast(lua, fep) )
        {
            std::cout << "Problems with parameters. Check configuration." << std::endl;
            return false;
        }

	    return true;
	}

public:
	hierarchical_eigval_solver()
	{}

	hierarchical_eigval_solver(const std::string& config_fn)
	{
		init(config_fn);
	}

	void
	load_config(const std::string& config_fn)
	{
		init(config_fn);
	}

    void
    run_computations()
    {
        using namespace revolution;

        auto ref_msh = mesh_hier.finest_mesh();
        auto sol_msh = *std::next(mesh_hier.meshes_begin(), which_level);

        std::string visit_output_filename = "hier_eigs.silo";
        silo_db.create(visit_output_filename);
        std::cout << green << "Reference solution, FEM" << reset << std::endl;
        std::cout << "Mesh size: " << mesh_h(ref_msh) << std::endl;
        cfem_eigenvalue_solver(ref_msh);
        std::cout << green << "Computed solution, HHO" << reset << std::endl;
        std::cout << "Mesh size: " << mesh_h(sol_msh) << std::endl;
        hho_eigenvalue_solver(sol_msh, hho_degree);



        std::vector<T> fem_grad_x, fem_grad_y, hho_grad_x, hho_grad_y;
        fem_grad_x.reserve( ref_msh.cells_size() );
        fem_grad_y.reserve( ref_msh.cells_size() );
        hho_grad_x.reserve( ref_msh.cells_size() );
        hho_grad_y.reserve( ref_msh.cells_size() );

        size_t num_cell_dofs = scalar_basis_size(hho_degree, 2);
        size_t num_face_dofs = scalar_basis_size(hho_degree, 1);

        auto exp_sol_size = num_cell_dofs * sol_msh.cells_size() +
                            num_face_dofs * sol_msh.faces_size();

        dynamic_vector<T> hho_exp_sol = dynamic_vector<T>::Zero(exp_sol_size);

        hho_exp_sol.block(0, 0, num_cell_dofs * sol_msh.cells_size(), 1) =
            hho_eigvecs.block(0, which_eigvec, num_cell_dofs * sol_msh.cells_size(), 1);

        auto num_int_faces = sol_msh.internal_faces_size();
        for (size_t face_i = 0; face_i < num_int_faces; face_i++)
        {
            auto face_base = num_cell_dofs * sol_msh.cells_size();
            auto face_exp_ofs = face_base + hho_expand_map.at(face_i) * num_face_dofs;
            auto face_comp_ofs = face_base + face_i * num_face_dofs;

            hho_exp_sol.block(face_exp_ofs, 0, num_face_dofs, 1) =
                hho_eigvecs.block(face_comp_ofs, which_eigvec, num_face_dofs, 1);
        }



        std::vector<dynamic_vector<T>> rec_sols;

        size_t sol_cell_i = 0;
        for (auto& sol_cl : sol_msh)
        {
            hho_degree_info hdi(hho_degree);
            auto gr = make_hho_scalar_laplacian(sol_msh, sol_cl, hdi);

            dynamic_vector<T> hho_evec = dynamic_vector<T>::Zero(num_cell_dofs + 3*num_face_dofs);
            hho_evec.head(num_cell_dofs) = hho_exp_sol.block(sol_cell_i*num_cell_dofs, 0, num_cell_dofs, 1);
            auto sol_fcs = faces(sol_msh, sol_cl);
            for (size_t i = 0; i < sol_fcs.size(); i++)
            {
                auto fc = sol_fcs[i];
                auto fc_num = sol_msh.lookup(fc);
                auto fc_ofs = num_cell_dofs * sol_msh.cells_size() + fc_num * num_face_dofs;
                hho_evec.block(num_cell_dofs+i*num_face_dofs, 0, num_face_dofs, 1) =
                    hho_exp_sol.block(fc_ofs, 0, num_face_dofs, 1);
            }

            rec_sols.push_back( gr.first * hho_evec );
            sol_cell_i++;
        }

        dynamic_vector<T> fem_evec = fem_eigvecs.block(0,which_eigvec,fem_eigvecs.rows(),1);
        dynamic_vector<T> e_fem_evec = dynamic_vector<T>::Zero(ref_msh.points_size());
        for (size_t i = 0; i < fem_evec.size(); i++)
            e_fem_evec( fem_expand_map.at(i) ) = fem_evec(i);

        T H1_error = 0.0;
        T H1_ierror = 0.0;
        size_t cell_i = 0;
        for (auto& ref_cl : ref_msh)
        {
            auto cfem_dphi = disk::cfem::eval_basis_grad(ref_msh, ref_cl);
            auto ptids = ref_cl.point_ids();

            Eigen::Matrix<T,1,2> fem_grad = Eigen::Matrix<T,1,2>::Zero();
            for (size_t i = 0; i < 3; i++)
                fem_grad += e_fem_evec(ptids[i]) * cfem_dphi.block(i,0,1,2);

            fem_grad_x.push_back( fem_grad(0) );
            fem_grad_y.push_back( fem_grad(1) );

            /* find parent cell */
            auto bar = barycenter(ref_msh, ref_cl);
            size_t sol_cl_ofs = mesh_hier.locate_point(bar, which_level);
            auto sol_cl = *std::next(sol_msh.cells_begin(), sol_cl_ofs);

            dynamic_vector<T> rec_sol = rec_sols.at(sol_cl_ofs);

            auto sol_cb = make_scalar_monomial_basis(sol_msh, sol_cl, hho_degree+1);
            auto qps = integrate(ref_msh, ref_cl, hho_degree+2);
            for (auto& qp : qps)
            {
                auto hho_dphi = sol_cb.eval_gradients(qp.point());

                Eigen::Matrix<T,1,2> hho_grad = Eigen::Matrix<T,1,2>::Zero();
                for (size_t i = 1; i < sol_cb.size(); i++)
                    hho_grad += rec_sol(i-1) * hho_dphi.block(i,0,1,2);

                Eigen::Matrix<T,1,2> diff = hho_grad - fem_grad;
                Eigen::Matrix<T,1,2> idiff = hho_grad + fem_grad;

                H1_error += qp.weight() * diff.dot(diff);
                H1_ierror += qp.weight() * idiff.dot(idiff);
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

        std::cout << "H1 error: " << std::sqrt(H1_error) << std::endl;
        std::cout << "H1 error: " << std::sqrt(H1_ierror) << std::endl;

        silo_db.close();
    }

};

int main(int argc, char **argv)
{
    using scalar_type = double;

    if (argc != 2)
    {
        std::cout << "Please specify configuration file" << std::endl;
        return 1;
    }

    hierarchical_eigval_solver<scalar_type> hes(argv[1]);

    hes.run_computations();

}
