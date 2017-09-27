/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
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
#include "output/silo.hpp"
#include "solvers/solver.hpp"

#include "contrib/sol2/sol.hpp"
#include "contrib/timecounter.h"

#include "core/output/hdf5_io.hpp"

template<typename Mesh>
bool
hho_solver(sol::state& lua, const Mesh& msh)
{
    typedef Mesh                                       mesh_type;
    typedef typename mesh_type::scalar_type            scalar_type;
    typedef typename mesh_type::cell                   cell_type;
    typedef typename mesh_type::face                   face_type;

    typedef disk::quadrature<mesh_type, cell_type>      cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>      face_quadrature_type;

    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

    typedef
    disk::basis_quadrature_data<mesh_type,
                                disk::scaled_monomial_scalar_basis,
                                disk::quadrature> bqdata_type;

    typedef disk::gradient_reconstruction_bq<bqdata_type>               gradrec_type;
    typedef disk::diffusion_like_stabilization_bq<bqdata_type>          stab_type;
    typedef disk::diffusion_like_static_condensation_bq<bqdata_type>    statcond_type;
    //typedef disk::assembler<mesh_type, face_basis_type, face_quadrature_type> assembler_type;
    typedef disk::assembler_homogeneus_dirichlet<mesh_type, face_basis_type, face_quadrature_type> assembler_type;
    typename assembler_type::sparse_matrix_type     system_matrix;
    typename assembler_type::vector_type            system_rhs, system_solution;
    typedef static_matrix<scalar_type, mesh_type::dimension, mesh_type::dimension> tensor_type;

    size_t          degree      = lua["config"]["degree"].get_or(1);

    auto bqd = bqdata_type(degree, degree);

    auto gradrec    = gradrec_type(bqd);
    auto stab       = stab_type(bqd);
    auto statcond   = statcond_type(bqd);
    auto assembler  = assembler_type(msh, degree);

    auto lua_rhs_fun = lua["right_hand_side"];
    if ( !lua_rhs_fun.valid() )
    {
        std::cout << "[config] right_hand_side() not defined" << std::endl;
        return false;
    }

    auto f = [&](const typename mesh_type::point_type& pt) -> scalar_type {
        return lua_rhs_fun(pt.x(), pt.y());
    };

    auto lua_matparam_fun = lua["material_parameter"];
    if ( !lua_matparam_fun.valid() )
    {
        std::cout << "[config] material_parameter() not defined" << std::endl;
        return false;
    }

    auto kappa = [&](const typename mesh_type::point_type& pt) -> tensor_type {
        tensor_type ret = tensor_type::Identity();
        scalar_type k = lua_matparam_fun(pt.x(), pt.y());
        return ret * k;
    };

    size_t elem_i = 0;
    for (auto& cl : msh)
    {
        elem_i++;

        gradrec.compute(msh, cl, kappa);
        stab.compute(msh, cl, gradrec.oper);
        auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, f, degree);
        dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
        auto scnp = statcond.compute(msh, cl, loc, cell_rhs);
        assembler.assemble(msh, cl, scnp);
        if (elem_i % 10000 == 0)
            std::cout << "\33[2K\rAssembly " << double(100*elem_i)/msh.cells_size() << "%" << std::flush;
    }

    std::cout << std::endl;

    //auto bcf = [](const typename mesh_type::point_type& pt) -> scalar_type {
    //    return 0;
    //};

    //assembler.impose_boundary_conditions(msh, bcf);
    assembler.finalize(system_matrix, system_rhs);

    size_t systsz = system_matrix.rows();
    size_t nnz = system_matrix.nonZeros();

    std::cout << "Starting linear solver..." << std::endl;
    std::cout << " * Solving for " << systsz << " unknowns." << std::endl;
    std::cout << " * Matrix fill: " << 100.0*double(nnz)/(systsz*systsz) << "%" << std::endl;

    dynamic_vector<scalar_type> sol = dynamic_vector<scalar_type>::Zero(systsz);
    disk::solvers::linear_solver(lua, system_matrix, system_rhs, sol);

    system_solution = assembler.expand_solution(msh, sol);
    //system_solution = sol;

    auto num_cell_dofs = bqd.cell_basis.size();
    auto num_face_dofs = bqd.face_basis.size();

    auto num_rec_dofs = bqd.cell_basis.computed_size();

    dynamic_vector<scalar_type> reconstructed_solution =
        dynamic_vector<scalar_type>::Zero(num_rec_dofs*msh.cells_size());

    std::vector<scalar_type> vals(msh.cells_size());
    std::vector<scalar_type> recvals(msh.cells_size());

    auto expected_solution = lua["expected_solution"];
    scalar_type error = 0.0;

    elem_i = 0;
    for (auto& cl : msh)
    {

        gradrec.compute(msh, cl, kappa);
        stab.compute(msh, cl, gradrec.oper);
        auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, f, degree);
        dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;


        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();
        auto num_local_face_dofs = num_face_dofs*num_faces;
        dynamic_vector<scalar_type> local_face_dofs =
            dynamic_vector<scalar_type>::Zero(num_local_face_dofs);

        size_t face_i = 0;
        for (auto& fc : fcs)
        {
            size_t face_ofs = msh.lookup(fc);
            auto local_face_ofs = face_i*num_face_dofs;
            auto global_face_ofs = face_ofs*num_face_dofs;
            local_face_dofs.segment(local_face_ofs, num_face_dofs) =
                system_solution.segment(global_face_ofs, num_face_dofs);

            face_i++;
        }

        dynamic_vector<scalar_type> local_dofs =
            statcond.recover(msh, cl, loc, cell_rhs, local_face_dofs);

        dynamic_vector<scalar_type> R = gradrec.oper*local_dofs;

        reconstructed_solution.segment(elem_i*num_rec_dofs+1, num_rec_dofs-1) = R;
        reconstructed_solution(elem_i*num_rec_dofs) = local_dofs(0);

        auto bar = barycenter(msh, cl);
        dynamic_vector<scalar_type> vT = local_dofs.segment(0, num_cell_dofs);
        auto phi = bqd.cell_basis.eval_functions(msh, cl, bar);
        vals.at(elem_i) = vT.dot(phi);

        auto phi1 = bqd.cell_basis.eval_functions(msh, cl, bar, 1, degree+1);
        recvals.at(elem_i) = R.dot(phi1) + vT(0);
        elem_i++;

        if ( expected_solution.valid() )
        {
            dynamic_matrix<scalar_type> cell_mm =
                dynamic_matrix<scalar_type>::Zero(num_cell_dofs, num_cell_dofs);

            dynamic_vector<scalar_type> cell_rhs =
                dynamic_vector<scalar_type>::Zero(num_cell_dofs);

            auto ep_f = [&](const typename mesh_type::point_type& pt) -> scalar_type {
                return expected_solution(pt.x(), pt.y());
            };

            auto cell_quadpoints = bqd.cell_quadrature.integrate(msh, cl);
            for (auto& qp : cell_quadpoints)
            {
                auto phi    = bqd.cell_basis.eval_functions(msh, cl, qp.point());
                auto wphi   = qp.weight() * phi;
                cell_mm     += wphi * phi.transpose();
                cell_rhs    += wphi * ep_f(qp.point());
            }

            dynamic_vector<scalar_type> vT_exact = cell_mm.llt().solve(cell_rhs);
            dynamic_vector<scalar_type> diff = vT_exact - vT;
            error += diff.dot(cell_mm*diff);
        }
    }

    if ( expected_solution.valid() )
        std::cout << "L2-norm error: " << sqrt(error) << std::endl;

    auto hdf5_output_filename = lua["config"]["hdf5_output"];
    if ( hdf5_output_filename.valid() )
    {
        std::string fn = hdf5_output_filename;

        hid_t file_id = H5Fcreate(fn.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        hid_t group_id = H5Gcreate(file_id, "/hho", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        hsize_t dim1_sz = reconstructed_solution.size();
        hid_t   dataspace_id = H5Screate_simple(1, &dim1_sz, NULL);

        hid_t   dataset_id = H5Dcreate(group_id, "solution", H5T_IEEE_F64LE, dataspace_id,
                                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, reconstructed_solution.data());


        hsize_t dims = 1;
        hid_t attr_dataspace_id = H5Screate_simple(1, &dims, NULL);

        hid_t attribute_id = H5Acreate2(dataset_id, "degree", H5T_STD_I32BE,
                                        attr_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

        int attr_data = bqd.cell_basis.degree()+1;
        H5Awrite(attribute_id, H5T_NATIVE_INT, &attr_data);

        H5Aclose(attribute_id);

        attribute_id = H5Acreate2(dataset_id, "elems", H5T_STD_I32BE,
                                        attr_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

        attr_data = msh.cells_size();
        H5Awrite(attribute_id, H5T_NATIVE_INT, &attr_data);

        H5Aclose(attribute_id);

        H5Sclose(attr_dataspace_id);

        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
    }

    auto visit_output_filename = lua["config"]["visit_output"];
    if ( visit_output_filename.valid() )
    {
        disk::silo_database silo_db;
        silo_db.create(visit_output_filename);
        silo_db.add_mesh(msh, "mesh");

        disk::silo_zonal_variable<scalar_type> u("u", vals);
        silo_db.add_variable("mesh", u);

        disk::silo_zonal_variable<scalar_type> Ru("Ru", recvals);
        silo_db.add_variable("mesh", Ru);

        silo_db.close();
    }

    return true;
}

template<typename T>
bool
cfem_eigenvalue_solver(sol::state& lua, const disk::simplicial_mesh<T, 2>& msh)
{
    feast_eigensolver_params<T> fep;

    if ( !setup_feast(lua, fep) )
    {
        std::cout << "Problems with parameters. Check configuration." << std::endl;
        return false;
    }

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

    std::vector<size_t> compress_map, expand_map;
    compress_map.resize( msh.points_size() );
    size_t system_size = std::count_if(dirichlet_nodes.begin(), dirichlet_nodes.end(), [](bool d) -> bool {return !d;});
    expand_map.resize( system_size );

    auto nnum = 0;
    for (size_t i = 0; i < msh.points_size(); i++)
    {
        if ( dirichlet_nodes.at(i) )
            continue;

        expand_map.at(nnum) = i;
        compress_map.at(i) = nnum++;
    }

    sparse_matrix_type  gK(system_size, system_size);
    sparse_matrix_type  gM(system_size, system_size);

    std::vector<triplet_type> gtK, gtM;

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

                gtK.push_back( triplet_type(compress_map.at(ptids[i]),
                                            compress_map.at(ptids[j]),
                                            K(i,j)) );

                gtM.push_back( triplet_type(compress_map.at(ptids[i]),
                                            compress_map.at(ptids[j]),
                                            M(i,j)) );
            }
        }
    }

    gK.setFromTriplets(gtK.begin(), gtK.end());
    gM.setFromTriplets(gtM.begin(), gtM.end());

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>   eigvecs;
    Eigen::Matrix<double, Eigen::Dynamic, 1>                eigvals;

    bool solver_ok = generalized_eigenvalue_solver(fep, gK, gM, eigvecs, eigvals);
    if (!solver_ok)
    {
        std::cout << "Problem during solution phase" << std::endl;
        return false;
    }

    auto eigenvalues_output_filename = lua["config"]["eigval_output"];
    if ( eigenvalues_output_filename.valid() )
    {
        std::string fname = eigenvalues_output_filename;

        std::ofstream ofs(fname);

        for (size_t i = 0; i < fep.eigvals_found; i++)
            ofs << std::setprecision(10) << eigvals(i) << std::endl;

        ofs.close();
    }

    auto visit_output_filename = lua["config"]["visit_output"];
    if ( visit_output_filename.valid() )
    {
        disk::silo_database silo_db;
        silo_db.create(visit_output_filename);
        silo_db.add_mesh(msh, "mesh");

        for (size_t i = 0; i < fep.eigvals_found; i++)
        {
            std::vector<T> solution_vals;
            solution_vals.resize(msh.points_size());

            dynamic_vector<T> gx = eigvecs.block(0,i,eigvecs.rows(),1);
            for (size_t i = 0; i < gx.size(); i++)
                solution_vals.at( expand_map.at(i) ) = gx(i);

            std::stringstream ss;
            ss << "u" << i;

            disk::silo_nodal_variable<T> u(ss.str(), solution_vals);
            silo_db.add_variable("mesh", u);
        }

        silo_db.close();
    }

    return true;
}





template<typename T>
bool
cfem_solver(sol::state& lua, const disk::simplicial_mesh<T, 2>& msh)
{
    auto lua_rhs_fun = lua["right_hand_side"];
    if ( !lua_rhs_fun.valid() )
    {
        std::cout << "[config] right_hand_side() not defined" << std::endl;
        return false;
    }

    auto f = [&](const typename disk::simplicial_mesh<T, 2>::point_type& pt) -> T {
        return lua_rhs_fun(pt.x(), pt.y());
    };

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

    std::vector<size_t> compress_map, expand_map;
    compress_map.resize( msh.points_size() );
    size_t system_size = std::count_if(dirichlet_nodes.begin(), dirichlet_nodes.end(), [](bool d) -> bool {return !d;});
    expand_map.resize( system_size );

    auto nnum = 0;
    for (size_t i = 0; i < msh.points_size(); i++)
    {
        if ( dirichlet_nodes.at(i) )
            continue;

        expand_map.at(nnum) = i;
        compress_map.at(i) = nnum++;
    }

    sparse_matrix_type              gA(system_size, system_size);
    dynamic_vector<T>               gb(system_size), gx(system_size);

    gb = dynamic_vector<T>::Zero(system_size);


    std::vector<triplet_type>       triplets;

    timecounter tc;
    tc.tic();

    auto lua_matparam_fun = lua["material_parameter"];
    if ( !lua_matparam_fun.valid() )
    {
        std::cout << "[config] material_parameter() not defined" << std::endl;
        return false;
    }

    for (auto& cl : msh)
    {
        static_matrix<T, 2, 2> kappa = static_matrix<T, 2, 2>::Zero();
        auto bar = barycenter(msh, cl);

        T eps = lua_matparam_fun(bar.x(), bar.y());

        kappa(0,0) = eps;
        kappa(1,1) = eps;

        auto A = disk::cfem::stiffness_matrix(msh, cl, kappa);
        auto b = disk::cfem::make_rhs(msh, cl, f);

        auto ptids = cl.point_ids();

        for (size_t i = 0; i < A.rows(); i++)
        {
            if ( dirichlet_nodes.at(ptids[i]) )
                continue;

            for (size_t j = 0; j < A.cols(); j++)
            {
                if ( dirichlet_nodes.at(ptids[j]) )
                    continue;

                triplets.push_back( triplet_type(compress_map.at(ptids[i]),
                                                 compress_map.at(ptids[j]),
                                                 A(i,j)) );
            }

            gb(compress_map.at(ptids[i])) += b(i);
        }
    }

    tc.toc();

    std::cout << "Assembly time: " << tc << " seconds" << std::endl;

    gA.setFromTriplets(triplets.begin(), triplets.end());

    size_t systsz = gA.rows();
    size_t nnz = gA.nonZeros();

    std::cout << "Starting linear solver..." << std::endl;
    std::cout << " * Solving for " << systsz << " unknowns." << std::endl;
    std::cout << " * Matrix fill: " << 100.0*double(nnz)/(systsz*systsz) << "%" << std::endl;

    disk::solvers::linear_solver(lua, gA, gb, gx);

    dynamic_vector<T> e_gx(msh.points_size());
    e_gx = dynamic_vector<T>::Zero(msh.points_size());

    for (size_t i = 0; i < gx.size(); i++)
        e_gx( expand_map.at(i) ) = gx(i);

    auto hdf5_output_filename = lua["config"]["hdf5_output"];
    if ( hdf5_output_filename.valid() )
    {
        //assert(std::is_same<T, double>::value);

        std::string fn = hdf5_output_filename;

        hid_t file_id = hdf5_open_or_create( fn.c_str() );

        herr_t status = H5Eset_auto(H5E_DEFAULT, NULL, NULL);
        status = H5Gget_objinfo (file_id, "/fem", 0, NULL);

        hid_t group_id;

        if (status != 0)
            group_id = H5Gcreate(file_id, "/fem", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        else
            group_id = H5Gopen(file_id, "/fem", H5P_DEFAULT);

        hsize_t dim1_sz = e_gx.size();
        hid_t   dataspace_id = H5Screate_simple(1, &dim1_sz, NULL);

        hid_t   dataset_id = H5Dcreate(group_id, "solution", H5T_IEEE_F64LE, dataspace_id,
                                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, e_gx.data());

        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
    }

    std::ofstream ofs("solution.dat");

    std::vector<T> solution_vals;
    solution_vals.reserve(msh.cells_size());

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
        solution_vals.push_back(val);
        cellnum++;
    }

    ofs.close();

    auto visit_output_filename = lua["config"]["visit_output"];
    if ( visit_output_filename.valid() )
    {
        disk::silo_database silo_db;
        silo_db.create(visit_output_filename);
        silo_db.add_mesh(msh, "mesh");

        disk::silo_zonal_variable<T> u("u", solution_vals);
        silo_db.add_variable("mesh", u);
        silo_db.close();
    }

    return true;
}






template<typename Mesh>
bool
eigval_solver(sol::state& lua, const Mesh& msh)
{
    typedef Mesh                                       mesh_type;
    typedef typename mesh_type::scalar_type            scalar_type;
    typedef typename mesh_type::cell                   cell_type;
    typedef typename mesh_type::face                   face_type;

    typedef disk::quadrature<mesh_type, cell_type>      cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>      face_quadrature_type;

    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

    typedef
    disk::basis_quadrature_data<mesh_type,
                                disk::scaled_monomial_scalar_basis,
                                disk::quadrature> bqdata_type;

    typedef disk::gradient_reconstruction_bq<bqdata_type>               gradrec_type;
    typedef disk::diffusion_like_stabilization_bq<bqdata_type>          stab_type;
    typedef disk::eigval_mass_matrix_bq<bqdata_type>                    mass_type;

    typedef Eigen::Triplet<scalar_type>                                 triplet_type;


    feast_eigensolver_params<scalar_type> fep;

    if ( !setup_feast(lua, fep) )
    {
        std::cout << "Problems with parameters. Check configuration." << std::endl;
        return false;
    }


    size_t          degree      = lua["config"]["degree"].get_or(1);

    auto bqd = bqdata_type(degree, degree);

    auto gradrec    = gradrec_type(bqd);
    auto stab       = stab_type(bqd);
    auto mass       = mass_type(bqd);

    std::cout << "Assembly" << std::endl;

    size_t num_cell_dofs = howmany_dofs(bqd.cell_basis);
    size_t num_face_dofs = howmany_dofs(bqd.face_basis);
    size_t num_global_cell_dofs = msh.cells_size() * num_cell_dofs;
    size_t num_global_face_dofs = msh.internal_faces_size() * num_face_dofs;
    size_t num_global_dofs = num_global_cell_dofs + num_global_face_dofs;

    /* Strong imposition of dirichlet bc */
    std::vector<size_t> compress_map, expand_map;
    compress_map.resize( msh.faces_size() );
    size_t num_non_dirichlet_faces = msh.internal_faces_size();
    size_t system_size = num_non_dirichlet_faces * num_face_dofs;
    expand_map.resize( system_size );

    size_t fnum = 0;
    size_t face_i = 0;
    for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
    {
        if ( msh.is_boundary(*itor) )
        {
            face_i++;
            continue;
        }

        expand_map.at(fnum) = face_i;
        compress_map.at(face_i) = fnum++;
        face_i++;
    }

    sparse_matrix<scalar_type> gK(num_global_dofs, num_global_dofs);
    sparse_matrix<scalar_type> gM(num_global_dofs, num_global_dofs);

    std::vector<triplet_type> gtK, gtM;

    scalar_type stab_weight = lua["config"]["stabilization_weight"].get_or(1.0);

    size_t elem_i = 0;
    for (auto& cl : msh)
    {
        gradrec.compute(msh, cl);
        stab.compute(msh, cl, gradrec.oper);
        mass.compute(msh, cl);

        auto K = gradrec.data + stab_weight * stab.data;
        auto M = mass.data;

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
                auto face_pos = compress_map.at(face_id);
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

    Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic> eigvecs;
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> eigvals;

    bool only_dump = lua["config"]["only_dump"].get_or(false);
    if (only_dump)
    {
        dump_sparse_matrix(gK, "stiff_matrix.txt");
        dump_sparse_matrix(gM, "mass_matrix.txt");
        return true;
    }

    generalized_eigenvalue_solver(fep, gK, gM, eigvecs, eigvals);

    std::cout << "Postprocessing" << std::endl;

    auto eigenvalues_output_filename = lua["config"]["eigval_output"];
    if ( eigenvalues_output_filename.valid() )
    {
        std::string fname = eigenvalues_output_filename;

        std::ofstream ofs(fname);

        for (size_t i = 0; i < fep.eigvals_found; i++)
            ofs << std::setprecision(10) << eigvals(i) << std::endl;

        ofs.close();
    }

    auto visit_output_filename = lua["config"]["visit_output"];
    if ( visit_output_filename.valid() )
    {
        disk::silo_database silo_db;
        silo_db.create(visit_output_filename);
        silo_db.add_mesh(msh, "mesh");

        for (size_t i = 0; i < fep.eigvals_found; i++)
        {
            std::vector<scalar_type> solution_vals;
            solution_vals.resize(msh.cells_size());

            dynamic_vector<scalar_type> gx = eigvecs.block(0,i,eigvecs.rows(),1);

            size_t cell_i = 0;
            for (auto& cl : msh)
            {
                auto bar = barycenter(msh, cl);
                auto local_dofs = gx.block(cell_i*num_cell_dofs, 0, num_cell_dofs, 1);
                auto phi = bqd.cell_basis.eval_functions(msh, cl, bar);
                auto val = local_dofs.dot(phi);
                solution_vals.at(cell_i) = val;
                cell_i++;
            }

            std::stringstream ss;
            ss << "u" << i;

            disk::silo_zonal_variable<scalar_type> u(ss.str(), solution_vals);
            silo_db.add_variable("mesh", u);
        }

        silo_db.close();
    }

    return true;
}

template<typename T>
bool
eigval_reference(sol::state& lua, const disk::simplicial_mesh<T, 2>& msh)
{
    disk::silo_database silo_db;
    auto visit_output_filename = lua["config"]["visit_output"];
    auto mode_x_max = lua["config"]["mode_x_max"].get_or(5);
    auto mode_y_max = lua["config"]["mode_y_max"].get_or(5);

    if ( visit_output_filename.valid() )
    {
        silo_db.create(visit_output_filename);
        silo_db.add_mesh(msh, "mesh");

        for (size_t i = 0; i < mode_x_max; i++)
        {
            for (size_t j = 0; j < mode_y_max; j++)
            {
                std::stringstream ss;
                ss << "ref_" << i << "_" << j;

                std::cout << i << " " << j << " " << M_PI*(i+1)*M_PI*(i+1) + M_PI*(j+1)*M_PI*(j+1) << std::endl;

                std::vector<T> solution_vals;
                solution_vals.reserve(msh.points_size());
                for (auto itor = msh.points_begin(); itor != msh.points_end(); itor++)
                {
                    auto pt = *itor;
                    auto ef = std::sin((i+1)*pt.x()*M_PI) * std::sin((j+1)*pt.y()*M_PI);
                    solution_vals.push_back(ef);
                }

                disk::silo_nodal_variable<T> u(ss.str(), solution_vals);
                silo_db.add_variable("mesh", u);
            }
        }

        silo_db.close();
    }

    return true;
}



class fem_solver
{
public:
    std::string     input_mesh;
    std::string     hdf5_output;
    std::string     visit_output;
    int             degree;
    bool            verbose;

    fem_solver() :
        input_mesh(""),
        degree(1),
        verbose(false)
    {}

    void run()
    {}

};

void
lua_register_fem_solver_usertype(sol::state& lua)
{
    lua.new_usertype<fem_solver>( "fem_solver",
        sol::constructors<fem_solver()>(),
        "run", &fem_solver::run,
        "input_mesh", &fem_solver::input_mesh,
        "degree", &fem_solver::degree,
        "verbose", &fem_solver::verbose
        );
}




int main(int argc, char **argv)
{
    using scalar_type = double;

    if (argc != 2)
    {
        std::cout << "Please specify configuration file" << std::endl;
        return 1;
    }

    sol::state lua;
    lua.open_libraries(sol::lib::math, sol::lib::base);
    lua["config"] = lua.create_table();
    disk::solvers::init_lua(lua);
    lua["solver"]["feast"] = lua.create_table();

    lua_register_fem_solver_usertype(lua);

    auto r = lua.do_file(argv[1]);

    if ( !r.valid() )
    {
        std::cout << "Problems opening configuration file" << std::endl;
        return 1;
    }


    std::string     method      = lua["config"]["method"].get_or(std::string("hho"));
    std::string     input_mesh  = lua["config"]["input_mesh"];
    std::string     hdf5_output = lua["config"]["hdf5_output"];
    bool            verbose     = lua["config"]["verbose"].get_or(false);

    if (std::regex_match(input_mesh, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;

        typedef disk::simplicial_mesh<scalar_type, 2>  mesh_type;

        mesh_type msh;
        disk::netgen_mesh_loader<scalar_type, 2> loader;
        loader.verbose(verbose);
        if (!loader.read_mesh(input_mesh))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }

        loader.populate_mesh(msh);

        std::cout << "Mesh avg. diameter: " << mesh_h(msh) << std::endl;

        if (method == "fem")
            cfem_solver(lua, msh);
        else if (method == "hho")
            hho_solver(lua, msh);
        else if (method == "hho_eigs")
            eigval_solver(lua, msh);
        else if (method == "fem_eigs")
            cfem_eigenvalue_solver(lua, msh);
        else if (method == "reference")
            eigval_reference(lua, msh);
    }

    if (std::regex_match(input_mesh, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;

        typedef disk::generic_mesh<scalar_type, 2>  mesh_type;

        mesh_type msh;
        disk::fvca5_mesh_loader<scalar_type, 2> loader;
        loader.verbose(verbose);
        if (!loader.read_mesh(input_mesh))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }

        loader.populate_mesh(msh);

        std::cout << "Mesh avg. diameter: " << mesh_h(msh) << std::endl;

        if (method == "hho")
            hho_solver(lua, msh);
        else if (method == "hho_eigs")
            eigval_solver(lua, msh);
    }

    if (std::regex_match(input_mesh, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: Cartesian 2D" << std::endl;

        typedef disk::cartesian_mesh<scalar_type, 2>  mesh_type;

        mesh_type msh;
        disk::cartesian_mesh_loader<scalar_type, 2> loader;
        loader.verbose(verbose);
        if (!loader.read_mesh(input_mesh))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }

        loader.populate_mesh(msh);

        std::cout << "Mesh avg. diameter: " << mesh_h(msh) << std::endl;

        if (method == "hho")
            hho_solver(lua, msh);
        else if (method == "hho_eigs")
            eigval_solver(lua, msh);
    }

    return 0;
}
