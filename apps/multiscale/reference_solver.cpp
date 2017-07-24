/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#include <iostream>
#include <regex>

#include "loaders/loader.hpp"
#include "cfem/cfem.hpp"
#include "output/silo.hpp"
#include "solvers/solver.hpp"

#include "contrib/sol2/sol.hpp"
#include "contrib/timecounter.h"


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

#ifdef HAVE_INTEL_MKL
        Eigen::PardisoLU<Eigen::SparseMatrix<T>>  solver;
        //solver.pardisoParameterArray()[59] = 0; //out-of-core
#else
        Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
#endif

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













int main(int argc, char **argv)
{
    using scalar_type = double;

    if (argc != 2)
    {
        std::cout << "Please specify configuration file" << std::endl;
        return 1;
    }

    sol::state lua;
    lua.open_libraries(sol::lib::math);
    lua["config"] = lua.create_table();
    disk::solvers::init_lua(lua);

    auto r = lua.do_file(argv[1]);

    if ( !r.valid() )
    {
        std::cout << "Problems opening configuration file" << std::endl;
        return 1;
    }

    std::string     method      = lua["config"]["method"].get_or(std::string("hho"));
    size_t          degree      = lua["config"]["degree"].get_or(1);
    std::string     input_mesh  = lua["config"]["input_mesh"];
    std::string     hdf5_output = lua["config"]["hdf5_output"];
    
    if (std::regex_match(input_mesh, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;

        typedef disk::simplicial_mesh<scalar_type, 2>  mesh_type;

        mesh_type msh;
        disk::netgen_mesh_loader<scalar_type, 2> loader;
        if (!loader.read_mesh(input_mesh))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        
        loader.populate_mesh(msh);

        cfem_solver(lua, msh);
    }
    
    return 0;
}