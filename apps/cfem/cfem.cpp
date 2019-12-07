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
#include <regex>

#include "loaders/loader.hpp"

#include "cfem/cfem.hpp"
#include "output/silo.hpp"

template<typename T>
void
cfem_solver(const disk::simplicial_mesh<T, 2>& msh)
{
    auto f = [](const typename disk::simplicial_mesh<T, 2>::point_type& pt) -> auto {
        return sin(pt.x()) * sin(pt.y());
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
    disk::dynamic_vector<T>               gb(system_size), gx(system_size);

    gb = disk::dynamic_vector<T>::Zero(system_size);


    std::vector<triplet_type>       triplets;

    //T integral = 0.0;

    for (auto& cl : msh)
    {
        disk::static_matrix<T, 2, 2> kappa = disk::static_matrix<T, 2, 2>::Zero();
        auto bar = barycenter(msh, cl);

        auto c = std::cos(M_PI * bar.x()/0.02);
        auto s = std::sin(M_PI * bar.y()/0.02);
        auto eps = 1 + 100*c*c*s*s;// + std::exp( (bar.x()*bar.x() + bar.y()*bar.y())/2. );
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

        solver.analyzePattern(gA);
        solver.factorize(gA);
        gx = solver.solve(gb);

        disk::dynamic_vector<T> e_gx(msh.points_size());
        e_gx = disk::dynamic_vector<T>::Zero(msh.points_size());

        for (size_t i = 0; i < gx.size(); i++)
            e_gx( expand_map.at(i) ) = gx(i);

        std::ofstream ofs("solution.dat");

        std::vector<double> solution_vals, solution_vals_nodes;
        solution_vals.reserve(msh.cells_size());
        solution_vals_nodes.resize(msh.points_size());

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


        disk::silo_database silo_db;
        silo_db.create("test.silo");
        silo_db.add_mesh(msh, "test");

        disk::silo_zonal_variable<double> u("u", solution_vals);
        silo_db.add_variable("mesh", u);

        for (size_t i = 0; i < e_gx.size(); i++)
            solution_vals_nodes[i] = e_gx(i);

        disk::silo_nodal_variable<double> u_nodal("u_nodal", solution_vals_nodes);
        silo_db.add_variable("mesh", u_nodal);

        silo_db.close();
}



int main(int argc, char **argv)
{
    using RealType = double;

    char    *filename       = nullptr;
    int     elems_1d        = 8;
    int ch;

    filename = argv[1];

    if (std::regex_match(filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;

        typedef disk::simplicial_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::netgen_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        cfem_solver(msh);
    }


    return 0;
}
