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
#include <regex>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/cfem/cfem.hpp"
#include "diskpp/output/silo.hpp"
#include "common.hpp"

template<typename T>
class cfem_fix_solver
{
    typedef disk::simplicial_mesh<T, 2>     mesh_type;
    typedef Eigen::SparseMatrix<T>          sparse_matrix_type;
    typedef Eigen::Triplet<T>               triplet_type;
    typedef disk::scalar_boundary_conditions<mesh_type> boundary_type;
    typedef static_vector<T, 3>             vector_type;
    typedef static_matrix<T, 3, 3>          matrix_type;

    std::vector<bool> dirichlet_nodes;
    std::vector<size_t> compress_map, expand_map;
    size_t system_size;
    disk::dynamic_vector<T>     e_gx_old;

public:
    cfem_fix_solver(const disk::simplicial_mesh<T, 2>& msh,
                const boundary_type& bnd)
    {

        dirichlet_nodes = std::vector<bool>( msh.points_size() );
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

        compress_map.resize( msh.points_size() );
        system_size = std::count_if(dirichlet_nodes.begin(), dirichlet_nodes.end(), [](bool d) -> bool {return !d;});
        expand_map.resize( system_size );

        auto nnum = 0;
        for (size_t i = 0; i < msh.points_size(); i++)
        {
            if ( dirichlet_nodes.at(i) )
                continue;

            expand_map.at(nnum) = i;
            compress_map.at(i) = nnum++;
        }
    }

    template<typename Function>
    auto
    solve(const mesh_type& msh,
        const Function & f,
        const boundary_type& bnd,
        const algorithm_parameters<T>& ap)
    {
        T tolerance = 5.e-8;
        T maxiter   = 3000;
        e_gx_old = disk::dynamic_vector<T>(msh.points_size());

        for(size_t iter= 0; iter < maxiter ; iter++)
        {
            sparse_matrix_type              gA(system_size, system_size);
            disk::dynamic_vector<T>               gb(system_size), gx(system_size);
            gb = disk::dynamic_vector<T>::Zero(system_size);
            std::vector<triplet_type>       triplets;

            auto is_contact_vector = make_is_contact_vector(msh, bnd);
            for (auto& cl : msh)
            {
                auto ptids = cl.point_ids();

                matrix_type Ah = disk::cfem::stiffness_matrix(msh, cl);
                vector_type Lh = disk::cfem::make_rhs(msh, cl, f);

                vector_type Bnegative  = vector_type::Zero();
                matrix_type Anitsche   = matrix_type::Zero();

                auto cl_id = msh.lookup(cl);
                if(is_contact_vector.at(cl_id))
                {
                    vector_type uloc;
                    for (size_t i = 0; i < 3; i++)
                        uloc(i) = e_gx_old(ptids[i]);

                    Anitsche  = make_fem_nitsche(msh, cl, bnd, ap.gamma_0, ap.theta);
                    Bnegative = make_fem_negative(msh, cl, bnd, ap.gamma_0, ap.theta, uloc);
                }

                matrix_type A = Ah - Anitsche;
                vector_type b = Lh - Bnegative;

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

            //std::cout << "Starting linear solver..." << std::endl;
            //std::cout << " * Solving for " << systsz << " unknowns." << std::endl;
            //std::cout << " * Matrix fill: " << 100.0*double(nnz)/(systsz*systsz) << "%" << std::endl;

            solver.analyzePattern(gA);
            solver.factorize(gA);
            gx = solver.solve(gb);

            disk::dynamic_vector<T> e_gx(msh.points_size());
            e_gx = disk::dynamic_vector<T>::Zero(msh.points_size());

            for (size_t i = 0; i < gx.size(); i++)
                e_gx( expand_map.at(i) ) = gx(i);

            T erroru (0);
            T errord (0);
            disk::dynamic_vector<T> diff_gx = e_gx - e_gx_old;

            for (auto& cl : msh)
            {
                matrix_type   mass_matrix = disk::cfem::mass_matrix(msh, cl);
                vector_type   diff_loc, uloc;

                auto ptids = cl.point_ids();

                for (size_t i = 0; i < 3; i++)
                {
                    uloc(i)     = e_gx_old(ptids[i]);
                    diff_loc(i) = diff_gx(ptids[i]);
                }

                erroru = std::sqrt( uloc.transpose() * mass_matrix * uloc);
                errord = std::sqrt( diff_loc.transpose() * mass_matrix * diff_loc);
            }

            e_gx_old = e_gx;

            std::cout << "error ( "<< iter<<") :    "<< errord << std::endl;
            //if(errord/erroru < tolerance)
            if(errord < tolerance)
                return 0;
        }
        return 1;
    }

    auto
    postprocess(const mesh_type& msh)
    {
        //dump_to_matlab(msh, "mesh.m");

        std::cout << "mesh size : "<< disk::average_diameter(msh) << std::endl;
        std::ofstream ofs("solution_ffem.dat");

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
                val += e_gx_old(ptids[i]) * phi(i);

            ofs << bar.x() << " " << bar.y() << " " << val << std::endl;
            solution_vals.push_back(val);
            cellnum++;
        }

        ofs.close();


        disk::silo_database silo_db;
        silo_db.create("test_diskpp.silo");
        silo_db.add_mesh(msh, "test");

        disk::silo_zonal_variable<double> u("u", solution_vals);
        silo_db.add_variable("mesh", u);

        for (size_t i = 0; i < e_gx_old.size(); i++)
            solution_vals_nodes[i] = e_gx_old(i);

        disk::silo_nodal_variable<double> u_nodal("u_nodal", solution_vals_nodes);
        silo_db.add_variable("mesh", u_nodal);

        silo_db.close();

        return;
    }
};


template<typename T>
void
run_signorini(  const disk::simplicial_mesh<T, 2>& msh,
                algorithm_parameters<T>& ap)
{
    typedef typename disk::simplicial_mesh<T, 2>::point_type point_type;

    auto force = [](const point_type& pt) -> T {
        return - 2. * M_PI * std::sin( 2. * M_PI * pt.x());
    };

    auto zero_fun = [](const point_type& p) -> T {
        return 0;
    };

    disk::scalar_boundary_conditions<disk::simplicial_mesh<T, 2>> bnd(msh);

    bnd.addDirichletBC(disk::DIRICHLET,1,zero_fun); //TOP
    bnd.addNeumannBC(disk::NEUMANN,2,zero_fun); //
    bnd.addNeumannBC(disk::NEUMANN,4,zero_fun); //
    bnd.addContactBC(disk::SIGNORINI,3); //BOTTOM


    cfem_fix_solver<T> ffem(msh, bnd); //, ap, rhs_fun);
    ffem.solve(msh, force, bnd, ap);
    ffem.postprocess(msh);

    return;
}



int main(int argc, char **argv)
{
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    char    *filename       = nullptr;
    using T = double;

    int ch;
    algorithm_parameters<T> ap;

    while ( (ch = getopt(argc, argv, "g:npz")) != -1 )
    {
        switch(ch)
        {
            case 'g':
                std::cout << "choosing gamma" << std::endl;
                ap.gamma_0 = atof(optarg);
                if (ap.gamma_0 <= 0)
                {
                    std::cout << "gamma_0 must be >0. Falling back to 0.1" << std::endl;
                    ap.gamma_0 = 0.1;
                }
                break;
            case 'n':
                ap.theta = -1.;
                std::cout << "theta negative chosen" << std::endl;
                break;
            case 'p':
                ap.theta = 1.;
                std::cout << "theta positive chosen" << std::endl;
                break;
            case 'z':
                ap.theta = 0.;
                std::cout << "theta zero chosen" << std::endl;
                break;
            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    filename = argv[0];

    /* Medit 2d*/
    if (std::regex_match(filename, std::regex(".*\\.mesh2d$") ))
    {
        typedef disk::simplicial_mesh<T, 2>  mesh_type;
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        mesh_type msh;
        disk::netgen_mesh_loader<T, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        run_signorini(msh, ap);
    }

    return 0;
}
