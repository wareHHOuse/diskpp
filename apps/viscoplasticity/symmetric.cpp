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

#include "viscoplasticity_solver.hpp"
template<typename Mesh>
auto
run_viscoplasticity(const Mesh& msh, size_t degree,
                    const typename Mesh::scalar_type & alpha,
                    const problem_type& problem)
{
    using T = typename Mesh::coordinate_type;
    T tolerance = 1.e-8, Ninf = 10.e+5;
    size_t max_iters = 100000;

    std::string name;
    switch (problem)
     {
        case DRIVEN:
            name = "driven";
            break;
        case COUETTE:
            name = "couette";
            break;
    }

    std::string info = name + "_k" + tostr(degree) + "_a" + tostr(alpha) + "_tri_h0025";
    std::ofstream ofs("errors_" + info + ".data");

    if (!ofs.is_open())
        std::cout << "Error opening errors "<<std::endl;

    typename revolution::hho_degree_info hdi(degree, degree);
    augmented_lagrangian_viscoplasticity<Mesh> als(msh, hdi, alpha);

    auto assembler = als.define_problem(msh, problem);
    als.initialize(msh, assembler);

    for(size_t i = 0; i < max_iters; i++)
    {
        als.run_stokes_like(msh, assembler, i);
        als.update_multiplier(msh, assembler);
        auto error = als.compute_errors(msh, assembler, false);

        T cvg_total (0.), cvg_stress(0.), cvg_gamma(0.);
        std::tie(cvg_total, cvg_stress, cvg_gamma) = als.convergence;
        ofs << std::sqrt(cvg_total)<< " "<< std::sqrt(cvg_stress) <<" ";
        ofs << std::sqrt(cvg_gamma)<< " "<< error.first << " " <<error.second <<std::endl;

        assert(std::sqrt(cvg_total) < Ninf);
        if( std::sqrt(cvg_total)  < tolerance)
            break;
        //if(error.first < 10.e-9 && i > 1)
        //    break;
    }
    ofs.close();

    auto final_error = als.compute_errors(msh, assembler, true);
    als.post_processing( msh, assembler, info, problem);

    return final_error;
}

//#if 0

int main(int argc, char **argv)
{
    using RealType = double;

    char    *filename       = nullptr;
    int ch;
    size_t degree = 0;
    RealType alpha = 1.;

    problem_type problem = DRIVEN;

    while ( (ch = getopt(argc, argv, "k:a:dc")) != -1 )
    {
        switch(ch)
        {
            case 'k':
                degree = atoi(optarg);
                break;

            case 'a':
                alpha = atof(optarg);
                if (alpha <= 0)
                {
                    std::cout << "alpha must be >=0. Falling back to 1." << std::endl;
                    alpha = 1.;
                }
            case 'd':
                problem = DRIVEN;
                break;

            case 'c':
                problem = COUETTE;
                break;
            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }


    argc -= optind;
    argv += optind;

    if (argc != 1)
    {
        std::cout << "Please specify a 2D mesh" << std::endl;

        return 0;
    }

    filename = argv[0];


    if (std::regex_match(filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;

        typedef disk::generic_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::fvca5_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        run_viscoplasticity(msh, degree, alpha, problem);
    }

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

        run_viscoplasticity(msh, degree, alpha, problem);
    }

    /* Medit 2d*/
    if (std::regex_match(filename, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        typedef disk::generic_mesh<RealType, 2>  mesh_type;
        mesh_type msh = disk::load_medit_2d_mesh<RealType>(filename);
        run_viscoplasticity(msh, degree, alpha, problem);
    }

    #if 0
    if (std::regex_match(filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: Cartesian 2D" << std::endl;

        typedef disk::cartesian_mesh<RealType, 2>   mesh_type;

        mesh_type msh;
        disk::cartesian_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        test_bases(msh);
    }
    #endif

    return 0;
}
//#endif
