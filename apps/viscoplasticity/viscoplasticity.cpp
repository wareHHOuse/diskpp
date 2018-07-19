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
#include "revolution/quadratures"
#include "methods/hho"

#include "core/loaders/loader.hpp"

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"
#include "output/postMesh.hpp"

#include "output/silo.hpp"

#include "viscoplasticity_solver.hpp"
template<typename T>
auto
run_viscoplasticity(size_t degree,
                    const T & alpha,
                    const problem_type& problem,
                    const std::string & other_info)
{
    /* Medit 2d*/
    std::cout << "Guessed mesh format: Medit format" << std::endl;
    typedef disk::generic_mesh<T, 2>  mesh_type;

    T tolerance = 1.e-10, Ninf = 10.e+5;
    size_t max_iters = 1; //50000;

    std::string name, filename;
    switch (problem)
    {
        case DRIVEN:
            name = "driven";
            filename = "../../../diskpp/meshes/2D_quads/medit/square_h0025.medit2d";
            break;
        case COUETTE:
            name = "couette";
            filename = "../../../diskpp/meshes/2D_triangles/medit/couronne_01.medit2d";
            break;
        default:
            std::cout << "wrong arguments" << std::endl;
            exit(1);
    }

    mesh_type                     msh;
    disk::medit_mesh_loader<T, 2> loader;
    loader.read_mesh(filename);
    loader.populate_mesh(msh);

    std::string info = name + "_k" + tostr(degree) + "_a" + tostr(alpha) + other_info;
    std::ofstream ofs("errors_" + info + ".data");

    if (!ofs.is_open())
        std::cout << "Error opening errors "<<std::endl;

    typename disk::hho_degree_info hdi(degree, degree);
    augmented_lagrangian_viscoplasticity<mesh_type> als(msh, hdi, alpha);

    auto assembler = als.define_problem(msh, problem);
    als.initialize(msh, assembler);

    size_t i;
    for(i = 0; i < max_iters; i++)
    {
        als.run_stokes_like(msh, assembler, i);

        //als.update_multiplier(msh, assembler);
        //auto error = als.compute_errors(msh, assembler, false);

        T cvg_total (0.), cvg_stress(0.), cvg_gamma(0.);
        std::tie(cvg_total, cvg_stress, cvg_gamma) = als.convergence;

        if(i % 1000 == 0)
        {
            std::cout << "  i : "<< i<<"  - " << std::sqrt(cvg_total)<<std::endl;
            als.post_processing( msh, assembler, info +"_i" + tostr(i), problem);
            std::cout << "done" << std::endl;
        }

        assert(std::sqrt(cvg_total) < Ninf);
        if( std::sqrt(cvg_total)  < tolerance)
            break;
    }
    ofs.close();

    std::cout << "Finish" << std::endl;
    auto final_error = als.compute_errors(msh, assembler, true);
    als.post_processing( msh, assembler, info +"_i" + tostr(i), problem);

    return final_error;
}


int main(int argc, char **argv)
{
    using RealType = double;

    char    *word      = nullptr;
    int ch;
    size_t degree = 0;
    RealType alpha = 1.;

    problem_type problem = DRIVEN;

    while ( (ch = getopt(argc, argv, "k:a:dcp")) != -1 )
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
                std::cout << "couette chosen" << std::endl;
                break;

            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    word = argv[0];

    if (word == nullptr)
    {
        std::cout << "no word specified" << std::endl;
        return 1;
    }

    run_viscoplasticity(degree, alpha, problem, word);

    return 0;
}
//#endif
