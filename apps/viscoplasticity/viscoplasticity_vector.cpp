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

#include "viscoplasticity_vector_solver.hpp"
template<typename T>
auto
run_viscoplasticity(size_t degree,
                    const T & alpha,
                    const problem_type& problem)
{
    /* Medit 2d*/
    std::cout << "Guessed mesh format: Medit format" << std::endl;
    typedef disk::generic_mesh<T, 2>  mesh_type;

    std::string name, filename;
    switch (problem)
    {
        case DRIVEN:
            name = "driven";
            filename = "../../../diskpp/meshes/2D_quads/medit/square_h005.medit2d";
            break;
        default:
            std::cout << "wrong arguments" << std::endl;
            exit(1);
    }

    mesh_type                     msh;
    disk::medit_mesh_loader<T, 2> loader;
    loader.read_mesh(filename);
    loader.populate_mesh(msh);

    std::string info = name + "_k" + tostr(degree) + "_a" + tostr(alpha);

    hho_degree_info hdi(degree , degree);
    augmented_lagrangian_viscoplasticity<mesh_type> als(msh, hdi, alpha);

    if(!als.run_alg(msh, info, problem))
        std::cout << "No convergence" << std::endl;

    //auto final_error = als.compute_errors(msh, assembler, true);
    return;
}


int main(int argc, char **argv)
{
    using RealType = double;

    char    *word      = nullptr;
    int ch;
    size_t degree = 0;
    RealType alpha = 1.;

    problem_type problem = DRIVEN;

    while ( (ch = getopt(argc, argv, "k:a:")) != -1 )
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
            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    run_viscoplasticity(degree, alpha, problem);

    return 0;
}
//#endif
