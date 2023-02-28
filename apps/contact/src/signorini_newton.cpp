/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
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
 #include <iomanip>
 #include <regex>
 #include <unistd.h>

 #include "diskpp/geometry/geometry.hpp"
 #include "diskpp/bases/bases.hpp"
 #include "diskpp/quadratures/quadratures.hpp"
 #include "diskpp/methods/hho"
 #include "diskpp/loaders/loader.hpp"
 #include "src/signorini_newton_solver.hpp"

int main(int argc, char **argv)
{
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    char    *filename       = nullptr;
    using T = double;

    int ch;
    T parameter = 1;
    size_t degree = 1;
    algorithm_parameters<T> ap;

    while ( (ch = getopt(argc, argv, "k:g:npzfcloe:")) != -1 )
    {
        switch(ch)
        {
            case 'k':
                degree = atoi(optarg);
                if (degree < 0)
                {
                    std::cout << "Degree must be positive. Falling back to 1." << std::endl;
                    degree = 1;
                }
                ap.degree = degree;
                break;
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
                break;
            case 'p':
                ap.theta = 1.;
                break;
            case 'z':
                ap.theta = 0.;
                break;
            case 'f':
                ap.solver = EVAL_ON_FACES;
                break;
            case 'c':
                ap.solver = EVAL_IN_CELLS;
                break;
            case 'l':
                ap.solver = EVAL_IN_CELLS_FULL;
                break;
            case 'o':
                ap.solver = EVAL_IN_CELLS_AS_FACES;
                break;
            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    std::cout << ap << std::endl;

    argc -= optind;
    argv += optind;

    filename = argv[0];
    auto run_exact = true;

    /* Netgen 2d*/
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

        run_signorini(msh, ap, run_exact);
        return 0;
    }

    /* FVCA5 2D */
    if (std::regex_match(filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        auto msh = disk::load_fvca5_2d_mesh<T>(filename);
        run_signorini(msh, ap, run_exact);
        return 0;
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        auto msh = disk::load_cartesian_2d_mesh<T>(filename);
        run_signorini(msh, ap, run_exact);
        return 0;
    }


    #if 0
    /* Medit 2d*/
    if (std::regex_match(filename, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        typedef disk::generic_mesh<T, 2>  mesh_type;
        mesh_type msh = disk::load_medit_2d_mesh<T>(filename);
        run_signorini(msh, ap);
    }
    #endif

    return 0;
}
