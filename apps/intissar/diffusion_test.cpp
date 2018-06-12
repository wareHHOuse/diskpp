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
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include <map>

#include "colormanip.h"

#include "config.h"

#include "loaders/loader.hpp"
#include "hho/hho.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "diffusion_hho.cpp"
enum problem_type
{
    TEST_DIRICHLET,
    TEST_NEUMANN,
    TEST_ROBIN
};


template<typename Mesh, typename LoaderType>
bool
verify_convergence(const std::vector<std::string>& paths,
                   size_t mindeg, size_t maxdeg,
                    problem_type problem)
{
    typedef typename Mesh::scalar_type scalar_type;

    bool success = true;

    for (size_t i = mindeg; i <= maxdeg; i++)
    {
        scalar_type expected_rate = i+1;

        std::vector<std::pair<scalar_type, scalar_type>> errdiams;

        std::cout << "Convergence rates for k = " << i << ":   " << std::flush;

        for (auto& tsp : paths)
        {
            Mesh    msh;
            LoaderType  loader;

            if (!loader.read_mesh(tsp))
            {
                std::cout << "Problem loading mesh." << std::endl;
                return false;
            }
            loader.populate_mesh(msh);

            scalar_type error = 0;

            switch(problem)
            {
                case TEST_DIRICHLET:
                    error = run_test_dirichlet(msh);
                    break;

                case TEST_NEUMANN:
                    error = run_test_neumann(msh);
                    break;
                case TEST_ROBIN:
                    error = run_test_robin(msh);
                    break;

                default:
                    error = run_test_dirichlet(msh);
            }

            auto diam = average_diameter(msh);
            errdiams.push_back( std::make_pair(diam, error) );
        }

        bool pass       = true;
        bool warning    = false;
        bool high, low, ok;

        for (size_t i = 0; i < errdiams.size(); i++)
        {
            if (i == 0)
            {
                std::cout << "    ";
                std::cout << std::scientific << std::setprecision(4) << errdiams.at(i).second << "    ";
                std::cout << std::scientific << std::setprecision(4) << errdiams.at(i).first << "    ";
                std::cout << "     -- " << std::endl;
            }
            else
            {
                auto d = log2(errdiams[i-1].first/errdiams[i].first);
                auto e = log2(errdiams[i-1].second/errdiams[i].second);
                auto rate = e/d;

                std::cout << "    ";
                std::cout << std::scientific << std::setprecision(4) << errdiams.at(i).second << "    ";
                std::cout << std::scientific << std::setprecision(4) << errdiams.at(i).first  << "    ";
                std::cout << std::fixed<< std::setprecision(2) << rate << std::endl;

                ok   = (std::abs(expected_rate - rate) < 0.4); /* Test passed */
                low  = ((expected_rate - rate) > 0.2); /* a bit too low, warn */
                high = ((rate - expected_rate) > 0.2); /* a bit too high, warn */

                if (low)    std::cout << magenta;
                if (high)   std::cout << cyan;
                std::cout << std::fixed << std::setprecision(3) << "rate" << rate << "  ";
                if (low or high)
                {
                    std::cout << reset;
                    warning = true;
                }
            }
        }



        std::string             okfail = "[\x1b[31;1mFAIL\x1b[0m]";
        if (ok && not warning)  okfail = "[\x1b[32;1m OK \x1b[0m]";
        if (ok && warning)      okfail = "[\x1b[33;1mWARN\x1b[0m]";

        std::cout << okfail << std::endl;

        success &= ok;
    }

    return success;
}


enum test_type
{
    TEST_VERIFY_CONVERGENCE,
    TEST_MEASURE_TIMES
};

void test_triangles_specialized(test_type tt, const problem_type problem)
{
    size_t runs = 2;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/2D_quads/medit/square_h01.medit2d");
    paths.push_back("../../../diskpp/meshes/2D_quads/medit/square_h005.medit2d");
    paths.push_back("../../../diskpp/meshes/2D_quads/medit/square_h0025.medit2d");
    paths.push_back("../../../diskpp/meshes/2D_quads/medit/square_h00125.medit2d");

    typedef disk::generic_mesh<double, 2>      MT;
    typedef disk::medit_mesh_loader<double, 2>   LT;

    switch(tt)
    {
        case TEST_VERIFY_CONVERGENCE:
            verify_convergence<MT, LT>(paths, 1, 1, problem);  // ce 3 était 1 avec karol et même le 0
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}


int main(int argc, char **argv)
{
    test_type tt = TEST_VERIFY_CONVERGENCE;
    problem_type problem = TEST_DIRICHLET;

    int ch;

    while ( (ch = getopt(argc, argv, "tc:dnr")) != -1 )
    {
        switch(ch)
        {
            case 't':
                tt = TEST_MEASURE_TIMES;
                break;

            case 'c':
                tt = TEST_VERIFY_CONVERGENCE;
                break;
            case 'd':
                problem = TEST_DIRICHLET;
                break;
            case 'n':
                problem = TEST_NEUMANN;
                break;
            case 'r':
                problem = TEST_ROBIN;
                break;
            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);

        }
    }

    argc -= optind;
    argv += optind;


    std::cout << bold << underline << "Triangles specialized" << reset << std::endl;
    test_triangles_specialized(tt, problem);

}

#if 0
int main(int argc, char **argv)
{
    using T = double;

    if (argc != 2)
    {
        std::cout << "Please specify file name." << std::endl;
        return 1;
    }

    char *mesh_filename = argv[1];

     /* Medit 2d*/

    if (std::regex_match(mesh_filename, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        auto msh = disk::load_medit_2d_mesh<T>(mesh_filename);
        run_test_dirichlet(msh);
        return 0;
    }
};
#endif
