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

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "diffusion_nitzche_solver.hpp"


template<typename MeshType, typename LoaderType>
void
test_mesh_format(const std::vector<std::string>& paths,
                 size_t runs, size_t mindeg, size_t maxdeg,
                 const std::string& output_basename)
{}


template<typename MeshType, typename LoaderType>
bool
verify_convergence(const std::vector<std::string>& paths,
                   size_t mindeg, size_t maxdeg,
                    algorithm_parameters<typename MeshType::scalar_type>& ap)
{
    typedef typename MeshType::scalar_type scalar_type;


    bool success = true;

    for (size_t i = mindeg; i <= maxdeg; i++)
    {
        scalar_type expected_rate = i+1;

        std::vector<std::pair<scalar_type, scalar_type>> errdiams;

        std::cout << "Convergence rates for k = " << i << ":   " << std::endl;

        ap.degree = i;

        for (auto& tsp : paths)
        {
            MeshType    msh;
            LoaderType  loader;

            if (!loader.read_mesh(tsp))
            {
                std::cout << "Problem loading mesh." << std::endl;
                return false;
            }
            loader.populate_mesh(msh);
            auto diam = average_diameter(msh);
            auto error = run_diffusion_solver(msh, ap);
            errdiams.push_back( std::make_pair(diam, error) );
        }

        bool pass       = true;
        bool warning    = false;
        bool high, low, ok;
        for (size_t i = 1; i < errdiams.size(); i++)
        {
            auto d = log2(errdiams[i-1].first/errdiams[i].first);
            auto e = log2(errdiams[i-1].second/errdiams[i].second);
            auto rate = e/d;

            ok   = (std::abs(expected_rate - rate) < 0.4); /* Test passed */
            low  = ((expected_rate - rate) > 0.2); /* a bit too low, warn */
            high = ((rate - expected_rate) > 0.2); /* a bit too high, warn */

            if (low)    std::cout << magenta;
            if (high)   std::cout << cyan;
            std::cout << std::fixed << std::setprecision(3) << errdiams[i].first;
            std::cout << "  " << rate << "  "<<std::endl;
            if (low or high)
            {
                std::cout << reset;
                warning = true;
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

void test_triangles_specialized(test_type tt, algorithm_parameters<double>& ap)
{
    size_t runs = 2;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri01.mesh2d");
    paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri02.mesh2d");
    paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri03.mesh2d");
    paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri04.mesh2d");

    typedef disk::simplicial_mesh<double, 2>      MT;
    typedef disk::netgen_mesh_loader<double, 2>   LT;

    switch(tt)
    {
        case TEST_MEASURE_TIMES:
            test_mesh_format<MT, LT>(paths, runs, 0, 3, "triangle_spec");
            break;

        case TEST_VERIFY_CONVERGENCE:
            verify_convergence<MT, LT>(paths, 0, 3, ap);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}

void test_triangles_generic(test_type tt, algorithm_parameters<double>& ap)
{
    size_t runs = 2;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_1.typ1");
    paths.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_2.typ1");
    paths.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_3.typ1");
    paths.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_4.typ1");

    typedef disk::generic_mesh<double, 2>       MT;
    typedef disk::fvca5_mesh_loader<double, 2>  LT;

    switch(tt)
    {
        case TEST_MEASURE_TIMES:
            test_mesh_format<MT, LT>(paths, runs, 0, 3, "triangle_gen");
            break;

        case TEST_VERIFY_CONVERGENCE:
            verify_convergence<MT, LT>(paths, 0, 3, ap);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}

void test_hexagons_generic(test_type tt, algorithm_parameters<double>& ap)
{
    size_t runs = 2;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_2.typ1");
    paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_3.typ1");
    paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_4.typ1");
    paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_5.typ1");

    typedef disk::generic_mesh<double, 2>       MT;
    typedef disk::fvca5_mesh_loader<double, 2>  LT;

    switch(tt)
    {
        case TEST_MEASURE_TIMES:
            test_mesh_format<MT, LT>(paths, runs, 0, 3, "hexagons_gen");
            break;

        case TEST_VERIFY_CONVERGENCE:
            verify_convergence<MT, LT>(paths, 0, 3, ap);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}

void test_hextri_generic(test_type tt, algorithm_parameters<double>& ap)
{
    size_t runs = 2;

    std::vector<std::string> paths;
    //paths.push_back("../hexagon_splitter/hextri1.typ1");
    paths.push_back("../hexagon_splitter/hextri2.typ1");
    paths.push_back("../hexagon_splitter/hextri3.typ1");
    paths.push_back("../hexagon_splitter/hextri4.typ1");
    paths.push_back("../hexagon_splitter/hextri5.typ1");


    typedef disk::generic_mesh<double, 2>       MT;
    typedef disk::fvca5_mesh_loader<double, 2>  LT;

    switch(tt)
    {
        case TEST_MEASURE_TIMES:
            test_mesh_format<MT, LT>(paths, runs, 0, 3, "hextri_gen");
            break;

        case TEST_VERIFY_CONVERGENCE:
            verify_convergence<MT, LT>(paths, 0, 3, ap);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}

void test_kershaw_2d(test_type tt, algorithm_parameters<double>& ap)
{
    size_t runs = 2;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_1.typ1");
    paths.push_back("../../../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_2.typ1");
    paths.push_back("../../../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_3.typ1");
    paths.push_back("../../../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_4.typ1");
    paths.push_back("../../../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_5.typ1");

    typedef disk::generic_mesh<double, 2>       MT;
    typedef disk::fvca5_mesh_loader<double, 2>  LT;

    switch(tt)
    {
        case TEST_MEASURE_TIMES:
            test_mesh_format<MT, LT>(paths, runs, 0, 3, "kershaw_2d");
            break;

        case TEST_VERIFY_CONVERGENCE:
            verify_convergence<MT, LT>(paths, 0, 3, ap);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}


int main(int argc, char **argv)
{
    test_type tt = TEST_VERIFY_CONVERGENCE;
    int ch;
    algorithm_parameters<double> ap;

    while ( (ch = getopt(argc, argv, "g:npzfcotv")) != -1 )
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
            case 'o':
                ap.solver = EVAL_IN_CELLS_AS_FACES;
                break;
            case 't':
                tt = TEST_MEASURE_TIMES;
                break;
            case 'v':
                tt = TEST_VERIFY_CONVERGENCE;
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
    test_triangles_specialized(tt, ap);

    std::cout << bold << underline << "Triangles generic" << reset << std::endl;
    test_triangles_generic(tt, ap);

    std::cout << bold << underline << "Hexagons" << reset << std::endl;
    test_hexagons_generic(tt, ap);

    std::cout << bold << underline << "Kershaw 2D" << reset << std::endl;
    test_kershaw_2d(tt, ap);

}
