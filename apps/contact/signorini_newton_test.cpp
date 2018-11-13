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

 #include "signorini_newton_solver.hpp"


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
                    algorithm_parameters<typename MeshType::coordinate_type>& ap,
                    const typename MeshType::coordinate_type& eta)
{
    typedef typename MeshType::coordinate_type coordinate_type;
    typedef std::pair<coordinate_type, coordinate_type> tuple_type;

    bool success = true;

    for (size_t i = mindeg; i <= maxdeg; i++)
    {
        coordinate_type expected_rate = i+1;


        std::vector<std::pair<coordinate_type, tuple_type>> errdiams;

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
            auto run_exact  = true;
            auto error_full = run_signorini(msh, ap, eta, run_exact);

            errdiams.push_back( std::make_pair(diam, error_full) );
        }

        bool pass       = true;
        bool warning    = false;
        bool high, low, ok;

        {
            auto error_full_i = errdiams[0].second;

            auto H1_error = error_full_i.first;
            auto L2_error = error_full_i.second;

            std::cout << std::fixed << std::setprecision(3) << errdiams[i].first <<" :    ";
            std::cout << " " <<  std::scientific<< std::setprecision(3)<< H1_error;
            std::cout << "  "<<  std::fixed << std::setprecision(3)<<"        -   ";
            std::cout << " " <<  std::scientific<< std::setprecision(3)<< L2_error <<std::endl;
        }

        for (size_t i = 1; i < errdiams.size(); i++)
        {

            auto error_full_i = errdiams[i].second;
            auto error_full_i_1 = errdiams[i-1].second;

            auto H1_e = log2(error_full_i_1.first/error_full_i.first);
            auto L2_e = log2(error_full_i_1.second/error_full_i.second);
            auto d = log2(errdiams[i-1].first/errdiams[i].first);

            auto H1_rate    = H1_e/d;
            auto L2_rate    = L2_e/d;

            ok   = (std::abs(expected_rate - H1_rate) < 0.4); /* Test passed */
            low  = ((expected_rate - H1_rate) > 0.2); /* a bit too low, warn */
            high = ((H1_rate - expected_rate) > 0.2); /* a bit too high, warn */

            if (low)    std::cout << magenta;
            if (high)   std::cout << cyan;


            auto H1_error = error_full_i.first;
            auto L2_error = error_full_i.second;

            std::cout << std::fixed << std::setprecision(3) << errdiams[i].first <<" :    ";
            std::cout << " " <<  std::scientific<< std::setprecision(3)<< H1_error;
            std::cout << "  "<<  std::fixed << std::setprecision(3)<< H1_rate << "   -   ";
            std::cout << " " <<  std::scientific<< std::setprecision(3)<< L2_error;
            std::cout << "  "<<  std::fixed << std::setprecision(3)<< L2_rate <<std::endl;
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

template<template<typename, size_t, typename> class Mesh,
      typename T, typename Storage>
void
renumber_boundaries(Mesh<T,2,Storage>& msh)
{
    auto storage = msh.backend_storage();

    for(size_t i = 0; i < msh.faces_size(); i++)
    {
        auto edge = *std::next(msh.faces_begin(), i);
        auto bar = barycenter(msh, edge);

        auto is_close_to = [](T val, T ref) -> bool {
            T eps = 1e-7;
            return std::abs(val - ref) < eps;
        };

        if (!storage->boundary_info.at(i).is_boundary)
            continue;

        size_t bid = 42;

        if ( is_close_to(bar.y(), 0.0) ) bid = 1;
        if ( is_close_to(bar.x(), 1.0) ) bid = 2;
        if ( is_close_to(bar.y(), -1.0) ) bid = 3;
        if ( is_close_to(bar.x(), -1.0) ) bid = 4;

        if (bid == 42)
            throw std::logic_error("Can not locate the edge");

        storage->boundary_info.at(i).boundary_id = bid;
    }
}
template<typename MeshType, typename LoaderType>
bool
verify_convergence_renum_bound(const std::vector<std::string>& paths,
                   size_t mindeg, size_t maxdeg,
                    algorithm_parameters<typename MeshType::coordinate_type>& ap,
                    const typename MeshType::coordinate_type& eta)
{
    typedef typename MeshType::coordinate_type coordinate_type;
    typedef std::pair<coordinate_type, coordinate_type> tuple_type;

    bool success = true;

    for (size_t i = mindeg; i <= maxdeg; i++)
    {
        coordinate_type expected_rate = i+1;


        std::vector<std::pair<coordinate_type, tuple_type>> errdiams;

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
            renumber_boundaries(msh);
            auto diam = average_diameter(msh);
            auto run_exact  = true;
            auto error_full = run_signorini(msh, ap, eta, run_exact);

            errdiams.push_back( std::make_pair(diam, error_full) );

        }

        bool pass       = true;
        bool warning    = false;
        bool high, low, ok;

        {
            auto error_full_i = errdiams[0].second;

            auto H1_error = error_full_i.first;
            auto L2_error = error_full_i.second;

            std::cout << std::fixed << std::setprecision(3) << errdiams[0].first <<" :    ";
            std::cout << " " <<  std::scientific<< std::setprecision(3)<< H1_error;
            std::cout << "  "<<  std::fixed << std::setprecision(3)<<"        -   ";
            std::cout << " " <<  std::scientific<< std::setprecision(3)<< L2_error <<std::endl;
        }

        for (size_t i = 1; i < errdiams.size(); i++)
        {

            auto error_full_i = errdiams[i].second;
            auto error_full_i_1 = errdiams[i-1].second;

            auto H1_e = log2(error_full_i_1.first/error_full_i.first);
            auto L2_e = log2(error_full_i_1.second/error_full_i.second);
            auto d = log2(errdiams[i-1].first/errdiams[i].first);

            auto H1_rate    = H1_e/d;
            auto L2_rate    = L2_e/d;

            ok   = (std::abs(expected_rate - H1_rate) < 0.4); /* Test passed */
            low  = ((expected_rate - H1_rate) > 0.2); /* a bit too low, warn */
            high = ((H1_rate - expected_rate) > 0.2); /* a bit too high, warn */

            if (low)    std::cout << magenta;
            if (high)   std::cout << cyan;


            auto H1_error = error_full_i.first;
            auto L2_error = error_full_i.second;

            std::cout << std::fixed << std::setprecision(3) << errdiams[i].first <<" :    ";
            std::cout << " " <<  std::scientific<< std::setprecision(3)<< H1_error;
            std::cout << "  "<<  std::fixed << std::setprecision(3)<< H1_rate << "   -   ";
            std::cout << " " <<  std::scientific<< std::setprecision(3)<< L2_error;
            std::cout << "  "<<  std::fixed << std::setprecision(3)<< L2_rate <<std::endl;
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

void test_triangles(test_type tt, algorithm_parameters<double>& ap, const double& eta)
{
    size_t runs = 2;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/square_tri1_big_neg.mesh2d");
    paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/square_tri2_big_neg.mesh2d");
    paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/square_tri3_big_neg.mesh2d");

    typedef disk::simplicial_mesh<double, 2>      MT;
    typedef disk::netgen_mesh_loader<double, 2>   LT;

    switch(tt)
    {
        case TEST_MEASURE_TIMES:
            test_mesh_format<MT, LT>(paths, runs, 0, 0, "triangle_spec");
            break;

        case TEST_VERIFY_CONVERGENCE:
            verify_convergence<MT, LT>(paths, 0, 3, ap, eta);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}

void test_triangles_specialized(test_type tt, algorithm_parameters<double>& ap, const double& eta)
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
            test_mesh_format<MT, LT>(paths, runs, 0, 0, "triangle_spec");
            break;

        case TEST_VERIFY_CONVERGENCE:
            verify_convergence<MT, LT>(paths, 0, 3, ap, eta);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}

void test_triangles_generic(test_type tt, algorithm_parameters<double>& ap, const double& eta)
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
            verify_convergence<MT, LT>(paths, 0, 3, ap, eta);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}

void test_quadrangles(test_type tt, algorithm_parameters<double>& ap, const double& eta)
{
    size_t runs = 2;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/2D_quads/diskpp/signorini_test_x16_y8.quad");
    paths.push_back("../../../diskpp/meshes/2D_quads/diskpp/signorini_test_x32_y16.quad");
    paths.push_back("../../../diskpp/meshes/2D_quads/diskpp/signorini_test_x64_y32.quad");
    paths.push_back("../../../diskpp/meshes/2D_quads/diskpp/signorini_test_x128_y64.quad");

    typedef disk::cartesian_mesh<double, 2>       MT;
    typedef disk::cartesian_mesh_loader<double, 2>  LT;

    switch(tt)
    {
        case TEST_MEASURE_TIMES:
            test_mesh_format<MT, LT>(paths, runs, 0, 3, "quadrangle_diskpp");
            break;

        case TEST_VERIFY_CONVERGENCE:
            verify_convergence_renum_bound<MT, LT>(paths, 0, 3, ap, eta);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}


void test_hexagons_generic(test_type tt, algorithm_parameters<double>& ap, const double& eta)
{
    size_t runs = 2;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_1_big_neg.typ1");
    paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_2_big_neg.typ1");
    paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_3_big_neg.typ1");
    paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_4_big_neg.typ1");
    //paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_5.typ1");

    typedef disk::generic_mesh<double, 2>       MT;
    typedef disk::fvca5_mesh_loader<double, 2>  LT;

    switch(tt)
    {
        case TEST_MEASURE_TIMES:
            test_mesh_format<MT, LT>(paths, runs, 0, 3, "hexagons_gen");
            break;

        case TEST_VERIFY_CONVERGENCE:
            verify_convergence_renum_bound<MT, LT>(paths, 0, 3, ap, eta);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}

void test_hextri_generic(test_type tt, algorithm_parameters<double>& ap, const double& eta)
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
            verify_convergence<MT, LT>(paths, 0, 3, ap, eta);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}

void test_kershaw_2d(test_type tt, algorithm_parameters<double>& ap, const double& eta)
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
            verify_convergence<MT, LT>(paths, 0, 3, ap, eta);
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
    double parameter = 1.;

    while ( (ch = getopt(argc, argv, "g:npzfcle:tv")) != -1 )
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
            case 'l':
                ap.solver = EVAL_IN_CELLS_FULL;
                break;
            case 'e':
                ap.solver = EVAL_WITH_PARAMETER;
                std::cout << "choosing parameter" << std::endl;
                parameter = atof(optarg);
                if (parameter < 0  || parameter > 1 )
                {
                    std::cout << "parameter must be in [0 1]. Falling back to 0.5" << std::endl;
                    parameter = 0.5;
                }
                std::cout << "parameter :"<< parameter << std::endl;
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

    std::cout << bold << underline << "Triangles for contact" << reset << std::endl;
    test_triangles(tt, ap, parameter);

    std::cout << bold << underline << "Quadrangles" << reset << std::endl;
    test_quadrangles(tt, ap, parameter);

    std::cout << bold << underline << "Hexagons" << reset << std::endl;
    test_hexagons_generic(tt, ap, parameter);


    #if 0
    std::cout << bold << underline << "Triangles specialized" << reset << std::endl;
    test_triangles_specialized(tt, ap, parameter);

    std::cout << bold << underline << "Hexagons" << reset << std::endl;
    test_hexagons_generic(tt, ap, parameter);

    std::cout << bold << underline << "Triangles generic" << reset << std::endl;
    test_triangles_generic(tt, ap, parameter);

    std::cout << bold << underline << "Kershaw 2D" << reset << std::endl;
    test_kershaw_2d(tt, ap, parameter);
    #endif
}
