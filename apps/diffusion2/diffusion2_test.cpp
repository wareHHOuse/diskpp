/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016, 2017 - matteo.cicuttin@enpc.fr
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

#include "diffusion2_solver.hpp"


template<typename MeshType, typename LoaderType>
void
test_mesh_format(const std::vector<std::string>& paths,
                 size_t runs, size_t mindeg, size_t maxdeg,
                 const std::string& output_basename)
{

    auto f = [](const point<typename MeshType::scalar_type, MeshType::dimension>& p) -> auto {
        return M_PI * M_PI * sin(p.x() * M_PI);
    };

    auto sf = [](const point<typename MeshType::scalar_type, MeshType::dimension>& p) -> auto {
        return sin(p.x() * M_PI);
    };


    for (size_t deg = mindeg; deg <= maxdeg; deg++)
    {
        std::stringstream ss;
        ss << output_basename << "_degree_" << deg << ".txt";

        std::cout << bold << red;
        std::cout << "Timing vs. # of DoFs, k = " << deg << ", runs = " << runs;
        std::cout << " [ Saved in " << ss.str() << " ]" << reset << std::endl;
        std::cout << bold << "   DoFs        ";
        std::cout << yellow << "Rec        Stab        StatC      " << nocolor;
        std::cout << green << "TotAsm        Sol" << nocolor;
        std::cout << reset << std::endl;

        std::ofstream ofs(ss.str());

        size_t mesh_index = 1;
        for (auto& tsp : paths)
        {
            MeshType    msh;
            LoaderType  loader;

            if (!loader.read_mesh(tsp))
            {
                std::cout << "Problem loading mesh." << std::endl;
                return;
            }
            loader.populate_mesh(msh);

            double time_gradrec     = 0.0;
            double time_stab        = 0.0;
            double time_statcond    = 0.0;
            double time_solver      = 0.0;

            size_t dofs;
            for (size_t run = 0; run < runs; run++)
            {
                diffusion_solver<MeshType> dp(msh, deg);

                assembly_info ai = dp.assemble(f, sf);
                time_gradrec    += ai.time_gradrec;
                time_stab       += ai.time_stab;
                time_statcond   += ai.time_statcond;
                dofs = ai.linear_system_size;

                solver_info si = dp.solve();
                time_solver += si.time_solver;
            }

            time_gradrec    /= double(runs);
            time_stab       /= double(runs);
            time_statcond   /= double(runs);
            time_solver     /= double(runs);

            ofs << dofs << " ";
            ofs << time_gradrec << " ";
            ofs << time_stab << " ";
            ofs << time_statcond << " ";
            ofs << time_solver << std::endl;

            std::cout << " " << std::setw(8) << dofs << " | ";

            if (time_gradrec > time_stab) std::cout << yellow;
            std::cout << sci3 << time_gradrec << "   " << reset;

            if (time_stab > time_gradrec) std::cout << yellow;
            std::cout << sci3 << time_stab << "   " << reset;

            std::cout << nocolor;

            std::cout << sci3 << time_statcond << " | ";

            auto time_assembly = time_gradrec + time_stab + time_statcond;

            if (time_assembly > time_solver) std::cout << green;
            std::cout << sci3 << time_assembly << reset;

            if (time_assembly < time_solver) std::cout << green;
            std::cout << "   " << sci3 << time_solver << reset;

            std::cout << std::endl;

        }

        ofs.close();
    }
}

template<typename MeshType, typename LoaderType>
bool
verify_convergence(const std::vector<std::string>& paths,
                   size_t mindeg, size_t maxdeg)
{
    typedef typename MeshType::scalar_type scalar_type;

    auto f = [](const point<typename MeshType::scalar_type, MeshType::dimension>& p) -> auto {
        return 2.0 * M_PI * M_PI * sin(p.x() * M_PI) * sin(p.y() * M_PI);
    };

    auto sf = [](const point<typename MeshType::scalar_type, MeshType::dimension>& p) -> auto {
        return sin(p.x() * M_PI) * sin(p.y() * M_PI);
    };

    bool success = true;

    for (size_t i = mindeg; i <= maxdeg; i++)
    {
        scalar_type expected_rate = i+2;

        std::vector<std::pair<scalar_type, scalar_type>> errdiams;

        std::cout << "Convergence rates for k = " << i << ":   " << std::flush;

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

            diffusion_solver<MeshType> dp(msh, i);

            dp.assemble(f, sf);
            dp.solve();
            dp.postprocess(f);
            auto error = dp.compute_l2_error(sf);
            auto diam = average_diameter(msh);

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
            std::cout << std::fixed << std::setprecision(3) << rate << "  ";
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

void test_triangles_specialized(test_type tt)
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
            verify_convergence<MT, LT>(paths, 0, 3);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}

void test_triangles_generic(test_type tt)
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
            verify_convergence<MT, LT>(paths, 0, 3);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}

void test_hexagons_generic(test_type tt)
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
            verify_convergence<MT, LT>(paths, 0, 3);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}

void test_hextri_generic(test_type tt)
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
            verify_convergence<MT, LT>(paths, 0, 3);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}

void test_kershaw_2d(test_type tt)
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
            verify_convergence<MT, LT>(paths, 0, 3);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}


void test_hexahedra_specialized(test_type tt)
{
    size_t runs = 2;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-2-2-2.hex");
    paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-4-4-4.hex");
    paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-8-8-8.hex");
    paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-16-16-16.hex");
    //paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-32-32-32.hex");

    typedef disk::cartesian_mesh<double, 3>         MT;
    typedef disk::cartesian_mesh_loader<double, 3>  LT;

    switch(tt)
    {
        case TEST_MEASURE_TIMES:
            test_mesh_format<MT, LT>(paths, runs, 0, 3, "hexahedra_spec");
            break;

        case TEST_VERIFY_CONVERGENCE:
            verify_convergence<MT, LT>(paths, 0, 3);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}

void test_hexahedra_generic(test_type tt)
{
    size_t runs = 2;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_2x2x2.msh");
    paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_4x4x4.msh");
    paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_8x8x8.msh");
    paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_16x16x16.msh");
    //paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_32x32x32.hex");

    typedef disk::generic_mesh<double, 3>       MT;
    typedef disk::fvca6_mesh_loader<double, 3>  LT;

    switch(tt)
    {
        case TEST_MEASURE_TIMES:
            test_mesh_format<MT, LT>(paths, runs, 0, 3, "hexahedra_gen");
            break;

        case TEST_VERIFY_CONVERGENCE:
            verify_convergence<MT, LT>(paths, 0, 3);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}

void test_tetrahedra_specialized(test_type tt)
{
    size_t runs = 2;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/3D_tetras/netgen/fvca6_tet1.mesh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/netgen/fvca6_tet2.mesh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/netgen/fvca6_tet3.mesh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/netgen/fvca6_tet4.mesh");

    typedef disk::simplicial_mesh<double, 3>    MT;
    typedef disk::netgen_mesh_loader<double, 3> LT;

    switch(tt)
    {
        case TEST_MEASURE_TIMES:
            test_mesh_format<MT, LT>(paths, runs, 0, 3, "tetrahedra_spec");
            break;

        case TEST_VERIFY_CONVERGENCE:
            verify_convergence<MT, LT>(paths, 0, 3);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}

void test_tetrahedra_generic(test_type tt)
{
    size_t runs = 2;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/3D_tetras/fvca6/tet.1.msh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/fvca6/tet.2.msh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/fvca6/tet.3.msh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/fvca6/tet.4.msh");

    typedef disk::generic_mesh<double, 3>       MT;
    typedef disk::fvca6_mesh_loader<double, 3>  LT;

    switch(tt)
    {
        case TEST_MEASURE_TIMES:
            test_mesh_format<MT, LT>(paths, runs, 0, 3, "tetrahedra_gen");
            break;

        case TEST_VERIFY_CONVERGENCE:
            verify_convergence<MT, LT>(paths, 0, 3);
            break;

        default:
            std::cout << "[ Unavailable Test ]" << std::endl;
            return;
    }
}

void test_polyhedra_generic(test_type tt)
{
    size_t runs = 2;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_10.msh");
    paths.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_20.msh");
    paths.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_30.msh");
    //paths.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_40.msh");

    typedef disk::generic_mesh<double, 3>       MT;
    typedef disk::fvca6_mesh_loader<double, 3>  LT;

    switch(tt)
    {
        case TEST_MEASURE_TIMES:
            test_mesh_format<MT, LT>(paths, runs, 0, 3, "polyhedra");
            break;

        case TEST_VERIFY_CONVERGENCE:
            verify_convergence<MT, LT>(paths, 0, 3);
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

    while ( (ch = getopt(argc, argv, "tc")) != -1 )
    {
        switch(ch)
        {
            case 't':
                tt = TEST_MEASURE_TIMES;
                break;

            case 'c':
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
    test_triangles_specialized(tt);

    std::cout << bold << underline << "Triangles generic" << reset << std::endl;
    test_triangles_generic(tt);

    std::cout << bold << underline << "Hexagons" << reset << std::endl;
    //test_hexagons_generic(tt);

    std::cout << bold << underline << "Kershaw 2D" << reset << std::endl;
    //test_kershaw_2d(tt);

    std::cout << bold << underline << "Hexahedra specialized" << reset << std::endl;
    //test_hexahedra_specialized(tt);

    std::cout << bold << underline << "Hexahedra generic" << reset << std::endl;
    //test_hexahedra_generic(tt);

    std::cout << bold << underline << "Tetrahedra specialized" << reset << std::endl;
    test_tetrahedra_specialized(tt);

    std::cout << bold << underline << "Tetrahedra generic" << reset << std::endl;
    test_tetrahedra_generic(tt);

    std::cout << bold << underline << "Polyhedra" << reset << std::endl;
    //test_polyhedra_generic(tt);

    std::cout << bold << underline << "Hextri generic" << reset << std::endl;
    //test_hextri_generic(tt);

}
