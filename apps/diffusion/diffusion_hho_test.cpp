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

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "geometry/geometry.hpp"
#include "loaders/loader.hpp"
#include "revolution/methods/hho"
#include "solvers/solver.hpp"

/***************************************************************************/
/* RHS definition */
template<typename Mesh>
struct rhs_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        return 2.0 * M_PI * M_PI * sin_px * sin_py;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        auto sin_pz = std::sin(M_PI * pt.z());
        return 3.0 * M_PI * M_PI * sin_px * sin_py * sin_pz;
    }
};

template<typename Mesh>
auto make_rhs_function(const Mesh& msh)
{
    return rhs_functor<Mesh>();
}

/***************************************************************************/
/* Expected solution definition */
template<typename Mesh>
struct solution_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        return sin_px * sin_py;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        auto sin_pz = std::sin(M_PI * pt.z());
        return sin_px * sin_py * sin_pz;
    }
};

template<typename Mesh>
auto make_solution_function(const Mesh& msh)
{
    return solution_functor<Mesh>();
}

using namespace revolution;

template<typename Mesh>
auto
run_hho_diffusion_solver(const Mesh& msh, const size_t degree)
{
    using T = typename Mesh::scalar_type;

    hho_degree_info hdi(degree);

    auto rhs_fun = make_rhs_function(msh);
    auto sol_fun = make_solution_function(msh);

    auto assembler = make_diffusion_assembler(msh, hdi);

    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
        auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
        auto rhs    = make_rhs(msh, cl, cb, rhs_fun);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = gr.second + stab;
        auto sc     = diffusion_static_condensation_compute(msh, cl, hdi, A, rhs);
        assembler.assemble(msh, cl, sc.first, sc.second, sol_fun);
    }

    assembler.finalize();

    size_t systsz = assembler.LHS.rows();
    size_t nnz = assembler.LHS.nonZeros();

    dynamic_vector<T> sol = dynamic_vector<T>::Zero(systsz);

    disk::solvers::pardiso_params<T> pparams;
    pparams.report_factorization_Mflops = true;
    mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

    T error = 0.0;

    std::ofstream ofs("sol.dat");

    for (auto& cl : msh)
    {
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
        auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
        auto rhs    = make_rhs(msh, cl, cb, rhs_fun);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = gr.second + stab;

        Eigen::Matrix<T, Eigen::Dynamic, 1> locsol =
            assembler.take_local_data(msh, cl, sol, sol_fun);

        Eigen::Matrix<T, Eigen::Dynamic, 1> fullsol =
            diffusion_static_condensation_recover(msh, cl, hdi, A, rhs, locsol);

        Eigen::Matrix<T, Eigen::Dynamic, 1> realsol = project_function(msh, cl, hdi, sol_fun);

        /// Just for residual
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
        typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              vector_type;

        matrix_type ATT = A.block(0,0, cbs, cbs);
        vector_type u_cell = fullsol.block(0,0, cbs, 1);
        vector_type res = ATT * u_cell - rhs;
        matrix_type mm  = make_mass_matrix(msh, cl, cb, hdi.cell_degree());
        error += res.dot(mm * res);

        //auto diff = realsol - fullsol;
        //error += diff.dot(A*diff);

        auto bar = barycenter(msh, cl);

        for (size_t i = 0; i < Mesh::dimension; i++)
            ofs << bar[i] << " ";
        ofs << fullsol(0) << std::endl;

    }

    std::cout << std::sqrt(error) << std::endl;

    ofs.close();
    return std::sqrt(error);
}



template<typename MeshType, typename LoaderType>
void
test_mesh_format(const std::vector<std::string>& paths,
                 size_t runs, size_t mindeg, size_t maxdeg,
                 const std::string& output_basename)
{

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

            auto error = run_hho_diffusion_solver(msh, i);
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
    test_hexagons_generic(tt);

    std::cout << bold << underline << "Kershaw 2D" << reset << std::endl;
    test_kershaw_2d(tt);

    std::cout << bold << underline << "Hexahedra specialized" << reset << std::endl;
    test_hexahedra_specialized(tt);

    std::cout << bold << underline << "Hexahedra generic" << reset << std::endl;
    test_hexahedra_generic(tt);

    std::cout << bold << underline << "Tetrahedra specialized" << reset << std::endl;
    test_tetrahedra_specialized(tt);

    std::cout << bold << underline << "Tetrahedra generic" << reset << std::endl;
    test_tetrahedra_generic(tt);

    std::cout << bold << underline << "Polyhedra" << reset << std::endl;
    test_polyhedra_generic(tt);

    std::cout << bold << underline << "Hextri generic" << reset << std::endl;
    test_hextri_generic(tt);

}
