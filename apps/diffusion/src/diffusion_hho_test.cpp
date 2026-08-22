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

#include "diskpp/common/colormanip.h"
#include "diskpp/common/timecounter.hpp"
#include "diskpp/methods/implementation_hho/methods_hho.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <sstream>

#include <unistd.h>

#define _USE_MATH_DEFINES
#include "diskpp/geometry/geometry.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/hho"
#include "diskpp/output/silo.hpp"
#include "diskpp/solvers/direct_solvers.hpp"

/***************************************************************************/
/* RHS definition */
template<typename Mesh>
struct rhs_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
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
    typedef typename mesh_type::coordinate_type scalar_type;
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
    typedef typename mesh_type::coordinate_type scalar_type;
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
    typedef typename mesh_type::coordinate_type scalar_type;
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

using namespace disk;

template<typename Mesh>
auto
run_hho_diffusion_solver(const Mesh& msh, const size_t degree)
{
    using T = typename Mesh::coordinate_type;

    hho_degree_info hdi(degree);

    auto rhs_fun = make_rhs_function(msh);
    auto sol_fun = make_solution_function(msh);

    auto assembler = make_diffusion_assembler(msh, hdi);

    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
        auto stab = make_scalar_hho_stabilization( msh, cl, gr.first, hdi );
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = gr.second + stab;
        auto rhst = make_rhs( msh, cl, cb, rhs_fun );
        dynamic_vector< T > rhs = dynamic_vector< T >::Zero( A.cols() );
        rhs.head( cb.size() ) = rhst;
        auto [lhsC, rhsC] = disk::static_condensation( A, rhs, cb.size() );
        assembler.assemble( msh, cl, lhsC, rhsC, sol_fun );
    }

    assembler.finalize();

    size_t systsz = assembler.LHS.rows();
    size_t nnz = assembler.LHS.nonZeros();

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(systsz);

    disk::solvers::sparse_lu( assembler.LHS, assembler.RHS, sol );

    T errorL2 = 0.0;
    T errorEnergy = 0.0;

    std::vector< T > u;

    for (auto& cl : msh)
    {
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
        auto stab = make_scalar_hho_stabilization( msh, cl, gr.first, hdi );
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = gr.second + stab;
        auto rhst = make_rhs( msh, cl, cb, rhs_fun );
        dynamic_vector< T > rhs = dynamic_vector< T >::Zero( A.cols() );
        rhs.head( cb.size() ) = rhst;

        Eigen::Matrix<T, Eigen::Dynamic, 1> locsol =
            assembler.take_local_data(msh, cl, sol, sol_fun);

        Eigen::Matrix< T, Eigen::Dynamic, 1 > fullsol = static_decondensation( A, rhs, locsol );

        Eigen::Matrix< T, Eigen::Dynamic, 1 > realsol =
            project_function( msh, cl, hdi, sol_fun, 2 );

        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
        typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              vector_type;

        vector_type diffT = fullsol.head( cbs ) - realsol.head( cbs );
        matrix_type mm = make_mass_matrix( msh, cl, cb );
        errorL2 += diffT.dot( mm * diffT );

        auto diff = realsol - fullsol;
        errorEnergy += diff.dot( A * diff );

        auto bar = barycenter(msh, cl);

        u.push_back( fullsol( 0 ) );
    }

    std::stringstream ss;
    ss << "diffusion_hho_test_" << degree << ".silo";

    disk::silo_database silo;
    silo.create( ss.str() );
    silo.add_mesh( msh, "mesh" );
    silo.add_variable( "mesh", "u", u, disk::zonal_variable_t );

    std::cout << "L2 error: " << std::sqrt( errorL2 ) << ", A error: " << std::sqrt( errorEnergy )
              << std::endl;

    return std::sqrt( errorL2 );
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
    typedef typename MeshType::coordinate_type scalar_type;

    bool success = true;

    for (size_t i = mindeg; i <= maxdeg; i++)
    {
        scalar_type expected_rate = i+2;

        std::vector<std::pair<scalar_type, scalar_type>> errdiams;

        std::cout << "Convergence rates for k = " << i << ":   " << std::endl;

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
            auto d = std::log( errdiams[i - 1].first / errdiams[i].first );
            auto e = std::log( errdiams[i - 1].second / errdiams[i].second );
            auto rate = e/d;

            ok   = (std::abs(expected_rate - rate) < 0.4); /* Test passed */
            low  = ((expected_rate - rate) > 0.2); /* a bit too low, warn */
            high = ((rate - expected_rate) > 0.2); /* a bit too high, warn */

            if (low)    std::cout << magenta;
            if (high)   std::cout << cyan;
            std::ios_base::fmtflags f( std::cout.flags() );
            std::cout << std::fixed << std::setprecision(3) << rate << "  ";
            std::cout.flags( f );
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

    char *env_mesh_base = getenv( "DISKPP_MESH_PATH" );
    std::string mesh_base = "../../../diskpp/meshes/";
    if ( env_mesh_base )
        mesh_base = env_mesh_base;

    std::vector< std::string > paths;
    paths.push_back( mesh_base + "/2D_triangles/netgen/tri01.mesh2d" );
    paths.push_back( mesh_base + "/2D_triangles/netgen/tri02.mesh2d" );
    paths.push_back( mesh_base + "/2D_triangles/netgen/tri03.mesh2d" );
    paths.push_back( mesh_base + "/2D_triangles/netgen/tri04.mesh2d" );

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

    char *env_mesh_base = getenv( "DISKPP_MESH_PATH" );
    std::string mesh_base = "../../../diskpp/meshes/";
    if ( env_mesh_base )
        mesh_base = env_mesh_base;

    std::vector< std::string > paths;
    paths.push_back( mesh_base + "/2D_triangles/fvca5/mesh1_1.typ1" );
    paths.push_back( mesh_base + "/2D_triangles/fvca5/mesh1_2.typ1" );
    paths.push_back( mesh_base + "/2D_triangles/fvca5/mesh1_3.typ1" );
    paths.push_back( mesh_base + "/2D_triangles/fvca5/mesh1_4.typ1" );

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

    char *env_mesh_base = getenv( "DISKPP_MESH_PATH" );
    std::string mesh_base = "../../../diskpp/meshes/";
    if ( env_mesh_base )
        mesh_base = env_mesh_base;

    std::vector< std::string > paths;
    paths.push_back( mesh_base + "/2D_hex/fvca5/hexagonal_1.typ1" );
    paths.push_back( mesh_base + "/2D_hex/fvca5/hexagonal_2.typ1" );
    paths.push_back( mesh_base + "/2D_hex/fvca5/hexagonal_3.typ1" );
    paths.push_back( mesh_base + "/2D_hex/fvca5/hexagonal_4.typ1" );
    paths.push_back( mesh_base + "/2D_hex/fvca5/hexagonal_5.typ1" );

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

void test_kershaw_2d(test_type tt)
{
    size_t runs = 2;

    char *env_mesh_base = getenv( "DISKPP_MESH_PATH" );
    std::string mesh_base = "../../../diskpp/meshes/";
    if ( env_mesh_base )
        mesh_base = env_mesh_base;

    std::vector< std::string > paths;
    paths.push_back( mesh_base + "/2D_kershaw/fvca5/mesh4_1_1.typ1" );
    paths.push_back( mesh_base + "/2D_kershaw/fvca5/mesh4_1_2.typ1" );
    paths.push_back( mesh_base + "/2D_kershaw/fvca5/mesh4_1_3.typ1" );
    paths.push_back( mesh_base + "/2D_kershaw/fvca5/mesh4_1_4.typ1" );

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

    char *env_mesh_base = getenv( "DISKPP_MESH_PATH" );
    std::string mesh_base = "../../../diskpp/meshes/";
    if ( env_mesh_base )
        mesh_base = env_mesh_base;

    std::vector< std::string > paths;
    paths.push_back( mesh_base + "/3D_hexa/diskpp/testmesh-2-2-2.hex" );
    paths.push_back( mesh_base + "/3D_hexa/diskpp/testmesh-4-4-4.hex" );
    paths.push_back( mesh_base + "/3D_hexa/diskpp/testmesh-8-8-8.hex" );
    paths.push_back( mesh_base + "/3D_hexa/diskpp/testmesh-16-16-16.hex" );

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

    char *env_mesh_base = getenv( "DISKPP_MESH_PATH" );
    std::string mesh_base = "../../../diskpp/meshes/";
    if ( env_mesh_base )
        mesh_base = env_mesh_base;

    std::vector< std::string > paths;
    paths.push_back( mesh_base + "/3D_hexa/fvca6/hexa_2x2x2.msh" );
    paths.push_back( mesh_base + "/3D_hexa/fvca6/hexa_4x4x4.msh" );
    paths.push_back( mesh_base + "/3D_hexa/fvca6/hexa_8x8x8.msh" );
    paths.push_back( mesh_base + "/3D_hexa/fvca6/hexa_16x16x16.msh" );

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

    char *env_mesh_base = getenv( "DISKPP_MESH_PATH" );
    std::string mesh_base = "../../../diskpp/meshes/";
    if ( env_mesh_base )
        mesh_base = env_mesh_base;

    std::vector< std::string > paths;
    paths.push_back( mesh_base + "/3D_tetras/netgen/fvca6_tet0.mesh" );
    paths.push_back( mesh_base + "/3D_tetras/netgen/fvca6_tet1.mesh" );
    paths.push_back( mesh_base + "/3D_tetras/netgen/fvca6_tet2.mesh" );
    paths.push_back( mesh_base + "/3D_tetras/netgen/fvca6_tet3.mesh" );
    paths.push_back( mesh_base + "/3D_tetras/netgen/fvca6_tet4.mesh" );

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

    char *env_mesh_base = getenv( "DISKPP_MESH_PATH" );
    std::string mesh_base = "../../../diskpp/meshes/";
    if ( env_mesh_base )
        mesh_base = env_mesh_base;

    std::vector< std::string > paths;
    paths.push_back( mesh_base + "/3D_tetras/netgen/fvca6_tet0.mesh" );
    paths.push_back( mesh_base + "/3D_tetras/netgen/fvca6_tet1.mesh" );
    paths.push_back( mesh_base + "/3D_tetras/netgen/fvca6_tet2.mesh" );
    paths.push_back( mesh_base + "/3D_tetras/netgen/fvca6_tet3.mesh" );
    paths.push_back( mesh_base + "/3D_tetras/netgen/fvca6_tet4.mesh" );

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

    char *env_mesh_base = getenv( "DISKPP_MESH_PATH" );
    std::string mesh_base = "../../../diskpp/meshes/";
    if ( env_mesh_base )
        mesh_base = env_mesh_base;

    std::vector< std::string > paths;
    paths.push_back( mesh_base + "/3D_general/fvca6/dbls_10.msh" );
    paths.push_back( mesh_base + "/3D_general/fvca6/dbls_20.msh" );
    paths.push_back( mesh_base + "/3D_general/fvca6/dbls_30.msh" );

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

    return 0;
}
