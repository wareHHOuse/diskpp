/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019, 2024                nicolas.pignet@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework.
 * M. Abbas, A. Ern, N. Pignet.
 * International Journal of Numerical Methods in Engineering (2019)
 * 120(3), 303-327
 * DOI: 10.1002/nme.6137
 */

#include "../share/tests_data.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/mechanics/NewtonSolver/NonLinearSolver.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>

#include <unistd.h>

template < template < typename, size_t, typename > class Mesh, typename T, size_t N,
           typename Storage >
void run_nl_solid_mechanics_solver( const Mesh< T, N, Storage > &msh,
                                    const disk::mechanics::NonLinearParameters< T > &rp,
                                    const STUDY &study ) {
    typedef Mesh< T, N, Storage > mesh_type;

    /* Get material parameters */
    const auto material_data = getMaterialData< T >( study );

    /* Get boundary conditions */
    const auto bnd = getBoundaryConditions( msh, material_data, study );

    /* Create nonlinear solver */
    disk::mechanics::NonLinearSolver< mesh_type > nl( msh, bnd, rp );

    /* Get external load */
    addExternalLoad( msh, material_data, study, nl );

    /* Add non linear option */
    addNonLinearOptions( msh, material_data, study, nl );

    if ( nl.verbose() ) {
        std::cout << "Solving the problem ..." << '\n';
    }

    disk::mechanics::SolverInfo solve_info = nl.compute();

    if ( nl.verbose() ) {
        solve_info.printInfo();
    }

    if ( nl.convergence() ) {
        std::cout << "average diameter h: " << average_diameter( msh ) << std::endl;
    }
}

int main( int argc, char **argv ) {
    using RealType = double;

    char *mesh_filename = nullptr;

    disk::mechanics::NonLinearParameters< RealType > rp;

    int ch;

    while ( ( ch = getopt( argc, argv, "r:" ) ) != -1 ) {
        switch ( ch ) {
        case 'r':
            if ( !rp.readParameters( optarg ) )
                exit( 1 );
            break;
        default:
            std::cout << "wrong arguments" << std::endl;
            exit( 1 );
        }
    }

    argc -= optind;
    argv += optind;

    if ( argc == 0 ) {
        std::cout << "Error" << std::endl;
        return 0;
    }

    mesh_filename = argv[0];

    /* Define study parameters to use */
    const STUDY study = STUDY::IMPACT_2D;

    addAdditionalParameters( study, rp );

    /* Poly 2d*/
    if ( std::regex_match( mesh_filename, std::regex( ".*\\.poly2d$" ) ) ) {
        std::cout << "Guessed mesh format: Poly2D format" << std::endl;
        disk::generic_mesh< RealType, 2 > msh;
        disk::load_mesh_poly< RealType >( mesh_filename, msh, rp.getVerbose() );
        run_nl_solid_mechanics_solver( msh, rp, study );
        return 0;
    }

    /* Poly 3d*/
    if ( std::regex_match( mesh_filename, std::regex( ".*\\.poly3d$" ) ) ) {
        std::cout << "Guessed mesh format: Poly3D format" << std::endl;
        disk::generic_mesh< RealType, 3 > msh;
        disk::load_mesh_poly< RealType >( mesh_filename, msh, rp.getVerbose() );
        run_nl_solid_mechanics_solver( msh, rp, study );
        return 0;
    }

    /* Medit 2d*/
    if ( std::regex_match( mesh_filename, std::regex( ".*\\.medit2d$" ) ) ) {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        disk::generic_mesh< RealType, 2 > msh;
        disk::load_mesh_medit< RealType >( mesh_filename, msh );
        run_nl_solid_mechanics_solver( msh, rp, study );
        return 0;
    }

    /* Medit 3d*/
    if ( std::regex_match( mesh_filename, std::regex( ".*\\.medit3d$" ) ) ) {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        disk::generic_mesh< RealType, 3 > msh;
        disk::load_mesh_medit< RealType >( mesh_filename, msh );
        run_nl_solid_mechanics_solver( msh, rp, study );
        return 0;
    }
}
