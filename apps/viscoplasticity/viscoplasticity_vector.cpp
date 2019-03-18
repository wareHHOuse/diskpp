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
#include "quadratures/quadratures.hpp"
#include "methods/hho"

#include "core/loaders/loader.hpp"

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"
#include "output/postMesh.hpp"

#include "output/silo.hpp"

#include "viscoplasticity_vector_solver.hpp"

template<typename Mesh>
bool
run_bingham(const Mesh& msh, const hho_degree_info& hdi,
            bingham_data<typename Mesh::coordinate_type, vector_problem> & vp)
{
    using T = typename Mesh::coordinate_type;
    typedef typename Mesh::point_type      point_type;
    typedef disk::mechanics::BoundaryConditions<Mesh> boundary_type;

    boundary_type   bnd(msh);
    std::string     name;

    auto wall = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
        return Matrix<T, Mesh::dimension, 1>::Zero();
    };
    auto movingWall = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
        return Matrix<T, Mesh::dimension, 1>{1,0};
    };
    auto symmetryPlane = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
        return Matrix<T, Mesh::dimension, 1>{0,0};
    };

    T omega = 1;
    auto rotation = [&](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
       return Matrix<T, Mesh::dimension, 1>{omega * p.y(), -omega * p.x()};
    };

    //to run with cartesian
    auto velocity  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
        if( std::abs(p.y() - 1.) < 1.e-8 )
            return Matrix<T, Mesh::dimension, 1>{1,0};
        else
            return Matrix<T, Mesh::dimension, 1>{0,0};
    };

    std::cout << "Im HERE  1 " << std::endl;

    switch (vp.problem)
    {
        case DRIVEN:
        std::cout << " I'm in DRIVEN" << std::endl;

            name = "driven";
            /*------------------------------------------------------------------
            *  Check boundary labels for the unitary square domain
            *          Netgen     _____          Medit     _____
            *                4   |     | 2                |     |
            *                    |_____|                  |_____|
            *                       3                        2
            *-------------------------------------------------------------------*/
            #if 0
            bnd.addDirichletBC(disk::mechanics::DIRICHLET,1, movingWall); //TOP
            bnd.addDirichletBC(disk::mechanics::DIRICHLET,2, wall); //TOP
            bnd.addDirichletBC(disk::mechanics::DIRICHLET,3, wall); //TOP
            bnd.addDirichletBC(disk::mechanics::DIRICHLET,4, wall); //TOP
            #endif
            //------------------------------------------------------------------
            bnd.addDirichletEverywhere(velocity);
            //------------------------------------------------------------------
            vp.yield = std::sqrt(2) * vp.Bn;

            break;
        case VANE:
            name = "vane";

            std::cout << " I'm in VANE" << std::endl;

            bnd.addDirichletBC( 0, 2, wall);
            bnd.addDirichletBC( 0, 1, rotation);
            bnd.addNeumannBC(10, 3, symmetryPlane);

            vp.yield = std::sqrt(2) * vp.Bn *  (vp.mu * omega);
            break;

        default:
            std::cout << "wrong arguments" << std::endl;
            exit(1);
            break;
    }

    ADMM<Mesh> admm_bingham(msh, hdi, vp);
    return admm_bingham.run(msh, bnd);
}


template<typename T, typename ProblemType>
std::tuple<std::string, hho_degree_info, bingham_data<T, ProblemType>>
read_data(const std::string& config_fn )
{
    bingham_data<T, ProblemType> vp_init, vp;

    sol::state lua;
    lua.open_libraries(sol::lib::math, sol::lib::base);
    lua["config"] = lua.create_table();
    lua["bi"]   = lua.create_table();

    auto r = lua.do_file(config_fn);
    if ( !r.valid() )
        throw std::invalid_argument("Problems opening configuration file");

    // Degree
    size_t hho_degree_cell          = lua["config"]["degree_cell"].get_or(0);
    size_t hho_degree_face          = lua["config"]["degree_face"].get_or(0);
    hho_degree_info hdi(hho_degree_cell, hho_degree_face);

    //Mesh
    std::string input_mesh  = lua["config"]["input_mesh"];
    std::string hname = lua["bi"]["hname"];

    std::string alpha = lua["bi"]["alpha"];
    std::string Bn    = lua["bi"]["Bn"];

    vp.alpha = atof(alpha.c_str());
    vp.Lref  = lua["bi"]["Lref"].get_or(vp_init.Lref);
    vp.Vref  = lua["bi"]["Vref"].get_or(vp_init.Vref);
    vp.mu    = lua["bi"]["mu"].get_or(vp_init.mu);
    vp.Bn    = atof(Bn.c_str());//lua["bi"]["Bn"].get_or(vp_init.Bn);
    vp.f     = lua["bi"]["f"].get_or(vp_init.f);

    std::string name = lua["bi"]["problem"];
    if( name == "DRIVEN")   { vp.problem  =  DRIVEN;}
    else if(name == "VANE") { vp.problem  =  VANE; }
    else
    {
        std::cout << "Problem must be {DRIVEN, VANE}//, COUETTE}.";
        std::cout << " Falling back to circular" << std::endl;
        vp.problem =  DRIVEN;
    }

    std::cout << vp << std::endl;
    std::cout << "DEGREE = ( "<< hdi.cell_degree() << " , ";
    std::cout << hdi.face_degree() <<" )" << std::endl;

    vp.info = name + "_k" + tostr(hdi.cell_degree())  + "_B" + tostr(vp.Bn)
                    + "_h"+ hname + "_a" + tostr(vp.alpha);

    return std::make_tuple(input_mesh, hdi, vp );
}

int main(int argc, char **argv)
{
    using RealType = double;

    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    if (argc != 2)
    {
        std::cout << "Please specify configuration file" << std::endl;
        return 1;
    }

    auto rd = read_data<RealType, vector_problem>(argv[1]);
    auto vp  = std::get<2>(rd);
    auto hdi = std::get<1>(rd);
    auto input_mesh =std::get<0>(rd);

    /* Medit 2d*/
    if (std::regex_match(input_mesh, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        typedef disk::generic_mesh<RealType, 2>  mesh_type;
        //mesh_type msh = disk::load_medit_2d_mesh<RealType>(input_mesh);

        //#if 0
        mesh_type  msh;
        disk::medit_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(input_mesh))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 0;
        }
        loader.populate_mesh(msh);
        //#endif

        if(!run_bingham(msh, hdi, vp))
            std::cout << "No convergence" << std::endl;

        return 1;
    }

    /* Netgen */
    if (std::regex_match(input_mesh, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;

        typedef disk::simplicial_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::netgen_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(input_mesh))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 0;
        }
        loader.populate_mesh(msh);

        if(!run_bingham(msh, hdi, vp))
            std::cout << "No convergence" << std::endl;

        return 1;
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(input_mesh, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;

        typedef disk::cartesian_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::cartesian_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(input_mesh))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 0;
        }
        loader.populate_mesh(msh);

        if(!run_bingham(msh, hdi, vp))
            std::cout << "No convergence" << std::endl;

        return 1;
    }
}
//#endif


#if 0
int main(int argc, char **argv)
{
    using RealType = double;

    char    *word      = nullptr;
    int ch;
    size_t degree = 1;
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
                break;
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
#endif
