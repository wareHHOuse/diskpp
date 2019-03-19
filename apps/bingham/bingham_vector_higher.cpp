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

#include "bingham_vector_higher_solver.hpp"

template<typename Mesh>
bool
run_bingham(Mesh& msh, const hho_degree_info& hdi,
            bingham_data<typename Mesh::coordinate_type, vector_problem> & vp)
{
    using T = typename Mesh::coordinate_type;
    typedef typename Mesh::point_type      point_type;
    typedef disk::BoundaryConditions<Mesh, false> boundary_type;

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

    switch (vp.problem)
    {
        case DRIVEN:

            name = "driven";
            //------------------------------------------------------------------
            renumber_boundaries(msh);
            bnd.addDirichletBC(disk::DIRICHLET,1, movingWall);
            bnd.addDirichletBC(disk::DIRICHLET,2, wall);
            bnd.addDirichletBC(disk::DIRICHLET,3, wall);
            bnd.addDirichletBC(disk::DIRICHLET,4, wall);
            //------------------------------------------------------------------
            vp.yield = std::sqrt(2) * vp.Bn;
            break;

        case VANE:

            name = "vane";
            //------------------------------------------------------------------
            bnd.addDirichletBC( 0, 2, wall);
            bnd.addDirichletBC( 0, 1, rotation);
            bnd.addNeumannBC(10, 3, symmetryPlane);
            //------------------------------------------------------------------
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
    using T = double;

    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    if (argc != 2)
    {
        std::cout << "Please specify configuration file" << std::endl;
        return 1;
    }

    auto rd = read_data<T, vector_problem>(argv[1]);
    auto vp  = std::get<2>(rd);
    auto hdi = std::get<1>(rd);
    auto input_mesh =std::get<0>(rd);

    /* Medit 2d*/
    if (std::regex_match(input_mesh, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        auto msh = disk::load_medit_2d_mesh<T>(input_mesh.c_str());
        if(!run_bingham(msh, hdi, vp))
            std::cout << "No convergence" << std::endl;

        return 1;
    }

    /* Netgen */
    if (std::regex_match(input_mesh, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<T>(input_mesh.c_str());
        if(!run_bingham(msh, hdi, vp))
            std::cout << "No convergence" << std::endl;

        return 1;
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(input_mesh, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        auto msh = disk::load_cartesian_2d_mesh<T>(input_mesh.c_str());
        if(!run_bingham(msh, hdi, vp))
            std::cout << "No convergence" << std::endl;

        return 1;
    }
}
//#endif
