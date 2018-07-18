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

#include "revolution/bases"
#include "revolution/quadratures"
#include "revolution/methods/hho"

#include "core/loaders/loader.hpp"

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"
#include "output/postMesh.hpp"

#include "output/silo.hpp"

#include "viscoplasticity_scalar_solver.hpp"


auto
run_viscoplasticity(const std::string& config_fn)
{
    using T = double;

    viscoplasticity_data<T> vp, vp_init;

    sol::state lua;
    lua.open_libraries(sol::lib::math, sol::lib::base);
    lua["config"] = lua.create_table();
    lua["vpst"]    = lua.create_table();

    auto r = lua.do_file(config_fn);
	if ( !r.valid() )
	{
    	std::cout << "Problems opening configuration file" << std::endl;
    	return false;
	}

    // Degree
    size_t hho_degree_cell          = lua["config"]["degree_cell"].get_or(0);
    size_t hho_degree_face          = lua["config"]["degree_face"].get_or(0);
    hho_degree_info hdi(hho_degree_cell, hho_degree_face);

    //Mesh
    std::string     input_mesh  = lua["config"]["input_mesh"];


    /* Medit 2d*/
    std::cout << "Guessed mesh format: Medit format" << std::endl;
    typedef disk::generic_mesh<T, 2>  mesh_type;


    vp.alpha = lua["vpst"]["alpha"].get_or(vp_init.alpha);
    vp.Lref  = lua["vpst"]["Lref"].get_or(vp_init.Lref);
    vp.Vref  = lua["vpst"]["Vref"].get_or(vp_init.Vref);
    vp.mu    = lua["vpst"]["mu"].get_or(vp_init.mu);
    vp.Bn    = lua["vpst"]["Bn"].get_or(vp_init.Bn);
    vp.f     = lua["vpst"]["f"].get_or(vp_init.f);

    std::string name = lua["vpst"]["problem"];
    if( name == "circular")    { vp.problem  =  CIRCULAR;}
    else if(name == "annulus") { vp.problem  =  ANNULUS; }
    else
    {
        std::cout << "Problem must be {circular, annulus}. Falling back to circular" << std::endl;
        vp.problem = CIRCULAR;
    }

    /* Medit 2d*/
    if (std::regex_match(input_mesh, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        typedef disk::generic_mesh<T, 2>  mesh_type;
        //mesh_type msh = disk::load_medit_2d_mesh<T>(input_mesh);

        //#if 0
        mesh_type                     msh;
        disk::medit_mesh_loader<T, 2> loader;
        loader.read_mesh(input_mesh);
        loader.populate_mesh(msh);
        //#endif

        augmented_lagrangian_viscoplasticity<mesh_type> als(msh, hdi, vp);
        if(!als.run_alg(msh))
            std::cout << "No convergence" << std::endl;

        return true;
    }

    return false;
}


int main(int argc, char **argv)
{
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    using scalar_type = double;

    if (argc != 2)
    {
        std::cout << "Please specify configuration file" << std::endl;
        return 1;
    }

    run_viscoplasticity(argv[1]);

    return 0;
}
//#endif
