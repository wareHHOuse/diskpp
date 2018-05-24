/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
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
#include <sstream>
#include <iomanip>
#include <regex>
#include <type_traits>
#include <cassert>

#include <Eigen/Eigenvalues>

#include "loaders/loader.hpp"
#include "cfem/cfem.hpp"
#include "revolution/methods/hho"
#include "mesh/mesh_hierarchy.hpp"

#include "output/silo.hpp"
#include "solvers/solver.hpp"

#include "contrib/sol2/sol.hpp"
#include "contrib/timecounter.h"

template<typename T>
class hierarchical_eigval_solver
{
	typedef disk::simplicial_mesh<T, 2>  mesh_type;

	mesh_type					initial_mesh;
	disk::mesh_hierarchy<T> 	mesh_hier;

	bool init(const std::string& config_fn)
	{
		sol::state lua;
    	lua.open_libraries(sol::lib::math, sol::lib::base);
    	lua["config"] = lua.create_table();
    	disk::solvers::init_lua(lua);
    	lua["solver"]["feast"] = lua.create_table();

    	auto r = lua.do_file(config_fn);

    	if ( !r.valid() )
    	{
        	std::cout << "Problems opening configuration file" << std::endl;
        	return false;
    	}

    	std::string     input_mesh  = lua["config"]["input_mesh"];
    	bool            verbose     = lua["config"]["verbose"].get_or(false);

	    if (std::regex_match(input_mesh, std::regex(".*\\.mesh2d$") ))
	    {
	        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;

	        

	        mesh_type msh;
	        disk::netgen_mesh_loader<T, 2> loader;
	        loader.verbose(verbose);
	        if (!loader.read_mesh(input_mesh))
	        {
	            std::cout << "Problem loading mesh." << std::endl;
	            return false;
	        }

	        loader.populate_mesh(msh);

	        mesh_hier.build_hierarchy(msh, 6);

	        std::cout << "Mesh avg. diameter: " << mesh_h(msh) << std::endl;
	    }

	    return true;
	}

public:
	hierarchical_eigval_solver()
	{}

	hierarchical_eigval_solver(const std::string& config_fn)
	{
		init(config_fn);
	}

	void
	load_config(const std::string& config_fn)
	{
		init(config_fn);
	}

};

int main(int argc, char **argv)
{
    using scalar_type = double;

    if (argc != 2)
    {
        std::cout << "Please specify configuration file" << std::endl;
        return 1;
    }

    hierarchical_eigval_solver<scalar_type> hes(argv[1]);

}