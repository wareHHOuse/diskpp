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
#include <optional>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/loaders/loader_gmsh.hpp"

#include "diskpp/cfem/cfem.hpp"
#include "diskpp/output/silo.hpp"
#include "diskpp/mesh/meshgen.hpp"

#include "cfem_heat.hpp"
#include "cfem_poisson.hpp"





int main_x(int argc, char **argv)
{
    using RealType = double;

    if (argc < 2) {
        return 1;
    }

    using mesh_type = disk::simplicial_mesh<RealType, 2>;

    disk::gmsh_geometry_loader<mesh_type> loader;
    loader.read_mesh(argv[1]);


    disk::cfem::heat::solver_state<RealType> state;
    state.dt = 1e-3;

    loader.populate_mesh(state.msh);

    std::cout << "init & asm\n"; 

    disk::cfem::heat::init(state);
    disk::cfem::heat::assemble(state);

    using dvec = disk::dynamic_vector<RealType>;
    dvec u_curr = dvec::Zero(state.msh.points_size());
    dvec u_next = dvec::Zero(state.msh.points_size());

    for (size_t i = 0; i < state.dirichlet_values.size(); i++) {
        if (state.dirichlet_values[i]) {
            u_curr[i] = *state.dirichlet_values[i];
        }
    }


    for (size_t ts = 0; ts < 10000; ts++) {
        if (ts % 10 == 0) {
            std::cout << "ts " << ts << "\n";
            disk::silo_database silo_db;
            std::string fname = "ts_" + std::to_string(ts) + ".silo";
            silo_db.create(fname);
            silo_db.add_mesh(state.msh, "mesh");
            silo_db.add_variable("mesh", "u", u_curr, disk::nodal_variable_t);
        }

        u_next = disk::cfem::heat::timestep(state, u_curr);

        u_curr = u_next;
    }


    return 0;
}



int main(int argc, char **argv)
{
    using RealType = double;

    if (argc < 2) {
        return 1;
    }

    using mesh_type = disk::simplicial_mesh<RealType, 2>;
    using point_type = typename mesh_type::point_type;

    


    disk::cfem::poisson::solver_state<RealType> state;
    
    //disk::gmsh_geometry_loader<mesh_type> loader;
    //loader.read_mesh(argv[1]);
    //loader.populate_mesh(state.msh);

    auto mesher = make_simple_mesher(state.msh);
    mesher.refine();
    mesher.refine();
    mesher.refine();
    mesher.refine();

    double theta = M_PI/2.0;

    state.msh.transform( [&](const point_type& pt){
        auto c = std::cos(theta);
        auto s = std::sin(theta);
        point_type newp{
            c*pt.x() - s*pt.y() + 1.0,
            s*pt.x() + c*pt.y() + 0.0
        };
        return newp;
    });


    std::cout << "init & asm\n"; 

    disk::cfem::poisson::init(state);

    /*
    disk::cfem::poisson::assemble(state);
    disk::cfem::poisson::solve(state);
    */

    disk::cfem::poisson::picard(state);

    disk::cfem::poisson::postpro(state);

    return 0;
}