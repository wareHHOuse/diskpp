/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
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
#include <list>

#include "geometry/geometry_generic.hpp"
#include "geometry/geometry_simplicial.hpp"
#include "geometry/geometry.hpp"
#include "loaders/loader.hpp"

int main(void)
{
    using T = double;

    using mesh_type = disk::simplicial_mesh<T,2>;

    mesh_type msh;
    disk::netgen_mesh_loader<T,2> loader;
    loader.read_mesh("../../../disk/meshes/square_submesh_test.mesh2d");
    loader.populate_mesh(msh);

    disk::submesher<mesh_type> sm;

    for (auto& cl : msh)
    {
        sm.generate_mesh(msh, cl);
    }

}
