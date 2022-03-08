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
#include <list>
#include <vector>
#include <sstream>

#include "geometry/geometry.hpp"
#include "loaders/loader.hpp"
#include "mesh/mesh_utility.hpp"
#include "mesh/mesh_hierarchy.hpp"


int main(int argc, char **argv)
{
    /*
    using T = double;
    using mesh_type = disk::generic_mesh<T,2>;

    mesh_type msh;
    disk::fvca5_mesh_loader<T,2> loader;
    loader.read_mesh("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_1.typ1");
    loader.populate_mesh(msh);

    dump_to_matlab(msh, "basemesh.m");

    disk::submesher<mesh_type> sm;
    size_t count = 0;
    for (auto& cl : msh)
    {
        auto submesh = sm.generate_mesh(msh, cl, 2);
        relax_mesh(submesh);

        std::stringstream ss;
        ss << "element_" << count << ".m";
        dump_to_matlab(submesh, ss.str());
        count++;
    }
    */

    using T = double;
    using mesh_type = disk::simplicial_mesh<T,2>;

    mesh_type msh;
    disk::netgen_mesh_loader<T,2> loader;

    loader.read_mesh("/Users/matteo/Desktop/lshape_1_1.mesh2d");
    loader.populate_mesh(msh);

    disk::mesh_hierarchy<T> mh(msh, 6);

    size_t i = 2;
    for (auto itor = mh.meshes_begin(); itor != mh.meshes_end(); itor++)
    {
        std::stringstream ss;
        ss << "/Users/matteo/Desktop/lshape_" << i << ".mesh2d";
        dump_netgen_format(*itor, ss.str());
        i++;
    }

    return 0;

}
