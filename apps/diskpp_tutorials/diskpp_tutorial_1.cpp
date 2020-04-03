/*
 *       /\        Matteo Cicuttin (C) 2016-2020
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

#include "core/mesh/mesh.hpp"
#include "core/loaders/loader.hpp"

/* DiSk++ tutorial: demonstrate how to iterate over mesh elements
 */

int main(int argc, const char **argv)
{
    using disk::cells;
    using disk::faces;

    if (argc < 2)
    {
        std::cout << "Please specify a DiSk++ cartesian mesh filename" << std::endl;
        return 1;
    }

    using T = double;
    disk::cartesian_mesh<T, 2> msh;

    /* Load the mesh */
    bool success = disk::load_cartesian_2d_mesh<T>(argv[1], msh);
    if (!success)
        return;

    /* Iterate on cells */
    for (auto& cl : cells(msh))
        std::cout << cl << std::endl;

    /* iterate on faces */
    for (auto& fc : faces(msh))
        std::cout << fc << std::endl;

    return 0;
}


