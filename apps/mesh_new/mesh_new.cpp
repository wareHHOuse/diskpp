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
#include <iomanip>
#include <regex>

#include "core/mesh/uniform_mesher.hpp"

#include "core/output/hdf5_io.hpp"


int main(int argc, char **argv)
{
    hdf5_context hctx("test.h5");

    typedef disk::mesh_v2::triangular_mesh<double>      mesh_type;
    typedef typename mesh_type::point_type              point_type;

    disk::uniform_mesher<mesh_type> um;



    mesh_type m = um.convert_to_user();

    dump_to_matlab(m, "user.m");

    save(hctx, m, "test");

}
