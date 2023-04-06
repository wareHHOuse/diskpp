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

#include "diskpp/mesh/meshgen.hpp"

using namespace disk;

int main(int argc, char **argv)
{
    using T = double;

    triangular_mesh<T>  msh;
    auto mesher = make_simple_mesher(msh);

    for (size_t i = 0; i < 7; i++)
    {
        std::cout << "*****************" << std::endl;
        mesher.refine();
        std::cout << msh.cells_size() << std::endl;

        T area = 0.;
        for(auto& cl : msh)
            area += measure(msh, cl);

        std::cout << area << std::endl;
    }
}