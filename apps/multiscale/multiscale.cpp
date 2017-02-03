/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016, 2017 - matteo.cicuttin@enpc.fr
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
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include <map>

#include "colormanip.h"

#include "../../config.h"

#include "loaders/loader.hpp"
#include "hho/hho_multiscale.hpp"






int main(int argc, char **argv)
{
    using RealType = double;

    if ( !std::regex_match(argv[1], std::regex(".*\\.mesh2d$") ) )
    {
        std::cout << "Mesh format must be Netgen 2D" << std::endl;
        return 1;
    }

    auto msh = disk::load_netgen_2d_mesh<RealType>(argv[1]);

    size_t k = 5;
    size_t rl = 2;
    disk::multiscale_local_problem<decltype(msh)> mlp(k, rl);
    for (auto& cl : msh)
    {
        mlp.assemble(msh, cl);
        break;
    }

    return 0;
}
