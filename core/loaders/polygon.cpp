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


#include "polygon.hpp"

int main(void)
{
    polygon_store ps;
    ps.add_polygon(std::vector<size_t>({10,29,35,43}));
    ps.add_polygon(std::vector<size_t>({3,29,35}));

    for (auto p : ps)
        std::cout << p << std::endl;

    return 0;
}
