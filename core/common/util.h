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

#pragma once

#include <cstddef>

unsigned int        fact(unsigned int);
unsigned int        binomial(unsigned int, unsigned int);

template<typename T>
T iexp_pow(T x, size_t n)
{
    if (n == 0)
        return 1;

    T y = 1;
    while (n > 1)
    {
        if (n % 2 == 0)
        {
            x = x * x;
            n = n / 2;
        }
        else
        {
            y = x * y;
            x = x * x;
            n = (n-1)/2;
        }
    }

    return x*y;
}
