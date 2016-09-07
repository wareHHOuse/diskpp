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
 
#include "util.h"

unsigned int fact(unsigned int n)
{
    return (n < 2) ? 1 : n*fact(n-1);
}

unsigned int binomial(unsigned int n, unsigned int k)
{
    if (n < k)
        return fact(n) / fact(k);

    return fact(n)/(fact(n-k)*fact(k));
}
