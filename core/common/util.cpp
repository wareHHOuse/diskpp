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


rusage_monitor::rusage_monitor()
{}

rusage_monitor::~rusage_monitor()
{
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);

    double u_secs = ru.ru_utime.tv_sec + ru.ru_utime.tv_usec/1e6;
    double s_secs = ru.ru_stime.tv_sec + ru.ru_stime.tv_usec/1e6;

    std::cout << "User CPU time:   " << u_secs << " seconds" << std::endl;
    std::cout << "System CPU time: " << s_secs << " seconds" << std::endl;
    std::cout << "Max RSS:         " << ru.ru_maxrss/1024 << " MB" << std::endl;
}