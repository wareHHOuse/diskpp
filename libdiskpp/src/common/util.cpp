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
#include <sys/select.h>

#include "diskpp/common/util.h"

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
    : rm_enabled(true)
{
    gettimeofday(&time_start, nullptr);
}

rusage_monitor::rusage_monitor(bool enable)
    : rm_enabled(enable)
{
    gettimeofday(&time_start, nullptr);
}

bool
rusage_monitor::enabled(void) const
{
    return rm_enabled;
}

void
rusage_monitor::enabled(bool enable)
{
    rm_enabled = enable;
}

rusage_monitor::~rusage_monitor()
{
    if (not enabled())
        return;

    gettimeofday(&time_end, nullptr);

    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);

    struct timeval diff;
    timersub(&time_end, &time_start, &diff);

    double u_secs = ru.ru_utime.tv_sec + ru.ru_utime.tv_usec/1e6;
    double s_secs = ru.ru_stime.tv_sec + ru.ru_stime.tv_usec/1e6;

    std::cout << "Resource usage report:" << std::endl;
    std::cout << "  Wall time:       " << diff.tv_sec << " seconds" << std::endl;
    std::cout << "  User CPU time:   " << u_secs << " seconds" << std::endl;
    std::cout << "  System CPU time: " << s_secs << " seconds" << std::endl;
    std::cout << "  Max RSS:         " << ru.ru_maxrss/1024 << " MB" << std::endl;
}