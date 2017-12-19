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

#pragma once

#include <thread>
#include <vector>

#include "loaders/strtot.hpp"

namespace disk {

namespace priv {

template<typename T>
point<T, 2>
read_2d_point_line(const char* str, char** endptr, T scalefactor)
{
   T t1, t2;

   t1 = strtot<T>(str, endptr);
   t2 = strtot<T>(*endptr, endptr);

   return point<T, 2>{t1 * scalefactor, t2 * scalefactor};
}

template<typename T>
point<T, 3>
read_3d_point_line(const char* str, char** endptr, T scalefactor)
{
   T t1, t2, t3;

   t1 = strtot<T>(str, endptr);
   t2 = strtot<T>(*endptr, endptr);
   t3 = strtot<T>(*endptr, endptr);

   return point<T, 3>{t1 * scalefactor, t2 * scalefactor, t3 * scalefactor};
}

template<typename T>
void
sort_uniq(std::vector<T>& v)
{
   std::sort(v.begin(), v.end());
   auto uniq_iter = std::unique(v.begin(), v.end());
   v.erase(uniq_iter, v.end());
}

} // namespace priv

} // namespace disk

#define THREADED
#ifdef THREADED
#define THREAD(name, body) std::thread name([&] { body })
#define WAIT_THREAD(name) name.join()
#else
#define THREAD(name, body)                                                                         \
   {                                                                                               \
      body                                                                                         \
   }
#define WAIT_THREAD(name)
#endif
