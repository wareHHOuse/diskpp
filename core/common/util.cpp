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

#include <ostream>
#include <iostream>
#include "util.h"

unsigned int fact(unsigned int n)
{
   return (n < 2) ? 1 : n*fact(n-1);
}

unsigned int binomial(unsigned int n, unsigned int k)
{
   if(k == 0)
      return 1;

   else if (n < k)
      return fact(n) / fact(k);

   if ( k > (n-k) )
   {
      unsigned int num(1);
      unsigned int dem(1);

      for(unsigned int i=1; i <= (n-k); i++)
      {
         num = num * (k+i);
         dem = dem * i;
      }

      return num/dem;
   }
   else
   {
      unsigned int num(1);
      unsigned int dem(1);

      for(unsigned int i=1; i <= k; i++)
      {
         num = num * (n-k+i);
         dem = dem * i;
      }

      return num/dem;
   }
}
