/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016-2017 - matteo.cicuttin@enpc.fr
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

/*
 * Copyright (C) 2013-2016, Matteo Cicuttin - matteo.cicuttin@uniud.it
 * Department of Electrical Engineering, University of Udine
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Udine nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR(s) ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(s) BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#ifndef _LOADER_HPP_WAS_INCLUDED_
#error "You must NOT include this file directly. Include loader.hpp"
#endif

#ifndef _LOADER_UNIFORM_HPP_
#define _LOADER_UNIFORM_HPP_


#include <vector>
#include <array>
#include <fstream>
#include <cassert>
#include <set>

#include "geometry/geometry.hpp"
#include "loader_utils.hpp"
#include "loader_template.hpp"

#include "mapped_file.h"
#include "strtot.hpp"

namespace disk {

   template<typename T, size_t N>
   class uniform_mesh_loader
   {
      static_assert(N == 1, "at the moment only 1D uniform meshes are available");
   };

   template<typename T>
   class uniform_mesh_loader<T,1> : public mesh_loader<generic_mesh<T, 1>>
   {
      typedef generic_mesh<T,1>                       mesh_type;
      typedef typename mesh_type::point_type          point_type;
      typedef typename mesh_type::node_type           node_type;
      typedef typename mesh_type::edge_type           edge_type;

      T       m_h, m_x_min, m_x_max;
      size_t  m_number_of_elements;

   public:
      uniform_mesh_loader()
      : m_x_min(T(0)), m_x_max(T(1)), m_number_of_elements(8)
      {
         m_h = fabs(m_x_max - m_x_min) / m_number_of_elements;
      }

      uniform_mesh_loader(T x_min, T x_max, size_t N)
      : m_x_min(x_min), m_x_max(x_max), m_number_of_elements(N)
      {
         m_h = fabs(m_x_max - m_x_min) / m_number_of_elements;
      }

      bool populate_mesh(mesh_type& msh)
      {
         if (this->verbose())
            std::cout << " *** POPULATING UNIFORM 1D MESH ***" << std::endl;

         auto storage = msh.backend_storage();

         auto num_edges = m_number_of_elements;
         auto num_nodes = m_number_of_elements + 1;

         storage->points.resize(num_nodes);
         storage->nodes.resize(num_nodes);
         storage->edges.resize(num_edges);

         for (size_t i = 0; i < num_edges; i++)
         {
            storage->points.at(i)   = point_type({m_x_min + (m_h * i)});
            storage->points.at(i+1) = point_type({m_x_min + (m_h * (i + 1))});

            auto n0 = typename node_type::id_type(i);
            auto n1 = typename node_type::id_type(i+1);
            auto e = edge_type{{n0, n1}};

            std::vector<point_identifier<1>> pts(2);
            pts[0] = point_identifier<1>(i);
            pts[1] = point_identifier<1>(i+1);
            e.set_point_ids(pts.begin(), pts.end());

            storage->edges.at(i) = e;
         }

         for (size_t i = 0; i < num_nodes; i++)
            storage->nodes.at(i) = node_type(point_identifier<1>(i));

         storage->boundary_info.resize(num_nodes);
         bnd_info bi{0, true};
         storage->boundary_info.at(0) = bi;
         storage->boundary_info.at(num_nodes - 1) = bi;

         return true;
      }
   };

} // namespace disk

#endif
