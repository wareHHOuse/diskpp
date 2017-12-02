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

#ifndef _LOADER_NETGEN_HPP_
#define _LOADER_NETGEN_HPP_


#include <vector>
#include <array>
#include <fstream>
#include <cassert>
#include <thread>
#include <set>

#include "geometry/geometry.hpp"
#include "loader_utils.hpp"
#include "loader_template.hpp"

#include "mapped_file.h"
#include "strtot.hpp"

namespace disk {

   template<typename T, size_t N>
   class netgen_mesh_loader
   {
      static_assert(N == 2 || N == 3, "Netgen supports only 2D or 3D");
   };

   namespace priv {

      template<typename T>
      point<T,2>
      read_2d_point_line(const char *str, char **endptr, T scalefactor)
      {
         T t1, t2;

         t1 = strtot<T>(str, endptr);
         t2 = strtot<T>(*endptr, endptr);

         return point<T,2>{t1*scalefactor, t2*scalefactor};
      }

      template<typename T>
      point<T,3>
      read_3d_point_line(const char *str, char **endptr, T scalefactor)
      {
         T t1, t2, t3;

         t1 = strtot<T>(str, endptr);
         t2 = strtot<T>(*endptr, endptr);
         t3 = strtot<T>(*endptr, endptr);

         return point<T,3>{t1*scalefactor, t2*scalefactor, t3*scalefactor};
      }

      template<typename T>
      std::tuple<T, T, T, T, T>
      read_tetrahedron_line(const char *str, char **endptr)
      {
         T t1, t2, t3, t4, t5;

         t1 = strtot<T>(str, endptr);
         t2 = strtot<T>(*endptr, endptr);
         t3 = strtot<T>(*endptr, endptr);
         t4 = strtot<T>(*endptr, endptr);
         t5 = strtot<T>(*endptr, endptr);

         return std::make_tuple(t1, t2-1, t3-1, t4-1, t5-1);
      }

      template<typename T>
      std::tuple<T, T, T, T>
      read_triangle_line(const char *str, char **endptr)
      {
         T t1, t2, t3, t4;

         t1 = strtot<T>(str, endptr);
         t2 = strtot<T>(*endptr, endptr);
         t3 = strtot<T>(*endptr, endptr);
         t4 = strtot<T>(*endptr, endptr);

         return std::make_tuple(t1, t2-1, t3-1, t4-1);
      }

      template<typename T>
      std::tuple<T, T, T>
      read_edge_line(const char *str, char **endptr)
      {
         T t1, t2, t3;

         t1 = strtot<T>(str, endptr);
         t2 = strtot<T>(*endptr, endptr);
         t3 = strtot<T>(*endptr, endptr);

         return std::make_tuple(t1, t2-1, t3-1);
      }

   } //namespace priv

   template<typename T>
   class netgen_mesh_loader<T,2> : public mesh_loader<simplicial_mesh<T,2>>
   {
      typedef simplicial_mesh<T,2>                    mesh_type;
      typedef typename mesh_type::point_type          point_type;
      typedef typename mesh_type::node_type           node_type;
      typedef typename mesh_type::edge_type           edge_type;
      typedef typename mesh_type::surface_type        surface_type;

      std::vector<point_type>                         points;
      std::vector<node_type>                          nodes;
      std::vector<edge_type>                          edges;
      std::vector<std::pair<edge_type, size_t> >      boundary_edges;
      std::vector<surface_type>                       surfaces;

      bool netgen_read(const std::string& filename)
      {
         /* Open file */
         if (filename.size() == 0)
         {
            std::cout << "Invalid mesh file name" << std::endl;
            return false;
         }

         size_t lines, linecount;

         mapped_file mf(filename);

         //std::cout << green << " * * * Reading NETGEN format mesh * * * ";
         //std::cout << nocolor << std::endl;

         /************************ Read points ************************/
         linecount = 0;

         const char *data = mf.mem();
         char *endptr;

         lines = strtot<size_t>(data, &endptr);

         points.reserve(lines);
         nodes.reserve(lines);

         while (linecount < lines)
         {
            if ( this->verbose() && ((linecount%100000) == 0) )
            {
               std::cout << "Reading points: " << linecount;
               std::cout << "/" << lines << "\r";
               std::cout.flush();
            }

            auto point = priv::read_2d_point_line<T>(endptr, &endptr, 1.0);

            points.push_back( point );

            auto point_id = point_identifier<2>( linecount );
            auto node = node_type( { point_id } );

            nodes.push_back(node);

            linecount++;
         }

         if (this->verbose())
         {
            std::cout << "Reading points: " << linecount;
            std::cout << "/" << lines << std::endl;
         }

         /************************ Read triangles ************************/
         linecount = 0;

         lines = strtot<size_t>(endptr, &endptr);

         edges.reserve(lines*3);
         surfaces.reserve(lines);

         while (linecount < lines)
         {
            if ( this->verbose() && ((linecount%100000) == 0) )
            {
               std::cout << "Reading triangles: " << linecount;
               std::cout << "/" << lines << "\r";
               std::cout.flush();
            }

            auto t = priv::read_triangle_line<size_t>(endptr, &endptr);

            point_identifier<2>     p0(std::get<1>(t));
            point_identifier<2>     p1(std::get<2>(t));
            point_identifier<2>     p2(std::get<3>(t));
            //domain_id_type      d(std::get<0>(t));

            edges.push_back( edge_type( { p0, p1 } ) );
            edges.push_back( edge_type( { p1, p2 } ) );
            edges.push_back( edge_type( { p0, p2 } ) );

            surfaces.push_back( surface_type( { p0, p1, p2 } ) );

            linecount++;
         }

         if (this->verbose())
         {
            std::cout << "Reading triangles: " << linecount;
            std::cout << "/" << lines << std::endl;
         }
         /************************ Read boundary edges ************************/
         linecount = 0;

         lines = strtot<size_t>(endptr, &endptr);

         boundary_edges.reserve(lines);

         while (linecount < lines)
         {
            if ( this->verbose() && ((linecount%50000) == 0) )
            {
               std::cout << "Reading edges: " << linecount;
               std::cout << "/" << lines << "\r";
               std::cout.flush();
            }

            auto t = priv::read_edge_line<size_t>(endptr, &endptr);

            point_identifier<2>     p0(std::get<1>(t));
            point_identifier<2>     p1(std::get<2>(t));

            edge_type   edge( { p0, p1 } );

            boundary_edges.push_back( std::make_pair(edge, std::get<0>(t)));

            linecount++;
         }

         if (this->verbose())
         {
            std::cout << "Reading edges: " << linecount;
            std::cout << "/" << lines << std::endl;
         }

         return true;
      }

   public:
      netgen_mesh_loader() = default;

      bool read_mesh(const std::string& s)
      {
         if (this->verbose())
            std::cout << " *** READING NETGEN 2D MESH ***" << std::endl;

         return netgen_read(s);
      }

      bool populate_mesh(mesh_type& msh)
      {
         auto storage = msh.backend_storage();

         if (this->verbose())
         {
            std::cout << "Sorting data...";
            std::cout.flush();
         }

         storage->points = std::move(points);
         storage->nodes = std::move(nodes);

         /* sort edges, make unique and move them in geometry */
         THREAD(edge_thread,
                priv::sort_uniq(edges);
                storage->edges = std::move(edges);
         );

         /* sort triangles, make unique and move them in geometry */
         THREAD(tri_thread,
                priv::sort_uniq(surfaces);
                storage->surfaces = std::move(surfaces);
         );

         /* wait for the threads */
         WAIT_THREAD(edge_thread);
         WAIT_THREAD(tri_thread);

         storage->boundary_info.resize(storage->edges.size());
         for (auto& be : boundary_edges)
         {
            auto position = find_element_id(storage->edges.begin(),
                                            storage->edges.end(), be.first);
            if (position.first == false)
            {
               std::cout << "Bad bug at " << __FILE__ << "("
               << __LINE__ << ")" << std::endl;
               return false;
            }
            bnd_info bi{be.second, true};
            storage->boundary_info.at(position.second) = bi;
         }

         if (this->verbose())
         {
            std::cout << "done." << std::endl;

            std::cout << "Nodes: " << storage->nodes.size() << std::endl;
            std::cout << "Edges: " << storage->edges.size() << std::endl;
            std::cout << "Faces: " << storage->surfaces.size() << std::endl;
         }

         boundary_edges.clear();

         return true;
      }
   };

   template<typename T>
   class netgen_mesh_loader<T,3> : public mesh_loader<simplicial_mesh<T,3>>
   {
      typedef simplicial_mesh<T,3>                    mesh_type;
      typedef typename mesh_type::point_type          point_type;
      typedef typename mesh_type::node_type           node_type;
      typedef typename mesh_type::edge_type           edge_type;
      typedef typename mesh_type::surface_type        surface_type;
      typedef typename mesh_type::volume_type         volume_type;

      std::vector<point_type>                         points;
      std::vector<node_type>                          nodes;
      std::vector<edge_type>                          edges;
      std::vector<surface_type>                       surfaces;
      std::vector<std::pair<surface_type, size_t> >   boundary_surfaces;
      std::vector<volume_type>                        volumes;


      bool netgen_read(const std::string& filename)
      {
         /* Open file */
         if (filename.size() == 0)
         {
            std::cout << "Invalid mesh file name" << std::endl;
            return false;
         }

         size_t lines, linecount;

         mapped_file mf(filename);

         //std::cout << green << " * * * Reading NETGEN format mesh * * * ";
         //std::cout << nocolor << std::endl;

         /************************ Read points ************************/
         linecount = 0;

         const char *data = mf.mem();
         char *endptr;

         lines = strtot<size_t>(data, &endptr);

         points.reserve(lines);
         nodes.reserve(lines);

         while (linecount < lines)
         {
            if ( this->verbose() && ((linecount%100000) == 0) )
            {
               std::cout << "Reading points: " << linecount;
               std::cout << "/" << lines << "\r";
               std::cout.flush();
            }

            auto point = priv::read_3d_point_line<T>(endptr, &endptr, 1.0);

            /*auto point = point_type( { std::get<0>(t),
             *                                     std::get<1>(t),
             *                                     std::get<2>(t) } );
             */
            points.push_back( point );

            auto point_id = point_identifier<3>( linecount );
            auto node = node_type( { point_id } );

            nodes.push_back(node);
            /* Do something with that point */

            linecount++;
         }

         if (this->verbose())
         {
            std::cout << "Reading points: " << linecount;
            std::cout << "/" << lines << std::endl;
         }

         /************************ Read tetrahedra ************************/
         linecount = 0;

         lines = strtot<size_t>(endptr, &endptr);

         edges.reserve(lines*6);
         surfaces.reserve(lines*4);
         volumes.reserve(lines);

         while (linecount < lines)
         {
            if ( this->verbose() && ((linecount%100000) == 0) )
            {
               std::cout << "Reading tetrahedra: " << linecount;
               std::cout << "/" << lines << "\r";
               std::cout.flush();
            }

            auto t = priv::read_tetrahedron_line<size_t>(endptr, &endptr);

            point_identifier<3>     p0(std::get<1>(t));
            point_identifier<3>     p1(std::get<2>(t));
            point_identifier<3>     p2(std::get<3>(t));
            point_identifier<3>     p3(std::get<4>(t));
            //domain_id_type      d(std::get<0>(t));

            edges.push_back( edge_type( { p0, p1 } ) );
            edges.push_back( edge_type( { p0, p2 } ) );
            edges.push_back( edge_type( { p0, p3 } ) );
            edges.push_back( edge_type( { p1, p2 } ) );
            edges.push_back( edge_type( { p1, p3 } ) );
            edges.push_back( edge_type( { p2, p3 } ) );

            surfaces.push_back( surface_type( { p0, p1, p2 } ) );
            surfaces.push_back( surface_type( { p0, p1, p3 } ) );
            surfaces.push_back( surface_type( { p0, p2, p3 } ) );
            surfaces.push_back( surface_type( { p1, p2, p3 } ) );

            //auto tuple = std::make_tuple(volume_type(p0, p1, p2, p3), d);
            //temp_tet.push_back( tuple );

            volumes.push_back( volume_type( { p0, p1, p2, p3 } ) );

            linecount++;
         }

         if (this->verbose())
         {
            std::cout << "Reading tetrahedra: " << linecount;
            std::cout << "/" << lines << std::endl;
         }

         /************************ Read boundary surfaces ************************/
         linecount = 0;

         lines = strtot<size_t>(endptr, &endptr);

         boundary_surfaces.reserve(lines);

         while (linecount < lines)
         {
            if ( this->verbose() && ((linecount%50000) == 0) )
            {
               std::cout << "Reading triangle: " << linecount;
               std::cout << "/" << lines << "\r";
               std::cout.flush();
            }

            auto t = priv::read_triangle_line<size_t>(endptr, &endptr);

            point_identifier<3>     p0(std::get<1>(t));
            point_identifier<3>     p1(std::get<2>(t));
            point_identifier<3>     p2(std::get<3>(t));

            surface_type   tri( { p0, p1, p2 } );

            boundary_surfaces.push_back(std::make_pair(tri, std::get<0>(t)) );

            linecount++;
         }

         if (this->verbose())
         {
            std::cout << "Reading triangle: " << linecount;
            std::cout << "/" << lines << std::endl;
         }

         return true;
      }

   public:
      netgen_mesh_loader() = default;

      bool read_mesh(const std::string& s)
      {
         if (this->verbose())
            std::cout << " *** READING NETGEN 3D MESH ***" << std::endl;
         return netgen_read(s);
      }

      bool populate_mesh(mesh_type& msh)
      {
         auto storage = msh.backend_storage();

         if (this->verbose())
         {
            std::cout << "Sorting data...";
            std::cout.flush();
         }

         storage->points = std::move(points);
         storage->nodes = std::move(nodes);

         /* sort edges, make unique and move them in geometry */
         THREAD(edge_thread,
                priv::sort_uniq(edges);
                storage->edges = std::move(edges);
         );

         /* sort triangles, make unique and move them in geometry */
         THREAD(tri_thread,
                priv::sort_uniq(surfaces);
                storage->surfaces = std::move(surfaces);
         );

         /* sort tetrahedra, make unique and move them in geometry */
         THREAD(tet_thread,
                std::sort(volumes.begin(), volumes.end());
                storage->volumes = std::move(volumes);
         );

         /* wait for the threads */
         WAIT_THREAD(edge_thread);
         WAIT_THREAD(tri_thread);
         WAIT_THREAD(tet_thread);

         storage->boundary_info.resize(storage->surfaces.size());
         for (auto& bs : boundary_surfaces)
         {
            auto position = find_element_id(storage->surfaces.begin(),
                                            storage->surfaces.end(), bs.first);
            if (position.first == false)
            {
               std::cout << "Bad bug at " << __FILE__ << "("
               << __LINE__ << ")" << std::endl;
               return false;
            }

            bnd_info bi{bs.second, true};
            storage->boundary_info.at(position.second) = bi;
         }

         if (this->verbose())
         {
            std::cout << "done." << std::endl;

            std::cout << "Nodes: " << storage->nodes.size() << std::endl;
            std::cout << "Edges: " << storage->edges.size() << std::endl;
            std::cout << "Faces: " << storage->surfaces.size() << std::endl;
            std::cout << "Volumes: " << storage->volumes.size() << std::endl;
         }

         boundary_surfaces.clear();

         return true;
      }

   };


} // namespace disk

#endif
