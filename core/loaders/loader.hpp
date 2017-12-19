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

#define _LOADER_HPP_WAS_INCLUDED_

#include "loader_cartesian.hpp"
#include "loader_fvca.hpp"
#include "loader_medit.hpp"
#include "loader_netgen.hpp"
#include "loader_uniform.hpp"

#include "geometry/geometry.hpp"

namespace disk {

/* Helper to load 2D meshes in DiSk++ format */
template<typename T>
disk::cartesian_mesh<T, 2>
load_cartesian_2d_mesh(const char* filename)
{
   typedef disk::cartesian_mesh<T, 2> mesh_type;

   mesh_type                         msh;
   disk::cartesian_mesh_loader<T, 2> loader;
   loader.read_mesh(filename);
   loader.populate_mesh(msh);

   return msh;
}

/* Helper to load 2D meshes in DiSk++ format */
template<typename T>
disk::cartesian_mesh<T, 2>
load_cartesian_2d_mesh2(const char* filename)
{
   typedef disk::cartesian_mesh<T, 2> mesh_type;

   mesh_type                          msh;
   disk::cartesian_mesh_loader2<T, 2> loader;
   loader.read_mesh(filename);
   loader.populate_mesh(msh);

   return msh;
}

/* Helper to load 3D meshes in DiSk++ format */
template<typename T>
disk::cartesian_mesh<T, 3>
load_cartesian_3d_mesh(const char* filename)
{
   typedef disk::cartesian_mesh<T, 3> mesh_type;

   mesh_type                         msh;
   disk::cartesian_mesh_loader<T, 3> loader;
   loader.read_mesh(filename);
   loader.populate_mesh(msh);

   return msh;
}

/* Helper to load 2D meshes in FVCA5 format */
template<typename T>
disk::generic_mesh<T, 2>
load_fvca5_2d_mesh(const char* filename)
{
   typedef disk::generic_mesh<T, 2> mesh_type;

   mesh_type                     msh;
   disk::fvca5_mesh_loader<T, 2> loader;
   loader.read_mesh(filename);
   loader.populate_mesh(msh);

   return msh;
}

/* Helper to load 3D meshes in FVCA6 format */
template<typename T>
disk::generic_mesh<T, 3>
load_fvca6_3d_mesh(const char* filename)
{
   typedef disk::generic_mesh<T, 3> mesh_type;

   mesh_type                     msh;
   disk::fvca6_mesh_loader<T, 3> loader;
   loader.read_mesh(filename);
   loader.populate_mesh(msh);

   return msh;
}

/* Helper to load 2D meshes in Medit format */
template<typename T>
disk::generic_mesh<T, 2>
load_medit_2d_mesh(const char* filename)
{
   typedef disk::generic_mesh<T, 2> mesh_type;

   mesh_type                     msh;
   disk::medit_mesh_loader<T, 2> loader;
   loader.read_mesh(filename);
   loader.populate_mesh(msh);

   return msh;
}

/* Helper to load 3D meshes in Medit format */
template<typename T>
disk::generic_mesh<T, 3>
load_medit_3d_mesh(const char* filename)
{
   typedef disk::generic_mesh<T, 3> mesh_type;

   mesh_type                     msh;
   disk::medit_mesh_loader<T, 3> loader;
   loader.read_mesh(filename);
   loader.populate_mesh(msh);

   return msh;
}

/* Helper to load 2D meshes in Netgen format */
template<typename T>
disk::simplicial_mesh<T, 2>
load_netgen_2d_mesh(const char* filename)
{
   typedef disk::simplicial_mesh<T, 2> mesh_type;

   mesh_type                      msh;
   disk::netgen_mesh_loader<T, 2> loader;
   loader.read_mesh(filename);
   loader.populate_mesh(msh);

   return msh;
}

/* Helper to load 3D meshes in Netgen format */
template<typename T>
disk::simplicial_mesh<T, 3>
load_netgen_3d_mesh(const char* filename)
{
   typedef disk::simplicial_mesh<T, 3> mesh_type;

   mesh_type                      msh;
   disk::netgen_mesh_loader<T, 3> loader;
   loader.read_mesh(filename);
   loader.populate_mesh(msh);

   return msh;
}

/* Helper to load uniform 1D meshes. */
template<typename T>
disk::generic_mesh<T, 1>
load_uniform_1d_mesh(T min, T max, size_t cells)
{
   typedef disk::generic_mesh<T, 1> mesh_type;

   mesh_type                       msh;
   disk::uniform_mesh_loader<T, 1> loader(min, max, cells);
   loader.populate_mesh(msh);

   return msh;
}

} // namespace disk

#undef _LOADER_HPP_WAS_INCLUDED_
