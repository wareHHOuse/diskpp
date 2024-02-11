/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */
/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2020, 2021
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group
 */
/*
 *       /\        Matteo Cicuttin (C) 2016-2019
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

#include <cassert>
#include <fstream>
#include <regex>

#include "diskpp/mesh/mesh.hpp"

#include "loader_cartesian.hpp"
#include "loader_fvca5.hpp"
#include "loader_fvca6.hpp"
#include "loader_medit.hpp"
#include "loader_netgen.hpp"
#include "loader_poly.hpp"
#include "loader_uniform.hpp"
#include "loader_utils.hpp"

namespace disk
{

/**************************************************************************/
// New mesh loader helpers
namespace priv
{

bool
check_filename_extension(const char* filename, const char* extension)
{
    std::stringstream ss;
    ss << ".*\\." << extension << "$";

    return std::regex_match(filename, std::regex(ss.str()));
}

template<typename LT, typename MT>
bool
load_mesh(const char* filename, LT& loader, MT& msh)
{
    if (!check_filename_extension(filename, LT::expected_extension))
    {
        std::cout << "Warning: unexpected filename extension for ";
        std::cout << "the required mesh type" << std::endl;
    }

    bool success = loader.read_mesh(filename);
    if (!success)
        return false;

    loader.populate_mesh(msh);
    return true;
}

} // namespace priv

template<typename T>
bool
load_mesh_netgen(const char* filename, disk::simplicial_mesh<T, 2>& msh)
{
    disk::netgen_mesh_loader<T, 2> loader;
    return priv::load_mesh(filename, loader, msh);
}

template<typename T>
bool
load_mesh_netgen(const char* filename, disk::simplicial_mesh<T, 3>& msh)
{
    disk::netgen_mesh_loader<T, 3> loader;
    return priv::load_mesh(filename, loader, msh);
}

template<typename T>
bool
load_mesh_diskpp_cartesian(const char* filename, disk::cartesian_mesh<T, 2>& msh)
{
    disk::cartesian_mesh_loader<T, 2> loader;
    return priv::load_mesh(filename, loader, msh);
}

template<typename T>
bool
load_mesh_diskpp_cartesian(const char* filename, disk::cartesian_mesh<T, 3>& msh)
{
    disk::cartesian_mesh_loader<T, 3> loader;
    return priv::load_mesh(filename, loader, msh);
}

template<typename T>
bool
load_mesh_fvca5_2d(const char* filename, disk::generic_mesh<T, 2>& msh)
{
    disk::fvca5_mesh_loader<T, 2> loader;
    return priv::load_mesh(filename, loader, msh);
}

template<typename T>
bool
load_mesh_fvca6_3d(const char* filename, disk::generic_mesh<T, 3>& msh)
{
    disk::fvca6_mesh_loader<T, 3> loader;
    return priv::load_mesh(filename, loader, msh);
}

template<typename T>
bool
load_mesh_poly(const char* filename, disk::generic_mesh<T, 2>& msh, bool verbose = false)
{
    disk::poly_mesh_loader<T, 2> loader;
    loader.verbose(verbose);
    return priv::load_mesh(filename, loader, msh);
}

template<typename T>
bool
load_mesh_poly(const char* filename, disk::generic_mesh<T, 3>& msh, bool verbose = false)
{
    disk::poly_mesh_loader<T, 3> loader;
    loader.verbose(verbose);
    return priv::load_mesh(filename, loader, msh);
}

template<typename T>
bool
load_mesh_medit(const char* filename, disk::generic_mesh<T, 2>& msh)
{
    disk::medit_mesh_loader<T, 2> loader;
    return priv::load_mesh(filename, loader, msh);
}

template<typename T>
bool
load_mesh_medit(const char* filename, disk::generic_mesh<T, 3>& msh)
{
    disk::medit_mesh_loader<T, 3> loader;
    return priv::load_mesh(filename, loader, msh);
}

} // namespace disk

#ifdef HAVE_GMSH
#include "loader_gmsh.hpp"
#endif

namespace disk
{

/**
 * @brief Deduce mesh type from filename extension, create appropriate mesh object
 *        and dispatch a function on it
 *
 * If you have the function `process_mesh(msh, p1, p2, ..., pn)` you can
 *
 *  ```
 *  disk::dispatch_all_meshes(mesh_filename,
 *        [](auto ...args) { process_mesh(args...); },
 *        p1, p2, ..., pn);
 *  ```
 *
 * @tparam Function
 * @tparam Args
 * @param mesh_filename Name of the file containing the mesh.
 * @param func Function to dispatch. The mesh object is passed as first parameter.
 * @param args Arguments to the function to dispatch.
 * @return int
 */

template<typename Function, typename... Args>
int
dispatch_all_meshes(const char* mesh_filename, Function func, Args&&... args)
{
    using T = double;
    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$")))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        disk::generic_mesh<T, 2> msh;
        if (load_mesh_fvca5_2d(mesh_filename, msh))
        {
            func(msh, args...);
            return 0;
        }
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$")))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        disk::simplicial_mesh<T, 2> msh;
        if (load_mesh_netgen(mesh_filename, msh))
        {
            func(msh, args...);
            return 0;
        }
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$")))
    {
        disk::cartesian_mesh<T, 2> msh;
        if (load_mesh_diskpp_cartesian(mesh_filename, msh))
        {
            func(msh, args...);
            return 0;
        }
    }

    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$")))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        disk::simplicial_mesh<T, 3> msh;
        if (load_mesh_netgen(mesh_filename, msh))
        {
            func(msh, args...);
            return 0;
        }
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$")))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        disk::cartesian_mesh<T, 3> msh;
        if (load_mesh_diskpp_cartesian(mesh_filename, msh))
        {
            func(msh, args...);
            return 0;
        }
    }

    /* FVCA6 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.msh$")))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
        disk::generic_mesh<T, 3> msh;
        if (load_mesh_fvca6_3d(mesh_filename, msh))
        {
            func(msh, args...);
            return 0;
        }
    }

#ifdef HAVE_GMSH
    /* GMSH 2D simplicials */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo2s$")))
    {
        std::cout << "Guessed mesh format: GMSH 2D simplicials" << std::endl;
        disk::simplicial_mesh<T, 2>                             msh;
        disk::gmsh_geometry_loader<disk::simplicial_mesh<T, 2>> loader;
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);
        func(msh, args...);
        return 0;
    }

    /* GMSH 3D simplicials */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3s$")))
    {
        std::cout << "Guessed mesh format: GMSH 3D simplicials" << std::endl;
        disk::simplicial_mesh<T, 3>                             msh;
        disk::gmsh_geometry_loader<disk::simplicial_mesh<T, 3>> loader;
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);
        func(msh, args...);
        return 0;
    }

    /* GMSH 3D generic */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3g$")))
    {
        std::cout << "Guessed mesh format: GMSH 3D general elements" << std::endl;
        disk::generic_mesh<T, 3>                             msh;
        disk::gmsh_geometry_loader<disk::generic_mesh<T, 3>> loader;
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);
        func(msh, args...);
        return 0;
    }
#endif

    return 1;
}

} // namespace disk
