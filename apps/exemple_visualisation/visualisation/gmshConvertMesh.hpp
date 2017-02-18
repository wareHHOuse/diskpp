/*
 *       /\
 *      /__\       Nicolas Pignet (C) 2017
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    Gmsh tools
 *  /_\/_\/_\/_\
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#ifndef _GmshConvertMesh_HPP
#define _GmshConvertMesh_HPP


#include<vector>
#include<string>
#include<utility>
#include <array>
#include <cassert>

#include "loaders/loader.hpp"
#include "geometry/geometry.hpp"
#include "mesh/point.hpp"
#include "gmshMesh.h"
#include "gmshElement.h"
#include "gmshData.h"
#include "gmshDisk.h"

namespace visu{
template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage, size_t DIM>
Gmesh convertMesh(const Mesh<T, DIM, Storage>& mesh)
{
    typedef Mesh<T, DIM, Storage>                 mesh_type;

    static_assert(sizeof(mesh_type) == -1, "convert: not suitable for the requested kind of mesh");
}

} //visu

#include "gmshConvertMeshSimplicial.hpp"
#include "gmshConvertMeshHexahedral.hpp"
#include "gmshConvertMeshGeneric.hpp"

#endif
