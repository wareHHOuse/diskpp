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

#pragma once

#define _GEOMETRY_HPP_WAS_INCLUDED_

namespace disk {

template<typename T, size_t DIM>
struct storage_class_trait
{
    static_assert(sizeof(T) == -1, "Undefined storage class");
};

}

#include "mesh/mesh.hpp"

#include "geometry_all.hpp"
#include "geometry_generic.hpp"
#include "geometry_simplicial.hpp"
#include "geometry_hexahedral.hpp"




#undef _GEOMETRY_HPP_WAS_INCLUDED_
