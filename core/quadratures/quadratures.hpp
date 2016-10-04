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

#define _QUADRATURES_HPP_WAS_INCLUDED_

#include "geometry/geometry.hpp"
#include "quadratures/quad_bones.hpp"

namespace disk {
template<typename MeshType, typename Element>
class quadrature
{
    typedef MeshType                            mesh_type;
    typedef Element                             element_type;
    typedef typename mesh_type::point_type      point_type;
    typedef typename mesh_type::scalar_type     scalar_type;

    static_assert(sizeof(MeshType) == -1, "quadrature: not suitable for the requested kind of mesh");
    static_assert(sizeof(Element) == -1, "quadrature: not suitable for the requested kind of element");
};
} // namespace disk

#include "quadratures/quad_simplicial.hpp"
#include "quadratures/quad_generic.hpp"

#undef _QUADRATURES_HPP_WAS_INCLUDED_
