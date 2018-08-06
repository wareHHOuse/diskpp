/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
 * Nicolas Pignet  (C) 2018                     nicolas.pignet@enpc.fr
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

#define _QUADRATURES_HPP_WAS_INCLUDED_

#include "geometry/geometry.hpp"
#include "quadratures/quad_bones.hpp"

namespace disk {

template<typename MeshType, typename Element>
std::vector<disk::quadrature_point<typename MeshType::scalar_type, MeshType::dimension>>
integrate(const MeshType& msh, const Element& elem, const size_t degree)
{
    static_assert(sizeof(MeshType) == -1, "quadrature: not suitable for the requested kind of mesh");
    static_assert(sizeof(Element) == -1, "quadrature: not suitable for the requested kind of element");
};

} // namespace disk

#include "quadratures/quad_simplicial.hpp"
#include "quadratures/quad_hexahedral.hpp"
#include "quadratures/quad_generic.hpp"

#undef _QUADRATURES_HPP_WAS_INCLUDED_
