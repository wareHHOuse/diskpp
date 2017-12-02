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

#define _QUADRATURES_HPP_WAS_INCLUDED_

#include "geometry/geometry.hpp"
#include "quadratures/quad_bones.hpp"

namespace disk
{
template <typename MeshType, typename Element>
class quadrature
{
    typedef MeshType                            mesh_type;
    typedef Element                             element_type;
    typedef typename mesh_type::point_type      point_type;
    typedef typename mesh_type::scalar_type     scalar_type;

    static_assert(sizeof(MeshType) == -1, "quadrature: not suitable for the requested kind of mesh");
    static_assert(sizeof(Element) == -1, "quadrature: not suitable for the requested kind of element");
};

template <typename Mesh>
std::pair<quadrature<Mesh, typename Mesh::cell>, quadrature<Mesh, typename Mesh::face>>
make_quadrature(const Mesh &msh, size_t order_cell, size_t order_face)
{
    const auto cq = quadrature<Mesh, typename Mesh::cell>(order_cell);
    const auto fq = quadrature<Mesh, typename Mesh::face>(order_face);

    return std::make_pair(cq, fq);
}

} // namespace disk

#include "quadratures/quad_simplicial.hpp"
#include "quadratures/quad_hexahedral.hpp"
#include "quadratures/quad_generic.hpp"

#undef _QUADRATURES_HPP_WAS_INCLUDED_
