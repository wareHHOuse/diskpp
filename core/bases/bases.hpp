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

#include <type_traits>

#define _BASES_HPP_WAS_INCLUDED_

namespace disk {

template<typename MeshType, typename Element>
struct scaled_monomial_scalar_basis
{
    static_assert(sizeof(MeshType) == -1, "scaled_monomial_scalar_basis: not suitable for the requested kind of mesh");
    static_assert(sizeof(Element) == -1, "scaled_monomial_scalar_basis: not suitable for the requested kind of element");
};

template<typename MeshType, typename Element>
struct scaled_monomial_vector_sg_basis
{
    static_assert(sizeof(MeshType) == -1, "scaled_monomial_vector_sg_basis: not suitable for the requested kind of mesh");
    static_assert(sizeof(Element) == -1, "scaled_monomial_vector_sg_basis: not suitable for the requested kind of element");
};

} //namespace disk

//#define POWER_CACHE

#include "geometry/geometry.hpp"
#include "bases/bases_ranges.hpp"
#include "bases/bases_bones.hpp"
#include "bases/bases_utils.hpp"

#include "bases/bases_all.hpp"
#include "bases/bases_simplicial.hpp"


namespace disk {

namespace priv {

template<typename MeshType>
using mct = typename MeshType::cell;

template<typename MeshType>
using mft = typename MeshType::face;

template<typename MeshType, typename Element>
using smsb = scaled_monomial_scalar_basis<MeshType, Element>;

} // namespace priv

template<typename MeshType>
std::pair<  priv::smsb<MeshType, priv::mct<MeshType>>,
            priv::smsb<MeshType, priv::mft<MeshType>>  >
make_scaled_monomial_scalar_basis(const MeshType& msh, size_t degree)
{
    auto cb = priv::smsb<MeshType, priv::mct<MeshType>>(degree);
    auto fb = priv::smsb<MeshType, priv::mft<MeshType>>(degree);

    return std::make_pair(cb, fb);
}

template<typename MeshType>
std::pair<  priv::smsb<MeshType, priv::mct<MeshType>>,
            priv::smsb<MeshType, priv::mft<MeshType>>  >
make_scaled_monomial_scalar_basis(const MeshType& msh, size_t cell_degree,
                                  size_t face_degree)
{
    auto cb = priv::smsb<MeshType, priv::mct<MeshType>>(cell_degree);
    auto fb = priv::smsb<MeshType, priv::mft<MeshType>>(face_degree);

    return std::make_pair(cb, fb);
}

} // namespace disk



#undef _BASES_HPP_WAS_INCLUDED_
