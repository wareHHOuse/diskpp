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

namespace disk {

template<typename MeshType, typename Element>
class scaled_monomial_scalar_basis
{
    static_assert(sizeof(MeshType) == -1, "scaled_monomial_scalar_basis: not suitable for the requested kind of mesh");
    static_assert(sizeof(Element) == -1, "scaled_monomial_scalar_basis: not suitable for the requested kind of element");
};

template<typename MeshType, typename Element>
class scaled_monomial_vector_sg_basis
{
    static_assert(sizeof(MeshType) == -1, "scaled_monomial_vector_sg_basis: not suitable for the requested kind of mesh");
    static_assert(sizeof(Element) == -1, "scaled_monomial_vector_sg_basis: not suitable for the requested kind of element");
};

} //namespace disk
