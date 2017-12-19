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

#include "bases/bases.hpp"

namespace disk {

// is_scalar basis ?
template<typename Basis>
struct is_scalar_basis
{
   static const bool value = false;
};

template<typename MeshType>
struct is_scalar_basis<scaled_monomial_scalar_basis<MeshType, typename MeshType::cell>>
{
   static const bool value = true;
};

template<typename MeshType>
struct is_scalar_basis<scaled_monomial_scalar_basis<MeshType, typename MeshType::face>>
{
   static const bool value = true;
};

// is_vector basis ?
template<typename Basis>
struct is_vector_basis
{
   static const bool value = false;
};

template<typename MeshType>
struct is_vector_basis<scaled_monomial_vector_basis<MeshType, typename MeshType::cell>>
{
   static const bool value = true;
};

template<typename MeshType>
struct is_vector_basis<scaled_monomial_vector_basis<MeshType, typename MeshType::face>>
{
   static const bool value = true;
};

template<typename MeshType>
struct is_vector_basis<Raviart_Thomas_vector_basis<MeshType, typename MeshType::cell>>
{
   static const bool value = true;
};

template<typename MeshType>
struct is_vector_basis<Raviart_Thomas_vector_basis<MeshType, typename MeshType::face>>
{
   static const bool value = true;
};

// is_matrix basis ?
template<typename Basis>
struct is_matrix_basis
{
   static const bool value = false;
};

template<typename MeshType>
struct is_matrix_basis<scaled_monomial_matrix_basis<MeshType, typename MeshType::cell>>
{
   static const bool value = true;
};

template<typename MeshType>
struct is_matrix_basis<scaled_monomial_matrix_basis<MeshType, typename MeshType::face>>
{
   static const bool value = true;
};

template<typename MeshType>
struct is_matrix_basis<scaled_monomial_sym_matrix_basis<MeshType, typename MeshType::cell>>
{
   static const bool value = true;
};

template<typename MeshType>
struct is_matrix_basis<scaled_monomial_sym_matrix_basis<MeshType, typename MeshType::face>>
{
   static const bool value = true;
};

template<typename MeshType>
struct is_matrix_basis<Raviart_Thomas_matrix_basis<MeshType, typename MeshType::cell>>
{
   static const bool value = true;
};

template<typename MeshType>
struct is_matrix_basis<Raviart_Thomas_matrix_basis<MeshType, typename MeshType::face>>
{
   static const bool value = true;
};

// use_vector_container basis ?
template<typename Basis>
struct use_vector_container
{
   static const bool value = false;
};

template<typename MeshType>
struct use_vector_container<scaled_monomial_vector_basis<MeshType, typename MeshType::cell>>
{
   static const bool value = true;
};

template<typename MeshType>
struct use_vector_container<scaled_monomial_vector_basis<MeshType, typename MeshType::face>>
{
   static const bool value = true;
};

template<typename MeshType>
struct use_vector_container<scaled_monomial_sym_matrix_basis<MeshType, typename MeshType::cell>>
{
   static const bool value = true;
};

template<typename MeshType>
struct use_vector_container<scaled_monomial_sym_matrix_basis<MeshType, typename MeshType::face>>
{
   static const bool value = true;
};

template<typename MeshType>
struct use_vector_container<Raviart_Thomas_vector_basis<MeshType, typename MeshType::cell>>
{
   static const bool value = true;
};

template<typename MeshType>
struct use_vector_container<Raviart_Thomas_vector_basis<MeshType, typename MeshType::face>>
{
   static const bool value = true;
};

template<typename MeshType>
struct use_vector_container<scaled_monomial_matrix_basis<MeshType, typename MeshType::cell>>
{
   static const bool value = true;
};

template<typename MeshType>
struct use_vector_container<scaled_monomial_matrix_basis<MeshType, typename MeshType::face>>
{
   static const bool value = true;
};

template<typename MeshType>
struct use_vector_container<Raviart_Thomas_matrix_basis<MeshType, typename MeshType::cell>>
{
   static const bool value = true;
};

template<typename MeshType>
struct use_vector_container<Raviart_Thomas_matrix_basis<MeshType, typename MeshType::face>>
{
   static const bool value = true;
};

} // end disk