/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
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

#include "Cavitation/Cavitation_qp.hpp"
#include "HenckyMises/HenckyMises_qp.hpp"
#include "LinearIsotropicAndKinematicHardening/LinearIsotropicAndKinematicHardening_qp.hpp"
#include "LinearLaw/LinearLaw_qp.hpp"
#include "Neohookean/Neohookean_qp.hpp"
#include "law_bones.hpp"

namespace disk
{

template<typename MeshType>
using Cavitation = LawTypeBones<MeshType, Cavitation_qp<typename MeshType::coordinate_type, MeshType::dimension>, false>;

template<typename MeshType>
using Neohookean =
  LawTypeBones<MeshType, Neohookean_qp<typename MeshType::coordinate_type, MeshType::dimension>, false>;

template<typename MeshType>
using HenckyMises =
  LawTypeBones<MeshType, HenckyMises_qp<typename MeshType::coordinate_type, MeshType::dimension>, false>;

template<typename MeshType>
using LinearLaw = LawTypeBones<MeshType, LinearLaw_qp<typename MeshType::coordinate_type, MeshType::dimension>, false>;

template<typename MeshType>
using LinearIsotropicAndKinematicHardening =
  LawTypeBones<MeshType,
               LinearIsotropicAndKinematicHardening_qp<typename MeshType::coordinate_type, MeshType::dimension>,
               true>;
}
