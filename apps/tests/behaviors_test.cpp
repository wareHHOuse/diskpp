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

#include "diskpp/mechanics/behaviors/laws/Cavitation/Cavitation_qp.hpp"
#include "diskpp/mechanics/behaviors/laws/HenckyMises/HenckyMises_qp.hpp"
#include "diskpp/mechanics/behaviors/laws/IsotropicHardeningVMis/IsotropicHardeningVMis_qp.hpp"
#include "diskpp/mechanics/behaviors/laws/LinearIsotropicAndKinematicHardening/LinearIsotropicAndKinematicHardening_qp.hpp"
#include "diskpp/mechanics/behaviors/laws/LinearLaw/LinearLaw_qp.hpp"
#include "diskpp/mechanics/behaviors/laws/Neohookean/Neohookean_qp.hpp"

// this is to test that the compilation of the different behaviors is ok

int
main(int argc, char const* argv[])
{
    using T = double;

    disk::mechanics::Cavitation_qp< T, 3 > cav3D;
    disk::mechanics::Cavitation_qp< T, 2 > cav2D;
    disk::mechanics::HenckyMises_qp< T, 3 > Hen3D;
    disk::mechanics::HenckyMises_qp< T, 2 > Hen2D;
    disk::mechanics::IsotropicHardeningVMis_qp< T, 3 > isot3D;
    disk::mechanics::IsotropicHardeningVMis_qp< T, 2 > isot2D;
    disk::mechanics::LinearIsotropicAndKinematicHardening_qp< T, 3 > mix3D;
    disk::mechanics::LinearIsotropicAndKinematicHardening_qp< T, 2 > mix2D;
    disk::mechanics::LinearLaw_qp< T, 3 > lin3D;
    disk::mechanics::LinearLaw_qp< T, 2 > lin2D;
    disk::mechanics::Neohookean_qp< T, 3 > neo3D;
    disk::mechanics::Neohookean_qp< T, 2 > neo2D;

    return 0;
}
