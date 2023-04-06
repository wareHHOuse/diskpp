/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2021                nicolas.pignet@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework.
 * M. Abbas, A. Ern, N. Pignet.
 * International Journal of Numerical Methods in Engineering (2019)
 * 120(3), 303-327
 * DOI: 10.1002/nme.6137
 */

#pragma once

#include <string>

namespace disk
{


enum DeformationMeasure : size_t
{
    SMALL_DEF       = 0,
    LOGARITHMIC_DEF = 1,
    F_DEF           = 2,
};

enum LawType : size_t
{
    ELASTIC             = 0,
    LINEAR_HARDENING    = 1,
    NONLINEAR_HARDENING = 2,
    HENCKY_MISES        = 3,
    NEOHOKEAN           = 4,
    CAVITATION          = 5,
    MFRONT              = 6
};

enum MfrontVariableType : size_t
{
    SCALAR  = 0,
    VECTOR  = 1,
    STENSOR = 2,
    TENSOR  = 3,
};

std::string
DeformationMeasureName(const size_t& def)
{
    switch (def)
    {
        case DeformationMeasure::SMALL_DEF: return "SMALL_DEF"; break;
        case DeformationMeasure::LOGARITHMIC_DEF: return "LOGARITHMIC_DEF"; break;
        case DeformationMeasure::F_DEF: return "F_DEF"; break;
        default:
            throw std::invalid_argument("Unknown deformation");
            break;
    }
}

std::string
LawTypeName(const size_t& law)
{
    switch (law)
    {
        case LawType::ELASTIC: return "ELASTIC"; break;
        case LawType::LINEAR_HARDENING: return "LINEAR_HARDENING"; break;
        case LawType::NONLINEAR_HARDENING: return "NONLINEAR_HARDENING"; break;
        case LawType::HENCKY_MISES: return "HENCKY_MISES"; break;
        case LawType::NEOHOKEAN: return "NEOHOKEAN"; break;
        case LawType::CAVITATION: return "CAVITATION"; break;
        case LawType::MFRONT: return "MFRONT"; break;

        default: throw std::invalid_argument("Not known law"); break;
    }
}

std::string
MfrontVariableTypeName(const size_t& var)
{
    switch (var)
    {
        case MfrontVariableType::SCALAR: return "SCALAR"; break;
        case MfrontVariableType::VECTOR: return "VECTOR"; break;
        case MfrontVariableType::STENSOR: return "STENSOR"; break;
        case MfrontVariableType::TENSOR: return "TENSOR"; break;

        default: throw std::invalid_argument("Not known variable"); break;
    }
}
};