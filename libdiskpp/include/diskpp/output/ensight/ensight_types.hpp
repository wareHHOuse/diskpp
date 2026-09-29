/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2026
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

#include <array>
#include <cstddef>
#include <cstdint>
#include <set>
#include <string>
#include <vector>

namespace disk::output::ensight {

enum class FieldType { scalar, vector, symmetric_tensor, tensor };
enum class FieldLocation { node, element };
enum class ElementBlock { bar2, tria3, quad4, nsided, nfaced };

inline const char *
case_keyword( const FieldType type, const FieldLocation location ) {
    if ( location == FieldLocation::node ) {
        switch ( type ) {
        case FieldType::scalar:
            return "scalar per node";
        case FieldType::vector:
            return "vector per node";
        case FieldType::symmetric_tensor:
            return "tensor symm per node";
        case FieldType::tensor:
            return "tensor asym per node";
        }
    } else {
        switch ( type ) {
        case FieldType::scalar:
            return "scalar per element";
        case FieldType::vector:
            return "vector per element";
        case FieldType::symmetric_tensor:
            return "tensor symm per element";
        case FieldType::tensor:
            return "tensor asym per element";
        }
    }
    return "";
}

struct ElementReference {
    ElementBlock block = ElementBlock::bar2;
    std::size_t index_in_block = 0;
};

struct NfacedCell {
    std::vector< std::vector< std::int64_t > > faces;
};

struct PointCloud {
    std::string name;
    std::int64_t part_id = 2;
    std::vector< std::array< double, 3 > > points;
    std::size_t
    number_of_points() const noexcept {
        return points.size();
    }
};

struct Geometry {
    std::string part_name = "HHO mesh";
    std::int64_t part_id = 1;
    std::vector< std::array< double, 3 > > points;
    std::vector< std::array< std::int64_t, 2 > > bars;
    std::vector< std::array< std::int64_t, 3 > > triangles;
    std::vector< std::array< std::int64_t, 4 > > quadrilaterals;
    std::vector< std::vector< std::int64_t > > polygons;
    std::vector< NfacedCell > polyhedra;
    std::vector< ElementReference > element_reference;

    std::size_t
    number_of_points() const noexcept {
        return points.size();
    }
    std::size_t
    number_of_elements() const noexcept {
        return element_reference.size();
    }
};

struct FieldDescription {
    std::string name;
    std::string safe_name;
    FieldType type = FieldType::scalar;
    FieldLocation location = FieldLocation::node;
    std::int64_t part_id = 1;
    std::set< std::size_t > written_steps;
};

} // namespace disk::output::ensight
