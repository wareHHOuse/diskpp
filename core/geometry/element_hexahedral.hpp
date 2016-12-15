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

#include <vector>
#include <array>
#include <cassert>

#include "mesh/mesh.hpp"
#include "mesh/mesh_storage.hpp"
#include "mesh/ident.hpp"
#include "mesh/point.hpp"

namespace disk {

namespace cartesian_priv {

template<size_t DIM, size_t CODIM>
struct howmany;

}

template<size_t DIM, size_t CODIM>
class cartesian_element
{
    static_assert(DIM == 2 or DIM == 3, "cartesian elements must be 2D or 3D");

    typedef point_identifier<DIM>       point_id_type;

    typedef std::array<point_id_type, cartesian_priv::howmany<DIM, CODIM>::nodes>
        node_array_type;

    node_array_type     m_pts_ptrs;

public:
    typedef identifier<cartesian_element, ident_impl_t, 0> id_type;

    cartesian_element() = default;

    cartesian_element(std::initializer_list<point_id_type> l)
    {
        std::copy(l.begin(), l.end(), m_pts_ptrs.begin());
    }

    node_array_type point_ids(void) const
    {
        return m_pts_ptrs;
    }

    bool operator<(const cartesian_element& other) const
    {
        return m_pts_ptrs < other.m_pts_ptrs;
    }

    bool operator==(const cartesian_element& other) const
    {
        return m_pts_ptrs == other.m_pts_ptrs;
    }

    size_t subelement_size() const
    {
        return cartesian_priv::howmany<DIM, CODIM>::subelements;
    }
};

namespace cartesian_priv {
    template<>
    struct howmany<3,0>
    {
        static const size_t nodes = 8;
        static const size_t edges = 12;
        static const size_t surfaces = 6;
        static const size_t volumes = 1;
        static const size_t subelements = surfaces;
    };

    template<>
    struct howmany<3,1>
    {
        static const size_t nodes = 4;
        static const size_t edges = 4;
        static const size_t surfaces = 1;
        static const size_t volumes = 0;
        static const size_t subelements = edges;
    };

    template<>
    struct howmany<3,2>
    {
        static const size_t nodes = 2;
        static const size_t edges = 1;
        static const size_t surfaces = 0;
        static const size_t volumes = 0;
        static const size_t subelements = nodes;
    };

    template<>
    struct howmany<3,3>
    {
        static const size_t nodes = 1;
        static const size_t edges = 0;
        static const size_t surfaces = 0;
        static const size_t volumes = 0;
    };

    template<>
    struct howmany<2,0>
    {
        static const size_t nodes = 4;
        static const size_t edges = 4;
        static const size_t surfaces = 4;
        static const size_t volumes = 0;
        static const size_t subelements = edges;
    };

    template<>
    struct howmany<2,1>
    {
        static const size_t nodes = 2;
        static const size_t edges = 1;
        static const size_t surfaces = 0;
        static const size_t volumes = 0;
        static const size_t subelements = nodes;
    };

    template<>
    struct howmany<2,2>
    {
        static const size_t nodes = 1;
        static const size_t edges = 0;
        static const size_t surfaces = 0;
        static const size_t volumes = 0;
    };
} // namespace priv

/* Output streaming operator for elements */
template<size_t DIM, size_t CODIM>
std::ostream&
operator<<(std::ostream& os, const cartesian_element<DIM, CODIM>& e)
{
    os << "cartesian_element<" << DIM << "," << CODIM << ">: ";
    auto pts = e.point_ids();
    for (auto itor = pts.begin(); itor != pts.end(); itor++)
        os << *itor << " ";

    return os;
}

} // namespace hho
