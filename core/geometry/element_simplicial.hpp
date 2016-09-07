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

namespace priv {

template<size_t DIM, size_t CODIM>
struct howmany;

}

template<size_t DIM, size_t CODIM>
class simplicial_element
{
    typedef point_identifier<DIM>       point_id_type;

    typedef std::array<point_id_type, priv::howmany<DIM, CODIM>::nodes> node_array_type;

    node_array_type     m_pts_ptrs;

public:
    typedef identifier<simplicial_element, ident_impl_t, 0> id_type;

    simplicial_element() = default;

    simplicial_element(std::initializer_list<point_id_type> l)
    {
        std::copy(l.begin(), l.end(), m_pts_ptrs.begin());
        std::sort(m_pts_ptrs.begin(), m_pts_ptrs.end());
    }

    node_array_type point_ids(void) const
    {
        return m_pts_ptrs;
    }

    bool operator<(const simplicial_element& other) const
    {
        return m_pts_ptrs < other.m_pts_ptrs;
    }

    bool operator==(const simplicial_element& other) const
    {
        return m_pts_ptrs == other.m_pts_ptrs;
    }

    size_t subelement_size() const
    {
        return priv::howmany<DIM, CODIM>::subelements;
    }
};

namespace priv {
    template<size_t DIM>
    struct howmany<DIM, DIM>
    {
        static const size_t nodes = 1;
        static const size_t edges = 0;
        static const size_t surfaces = 0;
        static const size_t volumes = 0;
    };

    template<>
    struct howmany<3,0>
    {
        static const size_t nodes = 4;
        static const size_t edges = 6;
        static const size_t surfaces = 4;
        static const size_t volumes = 1;
        static const size_t subelements = surfaces;
    };

    template<>
    struct howmany<3,1>
    {
        static const size_t nodes = 3;
        static const size_t edges = 3;
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
    struct howmany<2,0>
    {
        static const size_t nodes = 3;
        static const size_t edges = 3;
        static const size_t surfaces = 1;
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
    struct howmany<1,0>
    {
        static const size_t nodes = 2;
        static const size_t edges = 1;
        static const size_t surfaces = 0;
        static const size_t volumes = 0;
        static const size_t subelements = nodes;
    };
} // namespace priv

/* Output streaming operator for elements */
template<size_t DIM, size_t CODIM>
std::ostream&
operator<<(std::ostream& os, const simplicial_element<DIM, CODIM>& e)
{
    os << "simplicial_element<" << DIM << "," << CODIM << ">: ";
    auto pts = e.point_ids();
    for (auto itor = pts.begin(); itor != pts.end(); itor++)
        os << *itor << " ";

    return os;
}

} // namespace hho
