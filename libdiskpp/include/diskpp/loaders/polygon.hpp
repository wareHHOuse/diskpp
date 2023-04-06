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

#include <vector>
#include <iostream>

namespace disk{

class polygon
{
    std::vector<size_t>::const_iterator   m_my_position;

public:
    polygon (){}
    polygon(std::vector<size_t>::const_iterator my_position)
        : m_my_position(my_position)
    {}

    size_t
    number_of_edges(void) const
    {
        return *m_my_position;
    }

    std::vector<size_t>
    point_ids(void) const
    {
        auto begin = m_my_position;
        auto end = m_my_position;
        auto num_edges = number_of_edges();
        std::advance(begin, 1);
        std::advance(end, num_edges+1);
        return std::vector<size_t>(begin, end);
    }
};


namespace priv {

class polygon_iterator
{
    std::vector<size_t>::const_iterator   m_itor, m_end;

    void advance(void)
    {
        if (m_itor == m_end)
            return;

        size_t next_poly_ofs = *m_itor + 1;
        std::advance(m_itor, next_poly_ofs);
    }

public:

    typedef std::vector<size_t>::const_iterator const_iterator;

    polygon_iterator() {}

    polygon_iterator(const_iterator begin, const_iterator end)
        : m_itor(begin), m_end(end)
    {}

    polygon operator*() { return polygon(m_itor); }

    polygon_iterator& operator++()
    {
        if ( m_itor != m_end )
            advance();

        return *this;
    }

    polygon_iterator operator++(int)
    {
        auto it = *this;
        ++(*this);
        return it;
    }

    bool operator==(const polygon_iterator& other) const
    {
        return (m_itor == other.m_itor);
    }

    bool operator!=(const polygon_iterator& other) const
    {
        return (m_itor != other.m_itor);
    }

};

} // namespace priv


std::ostream&
operator<<(std::ostream& os, const polygon& poly)
{
    os << "[polygon] " << poly.number_of_edges() << " ( ";
    auto ptids = poly.point_ids();
    for (auto& pi : ptids)
        os << pi << " ";
    os << ")";
    return os;
}

class polygon_store
{
    std::vector<size_t>     m_data;

public:
    polygon_store() = default;

    void add_polygon(const std::vector<size_t>& point_ids)
    {
        m_data.push_back(point_ids.size());
        m_data.insert(m_data.end(), point_ids.begin(), point_ids.end());
    }

    priv::polygon_iterator begin() const
    {
        return priv::polygon_iterator(m_data.begin(), m_data.end());
    }

    priv::polygon_iterator end() const
    {
        return priv::polygon_iterator(m_data.end(), m_data.end());
    }
};

} // end disk
