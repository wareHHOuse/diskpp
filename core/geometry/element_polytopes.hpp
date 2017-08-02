/*
 *       /\        Matteo Cicuttin (C) 2016-2017
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
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#include <iostream>
#include <array>

//#include "diskpp_types.hpp"

namespace disk {

typedef uint32_t diskpp_index_t;
typedef uint32_t diskpp_unsigned_t;

template<typename UserData>
class userdata_stuff
{
public:
    UserData    user_data;
};

template<>
class userdata_stuff<void>
{};


template<size_t DIM, typename polytope_tag, typename UserData>
class polytope
{
    static_assert( sizeof(polytope_tag) == -1, "Unknown polytope" );
};


template<size_t DIM, size_t VTX, typename UserData>
class fixed_size_element : public userdata_stuff<UserData>
{
    std::array<diskpp_index_t, VTX>   m_point_ids;

public:
    fixed_size_element() = default;

    fixed_size_element(std::initializer_list<diskpp_index_t> l)
    {
        std::copy(l.begin(), l.end(), m_point_ids.begin());
        std::sort(m_point_ids.begin(), m_point_ids.end());
    }

    std::array<diskpp_index_t, VTX>
    point_ids(void) const
    {
        return m_point_ids;
    }

    bool operator<(const fixed_size_element& other) const
    {
        return m_point_ids < other.m_point_ids;
    }

    bool operator==(const fixed_size_element& other) const
    {
        return m_point_ids == other.m_pts_ptrs;
    }
};

/* per fare le conversioni uso std::transform */

struct polytope_simplex;

constexpr size_t simplex_vertices(size_t DIM) { return DIM+1; }

template<size_t DIM, typename UserData>
class polytope<DIM, polytope_simplex, UserData>
    : public fixed_size_element<DIM, simplex_vertices(DIM), UserData>
{
    typedef fixed_size_element<DIM, simplex_vertices(DIM), UserData> base;
    using base::base;
};

template<size_t DIM, typename UserData = void>
using simplicial_element = polytope<DIM, polytope_simplex, UserData>;

template<size_t DIM, typename UserData>
std::ostream&
operator<<(std::ostream& os, const polytope<DIM, polytope_simplex, UserData>& s)
{
    os << "simplex<" << DIM << ">: ";
    auto pts = s.point_ids();
    for (auto itor = pts.begin(); itor != pts.end(); itor++)
        os << *itor << " ";

    return os;
}

} // namespace disk
