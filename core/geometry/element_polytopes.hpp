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
#include <vector>
#include <array>
#include <cassert>
#include <type_traits>

//#include "core/mesh/point.hpp"

typedef size_t  domain_id_type;
typedef size_t  boundary_id_type;
typedef size_t  point_id_type;

namespace disk {

namespace priv {

/**
 The userdata class allows to add some user-specified data structure to
 the mesh elements. Normally an element is instantiated with UserData = void,
 however having some user-specified data can be useful in case of mesh
 refinement or similar stuff. In that case, one can attach flags and values
 to the elements.
 */
template<typename UserData>
class userdata
{
public:
    UserData    user_data;
};

template<>
class userdata<void>
{};


/**
 Here we have the element_info. The rationale here is that an element of
 codimension 0 is always a cell and an element of codimension 1 is always
 a face. Based on the codimension then, we are going to attach to each element
 information about the domain or the boundary it belongs to.
 */
template<size_t CODIM, typename UserData>
class element_info : public userdata<UserData>
{};

template<typename UserData>
class element_info<0, UserData> : public userdata<UserData>
{
    domain_id_type      m_domain_id;

public:
    element_info()
        : m_domain_id(0)
    {}

    domain_id_type domain_id() const
    {
        return m_domain_id;
    }

    void domain_id(domain_id_type d_id)
    {
        m_domain_id = d_id;
    }
};

template<typename UserData>
class element_info<1, UserData> : public userdata<UserData>
{
    boundary_id_type    m_boundary_id;
    bool                m_is_boundary;

public:
    element_info()
        : m_boundary_id(0), m_is_boundary(false)
    {}

    bool is_boundary() const
    {
        return m_is_boundary;
    }

    void is_boundary(bool b)
    {
        m_is_boundary = b;
    }

    std::pair<boundary_id_type, bool> boundary_id() const
    {
        return std::make_pair(m_boundary_id, m_is_boundary);
    }

    void boundary_id(boundary_id_type b_id)
    {
        m_boundary_id = b_id;
        m_is_boundary = true;
    }
};

template<typename UserDataS, typename UserDataD>
void
copy_element_info(const element_info<0, UserDataS>& src, element_info<0, UserDataD>& dst)
{
    dst.domain_id( src.domain_id() );
}

template<typename UserDataS, typename UserDataD>
void
copy_element_info(const element_info<1, UserDataS>& src, element_info<1, UserDataD>& dst)
{
    dst.boundary_id( src.boundary_id() );
    dst.is_boundary( src.is_boundary() );
}

template<size_t CODIM, typename UserDataS, typename UserDataD>
void
copy_element_info(const element_info<CODIM, UserDataS>&, element_info<CODIM, UserDataD>&)
{}

} //namespace priv

template<size_t N>                  /* Tag to identify fixed-size polytopes */
struct fixed_storage_polytope;

template<size_t N>                  /* Tag to identify variable-size fixed-storage polytopes */
struct limited_storage_polytope;

struct dynamic_storage_polytope;    /* Tag to identify fully dynamic storage polytopes */

namespace priv {

template<size_t DIM, size_t CODIM, typename UserData, typename StoragePolicy>
class polytope;

template<size_t DIM, size_t CODIM, typename UserData, size_t N>
class polytope<DIM, CODIM, UserData, fixed_storage_polytope<N>> : public element_info<CODIM, UserData>
{
    std::array<point_id_type, N>    m_pts;

public:
    polytope()
    {}

    polytope(std::initializer_list<point_id_type> l)
    {
        assert(l.size() == N);
        std::copy(l.begin(), l.end(), m_pts.begin());
    }

    /* Userdata-changing constructor */
    template<typename OtherUD>
    explicit polytope(const polytope<DIM, CODIM, OtherUD, fixed_storage_polytope<N>>& other)
    {
        priv::copy_element_info(other, *this);
        m_pts = other.point_identifiers();
    }

    std::array<point_id_type, N>
    point_identifiers() const
    {
        return m_pts;
    }

    bool operator<(const polytope& other)
    {
        return (m_pts < other.m_pts);
    }

    bool operator==(const polytope& other)
    {
        return (m_pts == other.m_pts);
    }
};

template<size_t DIM, size_t CODIM, typename UserData, size_t N>
class polytope<DIM, CODIM, UserData, limited_storage_polytope<N>> : public priv::element_info<CODIM, UserData>
{
    std::array<point_id_type, N>    m_pts;
    size_t                          m_actual_pts;

public:
    polytope()
    {}

    polytope(std::initializer_list<point_id_type> l)
    {
        assert(l.size() >= 3 && l.size() <= N);
        std::copy(l.begin(), l.end(), m_pts.begin());
        m_actual_pts = l.size();
    }

    /* Userdata-changing constructor */
    template<typename OtherUD>
    explicit polytope(const polytope<DIM, CODIM, OtherUD, limited_storage_polytope<N>>& other)
    {
        priv::copy_element_info(other, *this);
        auto other_pts = other.point_identifiers();
        assert(other_pts.size() <= m_pts.size());
        std::copy(other_pts.begin(), other_pts.end(), m_pts.begin());
        m_actual_pts = other_pts.size();
    }

    std::vector<point_id_type>
    point_identifiers() const
    {
        assert(m_actual_pts <= m_pts.size());
        return std::vector<point_id_type>(m_pts.begin(), m_pts.begin() + m_actual_pts);
    }

    bool operator<(const polytope& other)
    {
        return (m_pts < other.m_pts);
    }

    bool operator==(const polytope& other)
    {
        return (m_pts == other.m_pts);
    }
};

} // namespace priv


/**
 The actual definitions of the elements.
 */

/******************************************************************************/
template<size_t CODIM, typename UserData>
using node = priv::polytope<0, CODIM, UserData, fixed_storage_polytope<1>>;

template<size_t CODIM, typename UserData>
using edge = priv::polytope<1, CODIM, UserData, fixed_storage_polytope<2>>;

template<size_t CODIM, typename UserData, typename StoragePolicy>
using surface = priv::polytope<2, CODIM, UserData, StoragePolicy>;


/******************************************************************************/
template<size_t DIM, size_t CODIM, typename UserData, typename StorageScheme>
std::ostream& operator<<(std::ostream& os, const priv::polytope<DIM, CODIM, UserData, StorageScheme>& s)
{
    switch(DIM)
    {
        case 0: os << "node<" << CODIM << ">: "; break;
        case 1: os << "edge<" << CODIM << ">: "; break;
        case 2: os << "surface<" << CODIM << ">: "; break;
        case 3: os << "volume<" << CODIM << ">: "; break;
    }
    
    auto ptids = s.point_identifiers();
    for (auto& ptid : ptids)
        os << ptid << " ";

    if ( !std::is_same<UserData, void>::value )
        os << "[UD] ";
    
    return os;    
}
/*
template<typename UserData, typename StorageScheme>
std::ostream& operator<<(std::ostream& os, const surface<1, UserData, StorageScheme>& s)
{
    os << "surface<" << 1 << ">: ";
    auto ptids = s.point_identifiers();
    for (auto& ptid : ptids)
        os << ptid << " ";

    if ( !std::is_same<UserData, void>::value )
        os << "[UD] ";

    if ( s.is_boundary() )
    {
        auto bid = s.boundary_id();
        os << "[B " << bid.first << "] ";
    }
    
    return os;    
}

template<typename UserData, typename StorageScheme>
std::ostream& operator<<(std::ostream& os, const surface<0, UserData, StorageScheme>& s)
{
    os << "surface<" << 0 << ">: ";
    auto ptids = s.point_identifiers();
    for (auto& ptid : ptids)
        os << ptid << " ";

    if ( !std::is_same<UserData, void>::value )
        os << "[UD] ";

    os << "[D " << s.domain_id() << "] ";
    
    return os;    
}
*/

/******************************************************************************/


/******************************************************************************/
template<size_t CODIM, typename UserData>
class volume : public priv::element_info<CODIM, UserData>
{};




/******************************************************************************/
/*
template<size_t CODIM, typename UserData>
using triangle = surface<CODIM, UserData, fixed_polygon<3>>;

template<size_t CODIM, typename UserData>
using quadrilateral = surface<CODIM, UserData, fixed_polygon<4>>;

template<size_t CODIM, typename UserData>
using arbitrary_polygon = surface<CODIM, UserData, dynamic_polygon>;
*/


} // namespace disk
