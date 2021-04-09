/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet (C) 2019                      nicolas.pignet@enpc.fr
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

#include "mesh/mesh.hpp"
#include "mesh/mesh_storage.hpp"
#include "mesh/ident.hpp"
#include "mesh/point.hpp"

namespace disk {

template<size_t DIM, size_t CODIM>
class generic_element;

/* type traits for the elements */
template<typename T>
struct generic_element_traits;

template<size_t DIM, size_t CODIM>
struct generic_element_traits<generic_element<DIM, CODIM>>
{
    static_assert(CODIM < DIM, "generic_element_traits: CODIM must be less than DIM");

    typedef generic_element<DIM, CODIM+1>                           subelement_type;
    typedef identifier<generic_element<DIM,CODIM>, ident_raw_t, 0> id_type;
    static const size_t dimension = DIM;
    static const size_t codimension = CODIM;
};

template<size_t DIM>
struct generic_element_traits<generic_element<DIM, DIM>>
{
    typedef identifier<generic_element<DIM,DIM>, ident_raw_t, 0>   id_type;
    static const size_t dimension = DIM;
    static const size_t codimension = DIM;
};

template<typename To, typename From>
std::vector<To>
convert_to(const std::vector<From>& vec)
{
    std::vector<To> ret;
    ret.reserve(vec.size());
    for (auto& v : vec)
        ret.push_back( To(v) );

    return ret;
}


/* element base class (for subelements) */
template<typename T>
class generic_element_base
{
    typedef typename generic_element_traits<T>::subelement_type::id_type    sub_id_type;
    typedef point_identifier<generic_element_traits<T>::dimension>          point_id_type;

    void
    set_point_ids()
    {
       m_pts_ptrs.clear();
       m_pts_ptrs.reserve(m_sids_ptrs.size());

       for (auto& sid : m_sids_ptrs)
          m_pts_ptrs.push_back(point_id_type(sid));
    }

  protected:
    std::vector<sub_id_type>            m_sids_ptrs;
    std::vector<point_id_type>          m_pts_ptrs;

public:
    typedef typename std::vector<sub_id_type>::iterator         sid_iterator;
    typedef typename std::vector<sub_id_type>::const_iterator   const_sid_iterator;

    generic_element_base()
    {}

    generic_element_base(std::initializer_list<sub_id_type> l)
        : m_sids_ptrs(l)
    {
        set_point_ids();
    }

    explicit generic_element_base(const std::vector<sub_id_type>& sids_ptrs)
        : m_sids_ptrs(sids_ptrs)
    {
       set_point_ids();
    }

    explicit generic_element_base(std::vector<sub_id_type>&& sids_ptrs)
        : m_sids_ptrs( std::move(sids_ptrs) )
    {
       set_point_ids();
    }

    template<typename Itor>
    void set_subelement_ids(const Itor e_begin, const Itor e_end)
    {
        m_sids_ptrs = std::vector<sub_id_type>(e_begin, e_end);
    }

    /* All this point-related stuff is simply horrible, because points should
     * be automatically derived from lower-dimensional elements. However,
     * this is not possible because points must be stored in a specific order,
     * otherwise geometry-related computations fail.
     *
     * No, in 2D they won't fail! But what about 3D?
     */
    template<typename Itor>
    void set_point_ids(const Itor p_begin, const Itor p_end)
    {
        m_pts_ptrs = std::vector<point_id_type>(p_begin, p_end);
    }

    void set_point_ids(std::vector<point_id_type>&& pts)
    {
        m_pts_ptrs = std::move(pts);
    }

    std::vector<point_id_type>
    point_ids(void) const
    {
        assert( m_pts_ptrs.size() > 0 );
        return m_pts_ptrs;
    }

    std::vector<sub_id_type>
    faces_ids(void) const
    {
        assert(m_sids_ptrs.size() > 0);
        return m_sids_ptrs;
    }

    sid_iterator        subelement_id_begin() { return m_sids_ptrs.begin(); }
    sid_iterator        subelement_id_end()   { return m_sids_ptrs.end(); }
    const_sid_iterator  subelement_id_begin() const { return m_sids_ptrs.begin(); }
    const_sid_iterator  subelement_id_end()   const { return m_sids_ptrs.end(); }

    size_t              subelement_size() const { return m_sids_ptrs.size(); }
};

/* element class, CODIM < DIM case */
template<size_t DIM, size_t CODIM>
class generic_element : public generic_element_base<generic_element<DIM, CODIM>>
{
    static_assert(CODIM <= DIM, "generic_element: CODIM > DIM does not make sense");

public:
    typedef typename generic_element_traits<generic_element>::subelement_type  subelement_type;
    typedef typename generic_element_traits<generic_element>::id_type          id_type;

    generic_element()
    {}

    generic_element(std::initializer_list<typename subelement_type::id_type> l)
        : generic_element_base<generic_element<DIM, CODIM>>(l)
    {}

    explicit generic_element(const std::vector<typename subelement_type::id_type>& sids_ptrs)
        : generic_element_base<generic_element<DIM, CODIM>>(sids_ptrs)
    {}

    bool operator<(const generic_element& other) const
    {
        return (this->m_sids_ptrs < other.m_sids_ptrs);
    }

    bool operator==(const generic_element& other) const
    {
        return std::equal(this->m_sids_ptrs.begin(), this->m_sids_ptrs.end(),
                          other.m_sids_ptrs.begin());
    }
};

/* NODE */
template<size_t DIM>
class generic_element<DIM, DIM>
{
    point_identifier<DIM>    m_assoc_point;

public:
    typedef typename generic_element_traits<generic_element>::id_type          id_type;

    generic_element()
    {}

    explicit generic_element(point_identifier<DIM> assoc_point)
    {
        m_assoc_point = assoc_point;
    }

    std::vector<point_identifier<DIM>>
    point_ids(void) const
    {
        return std::vector<point_identifier<DIM>>{{m_assoc_point}};
    }

    bool operator<(const generic_element& other) const
    {
        return this->m_assoc_point < other.m_assoc_point;
    }

    bool operator==(const generic_element& other) const
    {
        return this->m_assoc_point == other.m_assoc_point;
    }
};

/* EDGE */
template<size_t DIM, typename Derived>
class generic_element_edge
{
public: 
    using sub_id_type = typename generic_element_traits<Derived>::subelement_type::id_type;
    using ptid_type = point_identifier<DIM>;

private:
    std::array<sub_id_type, 2>  sub_elem_ids;
    std::array<ptid_type, 2>    assoc_points;

public:
    typedef typename generic_element_traits<Derived>::id_type          id_type;

    generic_element_edge()
    {}

    generic_element_edge(const sub_id_type& n0, const sub_id_type& n1)
    {
        /* Assume node ids == point ids. Use set_point_ids() below if that
         * is not the case. */
        if (n0 < n1)
        {
            sub_elem_ids[0] = n0;
            sub_elem_ids[1] = n1;
        }
        else
        {
            sub_elem_ids[0] = n1;
            sub_elem_ids[1] = n0;
        }

        assoc_points[0] = ptid_type(n0);
        assoc_points[1] = ptid_type(n1);
    }

    
    template<typename Itor>
    void set_point_ids(const ptid_type& p0, const ptid_type& p1)
    {
        assoc_points[0] = p0;
        assoc_points[1] = p1;
    }

    auto point_ids(void) const
    {
        return assoc_points;
    }

    auto faces_ids(void) const
    {
        return sub_elem_ids;
    }

    bool operator<(const generic_element_edge& other) const
    {
        return (sub_elem_ids < other.sub_elem_ids);
    }

    bool operator==(const generic_element_edge& other) const
    {
        return std::equal(sub_elem_ids.begin(), sub_elem_ids.end(),
                          other.sub_elem_ids.begin());
    }

    auto subelement_id_begin() { return sub_elem_ids.begin(); }
    auto subelement_id_end()   { return sub_elem_ids.end(); }
    auto subelement_id_begin() const { return sub_elem_ids.begin(); }
    auto subelement_id_end()   const { return sub_elem_ids.end(); }

    size_t subelement_size() const { return sub_elem_ids.size(); }
};

template<>
class generic_element<3,2> : public generic_element_edge<3, generic_element<3,2>>
{
    using base = generic_element_edge<3, generic_element<3,2>>;

public:
    generic_element()
        : base()
    {}

    generic_element(const base::sub_id_type& n0, const base::sub_id_type& n1)
        : base(n0, n1)
    {}
};

template<>
class generic_element<2,1> : public generic_element_edge<2, generic_element<2,1>>
{
    using base = generic_element_edge<2, generic_element<2,1>>;

public:
    generic_element()
        : base()
    {}

    generic_element(const base::sub_id_type& n0, const base::sub_id_type& n1)
        : base(n0, n1)
    {}
};

template<>
class generic_element<1,0> : public generic_element_edge<1, generic_element<1,0>>
{
    using base = generic_element_edge<1, generic_element<1,0>>;

public:
    generic_element()
        : base()
    {}

    generic_element(const base::sub_id_type& n0, const base::sub_id_type& n1)
        : base(n0, n1)
    {}
};




/* Output streaming operator for elements */
template<size_t DIM, size_t CODIM>
std::ostream&
operator<<(std::ostream& os, const generic_element<DIM, CODIM>& e)
{
    os << "generic_element<" << DIM << "," << CODIM << ">: ";
    for (auto itor = e.subelement_id_begin();
         itor != e.subelement_id_end();
         itor++)
        os << *itor << " ";

    os << "( ";
    auto pts = e.point_ids();
    for (auto& pt : pts)
        os << pt << " ";

    std::cout << ")";

    return os;
}

template<size_t DIM>
std::ostream&
operator<<(std::ostream& os, const generic_element<DIM, DIM>& e)
{
    os << "element<" << DIM << "," << DIM << ">: ";
    return os;
}

} // namespace disk
