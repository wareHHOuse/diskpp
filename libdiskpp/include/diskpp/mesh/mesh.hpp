/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023, 2024
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */
/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2020, 2021
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group
 */
/*
 *       /\        Matteo Cicuttin (C) 2016-2019
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

 /*
  * This source file is part of EMT, the ElectroMagneticTool.
  *
  * Copyright (C) 2013-2015, Matteo Cicuttin - matteo.cicuttin@uniud.it
  * Department of Electrical Engineering, University of Udine
  * All rights reserved.
  *
  * Redistribution and use in source and binary forms, with or without
  * modification, are permitted provided that the following conditions are met:
  *     * Redistributions of source code must retain the above copyright
  *       notice, this list of conditions and the following disclaimer.
  *     * Redistributions in binary form must reproduce the above copyright
  *       notice, this list of conditions and the following disclaimer in the
  *       documentation and/or other materials provided with the distribution.
  *     * Neither the name of the University of Udine nor the
  *       names of its contributors may be used to endorse or promote products
  *       derived from this software without specific prior written permission.
  *
  * THIS SOFTWARE IS PROVIDED BY THE AUTHOR(s) ``AS IS'' AND ANY
  * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  * DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(s) BE LIABLE FOR ANY
  * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  */


#pragma once

#include <iostream>
#include <vector>
#include <memory>
#include <numeric>
#include <cassert>
#include <iterator>
#include <set>
#include <optional>

#include "ident.hpp"

#include "mesh_storage.hpp"

namespace disk {

namespace priv {
    template<typename Iterator, typename Predicate>
    class filter_iterator
    {
        Predicate       m_predicate;
        Iterator        m_itor{}, m_end{};

        void advance(void) { m_itor++; }

        void find_next(void)
        {
            while ( m_itor != m_end )
            {
                if ( m_predicate(*m_itor) )
                    return;
                advance();
            }
        }

    public:
        using value_type = typename std::iterator_traits<Iterator>::value_type;
        using reference = typename std::iterator_traits<Iterator>::reference;

        filter_iterator() = default;

        filter_iterator(Predicate pred, Iterator itor, Iterator end)
            : m_predicate(pred), m_itor(itor), m_end(end)
        {
            if (itor != end)
                find_next();
        }

        reference operator*() { return *m_itor; }
        value_type *operator->() const { return &(*m_itor); }

        filter_iterator& operator++() {
            if ( m_itor != m_end )
                advance();
            find_next();
            return *this;
        }

        filter_iterator operator++(int) {
            auto it = *this;
            ++(*this);
            return it;
        }

        bool operator==(const filter_iterator& other) const { return (m_itor == other.m_itor); }
        bool operator!=(const filter_iterator& other) const { return (m_itor != other.m_itor); }
    };

    template<typename mesh_type>
    class is_boundary_pred
    {
        const mesh_type&    msh_;
    public:
        is_boundary_pred(const mesh_type& msh) : msh_(msh) {}

        template<typename T>
        bool operator()(const T& elem) { return msh_.is_boundary(elem); }
    };

    template<typename mesh_type>
    class is_boundary_pred_with_id
    {
        const mesh_type&    msh_;
        size_t              m_element_id;
    public:
        is_boundary_pred_with_id(const mesh_type& msh, size_t id)
            : msh_(msh), m_element_id(id)
        {}

        template<typename T>
        bool operator()(const T& elem)
        {
            auto eid = msh_.lookup(elem);

            return msh_.is_boundary(eid) and (msh_.boundary_id(eid) == m_element_id);
        }
    };

    template<typename mesh_type>
    class is_internal_pred
    {
        const mesh_type&    msh_;
    public:
        is_internal_pred(const mesh_type& msh) : msh_(msh) {}

        template<typename T>
        bool operator()(const T& elem) { return !msh_.is_boundary(elem); }
    };

} //namespace priv

template<typename T, typename Iterator>
std::pair<bool, typename T::id_type>
find_element_id(const Iterator& begin, const Iterator& end, const T& element)
{
    auto itor = std::lower_bound(begin, end, element);

    if (itor != end && !(element < *itor))
    {
        typename T::id_type id(std::distance(begin, itor));
        return std::make_pair(true, id);
    }

    return std::make_pair(false, typename T::id_type());
}

/****************************************************************************/
namespace priv {

/* Mesh bones.
 *
 * @DIM     Space dimension
 * @T       Type that the mesh uses to represent points.
 */
template<typename T, size_t DIM, typename Storage>
class mesh_bones
{
    static_assert(DIM > 0 && DIM <= 3, "mesh: Allowed dimensions are 1, 2 and 3");

    std::shared_ptr<Storage>    m_storage;

public:
    mesh_bones() { m_storage = std::make_shared<Storage>(); }

    /* Return a shared_ptr to the backend storage. */
    std::shared_ptr<Storage>
    backend_storage(void)
    {
        return m_storage;
    }

    /* Return a shared_ptr to the backend storage. */
    const std::shared_ptr<Storage>
    backend_storage(void) const
    {
        return m_storage;
    }
};

/* \brif Generic template for a mesh.
 *
 * This template has to be specialized for the 1D, 2D and 3D cases and it
 * represents the actual interface between the user and the mesh. It is in
 * `priv` and `mesh` inherits from `mesh_base` to provide an additional
 * decoupling layer.
 * The user should interact with the mesh in terms of cells and faces only.
 *
 * \tparam DIM     Space dimension
 * \tparam T       Type that the mesh uses to represent points.
 * \tparam Storage Type of storage to represent the data
 */
template<typename T, size_t DIM, typename Storage>
class mesh_base
{
    static_assert(DIM > 0 && DIM <= 3, "mesh: Allowed dimensions are 1, 2 and 3");
};



/* \brief Template specialization for 3D meshes.
 *
 * \tparam T       Type that the mesh uses to represent points.
 * \tparam Storage Type of storage to represent the data
 */
template<typename T, typename Storage>
class mesh_base<T,3,Storage> : public mesh_bones<T,3,Storage>
{
public:
    typedef typename Storage::volume_type                       volume_type;
    typedef typename Storage::surface_type                      surface_type;
    typedef typename Storage::edge_type                         edge_type;
    typedef typename Storage::node_type                         node_type;
    typedef typename Storage::point_type                        point_type;

    typedef volume_type                                         cell;
    typedef surface_type                                        face;
    const static size_t dimension = 3;

    /* cell iterators */
    typedef typename std::vector<volume_type>::iterator         cell_iterator;
    typedef typename std::vector<volume_type>::const_iterator   const_cell_iterator;

    cell_iterator           cells_begin() { return this->backend_storage()->volumes.begin(); }
    cell_iterator           cells_end()   { return this->backend_storage()->volumes.end(); }
    const_cell_iterator     cells_begin() const { return this->backend_storage()->volumes.begin(); }
    const_cell_iterator     cells_end()   const { return this->backend_storage()->volumes.end(); }

    /* face iterators */
    typedef typename std::vector<surface_type>::iterator        face_iterator;
    typedef typename std::vector<surface_type>::const_iterator  const_face_iterator;

    face_iterator           faces_begin() { return this->backend_storage()->surfaces.begin(); }
    face_iterator           faces_end()   { return this->backend_storage()->surfaces.end(); }
    const_face_iterator     faces_begin() const { return this->backend_storage()->surfaces.begin(); }
    const_face_iterator     faces_end()   const { return this->backend_storage()->surfaces.end(); }

    size_t  cells_size() const { return this->backend_storage()->volumes.size(); }
    size_t  faces_size() const { return this->backend_storage()->surfaces.size(); }

};



/* \brief Template specialization for 2D meshes.
 *
 * \tparam T       Type that the mesh uses to represent points.
 * \tparam Storage Type of storage to represent the data
 */
template<typename T, typename Storage>
class mesh_base<T,2,Storage> : public mesh_bones<T,2,Storage>
{
public:
    typedef typename Storage::surface_type                      surface_type;
    typedef typename Storage::edge_type                         edge_type;
    typedef typename Storage::node_type                         node_type;
    typedef typename Storage::point_type                        point_type;

    typedef surface_type                                        cell;
    typedef edge_type                                           face;
    const static size_t dimension = 2;

    /* cell iterators */
    typedef typename std::vector<surface_type>::iterator        cell_iterator;
    typedef typename std::vector<surface_type>::const_iterator  const_cell_iterator;

    cell_iterator           cells_begin() { return this->backend_storage()->surfaces.begin(); }
    cell_iterator           cells_end()   { return this->backend_storage()->surfaces.end(); }
    const_cell_iterator     cells_begin() const { return this->backend_storage()->surfaces.begin(); }
    const_cell_iterator     cells_end()   const { return this->backend_storage()->surfaces.end(); }

    /* face iterators */
    typedef typename std::vector<edge_type>::iterator           face_iterator;
    typedef typename std::vector<edge_type>::const_iterator     const_face_iterator;

    face_iterator           faces_begin() { return this->backend_storage()->edges.begin(); }
    face_iterator           faces_end()   { return this->backend_storage()->edges.end(); }
    const_face_iterator     faces_begin() const { return this->backend_storage()->edges.begin(); }
    const_face_iterator     faces_end()   const { return this->backend_storage()->edges.end(); }

    size_t  cells_size() const { return this->backend_storage()->surfaces.size(); }
    size_t  faces_size() const { return this->backend_storage()->edges.size(); }
};

/* \brief Template specialization for 1D meshes.
 *
 * \tparam T       Type that the mesh uses to represent points.
 * \tparam Storage Type of storage to represent the data
 */
template<typename T, typename Storage>
class mesh_base<T,1,Storage> : public mesh_bones<T,1,Storage>
{
public:
    typedef typename Storage::edge_type                     edge_type;
    typedef typename Storage::node_type                     node_type;
    typedef typename Storage::point_type                    point_type;

    typedef edge_type                                       cell;
    typedef node_type                                       face;
    const static size_t dimension = 1;

    /* cell iterators */
    typedef typename std::vector<edge_type>::iterator       cell_iterator;
    typedef typename std::vector<edge_type>::const_iterator const_cell_iterator;

    cell_iterator           cells_begin() { return this->backend_storage()->edges.begin(); }
    cell_iterator           cells_end()   { return this->backend_storage()->edges.end(); }
    const_cell_iterator     cells_begin() const { return this->backend_storage()->edges.begin(); }
    const_cell_iterator     cells_end()   const { return this->backend_storage()->edges.end(); }

    /* face iterators */
    typedef typename std::vector<node_type>::iterator       face_iterator;
    typedef typename std::vector<node_type>::const_iterator const_face_iterator;

    face_iterator           faces_begin() { return this->backend_storage()->nodes.begin(); }
    face_iterator           faces_end()   { return this->backend_storage()->nodes.end(); }
    const_face_iterator     faces_begin() const { return this->backend_storage()->nodes.begin(); }
    const_face_iterator     faces_end()   const { return this->backend_storage()->nodes.end(); }

    size_t  cells_size() const { return this->backend_storage()->edges.size(); }
    size_t  faces_size() const { return this->backend_storage()->nodes.size(); }

};

} // namespace priv

template<typename T, size_t DIM, typename Storage>
class mesh : public priv::mesh_base<T,DIM,Storage>
{
    typedef priv::mesh_base<T,DIM,Storage> base_type;

public:
    static const size_t dimension = DIM;

    typedef T           coordinate_type;
    typedef Storage     storage_type;

    typedef typename priv::mesh_base<T, DIM, Storage>::point_type   point_type;

    typedef typename priv::mesh_base<T, DIM, Storage>::cell         cell;
    typedef typename priv::mesh_base<T, DIM, Storage>::face         face;
    typedef typename priv::mesh_base<T, DIM, Storage>::cell         cell_type;
    typedef typename priv::mesh_base<T, DIM, Storage>::face         face_type;

    /* point iterators */
    typedef typename std::vector<point_type>::iterator              point_iterator;
    typedef typename std::vector<point_type>::const_iterator        const_point_iterator;

    point_iterator          points_begin() { return this->backend_storage()->points.begin(); }
    point_iterator          points_end()   { return this->backend_storage()->points.end(); }
    const_point_iterator    points_begin() const { return this->backend_storage()->points.begin(); }
    const_point_iterator    points_end()   const { return this->backend_storage()->points.end(); }

    size_t  points_size() const { return this->backend_storage()->points.size(); }

    typedef priv::filter_iterator<typename mesh::face_iterator,
                                  priv::is_boundary_pred<mesh>>
                                  boundary_face_iterator;

    typedef priv::filter_iterator<typename mesh::const_face_iterator,
                                  priv::is_boundary_pred<mesh>>
                                  const_boundary_face_iterator;

    typedef priv::filter_iterator<typename mesh::face_iterator,
                                  priv::is_boundary_pred_with_id<mesh>>
                                  boundary_face_with_id_iterator;

    typedef priv::filter_iterator<typename mesh::face_iterator,
                                  priv::is_internal_pred<mesh>>
                                  internal_face_iterator;

    /* Return a vector with the ids of all the boundaries */
    std::vector<size_t>
    boundary_id_list(void) const
    {
        assert( this->backend_storage()->boundary_info.size() == this->faces_size() );
        std::set<size_t> boundary_ids;
        for (auto& bi : this->backend_storage()->boundary_info)
            if (bi.is_boundary())
                boundary_ids.insert(bi.id());

        return std::vector<size_t>(boundary_ids.begin(), boundary_ids.end());
    }

    /* Determine if the specified face is a boundary (via face id) */
    bool is_boundary(typename face::id_type id) const
    {
        assert( this->backend_storage()->boundary_info.size() == this->faces_size() );
        return this->backend_storage()->boundary_info.at(id).is_boundary();
    }

    /* Determine if the specified face is a boundary (via actual face) */
    bool is_boundary(const face& f) const
    {
        assert( this->backend_storage()->boundary_info.size() == this->faces_size() );
        auto fid = lookup(f);
        return this->backend_storage()->boundary_info.at(fid).is_boundary();
    }

    /* Determine if the specified face is a boundary (via iterator) */
    bool is_boundary(const typename base_type::face_iterator& itor) const
    {
        assert( this->backend_storage()->boundary_info.size() == this->faces_size() );
        auto ofs = std::distance(this->faces_begin(), itor);
        return this->backend_storage()->boundary_info.at(ofs).is_boundary();
    }

    size_t boundary_id(const typename face::id_type& fid) const
    {
        assert( this->backend_storage()->boundary_info.size() == this->faces_size() );
        if ( !is_boundary(fid) )
            throw std::invalid_argument("Face is not a boundary.");

        return this->backend_storage()->boundary_info.at(fid).id();
    }

    size_t boundary_id(const face& f) const
    {
        assert( this->backend_storage()->boundary_info.size() == this->faces_size() );
        auto fid = lookup(f);

        if ( !is_boundary(fid) )
            throw std::invalid_argument("Face is not a boundary.");

        return this->backend_storage()->boundary_info.at(fid).id();
    }

    /* Get the total number of boundary faces */
    size_t  boundary_faces_size() const
    {
        assert( this->backend_storage()->boundary_info.size() == this->faces_size() );

        auto count_lambda = [](const boundary_descriptor& bi) -> bool {
            return bi.is_boundary();
        };

        return std::count_if(this->backend_storage()->boundary_info.begin(),
                             this->backend_storage()->boundary_info.end(),
                             count_lambda);
    }

    /* Get the number of boundary faces in the boundary with specified id */
    size_t  boundary_faces_size(size_t id) const
    {
        assert( this->backend_storage()->boundary_info.size() == this->faces_size() );

        auto count_lambda = [&](const boundary_descriptor& bi) -> bool {
            return bi.is_boundary() and bi.id() == id;
        };

        return std::count_if(this->backend_storage()->boundary_info.begin(),
                             this->backend_storage()->boundary_info.end(),
                             count_lambda);
    }

    /* Get the total number of internal faces */
    size_t  internal_faces_size() const
    {
        assert( this->backend_storage()->boundary_info.size() == this->faces_size() );

        auto count_lambda = [](const boundary_descriptor& bi) -> bool {
            return not bi.is_boundary();
        };

        return std::count_if(this->backend_storage()->boundary_info.begin(),
                             this->backend_storage()->boundary_info.end(),
                             count_lambda);
    }

    /* Begin iterator to all the boundary faces */
    boundary_face_iterator boundary_faces_begin()
    {
        assert( this->backend_storage()->boundary_info.size() == this->faces_size() );
        typedef priv::is_boundary_pred<mesh> ibp;
        return boundary_face_iterator(ibp(*this), this->faces_begin(), this->faces_end());
    }

    const_boundary_face_iterator boundary_faces_begin() const
    {
        assert( this->backend_storage()->boundary_info.size() == this->faces_size() );
        typedef priv::is_boundary_pred<mesh> ibp;
        return const_boundary_face_iterator(ibp(*this), this->faces_begin(), this->faces_end());
    }

    /* End iterator to all the boundary faces */
    boundary_face_iterator boundary_faces_end()
    {
        assert( this->backend_storage()->boundary_info.size() == this->faces_size() );
        typedef priv::is_boundary_pred<mesh> ibp;
        return boundary_face_iterator(ibp(*this), this->faces_end(), this->faces_end());
    }

    const_boundary_face_iterator boundary_faces_end() const
    {
        assert( this->backend_storage()->boundary_info.size() == this->faces_size() );
        typedef priv::is_boundary_pred<mesh> ibp;
        return const_boundary_face_iterator(ibp(*this), this->faces_end(), this->faces_end());
    }

    /* Begin iterator to the boundary faces of boundary 'id' */
    boundary_face_with_id_iterator boundary_faces_begin(size_t id)
    {
        assert( this->backend_storage()->boundary_info.size() == this->faces_size() );
        typedef priv::is_boundary_pred_with_id<mesh> ibp;
        return boundary_face_with_id_iterator(ibp(*this, id), this->faces_begin(), this->faces_end());
    }

    /* End iterator to the boundary faces of boundary 'id' */
    boundary_face_with_id_iterator  boundary_faces_end(size_t id)
    {
        assert( this->backend_storage()->boundary_info.size() == this->faces_size() );
        typedef priv::is_boundary_pred_with_id<mesh> ibp;
        return boundary_face_with_id_iterator(ibp(*this, id), this->faces_end(), this->faces_end());
    }

    /* Begin iterator to the internal faces */
    internal_face_iterator  internal_faces_begin()
    {
        assert( this->backend_storage()->boundary_info.size() == this->faces_size() );
        typedef priv::is_internal_pred<mesh> iip;
        return internal_face_iterator(iip(*this), this->faces_begin(), this->faces_end());
    }

    /* End iterator to the internal faces */
    internal_face_iterator  internal_faces_end()
    {
        assert( this->backend_storage()->boundary_info.size() == this->faces_size() );
        typedef priv::is_internal_pred<mesh> iip;
        return internal_face_iterator(iip(*this), this->faces_end(), this->faces_end());
    }

    subdomain_descriptor domain_info(const cell& c) const
    {
        assert( this->backend_storage()->subdomain_info.size() == this->cells_size() );
        auto cid = lookup(c);
        return this->backend_storage()->subdomain_info.at(cid);
    }

    boundary_descriptor boundary_info(const face& f) const
    {
        assert( this->backend_storage()->boundary_info.size() == this->faces_size() );
        auto fid = lookup(f);
        return this->backend_storage()->boundary_info.at(fid);
    }

    /* Apply a transformation to the mesh. Transform should be a functor or
     * a lambda function of type
     *      mesh_type::point_type -> mesh_type::point_type
     */
    template<typename Transform>
    void transform(const Transform& tr)
    {
        std::transform(points_begin(), points_end(), points_begin(), tr);
    }

    /* Returns the numerial ID of a cell. */
    typename cell::id_type lookup(const cell& cl) const
    {
        auto ci = find_element_id(this->cells_begin(), this->cells_end(), cl);
        if (!ci.first)
            throw std::invalid_argument("Cell not present in mesh");

        return ci.second;
    }

    /* Returns the numerial ID of a face. */
    typename face::id_type lookup(const face& fc) const
    {
        auto fi = find_element_id(this->faces_begin(), this->faces_end(), fc);
        if (!fi.first)
        {
            std::stringstream ss;
            ss << fc << ": Face not present in mesh";
            throw std::invalid_argument(ss.str());
        }
        return fi.second;
    }

    auto operator[](size_t cell_num) const {
        assert(cell_num < this->cells_size());
        return *std::next(this->cells_begin(), cell_num);
    }

    void statistics(void) const
    {
        this->backend_storage()->statistics();
    }
};

template<typename Mesh>
concept mesh_1D = requires {
    typename Mesh::cell_type;
    typename Mesh::face_type;
    Mesh::dimension;
    requires Mesh::dimension == 1;
};

template<typename Mesh>
concept mesh_2D = requires {
    typename Mesh::cell_type;
    typename Mesh::face_type;
    Mesh::dimension;
    requires Mesh::dimension == 2;
};

template<typename Mesh>
concept mesh_3D = requires {
    typename Mesh::cell_type;
    typename Mesh::face_type;
    Mesh::dimension;
    requires Mesh::dimension == 3;
};

template<typename Mesh>
void mark_internal_faces(Mesh& msh)
{
    std::vector<size_t> neigh_count(msh.faces_size(), 0);

    for (auto& cl : cells(msh))
    {
        auto fcs = faces(msh, cl);
        for (auto& fc : fcs)
        {
            auto ofs = offset(msh, fc);
            neigh_count[ofs]++;
        }
    }

    auto storage = msh.backend_storage();
    for (size_t i = 0; i < msh.faces_size(); i++)
    {
        auto nc = neigh_count[i];
        assert( (nc == 0) or (nc == 1) or (nc == 2) );
        assert(storage->boundary_info.size() == msh.faces_size());
        auto bi = storage->boundary_info.at(i);
        if ( bi.is_boundary() and (nc == 2) )
            bi.is_internal(true);
        storage->boundary_info.at(i) = bi;
    }
}

template<typename Mesh>
class neighbour_connectivity
{
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using cit = typename cell_type::id_type;
    using cell_opt = std::optional<cit>;
    using face_owners_type = std::array<cell_opt, 2>;
    
    std::vector<face_owners_type> face_owners;

    void compute_connectivity(const Mesh& msh)
    {
        face_owners.resize( msh.faces_size() );

        size_t cell_i = 0;
        for (auto& cl : msh)
        {
            auto fcs = faces(msh, cl);

            for (auto& fc : fcs)
            {
                auto fi = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
                if (!fi.first)
                {
                    std::stringstream ss;
                    ss << fc << ": Face not present in mesh";
                    throw std::invalid_argument(ss.str());
                }
                auto fc_id = offset(msh, fc);

                auto& fo = face_owners.at(fc_id);
                if (not fo[0])
                    fo[0] = cit(cell_i);
                else if (not fo[1])
                    fo[1] = cit(cell_i);
                else
                    throw std::logic_error("BUG: a face has max 2 owners");
            }

            ++cell_i;
        }
    }

public:
    neighbour_connectivity(const Mesh& msh)
    {
        compute_connectivity(msh);
    }

    std::pair<cell_type, bool>
    neighbour_via(const Mesh& msh, const cell_type& cl, const face_type& fc)
    {
        if ( face_owners.size() != msh.faces_size() )
            throw std::logic_error("Inconsistent neighbour information.");

        auto cl_ofs = offset(msh, cl);
        auto fc_ofs = offset(msh, fc);

        auto fo = face_owners.at( fc_ofs );

        if ( fo[0] and (fo[0] != cl_ofs) )
            std::swap(fo[0], fo[1]);

        assert(fo[0].value() == cl_ofs);

        if (not fo[1])
            return std::make_pair(*msh.cells_begin(), false);

        return std::make_pair( *std::next(msh.cells_begin(), fo[1].value()), true);
    }
};

template<typename Mesh>
auto connectivity_via_faces(const Mesh& msh)
{
    return neighbour_connectivity<Mesh>(msh);
}

template<typename T, size_t DIM, typename Storage>
typename mesh<T, DIM, Storage>::cell_iterator
begin(mesh<T, DIM, Storage>& msh)
{
    return msh.cells_begin();
}

template<typename T, size_t DIM, typename Storage>
typename mesh<T, DIM, Storage>::cell_iterator
end(mesh<T, DIM, Storage>& msh)
{
    return msh.cells_end();
}

template<typename T, size_t DIM, typename Storage>
typename mesh<T, DIM, Storage>::const_cell_iterator
begin(const mesh<T, DIM, Storage>& msh)
{
    return msh.cells_begin();
}

template<typename T, size_t DIM, typename Storage>
typename mesh<T, DIM, Storage>::const_cell_iterator
end(const mesh<T, DIM, Storage>& msh)
{
    return msh.cells_end();
}

namespace priv {

template<typename Mesh>
class cells_iterproxy
{
    Mesh& m_msh;

public:
    cells_iterproxy(Mesh& msh)
        : m_msh(msh)
    {}

    auto begin()        { return m_msh.cells_begin(); }
    auto begin() const  { return m_msh.cells_begin(); }
    auto end()          { return m_msh.cells_end(); }
    auto end() const    { return m_msh.cells_end(); }
};

template<typename Mesh>
class faces_iterproxy
{
    const Mesh& m_msh;

public:
    faces_iterproxy(const Mesh& msh)
        : m_msh(msh)
    {}

    auto begin()        { return m_msh.faces_begin(); }
    auto begin() const  { return m_msh.faces_begin(); }
    auto end()          { return m_msh.faces_end(); }
    auto end() const    { return m_msh.faces_end(); }
};


template<typename Mesh>
class boundary_faces_iterproxy
{
    Mesh& m_msh;

public:
    boundary_faces_iterproxy(Mesh& msh)
        : m_msh(msh)
    {}
    auto begin()        { return m_msh.boundary_faces_begin(); }
    auto begin() const  { return m_msh.boundary_faces_begin(); }
    auto end()          { return m_msh.boundary_faces_end(); }
    auto end() const    { return m_msh.boundary_faces_end(); }
};

} // namespace priv

/**
 * Get iterator proxy to iterate on cells with
 *
 *      for (auto& cl : cells(msh)) {}
 */
template<typename T, size_t DIM, typename Storage>
auto cells(mesh<T, DIM, Storage>& msh)
{
    using Mesh = mesh<T, DIM, Storage>;
    return priv::cells_iterproxy<Mesh>(msh);
}

/**
 * Get iterator proxy to iterate on faces with
 *
 *      for (auto& fc : faces(msh)) {}
 */

template<typename T, size_t DIM, typename Storage>
auto faces(const mesh<T, DIM, Storage>& msh)
{
    using Mesh = mesh<T, DIM, Storage>;
    return priv::faces_iterproxy<Mesh>(msh);
}


/**
 * Get iterator proxy to iterate on boundry faces with
 *
 *      for (auto& fc : boundary_faces(msh)) {}
*/

template<typename T, size_t DIM, typename Storage>
auto boundary_faces(mesh<T, DIM, Storage>& msh)
{
    using Mesh = mesh<T, DIM, Storage>;
    return priv::boundary_faces_iterproxy<Mesh>(msh);
}

/*
template <typename T, size_t DIM, typename Storage>
auto boundary_faces(mesh<T, DIM, Storage>& msh, size_t boundary_id)
{
    using Mesh = mesh<T, DIM, Storage>;
    return priv::boundary_faces_iterproxy<Mesh>(msh);
}
*/


template<typename Mesh>
void cell_info(const Mesh& msh, const typename Mesh::cell& cl)
{
    std::cout << "** CELL INFORMATION BEGIN **" << std::endl;
    std::cout << cl << std::endl;
    auto fcs = faces(msh, cl);
    for (auto& fc : fcs)
    {
        std::cout << "  - " << fc;
        if ( msh.is_boundary(fc) )
            std::cout << " [B " << msh.boundary_id(fc) << "]";

        std::cout << std::endl;
    }
    std::cout << "** CELL INFORMATION END **" << std::endl;
}

template<typename Mesh>
size_t
offset(const Mesh& msh, const typename Mesh::cell_type& cl)
{
    auto itor = std::lower_bound(msh.cells_begin(), msh.cells_end(), cl);
    if ( itor == msh.cells_end() )
        throw std::logic_error("Cell not found: this is likely a bug.");

    return std::distance(msh.cells_begin(), itor);
}

template<typename Mesh>
size_t
offset(const Mesh& msh, const typename Mesh::face_type& fc)
{
    auto itor = std::lower_bound(msh.faces_begin(), msh.faces_end(), fc);
    if ( itor == msh.faces_end() )
        throw std::logic_error("Face not found: this is likely a bug.");

    return std::distance(msh.faces_begin(), itor);
}

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
void
dump_to_matlab(const Mesh<T, 2, Storage>& msh, const std::string& filename)
{
    std::ofstream ofs(filename);

    size_t elemnum = 0;
    for (auto cl : msh)
    {
        auto bar = barycenter(msh, cl);

        ofs << "text(" << bar.x() << "," << bar.y() << ",'" << elemnum << "');" << std::endl;

        auto fcs = faces(msh, cl);
        for (auto fc : fcs)
        {
            auto ptids = fc.point_ids();
            auto pts = points(msh, fc);
            assert(ptids.size() == pts.size());

            for (size_t i = 0; i < ptids.size(); i++)
            {
                ofs << "text(" << pts[i].x() << "," << pts[i].y() << ",'" << ptids[i] << "');" << std::endl;
            }

            if ( msh.is_boundary(fc) )
            {
                ofs << "line([" << pts[0].x() << " " << pts[1].x() << "], [";
                ofs << pts[0].y() << " " << pts[1].y() << "], 'Color', 'r');";
                ofs << std::endl;
            }
            else
            {
                ofs << "line([" << pts[0].x() << " " << pts[1].x() << "], [";
                ofs << pts[0].y() << " " << pts[1].y() << "], 'Color', 'g');";
                ofs << std::endl;
            }
        }
        elemnum++;
    }

    ofs.close();
}

template<typename Mesh>
class bounding_box
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::point_type      point_type;

    point_type                                  m_bb_min, m_bb_max;

public:
    bounding_box()
    {}

    bounding_box(const Mesh& msh)
    {
        compute(msh);
    }

    void compute(const Mesh& msh)
    {
        m_bb_min = *msh.points_begin();
        m_bb_max = *msh.points_begin();

        for (auto itor = std::next(msh.points_begin());
                  itor != msh.points_end();
                  itor++)
        {
            bool lt = false, gt = false;
            auto curpt = *itor;

            for (size_t i = 0; i < point_type::dimension; i++)
            {
                lt = lt or (curpt[i] < m_bb_min[i]);
                gt = gt or (curpt[i] > m_bb_max[i]);
            }
            assert(lt != gt); //lt xor gt
            if (lt) m_bb_min = curpt;
            if (gt) m_bb_max = curpt;
        }
    }

    point_type min() const
    {
        return m_bb_min;
    }

    point_type max() const
    {
        return m_bb_max;
    }
};

template<typename Mesh>
std::ostream& operator<<(std::ostream& os, const bounding_box<Mesh>& bb)
{
    os << "[" << bb.min() << ", " << bb.max() << "]";
    return os;
}

template<typename MeshType>
void
mesh_aspect(const MeshType& msh, const std::string& filename)
{
    using T = typename MeshType::coordinate_type;

    std::ofstream myfile;
    myfile.open(filename);
    if (myfile.is_open())
    {
        myfile << "# cell_id  meas_T  hT  hF_min  hF_max  hT / hF_min  hT / hF_max  hF_max / hF_min \n";
        for (auto& cl : msh)
        {
            const auto cell_id = msh.lookup(cl);
            const auto hT      = diameter(msh, cl);
            const auto mT      = measure(msh, cl);
            T          hF_min  = hT;
            T          hF_max  = 0;

            const auto fcs = faces(msh, cl);
            for (auto& fc : fcs)
            {
                const auto hF = diameter(msh, fc);
                hF_min        = std::min(hF_min, hF);
                hF_max        = std::max(hF_max, hF);
            }

            myfile << cell_id << " " << mT << " " << hT << " " << hF_min << " " << hF_max << " " << hT / hF_min << " "
                   << hT / hF_max << " " << hF_max / hF_min << "\n";
        }
    }

    myfile.close();
}

} // namespace disk
