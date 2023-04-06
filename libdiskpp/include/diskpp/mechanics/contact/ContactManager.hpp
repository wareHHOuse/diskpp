/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018, 2019
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2020                     nicolas.pignet@enpc.fr
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

#include <map>

#include "diskpp/common/eigen.hpp"

#include "diskpp/adaptivity/adaptivity.hpp"
#include "diskpp/boundary_conditions/boundary_conditions.hpp"

#include "diskpp/bases/bases.hpp"

namespace disk
{

namespace mechanics
{

template<typename MeshType>
class ContactManager
{
  private:
    typedef MeshType                            mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;

    typedef typename mesh_type::cell cell_type;
    typedef typename mesh_type::face face_type;

    typedef vector_boundary_conditions<mesh_type> bnd_type;

    std::map<size_t, size_t> mapOfFaces;
    std::map<size_t, size_t> mapOfMult;
    std::vector<size_t>      degree_mult;

  public:
    ContactManager(void)
    {
        mapOfFaces.clear();
        mapOfMult.clear();
        degree_mult.clear();
    }

    ContactManager(const mesh_type& msh, const bnd_type& bnd)
    {
        const auto num_contact_faces = bnd.nb_faces_contact();
        mapOfFaces.clear();
        mapOfMult.clear();
        degree_mult.clear();

        degree_mult.resize(num_contact_faces);
    }

    void
    addMapping(const size_t& face_id, const size_t& mult_id)
    {
        mapOfFaces.insert(std::make_pair(face_id, mult_id));
        mapOfMult.insert(std::make_pair(mult_id, face_id));
    }

    size_t
    getMappingFaceToMult(const size_t& face_id) const
    {
        return mapOfFaces.at(face_id);
    }

    size_t
    getMappingMultToFace(const size_t& mult_id) const
    {
        return mapOfMult.at(mult_id);
    }

    size_t
    numberOfMult(const mesh_type&                 msh,
                 const cell_type&                 cl,
                 const bnd_type&                  bnd,
                 const CellDegreeInfo<mesh_type>& cell_infos) const
    {

        return bnd.howmany_contact_faces(cl) *
               vector_basis_size(cell_infos.grad_degree(), mesh_type::dimension-1, mesh_type::dimension);
    }

    size_t
    numberOfMultFace(const size_t& mult_id) const
    {

        const auto degree = getDegreeMultFace(mult_id);
        return vector_basis_size(degree, mesh_type::dimension-1, mesh_type::dimension);
    }

    void
    setDegreeMultFace(const size_t& mult_id, const size_t& mult_degree)
    {
        degree_mult.at(mult_id) = mult_degree;
    }

    size_t
    getDegreeMultFace(const size_t& mult_id) const
    {
        return degree_mult.at(mult_id);
    }
};
}
}