/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
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

namespace disk
{

class DegreeInfo
{
  private:
    bool   m_has_unknowns;
    size_t m_degree;

  public:
    DegreeInfo(void) : m_has_unknowns(false), m_degree(0) {}
    DegreeInfo(const size_t degree) : m_has_unknowns(true), m_degree(degree) {}
    DegreeInfo(const bool unknowns) : m_has_unknowns(unknowns), m_degree(0) {}
    DegreeInfo(const bool unknowns, const size_t degree) : m_has_unknowns(unknowns), m_degree(degree) {}

    size_t
    degree(void) const
    {
        return m_degree;
    }

    void
    degree(const size_t degree)
    {
        m_degree = degree;
    }

    bool
    hasUnknowns(void) const
    {
        return m_has_unknowns;
    }

    void
    hasUnknowns(const bool unknowns)
    {
        m_has_unknowns = unknowns;
    }
};

template<typename MeshType>
class MeshDegree
{
  private:
    typedef typename MeshType::cell            cell_type;
    typedef typename MeshType::face            face_type;
    typedef typename MeshType::coordinate_type scalar_type;

    std::vector<DegreeInfo> m_cells_degree, m_faces_degree;

  public:
    MeshDegree(void)
    {
        m_cells_degree.clear();
        m_faces_degree.clear();
    }
    MeshDegree(const MeshType& msh)
    {
        m_cells_degree.resize(msh.cells_size());
        m_faces_degree.resize(msh.faces_size());
    }

    MeshDegree(const MeshType& msh, const size_t degree)
    {
        const DegreeInfo di(degree);
        m_cells_degree.assign(msh.cells_size(), di);
        m_faces_degree.assign(msh.faces_size(), di);
    }

    MeshDegree(const MeshType& msh, const size_t cell_degree, const size_t face_degree)
    {
        const DegreeInfo dic(cell_degree);
        const DegreeInfo dif(face_degree);

        m_cells_degree.assign(msh.cells_size(), dic);
        m_faces_degree.assign(msh.faces_size(), dif);
    }

    void
    degree(const MeshType& msh, const cell_type& cl, const size_t degree)
    {
        const auto cell_id = msh.lookup(cl);
        m_cells_degree[cell_id].degree(degree);
    }

    size_t
    degree(const MeshType& msh, const cell_type& cl) const
    {
        const auto cell_id = msh.lookup(cl);
        return m_cells_degree[cell_id].degree();
    }

    void
    degree(const MeshType& msh, const face_type& fc, const size_t degree)
    {
        const auto face_id = msh.lookup(fc);
        m_faces_degree[face_id].degree(degree);
    }

    size_t
    degree(const MeshType& msh, const face_type& fc) const
    {
        const auto face_id = msh.lookup(fc);
        return m_faces_degree[face_id].degree();
    }

    void
    hasUnknowns(const MeshType& msh, const cell_type& cl, const bool unknowns)
    {
        const auto cell_id = msh.lookup(cl);
        m_cells_degree[cell_id].hasUnknowns(unknowns);
    }

    bool
    hasUnknowns(const MeshType& msh, const cell_type& cl) const
    {
        const auto cell_id = msh.lookup(cl);
        return m_cells_degree[cell_id].hasUnknowns();
    }

    void
    hasUnknowns(const MeshType& msh, const face_type& fc, const bool unknowns)
    {
        const auto face_id = msh.lookup(fc);
        m_faces_degree[face_id].hasUnknowns(unknowns);
    }

    bool
    hasUnknowns(const MeshType& msh, const face_type& fc) const
    {
        const auto face_id = msh.lookup(fc);
        return m_faces_degree[face_id].hasUnknowns();
    }

    DegreeInfo
    degreeInfo(const MeshType& msh, const cell_type& cl) const
    {
        const auto cell_id = msh.lookup(cl);
        return m_cells_degree[cell_id];
    }

    void
    degreeInfo(const MeshType& msh, const cell_type& cl, const DegreeInfo& di)
    {
        const auto cell_id             = msh.lookup(cl);
        return m_cells_degree[cell_id] = di;
    }

    DegreeInfo
    degreeInfo(const MeshType& msh, const face_type& fc) const
    {
        const auto face_id = msh.lookup(fc);
        return m_faces_degree[face_id];
    }

    template<typename faces_type>
    std::vector<DegreeInfo>
    degreeInfo(const MeshType& msh, const faces_type& fcs) const
    {
        std::vector<DegreeInfo> ret(fcs.size());

        for (auto& fc : fcs)
        {
            const auto face_id = msh.lookup(fc);
            ret.push_back(m_faces_degree[face_id]);
        }

        return ret;
    }

    void
    degreeInfo(const MeshType& msh, const face_type& fc, const DegreeInfo& di)
    {
        const auto face_id             = msh.lookup(fc);
        return m_faces_degree[face_id] = di;
    }

    std::vector<DegreeInfo>
    faces_degree(void) const
    {
        return m_faces_degree;
    }

    std::vector<DegreeInfo>
    cells_degree(void) const
    {
        return m_cells_degree;
    }
};

} // end disk