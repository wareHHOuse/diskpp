/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
 * Karol Cascavita (C) 2018                     karol.cascavita@enpc.fr
 * Nicolas Pignet  (C) 2018, 2019               nicolas.pignet@enpc.fr
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

#include "adaptivity/adaptivity.hpp"
#include "bases/bases.hpp"

namespace disk
{

/**
 * @brief this class contains the different polynomial degrees which can used with HHO methods
 *
 */
class hho_degree_info
{
    /**
     * @brief
     *  cell_deg : polynomial degree used for the cells
     *  face_deg : polynomial degree used for the faces
     *  grad_deg : polynomial degree used for the gradient
     *
     */
    size_t cell_deg, face_deg, grad_deg;

  public:
    /**
     * @brief Construct a new hho degree info object
     * The default polynomial degree is 1 for the face and cell degrees
     *
     */
    hho_degree_info() : cell_deg(1), face_deg(1), grad_deg(1) {}

    /**
     * @brief Construct a new hho degree info object
     *
     * @param degree polynomial degree used for the cells, faces and gradient
     */
    explicit hho_degree_info(size_t degree) : cell_deg(degree), face_deg(degree), grad_deg(degree) {}

    /**
     * @brief Construct a new hho degree info object
     *
     * @param cd polynomial degree used for the cells
     * @param fd polynomial degree used for the faces
     *
     * @brief Note that
     * fd >= 0 and cd >=0
     * fd - 1 <= cd <= fd +1
     */
    hho_degree_info(size_t cd, size_t fd)
    {
        bool c1 = fd > 0 && (cd == fd - 1 || cd == fd || cd == fd + 1);
        bool c2 = fd == 0 && (cd == fd || cd == fd + 1);
        if (c1 || c2)
        {
            cell_deg = cd;
            face_deg = fd;
            grad_deg = fd;
        }
        else
        {
            std::cout << "Invalid cell degree. Reverting to equal-order" << std::endl;
            cell_deg = fd;
            face_deg = fd;
            grad_deg = fd;
        }
    }

    /**
     * @brief Construct a new hho degree info object
     *
     * @param cd polynomial degree used for the cells
     * @param fd polynomial degree used for the faces
     * @param gd polynomial degree used for the gradient
     *
     * @brief Note that
     * fd >= 0, cd >=0 and gd >=0
     * fd - 1 <= cd <= fd +1
     * gd >= fd
     */
    hho_degree_info(size_t cd, size_t fd, size_t gd)
    {
        bool c1 = fd > 0 && (cd == fd - 1 || cd == fd || cd == fd + 1);
        bool c2 = fd == 0 && (cd == fd || cd == fd + 1);
        bool c3 = gd >= fd;
        if (c1 || c2 || c3)
        {
            cell_deg = cd;
            face_deg = fd;
            grad_deg = gd;
        }
        else
        {
            std::cout << "Invalid cell degree. Reverting to equal-order" << std::endl;
            cell_deg = fd;
            face_deg = fd;
            grad_deg = fd;
        }
    }

    /**
     * @brief Return the polynomial degree used for the cells
     *
     * @return size_t polynomial degree used for the cells
     */
    size_t
    cell_degree() const
    {
        return cell_deg;
    }

    /**
     * @brief Return the polynomial degree used for the faces
     *
     * @return size_t polynomial degree used for the faces
     */
    size_t
    face_degree() const
    {
        return face_deg;
    }

    /**
     * @brief Return the polynomial degree used for the reconstruction operator
     *
     * @return size_t polynomial degree used for the reconstruction operator
     */
    size_t
    reconstruction_degree() const
    {
        return face_deg + 1;
    }

    /**
     * @brief Return the polynomial degree used for the gradient reconstruction
     *
     * @return size_t polynomial degree used for the gradient reconstruction
     */
    size_t
    grad_degree() const
    {
        return grad_deg;
    }

    /**
     * @brief Print the different degree
     *
     */
    void
    info_degree() const
    {
        std::cout << "cell degree: " << cell_deg << std::endl;
        std::cout << "face degree: " << face_deg << std::endl;
        std::cout << "gradient degree: " << grad_deg << std::endl;
        std::cout << "reconstruction degree: " << face_deg + 1 << std::endl;
    }
};

// const MeshDegreeInfo<Mesh>& degree_infos
template<typename Mesh>
size_t
scalar_face_dofs(const Mesh& msh, const DegreeInfo& di)
{
    if (di.hasUnknowns())
    {
        return scalar_basis_size(di.degree(), Mesh::dimension - 1);
    }

    return 0;
}

template<typename Mesh>
size_t
vector_face_dofs(const Mesh& msh, const DegreeInfo& di)
{
    if (di.hasUnknowns())
    {
        return vector_basis_size(di.degree(), Mesh::dimension - 1, Mesh::dimension);
    }

    return 0;
}

template<typename Mesh>
size_t
scalar_faces_dofs(const Mesh& msh, const typename Mesh::cell_type& cl, const MeshDegreeInfo<Mesh>& degree_infos)
{
    const auto fcs      = faces(msh, cl);
    size_t     num_dofs = 0;

    for (auto fc : fcs)
    {
        num_dofs += scalar_face_dofs(msh, degree_infos.degreeInfo(msh, fc));
    }

    return num_dofs;
}

template<typename Mesh>
size_t
scalar_faces_dofs(const Mesh& msh, const std::vector<DegreeInfo>& faces_degree_infos)
{
    size_t num_dofs = 0;

    for (auto& fdi : faces_degree_infos)
    {
        num_dofs += scalar_face_dofs(msh, fdi);
    }

    return num_dofs;
}

template<typename Mesh>
size_t
vector_faces_dofs(const Mesh& msh, const typename Mesh::cell_type& cl, const MeshDegreeInfo<Mesh>& degree_infos)
{
    const auto fcs      = faces(msh, cl);
    size_t     num_dofs = 0;

    for (auto fc : fcs)
    {
        num_dofs += vector_face_dofs(msh, degree_infos.degreeInfo(msh, fc));
    }

    return num_dofs;
}

template<typename Mesh>
size_t
vector_faces_dofs(const Mesh& msh, const std::vector<DegreeInfo>& faces_degree_infos)
{
    size_t num_dofs = 0;

    for (auto& fdi : faces_degree_infos)
    {
        num_dofs += vector_face_dofs(msh, fdi);
    }

    return num_dofs;
}

template<typename Mesh>
std::vector<size_t>
scalar_faces_offset(const Mesh& msh, const typename Mesh::cell_type& cl, const MeshDegreeInfo<Mesh>& degree_infos)
{
    const auto fcs      = faces(msh, cl);
    size_t     num_dofs = 0;

    std::vector<size_t> ret;
    ret.reserve(fcs.size());

    for (auto fc : fcs)
    {
        ret.push_back(num_dofs);
        num_dofs += scalar_face_dofs(msh, fc);
    }

    return ret;
}

template<typename Mesh>
std::vector<size_t>
vector_faces_offset(const Mesh& msh, const typename Mesh::cell_type& cl, const MeshDegreeInfo<Mesh>& degree_infos)
{
    const auto fcs      = faces(msh, cl);
    size_t     num_dofs = 0;

    std::vector<size_t> ret;
    ret.reserve(fcs.size());

    for (auto fc : fcs)
    {
        ret.push_back(num_dofs);
        num_dofs += vector_face_dofs(msh, fc);
    }

    return ret;
}

} // end diskpp