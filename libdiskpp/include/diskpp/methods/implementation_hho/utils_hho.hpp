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

#include <algorithm>
#include <vector>

#include "diskpp/adaptivity/adaptivity.hpp"
#include "diskpp/bases/bases.hpp"

namespace disk
{

namespace priv
{
    struct hdi_named_args
    {
        size_t rd;
        size_t cd;
        size_t fd;
    };
}

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
    size_t cell_deg, face_deg, grad_deg, rec_deg;

  public:
    /**
     * @brief Construct a new hho degree info object
     * The default polynomial degree is 1 for the face and cell degrees
     *
     */
    hho_degree_info()
        : cell_deg(1), face_deg(1), grad_deg(1), rec_deg(2)
    {}

    /**
     * @brief Construct a new hho degree info object
     *
     * @param degree polynomial degree used for the cells, faces and gradient
     */
    explicit hho_degree_info(size_t degree)
        : cell_deg(degree), face_deg(degree), grad_deg(degree), rec_deg(degree+1)
    {}

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
            rec_deg  = fd+1;
        }
        else
        {
            std::cout << "Invalid cell degree. Reverting to equal-order" << std::endl;
            cell_deg = fd;
            face_deg = fd;
            grad_deg = fd;
            rec_deg  = fd+1;
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
            rec_deg = fd+1;
        }
        else
        {
            std::cout << "Invalid cell degree. Reverting to equal-order" << std::endl;
            cell_deg = fd;
            face_deg = fd;
            grad_deg = fd;
            rec_deg = fd+1;
        }
    }

    hho_degree_info(const priv::hdi_named_args& hna)
    {
        rec_deg = hna.rd;
        cell_deg = hna.cd;
        face_deg = hna.fd;
        grad_deg = hna.fd;
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

    void cell_degree(size_t k)
    {
        cell_deg = k;
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

    void face_degree(size_t k)
    {
        face_deg = k;
    }

    /**
     * @brief Return the polynomial degree used for the reconstruction operator
     *
     * @return size_t polynomial degree used for the reconstruction operator
     */
    size_t
    reconstruction_degree() const
    {
        return rec_deg;
    }

    void reconstruction_degree(size_t k)
    {
        rec_deg = k;
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
        std::cout << "reconstruction degree: " << rec_deg << std::endl;
    }

    bool operator<(const hho_degree_info& other) const
    {
        std::array<size_t, 4> degs_this{cell_deg, face_deg, grad_deg, rec_deg};
        std::array<size_t, 4> degs_other{other.cell_deg, other.face_deg, other.grad_deg, other.rec_deg};
        return std::lexicographical_compare(degs_this.begin(), degs_this.end(), degs_other.begin(), degs_other.end());
    }
};

class assembly_index
{
    size_t idx;
    bool   assem;

public:
    assembly_index(size_t i, bool as) : idx(i), assem(as) {}

    operator size_t() const {
        if (!assem)
            throw std::logic_error("Invalid assembly_index");

        return idx;
    }

    bool assemble() const {
        return assem;
    }

    friend std::ostream&
    operator<<(std::ostream& os, const assembly_index& as) {
        os << "(" << as.idx << "," << as.assem << ")";
        return os;
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

template<typename Mesh>
size_t
scalar_diff_dofs(const Mesh& msh, const CellDegreeInfo<Mesh>& cell_infos)
{
    const auto faces_infos = cell_infos.facesDegreeInfo();

    const auto dbs = scalar_basis_size(cell_infos.cell_degree(), Mesh::dimension-1);

    size_t     num_dofs = 0;

    for (auto fdi : faces_infos)
    {
        const auto fbs = scalar_face_dofs(msh, fdi);
        num_dofs += std::max(dbs, fbs);
    }

    return num_dofs;
}

// define some optimization
namespace priv
{

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
compute_lhs_vector(const Mesh&                                           msh,
                   const typename Mesh::cell&                            cl,
                   const CellDegreeInfo<Mesh>&                           cell_infos,
                   const dynamic_matrix<typename Mesh::coordinate_type>& lhs_scalar)
{
    typedef typename Mesh::coordinate_type scalar_type;

    const int  dimension      = Mesh::dimension;
    const auto num_cell_dofs  = vector_basis_size(cell_infos.cell_degree(), dimension, dimension);
    const auto num_faces_dofs = vector_faces_dofs(msh, cell_infos.facesDegreeInfo());

    const auto total_dofs = num_cell_dofs + num_faces_dofs;

    dynamic_matrix<scalar_type> lhs = dynamic_matrix<scalar_type>::Zero(total_dofs, total_dofs);

    const auto scal_total_dofs = total_dofs / dimension;

    assert(lhs_scalar.rows() == scal_total_dofs);
    assert(lhs_scalar.cols() == scal_total_dofs);

    int row, col;
#ifdef FILL_COLMAJOR
    for (int j = 0; j < scal_total_dofs; j++)
    {
        col = j * dimension;
        for (int i = 0; i < scal_total_dofs; i++)
        {
            row = i * dimension;
            for (int k = 0; k < dimension; k++)
            {
                lhs(row + k, col + k) = lhs_scalar(i, j);
            }
        }
    }
#else
    for (int i = 0; i < scal_total_dofs; i++)
    {
        row = i * dimension;
        for (int j = 0; j < scal_total_dofs; j++)
        {
            col = j * dimension;
            for (int k = 0; k < dimension; k++)
            {
                lhs(row + k, col + k) = lhs_scalar(i, j);
            }
        }
    }
#endif

    return lhs;
}

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
compute_lhs_vector(const Mesh&                                           msh,
                   const typename Mesh::cell&                            cl,
                   const hho_degree_info&                                hdi,
                   const dynamic_matrix<typename Mesh::coordinate_type>& lhs_scalar)
{
    const CellDegreeInfo<Mesh> cell_infos(msh, cl, hdi.cell_degree(), hdi.face_degree(), hdi.grad_degree());

    return compute_lhs_vector(msh, cl, cell_infos, lhs_scalar);
}

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
compute_grad_vector(const Mesh&                                           msh,
                    const typename Mesh::cell&                            cl,
                    const CellDegreeInfo<Mesh>&                           cell_infos,
                    const dynamic_matrix<typename Mesh::coordinate_type>& grad_scalar)
{
    typedef typename Mesh::coordinate_type scalar_type;

    const int  dimension     = Mesh::dimension;
    const auto rbs           = vector_basis_size(cell_infos.reconstruction_degree(), dimension, dimension);
    const auto num_cell_dofs = vector_basis_size(cell_infos.cell_degree(), dimension, dimension);

    const auto num_faces_dofs = vector_faces_dofs(msh, cell_infos.facesDegreeInfo());

    const auto total_dofs = num_cell_dofs + num_faces_dofs;

    dynamic_matrix<scalar_type> grad = dynamic_matrix<scalar_type>::Zero(rbs - dimension, total_dofs);

    const auto scal_rbs        = rbs / dimension;
    const auto scal_total_dofs = total_dofs / dimension;

    assert(grad_scalar.rows() == scal_rbs - 1);
    assert(grad_scalar.cols() == scal_total_dofs);

    int row, col;
#ifdef FILL_COLMAJOR
    for (int j = 0; j < scal_total_dofs; j++)
    {
        col = j * dimension;
        for (int i = 0; i < scal_total_dofs; i++)
        {
            row = i * dimension;
            for (int k = 0; k < dimension; k++)
            {
                grad(row + k, col + k) = grad_scalar(i, j);
            }
        }
    }
#else
    for (int i = 0; i < scal_rbs - 1; i++)
    {
        row = i * dimension;
        for (int j = 0; j < scal_total_dofs; j++)
        {
            col = j * dimension;
            for (int k = 0; k < dimension; k++)
            {
                grad(row + k, col + k) = grad_scalar(i, j);
            }
        }
    }
#endif

    return grad;
}

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
compute_grad_vector(const Mesh&                                           msh,
                    const typename Mesh::cell&                            cl,
                    const hho_degree_info&                                hdi,
                    const dynamic_matrix<typename Mesh::coordinate_type>& grad_scalar)
{
    const CellDegreeInfo<Mesh> cell_infos(msh, cl, hdi.cell_degree(), hdi.face_degree(), hdi.grad_degree());

    return compute_grad_vector(msh, cl, cell_infos, grad_scalar);
}

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
compute_grad_matrix(const Mesh&                                           msh,
                    const typename Mesh::cell&                            cl,
                    const CellDegreeInfo<Mesh>&                           cell_infos,
                    const dynamic_matrix<typename Mesh::coordinate_type>& grad_vector)
{
    typedef typename Mesh::coordinate_type scalar_type;

    const int  dimension      = Mesh::dimension;
    const auto gbs            = matrix_basis_size(cell_infos.grad_degree(), dimension, dimension);
    const auto num_cell_dofs  = vector_basis_size(cell_infos.cell_degree(), dimension, dimension);
    const auto num_faces_dofs = vector_faces_dofs(msh, cell_infos.facesDegreeInfo());

    const auto total_dofs = num_cell_dofs + num_faces_dofs;

    dynamic_matrix<scalar_type> grad = dynamic_matrix<scalar_type>::Zero(gbs, total_dofs);

    const auto vec_gbs         = gbs / dimension;
    const auto scal_total_dofs = total_dofs / dimension;

    assert(grad_vector.rows() == vec_gbs);
    assert(grad_vector.cols() == scal_total_dofs);

    int row, col;
#ifdef FILL_COLMAJOR
    for (int j = 0; j < scal_total_dofs; j++)
    {
        col = j * dimension;
        for (int i = 0; i < vec_gbs; i++)
        {
            row = i * dimension;
            for (int k = 0; k < dimension; k++)
            {
                grad(row + k, col + k) = grad_vector(i, j);
            }
        }
    }
#else
    for (int i = 0; i < vec_gbs; i++)
    {
        row = i * dimension;
        for (int j = 0; j < scal_total_dofs; j++)
        {
            col = j * dimension;
            for (int k = 0; k < dimension; k++)
            {
                grad(row + k, col + k) = grad_vector(i, j);
            }
        }
    }
#endif

    return grad;
}

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
compute_grad_matrix(const Mesh&                                           msh,
                    const typename Mesh::cell&                            cl,
                    const hho_degree_info&                                hdi,
                    const dynamic_matrix<typename Mesh::coordinate_type>& grad_vector)
{
    const CellDegreeInfo<Mesh> cell_infos(msh, cl, hdi.cell_degree(), hdi.face_degree(), hdi.grad_degree());

    return compute_grad_matrix(msh, cl, cell_infos, grad_vector);
}
} // end priv

} // end diskpp