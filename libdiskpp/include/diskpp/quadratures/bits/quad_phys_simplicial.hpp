/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */
 
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

#include "diskpp/common/simplicial_formula.hpp"
#include "quad_raw_gauss.hpp"
#include "quad_raw_triangle.hpp"
#include "quad_raw_tetra.hpp"

namespace disk
{

/**
 * @brief Return quadrature points on the interior of a physical 2D triangle
 *
 * @param msh Reference to the mesh
 * @param cl Refererence to the cell
 * @param degree Order of the quadrature
 * @return Quadrature points
 */
template<typename T>
std::vector<disk::quadrature_point<T, 2>>
integrate(const disk::simplicial_mesh<T, 2>&                msh,
          const typename disk::simplicial_mesh<T, 2>::cell& cl,
          const size_t                                      degree)
{
    const auto pts = points(msh, cl);
    assert(pts.size() == 3);
    return disk::quadrature::triangle_gauss(degree, pts);
}

/**
 * @brief Return quadrature points on a face of a physical 2D triangle
 *
 * @param msh Reference to the mesh
 * @param fc Refererence to the face
 * @param degree Order of the quadrature
 * @return Quadrature points
 */
template<typename T>
std::vector<disk::quadrature_point<T, 2>>
integrate(const disk::simplicial_mesh<T, 2>&                msh,
          const typename disk::simplicial_mesh<T, 2>::face& fc,
          const size_t                                      degree)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 2);
    return disk::quadrature::gauss_legendre(degree, pts[0], pts[1]);
}

/**
 * @brief Return quadrature points on the interior of a physical 3D tetrahedron
 *
 * @param msh Reference to the mesh
 * @param cl Refererence to the cell
 * @param degree Order of the quadrature
 * @return Quadrature points
 */
template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate(const disk::simplicial_mesh<T, 3>& msh, const typename disk::simplicial_mesh<T, 3>::face& fc, size_t degree)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 3);
    return disk::quadrature::triangle_gauss(degree, pts[0], pts[1], pts[2]);
}

/**
 * @brief Return quadrature points on a face of a physical 3D tetrahedron
 *
 * @param msh Reference to the mesh
 * @param fc Refererence to the face
 * @param degree Order of the quadrature
 * @return Quadrature points
 */
template<typename T>
std::vector<disk::quadrature_point<T, 3>>
integrate(const disk::simplicial_mesh<T, 3>& msh,
    const typename disk::simplicial_mesh<T, 3>::cell_type& cl,
    size_t degree)
{
    const auto pts = points(msh, cl);
    assert(pts.size() == 4);
    return disk::quadrature::arbq(degree, pts[0], pts[1], pts[2], pts[3]);
}

} // namespace disk

