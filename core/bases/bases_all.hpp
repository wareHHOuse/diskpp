/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016 - 2020
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2020                     matteo.cicuttin@enpc.fr
 * Nicolas Pignet  (C) 2020                     nicolas.pignet@enpc.fr
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

#include "common/eigen.hpp"
#include "mesh/mesh.hpp"

using namespace Eigen;

namespace disk
{

/* Perform exponentiation by integer exponent. */
template<typename T>
T
iexp_pow(T x, size_t n)
{
    if (n == 0)
        return 1;

    T y = 1;
    while (n > 1)
    {
        if (n % 2 == 0)
        {
            x = x * x;
            n = n / 2;
        }
        else
        {
            y = x * y;
            x = x * x;
            n = (n - 1) / 2;
        }
    }

    return x * y;
}

template<typename Mesh, typename Element, typename ScalarType>
class scaled_monomial_abstract_basis;

template<typename Mesh, typename ScalarType>
class scaled_monomial_abstract_basis<Mesh, typename Mesh::cell, ScalarType>
{
  public:
    typedef Mesh                                mesh_type;
    typedef ScalarType                          scalar_type;
    typedef typename mesh_type::coordinate_type coordinate_type;
    typedef typename mesh_type::point_type      point_type;
    typedef typename mesh_type::cell            cell_type;

  private:
    typedef static_vector<coordinate_type, Mesh::dimension>                  vector_type;
    typedef static_matrix<coordinate_type, Mesh::dimension, Mesh::dimension> matrix_type;

    point_type                                   _bar;
    std::array<coordinate_type, Mesh::dimension> _length_box;
    vector_type                                  _scaling_factor;
    matrix_type                                  _axes;
    matrix_type                                  _passage;
    bool                                         _use_inertia_axes;

    void
    compute_scaling_factor()
    {
        for (size_t i = 0; i < Mesh::dimension; i++)
        {
            _scaling_factor(i) = coordinate_type(2) / _length_box[i];
        }
    }

    void
    compute_passage_matrix()
    {
        _passage = this->axes().transpose();
        for (size_t i = 0; i < Mesh::dimension; i++)
            _passage.row(i) *= _scaling_factor(i);
    }

  protected:
    /* This function maps a point to rescaled point axes * (pt - bar) / (2 * length_box) */
    point_type
    scaling_point(const point_type& pt) const
    {
        return point_type(this->passage_old2new() * (pt - _bar).to_vector());
    }

    // This matrix convert coordinates in the canonical basis to the user basis (rotation of axes)
    matrix_type
    passage_old2new() const
    {
        return _passage;
    }

    // This matrix convert coordinates in the user basis to the canonical basis (rotation of axes)
    matrix_type
    passage_new2old() const
    {
        return _passage.transpose();
    }

    matrix_type
    axes() const
    {
        return _axes;
    }

    std::array<coordinate_type, Mesh::dimension>
    length_box() const
    {
        return _length_box;
    }

    bool
    use_inertial_axes() const
    {
        return _use_inertia_axes;
    }

    vector_type
    scaling_factor() const
    {
        return _scaling_factor;
    }

  public:
    scaled_monomial_abstract_basis(const mesh_type& msh, const cell_type& cl, bool use_inertia_axes = true)
    {
        _bar              = barycenter(msh, cl);
        _use_inertia_axes = use_inertia_axes;
        if (_use_inertia_axes)
        {
            _axes = inertia_axes(msh, cl);
        }
        else
        {
            _axes = matrix_type::Identity();
        }
        _length_box = diameter_boundingbox(msh, cl, _axes);
        this->compute_scaling_factor();
        this->compute_passage_matrix();
    }

    bool
    is_orthonormal() const
    {
        return false;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename ScalarType>
class scaled_monomial_abstract_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::face, ScalarType>
{
  public:
    typedef Mesh<T, 2, Storage>                 mesh_type;
    typedef ScalarType                          scalar_type;
    typedef typename mesh_type::coordinate_type coordinate_type;
    typedef typename mesh_type::point_type      point_type;
    typedef typename mesh_type::face            face_type;

  private:
    typedef static_vector<coordinate_type, 2>    vector_type;
    typedef static_matrix<coordinate_type, 2, 2> matrix_type;

    point_type                           _bar;
    coordinate_type                      _length;
    vector_type                          _axes;
    bool                                 _use_inertia_axes;
    static_matrix<coordinate_type, 1, 2> _passage;

  protected:
    /* This function maps a 2D point on a face to a 1D reference system, to compute the
     * face basis. */
    point<coordinate_type, 1>
    map_face_point_2d_to_1d(const point_type& pt) const
    {
        return point<coordinate_type, 1>(_passage * (pt - _bar).to_vector());
    }

    /* This function maps a point to rescaled point axes * (pt - bar) / (2 * length_box) */
    point<coordinate_type, 1>
    scaling_point(const point_type& pt) const
    {
        return point<coordinate_type, 1>(_passage * (pt - _bar).to_vector());
    }

    vector_type
    axes() const
    {
        return _axes;
    }

    // This matrix convert coordinates in the canonical basis to the user basis (rotation of axes)
    static_matrix<coordinate_type, 1, 2>
    passage_old2new() const
    {
        return _passage;
    }

    // This matrix convert coordinates in the user basis to the canonical basis (rotation of axes)
    vector_type
    passage_new2old() const
    {
        return _passage.transpose();
    }

    coordinate_type
    length_box() const
    {
        return _length;
    }

    bool
    use_inertial_axes() const
    {
        return _use_inertia_axes;
    }

    coordinate_type
    scaling_factor() const
    {
        return 2.0 / _length;
    }

  public:
    scaled_monomial_abstract_basis(const mesh_type& msh, const face_type& fc, bool use_inertia_axes = true)
    {
        _bar              = barycenter(msh, fc);
        _use_inertia_axes = use_inertia_axes;
        const auto pts    = points(msh, fc);
        _axes             = (_bar - pts[0]).to_vector();
        _axes.normalize();
        _length  = diameter(msh, fc);
        _passage = _axes.transpose() * this->scaling_factor();
    }

    bool
    is_orthonormal() const
    {
        return false;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename ScalarType>
class scaled_monomial_abstract_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::face, ScalarType>
{
  public:
    typedef Mesh<T, 3, Storage>                 mesh_type;
    typedef ScalarType                          scalar_type;
    typedef typename mesh_type::coordinate_type coordinate_type;
    typedef typename mesh_type::point_type      point_type;
    typedef typename mesh_type::face            face_type;

  private:
    typedef static_vector<coordinate_type, 2>    vector_type;
    typedef static_matrix<coordinate_type, 2, 3> matrix_type;

    point_type                           _bar;
    std::array<coordinate_type, 2>       _length_box;
    vector_type                          _scaling_factor;
    static_matrix<coordinate_type, 3, 2> _axes;
    matrix_type                          _passage;
    bool                                 _use_inertia_axes;

    vector_type e0, e1;

    void
    compute_axes_box(const mesh_type& msh, const face_type& fc)
    {
        const auto axes_3d       = inertia_axes(msh, fc);
        const auto length_box_3d = diameter_boundingbox(msh, fc, axes_3d);

        // std::cout << "axes3D" << std::endl;
        // std::cout << axes_3d << std::endl;
        // std::cout << "box3D" << std::endl;
        // std::cout << length_box_3d[0] << ", " << length_box_3d[1] << ", " << length_box_3d[2] << std::endl;

        // search 0-value (min) of bounding box
        auto   min_val = length_box_3d[0];
        size_t min_ind = 0;

        for (size_t i = 1; i < 3; i++)
        {
            if (min_val > length_box_3d[i])
            {
                min_val = length_box_3d[i];
                min_ind = i;
            }
        }

        // remove axis i - orthogonal to the face
        size_t pos = 0;
        for (size_t i = 0; i < 3; i++)
        {
            if (i != min_ind)
            {
                _length_box[pos] = length_box_3d[i];
                _axes.col(pos)   = axes_3d.col(i);
                pos++;
            }
        }
        assert(pos == 2);

        // std::cout << "axes" << std::endl;
        // std::cout << _axes << std::endl;
        // std::cout << "box" << std::endl;
        // std::cout << _length_box[0] << ", " << _length_box[1] << std::endl;
    }

    void
    compute_scaling_factor()
    {
        _scaling_factor(0) = coordinate_type(2) / _length_box[0];
        _scaling_factor(1) = coordinate_type(2) / _length_box[1];
    }

    void
    compute_passage_matrix()
    {
        _passage = this->axes().transpose();
        _passage.row(0) *= _scaling_factor(0);
        _passage.row(1) *= _scaling_factor(1);
    }

  protected:
    /* This function maps a 3D point on a face to a 2D reference system, to compute the
     * face basis. */
    point<coordinate_type, 2>
    map_face_point_3d_to_2d(const point_type& pt) const
    {
        const vector_type map = (this->axes().transpose()) * (pt - _bar).to_vector();

        return point<coordinate_type, 2>({map(0), map(1)});
    }

    /* This function maps a point to rescaled point axes * (pt - bar) / (2 * length_box) */
    point<coordinate_type, 2>
    scaling_point(const point_type& pt) const
    {
        const vector_type map = (this->passage_old2new()) * (pt - _bar).to_vector();

        return point<coordinate_type, 2>({map(0), map(1)});
    }

    // This matrix convert coordinates in the canonical basis to the user basis (rotation of axes)
    static_matrix<coordinate_type, 2, 3>
    passage_old2new() const
    {
        return _passage;
    }

    // This matrix convert coordinates in the user basis to the canonical basis (rotation of axes)
    static_matrix<coordinate_type, 3, 2>
    passage_new2old() const
    {
        return _passage.transpose();
    }

    static_matrix<coordinate_type, 3, 2>
    axes() const
    {
        return _axes;
    }

    auto
    length_box() const
    {
        return _length_box;
    }

    bool
    use_inertial_axes() const
    {
        return _use_inertia_axes;
    }

    vector_type
    scaling_factor() const
    {
        return _scaling_factor;
    }

  public:
    scaled_monomial_abstract_basis(const mesh_type& msh, const face_type& fc, bool use_inertia_axes = true)
    {
        _bar              = barycenter(msh, fc);
        _use_inertia_axes = use_inertia_axes;
        this->compute_axes_box(msh, fc);
        this->compute_scaling_factor();
        this->compute_passage_matrix();
    }

    bool
    is_orthonormal() const
    {
        return false;
    }
};

template<typename Mesh, typename Element, typename ScalarType>
class scaled_monomial_abstract_face_basis;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, typename ScalarType>
class scaled_monomial_abstract_face_basis<Mesh<T, 3, Storage>, typename Mesh<T, 3, Storage>::face, ScalarType>
{
  public:
    typedef Mesh<T, 3, Storage>                 mesh_type;
    typedef ScalarType                          scalar_type;
    typedef typename mesh_type::coordinate_type coordinate_type;
    typedef typename mesh_type::point_type      point_type;
    typedef typename mesh_type::face            face_type;
    typedef Matrix<scalar_type, Dynamic, 1>     function_type;

  private:
    point_type      face_bar;
    coordinate_type face_h;

    /* Local reference frame */
    typedef static_vector<coordinate_type, 3> vector_type;
    vector_type                               e0, e1;

    /* It takes two edges of an element's face and uses them as the coordinate
     * axis of a 2D reference system. Those two edges are accepted only if they have an angle
     * between them greater than 8 degrees, then they are orthonormalized via G-S. */

    void
    compute_axis(const mesh_type& msh, const face_type& fc, vector_type& e0, vector_type& e1)
    {
        vector_type v0;
        vector_type v1;

        bool ok = false;

        const auto pts = points(msh, fc);

        const size_t npts = pts.size();
        for (size_t i = 1; i <= npts; i++)
        {
            const size_t i0 = (i + 1) % npts;
            const size_t i1 = (i - 1) % npts;
            v0              = (pts[i0] - pts[i]).to_vector();
            v1              = (pts[i1] - pts[i]).to_vector();

            const vector_type v0n = v0 / v0.norm();
            const vector_type v1n = v1 / v1.norm();

            if (v0n.dot(v1n) < 0.99) // we want at least 8 degrees angle
            {
                ok = true;
                break;
            }
        }

        if (!ok)
            throw std::invalid_argument("Degenerate polyhedron, cannot proceed");

        /* Don't normalize, in order to keep axes of the same order of lenght
         * of v in make_face_point_3d_to_2d() */
        e0 = v0; // / v0.norm();
        e1 = v1 - (v1.dot(v0) * v0) / (v0.dot(v0));
        e1 = e1; // / e1.norm();
    }

  protected:
    /* This function maps a 3D point on a face to a 2D reference system, to compute the
     * face basis. It takes two edges of an element's face and uses them as the coordinate
     * axis of a 2D reference system. Those two edges are accepted only if they have an angle
     * between them greater than 8 degrees, then they are orthonormalized via G-S. */
    point<T, 2>
    map_face_point_3d_to_2d(const point_type& pt) const
    {
        vector_type v = (pt - face_bar).to_vector();

        const auto eta = v.dot(e0);
        const auto xi  = v.dot(e1);

        return point<T, 2>({eta, xi});
    }

    auto
    face_barycenter() const
    {
        return face_bar;
    }

    auto
    face_diameter() const
    {
        return face_h;
    }

    auto
    reference_frame() const
    {
        return std::make_pair(e0, e1);
    }

  public:
    scaled_monomial_abstract_face_basis(const mesh_type& msh, const face_type& fc)
    {
        face_bar = barycenter(msh, fc);
        face_h   = diameter(msh, fc);
        compute_axis(msh, fc, e0, e1);
    }

    bool
    is_orthonormal() const
    {
        return false;
    }
};

}
