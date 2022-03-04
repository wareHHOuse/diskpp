/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018, 2019
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
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

#include <iostream>
#include <array>
#include <stdexcept>
#include <initializer_list>

#include "common/eigen.hpp"
#include "ident.hpp"

namespace disk{

/** \class point
  * \ingroup Core_Module
  *
  * \brief Point type which contains the coordinates of a point \f$ x \in R^d \f$ (d is the dimension) in the euclidean space
  *
  *
  * The two template parameters are required:
  * \tparam T Numeric type, e.g. float, double...
  * \tparam DIM Number of dimension (DIM=2 in 2D and DIM=3 in 3D)
  *
  *
  */

template<typename T, size_t DIM>
class point
{
    std::array<T, DIM>     m_coords;

public:
    /** \brief value_type typedef (Numeric type).
      */
    typedef T                                   value_type;

    /** \brief dimension of the point.
      */
    const static size_t                         dimension = DIM;

    /** \brief id_type typedef.
      * \sa identifier class
      */
    typedef identifier<point, ident_raw_t, 0>  id_type;

    /** \brief default constructor
     * create a std::array<T, DIM> of size DIM
     *
     * The coordinates are initialized to zero
     *
     */
    point()
    {
        for (size_t i = 0; i < DIM; i++)
            m_coords[i] = T(0);
    }

    /** \brief constructor by copy
     *
     * copy a given point in this*
     *
     */
    point(const point& other) : m_coords(other.m_coords) {}

    /** \brief constructor with a list of value
     *
     * The coordinates are copied from the list of values l
     *
     * \note The size of the list has to be equal to DIM
     */
    point(std::initializer_list<T> l)
    {
        if (l.size() != DIM)
            throw std::invalid_argument("Wrong initializer list size");

        for (size_t i = 0; i < DIM; i++)
            m_coords[i] = *(l.begin()+i);

    }

    point(const static_vector<T, DIM>& vec)
    {
        for (size_t i = 0; i < DIM; i++)
            m_coords[i] = vec(i);
    }

    template<typename U = T>
    point(const typename std::enable_if<DIM == 1, U>::type& x)
    {
        m_coords[0] = x;
    }

    template<typename U = T>
    point(const typename std::enable_if<DIM == 2, U>::type& x, const U& y)
    {
        m_coords[0] = x;
        m_coords[1] = y;
    }

    template<typename U = T>
    point(const typename std::enable_if<DIM == 3, U>::type& x, const U& y, const U& z)
    {
        m_coords[0] = x;
        m_coords[1] = y;
        m_coords[2] = z;
    }

    T   at(size_t pos) const { return m_coords.at(pos); }
    T&  at(size_t pos)       { return m_coords.at(pos); }

    T   operator[](size_t pos) const { return m_coords[pos]; }
    T&  operator[](size_t pos)       { return m_coords[pos]; }

    /**
      *
      *  Assigns points to each other.
      *
      */

    point operator=(const point& other)
    {
        m_coords = other.m_coords;
        return *this;
    }

    /**
      *
      * \return Return a point with opposite sign of this*
      *
      */

    point   operator-() const {
        auto ret = -1.0 * (*this);
        return ret;
    }

    /**
      * \return Return the x-coordinate of this* (first coordinate)
      *
      * \note DIM has to be \f$ DIM \geq 1 \f$
      */

    template<typename U = T>
    typename std::enable_if<DIM == 1 || DIM == 2 || DIM == 3, U>::type
    x() const { return m_coords[0]; }

    /**
      * \return Return a reference to the x-coordinate of this* (first coordinate)
      *
      * \note DIM has to be \f$ DIM \geq 1 \f$
      */
    template<typename U = T>
    typename std::enable_if<DIM == 1 || DIM == 2 || DIM == 3, U>::type&
    x() { return m_coords[0]; }

    /**
      * \return Return the x-coordinate of this* (first coordinate)
      *
      * \note DIM has to be \f$ DIM \geq 2 \f$
      */

    template<typename U = T>
    typename std::enable_if<DIM == 2 || DIM == 3, U>::type
    y() const { return m_coords[1]; }

    /**
      * \return Return a reference to the y-coordinate of this* (second coordinate)
      *
      * \note DIM has to be \f$ DIM \geq 2 \f$
      */

    template<typename U = T>
    typename std::enable_if<DIM == 2 || DIM == 3, U>::type&
    y() { return m_coords[1]; }

    /**
      * \return Return the z-coordinate of this* (third coordinate)
      *
      * \note DIM has to be \f$ DIM = 3 \f$
      */
    template<typename U = T>
    typename std::enable_if<DIM == 3, U>::type
    z() const { return m_coords[2]; }

    /**
      * \return Return a reference to the z-coordinate of this* (third coordinate)
      *
      * \note DIM has to be \f$ DIM = 3 \f$
      */
    template<typename U = T>
    typename std::enable_if<DIM == 3, U>::type&
    z() { return m_coords[2]; }


    /**
      * \brief Convert a point in a static Eigen vector.
      *
      * \return Return a static Eigen vector which contains the coordinates of the point
      *
      */

    static_vector<T, DIM> to_vector() const
    {
        static_vector<T, DIM> ret;
        for (size_t i = 0; i < DIM; i++)
            ret(i) = m_coords[i];
        return ret;
    }

    friend point operator+(const point& p1, const point& p2)
    {
        point ret;
        for (size_t i = 0; i < DIM; i++)
            ret.m_coords[i] = p1.m_coords[i] + p2.m_coords[i];

        return ret;
    }

    template<int N>
    friend point
    operator+(const point& p1, const static_vector<T, N>& vec)
    {
        static_assert(N == DIM, "wrong dimension");
        point ret;
        for (size_t i = 0; i < DIM; i++)
            ret.m_coords[i] = p1.m_coords[i] + vec(i);

        return ret;
    }

    template<int N>
    friend point
    operator+(const static_vector<T, N>& vec, const point& p1)
    {
        return p1 + vec;
    }

    friend point operator-(const point& p1, const point& p2)
    {
        point ret;
        for (size_t i = 0; i < DIM; i++)
            ret.m_coords[i] = p1.m_coords[i] - p2.m_coords[i];

        return ret;
    }

    friend point operator*(const point& p, T scalefactor)
    {
        point ret;
        for (size_t i = 0; i < DIM; i++)
            ret.m_coords[i] = p.m_coords[i] * scalefactor;

        return ret;
    }

    friend point operator*(T scalefactor, const point& p)
    {
        return p * scalefactor;
    }

    friend point operator/(const point& p, T scalefactor)
    {
        point ret;
        for (size_t i = 0; i < DIM; i++)
            ret.m_coords[i] = p.m_coords[i] / scalefactor;

        return ret;
    }
};

template<typename T, size_t DIM>
static_vector<T, DIM>
to_vector(const point<T,DIM>& pt)
{
    static_vector<T, DIM> ret;
    for (size_t i = 0; i < DIM; i++)
        ret(i) = pt[i];
    return ret;
}

template<typename T>
T
det(const point<T,2>& p1, const point<T,2>& p2)
{
    return p1.x() * p2.y() - p1.y() * p2.x();
}

template<typename T, size_t DIM>
std::ostream&
operator<<(std::ostream& os, const point<T, DIM>& pt)
{
    os << "( ";
    for (size_t i = 0; i < DIM; i++)
    {
        os << pt[i];

        if (i < DIM-1)
            os << ", ";
    }
    os << " )";
    return os;
}

struct dummy_point_type {};
template<size_t DIM>
using point_identifier = ident_raw_t;

} // end disk
