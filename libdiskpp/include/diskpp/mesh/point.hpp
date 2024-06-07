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

#include "diskpp/common/eigen.hpp"
#include "ident.hpp"

namespace disk{

template<typename T, size_t DIM>
class point
{
    std::array<T, DIM>     m_coords;

public:

    typedef T                                   value_type;
    const static size_t                         dimension = DIM;

    typedef identifier<point, ident_raw_t, 0>  id_type;

    point()
    {
        for (size_t i = 0; i < DIM; i++)
            m_coords[i] = T(0);
    }

    point(const point& other) : m_coords(other.m_coords) {}

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

    point operator=(const point& other)
    {
        m_coords = other.m_coords;
        return *this;
    }

    void set_all(const T& val)
    {
        for (auto& c : m_coords)
            c = val;
    }

    point   operator-() const {
        auto ret = -1.0 * (*this);
        return ret;
    }

    template<typename U = T>
    typename std::enable_if<DIM == 1 || DIM == 2 || DIM == 3, U>::type
    x() const { return m_coords[0]; }

    template<typename U = T>
    typename std::enable_if<DIM == 1 || DIM == 2 || DIM == 3, U>::type&
    x() { return m_coords[0]; }

    template<typename U = T>
    typename std::enable_if<DIM == 2 || DIM == 3, U>::type
    y() const { return m_coords[1]; }

    template<typename U = T>
    typename std::enable_if<DIM == 2 || DIM == 3, U>::type&
    y() { return m_coords[1]; }

    template<typename U = T>
    typename std::enable_if<DIM == 3, U>::type
    z() const { return m_coords[2]; }

    template<typename U = T>
    typename std::enable_if<DIM == 3, U>::type&
    z() { return m_coords[2]; }

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
T distance(const point<T,DIM>& a, const point<T,DIM>& b)
{
    return (b-a).to_vector().norm();
}

template<typename T, size_t DIM>
point<T,DIM>
midpoint(const point<T,DIM>& a, const point<T,DIM>& b)
{
    return (a+b)/2.;
}

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
