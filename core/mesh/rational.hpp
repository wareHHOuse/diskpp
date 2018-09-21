/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
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

#pragma once

#include <iostream>
#include <cassert>
#include <cmath>

enum class sign {
	POSITIVE = 0,
	NEGATIVE = 1
};

template<typename T>
enum sign
sgn(const T& n)
{
	if (n >= T(0))
		return sign::POSITIVE;
	else
		return sign::NEGATIVE;
}

template<typename I>
I gcd(I a, I b)
{
    if (b == 0)
       return a;
    
    return gcd(b, a%b);
}

template<typename I>
I lcm(I a, I b)
{
	return (std::abs(a)/gcd(a,b))*std::abs(b);
}

template<typename I>
class rational
{
	I	m_num, m_den;
	
	void adjust_signs(void)
	{
		if (m_den >= 0)
            return;
        
        m_num = -m_num;
        m_den = -m_den;
	}
	
public:
	rational()
		: m_num(0), m_den(1)
	{}
	
    rational(I n)
		: m_num(n), m_den(1)
	{}
	
	rational(I n, I d)
		: m_num(n), m_den(d)
	{
		adjust_signs();
		assert(m_den > 0);
	}
	
	explicit operator float() const {
		return float(m_num)/float(m_den);
	}
	
	explicit operator double() const {
		return double(m_num)/double(m_den);
	}
	
	I num() const {
		return m_num;
	}
	
	I den() const {
		return m_den;
	}
	
	void simplify(void)
	{
		I g = gcd(m_num, m_den);
		m_num /= g;
		m_den /= g;
		adjust_signs();
		assert(m_den > 0);
	}
	
	rational& operator+=(const rational& other)
	{
		auto a = m_num;
		auto b = m_den;
		auto c = other.m_num;
		auto d = other.m_den;
		
		m_num = a*d + b*c;
		m_den = b*d;
		
		simplify();
		
		return *this;
	}
	
	rational operator+(const rational& other) const
	{
		rational ret = *this;
		ret += other;
		return ret;
	}
	
	rational& operator-=(const rational& other)
	{
		auto a = m_num;
		auto b = m_den;
		auto c = other.m_num;
		auto d = other.m_den;
		
		m_num = a*d - b*c;
		m_den = b*d;
		
		simplify();
		
		return *this;
	}
	
	rational operator-(const rational& other) const
	{
		rational ret = *this;
		ret -= other;
		return ret;
	}
	
	rational operator-() const
	{
		return rational(-m_num, m_den);
	}
	
	rational& operator*=(const rational& other)
	{
		auto a = m_num;
		auto b = m_den;
		auto c = other.m_num;
		auto d = other.m_den;
		
		auto g1 = gcd(a,d);
		auto g2 = gcd(b,c);
		
		a /= g1;
		b /= g2;
		c /= g2;
		d /= g1;
		
		m_num = a*c;
		m_den = b*d;
		assert(m_den > 0);
		return *this;
	}
    
    rational& operator*=(const I& other)
    {
        m_num *= other;
        simplify();
        return *this;
    }
	
	rational operator*(const rational& other) const
	{
		rational ret = *this;
		ret *= other;
		return ret;
	}
    
    rational operator*(const I& other) const
    {
        rational ret = *this;
        ret *= other;
        return ret;
    }
    
	rational operator/=(const rational& other)
	{
		auto r = *this * inv(other);
		m_num = r.num();
		m_den = r.den();
		return *this;
	}
    
    rational operator/=(I other)
    {
        m_den *= other;
        assert(m_den != 0);
        return *this;
    }
	
	rational operator/(const rational& other) const
	{
		return *this * inv(other);
	}
    
    rational operator/(I other) const
    {
        rational ret = *this;
        ret /= other;
        return ret;
    }
	
	bool operator<(const rational& other) const
	{
		auto a = m_num;
		auto b = m_den;
		auto c = other.m_num;
		auto d = other.m_den;
		
		return a*d < b*c;
	}
	
	bool operator<=(const rational& other) const
	{
		auto a = m_num;
		auto b = m_den;
		auto c = other.m_num;
		auto d = other.m_den;
		
		return a*d <= b*c;
	}
	
	bool operator>(const rational& other) const
	{
		return ! ((*this) <= other);
	}
	
	bool operator>=(const rational& other) const
	{
		return ! ((*this) < other);
	}
};

template<typename I>
rational<I>
abs(rational<I> r)
{
	if (r < rational<I>(0))
		return -r;
	
	return r;
}

template<typename I>
rational<I> inv(const rational<I>& r)
{
	return rational<I>(r.den(), r.num());
}

template<typename T>
std::ostream&
operator<<(std::ostream& os, const rational<T>& r)
{
	os << r.num() << "/" << r.den();
	return os;
}

