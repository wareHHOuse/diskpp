/* Simple postscript library - Matteo Cicuttin (c) 2016 */

#pragma once

#include <fstream>
#include <list>

#include "mesh/point.hpp"

namespace sps {

typedef int             simpleps_int;
typedef unsigned int    simpleps_uint;
typedef double          simpleps_real;

/*
template<typename T>
class simpleps_point
{
    T  m_x, m_y;

public:
    simpleps_point()
    {}

    simpleps_point(T x, T y)
        : m_x(x), m_y(y)
    {}

    T x() const { return m_x; }
    T y() const { return m_y; }
};
*/

typedef point<simpleps_real,2> ps_point;

class postscript_object
{
public:

                    postscript_object() {}
    virtual         ~postscript_object() {}
    virtual void    emit_code(std::ostream&) = 0;
};

class path : public postscript_object
{
    std::list<ps_point>  m_points;

public:
    path();

    path(std::initializer_list<ps_point> l)
        : m_points(l)
    {}

    void emit_code(std::ostream&);
};

class text : public postscript_object
{
    ps_point        m_position;
    std::string     m_text;

public:
    text();

    text(const std::string& txt, const ps_point& pos)
        : m_position(pos), m_text(txt)
    {}

    void emit_code(std::ostream&);
};

class simpleps
{
    std::list<postscript_object *>  m_ps_objs;

public:
    simpleps();
    ~simpleps();
    void add(postscript_object *);

    bool write(const std::string&);
    void output(std::ostream&);
};

} // namespace simpleps
