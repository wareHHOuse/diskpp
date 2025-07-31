#pragma once

#include <sol/sol.hpp>
#include "point.hpp"

namespace disk {

template<typename T, size_t DIM>
void
register_point_usertype(sol::state& lua, const point<T,DIM>&)
{
    using point_type = point<T,DIM>;
    lua.new_usertype<point_type>("point",
        sol::constructors<point_type()
        >()
    );
}

template<typename T>
void
register_point_usertype(sol::state& lua, const point<T,1>&)
{
    using point_type = point<T,1>;
    lua.new_usertype<point_type>("point",
        sol::constructors<point_type(),
            point_type(const T&)
            >(),
        "x", &point_type::template x<T>
    );
}

template<typename T>
void
register_point_usertype(sol::state& lua, const point<T,2>&)
{
    using point_type = point<T,2>;
    lua.new_usertype<point_type>("point",
        sol::constructors<point_type(),
            point_type(const T&, const T&)
            >(),
        "x", sol::resolve<T(void) const>(&point_type::template x<T>),
        "y", sol::resolve<T(void) const>(&point_type::template y<T>)
    );
}

template<typename T>
void
register_point_usertype(sol::state& lua, const point<T,3>&)
{
    using point_type = point<T,3>;
    lua.new_usertype<point_type>("point",
        sol::constructors<point_type(),
            point_type(const T&, const T&, const T&)
            >(),
        "x", sol::resolve<T(void) const>(&point_type::template x<T>),
        "y", sol::resolve<T(void) const>(&point_type::template y<T>),
        "z", sol::resolve<T(void) const>(&point_type::template z<T>)
    );
}

} //namespace disk