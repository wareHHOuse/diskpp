/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#pragma once

#include <iostream>
#include <cstddef>

template<typename T, class impl, impl default_value>
struct identifier
{
    typedef impl    value_type;

    static identifier invalid() { return identifier(); }

    identifier() : id_val(default_value), valid(false)
    {}

    identifier(const identifier&) = default;

    explicit identifier(impl val) : id_val(val), valid(true)
    {}

    operator impl() const
    {
        if (!valid)
            throw std::logic_error("Invalid identifier used");

        return id_val;
    }

    bool operator==(const identifier& other) const
    {
        return id_val == other.id_val;
    }

    bool operator!=(const identifier& other) const
    {
        return id_val != other.id_val;
    }

    bool operator<(const identifier& other) const
    {
        return this->id_val < other.id_val;
    }

private:
    impl    id_val;
    bool    valid;  /* This could waste tons of memory. A specialization for
                     * size_t using the MSB of id_val as flag could be useful.
                     */
};

template<typename T, class impl, impl default_value>
std::ostream&
operator<<(std::ostream& os, const identifier<T, impl, default_value>& id)
{
    os << impl(id);
    return os;
}

typedef size_t  ident_impl_t;
