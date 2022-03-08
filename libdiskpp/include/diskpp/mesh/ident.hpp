/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
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
#include <cstddef>

namespace disk{

template<typename T, typename impl, impl default_value>
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

    /* This *actually* wastes tons of memory. The specialization below must
     * be used, but before there are some places where hardcoded size_t's had
     * to be replaced with ident_raw_t */
};

template<typename T, uint32_t default_value>
class identifier<T, uint32_t, default_value>
{
    uint32_t    id_val;

public:
    typedef uint32_t    value_type;

    static identifier invalid() { return identifier(); }

    identifier() : id_val(default_value & 0x7FFFFFFF)
    {}

    identifier(const identifier&) = default;

    explicit identifier(uint32_t val) : id_val( (val & 0x7FFFFFFF) | 0x80000000 )
    {}

    operator uint32_t() const
    {
        if ( !(id_val & 0x80000000) )
            throw std::logic_error("Invalid identifier used");

        return id_val & 0x7FFFFFFF;
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
        return (this->id_val & 0x7FFFFFFF) < (other.id_val & 0x7FFFFFFF);
    }
};

template<typename T, class impl, impl default_value>
std::ostream&
operator<<(std::ostream& os, const identifier<T, impl, default_value>& id)
{
    os << impl(id);
    return os;
}

typedef size_t      ident_raw_t;

} // end disk
