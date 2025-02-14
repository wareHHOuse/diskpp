/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#pragma once

namespace disk::lua {

struct hho_parameters {
    size_t order;
    double stabparam;

    hho_parameters()
        : order(0), stabparam(1.0)
    {}
};

void register_hho_usertypes(sol::state& lua)
{
    sol::usertype<hho_parameters> hpt = lua.new_usertype<hho_parameters>(
        "hho_parameters",
        sol::constructors<hho_parameters()>()
    );
    hpt["order"] = &hho_parameters::order;
    hpt["stabparam"] = &hho_parameters::stabparam;
}

} //namespace disk::lua