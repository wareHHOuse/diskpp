#pragma once

namespace disk::lua {

template<typename T>
void register_feast_usertypes(sol::state& lua)
{
    using fep_t = disk::feast_eigensolver_params<T>;

    sol::usertype<fep_t> fept = lua.new_usertype<fep_t>(
        "feast_eigensolver_params",
        sol::constructors<fep_t()>()
    );

    fept["verbose"] = &fep_t::verbose;
    fept["tolerance"] = &fep_t::tolerance;
    fept["min_eigval"] = &fep_t::min_eigval;
    fept["max_eigval"] = &fep_t::max_eigval;
    fept["subspace_size"] = &fep_t::subspace_size;
    fept["max_iter"] = &fep_t::max_iter;
}

} // namespace disk::lua
