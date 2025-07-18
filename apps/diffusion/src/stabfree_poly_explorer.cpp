#include <iostream>



#include "stabfree_explorer.hpp"

void testfunc(int a) {
    std::cout << "testfunc: " << a << std::endl;
}

int run(const config& cfg) {
    using T = double;
    using mesh_type = disk::generic_mesh<T,2>;
    using point_type = mesh_type::point_type;

    mesh_type msh_gen;
    double scale = 1.0;
    disk::make_single_element_mesh(msh_gen, scale, cfg.num_vertices);

    explore(msh_gen, cfg);
    
    return 0;
}

extern "C" int luaopen_sfexp(lua_State* L) {

    std::cout << "Stabfree polygon explorer" << std::endl;

    sol::state_view lua(L);

    lua["run"] = &run;

    return sol::stack::call_lua(L, 1, init_config );
}
