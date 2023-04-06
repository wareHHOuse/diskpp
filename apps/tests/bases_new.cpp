/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <iostream>

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/bases/bases_new.hpp"

int main(void)
{
    using namespace disk;
    using namespace disk::basis;

    using T = double;
    using mesh_type = generic_mesh<T,1>;

    using point_type = typename mesh_type::point_type;
    using cell_type = typename mesh_type::cell_type;

    mesh_type msh;
    cell_type cl = *msh.cells_begin();

    auto u = scaled_monomial_basis(msh, cl, 7);
    auto v = scaled_monomial_basis(msh, cl, 7);

    point_type x;
    auto n = 2.0;
    std::cout << u(x) << std::endl;

    //std::cout << grad(v).dot(n)(x) << std::endl;

    //auto a = scalar_product( u, u );
}