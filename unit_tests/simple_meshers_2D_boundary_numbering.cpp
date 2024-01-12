/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023, 2024
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <iostream>
#include <fstream>
#include <limits>

#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/geometry/geometry.hpp"
#include "diskpp/loaders/loader_utils.hpp"

template<typename Mesh>
bool test(Mesh& msh)
{
    // Check that the boundary number is the one specified by
    // loader_utils.hpp:boundary_number()
    for (const auto& cl : msh) {
        const auto fcs = faces(msh, cl);
        for (const auto& fc : fcs) {
            if ( not msh.is_boundary(fc) )
                continue;

            auto n = normal(msh, cl, fc);
            auto bi = msh.boundary_info(fc);
            auto bn = disk::boundary_number(n);
            if (not bn) {
                std::cout << "Can't compute boundary number from normal. ";
                std::cout << "Are faces parallel to xy, xz or yz planes?";
                std::cout << std::endl;
                return false;
            }
            if ( bi.id() != bn.value() ) {
                std::cout << "Wrong boundary numbering: bi.id() = " << bi.id();
                std::cout << ", bn.value() = " << bn.value() << ", normal = ";
                std::cout << n.transpose() << std::endl;
                return false;
            }
        }
    }

    return true;
}


int main(void)
{
    using T = double;

    std::cout << "Simplicial 2D" << std::endl;
    disk::simplicial_mesh<T, 2> msh_tri;
    auto mesher_tri = make_simple_mesher(msh_tri);
    mesher_tri.refine();

    if ( not test(msh_tri) ) {
        return 1;
    }

    std::cout << "Cartesian 2D" << std::endl;
    disk::cartesian_mesh<T, 2> msh_quad;
    auto mesher_quad = make_simple_mesher(msh_quad);
    mesher_quad.refine();

    if ( not test(msh_quad) ) {
        return 1;
    }

    std::cout << "FVCA5 hex-dominant 2D" << std::endl;
    disk::generic_mesh<T, 2> msh_poly;
    auto mesher_poly = make_fvca5_hex_mesher(msh_poly);
    mesher_poly.make_level(2);

    if ( not test(msh_poly) ) {
        return 1;
    }

    return 0;
}