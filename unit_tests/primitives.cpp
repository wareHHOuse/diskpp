/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

/* Test primitive operations like barycenter(), measure() etc. */

#include <iostream>
#include <fstream>

#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/geometry/geometry.hpp"

#define MEAS_THRESH 1e-15
#define BARY_THRESH 1e-15

template<typename Mesh>
bool test_primitives(const Mesh& msh)
{
    using T = typename Mesh::coordinate_type;
    using point_type = typename Mesh::point_type;

    T domain_measure = 0.0;
    point_type domain_barycenter;

    for (auto& cl : msh)
    {
        auto bar = barycenter(msh, cl);
        auto meas = measure(msh, cl);
        domain_barycenter = domain_barycenter + bar*meas;
        domain_measure += meas;
    }
    domain_barycenter = domain_barycenter/domain_measure;

    bool success = true;

    auto meas_error = std::abs(1.0 - domain_measure);
    std::cout << "  Measure error = " << meas_error << std::endl;
    success &= (meas_error < MEAS_THRESH);

    point_type exp_barycenter;
    exp_barycenter.set_all(0.5);
    auto bary_error = distance(domain_barycenter, exp_barycenter);
    std::cout << "  Barycenter error = " << bary_error << std::endl;
    success &= (bary_error < BARY_THRESH);

    return success;
}

int main(void)
{
    using T = double;

    bool success = true;

    {
        std::cout << "Simplicial 2D" << std::endl;
        disk::simplicial_mesh<T,2> msh;
        auto mesher = make_simple_mesher(msh);
        mesher.refine();
        mesher.refine();
        mesher.refine();
        mesher.refine();
        success &= test_primitives(msh);
    }

    {
        std::cout << "Polygonal 2D" << std::endl;
        disk::generic_mesh<T,2> msh;
        auto mesher = make_fvca5_hex_mesher(msh);
        mesher.make_level(4);
        success &= test_primitives(msh);
    }

    return (success == false);
}