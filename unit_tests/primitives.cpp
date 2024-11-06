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

#define MEAS_THRESH 5e-13
#define BARY_THRESH 5e-13

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

template<disk::mesh_2D Mesh>
bool test_area_and_normals(const Mesh& msh)
{
    double area_shoelace = 0.0;
    double area_divthm = 0.0;

    for (auto& cl : msh)
    {
        /* Use the shoelace formula to compute the area of the
         * polygon. This is sensible to the node ordering and
         * polygon orientation. It gives positive area if the
         * polygon vertices are oriented counterclockwise.
         */
        auto pts = points(msh, cl);
        for (size_t i = 0; i < pts.size(); i++) {
            auto p0 = pts[i];
            auto p1 = pts[(i+1)%pts.size()];
            area_shoelace += 0.5*(p0.x() - p1.x())*(p0.y() + p1.y());
        }
    
        /* Use the divergence theorem to compute the area of the
         * polygon. This is sensible to the direction of the normals.
         */
        auto fcs = faces(msh, cl);
        for (auto& fc : fcs) {
            auto bar = barycenter(msh, fc);
            auto meas = measure(msh, fc);
            auto n = normal(msh, cl, fc);
            area_divthm += 0.5 * meas * (bar.x()*n[0] + bar.y()*n[1]);
        }
    }
    auto s_error = std::abs(1.0 - area_shoelace);
    auto d_error = std::abs(1.0 - area_divthm);

    std::cout << "  Area via shoelace formula: " << area_shoelace << ", ";
    std::cout << "error: " << s_error << std::endl;
    std::cout << "  Area via divergence theorem: " << area_divthm << ", ";
    std::cout << "error: " << d_error << std::endl;

    

    return (s_error < MEAS_THRESH) and (d_error < MEAS_THRESH);
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

        success &= test_area_and_normals(msh);
    }

    {
        std::cout << "Simplicial 3D" << std::endl;
        disk::simplicial_mesh<T,3> msh;
        auto mesher = make_simple_mesher(msh);
        mesher.refine();
        mesher.refine();
        mesher.refine();
        mesher.refine();
        success &= test_primitives(msh);
    }

    return (success == false);
}