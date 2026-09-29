#include <iostream>

#include "diskpp/common/simplicial_formula.hpp"
#include "diskpp/geometry/geometry.hpp"
#include "diskpp/mesh/mesh.hpp"
#include "diskpp/mesh/meshgen.hpp"

#define THRESH 1e-14

bool test_kahan_2d(void)
{
    using T = double;
    disk::triangular_mesh<T> msh;
    auto mesher = disk::make_simple_mesher(msh);
    mesher.refine();
    mesher.refine();

    T area = 0.0;
    for (auto& cl : msh) {
        auto pts = points(msh, cl);
        area += disk::area_triangle_kahan(pts[0], pts[1], pts[2]);
    }

    auto numelemes = msh.cells_size();
    auto error = std::abs(area - 1.0);
    std::cout << "Elements: " << numelemes << ", area error: " << error << std::endl;

    return error < THRESH;
}

bool test_kahan_3d(void)
{
    using T = double;
    disk::tetrahedral_mesh<T> msh;
    auto mesher = disk::make_simple_mesher(msh);
    mesher.refine();
    mesher.refine();

    T volume = 0.0;
    for (auto& cl : msh) {
        auto pts = points(msh, cl);
        volume += disk::volume_tetrahedron_kahan(pts[0], pts[1], pts[2], pts[3]);
    }

    auto numelemes = msh.cells_size();
    auto error = std::abs(volume - 1.0);
    std::cout << "Elements: " << numelemes << ", volume error: " << error << std::endl;

    return error < THRESH;
}

int main(void)
{
    using T = double;

    bool success = true;

    success &= test_kahan_2d();
    success &= test_kahan_3d();

    return (success == false);
}