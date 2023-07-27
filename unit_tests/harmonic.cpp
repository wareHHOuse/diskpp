/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

 /* Test for harmonic basis implementation: for each element integrate
  * grad(phi).dot(n) on the faces. For the divergence theorem this is
  * equal to the integral of the laplacian on the cell. But as the
  * functions are harmonic, this must be zero. */

#include <iostream>
#include <fstream>
#include <functional>
#include <silo.h>

#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/quadratures/quadratures.hpp"
#include "diskpp/bases/bases.hpp"

template<typename Mesh>
bool
test_divergence_theorem(const Mesh& msh)
{
    using T = typename Mesh::coordinate_type;
    bool success = true;

    size_t degree = 5;

    for (const auto& cl : msh)
    {
        auto hb = make_scalar_harmonic_top_basis(msh, cl, degree);
        hb.maximum_polynomial_degree(0);

        disk::dynamic_vector<T> acc = disk::dynamic_vector<T>::Zero(hb.size());

        auto fcs = faces(msh, cl);
        for (const auto& fc : fcs)
        {
            auto n = normal(msh, cl, fc);
            auto qps = disk::integrate(msh, fc, degree);
            for (auto& qp : qps)
            {
                Eigen::Matrix<T, Eigen::Dynamic, Mesh::dimension> grad =
                    hb.eval_gradients(qp.point());

                acc += qp.weight() * (grad*n);
            }
        }

        if (acc.sum() > 1e-14) {
            std::cout << acc.transpose() << std::endl;
            success = false;
        }
    } 
   
    return success;
}

bool test_simplicial_2D(void)
{
    std::cout << "Simplicial 2D" << std::endl;
    using T = double;
    disk::simplicial_mesh<T,2> msh;
    auto mesher = make_simple_mesher(msh);
    return test_divergence_theorem(msh);
}

bool test_generic_2D(void)
{
    std::cout << "Generic 2D" << std::endl;
    using T = double;
    disk::generic_mesh<T,2> msh;
    auto mesher = make_fvca5_hex_mesher(msh);
    mesher.make_level(0);
    return test_divergence_theorem(msh);
}

bool test_simplicial_3D(void)
{
    std::cout << "Simplicial 3D" << std::endl;
    using T = double;
    disk::simplicial_mesh<T,3> msh;
    auto mesher = make_simple_mesher(msh);
    return test_divergence_theorem(msh);
}

int main(void)
{
    bool success = true;

    success &= test_simplicial_2D();
    success &= test_generic_2D();
    success &= test_simplicial_3D();
    
    return success ? 0 : 1;
}