/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

/* Compute the condition number of the local mass matrices using different
 * scalings in the computation of the basis functions. */

#include <iostream>
#include <fstream>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/bases/bases_new.hpp"
#include "diskpp/bases/bases_operations.hpp"
#include "diskpp/quadratures/quadratures.hpp"

using namespace disk;
using namespace disk::basis;

template<typename T>
T cond(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m)
{
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    Eigen::JacobiSVD<matrix_type> svd(m);
    return svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
}

template<template<typename, size_t, typename> class Mesh, typename T, size_t DIM, typename Storage>
void
eval_conditioning(const Mesh<T, DIM, Storage>& msh, size_t degree)
{}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
eval_conditioning(const Mesh<T, 2, Storage>& msh, size_t degree)
{
    using mesh_type = Mesh<T, 2, Storage>;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    std::vector<T> cond1, cond2, cond3;
    std::vector<T> aniso_ratio;

    T s1 = 0.0;
    T s2 = 0.0;
    T s3 = 0.0;

    for (auto& cl : msh)
    {
        auto phi1 = scaled_monomial_basis(msh, cl, degree, rescaling_strategy::none);
        matrix_type mass1 = integrate(msh, cl, phi1, phi1);
        cond1.push_back( cond(mass1) );
        s1 += mass1(0,0);

        auto phi2 = scaled_monomial_basis(msh, cl, degree, rescaling_strategy::gram_schmidt);
        matrix_type mass2 = integrate(msh, cl, phi2, phi2);
        cond2.push_back( cond(mass2) );
        s2 += mass2(0,0);

        auto phi3 = scaled_monomial_basis(msh, cl, degree, rescaling_strategy::inertial);
        matrix_type mass3 = integrate(msh, cl, phi3, phi3);
        cond3.push_back( cond(mass3) );
        s3 += mass3(0,0);

        auto em = scaled_inertia_axes(msh, cl);
        auto ar = std::max(em.col(0).norm(), em.col(1).norm()) / std::min(em.col(0).norm(), em.col(1).norm());
        aniso_ratio.push_back(ar);
    }

    std::cout << s1 << " " << s2 << " " << s3 << std::endl;

    std::ofstream ofs("conditioning.txt");

    for (size_t i = 0; i < msh.cells_size(); i++)
    {
        ofs << i << " " << cond1[i] << " " << cond2[i] << " " << cond3[i] << " " << aniso_ratio[i] << std::endl;
    }
}

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        std::cout << argv[0] << " <mesh file> <degree>" << std::endl;
        return 1;
    }

    auto mesh_filename = argv[1];
    size_t degree = std::stoi(argv[2]);

    return disk::dispatch_all_meshes(mesh_filename,
              [](auto ...args) { eval_conditioning(args...); },
            degree);
}
