/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

 /* Test for quadrature rules: verify divergence theorem on each element */

#include <iostream>
#include <fstream>
#include <functional>
#include <silo.h>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/quadratures/quadratures.hpp"

template<typename Mesh>
struct test_field {};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct test_field<Mesh<T,2,Storage>> {
    using mesh_type = Mesh<T,2,Storage>;
    using point_type = typename mesh_type::point_type;

    size_t m, n;

    test_field(size_t pm, size_t pn)
        : m(pm), n(pn)
    {}

    Eigen::Matrix<T,2,1>
    operator()(const point_type& pt)
    {
        Eigen::Matrix<T,2,1> ret;
        ret(0) = std::pow(pt.x(), m)*std::pow(pt.y(), n);
        ret(1) = std::pow(pt.x(), n)*std::pow(pt.y(), m);
        return ret;
    }

    T divergence(const point_type& pt)
    {
        auto dFdx = (m > 0) ? m*std::pow(pt.x(), m-1)*std::pow(pt.y(), n) : 0;
        auto dFdy = (m > 0) ? m*std::pow(pt.x(), n)*std::pow(pt.y(), m-1) : 0;
        return dFdx + dFdy;
    }
};


template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct test_field<Mesh<T,3,Storage>> {
    using mesh_type = Mesh<T,3,Storage>;
    using point_type = typename mesh_type::point_type;

    size_t m, n, p;

    test_field(size_t pm, size_t pn, size_t pp)
        : m(pm), n(pn), p(pp)
    {}

    Eigen::Matrix<T,3,1>
    operator()(const point_type& pt)
    {
        Eigen::Matrix<T,3,1> ret;
        ret(0) = std::pow(pt.x(), m)*std::pow(pt.y(), n)*std::pow(pt.z(), p);
        ret(1) = std::pow(pt.x(), p)*std::pow(pt.y(), m)*std::pow(pt.z(), n);
        ret(2) = std::pow(pt.x(), n)*std::pow(pt.y(), p)*std::pow(pt.z(), m);
        return ret;
    }

    T divergence(const point_type& pt)
    {
        auto dFdx = (m > 0) ? m*std::pow(pt.x(), m-1)*std::pow(pt.y(), n)*std::pow(pt.z(), p) : 0;
        auto dFdy = (m > 0) ? m*std::pow(pt.x(), p)*std::pow(pt.y(), m-1)*std::pow(pt.z(), n) : 0;
        auto dFdz = (m > 0) ? m*std::pow(pt.x(), n)*std::pow(pt.y(), p)*std::pow(pt.z(), m-1) : 0;
        return dFdx + dFdy + dFdz;
    }
};

template<disk::mesh_2D Mesh>
bool
test_divergence_theorem(const Mesh& msh)
{
    using T = typename Mesh::coordinate_type;
    bool success = true;

    
    for (size_t m = 0; m < 8; m++)
    {
        for (size_t n = 0; n < 8; n++)
        {
            test_field<Mesh> F(m, n);
    
            for (auto& cl : msh)
            {
                /* Volume integral of div(F) */
                T vol_int = 0.0;
                auto qps = disk::integrate(msh, cl, m+n);
                for (auto& qp : qps)
                    vol_int += qp.weight() * F.divergence(qp.point());
                
                /* Face integral of F.n */
                T surf_int = 0.0;
                auto fcs = faces(msh, cl);
                for (auto& fc : fcs)
                {
                    auto norm = normal(msh, cl, fc);
                    auto f_qps = disk::integrate(msh, fc, m+n);
                    for (auto& qp : f_qps)
                        surf_int += qp.weight() * F(qp.point()).dot(norm);
                }

                auto err = std::abs(vol_int - surf_int);
                if (err > 2e-15) {
                    std::cout << m << " " << n << ": " << err << std::endl;
                    success = false;
                }
            }
        }
    }

    double totarea = 0.0;
    for (auto& cl : msh)
        totarea += measure(msh, cl);
    std::cout << "Global area: " << totarea << std::endl;

    return success;
}

template<disk::mesh_3D Mesh>
bool
test_divergence_theorem(const Mesh& msh)
{
    using T = typename Mesh::coordinate_type;
    bool success = true;

    
    for (size_t m = 0; m < 5; m++)
    {
        for (size_t n = 0; n < 5; n++)
        {
            for (size_t p = 0; p < 5; p++)
            {
                test_field<Mesh> F(m, n, p);
        
                for (auto& cl : msh)
                {
                    /* Volume integral of div(F) */
                    T vol_int = 0.0;
                    auto qps = disk::integrate(msh, cl, m+n+p);
                    for (auto& qp : qps)
                        vol_int += qp.weight() * F.divergence(qp.point());
                    
                    /* Face integral of F.n */
                    T surf_int = 0.0;
                    auto fcs = faces(msh, cl);
                    for (auto& fc : fcs)
                    {
                        auto norm = normal(msh, cl, fc);
                        auto f_qps = disk::integrate(msh, fc, m+n+p);
                        for (auto& qp : f_qps)
                            surf_int += qp.weight() * F(qp.point()).dot(norm);
                    }

                    auto err = std::abs(vol_int - surf_int);
                    if (err > 2e-15) {
                        std::cout << m << " " << n << " " << p << ": " << err << std::endl;
                        success = false;
                    }
                }
            }
        }
    }
    
    double totvol = 0.0;
    for (auto& cl : msh)
        totvol += measure(msh, cl);
    std::cout << "Global volume: " << totvol << std::endl;

    return success;
}

bool test_simplicial_2D(void)
{
    std::cout << "Simplicial 2D" << std::endl;
    using T = double;
    disk::simplicial_mesh<T,2> msh;
    auto mesher = make_simple_mesher(msh);
    mesher.refine();
    mesher.refine();
    return test_divergence_theorem(msh);
}

bool test_generic_2D(void)
{
    std::cout << "Generic 2D" << std::endl;
    using T = double;
    disk::generic_mesh<T,2> msh;
    auto mesher = make_fvca5_hex_mesher(msh);
    mesher.make_level(2);
    return test_divergence_theorem(msh);
}


bool test_simplicial_3D(void)
{
    std::cout << "Simplicial 3D internal mesher" << std::endl;
    using T = double;
    disk::simplicial_mesh<T,3> msh;
    auto mesher = make_simple_mesher(msh);
    mesher.refine();
    return test_divergence_theorem(msh);
}

bool test_fvca6(const std::string& filename)
{
    using T = double;
    std::cout << "FVCA6 " + filename << std::endl;
    disk::generic_mesh<T,3> msh;
    disk::load_mesh_fvca6_3d<T>(filename.c_str(), msh);
    return test_divergence_theorem(msh);
}

int main(void)
{
    bool success = true;

    success &= test_simplicial_2D();
    success &= test_generic_2D();
    success &= test_simplicial_3D();
    success &= test_fvca6("test_data/tet.1.msh");
    success &= test_fvca6("test_data/dbls_10.msh");
    success &= test_fvca6("test_data/dbls_20.msh");

    return success ? 0 : 1;
}