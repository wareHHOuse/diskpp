/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2020
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group
 */

#include <iomanip>
#include <iostream>
#include <regex>

#include <unistd.h>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/hho"
#include "diskpp/methods/implementation_hho/curl.hpp"


#include "common.hpp"

template<typename Mesh, bool mixed>
struct test_functor_curl_reconstruction
{
    typename Mesh::coordinate_type
    operator()(const Mesh& msh, size_t degree) const
    {
        return 0;
    }

    size_t
    expected_rate(size_t k)
    {
        return k + 50;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage, bool mixed>
struct test_functor_curl_reconstruction<Mesh<T,3,Storage>, mixed>
{
    using mesh_type = Mesh<T,3,Storage>;
    /* Expect k+1 convergence (hho stabilization) */
    typename mesh_type::coordinate_type
    operator()(const mesh_type& msh, size_t degree) const
    {
        typedef typename mesh_type::cell            cell_type;
        typedef typename mesh_type::face            face_type;
        typedef typename mesh_type::coordinate_type scalar_type;
        typedef typename mesh_type::point_type      point_type;

        typedef Matrix<scalar_type, mesh_type::dimension, 1> ret_type;

        auto f = [](const point_type& pt) -> Matrix<T,3,1> {
            Matrix<T,3,1> ret;
            ret(0) = std::sin(M_PI*pt.y());
            ret(1) = std::sin(M_PI*pt.z());
            ret(2) = std::sin(M_PI*pt.x());
            return ret;
        };

        auto sol = [](const point_type& pt) -> Matrix<T,3,1> {
            Matrix<T,3,1> ret;
            ret(0) = -M_PI*std::cos(M_PI*pt.z());
            ret(1) = -M_PI*std::cos(M_PI*pt.x());
            ret(2) = -M_PI*std::cos(M_PI*pt.y());
            return ret;
        };

        size_t fd = degree;
        size_t cd = mixed ? degree+1 : degree;
        size_t rd = degree;

        disk::hho_degree_info hdi(disk::priv::hdi_named_args{.rd = rd, .cd = cd, .fd = fd});

        scalar_type error = 0.0;
        for (const auto& cl : msh)
        {
            auto CR = disk::curl_reconstruction_pk(msh, cl, hdi);
            auto proj = disk::project_tangent(msh, cl, hdi, f, 1);
            Matrix<T, Dynamic, 1> rf = CR.first * proj;

            auto rb = disk::make_vector_monomial_basis(msh, cl, rd);

            auto qps = integrate(msh, cl, 2*std::max(rd,cd)+1);
            for (const auto& qp : qps)
            {
                auto cphi = rb.eval_functions(qp.point());
                Matrix<T,3,1> diff = disk::eval(rf, cphi) - sol(qp.point());
                error += qp.weight() * diff.dot(diff);
            }

        }

        return std::sqrt(error);
    }

    size_t
    expected_rate(size_t k)
    {
        return k + 1;
    }
};

template<typename Mesh>
using test_functor_curl_reconstruction_eo = test_functor_curl_reconstruction<Mesh, false>;

template<typename Mesh>
using test_functor_curl_reconstruction_mo = test_functor_curl_reconstruction<Mesh, true>;

int
main(void)
{
    std::cout << red << "Test HHO curl reconstruction operator" << std::endl;
    // face order: k, cell order: k
    std::cout << cyan << "Face order: k and Cell order: k" << std::endl;
    tester<test_functor_curl_reconstruction_eo> tstr1;
    tstr1.run();

    std::cout << cyan << "Face order: k and Cell order: k+1" << std::endl;
    tester<test_functor_curl_reconstruction_mo> tstr2;
    tstr2.run();

    return 0;
}
