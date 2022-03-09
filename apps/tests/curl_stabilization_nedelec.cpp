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
struct test_functor_curl_stab
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
struct test_functor_curl_stab<Mesh<T,3,Storage>, mixed>
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
            //ret(0) = std::sin(M_PI*pt.x());
            //ret(1) = std::sin(M_PI*pt.y());
            //ret(2) = std::sin(M_PI*pt.z());
            ret(0) = 0.0;
            ret(1) = 0.0;
            ret(2) = pt.x()*pt.y()*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
            //ret(2) = pt.x()*pt.x(); /* Beware of NaNs in std::sqrt */
            return ret;
        };

        size_t fd = degree+1;
        size_t cd = degree+1;
        size_t rd = degree;

        disk::hho_degree_info hdi(disk::priv::hdi_named_args{.rd = rd, .cd = cd, .fd = fd});

        scalar_type error = 0.0;
        for (auto& cl : msh)
        {
            auto stab = disk::curl_hdg_stabilization_nedelec(msh, cl, hdi);
            Matrix<scalar_type, Dynamic, 1> proj = disk::project_tangent_nedelec(msh, cl, hdi, f, 1);

            error += proj.dot(stab * proj);
        }

        return std::sqrt( std::abs(error) );
    }

    size_t
    expected_rate(size_t k)
    {
        return k + (mixed ? 1 : 0);
    }
};

template<typename Mesh>
using test_functor_curl_stab_eo = test_functor_curl_stab<Mesh, false>;

template<typename Mesh>
using test_functor_curl_stab_mo = test_functor_curl_stab<Mesh, true>;

int
main(void)
{
    std::cout << red << "Test HHO curl stabilization operator" << std::endl;
    // face order: k, cell order: k
    std::cout << cyan << "Face order: k and Cell order: k" << std::endl;
    tester<test_functor_curl_stab_eo> tstr1;
    tstr1.run();

    // face order: k, cell order: k
    std::cout << cyan << "Face order: k and Cell order: k+1" << std::endl;
    tester<test_functor_curl_stab_mo> tstr2;
    tstr2.run();

    return 0;
}
