/*
 *       /\        Matteo Cicuttin (C) 2016-2020
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2020                     nicolas.pignet@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */

#include <iomanip>
#include <iostream>
#include <regex>

#include <unistd.h>

#include "diskpp/bases/bases.hpp"

#include "diskpp/loaders/loader.hpp"

#include "common.hpp"

template<typename Mesh, typename Function>
void
do_testing_basis(std::vector<Mesh>& meshes,
                 const Function&    run_test,
                 size_t             min_test_degree = MIN_TEST_DEGREE,
                 size_t             max_test_degree = MAX_TEST_DEGREE)
{
    using T = typename Mesh::coordinate_type;

    for (size_t k = min_test_degree; k <= max_test_degree; k++)
    {
        std::cout << "  Testing degree " << k << std::endl;
        std::cout << "    mesh_h         scaled_mono_cl scaled_lege_cl scaled_mono_fc scaled_lege_fc" << std::endl;

        std::vector<T> mesh_hs;
        std::vector<std::array<T,4>> conds;

        for (auto& msh : meshes)
        {
            auto error = run_test(msh, k);
            mesh_hs.push_back(disk::average_diameter(msh));
            conds.push_back(error);
        }

        for (size_t i = 0; i < mesh_hs.size(); i++)
        {
            const auto cond = conds.at(i);
            std::cout << "    ";
            std::cout << std::scientific << std::setprecision(5) << mesh_hs.at(i) << "    ";
            std::cout << std::scientific << std::setprecision(5) << cond[0] << "    ";
            std::cout << std::scientific << std::setprecision(5) << cond[1] << "    ";
            std::cout << std::scientific << std::setprecision(5) << cond[2] << "    ";
            std::cout << std::scientific << std::setprecision(5) << cond[3] << "    ";
            std::cout << "     -- " << std::endl;
        }
    }
}

template<template<typename> class TestFunctor>
class tester2
{

    template<typename Mesh>
    TestFunctor<Mesh>
    get_test_functor(const std::vector<Mesh>& meshes)
    {
        return TestFunctor<Mesh>();
    }

    void
    test_triangles_generic(size_t min_degree = MIN_TEST_DEGREE, size_t max_degree = MAX_TEST_DEGREE)
    {
        std::cout << yellow << "Mesh under test: triangles on generic mesh";
        std::cout << nocolor << std::endl;
        using T = double;

        auto meshes = get_triangle_generic_meshes<T>();
        auto tf     = get_test_functor(meshes);

        do_testing_basis(meshes, tf, min_degree, max_degree);
    }

    void
    test_polygonal_generic(size_t min_degree = MIN_TEST_DEGREE, size_t max_degree = MAX_TEST_DEGREE)
    {
        std::cout << yellow << "Mesh under test: polygons on generic mesh";
        std::cout << nocolor << std::endl;
        using T = double;

        auto meshes = get_polygonal_generic_meshes<T>();
        auto tf     = get_test_functor(meshes);

        do_testing_basis(meshes, tf, min_degree, max_degree);
    }

    void
    test_triangles_netgen(size_t min_degree = MIN_TEST_DEGREE, size_t max_degree = MAX_TEST_DEGREE)
    {
        std::cout << yellow << "Mesh under test: triangles on netgen mesh";
        std::cout << nocolor << std::endl;
        using T = double;

        auto meshes = get_triangle_netgen_meshes<T>();
        auto tf     = get_test_functor(meshes);

        do_testing_basis(meshes, tf, min_degree, max_degree);
    }

    void
    test_quads(size_t min_degree = MIN_TEST_DEGREE, size_t max_degree = MAX_TEST_DEGREE)
    {
        std::cout << yellow << "Mesh under test: quads on generic mesh";
        std::cout << nocolor << std::endl;
        using T = double;

        auto meshes = get_quad_generic_meshes<T>();
        auto tf     = get_test_functor(meshes);

        do_testing_basis(meshes, tf, min_degree, max_degree);
    }

    void
    test_cartesian_2d_diskpp(size_t min_degree = MIN_TEST_DEGREE, size_t max_degree = MAX_TEST_DEGREE)
    {
        std::cout << yellow << "Mesh under test: 2D cartesian mesh (DiSk++)";
        std::cout << nocolor << std::endl;
        using T = double;

        auto meshes = get_cartesian_2d_diskpp_meshes<T>();
        auto tf     = get_test_functor(meshes);

        do_testing_basis(meshes, tf, min_degree, max_degree);
    }

    void
    test_tetrahedra_netgen(size_t min_degree = MIN_TEST_DEGREE, size_t max_degree = MAX_TEST_DEGREE)
    {
        std::cout << yellow << "Mesh under test: tetrahedra on netgen mesh";
        std::cout << nocolor << std::endl;
        using T = double;

        auto meshes = get_tetrahedra_netgen_meshes<T>();
        auto tf     = get_test_functor(meshes);

        do_testing_basis(meshes, tf, min_degree, max_degree);
    }

    void
    test_cartesian_3d_diskpp(size_t min_degree = MIN_TEST_DEGREE, size_t max_degree = MAX_TEST_DEGREE)
    {
        std::cout << yellow << "Mesh under test: 3D cartesian mesh (DiSk++)";
        std::cout << nocolor << std::endl;
        using T = double;

        auto meshes = get_cartesian_3d_diskpp_meshes<T>();
        auto tf     = get_test_functor(meshes);

        do_testing_basis(meshes, tf, min_degree, max_degree);
    }

    void
    test_generic_fvca6(size_t min_degree = MIN_TEST_DEGREE, size_t max_degree = MAX_TEST_DEGREE)
    {
        std::cout << yellow << "Mesh under test: polyhedra on generic mesh";
        std::cout << nocolor << std::endl;
        using T = double;

        auto meshes = get_generic_fvca6_meshes<T>();
        auto tf     = get_test_functor(meshes);

        do_testing_basis(meshes, tf, min_degree, max_degree);
    }

  public:
    int
    run(size_t min_degree = MIN_TEST_DEGREE, size_t max_degree = MAX_TEST_DEGREE)
    {
        sol::state lua;

        bool crash_on_nan           = false;
        bool do_triangles_generic   = true;
        bool do_polygonal_generic   = true;
        bool do_triangles_netgen    = true;
        bool do_quads               = true;
        bool do_cartesian_2d_diskpp = true;
        bool do_tetrahedra_netgen   = true;
        bool do_cartesian_3d_diskpp = true;
        bool do_generic_fvca6       = true;

        auto r = lua.do_file("test_config.lua");
        if (r.valid())
        {
            crash_on_nan           = lua["crash_on_nan"].get_or(false);
            do_triangles_generic   = lua["do_triangles_generic"].get_or(false);
            do_polygonal_generic   = lua["do_polygonal_generic"].get_or(false);
            do_triangles_netgen    = lua["do_triangles_netgen"].get_or(false);
            do_quads               = lua["do_quads"].get_or(false);
            do_cartesian_2d_diskpp = lua["do_cartesian_2d_diskpp"].get_or(false);
            do_tetrahedra_netgen   = lua["do_tetrahedra_netgen"].get_or(false);
            do_cartesian_3d_diskpp = lua["do_cartesian_3d_diskpp"].get_or(false);
            do_generic_fvca6       = lua["do_generic_fvca6"].get_or(false);
        }

        if (crash_on_nan)
            _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

        if (do_triangles_generic)
            test_triangles_generic(min_degree, max_degree);

        if (do_triangles_netgen)
            test_triangles_netgen(min_degree, max_degree);

        if (do_polygonal_generic)
            test_polygonal_generic(min_degree, max_degree);

        if (do_quads)
            test_quads(min_degree, max_degree);

        if (do_cartesian_2d_diskpp)
            test_cartesian_2d_diskpp(min_degree, max_degree);

        if (do_tetrahedra_netgen)
            test_tetrahedra_netgen(min_degree, max_degree);

        if (do_cartesian_3d_diskpp)
            test_cartesian_3d_diskpp(min_degree, max_degree);

        if (do_generic_fvca6)
            test_generic_fvca6(min_degree, max_degree);

        return 0;
    }
};

template<typename Mesh>
struct test_basis_cond
{
    typename std::array<typename Mesh::coordinate_type, 4>
    operator()(const Mesh& msh, size_t degree) const
    {
        typedef Mesh                                mesh_type;
        typedef typename mesh_type::cell            cell_type;
        typedef typename mesh_type::face            face_type;
        typedef typename mesh_type::coordinate_type scalar_type;
        typedef typename mesh_type::point_type      point_type;

        using Mat = Matrix<scalar_type, Dynamic, Dynamic>;

        scalar_type average_cl_sc = 0.0, average_cl_lg = 0.0;
        scalar_type average_fc_sc = 0.0, average_fc_lg = 0.0;
        size_t      nb_cells = 0, nb_faces = 0;

        for (auto& cl : msh)
        {
            auto cb_CSC   = make_scalar_monomial_basis(msh, cl, degree);
            auto mass_CSC = make_mass_matrix(msh, cl, cb_CSC);
            // std::cout << "Mass_CSC: " << std::endl;
            // std::cout << mass_CSC << std::endl;

            const auto mass_CSC_ev   = SelfAdjointEigenSolver<Mat>(mass_CSC).eigenvalues();
            average_cl_sc            += mass_CSC_ev(mass_CSC_ev.size() - 1) / mass_CSC_ev(0);

            auto cb_CLG   = make_scalar_legendre_basis(msh, cl, degree);
            auto mass_CLG = make_mass_matrix(msh, cl, cb_CLG);
            // std::cout << "Mass_CLG: " << std::endl;
            // std::cout << mass_CLG << std::endl;

            const auto mass_CLG_ev   = SelfAdjointEigenSolver<Mat>(mass_CLG).eigenvalues();
            average_cl_lg            += mass_CLG_ev(mass_CLG_ev.size() - 1) / mass_CLG_ev(0);

            // std::cout << "COND CL: " << mass_CSC_cond << " vs " << mass_CLG_cond << std::endl;

            nb_cells++;
        }

        for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
        {
            const auto fc      = *itor;
            auto       fb_SC   = make_scalar_monomial_basis(msh, fc, degree);
            auto mass_SC = make_mass_matrix(msh, fc, fb_SC);

            auto fb_LG   = make_scalar_legendre_basis(msh, fc, degree);
            auto mass_LG = make_mass_matrix(msh, fc, fb_LG);

            const auto mass_SC_ev = SelfAdjointEigenSolver<Mat>(mass_SC).eigenvalues();
            average_fc_sc += mass_SC_ev(mass_SC_ev.size() - 1) / mass_SC_ev(0);

            const auto mass_LG_ev = SelfAdjointEigenSolver<Mat>(mass_LG).eigenvalues();
            average_fc_lg += mass_LG_ev(mass_LG_ev.size() - 1) / mass_LG_ev(0);

            nb_faces++;
        }

        average_cl_sc /= nb_cells;
        average_cl_lg /= nb_cells;

        average_fc_sc /= nb_faces;
        average_fc_lg /= nb_faces;

        // std::cout << "Average conditionning: scaled_monomials vs scaled_legendre" << std::endl;
        // std::cout << "Cells : " << average_cl_sc << " vs " << average_cl_lg << std::endl;
        // std::cout << "Faces : " << average_fc_sc << " vs " << average_fc_lg << std::endl;

        return {average_cl_sc, average_cl_lg, average_fc_sc, average_fc_lg};
    }

    size_t
    expected_rate(size_t k)
    {
        return k + 1;
    }
};

int
main(int argc, char** argv)
{
    tester2<test_basis_cond> tstr1;
    tstr1.run(1, 6);
}