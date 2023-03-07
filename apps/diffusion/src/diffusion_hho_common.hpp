/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */
/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2020, 2021
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group
 */
/*
 *       /\        Matteo Cicuttin (C) 2016-2019
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
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

#pragma once

#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/hho"

/***************************************************************************/
/* RHS definition */
template<typename Mesh>
struct rhs_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 1, Storage> >
{
    typedef Mesh<T,1,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        return M_PI * M_PI * sin_px;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        return 2.0 * M_PI * M_PI * sin_px * sin_py;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        auto sin_pz = std::sin(M_PI * pt.z());
        return 3.0 * M_PI * M_PI * sin_px * sin_py * sin_pz;
    }
};

template<typename Mesh>
auto make_rhs_function(const Mesh& msh)
{
    return rhs_functor<Mesh>();
}

/***************************************************************************/
/* Expected solution definition */
template<typename Mesh>
struct solution_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_functor< Mesh<T, 1, Storage> >
{
    typedef Mesh<T,1,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        return sin_px;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        return sin_px * sin_py;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        auto sin_pz = std::sin(M_PI * pt.z());
        return sin_px * sin_py * sin_pz;
    }
};

template<typename Mesh>
auto make_solution_function(const Mesh& msh)
{
    return solution_functor<Mesh>();
}

using namespace disk;

template<typename Mesh>
typename Mesh::coordinate_type
run_hho_diffusion_solver(const Mesh& msh, const hho_degree_info& hdi, const bool statcond, const bool stab_diam_F)
{
    using T = typename Mesh::coordinate_type;

    std::cout << "Running HHO(" << hdi.cell_degree() << ", " << hdi.face_degree();
    std::cout << ", " << hdi.reconstruction_degree() << ")" << std::endl;

    auto rhs_fun = make_rhs_function(msh);
    auto sol_fun = make_solution_function(msh);

    auto assembler = make_diffusion_assembler(msh, hdi);

    bool hdgstab = hdi.cell_degree() > hdi.face_degree();

    bool stabfree = hdi.reconstruction_degree() > hdi.face_degree()+1;
    std::cout << "Stabilization-free: " << std::boolalpha << stabfree << std::endl; 

    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr = make_scalar_hho_laplacian(msh, cl, hdi);

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = gr.second;

        if (not stabfree)
        {
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> stab;
            
            if (hdgstab)
                stab = make_scalar_hdg_stabilization(msh, cl, hdi, stab_diam_F);
            else
                stab = make_scalar_hho_stabilization(msh, cl, gr.first, hdi, stab_diam_F);

            A = A + stab;
        }

        Eigen::Matrix<T, Eigen::Dynamic, 1> rhs = make_rhs(msh, cl, cb, rhs_fun);
        auto sc     = make_scalar_static_condensation(msh, cl, hdi, A, rhs);
        assembler.assemble(msh, cl, sc.first, sc.second, sol_fun);
    }

    assembler.finalize();

    size_t systsz = assembler.LHS.rows();
    size_t nnz = assembler.LHS.nonZeros();

    std::cout << "Mesh has " << msh.cells_size() << " elements." << std::endl;
    std::cout << "System has " << assembler.LHS.rows() << " unknowns and ";
    std::cout << assembler.LHS.nonZeros() << " nonzeros." << std::endl;

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(systsz);

    std::cout << "Running MUMPS" << std::endl;
    sol = mumps_lu(assembler.LHS, assembler.RHS);

    T error = 0.0;

    //std::ofstream ofs("sol.dat");

    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());

        auto gr = make_scalar_hho_laplacian(msh, cl, hdi);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = gr.second;

        if (not stabfree)
        {
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> stab;
            
            if (hdgstab)
                stab = make_scalar_hdg_stabilization(msh, cl, hdi, stab_diam_F);
            else
                stab = make_scalar_hho_stabilization(msh, cl, gr.first, hdi, stab_diam_F);
            
            A = A + stab;
        }

        auto rhs = make_rhs(msh, cl, cb, rhs_fun);

        Eigen::Matrix<T, Eigen::Dynamic, 1> locsol =
            assembler.take_local_data(msh, cl, sol, sol_fun);

        Eigen::Matrix<T, Eigen::Dynamic, 1> fullsol = make_scalar_static_decondensation(msh, cl, hdi, A, rhs, locsol);

        Eigen::Matrix<T, Eigen::Dynamic, 1> realsol = project_function(msh, cl, hdi, sol_fun, 2);


        auto diff = realsol - fullsol;
        //error += diff.dot(A*diff);

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MM = make_mass_matrix(msh, cl, cb);

        error += diff.segment(0,cb.size()).dot(MM*diff.segment(0,cb.size()));

        //auto bar = barycenter(msh, cl);

        //for (size_t i = 0; i < Mesh::dimension; i++)
        //    ofs << bar[i] << " ";
        //ofs << fullsol(0) << std::endl;

    }

    std::cout << "h = " << disk::average_diameter(msh) << " ";
    std::cout << "err = " << std::sqrt(error) << std::endl;

    return std::sqrt(error);
}