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
#include "diskpp/bases/bases_utils.hpp"
#include "diskpp/methods/hho"
#include "diskpp/output/silo.hpp"

#define KXX 0.1
#define KXY 0.0
#define KYX 0.0
#define KYY 1.0

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

    disk::diffusion_tensor<mesh_type> diff_tens;

    rhs_functor() {
        diff_tens = disk::diffusion_tensor<mesh_type>::Identity();
    }

    rhs_functor(const disk::diffusion_tensor<mesh_type>& dt) {
        diff_tens = dt;
    }

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

    disk::diffusion_tensor<mesh_type> diff_tens;

    rhs_functor() {
        diff_tens = disk::diffusion_tensor<mesh_type>::Identity();
    }

    rhs_functor(const disk::diffusion_tensor<mesh_type>& dt) {
        diff_tens = dt;
    }

    /*
    scalar_type operator()(const point_type& pt) const
    {
        auto k00 = diff_tens(0,0);
        auto k01 = diff_tens(0,1);
        auto k10 = diff_tens(1,0);
        auto k11 = diff_tens(1,1);
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        auto cos_px = std::cos(M_PI * pt.x());
        auto cos_py = std::cos(M_PI * pt.y());

        return - M_PI * M_PI * ( (k01+k10) * cos_px * cos_py - (k00+k11) * sin_px * sin_py );
    }
    */

   scalar_type operator()(const point_type& pt) const
   {
        auto A = 0.01;
        auto kxx = diff_tens(0,0);
        auto kyy = diff_tens(1,1);
        auto x = pt.x();
        auto y = pt.y();
        using namespace std;
        auto ret =  -2*A*kyy*x*(-x + 1)*(exp(20*x) - 1)
            + kxx*(400*A*x*y*(-x + 1)*(-y + 1)*exp(20*x)
            - 40*A*x*y*(-y + 1)*exp(20*x)
            + 40*A*y*(-x + 1)*(-y + 1)*exp(20*x)
            - 2*A*y*(-y + 1)*(exp(20*x) - 1));
        return -ret;
   }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    disk::diffusion_tensor<mesh_type> diff_tens;

    rhs_functor() {
        diff_tens = disk::diffusion_tensor<mesh_type>::Identity();
    }

    rhs_functor(const disk::diffusion_tensor<mesh_type>& dt) {
        diff_tens = dt;
    }

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

template<typename Mesh>
auto make_rhs_function(const Mesh& msh, const disk::diffusion_tensor<Mesh>& dt)
{
    return rhs_functor<Mesh>(dt);
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

    /*
    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        return sin_px * sin_py;
    }
    */
    scalar_type operator()(const point_type& pt) const
    {
        auto A = 0.01;
        auto x = pt.x();
        auto y = pt.y();
        using namespace std;
        auto ret = A*x*y*(1-x)*(1-y)*(exp(20*x)-1);
        return ret;
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

/***************************************************************************/
/* Expected solution definition */
template<typename Mesh>
struct solution_gradient_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_gradient_functor< Mesh<T, 1, Storage> >
{
    typedef Mesh<T,1,Storage>                   mesh_type;
    typedef typename Eigen::Matrix<T,1,1>       gradient_type;
    typedef typename mesh_type::point_type      point_type;

    gradient_type operator()(const point_type& pt) const
    {
        gradient_type ret;
        ret(0,0) = M_PI * std::cos(M_PI * pt.x());
        return ret;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_gradient_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>                   mesh_type;
    typedef typename Eigen::Matrix<T,2,1>       gradient_type;
    typedef typename mesh_type::point_type      point_type;

    /*
    gradient_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        auto cos_px = std::cos(M_PI * pt.x());
        auto cos_py = std::cos(M_PI * pt.y());

        gradient_type ret;
        ret(0) = M_PI * cos_px * sin_py;
        ret(1) = M_PI * sin_px * cos_py;
        return ret;
    }
    */

    gradient_type operator()(const point_type& pt) const
    {
        auto A = 0.01;
        auto x = pt.x();
        auto y = pt.y();
        using namespace std;
        gradient_type ret;
        ret(0) = 20*A*x*y*(-x + 1)*(-y + 1)*exp(20*x)
            - A*x*y*(-y + 1)*(exp(20*x) - 1)
            + A*y*(-x + 1)*(-y + 1)*(exp(20*x) - 1);
        ret(1) = -A*x*y*(-x + 1)*(exp(20*x) - 1)
            + A*x*(-x + 1)*(-y + 1)*(exp(20*x) - 1);
        return ret;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_gradient_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>                   mesh_type;
    typedef typename Eigen::Matrix<T,3,1>       gradient_type;
    typedef typename mesh_type::point_type      point_type;

    gradient_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        auto sin_pz = std::sin(M_PI * pt.z());
        auto cos_px = std::cos(M_PI * pt.x());
        auto cos_py = std::cos(M_PI * pt.y());
        auto cos_pz = std::cos(M_PI * pt.z());

        Eigen::Matrix<T,3,1> ret;
        ret(0) = M_PI * cos_px * sin_py * sin_pz;
        ret(1) = M_PI * sin_px * cos_py * sin_pz;
        ret(2) = M_PI * sin_px * sin_py * cos_pz;
        return ret;
    }
};

template<typename Mesh>
auto make_solution_function(const Mesh& msh)
{
    return solution_functor<Mesh>();
}

template<typename Mesh>
auto make_solution_gradient_function(const Mesh& msh)
{
    return solution_gradient_functor<Mesh>();
}

using namespace disk;

template<typename T>
struct error_info {
    T   h;
    T   L2err;
    T   H1err;
    T   Aerr;
};

template<typename Mesh>
auto
run_hho_diffusion_solver(Mesh& msh, const hho_degree_info& hdi, const bool statcond,
    const bool stab_diam_F, const diffusion_tensor<Mesh>& diff_tens, const bool use_proj = false)
{
    using T = typename Mesh::coordinate_type;

    std::cout << "Running HHO(" << hdi.cell_degree() << ", " << hdi.face_degree();
    std::cout << ", " << hdi.reconstruction_degree() << ")" << std::endl;

    auto rhs_fun = make_rhs_function(msh, diff_tens);
    auto sol_fun = make_solution_function(msh);
    auto sol_grad = make_solution_gradient_function(msh);

    auto assembler = make_diffusion_assembler(msh, hdi);

    bool hdgstab = hdi.cell_degree() > hdi.face_degree();

    bool stabfree = hdi.reconstruction_degree() > hdi.face_degree()+1;
    std::cout << "Stabilization-free: " << std::boolalpha << stabfree << std::endl;
    std::cout << "Use projection in reconstruction: " << std::boolalpha << use_proj << std::endl;

    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> GR;
        if (use_proj) {
            auto oper = make_shlp(msh, cl, hdi, diff_tens);
            A = oper.second;
            GR = oper.first;
        }
        else {
            auto oper = make_scalar_hho_laplacian(msh, cl, hdi, diff_tens);
            A = oper.second;
            GR = oper.first;
        }

        if (not stabfree)
        {   
            if (hdgstab)
                A = A + make_scalar_hdg_stabilization(msh, cl, hdi, stab_diam_F);
            else
                A = A + make_scalar_hho_stabilization(msh, cl, GR, hdi, stab_diam_F);
        }

        Eigen::Matrix<T, Eigen::Dynamic, 1> rhs = make_rhs(msh, cl, cb, rhs_fun);
        auto sc = make_scalar_static_condensation(msh, cl, hdi, A, rhs);
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

    T Merr = 0.0;
    T L2err = 0.0;
    T L2norm = 0.0;
    T H1err = 0.0;
    T H1norm = 0.0;
    T Aerr = 0.0;
    T Anorm = 0.0;

    std::vector<T> silo_u, silo_uh, silo_f;

    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto rb = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
        auto rbs = rb.size();

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> GR;
        if (use_proj) {
            auto oper = make_shlp(msh, cl, hdi, diff_tens);
            A = oper.second;
            GR = oper.first;
        }
        else {
            auto oper = make_scalar_hho_laplacian(msh, cl, hdi, diff_tens);
            A = oper.second;
            GR = oper.first;
        }

        if (not stabfree)
        {   
            if (hdgstab)
                A = A + make_scalar_hdg_stabilization(msh, cl, hdi, stab_diam_F);
            else
                A = A + make_scalar_hho_stabilization(msh, cl, GR, hdi, stab_diam_F);
        }

        auto rhs = make_rhs(msh, cl, cb, rhs_fun);

        Eigen::Matrix<T, Eigen::Dynamic, 1> locsol =
            assembler.take_local_data(msh, cl, sol, sol_fun);

        Eigen::Matrix<T, Eigen::Dynamic, 1> fullsol =
            make_scalar_static_decondensation(msh, cl, hdi, A, rhs, locsol);

        Eigen::Matrix<T, Eigen::Dynamic, 1> realsol =
            project_function(msh, cl, hdi, sol_fun, 2);

        {
            Eigen::Matrix<T, Eigen::Dynamic, 1> dofs = fullsol.head( cb.size() );
            auto bar = barycenter(msh, cl);
            auto phi = cb.eval_functions( bar );
            auto uh = disk::eval(dofs, phi);
            auto u = sol_fun(bar);
            auto f = rhs_fun(bar);
            silo_u.push_back(u);
            silo_uh.push_back(uh);
            silo_f.push_back(f);
        }

        Matrix<T, Dynamic, 1> Rsol = GR*fullsol;

        auto qps = disk::integrate(msh, cl, 2*hdi.cell_degree()+2);
        for (auto& qp : qps)
        {
            Eigen::Matrix<T, Eigen::Dynamic, 1> dofs = fullsol.head( cb.size() );

            auto u = sol_fun( qp.point() );
            //auto phi = cb.eval_functions( qp.point() );
            //auto uh = disk::eval(dofs, phi);
            L2norm += qp.weight() * u * u;

            auto grad_u = sol_grad( qp.point() );
            auto dphi_tmp = rb.eval_gradients( qp.point() );
            Matrix<T, Dynamic, Mesh::dimension> dphi = dphi_tmp.block(1, 0, rbs-1, Mesh::dimension);
            auto grad_Ruh = disk::eval(Rsol, dphi);
            H1err += qp.weight() * (grad_u - grad_Ruh).dot( diff_tens*(grad_u - grad_Ruh) );
            H1norm += qp.weight() * grad_u.dot(diff_tens*grad_u);
        }

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MM = disk::make_mass_matrix(msh, cl, cb);
        auto diff = realsol - fullsol;
        Aerr += diff.dot(A*diff);
        Anorm += realsol.dot(A*realsol);
        auto diffT = diff.head(cb.size());
        Merr += diffT.dot(MM*diffT);
    }

    /*
    std::stringstream ss;
    ss << "stabfree_" << hdi.cell_degree() << "_" << hdi.face_degree() << "_" << hdi.reconstruction_degree();
    ss << ".silo";
    disk::silo_database db;
    db.create( ss.str() );

    disk::silo_zonal_variable<T> silo_uh_var("uh", silo_uh);
    disk::silo_zonal_variable<T> silo_u_var("u", silo_u);
    disk::silo_zonal_variable<T> silo_f_var("f", silo_f);

    db.add_mesh(msh, "mesh");
    db.add_variable("mesh", silo_uh_var);
    db.add_variable("mesh", silo_u_var);
    db.add_variable("mesh", silo_f_var);
    db.close();
    */
    

    error_info<T> ei;
    ei.h = disk::average_diameter(msh); 
    ei.L2err = std::sqrt(Merr)/std::sqrt(L2norm);
    ei.H1err = std::sqrt(H1err)/std::sqrt(H1norm);
    ei.Aerr = std::sqrt(Aerr)/std::sqrt(Anorm);

    return ei;
}

template<typename Mesh>
auto
run_hho_diffusion_solver(Mesh& msh, const hho_degree_info& hdi, const bool statcond,
    const bool stab_diam_F)
{
    diffusion_tensor<Mesh> dt = diffusion_tensor<Mesh>::Identity();
    return run_hho_diffusion_solver(msh, hdi, statcond, stab_diam_F, dt, false);
}

