/*
*       /\         DISK++, a template library for DIscontinuous SKeletal
*      /__\        methods.
*     /_\/_\
*    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
*   /__\  /__\     matteo.cicuttin@enpc.fr
*  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
*
* This file is copyright of the following authors:
* Karol Cascavita (C) 2018                     karol.cascavita@enpc.fr
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
#include "geometry/geometry.hpp"
#include "loaders/loader.hpp"
#include "revolution/methods/hho"
#include "solvers/solver.hpp"
 #include "common.hpp"
/***************************************************************************/
/* RHS definition */
template<typename Mesh>
struct rhs_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto cos_px = std::cos(M_PI * pt.x());
        auto cos_py = std::cos(M_PI * pt.y());

        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());

        //return 2.0 * M_PI * M_PI * cos_px * cos_py;
        return 2.0 * M_PI * M_PI * sin_px * sin_py;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
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
struct solution_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto cos_px = std::cos(M_PI * pt.x());
        auto cos_py = std::cos(M_PI * pt.y());

        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());

        //return cos_px * cos_py;
        return sin_px * sin_py;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
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

using namespace revolution;

#if 0
    template<typename Function>
    dynamic_vector<T>
    diffusion(const Mesh&  msh,
            const Function& rhs_fun,
            const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd)
    {
        hho_degree_info     hdi(ap.degree +1, ap.degree); //Also allow (degree + 1, degree)
        std::cout << " * cell degree :"<< hdi.cell_degree() << std::endl;
        std::cout << " * face degree :"<< hdi.face_degree() << std::endl;

        auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
        auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

        auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                        - fbs*msh.boundary_faces_size() ;

        auto offset_vector = full_offset(msh, hdi);

        dynamic_vector<T>  diff_sol = dynamic_vector<T>::Zero(num_full_dofs);

        auto cl_count = 0;
        auto assembler = make_diffusion_assembler2(msh, hdi, bnd);

        for (auto& cl : msh)
        {
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);

            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            matrix_type Ah  = gr.second + stab;

            auto sc = diffusion_static_condensation_compute(msh, cl, hdi, Ah, Lh);
            assembler.assemble(msh, cl, sc.first, sc.second);
            cl_count++;
        }

        assembler.impose_neumann_boundary_conditions(msh, bnd);
        assembler.finalize();

        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        //std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
        //std::cout << "Mesh faces: " << msh.faces_size() << std::endl;
        //std::cout << "Dofs: " << systsz << std::endl;

        dynamic_vector<T> sol = dynamic_vector<T>::Zero(systsz);

        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = true;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

        dynamic_vector<T> diffusion_sol = dynamic_vector<T>::Zero(num_full_dofs);

        cl_count = 0;
        for (auto& cl : msh)
        {
            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;

            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);

            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            matrix_type Ah  = gr.second + stab;

            vector_type u_faces = assembler.take_local_data(msh, cl, sol);
            vector_type u_full  =
            diffusion_static_condensation_recover(msh, cl, hdi, Ah, Lh, u_faces);

            diffusion_sol.block(cell_ofs, 0, num_total_dofs ,1) = u_full;
            cl_count++;
        }
        return diffusion_sol;
    }

#endif

template<typename Mesh>
auto
run_hho_diffusion_nitsche_par(const Mesh& msh,
    const algorithm_parameters<typename Mesh::scalar_type>& ap,
    const typename Mesh::scalar_type & eta)
{
    using T =  typename Mesh::scalar_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;


    hho_degree_info hdi(ap.degree + 1, ap.degree );
    //hho_degree_info hdi(ap.degree);

    //std::cout << " * CELL DEGREE =  "<<hdi.cell_degree()<< std::endl;
    //std::cout << " * FACE DEGREE =  "<<hdi.face_degree()<< std::endl;

    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    auto rhs_fun = make_rhs_function(msh);
    auto sol_fun = make_solution_function(msh);
    auto assembler = make_diffusion_assembler_nitsche(msh, hdi);

    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
        auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
        auto Lh     = make_rhs(msh, cl, cb, rhs_fun);
        matrix_type Ah = gr.second + stab;

        matrix_type Aconsist   = matrix_type::Zero(Ah.rows(), Ah.cols());
        matrix_type Anitsche  = matrix_type::Zero(Ah.rows(), Ah.cols());
        vector_type Bnitsche  = vector_type::Zero(Ah.cols());

        bool has_a_boundary_face = false;
        auto fcs = faces(msh, cl);
        for(auto fc : fcs)
        {
            if(msh.is_boundary(fc))
            {
                has_a_boundary_face = true;
                break;
            }
        }

        if (has_a_boundary_face)
        {
            Aconsist   = make_hho_consist_diff_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, eta);
            auto ntz   = make_hho_nitsche_diff_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, sol_fun, eta );
            Anitsche  = ntz.first;
            Bnitsche  = ntz.second;
        }

        matrix_type A = Ah - Anitsche - Aconsist;
        vector_type rhs = -Bnitsche;
        rhs.block(0, 0, cbs, 1) += Lh;
        auto sc = diffusion_static_condensation_compute_full(msh, cl, hdi, A, rhs);
        assembler.assemble(msh, cl, sc.first, sc.second);
    }

    assembler.finalize();

    size_t systsz = assembler.LHS.rows();
    size_t nnz = assembler.LHS.nonZeros();

    #if 0
    std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
    std::cout << "Mesh faces: " << msh.faces_size() << std::endl;
    std::cout << "Dofs: " << systsz << std::endl;
    #endif

    dynamic_vector<T> sol = dynamic_vector<T>::Zero(systsz);

    disk::solvers::pardiso_params<T> pparams;
    pparams.report_factorization_Mflops = true;
    mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

    dump_sparse_matrix(assembler.LHS , "Amat.dat");

    T H1_error = 0.0;
    T L2_error = 0.0;

    std::ofstream ofs("sol.dat");

    std::vector<T> inf_errors;

    for (auto& cl : msh)
    {
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
        auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
        auto Lh     = make_rhs(msh, cl, cb, rhs_fun);
        matrix_type Ah = gr.second + stab;

        matrix_type Aconsist  = matrix_type::Zero(Ah.rows(), Ah.cols());
        matrix_type Anitsche  = matrix_type::Zero(Ah.rows(), Ah.cols());
        vector_type Bnitsche  = vector_type::Zero(Ah.rows());

        bool has_a_boundary_face = false;
        auto fcs = faces(msh, cl);
        for(auto fc : fcs)
        {
            if(msh.is_boundary(fc))
            {
                has_a_boundary_face = true;
                break;
            }
        }

        if (has_a_boundary_face)
        {
            Aconsist   = make_hho_consist_diff_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, eta);
            auto ntz   = make_hho_nitsche_diff_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, sol_fun, eta );
            Anitsche  = ntz.first;
            Bnitsche  = ntz.second;
        }

        matrix_type A   = Ah - Anitsche - Aconsist;
        vector_type rhs = Lh - Bnitsche.block(0, 0, cbs, 1);

        vector_type locsol = assembler.take_local_data(msh, cl, sol);

        vector_type fullsol =
            diffusion_static_condensation_recover(msh, cl, hdi, A, rhs, locsol);

        vector_type realsol = project_function(msh, cl, hdi, sol_fun);

        auto diff = realsol - fullsol;
        H1_error += diff.dot(A*diff);

        matrix_type mass  = make_mass_matrix(msh, cl, cb, hdi.cell_degree());
        vector_type u_diff = diff.block(0, 0, cbs, 1);
        L2_error += u_diff.dot(mass * u_diff);

        vector_type ucell = fullsol.block(0, 0, cbs, 1);

        auto qps = integrate(msh, cl, 2 * hdi.cell_degree());

        for(auto qp : qps)
        {
            auto c_phi = cb.eval_functions(qp.point());
            auto realeval = sol_fun(qp.point());
            auto hho_eval = ucell.dot(c_phi);
            inf_errors.push_back( std::abs(hho_eval -realeval));
        }


        auto bar = barycenter(msh, cl);

        for (size_t i = 0; i < Mesh::dimension; i++)
            ofs << bar[i] << " ";
        ofs << fullsol(0) << std::endl;

    }

    auto itor = std::max_element(inf_errors.begin(), inf_errors.end());
    size_t distance = std::distance(inf_errors.begin(), itor);
    auto Linf_error = *std::next(inf_errors.begin(), distance);

    //std::cout << std::sqrt(error) << std::endl;

    ofs.close();

    H1_error = std::sqrt(H1_error);
    L2_error = std::sqrt(L2_error);

    auto error = std::make_tuple(H1_error, L2_error, Linf_error);
    return error;
}


#if 0
template<typename Mesh>
auto
run_hho_diffusion_nitsche(const Mesh& msh,
    const algorithm_parameters<typename Mesh::scalar_type>& ap)
{
    using T =  typename Mesh::scalar_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    hho_degree_info hdi(ap.degree + 1, ap.degree);

    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    auto rhs_fun = make_rhs_function(msh);
    auto sol_fun = make_solution_function(msh);
    auto assembler = make_diffusion_assembler_nitsche(msh, hdi);

    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
        //auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
        auto stab   = make_hdg_scalar_stabilization(msh, cl, hdi);
        auto Lh     = make_rhs(msh, cl, cb, rhs_fun);
        matrix_type Ah = gr.second + stab;

        matrix_type Aconsist  = matrix_type::Zero(Ah.rows(), Ah.cols());
        matrix_type Anitsche  = matrix_type::Zero(Ah.rows(), Ah.cols());
        vector_type Bnitsche  = vector_type::Zero(Ah.cols());

        bool has_a_boundary_face = false;
        auto fcs = faces(msh, cl);
        for(auto fc : fcs)
        {
            if(msh.is_boundary(fc))
            {
                has_a_boundary_face = true;
                break;
            }
        }

        if (has_a_boundary_face)
        {
            Aconsist   = make_hho_consist_diff(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta);
            auto ntz   = make_hho_nitsche_diff(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, sol_fun );
            Anitsche  = ntz.first;
            Bnitsche  = ntz.second;
        }

        matrix_type A = Ah - Anitsche - Aconsist;
        vector_type rhs = -Bnitsche;
        rhs.block(0, 0, cbs, 1) += Lh;
        auto sc = diffusion_static_condensation_compute_full(msh, cl, hdi, A, rhs);
        assembler.assemble(msh, cl, sc.first, sc.second);
    }

    assembler.finalize();

    size_t systsz = assembler.LHS.rows();
    size_t nnz = assembler.LHS.nonZeros();

    #if 0
    std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
    std::cout << "Mesh faces: " << msh.faces_size() << std::endl;
    std::cout << "Dofs: " << systsz << std::endl;
    #endif

    dynamic_vector<T> sol = dynamic_vector<T>::Zero(systsz);

    disk::solvers::pardiso_params<T> pparams;
    pparams.report_factorization_Mflops = true;
    mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

    T error = 0.0;

    std::ofstream ofs("sol.dat");

    for (auto& cl : msh)
    {
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
        //auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
        auto stab   = make_hdg_scalar_stabilization(msh, cl, hdi);
        auto Lh     = make_rhs(msh, cl, cb, rhs_fun);
        matrix_type Ah = gr.second + stab;

        matrix_type Aconsist   = matrix_type::Zero(Ah.rows(), Ah.cols());
        matrix_type Anitsche  = matrix_type::Zero(Ah.rows(), Ah.cols());
        vector_type Bnitsche  = vector_type::Zero(Ah.rows());

        bool has_a_boundary_face = false;
        auto fcs = faces(msh, cl);
        for(auto fc : fcs)
        {
            if(msh.is_boundary(fc))
            {
                has_a_boundary_face = true;
                break;
            }
        }

        if (has_a_boundary_face)
        {
            Aconsist   = make_hho_consist_diff(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta);
            auto ntz   = make_hho_nitsche_diff(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, sol_fun );
            Anitsche  = ntz.first;
            Bnitsche  = ntz.second;
        }

        matrix_type A   = Ah - Anitsche - Aconsist;
        vector_type rhs = Lh - Bnitsche.block(0, 0, cbs, 1);

        vector_type locsol = assembler.take_local_data(msh, cl, sol);

        vector_type fullsol =
            diffusion_static_condensation_recover(msh, cl, hdi, A, rhs, locsol);

        vector_type realsol = project_function(msh, cl, hdi, sol_fun);

        auto diff = realsol - fullsol;
        error += diff.dot(A*diff);

        auto bar = barycenter(msh, cl);

        for (size_t i = 0; i < Mesh::dimension; i++)
            ofs << bar[i] << " ";
        ofs << fullsol(0) << std::endl;

    }

    ofs.close();

    //std::cout << std::sqrt(error) << std::endl;
    return std::sqrt(error);
}
#endif


template<typename Mesh>
auto
run_hho_diffusion_nitsche_faces(const Mesh& msh,
    const algorithm_parameters<typename Mesh::scalar_type>& ap)
{
    using T =  typename Mesh::scalar_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    hho_degree_info hdi(ap.degree);

    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    auto rhs_fun = make_rhs_function(msh);
    auto sol_fun = make_solution_function(msh);
    auto assembler = make_diffusion_assembler_nitsche(msh, hdi);

    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
        auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
        auto Lh     = make_rhs(msh, cl, cb, rhs_fun);
        matrix_type Ah = gr.second + stab;

        matrix_type Aconsist   = matrix_type::Zero(Ah.rows(), Ah.cols());
        matrix_type Anitsche  = matrix_type::Zero(Ah.rows(), Ah.cols());
        vector_type Bnitsche  = vector_type::Zero(Ah.cols());

        bool has_a_boundary_face = false;
        auto fcs = faces(msh, cl);
        for(auto fc : fcs)
        {
            if(msh.is_boundary(fc))
            {
                has_a_boundary_face = true;
                break;
            }
        }

        if (has_a_boundary_face)
        {
            Aconsist   = make_hho_consist_diff_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta);
            auto ntz   = make_hho_nitsche_diff_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, sol_fun );
            Anitsche  = ntz.first;
            Bnitsche  = ntz.second;
        }

        matrix_type A = Ah - Anitsche - Aconsist;
        vector_type rhs = -Bnitsche;
        rhs.block(0, 0, cbs, 1) += Lh;
        auto sc = diffusion_static_condensation_compute_full(msh, cl, hdi, A, rhs);
        assembler.assemble(msh, cl, sc.first, sc.second);
    }

    assembler.finalize();

    size_t systsz = assembler.LHS.rows();
    size_t nnz = assembler.LHS.nonZeros();

    #if 0
    std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
    std::cout << "Mesh faces: " << msh.faces_size() << std::endl;
    std::cout << "Dofs: " << systsz << std::endl;
    #endif

    dynamic_vector<T> sol = dynamic_vector<T>::Zero(systsz);

    disk::solvers::pardiso_params<T> pparams;
    pparams.report_factorization_Mflops = true;
    mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

    T H1_error = 0.0;
    T L2_error = 0.0;

    std::ofstream ofs("sol.dat");

    std::vector<T> inf_errors;

    for (auto& cl : msh)
    {
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
        auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
        auto Lh     = make_rhs(msh, cl, cb, rhs_fun);
        matrix_type Ah = gr.second + stab;

        matrix_type Aconsist  = matrix_type::Zero(Ah.rows(), Ah.cols());
        matrix_type Anitsche  = matrix_type::Zero(Ah.rows(), Ah.cols());
        vector_type Bnitsche  = vector_type::Zero(Ah.rows());

        bool has_a_boundary_face = false;
        auto fcs = faces(msh, cl);
        for(auto fc : fcs)
        {
            if(msh.is_boundary(fc))
            {
                has_a_boundary_face = true;
                break;
            }
        }

        if (has_a_boundary_face)
        {
            Aconsist   = make_hho_consist_diff_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta);
            auto ntz   = make_hho_nitsche_diff_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, sol_fun );
            Anitsche  = ntz.first;
            Bnitsche  = ntz.second;
        }

        matrix_type A   = Ah - Anitsche - Aconsist;
        vector_type rhs = Lh - Bnitsche.block(0, 0, cbs, 1);

        vector_type locsol = assembler.take_local_data(msh, cl, sol);

        vector_type fullsol =
            diffusion_static_condensation_recover(msh, cl, hdi, A, rhs, locsol);

        vector_type realsol = project_function(msh, cl, hdi, sol_fun);

        auto diff = realsol - fullsol;
        H1_error += diff.dot(A*diff);

        matrix_type mass  = make_mass_matrix(msh, cl, cb, hdi.cell_degree());
        vector_type u_diff = diff.block(0, 0, cbs, 1);
        L2_error += u_diff.dot(mass * u_diff);

        vector_type ucell= fullsol.block(0,0, cbs,1);
        auto qps = integrate(msh, cl, 2 * hdi.cell_degree());

        for(auto qp : qps)
        {
            auto c_phi = cb.eval_functions(qp.point());
            auto realeval = sol_fun(qp.point());
            auto hho_eval = ucell.dot(c_phi);
            inf_errors.push_back( std::abs(hho_eval -realeval));
        }


        auto bar = barycenter(msh, cl);

        for (size_t i = 0; i < Mesh::dimension; i++)
            ofs << bar[i] << " ";
        ofs << fullsol(0) << std::endl;

    }

    auto itor = std::max_element(inf_errors.begin(), inf_errors.end());
    size_t distance = std::distance(inf_errors.begin(), itor);
    auto Linf_error = *std::next(inf_errors.begin(), distance);

    //std::cout << std::sqrt(error) << std::endl;

    ofs.close();

    H1_error = std::sqrt(H1_error);
    L2_error = std::sqrt(L2_error);

    auto error = std::make_tuple(H1_error, L2_error, Linf_error);
    return error;
}
template<typename Mesh, typename T>
auto
run_diffusion_solver(const Mesh& msh, const algorithm_parameters<T>& ap, const T& eta)
{
    std::tuple<T, T, T> error;
    switch (ap.solver)
    {
        case EVAL_IN_CELLS:
            //error = run_hho_diffusion_nitsche(msh, ap);
            break;
        case EVAL_ON_FACES:
            error = run_hho_diffusion_nitsche_faces(msh, ap);
            break;
        case EVAL_WITH_PARAMETER:
            error = run_hho_diffusion_nitsche_par(msh, ap, eta);
            break;
        case EVAL_IN_CELLS_AS_FACES:
            std::cout << "Not valid in this case. Choose faces (f) or cells( c)" << std::endl;
            break;


        default:
            throw std::invalid_argument("Invalid solver");
    }
    return error;
}
