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
#include "methods/hho"
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
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto cos_px = std::cos(M_PI * pt.x());
        auto cos_py = std::cos(M_PI * pt.y());

        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());

        return 2.0 * M_PI * M_PI * cos_px * cos_py;
        //return 2.0 * M_PI * M_PI * sin_px * sin_py;
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
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto cos_px = std::cos(M_PI * pt.x());
        auto cos_py = std::cos(M_PI * pt.y());

        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());

        return cos_px * cos_py;
        //return sin_px * sin_py;
    }
};

using namespace disk;

template<typename Mesh>
auto make_solution_function(const Mesh& msh)
{
    return solution_functor<Mesh>();
}

template<typename Mesh>
auto
run_hho_diffusion_nitsche_faces(const Mesh& msh,
    const algorithm_parameters<typename Mesh::coordinate_type>& ap)
{
    using T =  typename Mesh::coordinate_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    hho_degree_info hdi(ap.degree, ap.degree );

    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    auto rhs_fun = make_rhs_function(msh);
    auto sol_fun = make_solution_function(msh);
    auto assembler = make_diffusion_assembler_nitsche_faces(msh, hdi);

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


    dynamic_vector<T> sol = dynamic_vector<T>::Zero(systsz);

    disk::solvers::pardiso_params<T> pparams;
    pparams.report_factorization_Mflops = true;
    mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

    T H1_error = 0.0;
    T L2_error = 0.0;

    std::ofstream ofs("sol_faces_solver.dat");

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

        vector_type realsol = project_function(msh, cl, hdi, sol_fun, 2);

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

    ofs.close();

    H1_error = std::sqrt(H1_error);
    L2_error = std::sqrt(L2_error);

    auto error = std::make_tuple(H1_error, L2_error, Linf_error);
    return error;
}


template<typename Mesh>
auto
run_hho_diffusion_nitsche_cells_full(const Mesh& msh,
    const algorithm_parameters<typename Mesh::coordinate_type>& ap,
    const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd)
{
    using T =  typename Mesh::coordinate_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    hho_degree_info hdi(ap.degree + 1, ap.degree);

    auto is_contact_vector = make_is_contact_vector(msh, bnd);

    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension - 1);

    auto rhs_fun = make_rhs_function(msh);
    auto sol_fun = make_solution_function(msh);

    // Mixed assembler: to allow Dirichlet with Nitsche's method and by imposition
    //                  of the values on the boundary faces
    auto assembler = make_mix_full_assembler(msh, hdi, bnd);

    auto cl_count = 0;

    for (auto& cl : msh)
    {
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        const auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;

        matrix_type A  = matrix_type::Zero(num_total_dofs, num_total_dofs);
        vector_type fullsol  = vector_type::Zero(num_total_dofs);

        if (is_contact_vector.at(cl_count) == 1)
        {
            auto gr   = make_hho_contact_scalar_laplacian(msh, cl, hdi, bnd);
            auto stab = make_hdg_contact_stabilization(msh, cl, hdi, bnd);

            matrix_type Ah  = gr.second + stab;
            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun);//, hdi.cell_degree());

            matrix_type Aconsist = make_hho_consist_mix(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd);
            auto ntz  = make_hho_nitsche_mix(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, sol_fun, bnd);
            matrix_type Anitsche  = ntz.first;
            vector_type Bnitsche  = ntz.second;

            matrix_type A = Ah - Anitsche - Aconsist;
            vector_type rhs = -Bnitsche;
            rhs.block(0, 0, cbs, 1) += Lh;
            assembler.assemble(msh, cl, A, rhs);
        }
        else
        {
            auto gr   = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab = make_hdg_scalar_stabilization(msh, cl, hdi);

            vector_type Lh = make_rhs(msh, cl, cb, rhs_fun);//, hdi.cell_degree());
            matrix_type Ah = gr.second + stab;

            assembler.assemble(msh, cl, Ah, Lh);
        }
        cl_count++;

    }
    assembler.impose_neumann_boundary_conditions(msh, bnd);
    assembler.finalize();

    size_t systsz = assembler.LHS.rows();
    size_t nnz = assembler.LHS.nonZeros();

    dynamic_vector<T> sol = dynamic_vector<T>::Zero(systsz);

    disk::solvers::pardiso_params<T> pparams;
    pparams.report_factorization_Mflops = true;
    mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

    T H1_error = 0.0;
    T L2_error = 0.0;
    cl_count = 0;

    std::ofstream ofs("sol_cfull_solver.dat");

    std::vector<T> inf_errors;

    for (auto& cl : msh)
    {

        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        const auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;

        matrix_type A  = matrix_type::Zero(num_total_dofs, num_total_dofs);
        vector_type fullsol  = vector_type::Zero(num_total_dofs);

        if (is_contact_vector.at(cl_count) == 1)
        {
            auto gr   = make_hho_contact_scalar_laplacian(msh, cl, hdi, bnd);
            auto stab = make_hdg_contact_stabilization(msh, cl, hdi, bnd);

            matrix_type Ah  = gr.second + stab;
            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun);//, hdi.cell_degree());

            matrix_type Aconsist = make_hho_consist_mix(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta,  bnd);
            auto ntz  = make_hho_nitsche_mix(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, sol_fun, bnd );
            matrix_type Anitsche  = ntz.first;
            vector_type Bnitsche  = ntz.second;

            A = Ah - Anitsche - Aconsist;
            vector_type rhs = Lh - Bnitsche.block(0, 0, cbs, 1);

            fullsol = assembler.take_local_data(msh, cl, sol);
        }
        else
        {
            auto gr   = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab = make_hdg_scalar_stabilization(msh, cl, hdi);

            vector_type rhs = make_rhs(msh, cl, cb, rhs_fun);//, hdi.cell_degree());
            A = gr.second + stab;

            vector_type sol_faces = assembler.take_local_data(msh, cl, sol);
            fullsol = assembler.take_local_data(msh, cl, sol);
        }

        vector_type realsol = project_function(msh, cl, hdi, sol_fun, 2);

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

        auto fcs = faces(msh, cl);
        for(size_t fi = 0; fi < fcs.size(); fi++)
        {
            auto fc = fcs[fi];
            auto face_id = msh.lookup(fc);

            if(bnd.is_contact_face(face_id))
                continue;

            auto fbar = barycenter(msh, fc);
            for (size_t i = 0; i < Mesh::dimension; i++)
                ofs << fbar[i] << " ";
            auto idx = cbs + fbs * fi;
            ofs << fullsol(idx) << std::endl;
        }

        cl_count++;
    }
    ofs.close();

    auto itor = std::max_element(inf_errors.begin(), inf_errors.end());
    size_t distance = std::distance(inf_errors.begin(), itor);
    auto Linf_error = *std::next(inf_errors.begin(), distance);

    H1_error = std::sqrt(H1_error);
    L2_error = std::sqrt(L2_error);

    auto error = std::make_tuple(H1_error, L2_error, Linf_error);

    return error;
}


template<typename Mesh, typename T>
auto
run_diffusion_solver(const Mesh& msh, const algorithm_parameters<T>& ap)
{

    dump_to_matlab(msh,"mesh.m");

    typedef disk::mechanics::BoundaryConditionsScalar<Mesh> boundary_type;
    boundary_type  bnd(msh);

    auto sol_fun = make_solution_function(msh);

    auto zero_fun = [](const typename Mesh::point_type& p) -> T {
        return 0.;
    };

    auto neu_left = [](const typename Mesh::point_type& p) -> T {
        auto cos_px = std::cos(M_PI * p.x());
        auto sin_py = std::sin(M_PI * p.y());

        auto sin_px = std::sin(M_PI * p.x());
        auto cos_py = std::cos(M_PI * p.y());

        return  M_PI * sin_px * cos_py;
        //return  M_PI * cos_px * sin_py;
    };

    auto neu_right = [](const typename Mesh::point_type& p) -> T {
        auto cos_px = std::cos(M_PI * p.x());
        auto sin_py = std::sin(M_PI * p.y());

        auto sin_px = std::sin(M_PI * p.x());
        auto cos_py = std::cos(M_PI * p.y());

        return  - M_PI * sin_px * cos_py;
        //return  -M_PI * cos_px * sin_py;
    };

    /* DIRICHLET conditions can be enforced in strong way using addDirichletBC
    *  or in a weak form with Nitsche's method using addContactBC
    */
    //test 1:
    #if 0
    bnd.addDirichletBC(disk::mechanics::DIRICHLET,1, sol_fun);
    bnd.addContactBC(disk::mechanics::SIGNORINI, 3);
    bnd.addNeumannBC(disk::mechanics::NEUMANN, 4, neu_left);
    bnd.addNeumannBC(disk::mechanics::NEUMANN, 2, neu_right);
    #endif

    //test 2 : All boundaries with DIRICHLET enforced via Nitsche's method
    bnd.addContactBC(disk::mechanics::SIGNORINI, 1);
    bnd.addContactBC(disk::mechanics::SIGNORINI, 2);
    bnd.addContactBC(disk::mechanics::SIGNORINI, 3);
    bnd.addContactBC(disk::mechanics::SIGNORINI, 4);

    std::tuple<T, T, T> error;
    switch (ap.solver)
    {
        case EVAL_ON_FACES:
            error = run_hho_diffusion_nitsche_faces(msh,  ap); // Just for test 1.
            break;
        case EVAL_IN_CELLS_FULL:
            error = run_hho_diffusion_nitsche_cells_full(msh, ap, bnd);
            break;
        case EVAL_IN_CELLS:
        case EVAL_IN_CELLS_AS_FACES:
        default:
            std::cout << "Not valid in this case. Choose only full cells (l)" << std::endl;
            throw std::invalid_argument("Invalid solver");
    }
    return error;
}
