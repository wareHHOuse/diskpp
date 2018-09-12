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

        return 2.0 * M_PI * M_PI * cos_px * cos_py;
        //return 2.0 * M_PI * M_PI * sin_px * sin_py;
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

        return cos_px * cos_py;
        //return sin_px * sin_py;
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


template<typename Mesh>
auto
run_hho_diffusion_nitsche_cells_full(const Mesh& msh,
    const algorithm_parameters<typename Mesh::scalar_type>& ap,
    const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
    const typename Mesh::scalar_type& eta)
{
    using T =  typename Mesh::scalar_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    hho_degree_info hdi(ap.degree + 1, ap.degree);

    auto is_contact_vector = make_is_contact_vector(msh, bnd);

    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension - 1);

    auto rhs_fun = make_rhs_function(msh);
    auto sol_fun = make_solution_function(msh);
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

            matrix_type Aconsist   = make_hho_consist_mix_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, eta, bnd);
            auto ntz   = make_hho_nitsche_mix_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, sol_fun, eta, bnd );
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
            //auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);

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

    std::ofstream ofs("solfull.dat");

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

            matrix_type Aconsist   = make_hho_consist_mix_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, eta, bnd);
            auto ntz   = make_hho_nitsche_mix_par(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, sol_fun, eta, bnd );
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
            //auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);

            vector_type rhs = make_rhs(msh, cl, cb, rhs_fun);//, hdi.cell_degree());
            A = gr.second + stab;

            vector_type sol_faces = assembler.take_local_data(msh, cl, sol);
            fullsol = assembler.take_local_data(msh, cl, sol);
        }

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

    auto itor = std::max_element(inf_errors.begin(), inf_errors.end());
    size_t distance = std::distance(inf_errors.begin(), itor);
    auto Linf_error = *std::next(inf_errors.begin(), distance);


    ofs.close();

    H1_error = std::sqrt(H1_error);
    L2_error = std::sqrt(L2_error);

    auto error = std::make_tuple(H1_error, L2_error, Linf_error);
    //std::cout << std::sqrt(error) << std::endl;

    return error;
}


template<typename Mesh, typename T>
auto
run_diffusion_solver(const Mesh& msh, const algorithm_parameters<T>& ap,
                     const T& eta)
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

    /*--------------------------------------------------------------------------
    *  Check boundary labels for the unitary square domain
    *          Netgen     __1__          Medit     __1__
    *                4   |     | 2                |     |
    *                    |_____|                  |_____|
    *                       3                        2
    *-------------------------------------------------------------------------*/
    //test 1
    //#if 0
    bnd.addDirichletBC(disk::mechanics::DIRICHLET,1,sol_fun); //TOP
    bnd.addContactBC(disk::mechanics::SIGNORINI, 3); //
    bnd.addNeumannBC(disk::mechanics::NEUMANN, 4, neu_left); //
    bnd.addNeumannBC(disk::mechanics::NEUMANN, 2, neu_right); //Bottom
    //bnd.addContactBC(disk::mechanics::SIGNORINI,3); //BOTTOM
    //#endif

    //test 2
    #if 0
    bnd.addContactBC(disk::mechanics::SIGNORINI, 1); //
    bnd.addContactBC(disk::mechanics::SIGNORINI, 2); //
    bnd.addContactBC(disk::mechanics::SIGNORINI, 3); //
    bnd.addContactBC(disk::mechanics::SIGNORINI, 4); //
    #endif



    std::tuple<T, T, T> error;
    switch (ap.solver)
    {
        //case EVAL_IN_CELLS:
        //    error = run_hho_diffusion_nitsche_cells(msh, ap, bnd, eta);
        //    break;
        case EVAL_IN_CELLS_FULL: //Temporal name
            error = run_hho_diffusion_nitsche_cells_full(msh, ap, bnd, eta);
            break;
        #if 0
        case EVAL_ON_FACES:
            error = run_hho_diffusion_nitsche_faces(msh, ap);
            break;
        case EVAL_WITH_PARAMETER:
            error = run_hho_diffusion_nitsche_par(msh, ap, bnd, eta);
            break;
        case EVAL_IN_CELLS_AS_FACES:
            std::cout << "Not valid in this case. Choose faces (f) or cells( c)" << std::endl;
            break;
        #endif
        default:
            throw std::invalid_argument("Invalid solver");
    }
    return error;
}
