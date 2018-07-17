/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
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

#include <iostream>
#include <iomanip>
#include <regex>

#include <unistd.h>

#include "revolution/bases"
#include "revolution/quadratures"
#include "revolution/methods/hho"

#include "core/loaders/loader.hpp"

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"
#include "output/postMesh.hpp"

#include "output/silo.hpp"
#include "solvers/solver.hpp"
#include "viscoplasticity_utils.hpp"

enum problem_type
{
    DRIVEN
};

template<typename Mesh>
class augmented_lagrangian_viscoplasticity
{
    typedef Mesh mesh_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;
    typedef typename mesh_type::scalar_type T;

    typedef disk::mechanics::BoundaryConditions<mesh_type> boundary_type;

    using point_type = typename mesh_type::point_type;

    typedef Matrix<T, Mesh::dimension, Mesh::dimension>         tensor_type;
    typedef Matrix<T, Mesh::dimension, 1>                       vector2d_type;

    typedef Matrix<T, Dynamic, Dynamic>                 matrix_type;
    typedef Matrix<T, Dynamic, 1>                       vector_type;

    typedef std::function<vector2d_type (const point_type &)>   vector_funtion_type;
    typedef std::function<T   (const point_type &)>             scalar_funtion_type;

    hho_degree_info di;
    T             viscosity;
    T             alpha;
    T             yield;
    size_t        cbs, fbs, pbs, sbs;

public:
    vector_type             sol_old;
    std::tuple<T, T, T>     convergence;
    boundary_type           bnd;
    tensors_at_quad_pts_utils   tsr_utils;

    augmented_lagrangian_viscoplasticity(const Mesh& msh,
                            const typename hho_degree_info & hdi,
                            const T& alpha_ext):
                            di(hdi), alpha(alpha_ext)
    {
        viscosity = 1.;
        T f = 1;
        T Lref = 1.;
        T Bn  =  2;
        yield =  Bn * f * Lref; // * viscosity;// * omegaExt; //* f * Lref;

        const auto dim =  Mesh::dimension;

        cbs = revolution::vector_basis_size(di.cell_degree(), dim, dim);
        fbs = revolution::vector_basis_size(di.face_degree(), dim - 1, dim);
        pbs = revolution::scalar_basis_size(di.face_degree(), dim);
        sbs = revolution::sym_matrix_basis_size(di.face_degree(), dim, dim);

        size_t quad_degree = 2. * di.face_degree();
        tsr_utils = tensors_at_quad_pts_utils<mesh_type>(msh, quad_degree);
    };

    auto
    define_problem(const mesh_type& msh, const problem_type& problem )
    {
        bnd = boundary_type(msh);

        auto wall = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
            return Matrix<T, Mesh::dimension, 1>::Zero();
        };
        auto movingWall = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
            return Matrix<T, Mesh::dimension, 1>{1,0};
        };

        switch (problem)
		{
            case DRIVEN:

                rhs_fun  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    return Matrix<T, Mesh::dimension, 1>::Zero();
                };

                bnd.addDirichletBC(0, 1, movingWall);
                bnd.addDirichletBC(0, 2, wall);
                bnd.addDirichletBC(0, 3, wall);
                bnd.addDirichletBC(0, 4, wall);
               break;
            default:
                throw std::invalid_argument("Invalid problem");
        }

    }


    template<typename Assembler>
    auto
    make_global_matrix(const mesh_type& msh, Assembler& assembler)
    {
        assembler.initialize_lhs();

        for (auto cl : msh)
        {
            auto G  = revolution::make_hlow_stokes(msh, cl, di, true);
            auto gr = revolution::make_hho_stokes(msh, cl, di, true);
            Matrix<T, Dynamic, Dynamic> stab;
            stab = make_hho_vector_stabilization(msh, cl, gr.first, di);
            auto dr = make_hho_divergence_reconstruction_stokes_rhs(msh, cl, di);

            Matrix<T, Dynamic, Dynamic> A = 2.*(alpha * G.second + viscosity * stab);

            assembler.assemble_lhs(msh, cl, A, -dr);
        }

        assembler.finalize_lhs();

        return;
    }

    run_alg(const mesh_type& msh)
    {
        auto assembler = revolution::make_stokes_assembler_alg(msh, di, bnd);
        auto systsz = assembler.global_system_size();
        sol_old = vector_type::Zero(systsz);
        dynamic_vector<T> sol =  dynamic_vector<T>::Zero(systsz);

        auto num_total_quad_dofs = tsr_utils.num_total_quad_points(msh);

        multiplier   = matrix_type::Zero(sbs, num_total_quads);
        auxiliar     = matrix_type::Zero(sbs, num_total_quads);
        auxiliar_old = matrix_type::Zero(sbs, num_total_quads);

        for(size_t i = 0; i < max_iters; i++)
        {
            make_global_rhs(msh,assembler);
            if(iter == 0)
                make_global_matrix(msh, assembler);

            //dump_sparse_matrix(assembler.LHS, "stokes.txt");
            size_t systsz = assembler.LHS.rows();
            size_t nnz = assembler.LHS.nonZeros();

            dynamic_vector<T>::Zero sol = dynamic_vector<T>::Zero(systsz);
            disk::solvers::pardiso_params<T> pparams;
            pparams.report_factorization_Mflops = true;
            mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

            //---------------------------------------------------------------------
            T cvg_total (0.), cvg_stress(0.), cvg_gamma(0.);
            std::tie(cvg_total, cvg_stress, cvg_gamma) = als.convergence;

            if(i % 1000 == 0)
                std::cout << "  i : "<< i<<"  - " << std::sqrt(cvg_total)<<std::endl;

            assert(std::sqrt(cvg_total) < Ninf);
            if( std::sqrt(cvg_total)  < tolerance)
            {
                std::cout << "  i : "<< i<<"  - " << std::sqrt(cvg_total)<<std::endl;
                break;
            }

            sol_old = sol;
        }

    }
}
