/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#include <iostream>

#include "loaders/loader.hpp"
#include "elasticity.hpp"
#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

template<typename MeshType, typename Function, typename Solution>
void test_new_elasticity(MeshType& msh, const Function& load, const Solution& solution, size_t degree)
{
    typedef MeshType                            mesh_type;

    typedef typename mesh_type::scalar_type     scalar_type;

    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    typedef disk::quadrature<mesh_type, cell_type>   cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>   face_quadrature_type;
    typedef disk::scaled_monomial_vector_sg_basis<mesh_type, cell_type>  cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>     div_cell_basis_type;
    typedef disk::scaled_monomial_vector_sg_basis<mesh_type, face_type>  face_basis_type;


    if (degree < 1)
    {
        std::cout << "Only K > 0 for this problem" << std::endl;
        return;
    }


    disk::gradient_reconstruction_nopre<mesh_type,
                                        cell_basis_type,
                                        cell_quadrature_type,
                                        face_basis_type,
                                        face_quadrature_type> gradrec(degree);

    disk::divergence_reconstruction_nopre<mesh_type,
                                          cell_basis_type,
                                          cell_quadrature_type,
                                          face_basis_type,
                                          face_quadrature_type,
                                          div_cell_basis_type,
                                          cell_quadrature_type> divrec(degree);


    disk::diffusion_like_stabilization_nopre<mesh_type,
                                             cell_basis_type,
                                             cell_quadrature_type,
                                             face_basis_type,
                                             face_quadrature_type> stab(degree);

    disk::diffusion_like_static_condensation_nopre<mesh_type,
                                                   cell_basis_type,
                                                   cell_quadrature_type,
                                                   face_basis_type,
                                                   face_quadrature_type> statcond(degree);

    disk::assembler_nopre<mesh_type,
                          face_basis_type,
                          face_quadrature_type> assembler(msh, degree);


    scalar_type mu      = 1.0;
    scalar_type lambda  = 1.0;

    timecounter tc;

    tc.tic();
    for (auto& cl : msh)
    {
        gradrec.compute(msh, cl);
        divrec.compute(msh, cl);
        stab.compute(msh, cl, gradrec.oper);
        auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, load, degree);
        dynamic_matrix<scalar_type> loc = 2 * mu * gradrec.data +
                                          lambda * divrec.data +
                                          2 * mu * stab.data;
        auto sc = statcond.compute(msh, cl, loc, cell_rhs);
        assembler.assemble(msh, cl, sc);
    }

    assembler.impose_boundary_conditions(msh, solution);
    assembler.finalize();
    tc.toc();

    std::cout << "Assembly time: " << tc << " seconds." << std::endl;

    tc.tic();

#ifdef HAVE_SOLVER_WRAPPERS
    agmg_solver<scalar_type> solver;
    dynamic_vector<scalar_type> X = solver.solve(assembler.matrix, assembler.rhs);
#else

#ifdef HAVE_INTEL_MKL
    Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
#else
    Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
#endif
    solver.analyzePattern(assembler.matrix);
    solver.factorize(assembler.matrix);
    dynamic_vector<scalar_type> X = solver.solve(assembler.rhs);
#endif

    tc.toc();
    std::cout << "Solver time: " << tc << " seconds." << std::endl;

    face_basis_type face_basis(degree);
    auto fbs = face_basis.size();

    scalar_type diam = 0.0;

    for (auto& cl : msh)
    {
        diam = std::max(diameter(msh, cl), diam);
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        dynamic_vector<scalar_type> xFs = dynamic_vector<scalar_type>::Zero(num_faces*fbs);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            dynamic_vector<scalar_type> xF = dynamic_vector<scalar_type>::Zero(fbs);
            xF = X.block(face_id * fbs, 0, fbs, 1);
            xFs.block(face_i * fbs, 0, fbs, 1) = xF;
        }

        gradrec.compute(msh, cl);
        divrec.compute(msh, cl);
        stab.compute(msh, cl, gradrec.oper);
        auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, load, degree);
        dynamic_matrix<scalar_type> loc = 2 * mu * gradrec.data +
                                          lambda * divrec.data +
                                          2 * mu * stab.data;

        dynamic_vector<scalar_type> x = statcond.recover(msh, cl, loc, cell_rhs, xFs);
    }

}

template<typename MeshType>
void test_elasticity(MeshType& msh)
{
    typedef MeshType                            mesh_type;

    typedef typename mesh_type::scalar_type     scalar_type;

    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    typedef static_vector<scalar_type, 3> result_type;

    typedef disk::quadrature<mesh_type, face_type>   face_quadrature_type;
    typedef disk::scaled_monomial_vector_sg_basis<mesh_type, cell_type>  cell_basis_type;
    typedef disk::scaled_monomial_vector_sg_basis<mesh_type, face_type>  face_basis_type;

    size_t degree = 1;

    typedef disk::elasticity_template<mesh_type> elasticity;

    auto f = [](const point<scalar_type,3>& p) -> result_type {
        scalar_type fx = M_PI*M_PI*sin(M_PI*p.x());
        scalar_type fy = M_PI*M_PI*sin(M_PI*p.y());
        scalar_type fz = M_PI*M_PI*sin(M_PI*p.z());

        return 2*result_type{fx,fy,fz} + result_type{fx,fy,fz};
    };

    auto sf = [](const point<scalar_type,3>& p) -> result_type {
        scalar_type fx = sin(M_PI*p.x());
        scalar_type fy = sin(M_PI*p.y());
        scalar_type fz = sin(M_PI*p.z());
        return result_type{fx,fy,fz};
    };

    typedef Eigen::Triplet<scalar_type> triplet_type;

    elasticity elast(degree);
    size_t face_basis_size = elast.face_basis_size();

    size_t nunkw = face_basis_size * (msh.faces_size() + msh.boundary_faces_size());
    Eigen::SparseMatrix<scalar_type> A(nunkw, nunkw);
    std::vector<triplet_type> triplets;
    dynamic_vector<scalar_type> rhs = dynamic_vector<scalar_type>::Zero(nunkw);

/*
    for (auto& cl : msh)
    {
        elast.test_operators(msh, cl);
    }
    return;
*/

    timecounter tc;

    tc.tic();

    for (auto& cl : msh)
    {
        auto LC = elast.build_local_contrib(msh, cl, f);
        return;
        /* PROBLEM ASSEMBLY */
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        std::vector<size_t> l2g(face_basis_size * num_faces);
        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            auto face_offset = face_id * face_basis_size;

            auto begin = face_i * face_basis_size;
            auto end = (face_i+1) * face_basis_size;
            for (size_t i = begin; i < end; i++)
                l2g.at(i) = face_offset - (face_i*face_basis_size);

        }

        auto AC = LC.first;
        auto bC = LC.second;

        for (size_t i = 0; i < AC.rows(); i++)
        {
            for (size_t j = 0; j < AC.cols(); j++)
                triplets.push_back( triplet_type( l2g.at(i)+i, l2g.at(j)+j, AC(i,j) ) );

            rhs(l2g.at(i)+i) += bC(i);
        }

    }

    /* BOUNDARY CONDITIONS */
    auto fbs = face_basis_size;
    size_t face_i = 0;
    for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
    {
        auto bfc = *itor;

        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
        if (!eid.first)
            throw std::invalid_argument("This is a bug: face not found");

        auto face_id = eid.second;

        auto face_offset = face_id * fbs;
        auto face_offset_lagrange = (msh.faces_size() + face_i) * fbs;


        face_quadrature_type fq(2*degree);
        auto fqd = fq.integrate(msh, bfc);

        face_basis_type fb(degree);

        dynamic_matrix<scalar_type> MFF = dynamic_matrix<scalar_type>::Zero(fbs, fbs);
        dynamic_vector<scalar_type> rhs_f = dynamic_vector<scalar_type>::Zero(fbs);

        for (auto& qp : fqd)
        {
            auto f_phi = fb.eval_functions(msh, bfc, qp.point());

            for (size_t i = 0; i < f_phi.size(); i++)
            {
                for (size_t j = 0; j < f_phi.size(); j++)
                    MFF(i,j) += qp.weight() * disk::mm_prod(f_phi[i], f_phi[j]);

                rhs_f(i) += qp.weight() * disk::mm_prod(f_phi[i], sf(qp.point()));
            }
        }


        for (size_t i = 0; i < MFF.rows(); i++)
        {
            for (size_t j = 0; j < MFF.cols(); j++)
            {
                triplets.push_back( triplet_type(face_offset + i, face_offset_lagrange + j, MFF(i,j)) );
                triplets.push_back( triplet_type(face_offset_lagrange + j, face_offset + i, MFF(i,j)) );
            }
            rhs(face_offset_lagrange+i) = rhs_f(i);
        }

        face_i++;
    }

    A.setFromTriplets(triplets.begin(), triplets.end());

    tc.toc();

    std::cout << "Assembly time: " << tc << " seconds." << std::endl;

    //saveMarket(A, "A.mtx");

    /* SOLVE */
    tc.tic();

#ifdef HAVE_INTEL_MKL
    Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>> solver;
#else
    Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>> solver;
#endif
    solver.analyzePattern(A);
    solver.factorize(A);
    dynamic_vector<scalar_type> X = solver.solve(rhs);

    tc.toc();

    std::cout << "Solver time: " << tc << " seconds." << std::endl;


    scalar_type diam = 0.0;
    scalar_type err = 0.0;

    for (auto& cl : msh)
    {
        diam = std::max(diameter(msh, cl), diam);
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        dynamic_vector<scalar_type> xFs = dynamic_vector<scalar_type>::Zero(num_faces*fbs);

        for (face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            dynamic_vector<scalar_type> xF = dynamic_vector<scalar_type>::Zero(fbs);
            xF = X.block(face_id * fbs, 0, fbs, 1);
            xFs.block(face_i * fbs, 0, fbs, 1) = xF;
        }

        dynamic_vector<scalar_type> x = elast.recover_full_solution(msh, cl, xFs, f);
        dynamic_vector<scalar_type> rec = elast.high_order_reconstruction(msh, cl, x);


        /*
        auto test_points = make_test_points(msh, cl);
        for (size_t itp = 0; itp < test_points.size(); itp++)
        {
            auto tp = test_points[itp];

            auto pot = elast.evaluate_at_point(msh, cl, rec, tp);
            std::cout << "***" << std::endl;
            std::cout << pot.transpose() << std::endl;
            std::cout << sf(tp).transpose() << std::endl;
        }
        */

        err += elast.compute_cell_error(msh, cl, rec, sf);
    }

    std::cout << "Max diam: " << diam << std::endl;
    std::cout << "L2 error: " << std::sqrt(err) << std::endl;
}

int main(void)
{
    typedef double                                      RealType;
    typedef disk::simplicial_mesh<RealType, 3>           mesh_type;


    mesh_type msh;
    disk::netgen_mesh_loader<RealType, 3> loader;
    if (!loader.read_mesh("/Users/matteo/mroot/matteo/workingcopies/git-mine/hho/hho/hho3d/tests/data/convt02.mesh"))
    {
        std::cout << "Problem loading mesh." << std::endl;
        return 1;
    }
    loader.populate_mesh(msh);

    //test_elasticity(msh);

    typedef typename mesh_type::scalar_type     scalar_type;

    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    typedef static_vector<scalar_type, 3> result_type;

    auto f = [](const point<scalar_type,3>& p) -> result_type {
        scalar_type fx = M_PI*M_PI*sin(M_PI*p.x());
        scalar_type fy = M_PI*M_PI*sin(M_PI*p.y());
        scalar_type fz = M_PI*M_PI*sin(M_PI*p.z());

        return 2*result_type{fx,fy,fz} + result_type{fx,fy,fz};
    };

    auto sf = [](const point<scalar_type,3>& p) -> result_type {
        scalar_type fx = sin(M_PI*p.x());
        scalar_type fy = sin(M_PI*p.y());
        scalar_type fz = sin(M_PI*p.z());
        return result_type{fx,fy,fz};
    };

    test_elasticity(msh);
    test_new_elasticity(msh, f, sf, 1);

}
