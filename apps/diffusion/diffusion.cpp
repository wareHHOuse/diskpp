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
#include <regex>
#include <unistd.h>

#include <map>

#include "config.h"

#ifdef HAVE_SOLVER_WRAPPERS
    #include "agmg/agmg.hpp"
#endif

#include "diffusion.hpp"

#include "loaders/loader.hpp"
#include "hho/hho.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

template<typename MeshType, typename Function, typename Solution>
void
test_diffusion(MeshType& msh,               /* handle to the mesh */
               const Function& load,        /* rhs */
               const Solution& solution,    /* solution of the problem */
               size_t degree)               /* degree of the method */
{
    typedef MeshType                                   mesh_type;
    typedef typename mesh_type::scalar_type            scalar_type;
    typedef typename mesh_type::cell                   cell_type;
    typedef typename mesh_type::face                   face_type;

    typedef disk::quadrature<mesh_type, cell_type>      cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>      face_quadrature_type;

    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;


    disk::gradient_reconstruction_nopre<mesh_type,
                                        cell_basis_type,
                                        cell_quadrature_type,
                                        face_basis_type,
                                        face_quadrature_type> gradrec_nopre(degree);


    disk::diffusion_like_stabilization_nopre<mesh_type,
                                             cell_basis_type,
                                             cell_quadrature_type,
                                             face_basis_type,
                                             face_quadrature_type> stab_nopre(degree);

    disk::diffusion_like_static_condensation_nopre<mesh_type,
                                                   cell_basis_type,
                                                   cell_quadrature_type,
                                                   face_basis_type,
                                                   face_quadrature_type> statcond_nopre(degree);

    disk::assembler_nopre<mesh_type,
                          face_basis_type,
                          face_quadrature_type> assembler_nopre(msh, degree);


    cell_basis_type         cb(degree+1);
    cell_quadrature_type    cq(2*degree+2);

    timecounter tc;
    std::map<std::string, double> timings;

    /* ASSEMBLE PROBLEM */
    std::cout << "Assembling" << std::endl;
    tc.tic();
    for (auto& cl : msh)
    {
        timecounter tc_detail;

        tc_detail.tic();
        gradrec_nopre.compute(msh, cl);
        tc_detail.toc();
        timings["gr_np"] += tc_detail.to_double();

        tc_detail.tic();
        stab_nopre.compute(msh, cl, gradrec_nopre.oper);
        tc_detail.toc();
        timings["stab_np"] += tc_detail.to_double();

        tc_detail.tic();
        auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, load, degree);
        dynamic_matrix<scalar_type> loc = gradrec_nopre.data + stab_nopre.data;
        auto scnp = statcond_nopre.compute(msh, cl, loc, cell_rhs);
        tc_detail.toc();
        timings["sc_np"] += tc_detail.to_double();

        tc_detail.tic();
        assembler_nopre.assemble(msh, cl, scnp);
        tc_detail.toc();
        timings["asm_np"] += tc_detail.to_double();

    }

    for (auto& t : timings)
        std::cout << t.first << ": " << t.second << std::endl;

    assembler_nopre.impose_boundary_conditions(msh, solution);
    assembler_nopre.finalize();
    tc.toc();
    std::cout << "Assembly time: " << tc << " seconds." << std::endl;

    /* SOLVE */
    tc.tic();

#ifdef HAVE_SOLVER_WRAPPERS
    agmg_solver<scalar_type> solver;
    dynamic_vector<scalar_type> X = solver.solve(assembler_nopre.matrix, assembler_nopre.rhs);
#else

#ifdef HAVE_INTEL_MKL
    Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
#else
    Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
#endif
    solver.analyzePattern(assembler_nopre.matrix);
    solver.factorize(assembler_nopre.matrix);
    dynamic_vector<scalar_type> X = solver.solve(assembler_nopre.rhs);
#endif

    tc.toc();
    std::cout << "Solver time: " << tc << " seconds." << std::endl;

    face_basis_type face_basis(degree);
    size_t fbs = face_basis.size();

    scalar_type diam = 0.0;
    scalar_type err_dof_k = 0.0;
    scalar_type err_fun_k = 0.0;
    scalar_type err_dof_kp = 0.0;
    scalar_type err_fun_kp = 0.0;
    std::ofstream ofs("plotnew.dat");
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

        gradrec_nopre.compute(msh, cl);
        stab_nopre.compute(msh, cl, gradrec_nopre.oper);
        dynamic_matrix<scalar_type> loc = gradrec_nopre.data + stab_nopre.data;
        auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, load, degree);
        dynamic_vector<scalar_type> x = statcond_nopre.recover(msh, cl, loc, cell_rhs, xFs);

        dynamic_vector<scalar_type> rec(cb.size());
        rec.tail(cb.size()-1) = gradrec_nopre.oper * x;
        rec(0) = x(0);

        auto test_points = make_test_points(msh, cl);
        for (size_t itp = 0; itp < test_points.size(); itp++)
        //for (auto& qp : qps)
        {
            auto tp = test_points[itp];
            //auto tp = qp.point();
            auto pot = 0.;

            auto phi = cb.eval_functions(msh, cl, tp);

            for (size_t i = 0; i < cb.size(); i++)
                pot += phi[i] * rec(i);

            for (size_t i = 0; i < MeshType::dimension; i++)
                ofs << tp[i] << " ";
            ofs << pot << " " << std::abs(pot - solution(tp)) << std::endl;
        }

        auto qps = cq.integrate(msh, cl);

        for (auto& qp : qps)
        {
            auto phi = cb.eval_functions(msh, cl, qp.point());

            scalar_type potk = 0.0;
            for (size_t i = 0; i < cb.range(0,degree).size(); i++)
                potk += phi[i] * rec(i);

            scalar_type potkp = potk;
            for (size_t i = cb.range(0,degree).size(); i < cb.size(); i++)
                potkp += phi[i] * rec(i);

            auto potr = solution(qp.point());

            scalar_type diffk = 0.0;
            scalar_type diffkp = 0.0;
            diffk = (potk - potr) * (potk - potr) * qp.weight();
            diffkp = (potkp - potr) * (potkp - potr) * qp.weight();

            err_fun_k += diffk;
            err_fun_kp += diffkp;
        }

        //err_dof_k += compute_L2_error(dld, degree, solution, x);
        //err_dof_kp += compute_L2_error(dld, degree+1, solution, x);

    }

    ofs.close();

    std::cout << "Mesh diameter: " << diam << std::endl;
    std::cout << "L2-norm error, dof, K:   " << std::sqrt(err_dof_k) << std::endl;
    std::cout << "L2-norm error, fun, K:   " << std::sqrt(err_fun_k) << std::endl;
    std::cout << "L2-norm error, dof, K+1: " << std::sqrt(err_dof_kp) << std::endl;
    std::cout << "L2-norm error, fun, K+1: " << std::sqrt(err_fun_kp) << std::endl;
}

int main(int argc, char **argv)
{
    using RealType = double;

    char    *filename   = nullptr;
    int     degree      = 1;
    int     elems_1d    = 8;
    int ch;

    while ( (ch = getopt(argc, argv, "k:n:")) != -1 )
    {
        switch(ch)
        {
            case 'k':
                degree = atoi(optarg);
                if (degree < 0)
                {
                    std::cout << "Degree must be positive. Falling back to 1." << std::endl;
                    degree = 1;
                }
                break;

            case 'n':
                elems_1d = atoi(optarg);
                if (elems_1d < 0)
                {
                    std::cout << "Num of elems must be positive. Falling back to 8." << std::endl;
                    elems_1d = 8;
                }
                break;

            case 'h':
            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;


    if (argc == 0)
    {
        std::cout << "Running 1D test simulation" << std::endl;

        typedef disk::generic_mesh<RealType, 1>  mesh_type;

        mesh_type msh;
        disk::uniform_mesh_loader<RealType, 1> loader(0,1,elems_1d);
        loader.populate_mesh(msh);

        auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            //return 2.;
            return M_PI * M_PI * sin(p.x() * M_PI);
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            //return - p.x() * p.x();
            return sin(p.x() * M_PI);
        };

        //test_diffusion(msh, f, sf, degree);
        test_diffusion(msh, f, sf, degree);

        return 0;
    }

    filename = argv[0];

    if (std::regex_match(filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;

        typedef disk::generic_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::fvca5_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return M_PI * M_PI * sin(p.x() * M_PI);
            //return 6. * p.x();
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return sin(p.x() * M_PI);
            //return - p.x() * p.x() * p.x();
        };

        //test_diffusion(msh, f, sf, degree);
        test_diffusion(msh, f, sf, degree);
    }

    if (std::regex_match(filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;

        typedef disk::simplicial_mesh<RealType, 3>   mesh_type;

        mesh_type msh;
        disk::netgen_mesh_loader<RealType, 3> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return M_PI * M_PI * sin(p.x() * M_PI);
            //return 12*p.y()*p.y();
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return sin(p.x() * M_PI);
            //return -p.y() * p.y() * p.y() * p.y();
        };

        //test_diffusion(msh, f, sf, degree);
        test_diffusion(msh, f, sf, degree);
    }

    return 0;
}
