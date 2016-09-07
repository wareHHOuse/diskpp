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

#include "diffusion.hpp"

#include "loaders/loader.hpp"
#include "hho/hho.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

template<typename MeshType, typename Function, typename Solution>
void
test_new_diffusion(MeshType& msh, const Function& load, const Solution& solution, size_t degree)
{
    typedef MeshType                                                    mesh_type;

    typedef typename mesh_type::scalar_type                             scalar_type;

    typedef typename mesh_type::cell                                    cell_type;
    typedef typename mesh_type::face                                    face_type;

    typedef disk::quadrature<mesh_type, cell_type>                       cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>                       face_quadrature_type;

    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>     cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>     face_basis_type;

    typedef disk::diffusion_local_data<mesh_type,
                                      cell_quadrature_type,
                                      cell_basis_type,
                                      face_quadrature_type,
                                      face_basis_type>                  localdata_type;


    localdata_type                                              dld(msh, degree);
    disk::gradient_reconstruction<localdata_type>                gradrec;
    disk::diffusion_like_stabilization<localdata_type>           stab;
    disk::diffusion_like_static_condensation<localdata_type>     statcond;
    disk::assembler<localdata_type>                              assembler(dld);

    cell_basis_type         cb(degree+1);
    cell_quadrature_type    cq(2*degree+2);

    timecounter tc;

    /* ASSEMBLE PROBLEM */
    tc.tic();

    std::cout << dld.num_cell_dofs() << std::endl;
    std::cout << dld.num_face_dofs() << std::endl;

    std::map<std::string, double> timings;

    //scalar_type proj_err = 0.0;
    //std::ofstream ofsd("debug.dat");
    for (auto& cl : msh)
    {
        timecounter tc_detail;

        tc_detail.tic();
        dld.compute(cl, load);
        tc_detail.toc();
        timings["pre"] += tc.to_double();

        tc_detail.tic();
        gradrec.compute(dld);
        tc_detail.toc();
        timings["gr"] += tc.to_double();

        tc_detail.tic();
        stab.compute(dld, gradrec.oper);
        tc_detail.toc();
        timings["stab"] += tc.to_double();

        tc_detail.tic();
        dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
        auto sc = statcond.compute(dld, loc);
        tc_detail.toc();
        timings["sc"] += tc.to_double();

        tc_detail.tic();
        assembler.assemble(dld, sc);
        tc_detail.toc();
        timings["asm"] += tc.to_double();

#if 0
        auto proj = project(dld, cl, load);
        auto gr = gradrec.oper * proj;
        auto gr2 = dynamic_vector<scalar_type>(cb.size());
        gr2.tail(cb.size()-1) = gr;
        gr2(0) = proj(0);


        //std::cout << " dump " << std::endl;
        //std::cout << gradrec.data << std::endl;

        //std::cout << std::endl;
        //std::cout << stab.data << std::endl;

        //std::cout << proj.transpose();
        //std::cout << std::endl;
        //std::cout << gr2.transpose() << std::endl;
        //std::cout << std::endl;

        auto test_points = make_test_points(msh, cl);
        for (size_t itp = 0; itp < test_points.size(); itp++)
        {
            auto tp = test_points[itp];
            auto pot = 0.;

            auto phi = cb.eval_functions(msh, cl, tp);

            for (size_t i = 0; i < cb.size(); i++)
                pot += phi[i] * gr2(i);

            for (size_t i = 0; i < MeshType::dimension; i++)
                ofsd << tp[i] << " ";
            ofsd << pot << std::endl;
        }

        proj_err += compute_L2_error(dld, degree, load, gr2);
#endif
    }

    //dld.dumptimings();
    //ofsd.close();
    //std::cout << sqrt(proj_err) << std::endl;

    for (auto& t : timings)
        std::cout << t.first << ": " << t.second << std::endl;

    assembler.impose_boundary_conditions(dld, solution);
    assembler.finalize();
    tc.toc();
    std::cout << "Assembly time: " << tc << " seconds." << std::endl;

    /* SOLVE */
    tc.tic();

#ifdef HAVE_INTEL_MKL
    Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
#else
    Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
#endif
    solver.analyzePattern(assembler.matrix);
    solver.factorize(assembler.matrix);
    dynamic_vector<scalar_type> X = solver.solve(assembler.rhs);

    tc.toc();
    std::cout << "Solver time: " << tc << " seconds." << std::endl;


    /* POSTPROCESS */
    size_t fbs = dld.num_face_dofs();
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

        dld.compute(cl, load);
        gradrec.compute(dld);
        stab.compute(dld, gradrec.oper);
        dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
        dynamic_vector<scalar_type> x = statcond.recover(dld, loc, xFs);

        dynamic_vector<scalar_type> rec(cb.size());
        rec.tail(cb.size()-1) = gradrec.oper * x;
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

            //std::cout << potk << " " << potkp << " " << potr << std::endl;

            scalar_type diffk = 0.0;
            scalar_type diffkp = 0.0;
            diffk = (potk - potr) * (potk - potr) * qp.weight();
            diffkp = (potkp - potr) * (potkp - potr) * qp.weight();

            err_fun_k += diffk;
            err_fun_kp += diffkp;
        }

        err_dof_k += compute_L2_error(dld, degree, solution, x);
        err_dof_kp += compute_L2_error(dld, degree+1, solution, x);

    }

    ofs.close();

    std::cout << "Mesh diameter: " << diam << std::endl;
    std::cout << "L2-norm error, dof, K:   " << std::sqrt(err_dof_k) << std::endl;
    std::cout << "L2-norm error, fun, K:   " << std::sqrt(err_fun_k) << std::endl;
    std::cout << "L2-norm error, dof, K+1: " << std::sqrt(err_dof_kp) << std::endl;
    std::cout << "L2-norm error, fun, K+1: " << std::sqrt(err_fun_kp) << std::endl;
}

template<typename MeshType, typename Function, typename Solution>
auto
test_diffusion(MeshType& msh, const Function& load, const Solution& solution, size_t degree)
{
    typedef MeshType                            mesh_type;

    typedef typename mesh_type::scalar_type     scalar_type;

    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    typedef disk::quadrature<mesh_type, cell_type>   cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>   face_quadrature_type;

    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>  cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>  face_basis_type;

    typedef disk::diffusion_template<mesh_type, cell_quadrature_type, cell_basis_type, face_quadrature_type, face_basis_type> diffusion;

    typedef Eigen::Triplet<scalar_type> triplet_type;

    diffusion diff(degree);
    size_t face_basis_size = diff.face_basis_size();

    size_t nunkw = face_basis_size * (msh.faces_size() + msh.boundary_faces_size());
    Eigen::SparseMatrix<scalar_type> A(nunkw, nunkw);
    std::vector<triplet_type> triplets;
    dynamic_vector<scalar_type> rhs = dynamic_vector<scalar_type>::Zero(nunkw);

    /*
    for (auto& cl : msh)
    {
        std::cout << cl << std::endl;
        diff.test_operators(msh, cl);
    }
    return;
    */

    timecounter tc;

    tc.tic();

    for (auto& cl : msh)
    {
        auto LC = diff.build_local_contrib(msh, cl, load);

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
                    MFF(i,j) += qp.weight() * f_phi[i] * f_phi[j];

                rhs_f(i) += qp.weight() * f_phi[i] * solution(qp.point());
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
    std::ofstream ofs("plot.dat");
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

        dynamic_vector<scalar_type> x = diff.recover_full_solution(msh, cl, xFs, load);
        dynamic_vector<scalar_type> rec = diff.high_order_reconstruction(msh, cl, x);

        auto test_points = make_test_points(msh, cl);
        for (size_t itp = 0; itp < test_points.size(); itp++)
        {
            auto tp = test_points[itp];

            auto pot = diff.evaluate_at_point(msh, cl, rec, tp);

            for (size_t i = 0; i < MeshType::dimension; i++)
                ofs << tp[i] << " ";
            ofs << pot << " " << std::abs(pot - solution(tp)) << std::endl;
        }

        err += diff.compute_cell_error(msh, cl, rec, solution);
    }

    ofs.close();

    std::cout << "Mesh diameter: " << diam << std::endl;
    std::cout << "L2-norm error: " << std::sqrt(err) << std::endl;

    return std::make_tuple(diam, err);
}

#if 0
int run_convtest_3d(void)
{
    std::vector<std::string> meshes_3d;
    meshes_3d.push_back("../../../meshes/convt01.mesh");
    meshes_3d.push_back("../../../meshes/convt02.mesh");
    meshes_3d.push_back("../../../meshes/convt03.mesh");
    meshes_3d.push_back("../../../meshes/convt04.mesh");
}

int run_convtest_2d(void)
{
    using RealType = double;
    std::vector<std::string> meshes_2d;
    meshes_2d.push_back("../../../meshes/mesh2_1.typ1");
    meshes_2d.push_back("../../../meshes/mesh2_2.typ1");
    meshes_2d.push_back("../../../meshes/mesh2_3.typ1");
    meshes_2d.push_back("../../../meshes/mesh2_4.typ1");

    std::ofstream ofs("conv2d.txt");

    for (auto& mn : meshes_2d)
    {
        typedef hho::generic_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        hho::fvca5_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(mn))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return M_PI * M_PI * sin(p.y() * M_PI);
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return sin(p.y() * M_PI);
        };

        for (size_t i = 0; i < 4; i++)
        {
            auto t1 = test_diffusion(msh, f, sf, i);
            auto t2 = test_new_diffusion(msh, f, sf, i);

            if (i == 0)
                ofs << std::get<0>(t1) << " ";
            ofs << std::get<1>(t1) << " ";
            ofs << std::get<1>(t2) << " ";
            ofs << std::get<2>(t2) << " ";
        }

        ofs << std::endl;
    }

    ofs.close();
}
#endif

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
        test_new_diffusion(msh, f, sf, degree);

        return 0;
    }

    filename = argv[0];

    if (std::regex_match (filename, std::regex(".*\\.typ1$") ))
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
        test_new_diffusion(msh, f, sf, degree);
    }

    if (std::regex_match (filename, std::regex(".*\\.mesh$") ))
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
        test_new_diffusion(msh, f, sf, degree);
    }

    return 0;
}

/*
int main(void)
{
    typedef double                                      RealType;
    typedef hho::simplicial_mesh<RealType, 3>           mesh_type;


    mesh_type msh;
    hho::netgen_mesh_loader<RealType, 3> loader;
    if (!loader.read_mesh("../../../meshes/convt02.mesh"))
    {
        std::cout << "Problem loading mesh." << std::endl;
        return 1;
    }
    loader.populate_mesh(msh);

    test_diffusion(msh);
}
*/

/*
int main(void)
{
    typedef double                                      RealType;
    typedef hho::generic_mesh<RealType, 2>              mesh_type;


    mesh_type msh;
    hho::fvca5_mesh_loader<RealType, 2> loader;
    if (!loader.read_mesh("../../../meshes/mesh2_4.typ1"))
    {
        std::cout << "Problem loading mesh." << std::endl;
        return 1;
    }
    loader.populate_mesh(msh);

    test_diffusion(msh);
}
*/
