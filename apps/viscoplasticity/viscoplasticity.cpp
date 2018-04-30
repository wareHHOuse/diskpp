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


#include "viscoplasticity_solver.hpp"
template<typename Mesh>
auto
run_viscoplasticity(const Mesh& msh, size_t degree,
                    const typename Mesh::scalar_type & alpha,
                    std::ofstream & ofs)
{
    using T = typename Mesh::coordinate_type;
    T tolerance = 1.e-8, Ninf = 10.e+5;
    size_t max_iters = 50000;

    typename revolution::hho_degree_info hdi(degree, degree);
    augmented_lagrangian_viscoplasticity<Mesh> als(msh, hdi, alpha);
    auto assembler = als.define_problem(msh);
    als.initialize(msh, assembler);

    T error_old = 0.;

    for(size_t i = 0; i < max_iters; i++)
    {
        als.run_stokes_like(msh, assembler, i);
        als.update_multiplier(msh, assembler);
        auto error = als.compute_errors(msh, assembler, false);

        T cvg_total (0.), cvg_stress(0.), cvg_gamma(0.);
        std::tie(cvg_total, cvg_stress, cvg_gamma) = als.convergence;

        ofs << std::sqrt(cvg_total)<< " "<< std::sqrt(cvg_stress) <<" ";
        ofs << std::sqrt(cvg_gamma)<< " "<< error.first << " " <<error.second <<std::endl;
        assert(std::sqrt(cvg_total) < Ninf);
        if( std::sqrt(cvg_total)  < tolerance)
            break;
        //if(error.first < 10.e-9 && i > 1)
        //    break;
    }

    auto final_error = als.compute_errors(msh, assembler, true);

    //post-processing
    auto dim = Mesh::dimension;
    auto rbs   = revolution::vector_basis_size(hdi.reconstruction_degree(), dim, dim);
    auto cbs   = revolution::vector_basis_size(hdi.cell_degree(), dim, dim);
    Matrix< T, Dynamic, 1> cell_sol(cbs * msh.cells_size());
    Matrix< T, Dynamic, 1> cell_rec_sol(rbs * msh.cells_size());
    Matrix< T, Dynamic, 1> press_vec(msh.cells_size());

    size_t cl_count = 0;
    for(auto cl : msh)
    {
        auto gr  = revolution::make_hho_stokes(msh, cl, hdi, als.use_sym_grad);
        auto cell_ofs = revolution::priv::offset(msh, cl);
        Matrix<T, Dynamic, 1> svel =  assembler.take_velocity(msh, cl, als.sol);
        cell_rec_sol.block(cell_ofs * rbs + dim, 0, rbs - dim, 1) = gr.first * svel;
        cell_rec_sol.block(cell_ofs * rbs, 0, dim, 1) = svel.block(0,0, dim, 1);
        cell_sol.block(cell_ofs * cbs, 0, cbs, 1) = svel.block(0,0, cbs, 1);

        //this is only for k = 0, since there is only one velocity;
        auto bar = barycenter(msh, cl);
        Matrix<T, Dynamic, 1> spress =  assembler.take_pressure(msh, cl, als.sol);
        auto pb  = revolution::make_scalar_monomial_basis(msh, cl, hdi.face_degree());
        auto p_phi = pb.eval_functions(bar);
        press_vec(cl_count++) =  p_phi.dot(spress);
    }

    typedef point<T,2>             point_type;

    std::string info = "_k" + tostr(hdi.face_degree()) + "_a"+ tostr(alpha);
    compute_discontinuous_velocity( msh, cell_sol, hdi, "driven_veclocity2d_quad_B0_h01" + info +".msh");

    std::pair<point_type, point_type> p_x, p_y;

    switch(als.problem_num)  //Do this with enum "driven")
    {
        case 1:
            p_x = std::make_pair(point_type({0.0, 0.5}), point_type({1.0, 0.5}));
            p_y = std::make_pair(point_type({0.5, 0.0}), point_type({0.5, 1.0}));
            plot_over_line(msh, p_x, cell_rec_sol, hdi.reconstruction_degree(), "plot_over_line_x_alg_quad_B0_h01.data");
            plot_over_line(msh, p_y, cell_rec_sol, hdi.reconstruction_degree(), "plot_over_line_y_alg_quad_B0_h01.data");
            break;
        case 2:
            p_x = std::make_pair(point_type({0., 0.}), point_type({0.866, 0.5}));
            plot_over_line(msh, p_x, cell_sol, hdi.cell_degree(), "plot_over_line_x_alg"+ info +".data");
            break;
    }
    save_coords(msh, "Coords.data");
    save_data(press_vec, "pressure.data");

    //quiver( msh, als.sol, assembler, hdi);
    return final_error;
}

#if 0
void convergence_test_typ1(void)
{
    using T = double;
    bool use_sym_grad = false;
    std::vector<std::string> meshfiles;

    //meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_1.typ1");
    //meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_2.typ1");
    //meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_3.typ1");
    //meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_4.typ1");
    //meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_5.typ1");
    //meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_6.typ1");

    /*
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_1.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_2.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_3.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_4.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_5.typ1");
    */
    /*
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_1.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_2.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_3.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_4.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_5.typ1");
    */
    std::cout << "                   velocity H1-error";
    std::cout << "    -     pressure L2-error "<< std::endl;

    for (size_t k = 2; k < 3; k++)
    {
        std::cout << "DEGREE " << k << std::endl;

        std::ofstream ofs("errors_k" + tostr(k) + ".data");
        if (!ofs.is_open())
            std::cout << "Error opening errors "<<std::endl;

        std::vector<T> mesh_hs;
        std::vector<std::pair<T,T>> errors;

        //for (size_t i = 0; i < 1; i++)
        for (size_t i = 0; i < meshfiles.size(); i++)
        {
            typedef disk::generic_mesh<T, 2>  mesh_type;

            std::cout << " Mesh : "<< i << std::endl;
            mesh_type msh;
            disk::fvca5_mesh_loader<T, 2> loader;
            if (!loader.read_mesh(meshfiles.at(i)))
            {
                std::cout << "Problem loading mesh." << std::endl;
                continue;
            }
            loader.populate_mesh(msh);

            auto error = run_alg_stokes(msh, k, ofs, use_sym_grad);

            mesh_hs.push_back( disk::mesh_h(msh) );
            errors.push_back(error);
            ofs << " " << std::endl;
        }
        ofs.close();

        for (size_t i = 0; i < mesh_hs.size(); i++)
        {
            if (i == 0)
            {
                std::cout << "    ";
                std::cout << std::scientific << std::setprecision(4) << mesh_hs.at(i) << "    ";
                std::cout << std::scientific << std::setprecision(4) << errors.at(i).first;
                std::cout << "     -- " << "          ";
                std::cout << std::scientific << std::setprecision(4) << errors.at(i).second;
                std::cout << "     -- " << std::endl;
            }
            else
            {
                auto rate = std::log( errors.at(i).first/errors.at(i-1).first ) /
                            std::log( mesh_hs.at(i)/mesh_hs.at(i-1) );
                std::cout << "    ";
                std::cout << std::scientific  << std::setprecision(4) << mesh_hs.at(i) << "    ";
                std::cout << std::scientific  << std::setprecision(4) << errors.at(i).first << "    ";
                std::cout << std::fixed<< std::setprecision(2) << rate << "          ";

                auto pres_rate = std::log( errors.at(i).second/errors.at(i-1).second ) /
                            std::log( mesh_hs.at(i)/mesh_hs.at(i-1) );
                std::cout << std::scientific  << std::setprecision(4) << errors.at(i).second << "    ";
                std::cout << std::fixed << std::setprecision(2) << pres_rate << std::endl;
            }
        }
    }
}

int main(void)
{
    convergence_test_typ1();
}
#endif

//#if 0

int main(int argc, char **argv)
{
    using RealType = double;

    char    *filename       = nullptr;
    int ch;
    size_t degree = 0;
    RealType alpha = 1.;

    while ( (ch = getopt(argc, argv, "k:a:")) != -1 )
    {
        switch(ch)
        {
            case 'k':
                degree = atoi(optarg);
                break;

            case 'a':
                alpha = atof(optarg);
                if (alpha <= 0)
                {
                    std::cout << "alpha must be >=0. Falling back to 1." << std::endl;
                    alpha = 1.;
                }
                break;
            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    if (argc != 1)
    {
        std::cout << "Please specify a 2D mesh" << std::endl;

        return 0;
    }

    filename = argv[0];

    std::ofstream ofs("errors_k" + tostr(degree) + "_a" + tostr(alpha)+"_quad_B0_h01.data");
    if (!ofs.is_open())
        std::cout << "Error opening errors "<<std::endl;

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

        run_viscoplasticity(msh, degree, alpha, ofs);
    }

    if (std::regex_match(filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;

        typedef disk::simplicial_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::netgen_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        run_viscoplasticity(msh, degree, alpha, ofs);
    }

    /* Medit 2d*/
    if (std::regex_match(filename, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        typedef disk::generic_mesh<RealType, 2>  mesh_type;
        mesh_type msh = disk::load_medit_2d_mesh<RealType>(filename);
        run_viscoplasticity(msh, degree, alpha, ofs);
    }

    #if 0
    if (std::regex_match(filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: Cartesian 2D" << std::endl;

        typedef disk::cartesian_mesh<RealType, 2>   mesh_type;

        mesh_type msh;
        disk::cartesian_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        test_bases(msh);
    }
    #endif
    ofs.close();
    return 0;
}
//#endif
