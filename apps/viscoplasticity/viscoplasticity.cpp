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
                    const problem_type& problem,
                    const std::string & other_info)
{
    using T = typename Mesh::coordinate_type;
    T tolerance = 1.e-10, Ninf = 10.e+5;
    size_t max_iters = 50000;

    std::string name;
    switch (problem)
    {
        case DRIVEN:
            name = "driven";
            break;
        case COUETTE:
            name = "couette";
            break;
        case POISEUILLE:
            name = "poiseuille";
            break;
        default:
            std::cout << "wrong arguments" << std::endl;
            exit(1);
    }

    std::string info = name + "_k" + tostr(degree) + "_a" + tostr(alpha) + other_info;
    std::ofstream ofs("errors_" + info + ".data");

    if (!ofs.is_open())
        std::cout << "Error opening errors "<<std::endl;

    typename revolution::hho_degree_info hdi(degree, degree);
    augmented_lagrangian_viscoplasticity<Mesh> als(msh, hdi, alpha);

    auto assembler = als.define_problem(msh, problem);
    als.initialize(msh, assembler);

    for(size_t i = 0; i < max_iters; i++)
    {
        als.run_stokes_like(msh, assembler, i);
        als.update_multiplier(msh, assembler);
        auto error = als.compute_errors(msh, assembler, false);

        T cvg_total (0.), cvg_stress(0.), cvg_gamma(0.);
        std::tie(cvg_total, cvg_stress, cvg_gamma) = als.convergence;
        //ofs << std::sqrt(cvg_total)<< " "<< std::sqrt(cvg_stress) <<" ";
        //ofs << std::sqrt(cvg_gamma)<< " "<< error.first << " " <<error.second <<std::endl;

        if(i % 100 == 0)
            std::cout << "  i : "<< i<<"  - " << std::sqrt(cvg_total)<<std::endl;

        assert(std::sqrt(cvg_total) < Ninf);
        if( std::sqrt(cvg_total)  < tolerance)
            break;
        //if(error.first < 10.e-9 && i > 1)
        //    break;
    }
    ofs.close();

    auto final_error = als.compute_errors(msh, assembler, true);
    als.post_processing( msh, assembler, info, problem);

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

    problem_type problem = DRIVEN;

    while ( (ch = getopt(argc, argv, "k:a:dcp")) != -1 )
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
            case 'd':
                problem = DRIVEN;
                break;

            case 'c':
                problem = COUETTE;
                break;

            case 'p':
                problem = POISEUILLE;
                break;
            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    filename = argv[0];

    if (filename == nullptr)
    {
        std::cout << "Please specified a 2D mesh" << std::endl;
        return 1;
    }

    char *word   = nullptr;
    word = argv[1];

    if (word == nullptr)
    {
        std::cout << "no word specified" << std::endl;
        return 1;
    }

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

        run_viscoplasticity(msh, degree, alpha, problem, word);
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

        run_viscoplasticity(msh, degree, alpha, problem, word);
    }

    /* Medit 2d*/
    if (std::regex_match(filename, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        typedef disk::generic_mesh<RealType, 2>  mesh_type;
        mesh_type msh = disk::load_medit_2d_mesh<RealType>(filename);
        run_viscoplasticity(msh, degree, alpha, problem, word);
    }

    //#if 0
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

        run_viscoplasticity(msh, degree, alpha, problem, word);
    }
    //#endif

    return 0;
}
//#endif
