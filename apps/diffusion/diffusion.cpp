/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016, 2017 - matteo.cicuttin@enpc.fr
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
#include <sstream>
#include <iomanip>

#include <map>

#include "colormanip.h"

#include "../../config.h"

#ifdef HAVE_SOLVER_WRAPPERS
    #include "agmg/agmg.hpp"
#endif

#include "loaders/loader.hpp"
#include "hho/hho.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "diffusion_solver.hpp"

struct run_params
{
    size_t  degree;
    int     l;
    bool    verbose;
};

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
void
run_diffusion_solver(const Mesh<T, 1, Storage>& msh, run_params& rp)
{
    typedef Mesh<T, 1, Storage> mesh_type;

    auto load = [](const point<T, 1>& p) -> auto {
        return M_PI * M_PI * sin(p.x() * M_PI);
    };

    auto solution = [](const point<T, 1>& p) -> auto {
        return sin(p.x() * M_PI);
    };

    diffusion_solver<mesh_type> dp(msh, rp.degree);
    dp.verbose(rp.verbose);

    dp.assemble(load, solution);
    dp.solve();
    dp.postprocess(load);
    //dp.plot_solution("plot.dat");
    std::cout << dp.compute_l2_error(solution) << std::endl;
}

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
void
run_diffusion_solver(const Mesh<T, 2, Storage>& msh, run_params& rp)
{
    typedef Mesh<T, 2, Storage> mesh_type;

    auto load = [](const point<T, 2>& p) -> auto {

        //return 2.0 * M_PI * M_PI * sin(p.x() * M_PI) * sin(p.y() * M_PI);
        return sin(p.x()) * sin(p.y());

    };

    auto solution = [](const point<T, 2>& p) -> auto {
        return sin(p.x() * M_PI) * sin(p.y() * M_PI);
    };

    diffusion_solver<mesh_type> dp(msh, rp.degree);
    dp.verbose(rp.verbose);

    std::cout << "ASM" << std::endl;
    dp.assemble(load, solution);
    std::cout << "solve" << std::endl;
    dp.solve();
    std::cout << "post" << std::endl;
    dp.postprocess(load);
    //dp.plot_solution("plot.dat");
    std::cout << dp.compute_l2_error(solution) << std::endl;
}

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
void
run_diffusion_solver(const Mesh<T, 3, Storage>& msh, run_params& rp)
{
    typedef Mesh<T, 3, Storage> mesh_type;

    auto load = [](const point<T, 3>& p) -> auto {
        return 3.0 * M_PI * M_PI * sin(p.x() * M_PI) * sin(p.y() * M_PI) * sin(p.z() * M_PI);
    };

    auto solution = [](const point<T, 3>& p) -> auto {
        return sin(p.x() * M_PI) * sin(p.y() * M_PI) * sin(p.z() * M_PI);
    };

    diffusion_solver<mesh_type> dp(msh, rp.degree);
    dp.verbose(rp.verbose);

    dp.assemble(load, solution);
    dp.solve();
    dp.postprocess(load);
    dp.plot_solution("plot.dat");
    std::cout << dp.compute_l2_error(solution) << std::endl;
}

int main(int argc, char **argv)
{
    using RealType = double;

    char    *mesh_filename  = nullptr;
    char    *plot_filename  = nullptr;
    int     degree          = 1;
    int     l               = 0;
    int     elems_1d        = 8;

    run_params rp;
    rp.degree   = 1;
    rp.l        = 0;
    rp.verbose  = false;

    int ch;

    while ( (ch = getopt(argc, argv, "k:l:n:p:v")) != -1 )
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
                rp.degree = degree;
                break;

            case 'l':
                rp.l = atoi(optarg);
                if (l < -1 or l > 1)
                {
                    std::cout << "l can be -1, 0 or 1. Falling back to 0." << std::endl;
                    rp.l = 0;
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

            case 'p':
                plot_filename = optarg;
                break;

            case 'v':
                rp.verbose = true;
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

    mesh_filename = argv[0];
    
#if 0
    if (argc == 0)
    {
        std::cout << "Mesh format: 1D uniform" << std::endl;
        auto msh = disk::load_uniform_1d_mesh<RealType>(0, 1, elems_1d);
        run_diffusion_solver(msh, rp);
        return 0;
    }

    

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        auto msh = disk::load_fvca5_2d_mesh<RealType>(mesh_filename);
        run_diffusion_solver(msh, rp);
        return 0;
    }
#endif
    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);
        run_diffusion_solver(msh, rp);
        return 0;
    }
#if 0
    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        auto msh = disk::load_cartesian_2d_mesh<RealType>(mesh_filename);
        run_diffusion_solver(msh, rp);
        return 0;
    }

    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        auto msh = disk::load_netgen_3d_mesh<RealType>(mesh_filename);
        run_diffusion_solver(msh, rp);
        return 0;
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        auto msh = disk::load_cartesian_3d_mesh<RealType>(mesh_filename);
        run_diffusion_solver(msh, rp);
        return 0;
    }
#endif

}




#if 0



int main(int argc, char **argv)
{
    using RealType = double;

    char    *filename       = nullptr;
    int     degree          = 1;
    int     l               = 0;
    int     elems_1d        = 8;
    bool    submesh_flag    = false;
    int ch;

    while ( (ch = getopt(argc, argv, "k:n:sl:")) != -1 )
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

            case 's':
                submesh_flag = true;
                break;

            case 'l':
                l = atoi(optarg);
                if (l < -1 or l > 1)
                {
                    std::cout << "l can be -1, 0 or 1. Falling back to 0." << std::endl;
                    l = 0;
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

        //test_gradrec(msh, degree);
        test_diffusion(msh, f, sf, degree, "plot.dat");

        return 0;
    }

    filename = argv[0];


    if (std::regex_match(filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;

        typedef disk::generic_mesh<RealType, 3>   mesh_type;

        mesh_type msh;
        disk::fvca6_mesh_loader<RealType, 3> loader;


        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }

        loader.populate_mesh(msh);

        auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return M_PI * M_PI * sin(p.x() * M_PI);
            //return 1.0;
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return sin(p.x() * M_PI);
            //return -p.x() * p.x() * 0.5;
        };

        double m = 0.0;
        for (auto& cl : msh)
            m += measure(msh, cl);
        std::cout << m << std::endl;

        test_diffusion(msh, f, sf, degree, "plot.dat");
        //test_gradrec(msh, degree);
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
            //return 1.0;
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return sin(p.x() * M_PI);
            //return -p.x() * p.x() * 0.5;
        };

        test_diffusion(msh, f, sf, degree, "plot.dat");
        //test_gradrec(msh, degree);
    }

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

        auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return M_PI * M_PI * sin(p.x() * M_PI);
            //return 1.0;
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return sin(p.x() * M_PI);
            //return -p.x() * p.x() * 0.5;
        };

        test_diffusion(msh, f, sf, degree, "plot.dat");
        //test_gradrec(msh, degree);
    }

    if (std::regex_match(filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: Cartesian 3D" << std::endl;

        typedef disk::cartesian_mesh<RealType, 3>   mesh_type;

        mesh_type msh;
        disk::cartesian_mesh_loader<RealType, 3> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return M_PI * M_PI * sin(p.x() * M_PI);
            //return 1.0;
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return sin(p.x() * M_PI);
            //return -p.x() * p.x() * 0.5;
        };

        test_diffusion(msh, f, sf, degree, "plot.dat");
        //test_gradrec(msh, degree);
    }

    return 0;
}

#endif
