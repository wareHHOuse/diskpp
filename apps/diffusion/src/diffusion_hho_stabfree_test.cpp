/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <iostream>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include <map>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/hho"

#include "mumps.hpp"

#include "diffusion_hho_common.hpp"

template<typename Mesh>
void
test_stabfree_hho(Mesh& msh, std::ofstream& ofs)
{
    size_t kmin = 0;
    size_t kmax = 2;
    size_t rincrmin = 1;
    size_t rincrmax = 3;

    
    ofs << "\t";
    for (size_t k = kmin; k <= kmax; k++)
        for (size_t l = k; l <= k+1; l++)
            for (size_t r = k+rincrmin; r <= k+rincrmax; r++)
                ofs << "HHO(" << l << "," << k << "," << r << ")\t";
    ofs << std::endl;

    auto h = disk::average_diameter(msh);

    ofs << h << "\t";

    for (size_t k = kmin; k <= kmax; k++)
    {
        for (size_t l = k; l <= k+1; l++)
        {
            for (size_t r = k+rincrmin; r <= k+rincrmax; r++)
            {
                hho_degree_info hdi;
                hdi.cell_degree(l);
                hdi.face_degree(k);
                hdi.reconstruction_degree(r);
                try {
                    auto error = run_hho_diffusion_solver(msh, hdi, true, true);
                    ofs << error << "\t";
                }
                catch (...) {
                    ofs << "FAIL\t";
                }
            }
        }
    }

    ofs << std::endl;
}

int
test_stabfree_hho_dispatch(const char *mesh_filename, std::ofstream& ofs)
{
    using T = double;

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        auto msh = disk::load_fvca5_2d_mesh<T>(mesh_filename);
        test_stabfree_hho(msh, ofs);
        return 0;
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<T>(mesh_filename);
        test_stabfree_hho(msh, ofs);
        return 0;
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        auto msh = disk::load_cartesian_2d_mesh<T>(mesh_filename);
        test_stabfree_hho(msh, ofs);
        return 0;
    }


    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        auto msh = disk::load_netgen_3d_mesh<T>(mesh_filename);
        test_stabfree_hho(msh, ofs);
        return 0;
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        auto msh = disk::load_cartesian_3d_mesh<T>(mesh_filename);
        test_stabfree_hho(msh, ofs);
        return 0;
    }

    /* FVCA6 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
        disk::generic_mesh<T,3> msh;
        disk::load_mesh_fvca6_3d<T>(mesh_filename, msh);
        test_stabfree_hho(msh, ofs);
        return 0;
    }

    return 1;
}

int main(int argc, char **argv)
{
    rusage_monitor rm;

    std::vector<char *> meshes;

    int ch;
    while ( (ch = getopt(argc, argv, "m:")) != -1 )
    {
        switch(ch)
        {
            case 'm':
                meshes.push_back(optarg);
                break;
        }
    }

    std::ofstream ofs("convtest_stabfree.txt");

    for (const auto& mesh : meshes)
    {
        test_stabfree_hho_dispatch(mesh, ofs);
    }
}

