#include <iostream>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include <map>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/hho"


template<typename Mesh>
void compute_eigs(const Mesh& msh, size_t sk)
{
    disk::hho_degree_info hdi;
    
    hdi.cell_degree(sk);
    hdi.face_degree(sk);

    for (size_t rk = sk+1; rk < sk+5; rk++)
    {
        std::cout << "HHO space degree: " << sk << ", reconstruction degree: " << rk << std::endl;
        hdi.reconstruction_degree(rk);

        for (auto& cl : msh)
        {
            auto gr = make_scalar_hho_laplacian(msh, cl, hdi);
            std::cout << gr.second.eigenvalues().transpose() << std::endl;
            break;
        }
    }
}


int main(int argc, char **argv)
{
    rusage_monitor rm;

    using T = double;

    size_t      num_elems = 16;
    size_t      degree = 1;
    char *      mesh_filename = nullptr;
    bool        stat_cond = true, stab_diam_F = true;

    int ch;
    while ( (ch = getopt(argc, argv, "k:m:twN:")) != -1 )
    {
        switch(ch)
        {
            case 't': stab_diam_F = false; break;

            case 'w': stat_cond = false; break;

            case 'k':
                degree = std::stoi(optarg);
                break;

            case 'N':
                /* Only for 1D */
                num_elems = std::stoi(optarg);
                break;

            case 'm':
                mesh_filename = optarg;
                break;

            case '?':
            default:
                std::cout << "Invalid option" << std::endl;
                return 1;
        }
    }

    if (mesh_filename == nullptr)
    {
        std::cout << "Please specify the mesh (-m)" << std::endl;
        return 1;
    }

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        auto msh = disk::load_fvca5_2d_mesh<T>(mesh_filename);
        compute_eigs(msh, degree);
        return 0;
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<T>(mesh_filename);
        std::cout << msh.faces_size() << std::endl;
        compute_eigs(msh, degree);
        return 0;
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        auto msh = disk::load_cartesian_2d_mesh<T>(mesh_filename);
        compute_eigs(msh, degree);
        return 0;
    }


    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        auto msh = disk::load_netgen_3d_mesh<T>(mesh_filename);
        compute_eigs(msh, degree);
        return 0;
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        auto msh = disk::load_cartesian_3d_mesh<T>(mesh_filename);
        compute_eigs(msh, degree);
        return 0;
    }

    /* FVCA6 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
        disk::generic_mesh<T,3> msh;
        disk::load_mesh_fvca6_3d<T>(mesh_filename, msh);
        compute_eigs(msh, degree);
        return 0;
    }

    std::cout << "Unable to detect mesh format." << std::endl;
    return 1;
}