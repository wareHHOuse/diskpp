#include <iostream>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include <map>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/hho"
#include "diskpp/mesh/meshgen.hpp"


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

    if (mesh_filename)
    {
        disk::dispatch_all_meshes(mesh_filename,
            [](auto& ...args) { compute_eigs(args...); },
            degree );
    }
    else
    {
        using T = double;
        disk::triangular_mesh<T> msh;
        auto mesher = make_simple_mesher(msh);
        
        for (size_t i = 0; i < 4; i++)
        {
            mesher.refine();
            compute_eigs(msh, degree);
        }
    }
}