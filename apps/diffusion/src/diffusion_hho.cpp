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

#if 0
template<typename Mesh>
struct test_functor
{
    /* Expect k+1 convergence (hho stabilization, energy norm) */
    typename Mesh::coordinate_type
    operator()(const Mesh& msh, size_t degree) const
    {
        return run_hho_diffusion_solver(msh, degree);
    }

    size_t
    expected_rate(size_t k)
    {
        return k+1;
    }
};

template<typename Mesh>
void
run_diffusion_solver(const Mesh& msh)
{
    run_hho_diffusion_solver(msh, 0);
}
#endif




int main(int argc, char **argv)
{
    rusage_monitor rm;

    using T = double;

    size_t      num_elems = 16;
    hho_degree_info hdi(1);
    char *      mesh_filename = nullptr;
    bool        stat_cond = true, stab_diam_F = true;

    size_t tmpdeg;
    bool rflag = false;
    bool lflag = false;

    int ch;
    while ( (ch = getopt(argc, argv, "k:m:N:r:l:")) != -1 )
    {
        switch(ch)
        {
            case 'k':
                tmpdeg = std::stoi(optarg);
                if (not lflag)
                    hdi.cell_degree(tmpdeg);
                
                hdi.face_degree(tmpdeg);
                
                if (not rflag)
                    hdi.reconstruction_degree(tmpdeg+1);
                break;

            case 'r':
                tmpdeg = std::stoi(optarg);
                rflag = true;
                hdi.reconstruction_degree(tmpdeg);
                break;

            case 'l':
                tmpdeg = std::stoi(optarg);
                lflag = true;
                hdi.cell_degree(tmpdeg);
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
        std::cout << "Mesh format: 1D uniform" << std::endl;

        typedef disk::generic_mesh<T, 1>  mesh_type;

        mesh_type msh;
        disk::uniform_mesh_loader<T, 1> loader(0, 1, num_elems);
        loader.populate_mesh(msh);

        stab_diam_F = false;
        run_hho_diffusion_solver(msh, hdi, stat_cond, stab_diam_F);

        return 0;
    }

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        auto msh = disk::load_fvca5_2d_mesh<T>(mesh_filename);
        run_hho_diffusion_solver(msh, hdi, stat_cond, stab_diam_F);
        return 0;
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<T>(mesh_filename);

        std::cout << msh.faces_size() << std::endl;

        run_hho_diffusion_solver(msh, hdi, stat_cond, stab_diam_F);
        return 0;
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        auto msh = disk::load_cartesian_2d_mesh<T>(mesh_filename);
        run_hho_diffusion_solver(msh, hdi, stat_cond, stab_diam_F);
        return 0;
    }


    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        auto msh = disk::load_netgen_3d_mesh<T>(mesh_filename);
        run_hho_diffusion_solver(msh, hdi, stat_cond, stab_diam_F);
        return 0;
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        auto msh = disk::load_cartesian_3d_mesh<T>(mesh_filename);
        run_hho_diffusion_solver(msh, hdi, stat_cond, stab_diam_F);
        return 0;
    }

    /* FVCA6 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
        disk::generic_mesh<T,3> msh;

        disk::load_mesh_fvca6_3d<T>(mesh_filename, msh);

        run_hho_diffusion_solver(msh, hdi, stat_cond, stab_diam_F);

        return 0;
    }
}


#if 0
int main(void)
{
    tester<test_functor> tstr;
    tstr.run();
    return 0;
}
#endif
