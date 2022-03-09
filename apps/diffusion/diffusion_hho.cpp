#include <iostream>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include <map>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/hho"

#include "mumps.hpp"

/***************************************************************************/
/* RHS definition */
template<typename Mesh>
struct rhs_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 1, Storage> >
{
    typedef Mesh<T,1,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        return M_PI * M_PI * sin_px;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        return 2.0 * M_PI * M_PI * sin_px * sin_py;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        auto sin_pz = std::sin(M_PI * pt.z());
        return 3.0 * M_PI * M_PI * sin_px * sin_py * sin_pz;
    }
};

template<typename Mesh>
auto make_rhs_function(const Mesh& msh)
{
    return rhs_functor<Mesh>();
}

/***************************************************************************/
/* Expected solution definition */
template<typename Mesh>
struct solution_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_functor< Mesh<T, 1, Storage> >
{
    typedef Mesh<T,1,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        return sin_px;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        return sin_px * sin_py;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        auto sin_pz = std::sin(M_PI * pt.z());
        return sin_px * sin_py * sin_pz;
    }
};

template<typename Mesh>
auto make_solution_function(const Mesh& msh)
{
    return solution_functor<Mesh>();
}

using namespace disk;

template<typename Mesh>
typename Mesh::coordinate_type
run_hho_diffusion_solver(const Mesh& msh, size_t degree, const bool statcond, const bool stab_diam_F)
{
    using T = typename Mesh::coordinate_type;

    hho_degree_info hdi(degree);

    auto rhs_fun = make_rhs_function(msh);
    auto sol_fun = make_solution_function(msh);

    auto assembler = make_diffusion_assembler(msh, hdi);

    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
        auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi, stab_diam_F);
        auto rhs    = make_rhs(msh, cl, cb, rhs_fun);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = gr.second + stab;
        auto sc     = make_scalar_static_condensation(msh, cl, hdi, A, rhs);
        assembler.assemble(msh, cl, sc.first, sc.second, sol_fun);
    }

    assembler.finalize();

    size_t systsz = assembler.LHS.rows();
    size_t nnz = assembler.LHS.nonZeros();

    std::cout << "Mesh has " << msh.cells_size() << " elements." << std::endl;
    std::cout << "System has " << assembler.LHS.rows() << " unknowns and ";
    std::cout << assembler.LHS.nonZeros() << " nonzeros." << std::endl;

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(systsz);

    std::cout << "Running MUMPS" << std::endl;
    sol = mumps_lu(assembler.LHS, assembler.RHS);

    T error = 0.0;

    //std::ofstream ofs("sol.dat");

    for (auto& cl : msh)
    {
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
        auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi, stab_diam_F);
        auto rhs    = make_rhs(msh, cl, cb, rhs_fun);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = gr.second + stab;

        Eigen::Matrix<T, Eigen::Dynamic, 1> locsol =
            assembler.take_local_data(msh, cl, sol, sol_fun);

        Eigen::Matrix<T, Eigen::Dynamic, 1> fullsol = make_scalar_static_decondensation(msh, cl, hdi, A, rhs, locsol);

        Eigen::Matrix<T, Eigen::Dynamic, 1> realsol = project_function(msh, cl, hdi, sol_fun, 2);


        auto diff = realsol - fullsol;
        //error += diff.dot(A*diff);

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MM = make_mass_matrix(msh, cl, cb);

        error += diff.segment(0,cb.size()).dot(MM*diff.segment(0,cb.size()));

        //auto bar = barycenter(msh, cl);

        //for (size_t i = 0; i < Mesh::dimension; i++)
        //    ofs << bar[i] << " ";
        //ofs << fullsol(0) << std::endl;

    }

    std::cout << "h = " << disk::average_diameter(msh) << " ";
    std::cout << "err = " << std::sqrt(error) << std::endl;

    return std::sqrt(error);
}


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
        std::cout << "Mesh format: 1D uniform" << std::endl;

        typedef disk::generic_mesh<T, 1>  mesh_type;

        mesh_type msh;
        disk::uniform_mesh_loader<T, 1> loader(0, 1, num_elems);
        loader.populate_mesh(msh);

        stab_diam_F = false;
        run_hho_diffusion_solver(msh, degree, stat_cond, stab_diam_F);

        return 0;
    }

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        auto msh = disk::load_fvca5_2d_mesh<T>(mesh_filename);
        run_hho_diffusion_solver(msh, degree, stat_cond, stab_diam_F);
        return 0;
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<T>(mesh_filename);

        std::cout << msh.faces_size() << std::endl;

        run_hho_diffusion_solver(msh, degree, stat_cond, stab_diam_F);
        return 0;
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        auto msh = disk::load_cartesian_2d_mesh<T>(mesh_filename);
        run_hho_diffusion_solver(msh, degree, stat_cond, stab_diam_F);
        return 0;
    }


    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        auto msh = disk::load_netgen_3d_mesh<T>(mesh_filename);
        run_hho_diffusion_solver(msh, degree, stat_cond, stab_diam_F);
        return 0;
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        auto msh = disk::load_cartesian_3d_mesh<T>(mesh_filename);
        run_hho_diffusion_solver(msh, degree, stat_cond, stab_diam_F);
        return 0;
    }

    /* FVCA6 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
        disk::generic_mesh<T,3> msh;

        disk::load_mesh_fvca6_3d<T>(mesh_filename, msh);

        run_hho_diffusion_solver(msh, degree, stat_cond, stab_diam_F);

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
