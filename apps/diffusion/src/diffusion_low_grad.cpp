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
auto
run_hho_diffusion_solver(const Mesh& msh, const size_t degree)
{
    using T = typename Mesh::coordinate_type;

    hho_degree_info hdi(degree+1, degree);

    auto rhs_fun = make_rhs_function(msh);
    auto sol_fun = make_solution_function(msh);

    auto assembler = make_diffusion_assembler(msh, hdi);

    for (auto& cl : msh)
    {
        auto cb   = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto G    = make_vector_hho_gradrec(msh, cl, hdi);
        auto gr   = make_scalar_hho_laplacian(msh, cl, hdi);
        auto stab = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);
        auto rhs  = make_rhs(msh, cl, cb, rhs_fun);

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = G.second + stab;
        const auto sc  = make_scalar_static_condensation(msh, cl, hdi, A, rhs);
        assembler.assemble(msh, cl, sc.first, sc.second, sol_fun);
    }

    assembler.finalize();

    size_t systsz = assembler.LHS.rows();
    size_t nnz = assembler.LHS.nonZeros();

    #if 0
    std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
    std::cout << "Mesh faces: " << msh.faces_size() << std::endl;
    std::cout << "Dofs: " << systsz << std::endl;
    #endif
    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(systsz);

    std::cout << "Running MUMPS" << std::endl;
    sol = mumps_lu(assembler.LHS, assembler.RHS);

    T error = 0.0;

    std::ofstream ofs("sol.dat");

    for (auto& cl : msh)
    {
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto G      = make_vector_hho_gradrec(msh, cl, hdi);
        auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
        auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);
        auto rhs    = make_rhs(msh, cl, cb, rhs_fun);

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = G.second + stab;

        Eigen::Matrix<T, Eigen::Dynamic, 1> locsol =
            assembler.take_local_data(msh, cl, sol, sol_fun);

        Eigen::Matrix<T, Eigen::Dynamic, 1> fullsol =
            make_scalar_static_decondensation(msh, cl, hdi, A, rhs, locsol);

        Eigen::Matrix<T, Eigen::Dynamic, 1> realsol = project_function(msh, cl, hdi, sol_fun);

        auto diff = realsol - fullsol;
        error += diff.dot(A*diff);

        auto bar = barycenter(msh, cl);

        for (size_t i = 0; i < Mesh::dimension; i++)
            ofs << bar[i] << " ";
        ofs << fullsol(0) << std::endl;

    }

    ofs.close();
    return std::sqrt(error);
}

template<typename Mesh>
auto
run_diffusion_solver(const Mesh& msh, const size_t k)
{
    using T = typename Mesh::coordinate_type;

    T error = run_hho_diffusion_solver(msh, k);

    return error;
}


void convergence_test_typ1(void)
{
    using T = double;
    bool use_sym_grad = false;
    std::vector<std::string> meshfiles;
    /*
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_1.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_2.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_3.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_4.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_5.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_6.typ1");
    */

    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_1.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_2.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_3.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_4.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_5.typ1");

    /*
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_1.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_2.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_3.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_4.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_5.typ1");
    */
    std::cout << "                   velocity H1-error";
    std::cout << "    -     pressure L2-error "<< std::endl;

    for (size_t k = 0; k < 3; k++)
    {
        std::cout << "DEGREE " << k << std::endl;

        std::vector<T> mesh_hs;
        std::vector<T> errors;

        //for (size_t i = 0; i < 2; i++)
        for (size_t i = 0; i < meshfiles.size(); i++)
        {
            typedef disk::generic_mesh<T, 2>  mesh_type;

            //std::cout << " Mesh : "<< i << std::endl;
            mesh_type msh;
            disk::fvca5_mesh_loader<T, 2> loader;
            if (!loader.read_mesh(meshfiles.at(i)))
            {
                std::cout << "Problem loading mesh." << std::endl;
                continue;
            }
            loader.populate_mesh(msh);

            auto error = run_diffusion_solver(msh, k);

            mesh_hs.push_back( disk::average_diameter(msh) );
            errors.push_back(error);
        }

        for (size_t i = 0; i < mesh_hs.size(); i++)
        {
            if (i == 0)
            {
                std::cout << "    ";
                std::cout << std::scientific << std::setprecision(4) << mesh_hs.at(i) << "    ";
                std::cout << std::scientific << std::setprecision(4) << errors.at(i) << std::endl;
            }
            else
            {
                auto rate = std::log( errors.at(i)/errors.at(i-1) ) /
                            std::log( mesh_hs.at(i)/mesh_hs.at(i-1) );
                std::cout << "    ";
                std::cout << std::scientific  << std::setprecision(4) << mesh_hs.at(i) << "    ";
                std::cout << std::scientific  << std::setprecision(4) << errors.at(i)  << "    ";
                std::cout << std::fixed<< std::setprecision(2) << rate <<  std::endl;
            }
        }
    }
}

int main(void)
{
    convergence_test_typ1();
}

#if 0
int main(int argc, char **argv)
{
    using T = double;

    if (argc != 2)
    {
        std::cout << "Please specify file name." << std::endl;
        return 1;
    }

    char *mesh_filename = argv[1];

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        auto msh = disk::load_fvca5_2d_mesh<T>(mesh_filename);
        run_diffusion_solver(msh);
        return 0;
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<T>(mesh_filename);

        std::cout << msh.faces_size() << std::endl;

        run_diffusion_solver(msh);
        return 0;
    }

    /* DiSk++ cartesian 2D */
    /*
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        auto msh = disk::load_cartesian_2d_mesh<T>(mesh_filename);
        run_diffusion_solver(msh);
        return 0;
    }
    */

    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        auto msh = disk::load_netgen_3d_mesh<T>(mesh_filename);
        run_diffusion_solver(msh);
        return 0;
    }

    /* DiSk++ cartesian 3D */
    /*
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        auto msh = disk::load_cartesian_3d_mesh<T>(mesh_filename);
        run_diffusion_solver(msh);
        return 0;
    }
    */

}
#endif
