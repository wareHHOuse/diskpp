#include <iostream>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include <map>

#include "colormanip.h"

#include "geometry/geometry.hpp"
#include "loaders/loader.hpp"
#include "revolution/methods/hho"
#include "apps/intissar/methods_hho_scalar_neumann.hpp"
#include "solvers/solver.hpp"





/***************************************************************************/
/* RHS definition */
template<typename Mesh>
struct rhs_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto cos_py = std::cos(M_PI * pt.y());
        auto sin_py = std::sin(M_PI * pt.y());
        return 2.0 * M_PI * M_PI * sin_px * cos_py;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        auto sin_pz = std::sin(M_PI * pt.z());
        auto cos_py = std::cos(M_PI * pt.y());
        return 3.0 * M_PI * M_PI * sin_px * cos_py * sin_pz;
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
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto cos_py = std::cos(M_PI * pt.y());
        auto sin_py = std::sin(M_PI * pt.y());
        return sin_px * cos_py;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        auto sin_pz = std::sin(M_PI * pt.z());
        auto cos_py = std::cos(M_PI * pt.y());
        return sin_px * cos_py * sin_pz;
    }
};

template<typename Mesh>
auto make_solution_function(const Mesh& msh)
{
    return solution_functor<Mesh>();
}

/***************************************************************************/
/* NEUMANN function */
template<typename Mesh>
struct neumann_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct neumann_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        //auto sin_px = std::sin(M_PI * pt.x());
        //auto sin_py = std::sin(M_PI * pt.y());
        auto cos_px = std::cos(M_PI * pt.x());
        auto cos_py = std::cos(M_PI * pt.y());
        return   M_PI * cos_px * cos_py;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct neumann_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto cos_px = std::cos(M_PI * pt.x());
        auto cos_py = std::cos(M_PI * pt.y());
        auto cos_pz = std::cos(M_PI * pt.z());
        return M_PI * cos_px * cos_py * cos_pz;
    }
};

template<typename Mesh>
 auto make_neumann_function(const Mesh& msh)
{
    return neumann_functor<Mesh>();
}



/***************************************************************************/
/* ROBIN function */
template<typename Mesh>
struct robin_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct robin_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        auto cos_py = std::cos(M_PI * pt.y());
        auto cos_px = std::cos(M_PI * pt.x());
        return M_PI * cos_px * cos_py + sin_px * cos_py;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct robin_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto cos_px = std::cos(M_PI * pt.x());
        auto cos_py = std::cos(M_PI * pt.y());
        auto cos_pz = std::cos(M_PI * pt.z());
        return M_PI * cos_px * cos_py * cos_pz;
    }
};

template<typename Mesh>
 auto make_robin_function(const Mesh& msh)
{
    return robin_functor<Mesh>();
}

using namespace revolution;

template<typename Mesh>
auto
run_test_neumann(const Mesh& msh)
{
    using T = typename Mesh::coordinate_type;
    const size_t degree = 0 ;
    hho_degree_info hdi(degree);
    auto rhs_fun = make_rhs_function(msh);
    auto sol_fun = make_solution_function(msh);

    typedef disk::mechanics::BoundaryConditionsScalar2<Mesh> boundary_type;
    boundary_type  m_bnd(msh);
    auto g = make_neumann_function(msh);

    for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
    {
        const auto bfc = *itor;

        const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

        const auto face_id = eid.second;
        auto b_id= msh.boundary_id(face_id);

        m_bnd.addDirichletBC(disk::mechanics::DIRICHLET,3,sol_fun);
        m_bnd.addDirichletBC(disk::mechanics::DIRICHLET,4, sol_fun);
        m_bnd.addNeumannBC(disk::mechanics::NEUMANN,1, g);
        m_bnd.addNeumannBC(disk::mechanics::NEUMANN,2, g);

    }
    auto assembler = make_diffusion_assembler2(msh, hdi, m_bnd);


    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
        auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
        auto rhs    = make_rhs(msh, cl, cb, rhs_fun);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = gr.second + stab;    // if ROBIN: see methods_hho.hpp line 150 and line 285
        auto sc     = diffusion_static_condensation_compute(msh, cl, hdi, A, rhs);
        assembler.assemble(msh, cl, sc.first, sc.second);
    }

    assembler.neumann_boundary_function(msh, m_bnd);
    assembler.finalize();

    size_t systsz = assembler.LHS.rows();
    size_t nnz = assembler.LHS.nonZeros();

    std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
    std::cout << "Mesh faces: " << msh.faces_size() << std::endl;
    std::cout << "Dofs: " << systsz << std::endl;

    dynamic_vector<T> sol = dynamic_vector<T>::Zero(systsz);

    disk::solvers::pardiso_params<T> pparams;
    pparams.report_factorization_Mflops = true;
    mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

    T error = 0.0 ;
    std::ofstream ofs("sol.dat");

    if(!ofs.is_open())
    {
        std::cout<< "Error opening file"<<std::endl;
    }


    for (auto& cl : msh)
    {
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
        auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
        auto rhs    = make_rhs(msh, cl, cb, rhs_fun);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = gr.second + stab;

        Eigen::Matrix<T, Eigen::Dynamic, 1> locsol =
            assembler.take_local_data(msh, cl, sol);

        Eigen::Matrix<T, Eigen::Dynamic, 1> fullsol =
            diffusion_static_condensation_recover(msh, cl, hdi, A, rhs, locsol);

        Eigen::Matrix<T, Eigen::Dynamic, 1> realsol = project_function(msh, cl, hdi, sol_fun);

        auto diff = realsol - fullsol;
        error += diff.dot(A*diff);

        auto bar = barycenter(msh, cl);

        for (size_t i = 0; i < Mesh::dimension; i++)
            ofs << bar[i] << " ";
        ofs << fullsol(0) << std::endl;

    }
    std::cout<<"L'erreur: "<<std::endl;
    std::cout << std::sqrt(error) << std::endl;

    ofs.close();
    return std::sqrt(error);
}

template<typename Mesh>
auto
run_test_dirichlet(const Mesh& msh)
{
    using T = typename Mesh::coordinate_type;
    const size_t degree = 0 ;
    hho_degree_info hdi(degree);
    auto rhs_fun = make_rhs_function(msh);
    auto sol_fun = make_solution_function(msh);

    typedef disk::mechanics::BoundaryConditionsScalar2<Mesh> boundary_type;
    boundary_type  m_bnd(msh);
    m_bnd.addDirichletEverywhere(sol_fun);
    auto assembler = make_diffusion_assembler2(msh, hdi, m_bnd);

    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
        auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
        auto rhs    = make_rhs(msh, cl, cb, rhs_fun);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = gr.second + stab;    // if ROBIN: see methods_hho.hpp line 150 and line 285
        auto sc     = diffusion_static_condensation_compute(msh, cl, hdi, A, rhs);
        assembler.assemble(msh, cl, sc.first, sc.second);
    }
    assembler.finalize();

    size_t systsz = assembler.LHS.rows();
    size_t nnz = assembler.LHS.nonZeros();

    std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
    std::cout << "Mesh faces: " << msh.faces_size() << std::endl;
    std::cout << "Dofs: " << systsz << std::endl;

    dynamic_vector<T> sol = dynamic_vector<T>::Zero(systsz);

    disk::solvers::pardiso_params<T> pparams;
    pparams.report_factorization_Mflops = true;
    mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

    T error = 0.0 ;
    std::ofstream ofs("sol.dat");

    if(!ofs.is_open())
    {
        std::cout<< "Error opening file"<<std::endl;
    }


    for (auto& cl : msh)
    {
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
        auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
        auto rhs    = make_rhs(msh, cl, cb, rhs_fun);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = gr.second + stab;

        Eigen::Matrix<T, Eigen::Dynamic, 1> locsol =
            assembler.take_local_data(msh, cl, sol);

        Eigen::Matrix<T, Eigen::Dynamic, 1> fullsol =
            diffusion_static_condensation_recover(msh, cl, hdi, A, rhs, locsol);

        Eigen::Matrix<T, Eigen::Dynamic, 1> realsol = project_function(msh, cl, hdi, sol_fun);

        auto diff = realsol - fullsol;
        error += diff.dot(A*diff);

        auto bar = barycenter(msh, cl);

        for (size_t i = 0; i < Mesh::dimension; i++)
            ofs << bar[i] << " ";
        ofs << fullsol(0) << std::endl;

    }
    std::cout<<"L'erreur: "<<std::endl;
    std::cout << std::sqrt(error) << std::endl;

    ofs.close();
    return std::sqrt(error);
}


template<typename Mesh>
auto
run_test_robin(const Mesh& msh)
{
    using T = typename Mesh::coordinate_type;
    const size_t degree = 0 ;
    hho_degree_info hdi(degree);
    auto rhs_fun = make_rhs_function(msh);
    auto sol_fun = make_solution_function(msh);

    typedef disk::mechanics::BoundaryConditionsScalar2<Mesh> boundary_type;
    boundary_type  m_bnd(msh);

    auto g = make_robin_function (msh);
    for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
    {
        const auto bfc = *itor;

        const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

        const auto face_id = eid.second;
        auto b_id= msh.boundary_id(face_id);
        m_bnd.addDirichletBC(disk::mechanics::DIRICHLET, 1, sol_fun);
        m_bnd.addDirichletBC(disk::mechanics::DIRICHLET, 0, sol_fun);
        m_bnd.addRobinBC(disk::mechanics::ROBIN, 3, g);
        m_bnd.addRobinBC(disk::mechanics::ROBIN, 2, g);
    }

    auto assembler = make_diffusion_assembler2(msh, hdi, m_bnd);

    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
        auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
        auto rhs    = make_rhs(msh, cl, cb, rhs_fun);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = gr.second + stab;    // if ROBIN: see methods_hho.hpp line 150 and line 285
        auto sc     = diffusion_static_condensation_compute(msh, cl, hdi, A, rhs);
        assembler.assemble(msh, cl, sc.first, sc.second);
    }

    assembler.robin_boundary_function(msh, m_bnd);
    assembler.finalize();

    size_t systsz = assembler.LHS.rows();
    size_t nnz = assembler.LHS.nonZeros();

    std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
    std::cout << "Mesh faces: " << msh.faces_size() << std::endl;
    std::cout << "Dofs: " << systsz << std::endl;

    dynamic_vector<T> sol = dynamic_vector<T>::Zero(systsz);

    disk::solvers::pardiso_params<T> pparams;
    pparams.report_factorization_Mflops = true;
    mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

    T error = 0.0 ;
    std::ofstream ofs("sol.dat");

    if(!ofs.is_open())
    {
        std::cout<< "Error opening file"<<std::endl;
    }


    for (auto& cl : msh)
    {
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
        auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
        auto rhs    = make_rhs(msh, cl, cb, rhs_fun);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = gr.second + stab;

        Eigen::Matrix<T, Eigen::Dynamic, 1> locsol =
            assembler.take_local_data(msh, cl, sol);

        Eigen::Matrix<T, Eigen::Dynamic, 1> fullsol =
            diffusion_static_condensation_recover(msh, cl, hdi, A, rhs, locsol);

        Eigen::Matrix<T, Eigen::Dynamic, 1> realsol = project_function(msh, cl, hdi, sol_fun);

        auto diff = realsol - fullsol;
        error += diff.dot(A*diff);

        auto bar = barycenter(msh, cl);

        for (size_t i = 0; i < Mesh::dimension; i++)
            ofs << bar[i] << " ";
        ofs << fullsol(0) << std::endl;

    }
    std::cout<<"L'erreur: "<<std::endl;
    std::cout << std::sqrt(error) << std::endl;

    ofs.close();
    return std::sqrt(error);
}
