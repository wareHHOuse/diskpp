#include <iostream>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include <map>

#include "colormanip.h"

#include "geometry/geometry.hpp"
#include "loaders/loader.hpp"
//#include "revolution/methods/hho"
#include "solvers/solver.hpp"
#include "methods_hho.hpp"


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
        return  -2.0 * M_PI * std::sin(2. * M_PI * pt.x());
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
        return 0.;
    }
};


template<typename Mesh>
auto make_rhs_function(const Mesh& msh)
{
    return rhs_functor<Mesh>();
}

/* Expected solution definition */
template<typename Mesh>
struct dirichlet_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct dirichlet_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T,2,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        return 0.;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct dirichlet_functor< Mesh<T, 3, Storage> >
{
    typedef Mesh<T,3,Storage>               mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        return 0.;
    }
};

template<typename Mesh>
auto make_dirichlet_function(const Mesh& msh)
{
    return dirichlet_functor<Mesh>();
}

/* NEUMANN */
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
        return 0.;
    }
};

template<typename Mesh>
 auto make_neumann_function(const Mesh& msh)
{
    return neumann_functor<Mesh>();
}

using namespace revolution;

template<typename Mesh>
auto
full_offset(const Mesh& msh, const hho_degree_info& hdi)
{
    size_t  i = 1, dofs = 0;
    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    std::vector<size_t> vec(msh.cells_size());
    for (auto itor = msh.cells_begin(); itor != msh.cells_end() -1; itor++)
    {
        auto cl = *itor;
        auto fcs = faces(msh, cl);
        dofs += fcs.size() * fbs + cbs;
        vec.at(i++) = dofs;
    }
    return vec;
}

template<typename Mesh>
auto
fix_point_solver(const Mesh& msh, const hho_degree_info& hdi,
                const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
                const std::vector<size_t>& is_contact_vector)
{
    using T = typename Mesh::coordinate_type;

    auto gamma_0 = 0.1;
    auto max_iter = 1000;
    auto tol = 1.e-6;
    auto rhs_fun = make_rhs_function(msh);

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_total_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                        - fbs*msh.boundary_faces_size() ;
    dynamic_vector<T> full_sol_old = dynamic_vector<T>::Zero(num_total_dofs);
    dynamic_vector<T> full_sol = dynamic_vector<T>::Zero(num_total_dofs);


    auto offset_vector = full_offset(msh, hdi);

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        auto cl_count = 0;
        full_sol_old = full_sol;
        auto assembler = make_diffusion_assembler2(msh, hdi, bnd);
        for (auto& cl : msh)
        {
            auto fcs = faces(msh, cl);
            auto num_dofs = cbs + fcs.size() *fbs;
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
            auto rhs_f  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            Matrix<T, Dynamic, Dynamic> A = gr.second + stab;
            Matrix<T, Dynamic, 1> rhs = Matrix<T, Dynamic, 1>(num_dofs);
            rhs.block(0, 0, cbs, 1) = rhs_f;

            if (is_contact_vector.at(cl_count))
            {
                auto cell_ofs = offset_vector.at(cl_count);
                Matrix<T, Dynamic,1> u_full = full_sol_old.block(cell_ofs, 0, num_dofs, 1);

            	Matrix<T, Dynamic, Dynamic> A_cont = make_contact_unnamed(msh, cl, hdi, gr.first, gamma_0, bnd);
                Matrix<T, Dynamic,1> rhs_cont = make_contact_negative(msh, cl, hdi, gr.first, gamma_0, bnd, u_full);

                #if 0
                std::cout << "A : " << std::endl;
                for(size_t row = 0; row < A_cont.rows(); row++)
                {
                    for(size_t col = 0; col < A_cont.cols(); col++)
                        std::cout << A_cont(row, col)<< " ";
                    std::cout << std::endl;
                }

                std::cout << "rhs : " << std::endl;
                for(size_t row = 0; row < A_cont.rows(); row++)
                    std::cout << rhs_cont(row, 0)<< " ";
                std::cout << std::endl;
                #endif
                A -= A_cont;
                rhs -=  rhs_cont;
            }

            auto sc = diffusion_static_condensation_compute_full(msh, cl, hdi, A, rhs);
            assembler.assemble(msh, cl, sc.first, sc.second);
            cl_count++;
        }

        assembler.neumann_boundary_function(msh, bnd);
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
            std::cout<< "Error opening file"<<std::endl;

        cl_count = 0;
        for (auto& cl : msh)
        {
            auto fcs = faces(msh, cl);
            auto num_dofs = cbs + fcs.size() * fbs;
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
            Matrix<T, Dynamic, 1> rhs  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            Matrix<T, Dynamic, Dynamic> A = gr.second + stab;

            if (is_contact_vector.at(cl_count))
            {
                auto cell_ofs = offset_vector.at(cl_count);
                Matrix<T, Dynamic,1> u_full = full_sol_old.block(cell_ofs, 0, num_dofs, 1);

            	Matrix<T, Dynamic, Dynamic> A_cont = make_contact_unnamed(msh, cl, hdi, gr.first, gamma_0, bnd );
                Matrix<T, Dynamic, 1> rhs_cont = make_contact_negative(msh, cl, hdi, gr.first, gamma_0, bnd, u_full);

                A -= A_cont;
                rhs -=  rhs_cont.block(0, 0, cbs, 1);
            }

            Matrix<T, Dynamic, 1> locsol =
                                        assembler.take_local_data(msh, cl, sol);
            Matrix<T, Dynamic, 1> u_full =
                diffusion_static_condensation_recover(msh, cl, hdi, A, rhs, locsol);

            auto cell_ofs = offset_vector.at(cl_count);
            Matrix<T, Dynamic, 1> u_full_old = full_sol_old.block(cell_ofs, 0, num_dofs ,1);
            full_sol.block(cell_ofs, 0, num_dofs ,1) = u_full;

            Matrix<T, Dynamic, 1>  diff = u_full - u_full_old;
            error += diff.dot(A*diff);

            auto bar = barycenter(msh, cl);
            ofs << bar.x() << " " << bar.y() <<" "<< u_full(0) << std::endl;
            cl_count++;
        }
        ofs.close();

        std::cout << "iter : "<< iter << "   ; error = "<< std::sqrt(error)<< std::endl;
        if( std::sqrt(error)  < tol)
        {
            //for (size_t i = 0; i < Mesh::dimension; i++)
            //    ofs << bar[i] << " ";
            //ofs << full_sol(0) << std::endl;
            return 0;

        }
    }

}

//template<typename Mesh>
auto
run_signorini()//const Mesh& msh)
{
    using T = double;
    const size_t degree = 0 ;
    hho_degree_info hdi(degree);
    typedef disk::generic_mesh<T, 2>  Mesh;

    std::cout << "Mesh format: Medit format" << std::endl;
    auto filename = "../../../diskpp/meshes/2D_quads/medit/square_h0025.medit2d";
    auto msh = disk::load_medit_2d_mesh<T>(filename);

    auto g_N = make_neumann_function(msh);
    auto u_D = make_dirichlet_function(msh);

    typedef disk::mechanics::BoundaryConditionsScalar<Mesh> boundary_type;
    boundary_type  m_bnd(msh);

    m_bnd.addDirichletBC(disk::mechanics::DIRICHLET,1,u_D);
    m_bnd.addNeumannBC(disk::mechanics::NEUMANN,3,g_N);
    m_bnd.addNeumannBC(disk::mechanics::NEUMANN,4,g_N);
    m_bnd.addContactBC(disk::mechanics::SIGNORINI,2);

    //cells with contact faces
    auto num_cells = msh.cells_size();
    std::vector<size_t> is_contact_vector = std::vector<size_t>(num_cells);
    size_t i =0;
    for (auto& cl : msh)
    {
        auto fcs=faces(msh, cl);
	    for (auto& fc:fcs)
	    {
	    	auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
            const auto face_id=eid.second;

		    if (m_bnd.is_contact_face(face_id))
		    {
			    is_contact_vector.at(i) = 1;
			    continue;
		    }
	    }
	    i++;
    }

    fix_point_solver( msh, hdi, m_bnd, is_contact_vector);

    return;
}

int main(void)
{
    run_signorini();
}
