/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Karol Cascavita (C) 2018                     karol.cascavita@enpc.fr
 * Intissar Addali (C) 2018                     intissar.addali@inria.fr
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

 #include "output/silo.hpp"
#include "signorini_solver.hpp"


enum method_type
{
    POSITIVE,
    NEGATIVE,
    ZERO
};

enum solver_type
{
    NEWTON,
    FIX_POINT,
    MIXED
};

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
        return - 2.* M_PI *  std::sin(2. * M_PI * pt.x());
    }
};

template<typename Mesh>
auto make_rhs_function(const Mesh& msh)
{
    return rhs_functor<Mesh>();
}


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


template<typename Mesh, typename T>
void
run_signorini(  const Mesh& msh,
                algorithm_parameters<T>& ap,
                const method_type& tt,
                const solver_type& ss)
{
    const size_t degree = 0 ;
    hho_degree_info hdi(degree+1, degree);

    auto g_N = make_neumann_function(msh);
    auto u_D = make_dirichlet_function(msh);
    auto rhs_fun = make_rhs_function(msh);

    typedef disk::mechanics::BoundaryConditionsScalar<Mesh> boundary_type;
    boundary_type  m_bnd(msh);

    m_bnd.addDirichletBC(disk::mechanics::DIRICHLET,1,u_D); //TOP
    m_bnd.addNeumannBC(disk::mechanics::NEUMANN,3,g_N); //
    m_bnd.addNeumannBC(disk::mechanics::NEUMANN,4,g_N); //
    m_bnd.addContactBC(disk::mechanics::SIGNORINI,2); //BOTTOM

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                        - fbs*msh.boundary_faces_size() ;

    dynamic_vector<T> full_sol = dynamic_vector<T>::Zero(num_full_dofs);

    if(ss == FIX_POINT)
        full_sol = fix_point_solver( msh, hdi, m_bnd, ap, rhs_fun);
    else if(ss == NEWTON)
        newton_solver( msh, hdi, m_bnd, ap, rhs_fun, full_sol);
    else if(ss == MIXED)
    {
        full_sol = fix_point_solver( msh, hdi, m_bnd, ap, rhs_fun);
        newton_solver( msh, hdi, m_bnd, ap, rhs_fun, full_sol);
    }
    else
        throw std::invalid_argument("No solver chosen");

    return;
}

int main(int argc, char **argv)
{
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    char    *filename       = nullptr;
    using T = double;

    int ch;
    algorithm_parameters<T> ap;
    method_type tt = POSITIVE;
    solver_type ss = FIX_POINT;

    while ( (ch = getopt(argc, argv, "g:npzfwm")) != -1 )
    {
        switch(ch)
        {
            case 'g':
                std::cout << "choosing gamma" << std::endl;
                ap.gamma_0 = atof(optarg);
                if (ap.gamma_0 <= 0)
                {
                    std::cout << "gamma_0 must be >0. Falling back to 0.1" << std::endl;
                    ap.gamma_0 = 0.1;
                }
                break;
            #if 0
            case 't':
                if(atoi(optarg) == 1 || atoi(optarg) == 0 || atoi(optarg) == -1)
                {
                    ap.theta = T(atoi(optarg));
                    std::cout << "theta chosen as "<< ap.theta << std::endl;
                }
                else
                    std::cout << "theta must be in {-1, 0,1}. Falling back to 1." << std::endl;
            #endif
            //#if 0
            case 'n':
                ap.theta = -1.;
                std::cout << "theta negative chosen" << std::endl;
                break;
            case 'p':
                ap.theta = 1.;
                std::cout << "theta positive chosen" << std::endl;
                break;
            case 'z':
                ap.theta = 0.;
                std::cout << "theta zero chosen" << std::endl;
                break;
            case 'f':
                ss = FIX_POINT;
                std::cout << "Fix point Solver chosen" << std::endl;
                break;
            case 'w':
                ss = NEWTON;
                std::cout << "Newton Solver chosen" << std::endl;
                break;
            case 'm':
                ss = MIXED;
                std::cout << "Mixed Solver chosen" << std::endl;
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

    /* Medit 2d*/
    if (std::regex_match(filename, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        typedef disk::generic_mesh<T, 2>  mesh_type;
        mesh_type msh = disk::load_medit_2d_mesh<T>(filename);
        run_signorini(msh, ap, tt, ss);
    }

    return 0;
}
