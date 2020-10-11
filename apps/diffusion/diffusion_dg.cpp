/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *  
 * Matteo Cicuttin (C) 2020
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group  
 */

/*
 *       /\        Matteo Cicuttin (C) 2016-2019
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
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
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include <map>

#include "colormanip.h"

#include "geometry/geometry.hpp"
#include "loaders/loader.hpp"
#include "methods/hho"
#include "solvers/solver.hpp"

#include "../tests/common.hpp"

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

template<typename Mesh>
void
run_diffusion_solver(Mesh& msh)
{
    size_t degree = 1;
    
    auto cvf = connectivity_via_faces(msh);
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;
    
    for (auto& tcl : msh)
    {
        auto tbasis = disk::make_scalar_monomial_basis(msh, tcl, degree);
        auto qps = disk::integrate(msh, tcl, 2*degree);
        
        matrix_type K = matrix_type::Zero(tbasis.size(), tbasis.size());
        vector_type loc_rhs = vector_type::Zero(tbasis.size());
        for (auto& qp : qps)
        {
            auto ep = qp.point();
            auto phi = tbasis.eval_functions(ep);
            auto dphi = tbasis.eval_gradients(ep);
            
            K += qp.weight() * dphi * dphi.transpose();
            loc_rhs += qp.weight() * phi;
        }
        
        auto fcs = faces(msh, tcl);
        for (auto& fc : fcs)
        {
            matrix_type Att = matrix_type::Zero(tbasis.size(), tbasis.size());
            matrix_type Atn = matrix_type::Zero(tbasis.size(), tbasis.size());
            
            auto nv = cvf.neighbour_via(msh, tcl, fc);
            auto ncl = nv.first;
            auto nbasis = disk::make_scalar_monomial_basis(msh, ncl, degree);
            assert(tbasis.size() == nbasis.size());
            
            auto n     = normal(msh, tcl, fc);
            auto eta_l = 1.0;//eta / diameter(msh, fc);
            auto f_qps = disk::integrate(msh, fc, 2*degree);
            
            for (auto& fqp : f_qps)
            {
                auto ep     = fqp.point();
                auto tphi   = tbasis.eval_functions(ep);
                auto tdphi  = tbasis.eval_gradients(ep);
                
                if (nv.second)
                {   /* NOT on a boundary */
                    Att += + fqp.weight() * eta_l * tphi * tphi.transpose();
                    Att += - fqp.weight() * 0.5 * tphi * (tdphi*n).transpose();
                    Att += - fqp.weight() * 0.5 * (tdphi*n) * tphi.transpose();
                }
                else
                {   /* On a boundary*/
                    Att += + fqp.weight() * eta_l * tphi * tphi.transpose();
                    Att += - fqp.weight() * tphi * (tdphi*n).transpose();
                    Att += - fqp.weight() * (tdphi*n) * tphi.transpose();
                    
                    loc_rhs -= fqp.weight() * (tdphi*n);
                    loc_rhs += fqp.weight() * eta_l * tphi;
                    continue;
                }
                
                auto nphi   = nbasis.eval_functions(ep);
                auto ndphi  = nbasis.eval_gradients(ep);
                
                Atn += - fqp.weight() * eta_l * tphi * nphi.transpose();
                Atn += - fqp.weight() * 0.5 * tphi * (ndphi*n).transpose();
                Atn += + fqp.weight() * 0.5 * (tdphi*n) * nphi.transpose();
            }
            
            //assm.assemble(msh, tcl, tcl, Att);
            //if (nv.second)
            //    assm.assemble(msh, tcl, ncl, Atn);
        }
        
        //assm.assemble(msh, tcl, K, loc_rhs);
    }
}

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
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        auto msh = disk::load_cartesian_2d_mesh<T>(mesh_filename);
        run_diffusion_solver(msh);
        return 0;
    }


    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        auto msh = disk::load_netgen_3d_mesh<T>(mesh_filename);
        run_diffusion_solver(msh);
        return 0;
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        auto msh = disk::load_cartesian_3d_mesh<T>(mesh_filename);
        run_diffusion_solver(msh);
        return 0;
    }
}

