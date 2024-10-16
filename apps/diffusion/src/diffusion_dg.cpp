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

#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/hho"
#include "diskpp/methods/dg"
#include "diskpp/output/silo.hpp"

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

template<typename Mesh>
void
run_diffusion_solver(Mesh& msh, size_t degree, const typename Mesh::coordinate_type eta)
{   
    auto cvf = connectivity_via_faces(msh);
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    auto f = make_rhs_function(msh);

    auto cbs = disk::scalar_basis_size(degree, Mesh::dimension);
    auto assm = make_discontinuous_galerkin_assembler(msh, cbs);
    
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
            loc_rhs += qp.weight() * phi * f(qp.point());
        }
        
        auto fcs = faces(msh, tcl);
        for (auto& fc : fcs)
        {
            auto nv = cvf.neighbour_via(msh, tcl, fc);
            
            auto n     = normal(msh, tcl, fc);
            auto eta_l = eta / diameter(msh, fc);
            auto f_qps = disk::integrate(msh, fc, 2*degree);
            
            if (nv) {
                matrix_type Att = matrix_type::Zero(tbasis.size(), tbasis.size());
                matrix_type Atn = matrix_type::Zero(tbasis.size(), tbasis.size());
                
                auto ncl = nv.value();
                auto nbasis = disk::make_scalar_monomial_basis(msh, ncl, degree);
                assert(tbasis.size() == nbasis.size());
                
                for (auto& fqp : f_qps) {
                    auto ep     = fqp.point();
                    auto tphi   = tbasis.eval_functions(ep);
                    auto tdphi  = tbasis.eval_gradients(ep);
                
                    /* NOT on a boundary */
                    Att += + fqp.weight() * eta_l * tphi * tphi.transpose();
                    Att += - fqp.weight() * 0.5 * tphi * (tdphi*n).transpose();
                    Att += - fqp.weight() * 0.5 * (tdphi*n) * tphi.transpose();
                
                    auto nphi   = nbasis.eval_functions(ep);
                    auto ndphi  = nbasis.eval_gradients(ep);
                
                    Atn += - fqp.weight() * eta_l * tphi * nphi.transpose();
                    Atn += - fqp.weight() * 0.5 * tphi * (ndphi*n).transpose();
                    Atn += + fqp.weight() * 0.5 * (tdphi*n) * nphi.transpose();
                }
                assm.assemble(msh, tcl, tcl, Att);
                assm.assemble(msh, tcl, ncl, Atn);
            }
            else {
                matrix_type Att = matrix_type::Zero(tbasis.size(), tbasis.size());
                for (auto& fqp : f_qps) {
                    auto ep     = fqp.point();
                    auto tphi   = tbasis.eval_functions(ep);
                    auto tdphi  = tbasis.eval_gradients(ep);
                    
                    /* On a boundary*/
                    Att += + fqp.weight() * eta_l * tphi * tphi.transpose();
                    Att += - fqp.weight() * tphi * (tdphi*n).transpose();
                    Att += - fqp.weight() * (tdphi*n) * tphi.transpose();
                    
                    //loc_rhs -= fqp.weight() * (tdphi*n);
                    //loc_rhs += fqp.weight() * eta_l * tphi;
                }
                assm.assemble(msh, tcl, tcl, Att);
            }   
        }
        
        assm.assemble(msh, tcl, K, loc_rhs);
    }

    assm.finalize();

    std::cout << "Mesh has " << msh.cells_size() << " elements." << std::endl;
    std::cout << "System has " << assm.LHS.rows() << " unknowns and ";
    std::cout << assm.LHS.nonZeros() << " nonzeros." << std::endl;

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.syssz);

    std::cout << "Running MUMPS" << std::endl;
    sol = mumps_lu(assm.LHS, assm.RHS);

    auto sol_fun = make_solution_function(msh);

    std::vector<double> data;

    T err = 0.0; size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, degree);

        Matrix<T, Dynamic, Dynamic> MMe = disk::make_mass_matrix(msh, cl, cb);
        Matrix<T, Dynamic, 1> arhs = disk::make_rhs(msh, cl, cb, sol_fun);
        Matrix<T, Dynamic, 1> asol = MMe.llt().solve(arhs);
        Matrix<T, Dynamic, 1> lsol = sol.segment(cell_i*cb.size(), cb.size());


        Matrix<T, Dynamic, 1> diff = lsol - asol;

        err += diff.dot(MMe*diff);

        data.push_back( lsol(0) );

        cell_i++;
    }

    disk::silo_database silo_db;
    silo_db.create("diffusion.silo");
    silo_db.add_mesh(msh, "mesh");

    disk::silo_zonal_variable<T> u("u", data);
    silo_db.add_variable("mesh", u);

    silo_db.close();

    std::cout << "h = " << disk::average_diameter(msh) << " ";
    std::cout << "err = " << std::sqrt(err) << std::endl;
}



int main(int argc, char **argv)
{
    rusage_monitor rm;

    double      stab_param = 1.0;
    size_t      degree = 1;
    char *      mesh_filename = nullptr;

    int ch;
    while ( (ch = getopt(argc, argv, "a:k:m:")) != -1 )
    {
        switch(ch)
        {
            case 'a':
                stab_param = std::stod(optarg);
                break;

            case 'k':
                degree = std::stoi(optarg);
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

    return disk::dispatch_all_meshes(mesh_filename,
              [](auto ...args) { run_diffusion_solver(args...); },
              degree, stab_param );
}

