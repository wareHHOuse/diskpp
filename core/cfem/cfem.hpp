/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#pragma once

#include "common/eigen.hpp"

namespace disk { namespace cfem {

template<typename T>
T cfem_measure(const disk::simplicial_mesh<T, 2>& msh,
               const typename disk::simplicial_mesh<T, 2>::cell& cl)
{
    auto pts = points(msh, cl);

    auto x0 = pts[0].x(); auto y0 = pts[0].y();
    auto x1 = pts[1].x(); auto y1 = pts[1].y();
    auto x2 = pts[2].x(); auto y2 = pts[2].y();

    return (x1*y2 - y1*x2 - x0*(y2 - y1) + y0*(x2 - x1))/2.0;
}

template<typename T>
static_vector<T, 3>
eval_basis(const disk::simplicial_mesh<T, 2>& msh,
           const typename disk::simplicial_mesh<T, 2>::cell& cl,
           const typename disk::simplicial_mesh<T, 2>::point_type& pt)
{
    static_vector<T, 3> ret;
    
    auto meas = cfem_measure(msh, cl);

    auto pts = points(msh, cl);

    auto C = 1./(2.*meas);
    auto x0 = pts[0].x(); auto y0 = pts[0].y();
    auto x1 = pts[1].x(); auto y1 = pts[1].y();
    auto x2 = pts[2].x(); auto y2 = pts[2].y();

    ret(0) = (x1*y2 - y1*x2 - pt.x() * (y2 - y1) + pt.y() * (x2 - x1)) * C;
    ret(1) = (x2*y0 - y2*x0 + pt.x() * (y2 - y0) - pt.y() * (x2 - x0)) * C;
    ret(2) = (x0*y1 - y0*x1 - pt.x() * (y1 - y0) + pt.y() * (x1 - x0)) * C;

    return ret;
}

template<typename T>
static_matrix<T, 2, 3>
eval_basis_grad(const disk::simplicial_mesh<T, 2>& msh,
                const typename disk::simplicial_mesh<T, 2>::cell& cl)
{
    static_matrix<T, 2, 3> ret;
    
    auto meas = cfem_measure(msh, cl);

    auto pts = points(msh, cl);

    auto C = 1./(2.*meas);
    auto x0 = pts[0].x(); auto y0 = pts[0].y();
    auto x1 = pts[1].x(); auto y1 = pts[1].y();
    auto x2 = pts[2].x(); auto y2 = pts[2].y();
    
    ret(0,0) = (y1 - y2) * C;
    ret(0,1) = (y2 - y0) * C;
    ret(0,2) = (y0 - y1) * C;
    ret(1,0) = (x2 - x1) * C;
    ret(1,1) = (x0 - x2) * C;
    ret(1,2) = (x1 - x0) * C;

    return ret;
}

template<typename T>
static_matrix<T, 3, 3>
mass_matrix(const disk::simplicial_mesh<T, 2>& msh,
            const typename disk::simplicial_mesh<T, 2>::cell& cl)
{
    static_matrix<T, 3, 3> ret;
    auto pts = points(msh, cl);
    
    auto meas = measure(msh, cl);
    auto dphi = eval_basis_grad(msh, cl);
    auto stiff = meas * dphi.transpose() * dphi;

    return stiff;
}

template<typename T>
static_matrix<T, 3, 3>
stiffness_matrix(const disk::simplicial_mesh<T, 2>& msh,
                 const typename disk::simplicial_mesh<T, 2>::cell& cl)
{
    static_matrix<T, 3, 3> ret;
    
    auto meas = measure(msh, cl);
    auto dphi = eval_basis_grad(msh, cl);
    auto stiff = meas * dphi.transpose() * dphi;

    return stiff;
}

template<typename T>
static_matrix<T, 3, 3>
stiffness_matrix(const disk::simplicial_mesh<T, 2>& msh,
                 const typename disk::simplicial_mesh<T, 2>::cell& cl,
                 const static_matrix<T, 2, 2>& kappa)
{
    static_matrix<T, 3, 3> ret;
    
    auto meas = measure(msh, cl);

    auto dphi = eval_basis_grad(msh, cl);
    auto stiff = meas * dphi.transpose() * (kappa*dphi);

    return stiff;
}
/*
template<typename T, typename Function>
static_matrix<T, 3, 3>
stiffness_matrix(const disk::simplicial_mesh<T, 2>& msh,
                 const typename disk::simplicial_mesh<T, 2>::cell& cl,
                 Function& kappa)
{
    static_matrix<T, 3, 3> ret;
    
    auto meas = measure(msh, cl);
    auto bar = barycenter(msh, cl);

    auto dphi = eval_basis_grad(msh, cl);
    auto stiff = meas * dphi.transpose() * (kappa(bar)*dphi);

    return stiff;
}
*/

template<typename T, typename Function>
static_vector<T, 3>
make_rhs(const disk::simplicial_mesh<T, 2>& msh,
         const typename disk::simplicial_mesh<T, 2>::cell& cl,
         const Function& f)
{
    static_vector<T, 3> ret;

    auto meas = cfem_measure(msh, cl);
    auto bar = barycenter(msh, cl);

    static_vector<T,3> r = static_vector<T,3>::Ones();

    ret = f(bar) * std::abs(meas) * r / 3.0;

    return ret;
}

} //namespace cfem
} //namespace disk