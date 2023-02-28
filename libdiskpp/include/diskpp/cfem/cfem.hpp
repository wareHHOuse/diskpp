/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
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

#pragma once

#include "diskpp/common/eigen.hpp"

namespace disk { namespace cfem {

template<typename T>
static_vector<T, 3>
eval_basis(const disk::simplicial_mesh<T, 2>& msh,
           const typename disk::simplicial_mesh<T, 2>::cell& cl,
           const typename disk::simplicial_mesh<T, 2>::point_type& pt)
{
    static_vector<T, 3> ret;

    auto pts = points(msh, cl);
    auto x0 = pts[0].x(); auto y0 = pts[0].y();
    auto x1 = pts[1].x(); auto y1 = pts[1].y();
    auto x2 = pts[2].x(); auto y2 = pts[2].y();

    auto m = (x1*y2 - y1*x2 - x0*(y2 - y1) + y0*(x2 - x1));

    ret(0) = (x1*y2 - y1*x2 - pt.x() * (y2 - y1) + pt.y() * (x2 - x1)) / m;
    ret(1) = (x2*y0 - y2*x0 + pt.x() * (y2 - y0) - pt.y() * (x2 - x0)) / m;
    ret(2) = (x0*y1 - y0*x1 - pt.x() * (y1 - y0) + pt.y() * (x1 - x0)) / m;

    return ret;
}

template<typename T>
static_matrix<T, 3, 2>
eval_basis_grad(const disk::simplicial_mesh<T, 2>& msh,
                const typename disk::simplicial_mesh<T, 2>::cell& cl)
{
    static_matrix<T, 3, 2> ret;

    auto pts = points(msh, cl);
    auto x0 = pts[0].x(); auto y0 = pts[0].y();
    auto x1 = pts[1].x(); auto y1 = pts[1].y();
    auto x2 = pts[2].x(); auto y2 = pts[2].y();

    auto m = (x1*y2 - y1*x2 - x0*(y2 - y1) + y0*(x2 - x1));

    ret(0,0) = (y1 - y2) / m;
    ret(1,0) = (y2 - y0) / m;
    ret(2,0) = (y0 - y1) / m;
    ret(0,1) = (x2 - x1) / m;
    ret(1,1) = (x0 - x2) / m;
    ret(2,1) = (x1 - x0) / m;

    return ret;
}

template<typename T>
static_matrix<T, 3, 3>
mass_matrix(const disk::simplicial_mesh<T, 2>& msh,
            const typename disk::simplicial_mesh<T, 2>::cell& cl)
{
    static_matrix<T, 3, 3> ret;
    auto pts = points(msh, cl);

    auto p0 = (pts[0] + pts[1])/2;
    auto p1 = (pts[1] + pts[2])/2;
    auto p2 = (pts[0] + pts[2])/2;

    auto meas = measure(msh, cl);
    auto phi = eval_basis(msh, cl, p0);

    ret = phi * phi.transpose();

    phi = eval_basis(msh, cl, p1);
    ret = ret + phi * phi.transpose();

    phi = eval_basis(msh, cl, p2);
    ret = ret + phi * phi.transpose();

    ret *= meas/3.0;

    return ret;
}

template<typename T>
static_matrix<T, 3, 3>
stiffness_matrix(const disk::simplicial_mesh<T, 2>& msh,
                 const typename disk::simplicial_mesh<T, 2>::cell& cl)
{
    static_matrix<T, 3, 3> ret;

    auto meas = measure(msh, cl);
    auto dphi = eval_basis_grad(msh, cl);
    auto stiff = meas * dphi * dphi.transpose();

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
    auto stiff = meas * dphi * (kappa * dphi.transpose());

    return stiff;
}

template<typename T, typename Function>
static_vector<T, 3>
make_rhs(const disk::simplicial_mesh<T, 2>& msh,
         const typename disk::simplicial_mesh<T, 2>::cell& cl,
         const Function& f)
{
    static_vector<T, 3> ret;

    auto meas = measure(msh, cl);
    auto bar = barycenter(msh, cl);

    static_vector<T,3> r = static_vector<T,3>::Ones();

    ret = f(bar) * meas * r / 3.0;

    return ret;
}

} //namespace cfem
} //namespace disk
