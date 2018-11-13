    /*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
 * Karol Cascavita (C) 2018                     klcascavitam@unal.edu.co
 * Nicolas Pignet  (C) 2018                     nicolas.pignet@enpc.fr
 * Intissar Addali (C) 2018                     intissar.addali@inria.fr
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

#include <iterator>

#include "common/eigen.hpp"
#include "bases/bases.hpp"
#include "quadratures/quadratures.hpp"
#include "mechanics/BoundaryConditions.hpp"

using namespace Eigen;

namespace disk
{

class hho_degree_info
{
    size_t  cell_deg, face_deg, reconstruction_deg, grad_deg;

public:
    hho_degree_info()
        : cell_deg(1), face_deg(1), reconstruction_deg(2), grad_deg(1)
    {}

    explicit hho_degree_info(size_t degree)
        : cell_deg(degree), face_deg(degree), reconstruction_deg(degree+1), grad_deg(degree)
    {}

    hho_degree_info(size_t cd, size_t fd)
    {
        bool c1 = fd > 0  && (cd == fd-1 || cd == fd || cd == fd+1);
        bool c2 = fd == 0 && (cd == fd || cd == fd+1);
        if ( c1 || c2 )
        {
            cell_deg            = cd;
            face_deg            = fd;
            reconstruction_deg  = fd+1;
            grad_deg            = fd;
        }
        else
        {
            std::cout << "Invalid cell degree. Reverting to equal-order" << std::endl;
            cell_deg            = fd;
            face_deg            = fd;
            reconstruction_deg  = fd+1;
            grad_deg            = fd;
        }
    }

    hho_degree_info(size_t cd, size_t fd, size_t gd)
    {
       bool c1 = fd > 0 && (cd == fd - 1 || cd == fd || cd == fd + 1);
       bool c2 = fd == 0 && (cd == fd || cd == fd + 1);
       bool c3 = gd >= fd;
       if (c1 || c2 || c3) {
          cell_deg           = cd;
          face_deg           = fd;
          reconstruction_deg = fd + 1;
          grad_deg           = gd;
       } else {
          std::cout << "Invalid cell degree. Reverting to equal-order" << std::endl;
          cell_deg           = fd;
          face_deg           = fd;
          reconstruction_deg = fd + 1;
          grad_deg           = fd;
       }
    }

    size_t cell_degree() const
    {
        return cell_deg;
    }

    size_t face_degree() const
    {
        return face_deg;
    }

    size_t reconstruction_degree() const
    {
        return reconstruction_deg;
    }

    size_t
    grad_degree() const
    {
       return grad_deg;
    }

    void
    info_degree() const
    {
       std::cout << cell_deg << " " << face_deg << " " << reconstruction_deg << " " << grad_deg
                 << std::endl;
    }
};


template<typename Mesh>
std::pair<   Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>  >
make_hho_scalar_laplacian(const Mesh& msh, const typename Mesh::cell_type& cl,
                          const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    const size_t DIM = Mesh::dimension;

    const auto recdeg = di.reconstruction_degree();
    const auto celdeg = di.cell_degree();
    const auto facdeg = di.face_degree();

    auto cb = make_scalar_monomial_basis(msh, cl, recdeg);

    const auto rbs = scalar_basis_size(recdeg, Mesh::dimension);
    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    const auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);

    const auto num_faces = howmany_faces(msh, cl);

    matrix_type stiff = matrix_type::Zero(rbs, rbs);
    matrix_type gr_lhs = matrix_type::Zero(rbs-1, rbs-1);
    matrix_type gr_rhs = matrix_type::Zero(rbs-1, cbs + num_faces*fbs);

    auto qps = integrate(msh, cl, 2 * (recdeg-1));
    for (auto& qp : qps)
    {
        const auto dphi = cb.eval_gradients(qp.point());
        stiff += qp.weight() * dphi * dphi.transpose();
    }

    gr_lhs = stiff.block(1, 1, rbs-1, rbs-1);
    gr_rhs.block(0, 0, rbs-1, cbs) = stiff.block(1, 0, rbs-1, cbs);

    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = normal(msh, cl, fc);
        auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

        auto qps_f = integrate(msh, fc, recdeg - 1 + std::max(facdeg,celdeg));
        for (auto& qp : qps_f)
        {
            vector_type c_phi_tmp = cb.eval_functions(qp.point());
            vector_type c_phi = c_phi_tmp.head(cbs);
            Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
            Matrix<T, Dynamic, DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
            vector_type f_phi = fb.eval_functions(qp.point());
            gr_rhs.block(0, cbs+i*fbs, rbs-1, fbs) += qp.weight() * (c_dphi * n) * f_phi.transpose();
            gr_rhs.block(0, 0, rbs-1, cbs) -= qp.weight() * (c_dphi * n) * c_phi.transpose();
        }
    }

    matrix_type oper = gr_lhs.ldlt().solve(gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hho_eigval_mass_matrix(const Mesh& msh, const typename Mesh::cell_type& cl,
                            size_t degree)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    auto cb = make_scalar_monomial_basis(msh, cl, degree);
    auto cbs = scalar_basis_size(degree, Mesh::dimension);

    matrix_type data = matrix_type::Zero(cbs, cbs);

    const auto qps = integrate(msh, cl, 2*degree);
    for (auto& qp : qps)
    {
        const auto phi = cb.eval_functions(qp.point());
        data += qp.weight() * phi * phi.transpose();
    }

    return data;
}

template<typename Mesh>
std::pair<Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
          Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>>
make_hho_gradrec_vector(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    const auto celdeg  = di.cell_degree();
    const auto facdeg  = di.face_degree();
    const auto graddeg = di.grad_degree();

    auto cb = make_scalar_monomial_basis(msh, cl, celdeg);
    auto gb = make_vector_monomial_basis(msh, cl, graddeg);

    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    const auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);
    const auto gbs = vector_basis_size(graddeg, Mesh::dimension, Mesh::dimension);

    const auto num_faces = howmany_faces(msh, cl);

    const matrix_type  gr_lhs = make_mass_matrix(msh, cl, gb);
    matrix_type        gr_rhs = matrix_type::Zero(gbs, cbs + num_faces * fbs);

    if(celdeg > 0)
    {
        const auto qps = integrate(msh, cl, celdeg - 1 + facdeg);
        for (auto& qp : qps)
        {
            const auto c_dphi = cb.eval_gradients(qp.point());
            const auto g_phi  = gb.eval_functions(qp.point());

            gr_rhs.block(0, 0, gbs, cbs) += qp.weight() * priv::outer_product(c_dphi, g_phi);
        }
    }

    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = normal(msh, cl, fc);
        auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

        const auto qps_f = integrate(msh, fc, facdeg + std::max(facdeg, celdeg));
        for (auto& qp : qps_f)
        {
            const vector_type c_phi      = cb.eval_functions(qp.point());
            const vector_type f_phi      = fb.eval_functions(qp.point());
            const auto        g_phi      = gb.eval_functions(qp.point());
            const vector_type qp_g_phi_n = qp.weight() * g_phi * n;

            gr_rhs.block(0, cbs + i * fbs, gbs, fbs) += qp_g_phi_n * f_phi.transpose();
            gr_rhs.block(0, 0, gbs, cbs) -= qp_g_phi_n * c_phi.transpose();
        }
    }

    matrix_type oper = gr_lhs.ldlt().solve(gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
std::pair<   Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>  >
make_hho_vector_laplacian(const Mesh& msh, const typename Mesh::cell_type& cl,
                          const hho_degree_info& di)
{
    using  T = typename Mesh::coordinate_type;
    const size_t N = Mesh::dimension;

    const auto recdeg = di.reconstruction_degree();
    const auto celdeg = di.cell_degree();
    const auto facdeg = di.face_degree();

    auto cb = make_vector_monomial_basis(msh, cl, recdeg);

    const auto rbs = vector_basis_size(recdeg, Mesh::dimension, Mesh::dimension);
    const auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
    const auto fbs = vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension);

    const auto num_faces = howmany_faces(msh, cl);

    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;
    typedef Matrix<T, N, N>             gradient_type;
    typedef Matrix<T, Dynamic, N>       function_type;

    matrix_type stiff  = matrix_type::Zero(rbs, rbs);
    matrix_type gr_lhs = matrix_type::Zero(rbs-N, rbs-N);
    matrix_type gr_rhs = matrix_type::Zero(rbs-N, cbs + num_faces*fbs);

    const auto qps = integrate(msh, cl, 2 * (recdeg - 1));
    for (auto& qp : qps)
    {
        const auto dphi = cb.eval_gradients(qp.point());
        stiff += qp.weight() * priv::outer_product(dphi, dphi);
    }

    gr_lhs = stiff.block(N, N, rbs-N, rbs-N);
    gr_rhs.block(0, 0, rbs-N, cbs) = stiff.block(N, 0, rbs-N, cbs);

    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = normal(msh, cl, fc);
        const auto fb = make_vector_monomial_basis(msh, fc, facdeg);

        const auto qps_f = integrate(msh, fc, std::max(facdeg, celdeg) + recdeg - 1);
        for (auto& qp : qps_f)
        {
            function_type   c_phi_tmp = cb.eval_functions(qp.point());
            function_type   c_phi = c_phi_tmp.block(0, 0, cbs, N);

            eigen_compatible_stdvector<gradient_type> c_dphi_tmp = cb.eval_gradients(qp.point());

            auto begin_iter = std::next(c_dphi_tmp.begin(), N);
            eigen_compatible_stdvector<gradient_type> c_dphi;
            c_dphi.resize(rbs - N);
            assert( std::distance(begin_iter, c_dphi_tmp.end()) == c_dphi.size() );
            std::copy(begin_iter, c_dphi_tmp.end(), c_dphi.begin());

            function_type   f_phi = fb.eval_functions(qp.point());

            Matrix<T, Dynamic, N> c_dphi_n = priv::outer_product(c_dphi, n);
            gr_rhs.block(0, cbs + i*fbs, rbs-N, fbs) +=
                    qp.weight() * priv::outer_product(f_phi, c_dphi_n);
            gr_rhs.block(0, 0, rbs-N, cbs) -=
                    qp.weight() * priv::outer_product(c_phi, c_dphi_n);
        }
    }

    matrix_type oper = gr_lhs.ldlt().solve(gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

namespace priv{
    size_t
    nb_lag(const size_t dim)
    {
       size_t lag = 1;
       if (dim == 3) lag = 3;
       return lag;
    }
}
template<typename Mesh>
std::pair<   Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>  >
make_hho_vector_symmetric_laplacian(const Mesh& msh,
                          const typename Mesh::cell_type& cl,
                          const hho_degree_info& di)
{
    using  T = typename Mesh::coordinate_type;
    const size_t N = Mesh::dimension;

    const auto recdeg = di.reconstruction_degree();
    const auto celdeg = di.cell_degree();
    const auto facdeg = di.face_degree();

    const auto cb = make_vector_monomial_basis(msh, cl, recdeg);

    const auto rbs = vector_basis_size(recdeg, N, N);
    const auto cbs = vector_basis_size(celdeg, N, N);
    const auto fbs = vector_basis_size(facdeg, N - 1, N);

    const auto num_faces = howmany_faces(msh, cl);

    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;
    typedef Matrix<T, N, N>             gradient_type;
    typedef Matrix<T, Dynamic, N>       function_type;

    const size_t rbs_ho         = rbs - N;
    const size_t num_total_dofs = cbs + num_faces * fbs;
    const size_t nb_lag         = priv::nb_lag(N);

    matrix_type stiff  = matrix_type::Zero(rbs, rbs);
    matrix_type gr_lhs = matrix_type::Zero(rbs_ho + nb_lag, rbs_ho + nb_lag);
    matrix_type gr_rhs = matrix_type::Zero(rbs_ho + nb_lag, num_total_dofs);

    auto qps = integrate(msh, cl, 2 * (recdeg - 1));
    for (auto& qp : qps)
    {
        const auto dphi = cb.eval_sgradients(qp.point());
        stiff += qp.weight() * priv::outer_product(dphi, dphi);
    }

    gr_lhs.block(0, 0, rbs_ho, rbs_ho) = stiff.block(N, N, rbs_ho, rbs_ho);
    gr_rhs.block(0, 0, rbs_ho, cbs)    = stiff.block(N, 0, rbs_ho, cbs);

    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = normal(msh, cl, fc);
        auto fb = make_vector_monomial_basis(msh, fc, facdeg);

        auto qps_f = integrate(msh, fc, std::max(facdeg, celdeg) + recdeg - 1);
        for (auto& qp : qps_f)
        {
            const function_type c_phi_tmp = cb.eval_functions(qp.point());
            const function_type c_phi     = c_phi_tmp.block(0, 0, cbs, N);

            eigen_compatible_stdvector<gradient_type> c_dphi_tmp = cb.eval_sgradients(qp.point());

            auto begin_iter = std::next(c_dphi_tmp.begin(), N);
            eigen_compatible_stdvector<gradient_type> c_dphi(rbs_ho);
            std::copy(begin_iter, c_dphi_tmp.end(), c_dphi.begin());

            function_type    f_phi = fb.eval_functions(qp.point());
            function_type c_dphi_n = priv::outer_product(c_dphi, n);
            gr_rhs.block(0, cbs + i*fbs, rbs_ho, fbs) +=
                    qp.weight() * priv::outer_product(f_phi, c_dphi_n);
            gr_rhs.block(0, 0, rbs_ho, cbs) -=
                    qp.weight() * priv::outer_product(c_phi, c_dphi_n);
        }
    }

    qps = integrate(msh, cl, recdeg);

    matrix_type rot = matrix_type::Zero(rbs, nb_lag);
    for (auto& qp : qps)
    {
        auto rphi = cb.eval_curls(qp.point());
        rot += qp.weight() * rphi;
    }
    gr_lhs.block(0, rbs_ho, rbs_ho, nb_lag) += rot.bottomLeftCorner(rbs_ho, nb_lag);
    gr_lhs.block(rbs_ho, 0, nb_lag, rbs_ho) += rot.bottomLeftCorner(rbs_ho, nb_lag).transpose();

    // use LU solver because lhs is only symmetric and positive
    matrix_type sol  = gr_lhs.lu().solve(gr_rhs);
    matrix_type oper = sol.block(0,0, rbs_ho, num_total_dofs);
    matrix_type gr   = gr_rhs.block(0,0, rbs_ho, num_total_dofs);
    matrix_type data = gr.transpose() * oper;

    return std::make_pair(oper, data);
}
//#endif

template<typename Mesh>
std::pair<Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
          Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>>
make_hho_sym_gradrec_matrix(const Mesh&                     msh,
                        const typename Mesh::cell_type& cl,
                        const hho_degree_info&          di)
{
   using T        = typename Mesh::coordinate_type;
   typedef Matrix<T, Dynamic, Dynamic> matrix_type;
   typedef Matrix<T, Dynamic, 1>       vector_type;

   const size_t N = Mesh::dimension;

   const auto graddeg = di.grad_degree();
   const auto celdeg  = di.cell_degree();
   const auto facdeg  = di.face_degree();

   const auto gb = make_sym_matrix_monomial_basis(msh, cl, graddeg);
   const auto cb = make_vector_monomial_basis(msh, cl, celdeg);

   const auto gbs = sym_matrix_basis_size(graddeg, Mesh::dimension, Mesh::dimension);
   const auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
   const auto fbs = vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension);

   const auto num_faces = howmany_faces(msh, cl);

   matrix_type gr_lhs = matrix_type::Zero(gbs, gbs);
   matrix_type gr_rhs = matrix_type::Zero(gbs, cbs + num_faces * fbs);

   const size_t dim2 = N * N;

   // this is very costly to build it
   const auto qps = integrate(msh, cl, 2 * graddeg);

   size_t dec = 0;
    if (N == 3)
        dec = 6;
    else if (N == 2)
        dec = 3;
    else
        std::logic_error("Expected 3 >= dim > 1");

    for (auto& qp : qps)
    {
        const auto gphi = gb.eval_functions(qp.point());

        for (size_t j = 0; j < gbs; j++)
        {
            const auto qp_gphi_j = priv::inner_product(qp.weight(), gphi[j]);
            for (size_t i = j; i < gbs; i += dec)
                gr_lhs(i, j) += priv::inner_product(gphi[i], qp_gphi_j);
        }
    }

   // upper part
    for (size_t j = 0; j < gbs; j++)
        for (size_t i = 0; i < j; i++)
            gr_lhs(i, j) = gr_lhs(j, i);

    // compute rhs
    if(celdeg)
    {
        const auto qpc = integrate(msh, cl, graddeg + celdeg - 1);
        for (auto& qp : qpc)
        {
            const auto gphi = gb.eval_functions(qp.point());
            const auto dphi = cb.eval_sgradients(qp.point());

            for (size_t j = 0; j < cbs; j++)
            {
                const auto qp_dphi_j = priv::inner_product(qp.weight(), dphi[j]);
                for (size_t i = 0; i < gbs; i++)
                    gr_rhs(i, j) += priv::inner_product(gphi[i], qp_dphi_j);
            }
        } // end qp
    }


    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = normal(msh, cl, fc);
        const auto fb = make_vector_monomial_basis(msh, fc, facdeg);

        const auto qps_f = integrate(msh, fc, graddeg + std::max(celdeg, facdeg));
        for (auto& qp : qps_f)
        {
            const auto cphi = cb.eval_functions(qp.point());
            const auto gphi = gb.eval_functions(qp.point());
            const auto fphi = fb.eval_functions(qp.point());

            const auto gphi_n = priv::outer_product(gphi, n);
            gr_rhs.block(0, cbs + i * fbs, gbs, fbs) +=
                                qp.weight() * priv::outer_product(fphi, gphi_n);
            gr_rhs.block(0, 0, gbs, cbs) -= qp.weight() * priv::outer_product(cphi, gphi_n);
        }
    }

    matrix_type oper = gr_lhs.ldlt().solve(gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
std::pair<   Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>  >
make_hho_divergence_reconstruction(const Mesh& msh, const typename Mesh::cell_type& cl,
                                   const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto celdeg = di.cell_degree();
    const auto facdeg = di.face_degree();
    const auto recdeg = di.face_degree();

    auto cbas_s = make_scalar_monomial_basis(msh, cl, recdeg);

    const auto rbs = scalar_basis_size(recdeg, Mesh::dimension);
    const auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
    const auto fbs = vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension);

    const auto num_faces = howmany_faces(msh, cl);

    const auto dr_lhs = make_mass_matrix(msh, cl, cbas_s);
    const auto dr_rhs = make_hho_divergence_reconstruction_stokes_rhs(msh, cl, di);

    assert(dr_lhs.rows() == rbs && dr_lhs.cols() == rbs);
    assert(dr_rhs.rows() == rbs && dr_rhs.cols() == cbs + num_faces * fbs);

    matrix_type oper = dr_lhs.ldlt().solve(dr_rhs);
    matrix_type data = dr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hho_divergence_reconstruction_stokes_rhs(const Mesh& msh, const typename Mesh::cell_type& cl,
                                   const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto celdeg = di.cell_degree();
    const auto facdeg = di.face_degree();
    const auto recdeg = di.face_degree();

    const auto cbas_v = make_vector_monomial_basis(msh, cl, celdeg);
    const auto cbas_s = make_scalar_monomial_basis(msh, cl, recdeg);

    const auto rbs = scalar_basis_size(recdeg, Mesh::dimension);
    const auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
    const auto fbs = vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension);

    const auto num_faces = howmany_faces(msh, cl);

    matrix_type dr_rhs = matrix_type::Zero(rbs, cbs + num_faces*fbs);

    if(recdeg > 0)
    {
        const auto qps = integrate(msh, cl, celdeg + recdeg - 1);
        for (auto& qp : qps)
        {
            const auto s_dphi = cbas_s.eval_gradients(qp.point());
            const auto v_phi  = cbas_v.eval_functions(qp.point());

            dr_rhs.block(0, 0, rbs, cbs) -= qp.weight() * priv::outer_product(v_phi, s_dphi);
        }
    }

    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc     = fcs[i];
        const auto n      = normal(msh, cl, fc);
        const auto fbas_v = make_vector_monomial_basis(msh, fc, facdeg);

        const auto qps_f = integrate(msh, fc, facdeg + recdeg);
        for (auto& qp : qps_f)
        {
            const auto s_phi = cbas_s.eval_functions(qp.point());
            const auto f_phi = fbas_v.eval_functions(qp.point());

            const Matrix<T, Dynamic, Mesh::dimension> s_phi_n = (s_phi * n.transpose());
            dr_rhs.block(0, cbs + i * fbs, rbs, fbs) += qp.weight() * priv::outer_product(f_phi, s_phi_n);
        }
    }

    return dr_rhs;
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hdg_scalar_stabilization(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    const auto celdeg = di.cell_degree();
    const auto facdeg = di.face_degree();

    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    const auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);

    const auto num_faces = howmany_faces(msh, cl);
    const auto total_dofs = cbs + num_faces * fbs;

    matrix_type       data = matrix_type::Zero(total_dofs, total_dofs);
    const matrix_type If   = matrix_type::Identity(fbs, fbs);

    auto cb = make_scalar_monomial_basis(msh, cl, celdeg);
    const auto fcs = faces(msh, cl);

    for (size_t i = 0; i < num_faces; i++)
    {
        const auto fc = fcs[i];
        const auto h  = diameter(msh, fc);
        auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

        matrix_type oper  = matrix_type::Zero(fbs, total_dofs);
        matrix_type tr    = matrix_type::Zero(fbs, total_dofs);
        matrix_type mass  = matrix_type::Zero(fbs, fbs);
        matrix_type trace = matrix_type::Zero(fbs, cbs);

        oper.block(0, cbs + i  * fbs, fbs, fbs) = -If;

        const auto qps = integrate(msh, fc, facdeg + celdeg);
        for (auto& qp : qps)
        {
            const auto c_phi = cb.eval_functions(qp.point());
            const auto f_phi = fb.eval_functions(qp.point());

            assert(c_phi.rows() == cbs);
            assert(f_phi.rows() == fbs);
            assert(c_phi.cols() == f_phi.cols());

            mass += qp.weight() * f_phi * f_phi.transpose();
            trace += qp.weight() * f_phi * c_phi.transpose();
        }

        tr.block(0, cbs + i * fbs, fbs, fbs) = -mass;
        tr.block(0, 0, fbs, cbs) = trace;

        oper.block(0, 0, fbs, cbs) = mass.ldlt().solve(trace);
        data += oper.transpose() * tr * (1./h);
    }

    return data;
}

template<typename T, int M, int N>
T estimate_conditioning(const Matrix<T,M,N>& m)
{
    Eigen::JacobiSVD< Matrix<T,M,N> > svd(m);
    T sigma_max = svd.singularValues()(0);
    T sigma_min = svd.singularValues()(svd.singularValues().size()-1);
    T cond =  sigma_max / sigma_min;
    return cond;
}


template<typename Mesh, typename T>
auto
diffusion_static_condensation_compute(const Mesh& msh,
        const typename Mesh::cell_type& cl, const hho_degree_info hdi,
        const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& local_mat,
        const typename Eigen::Matrix<T, Eigen::Dynamic, 1>& cell_rhs)
{
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    auto cell_degree = hdi.cell_degree();
    auto face_degree = hdi.face_degree();
    auto num_cell_dofs = scalar_basis_size(cell_degree, Mesh::dimension);
    auto num_face_dofs = scalar_basis_size(face_degree, Mesh::dimension-1);

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    assert(local_mat.rows() == local_mat.cols());
    assert(local_mat.cols() == num_cell_dofs + num_faces*num_face_dofs);
    assert(cell_rhs.rows() == num_cell_dofs);

    size_t cell_size = num_cell_dofs;
    size_t face_size = num_face_dofs * num_faces;

    matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
    matrix_type K_TF = local_mat.topRightCorner(cell_size, face_size);
    matrix_type K_FT = local_mat.bottomLeftCorner(face_size, cell_size);
    matrix_type K_FF = local_mat.bottomRightCorner(face_size, face_size);

    assert(K_TT.cols() == cell_size);
    assert(K_TT.cols() + K_TF.cols() == local_mat.cols());
    assert(K_TT.rows() + K_FT.rows() == local_mat.rows());
    assert(K_TF.rows() + K_FF.rows() == local_mat.rows());
    assert(K_FT.cols() + K_FF.cols() == local_mat.cols());

    auto K_TT_ldlt = K_TT.ldlt();
    matrix_type AL = K_TT_ldlt.solve(K_TF);
    vector_type bL = K_TT_ldlt.solve(cell_rhs);

    matrix_type AC = K_FF - K_FT * AL;
    vector_type bC = /* no projection on faces, eqn. 26*/ - K_FT * bL;

    return std::make_pair(AC, bC);
}

template<typename Mesh, typename T>
auto
diffusion_static_condensation_compute_full(const Mesh& msh,
        const typename Mesh::cell_type& cl, const hho_degree_info hdi,
        const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& local_mat,
        const typename Eigen::Matrix<T, Eigen::Dynamic, 1>& rhs)
{
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    auto cell_degree = hdi.cell_degree();
    auto face_degree = hdi.face_degree();
    auto num_cell_dofs = scalar_basis_size(cell_degree, Mesh::dimension);
    auto num_face_dofs = scalar_basis_size(face_degree, Mesh::dimension-1);

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    assert(local_mat.rows() == local_mat.cols());
    assert(local_mat.cols() == num_cell_dofs + num_faces*num_face_dofs);
    assert(rhs.rows() == num_cell_dofs + num_faces*num_face_dofs);

    size_t cell_size = num_cell_dofs;
    size_t face_size = num_face_dofs * num_faces;

    matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
    matrix_type K_TF = local_mat.topRightCorner(cell_size, face_size);
    matrix_type K_FT = local_mat.bottomLeftCorner(face_size, cell_size);
    matrix_type K_FF = local_mat.bottomRightCorner(face_size, face_size);

    assert(K_TT.cols() == cell_size);
    assert(K_TT.cols() + K_TF.cols() == local_mat.cols());
    assert(K_TT.rows() + K_FT.rows() == local_mat.rows());
    assert(K_TF.rows() + K_FF.rows() == local_mat.rows());
    assert(K_FT.cols() + K_FF.cols() == local_mat.cols());

    vector_type cell_rhs = rhs.block(0, 0, num_cell_dofs, 1);
    vector_type face_rhs = rhs.block(num_cell_dofs, 0, num_faces*num_face_dofs, 1);

    auto K_TT_ldlt = K_TT.ldlt();
    matrix_type AL = K_TT_ldlt.solve(K_TF);
    vector_type bL = K_TT_ldlt.solve(cell_rhs);

    matrix_type AC = K_FF - K_FT * AL;
    vector_type bC = /* no projection on faces, eqn. 26*/face_rhs - K_FT * bL;

    return std::make_pair(AC, bC);
}

template<typename Mesh, typename T>
auto
diffusion_static_condensation_compute_vector(const Mesh& msh,
        const typename Mesh::cell_type& cl, const hho_degree_info hdi,
        const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& local_mat,
        const typename Eigen::Matrix<T, Eigen::Dynamic, 1>& cell_rhs)
{
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    auto cell_degree = hdi.cell_degree();
    auto face_degree = hdi.face_degree();
    auto num_cell_dofs = vector_basis_size(cell_degree, Mesh::dimension, Mesh::dimension);
    auto num_face_dofs = vector_basis_size(face_degree, Mesh::dimension-1, Mesh::dimension);

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    assert(local_mat.rows() == local_mat.cols());
    assert(local_mat.cols() == num_cell_dofs + num_faces*num_face_dofs);
    assert(cell_rhs.rows() == num_cell_dofs);

    size_t cell_size = num_cell_dofs;
    size_t face_size = num_face_dofs * num_faces;

    matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
    matrix_type K_TF = local_mat.topRightCorner(cell_size, face_size);
    matrix_type K_FT = local_mat.bottomLeftCorner(face_size, cell_size);
    matrix_type K_FF = local_mat.bottomRightCorner(face_size, face_size);

    assert(K_TT.cols() == cell_size);
    assert(K_TT.cols() + K_TF.cols() == local_mat.cols());
    assert(K_TT.rows() + K_FT.rows() == local_mat.rows());
    assert(K_TF.rows() + K_FF.rows() == local_mat.rows());
    assert(K_FT.cols() + K_FF.cols() == local_mat.cols());

    auto K_TT_ldlt = K_TT.ldlt();
    matrix_type AL = K_TT_ldlt.solve(K_TF);
    vector_type bL = K_TT_ldlt.solve(cell_rhs);

    matrix_type AC = K_FF - K_FT * AL;
    vector_type bC = /* no projection on faces, eqn. 26*/ - K_FT * bL;

    return std::make_pair(AC, bC);
}

template<typename Mesh, typename T>
auto
diffusion_static_condensation_compute_vector_full(const Mesh& msh,
        const typename Mesh::cell_type& cl, const hho_degree_info hdi,
        const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& local_mat,
        const typename Eigen::Matrix<T, Eigen::Dynamic, 1>& rhs)
{
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    auto cell_degree = hdi.cell_degree();
    auto face_degree = hdi.face_degree();
    auto num_cell_dofs = vector_basis_size(cell_degree, Mesh::dimension, Mesh::dimension);
    auto num_face_dofs = vector_basis_size(face_degree, Mesh::dimension-1, Mesh::dimension);

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    assert(local_mat.rows() == local_mat.cols());
    assert(local_mat.cols() == num_cell_dofs + num_faces*num_face_dofs);

    size_t cell_size = num_cell_dofs;
    size_t face_size = num_face_dofs * num_faces;

    matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
    matrix_type K_TF = local_mat.topRightCorner(cell_size, face_size);
    matrix_type K_FT = local_mat.bottomLeftCorner(face_size, cell_size);
    matrix_type K_FF = local_mat.bottomRightCorner(face_size, face_size);

    assert(K_TT.cols() == cell_size);
    assert(K_TT.cols() + K_TF.cols() == local_mat.cols());
    assert(K_TT.rows() + K_FT.rows() == local_mat.rows());
    assert(K_TF.rows() + K_FF.rows() == local_mat.rows());
    assert(K_FT.cols() + K_FF.cols() == local_mat.cols());
    assert((rhs.rows() == cell_size) || (rhs.rows() == cell_size + num_faces*num_face_dofs));

    vector_type cell_rhs = rhs.block(0, 0, num_cell_dofs, 1);
    vector_type face_rhs = vector_type::Zero(face_size);
    if(rhs.cols() ==  num_cell_dofs + num_faces*num_face_dofs)
        face_rhs = rhs.block(num_cell_dofs, 0, num_faces*num_face_dofs, 1);

    auto K_TT_ldlt = K_TT.ldlt();
    matrix_type AL = K_TT_ldlt.solve(K_TF);
    vector_type bL = K_TT_ldlt.solve(cell_rhs);

    matrix_type AC = K_FF - K_FT * AL;
    vector_type bC = /* no projection on faces, eqn. 26*/face_rhs - K_FT * bL;

    return std::make_pair(AC, bC);
}

template<typename Mesh>
std::pair<Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>,
Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, 1>>
diffusion_static_condensation_compute_alg(const Mesh& msh,
        const typename Mesh::cell_type& cl, const hho_degree_info hdi,
        const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>& local_mat,
        const Matrix<typename Mesh::coordinate_type, Dynamic, 1>& rhs)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    auto cell_degree = hdi.cell_degree();
    auto face_degree = hdi.face_degree();
    auto num_cell_dofs = scalar_basis_size(cell_degree, Mesh::dimension);
    auto num_face_dofs = scalar_basis_size(face_degree, Mesh::dimension-1);

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    assert(local_mat.rows() == local_mat.cols());
    assert(local_mat.cols() == num_cell_dofs + num_faces*num_face_dofs);
    assert(rhs.rows() == num_cell_dofs  + num_faces*num_face_dofs);

    size_t cell_size = num_cell_dofs;
    size_t face_size = num_face_dofs * num_faces;

    matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
    matrix_type K_TF = local_mat.topRightCorner(cell_size, face_size);
    matrix_type K_FT = local_mat.bottomLeftCorner(face_size, cell_size);
    matrix_type K_FF = local_mat.bottomRightCorner(face_size, face_size);

    assert(K_TT.cols() == cell_size);
    assert(K_TT.cols() + K_TF.cols() == local_mat.cols());
    assert(K_TT.rows() + K_FT.rows() == local_mat.rows());
    assert(K_TF.rows() + K_FF.rows() == local_mat.rows());
    assert(K_FT.cols() + K_FF.cols() == local_mat.cols());

    vector_type cell_rhs = rhs.block(0, 0, num_cell_dofs, 1);
    vector_type face_rhs = rhs.block(num_cell_dofs, 0, num_faces*num_face_dofs, 1);

    auto K_TT_ldlt = K_TT.ldlt();
    matrix_type AL = K_TT_ldlt.solve(K_TF);
    vector_type bL = K_TT_ldlt.solve(cell_rhs);

    matrix_type AC = K_FF - K_FT * AL;
    vector_type bC = /* no projection on faces, eqn. 26*/ face_rhs - K_FT * bL ;

    return std::make_pair(AC, bC);
}

template<typename Mesh, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
diffusion_static_condensation_recover(const Mesh& msh,
    const typename Mesh::cell_type& cl, const hho_degree_info hdi,
    const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& local_mat,
    const typename Eigen::Matrix<T, Eigen::Dynamic, 1>& cell_rhs,
    const typename Eigen::Matrix<T, Eigen::Dynamic, 1>& solF)
{
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    auto cell_degree = hdi.cell_degree();
    auto face_degree = hdi.face_degree();
    auto num_cell_dofs = scalar_basis_size(cell_degree, Mesh::dimension);
    auto num_face_dofs = scalar_basis_size(face_degree, Mesh::dimension-1);

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    size_t cell_size        = num_cell_dofs;
    size_t all_faces_size   = num_face_dofs * num_faces;

    vector_type ret( cell_size + all_faces_size );

    matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
    matrix_type K_TF = local_mat.topRightCorner(cell_size, all_faces_size);

    vector_type solT = K_TT.ldlt().solve(cell_rhs - K_TF*solF);

    ret.head(cell_size)         = solT;
    ret.tail(all_faces_size)    = solF;

    return ret;
}


template<typename Mesh, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
diffusion_static_condensation_recover_vector(const Mesh& msh,
    const typename Mesh::cell_type& cl, const hho_degree_info hdi,
    const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& local_mat,
    const typename Eigen::Matrix<T, Eigen::Dynamic, 1>& cell_rhs,
    const typename Eigen::Matrix<T, Eigen::Dynamic, 1>& solF)
{
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    auto cell_degree = hdi.cell_degree();
    auto face_degree = hdi.face_degree();
    auto num_cell_dofs = vector_basis_size(cell_degree, Mesh::dimension, Mesh::dimension);
    auto num_face_dofs = vector_basis_size(face_degree, Mesh::dimension-1, Mesh::dimension);

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    size_t cell_size        = num_cell_dofs;
    size_t all_faces_size   = num_face_dofs * num_faces;

    vector_type ret( cell_size + all_faces_size );

    matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
    matrix_type K_TF = local_mat.topRightCorner(cell_size, all_faces_size);

    vector_type solT = K_TT.ldlt().solve(cell_rhs - K_TF*solF);

    ret.head(cell_size)         = solT;
    ret.tail(all_faces_size)    = solF;

    return ret;
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hho_scalar_stabilization(const Mesh&                                                     msh,
                              const typename Mesh::cell_type&                                 cl,
                              const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>& reconstruction,
                              const hho_degree_info&                                          di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto recdeg = di.reconstruction_degree();
    const auto celdeg = di.cell_degree();
    const auto facdeg = di.face_degree();

    const auto rbs = scalar_basis_size(recdeg, Mesh::dimension);
    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    const auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);

    auto cb = make_scalar_monomial_basis(msh, cl, recdeg);

    const matrix_type mass_mat = make_mass_matrix(msh, cl, cb);

    // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

    // Step 1: compute \pi_T^k p_T^k v (third term).
    const matrix_type M1    = mass_mat.block(0, 0, cbs, cbs);
    const matrix_type M2    = mass_mat.block(0, 1, cbs, rbs - 1);
    matrix_type       proj1 = -M1.ldlt().solve(M2 * reconstruction);

    assert(M2.cols() == reconstruction.rows());

    // Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
    proj1.block(0, 0, cbs, cbs) += matrix_type::Identity(cbs, cbs);

    const auto fcs       = faces(msh, cl);
    const auto num_faces = fcs.size();

    matrix_type data = matrix_type::Zero(cbs + num_faces * fbs, cbs + num_faces * fbs);

    // Step 3: project on faces (eqn. 21)
    for (size_t face_i = 0; face_i < num_faces; face_i++)
    {
        const auto fc = fcs[face_i];
        const auto hf = diameter(msh, fc);
        auto       fb = make_scalar_monomial_basis(msh, fc, facdeg);

        matrix_type face_mass_matrix  = matrix_type::Zero(fbs, fbs);
        matrix_type face_trace_matrix = matrix_type::Zero(fbs, rbs);

        const auto face_quadpoints = integrate(msh, fc, recdeg + facdeg);
        for (auto& qp : face_quadpoints)
        {
            const auto        f_phi   = fb.eval_functions(qp.point());
            const auto        c_phi   = cb.eval_functions(qp.point());
            const matrix_type q_f_phi = qp.weight() * f_phi;
            face_mass_matrix += q_f_phi * f_phi.transpose();
            face_trace_matrix += q_f_phi * c_phi.transpose();
        }

        LLT<matrix_type> piKF;
        piKF.compute(face_mass_matrix);

        // Step 3a: \pi_F^k( v_F - p_T^k v )
        const matrix_type MR1 = face_trace_matrix.block(0, 1, fbs, rbs - 1);

        matrix_type       proj2 = piKF.solve(MR1 * reconstruction);
        proj2.block(0, cbs + face_i * fbs, fbs, fbs) -= matrix_type::Identity(fbs, fbs);

        // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
        const matrix_type MR2   = face_trace_matrix.block(0, 0, fbs, cbs);
        const matrix_type proj3 = piKF.solve(MR2 * proj1);
        const matrix_type BRF   = proj2 + proj3;

        data += BRF.transpose() * face_mass_matrix * BRF / hf;
    }

    return data;
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hho_scalar_stabilization_2(const Mesh& msh, const typename Mesh::cell_type& cl,
                                const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic> reconstruction,
                                const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    auto recdeg = di.reconstruction_degree();
    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    auto rbs = scalar_basis_size(recdeg, Mesh::dimension);
    auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);

    auto cb = make_scalar_monomial_basis(msh, cl, recdeg);

    matrix_type mass_mat = matrix_type::Zero(rbs, rbs);
    auto cell_quadpoints = integrate(msh, cl, 2*recdeg);
    for (auto& qp : cell_quadpoints)
    {
        auto c_phi = cb.eval_functions(qp.point());
        mass_mat += qp.weight() * c_phi * c_phi.transpose();
    }

    // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

    //Step 1: compute \pi_T^k p_T^k v (third term).
    matrix_type M1 = mass_mat.block(0, 0, cbs, cbs);
    matrix_type M2 = mass_mat.block(0, 1, cbs, rbs-1);
    matrix_type proj1 = -M1.ldlt().solve(M2*reconstruction);

    //Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
    matrix_type I_T = matrix_type::Identity(cbs, cbs);
    proj1.block(0, 0, cbs, cbs) += I_T;

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    matrix_type data = matrix_type::Zero(cbs+num_faces*fbs, cbs+num_faces*fbs);

    auto h = diameter(msh, cl);

    // Step 3: project on faces (eqn. 21)
    for (size_t face_i = 0; face_i < num_faces; face_i++)
    {
        auto fc = fcs[face_i];
        auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

        matrix_type face_mass_matrix    = matrix_type::Zero(fbs, fbs);
        matrix_type face_trace_matrix   = matrix_type::Zero(fbs, rbs);

        auto face_quadpoints = integrate(msh, fc, 2*recdeg);
        for (auto& qp : face_quadpoints)
        {
            auto f_phi = fb.eval_functions(qp.point());
            auto c_phi = cb.eval_functions(qp.point());
            vector_type q_f_phi = qp.weight() * f_phi;
            face_mass_matrix += q_f_phi * f_phi.transpose();
            face_trace_matrix += q_f_phi * c_phi.transpose();
        }

        LLT<matrix_type> piKF;
        piKF.compute(face_mass_matrix);

        // Step 3a: \pi_F^k( v_F - p_T^k v )
        matrix_type MR1 = face_trace_matrix.block(0, 1, fbs, rbs-1);
        matrix_type A = MR1 * reconstruction;

        matrix_type proj2 = piKF.solve(A);
        matrix_type I_F = matrix_type::Identity(fbs, fbs);
        proj2.block(0, cbs+face_i*fbs, fbs, fbs) -= I_F;
        A.block(0, cbs+face_i*fbs, fbs, fbs) -= face_mass_matrix;

        // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
        matrix_type MR2 = face_trace_matrix.block(0, 0, fbs, cbs);
        matrix_type B = MR2 * proj1;
        matrix_type proj3 = piKF.solve(B);
        matrix_type BRF = proj2 + proj3;

        data += BRF.transpose() * (A+B) / h;
    }

    return data;
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hho_vector_stabilization(const Mesh&                                                     msh,
                              const typename Mesh::cell_type&                                 cl,
                              const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>& reconstruction,
                              const hho_degree_info&                                          di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    const auto recdeg = di.reconstruction_degree();
    const auto celdeg = di.cell_degree();
    const auto facdeg = di.face_degree();

    const auto rbs = vector_basis_size(recdeg, Mesh::dimension, Mesh::dimension);
    const auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
    const auto fbs = vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension);

    size_t N = Mesh::dimension;

    auto cb = make_vector_monomial_basis(msh, cl, recdeg);

    const matrix_type mass_mat = make_mass_matrix(msh, cl, cb);

    // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

    // Step 1: compute \pi_T^k p_T^k v (third term).
    const matrix_type M1    = mass_mat.block(0, 0, cbs, cbs);
    const matrix_type M2    = mass_mat.block(0, N, cbs, rbs - N);
    matrix_type       proj1 = -M1.ldlt().solve(M2 * reconstruction);

    assert(M2.cols() == reconstruction.rows());

    // Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
    proj1.block(0, 0, cbs, cbs) += matrix_type::Identity(cbs, cbs);

    const auto fcs       = faces(msh, cl);
    const auto num_faces = fcs.size();

    matrix_type data = matrix_type::Zero(cbs + num_faces * fbs, cbs + num_faces * fbs);

    // Step 3: project on faces (eqn. 21)
    for (size_t face_i = 0; face_i < num_faces; face_i++)
    {
        const auto fc = fcs[face_i];
        const auto hf = diameter(msh, fc);
        auto       fb = make_vector_monomial_basis(msh, fc, facdeg);

        matrix_type face_mass_matrix  = matrix_type::Zero(fbs, fbs);
        matrix_type face_trace_matrix = matrix_type::Zero(fbs, rbs);

        const auto face_quadpoints = integrate(msh, fc, recdeg + facdeg);
        for (auto& qp : face_quadpoints)
        {
            const auto        f_phi   = fb.eval_functions(qp.point());
            const auto        c_phi   = cb.eval_functions(qp.point());
            const Matrix<T, Dynamic, Mesh::dimension> q_f_phi = qp.weight() * f_phi;
            face_mass_matrix += priv::outer_product(f_phi, q_f_phi);
            face_trace_matrix += priv::outer_product(c_phi, q_f_phi);
        }

        LLT<matrix_type> piKF;
        piKF.compute(face_mass_matrix);

        // Step 3a: \pi_F^k( v_F - p_T^k v )
        const matrix_type MR1 = face_trace_matrix.block(0, N, fbs, rbs - N);

        matrix_type       proj2 = piKF.solve(MR1 * reconstruction);
        proj2.block(0, cbs + face_i * fbs, fbs, fbs) -= matrix_type::Identity(fbs, fbs);

        // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
        const matrix_type MR2   = face_trace_matrix.block(0, 0, fbs, cbs);
        const matrix_type proj3 = piKF.solve(MR2 * proj1);
        const matrix_type BRF   = proj2 + proj3;

        data += BRF.transpose() * face_mass_matrix * BRF / hf;
    }

    return data;
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
project_function(const Mesh& msh, const typename Mesh::cell_type& cl,
                 const hho_degree_info& hdi,
                 const scalar_rhs_function<Mesh>& f, size_t di = 0)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto num_faces = howmany_faces(msh, cl);

    vector_type ret = vector_type::Zero(cbs+num_faces*fbs);

    auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
    matrix_type cell_mm = make_mass_matrix(msh, cl, cb, di);
    vector_type cell_rhs = make_rhs(msh, cl, cb, f, di);
    ret.block(0, 0, cbs, 1) = cell_mm.ldlt().solve(cell_rhs);

    auto fcs = faces(msh, cl);
    for (size_t i = 0; i < num_faces; i++)
    {
        auto fc = fcs[i];
        auto fb = make_scalar_monomial_basis(msh, fc, hdi.face_degree());
        matrix_type face_mm = make_mass_matrix(msh, fc, fb, di);
        vector_type face_rhs = make_rhs(msh, fc, fb, f, di);
        ret.block(cbs+i*fbs, 0, fbs, 1) = face_mm.ldlt().solve(face_rhs);
    }

    return ret;
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
project_function(const Mesh& msh, const typename Mesh::cell_type& cl,
                 const hho_degree_info& hdi,
                 const vector_rhs_function<Mesh>& f, size_t di = 0)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    auto cbs = vector_basis_size(hdi.cell_degree(), Mesh::dimension, Mesh::dimension);
    auto fbs = vector_basis_size(hdi.face_degree(), Mesh::dimension-1, Mesh::dimension);
    auto num_faces = howmany_faces(msh, cl);

    vector_type ret = vector_type::Zero(cbs+num_faces*fbs);

    auto cb = make_vector_monomial_basis(msh, cl, hdi.cell_degree());
    matrix_type cell_mm = make_mass_matrix(msh, cl, cb, di);
    vector_type cell_rhs = make_rhs(msh, cl, cb, f, di);
    ret.block(0, 0, cbs, 1) = cell_mm.ldlt().solve(cell_rhs);

    auto fcs = faces(msh, cl);
    for (size_t i = 0; i < num_faces; i++)
    {
        ret.segment(cbs + i * fbs, fbs) = project_function(msh, fcs[i], hdi, f, di);
    }

    return ret;
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
project_function(const Mesh&                      msh,
                 const typename Mesh::face_type&  fc,
                 const hho_degree_info&           hdi,
                 const vector_rhs_function<Mesh>& f,
                 size_t                           di = 0)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1> vector_type;

    auto fbs       = vector_basis_size(hdi.face_degree(), Mesh::dimension - 1, Mesh::dimension);

    vector_type ret = vector_type::Zero(fbs);

    auto  fb             = make_vector_monomial_basis(msh, fc, hdi.face_degree());
    matrix_type face_mm  = make_mass_matrix(msh, fc, fb, di);
    vector_type face_rhs = make_rhs(msh, fc, fb, f, di);
    ret.block(0, 0, fbs, 1)  = face_mm.ldlt().solve(face_rhs);

   return ret;
}

#if 0
template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
project_tensor( const Mesh& msh, const typename Mesh::cell_type& cl,
                const hho_degree_info&           hdi,
                const vector_rhs_function<Mesh>& f, size_t di = 0)
{
    using T = typename Mesh::coordinate_type;

    auto cb  =  make_vector_monomial_basis(msh, cl, hdi.cell_degree());
    auto cbs =  vector_basis_size(hdi.cell_degree(), Mesh::dimension, Mesh::dimension);

    matrix_type cell_mm = make_mass_matrix(msh, cl, cb, di);
    vector_type rhs = vector_type::Zero(cbs);
    //vector_type cell_rhs = make_rhs(msh, cl, cb, f, di);

    auto qps = integrate(msh, cl, 2*(hdi.cell_degree()+di));

    for (auto& qp : qps)
    {
        auto phi = cb.eval_functions(qp.point());
        auto rhs_func = f(qp.point());
        for(auto i = 0; i < cbs; i++)
            rhs(i) += qp.weight() * phi(i) * rhs_func;
    }

    return cell_mm.llt().solve(rhs);
}
template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
project_tensor( const Mesh& msh, const typename Mesh::cell_type& cl,
                const hho_degree_info&           hdi,
                const matrix_rhs_function<Mesh>& f, size_t di = 0)
{
    using T = typename Mesh::coordinate_type;

    auto cb  =  make_matrix_monomial_basis(msh, cl, hdi.cell_degree());
    auto cbs =  matrix_basis_size(hdi.cell_degree(), Mesh::dimension, Mesh::dimension);

    matrix_type cell_mm = make_mass_matrix(msh, cl, cb, di);
    vector_type rhs = vector_type::Zero(cbs);
    //vector_type cell_rhs = make_rhs(msh, cl, cb, f, di);

    auto qps = integrate(msh, cl, 2*(hdi.cell_degree()+di));

    for (auto& qp : qps)
    {
        auto phi = cb.eval_functions(qp.point());
        auto rhs_func = f(qp.point());
        rhs += qp.weight() * priv::outer_product(phi, rhs_func);
    }

    return cell_mm.llt().solve(rhs);
}
#endif
namespace priv
{

template<typename Mesh>
size_t
offset(const Mesh& msh, const typename Mesh::cell_type& cl)
{
    auto itor = std::lower_bound(msh.cells_begin(), msh.cells_end(), cl);
    if ( itor == msh.cells_end() )
        throw std::logic_error("Cell not found: this is likely a bug.");

    return std::distance(msh.cells_begin(), itor);
}

template<typename Mesh>
size_t
offset(const Mesh& msh, const typename Mesh::face_type& fc)
{
    auto itor = std::lower_bound(msh.faces_begin(), msh.faces_end(), fc);
    if ( itor == msh.faces_end() )
        throw std::logic_error("Face not found: this is likely a bug.");

    return std::distance(msh.faces_begin(), itor);
}


} // priv


template<typename Mesh>
class diffusion_condensed_assembler
{
    using T = typename Mesh::coordinate_type;
    std::vector<size_t>                 compress_table;
    std::vector<size_t>                 expand_table;

    hho_degree_info                     di;

    std::vector< Triplet<T> >           triplets;

    size_t      num_all_faces, num_dirichlet_faces, num_other_faces;

    class assembly_index
    {
        size_t  idx;
        bool    assem;

    public:
        assembly_index(size_t i, bool as)
            : idx(i), assem(as)
        {}

        operator size_t() const
        {
            if (!assem)
                throw std::logic_error("Invalid assembly_index");

            return idx;
        }

        bool assemble() const
        {
            return assem;
        }

        friend std::ostream& operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };

public:
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    SparseMatrix<T>         LHS;
    vector_type   RHS;

    diffusion_condensed_assembler(const Mesh& msh, hho_degree_info hdi)
        : di(hdi)
    {
        auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
            return msh.is_boundary(fc);
        };

        num_all_faces = msh.faces_size();
        num_dirichlet_faces = std::count_if(msh.faces_begin(), msh.faces_end(), is_dirichlet);
        num_other_faces = num_all_faces - num_dirichlet_faces;

        compress_table.resize( num_all_faces );
        expand_table.resize( num_other_faces );

        size_t compressed_offset = 0;
        for (size_t i = 0; i < num_all_faces; i++)
        {
            auto fc = *std::next(msh.faces_begin(), i);
            if ( !is_dirichlet(fc) )
            {
                compress_table.at(i) = compressed_offset;
                expand_table.at(compressed_offset) = i;
                compressed_offset++;
            }
        }

        auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);

        auto system_size = fbs * num_other_faces;

        LHS = SparseMatrix<T>( system_size, system_size );
        RHS = vector_type::Zero( system_size );
    }

#if 0
    void dump_tables() const
    {
        std::cout << "Compress table: " << std::endl;
        for (size_t i = 0; i < compress_table.size(); i++)
            std::cout << i << " -> " << compress_table.at(i) << std::endl;
    }
#endif

    template<typename Function>
    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const matrix_type& lhs,
             const vector_type& rhs,
             const Function& dirichlet_bf)
    {
        auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
            return msh.is_boundary(fc);
        };

        auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);

        auto fcs = faces(msh, cl);

        std::vector<assembly_index> asm_map;
        asm_map.reserve(fcs.size()*fbs);

        vector_type dirichlet_data = vector_type::Zero(fcs.size()*fbs);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto face_offset = priv::offset(msh, fc);
            auto face_LHS_offset = compress_table.at(face_offset)*fbs;

            bool dirichlet = is_dirichlet(fc);

            for (size_t i = 0; i < fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );

            if (dirichlet)
            {
                auto fb = make_scalar_monomial_basis(msh, fc, di.face_degree());

                matrix_type mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                vector_type rhs = make_rhs(msh, fc, fb, dirichlet_bf, di.face_degree());
                dirichlet_data.block(face_i*fbs, 0, fbs, 1) = mass.ldlt().solve(rhs);
            }
        }

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if ( asm_map[j].assemble() )
                    triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
                else
                    RHS(asm_map[i]) -= lhs(i,j) * dirichlet_data(j);
            }

            RHS(asm_map[i]) += rhs(i);
        }
    } // assemble()

    template<typename Function>
    vector_type
    take_local_data(const Mesh& msh, const typename Mesh::cell_type& cl,
    const vector_type& solution, const Function& dirichlet_bf)
    {
        auto facdeg = di.face_degree();
        auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);
        auto fcs = faces(msh, cl);

        auto num_faces = fcs.size();

        vector_type ret = vector_type::Zero(num_faces*fbs);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];

            auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
                return msh.is_boundary(fc);
            };

            bool dirichlet = is_dirichlet(fc);

            if (dirichlet)
            {
                auto fb = make_scalar_monomial_basis(msh, fc, di.face_degree());

                matrix_type mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                vector_type rhs = make_rhs(msh, fc, fb, dirichlet_bf, di.face_degree());
                ret.block(face_i*fbs, 0, fbs, 1) = mass.ldlt().solve(rhs);
            }
            else
            {
                auto face_offset = priv::offset(msh, fc);
                auto face_SOL_offset = compress_table.at(face_offset)*fbs;
                ret.block(face_i*fbs, 0, fbs, 1) = solution.block(face_SOL_offset, 0, fbs, 1);
            }
        }

        return ret;
    }

    void finalize(void)
    {
        LHS.setFromTriplets( triplets.begin(), triplets.end() );
        triplets.clear();

        dump_sparse_matrix(LHS, "diff.dat");
    }

    size_t num_assembled_faces() const
    {
        return num_other_faces;
    }

};
template<typename Mesh>
auto make_diffusion_assembler(const Mesh& msh, hho_degree_info hdi)
{
    return diffusion_condensed_assembler<Mesh>(msh, hdi);
}

template<typename Mesh>
class diffusion_condensed_assembler2
{
    using T = typename Mesh::coordinate_type;
    std::vector<size_t>                 compress_table;
    std::vector<size_t>                 expand_table;

    hho_degree_info                     di;
    std::vector< Triplet<T> >           triplets;
    size_t      num_all_faces, num_dirichlet_faces, num_other_faces, system_size;

    typedef disk::mechanics::BoundaryConditionsScalar<Mesh> boundary_type;
    boundary_type m_bnd;

    class assembly_index
    {
        size_t  idx;
        bool    assem;

    public:
        assembly_index(size_t i, bool as)
            : idx(i), assem(as)
        {}

        operator size_t() const
        {
            if (!assem)
                throw std::logic_error("Invalid assembly_index");

            return idx;
        }

        bool assemble() const
        {
            return assem;
        }

        friend std::ostream& operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };

public:
  typedef Matrix<T, Dynamic, Dynamic> matrix_type;
  typedef Matrix<T, Dynamic, 1>       vector_type;

  SparseMatrix<T> LHS;
  vector_type     RHS;

  diffusion_condensed_assembler2(const Mesh& msh, hho_degree_info hdi, const boundary_type& bnd) : di(hdi), m_bnd(bnd)
  {
      auto is_dirichlet = [&](const typename Mesh::face& fc) -> bool {
          auto fc_id = msh.lookup(fc);
          return bnd.is_dirichlet_face(fc_id);
      };

      num_all_faces       = msh.faces_size();
      num_dirichlet_faces = std::count_if(msh.faces_begin(), msh.faces_end(), is_dirichlet);
      num_other_faces     = num_all_faces - num_dirichlet_faces;

      compress_table.resize(num_all_faces);
      expand_table.resize(num_other_faces);

      size_t compressed_offset = 0;
      for (size_t i = 0; i < num_all_faces; i++)
      {
          auto fc = *std::next(msh.faces_begin(), i);
          if (!is_dirichlet(fc))
          {
              compress_table.at(i)               = compressed_offset;
              expand_table.at(compressed_offset) = i;
              compressed_offset++;
          }
      }
      auto fbs         = scalar_basis_size(hdi.face_degree(), Mesh::dimension - 1);
      auto system_size = fbs * num_other_faces;

      LHS = SparseMatrix<T>(system_size, system_size);
      RHS = vector_type::Zero(system_size);
    }

#if 0
    void dump_tables() const
    {
        std::cout << "Compress table: " << std::endl;
        for (size_t i = 0; i < compress_table.size(); i++)
            std::cout << i << " -> " << compress_table.at(i) << std::endl;
    }
#endif

    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const matrix_type& lhs,
             const vector_type& rhs)
    {

        auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);
        auto fcs = faces(msh, cl);

        std::vector<assembly_index> asm_map;
        asm_map.reserve(fcs.size()*fbs);

        vector_type dirichlet_data = vector_type::Zero(fcs.size()*fbs);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto face_offset = priv::offset(msh, fc);
            auto face_LHS_offset = compress_table.at(face_offset)*fbs;

            auto face_id = msh.lookup(fc);
            bool dirichlet = m_bnd.is_dirichlet_face(face_id);

            for (size_t i = 0; i < fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );

            if (dirichlet)
            {
                auto fb = make_scalar_monomial_basis(msh, fc, di.face_degree());
                auto dirichlet_fun = m_bnd.dirichlet_boundary_func(face_id);
                matrix_type mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                vector_type rhs = make_rhs(msh, fc, fb, dirichlet_fun, di.face_degree());
                dirichlet_data.block(face_i*fbs, 0, fbs, 1) = mass.ldlt().solve(rhs);
            }
        }

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if ( asm_map[j].assemble() )
                    triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
                else
                    RHS(asm_map[i]) -= lhs(i,j) * dirichlet_data(j);
            }

            RHS(asm_map[i]) += rhs(i);
        }
    } // assemble()

    void
    impose_neumann_boundary_conditions(const Mesh& msh, const boundary_type& bnd)
    {

        if (bnd.nb_faces_neumann() > 0)
        {
            for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
            {
                const auto bfc = *itor;
                const auto face_id = msh.lookup(bfc);

                if (bnd.is_neumann_face(face_id))
                {
                    if (bnd.is_dirichlet_face(face_id))
                    {
                            throw std::invalid_argument("You tried to impose"
                                "both Dirichlet and Neumann conditions on the same face");
                    }
                    else if (bnd.is_robin_face(face_id))
                    {
                            throw std::invalid_argument("You tried to impose"
                                "both Robin and Neumann conditions on the same face");
		    }
                    else
                    {
                        const size_t face_degree   = di.face_degree();
                        const size_t num_face_dofs = scalar_basis_size(face_degree, Mesh::dimension - 1);
                        std::vector<assembly_index> asm_map;
                        asm_map.reserve(num_face_dofs);

                        auto face_offset = face_id;
                        auto face_LHS_offset =compress_table.at( face_offset)*num_face_dofs;

                        for (size_t i = 0; i < num_face_dofs; i++)
                        {
                            asm_map.push_back( assembly_index(face_LHS_offset+i, true) );
                        }

                        auto fb = make_scalar_monomial_basis(msh, bfc, face_degree);
                        vector_type neumann = make_rhs(msh, bfc, fb, bnd.neumann_boundary_func(face_id), face_degree);

                        assert(neumann.size() == num_face_dofs);
                        for (size_t i = 0; i < neumann.size() ; i++)
                        {
                            RHS(asm_map[i]) += neumann[i];
                        }
                    }
                }
            }
        }
        else
            throw std::invalid_argument("There are no Neumann faces");
    }

    void impose_robin_boundary_conditions(const Mesh& msh, const boundary_type& bnd)
    {

        if (bnd.nb_faces_robin() > 0)
        {
            for (auto itor = msh.boundary_faces_begin(); itor!= msh.boundary_faces_end(); itor++)
            {
                const auto bfc =*itor ;
                const auto eid = find_element_id (msh.faces_begin(), msh.faces_end(), bfc);
                if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
                const auto face_id = eid.second;

                if (bnd.is_robin_face(face_id))
                {
                    if (bnd.is_neumann_face(face_id))
                    {
                        switch (bnd.neumann_boundary_type(face_id))
                        {
                            case disk::mechanics::NEUMANN:
                                throw std::invalid_argument("You tried to impose both Neumann and Robin conditions on the same face");
                                break;
                            default:
                                throw std::logic_error("Unknown Neumann Conditions");
                            break;

                        }
                    }
                    else if (bnd.is_dirichlet_face(face_id))
                    {
                        switch (bnd.dirichlet_boundary_type(face_id))
                        {
                            case disk::mechanics::DIRICHLET:
                                throw std::invalid_argument("You tried to impose both Dirichlet and Robin conditions on the same face");
                                break;
                            default:
                                throw std::logic_error("Unknown Dirichlet Conditions");
                            break;
                        }
                    }
                    else
                    {
                        const size_t face_degree = di.face_degree();
                        const size_t num_face_dofs = scalar_basis_size(face_degree, Mesh::dimension -1);
                        std::vector<assembly_index> asm_map;
                        asm_map.reserve(num_face_dofs);

                        auto face_offset = face_id;
                        auto face_LHS_offset = compress_table.at(face_offset)*num_face_dofs;

                        for (size_t i = 0; i < num_face_dofs; i++)
                            asm_map.push_back( assembly_index(face_LHS_offset+i,true) );

                        auto fb = make_scalar_monomial_basis(msh, bfc, face_degree);
                        vector_type robin = make_rhs(msh,bfc,fb,bnd.robin_boundary_func(face_id), face_degree);
                        assert (robin.size() == num_face_dofs);

                        matrix_type mass = make_mass_matrix(msh, bfc, fb, face_degree);

                        for (size_t i = 0; i < num_face_dofs; i++)
                        {
                            RHS(asm_map[i]) += robin[i];

                            for (size_t j = 0; j < num_face_dofs; j++)
                            {
                                if ( asm_map[j].assemble() )
                                    triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], mass(i,j)) );
                            }
                        }
                    }
                }
            }
        }
    }



    vector_type
    take_local_data(const Mesh& msh, const typename Mesh::cell_type& cl,
    const vector_type& solution)
    {
        auto facdeg = di.face_degree();
        auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);
        auto fcs = faces(msh, cl);

        auto num_faces = fcs.size();

        vector_type ret = vector_type::Zero(num_faces*fbs);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];

            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
            const auto face_id=eid.second;

            bool dirichlet = m_bnd.is_dirichlet_face(face_id);

            if (dirichlet)
            {
                auto fb = make_scalar_monomial_basis(msh, fc, di.face_degree());
                auto dirichlet_bf = m_bnd.dirichlet_boundary_func(face_id);
                matrix_type mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                vector_type rhs = make_rhs(msh, fc, fb, dirichlet_bf, di.face_degree());
                ret.block(face_i*fbs, 0, fbs, 1) = mass.ldlt().solve(rhs);
            }
            else
            {
                auto face_offset = priv::offset(msh, fc);
                auto face_SOL_offset = compress_table.at(face_offset)*fbs;
                ret.block(face_i*fbs, 0, fbs, 1) = solution.block(face_SOL_offset, 0, fbs, 1);
            }
        }

        return ret;
    }

    void finalize(void)
    {
        LHS.setFromTriplets( triplets.begin(), triplets.end() );
        triplets.clear();

        dump_sparse_matrix(LHS, "diff.dat");
    }

    size_t num_assembled_faces() const
    {
        return num_other_faces;
    }

};

template<typename Mesh, typename BoundaryType>
auto make_diffusion_assembler2(const Mesh& msh, hho_degree_info hdi,
                                    const BoundaryType& bnd)
{
    return diffusion_condensed_assembler2<Mesh>(msh, hdi, bnd);
}


template<typename Mesh>
class diffusion_condensed_assembler_vector
{
    using T = typename Mesh::coordinate_type;
    std::vector<size_t>                 compress_table;
    std::vector<size_t>                 expand_table;

    hho_degree_info                     di;

    std::vector< Triplet<T> >           triplets;

    size_t   num_all_faces, num_dirichlet_faces, num_other_faces;

    class assembly_index
    {
        size_t  idx;
        bool    assem;

    public:
        assembly_index(size_t i, bool as)
            : idx(i), assem(as)
        {}

        operator size_t() const
        {
            if (!assem)
                throw std::logic_error("Invalid assembly_index");

            return idx;
        }

        bool assemble() const
        {
            return assem;
        }

        friend std::ostream& operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };

public:
  typedef Matrix<T, Dynamic, Dynamic> matrix_type;
  typedef Matrix<T, Dynamic, 1>       vector_type;

  SparseMatrix<T> LHS;
  vector_type     RHS;

  diffusion_condensed_assembler_vector(const Mesh& msh, hho_degree_info hdi) : di(hdi)
  {
      auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool { return msh.is_boundary(fc); };

      num_all_faces       = msh.faces_size();
      num_dirichlet_faces = std::count_if(msh.faces_begin(), msh.faces_end(), is_dirichlet);
      num_other_faces     = num_all_faces - num_dirichlet_faces;

      compress_table.resize(num_all_faces);
      expand_table.resize(num_other_faces);

      size_t compressed_offset = 0;
      for (size_t i = 0; i < num_all_faces; i++)
      {
          auto fc = *std::next(msh.faces_begin(), i);
          if (!is_dirichlet(fc))
          {
              compress_table.at(i)               = compressed_offset;
              expand_table.at(compressed_offset) = i;
              compressed_offset++;
          }
      }

      auto fbs = vector_basis_size(hdi.face_degree(), Mesh::dimension - 1, Mesh::dimension);

      auto system_size = fbs * num_other_faces;

      LHS = SparseMatrix<T>(system_size, system_size);
      RHS = vector_type::Zero(system_size);
    }

#if 0
    void dump_tables() const
    {
        std::cout << "Compress table: " << std::endl;
        for (size_t i = 0; i < compress_table.size(); i++)
            std::cout << i << " -> " << compress_table.at(i) << std::endl;
    }
#endif

    template<typename Function>
    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const matrix_type& lhs,
             const vector_type& rhs,
             const Function& dirichlet_bf)
    {
        auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
            return msh.is_boundary(fc);
        };

        auto fbs = vector_basis_size(di.face_degree(), Mesh::dimension-1, Mesh::dimension);

        auto fcs = faces(msh, cl);

        std::vector<assembly_index> asm_map;
        asm_map.reserve(fcs.size()*fbs);

        vector_type dirichlet_data = vector_type::Zero(fcs.size()*fbs);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto face_offset = priv::offset(msh, fc);
            auto face_LHS_offset = compress_table.at(face_offset)*fbs;

            bool dirichlet = is_dirichlet(fc);

            for (size_t i = 0; i < fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );

            if (dirichlet)
            {
                auto fb = make_vector_monomial_basis(msh, fc, di.face_degree());

                matrix_type mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                vector_type rhs = make_rhs(msh, fc, fb, dirichlet_bf, di.face_degree());
                dirichlet_data.block(face_i*fbs, 0, fbs, 1) = mass.ldlt().solve(rhs);
            }
        }

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if ( asm_map[j].assemble() )
                    triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
                else
                    RHS(asm_map[i]) -= lhs(i,j) * dirichlet_data(j);
            }

            RHS(asm_map[i]) += rhs(i);
        }
    } // assemble()

    template<typename Function>
    vector_type
    take_local_data(const Mesh& msh, const typename Mesh::cell_type& cl,
    const vector_type& solution, const Function& dirichlet_bf)
    {
        auto facdeg = di.face_degree();
        auto fbs = vector_basis_size(di.face_degree(), Mesh::dimension-1,Mesh::dimension);
        auto fcs = faces(msh, cl);

        auto num_faces = fcs.size();

        vector_type ret = vector_type::Zero(num_faces*fbs);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];

            auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
                return msh.is_boundary(fc);
            };

            bool dirichlet = is_dirichlet(fc);

            if (dirichlet)
            {
                auto fb = make_vector_monomial_basis(msh, fc, di.face_degree());

                matrix_type mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                vector_type rhs = make_rhs(msh, fc, fb, dirichlet_bf, di.face_degree());
                ret.block(face_i*fbs, 0, fbs, 1) = mass.ldlt().solve(rhs);
            }
            else
            {
                auto face_offset = priv::offset(msh, fc);
                auto face_SOL_offset = compress_table.at(face_offset)*fbs;
                ret.block(face_i*fbs, 0, fbs, 1) = solution.block(face_SOL_offset, 0, fbs, 1);
            }
        }

        return ret;
    }

    void finalize(void)
    {
        LHS.setFromTriplets( triplets.begin(), triplets.end() );
        triplets.clear();

        dump_sparse_matrix(LHS, "diff.dat");
    }

    size_t num_assembled_faces() const
    {
        return num_other_faces;
    }

};

template<typename Mesh>
auto make_diffusion_assembler_vector(const Mesh& msh, hho_degree_info hdi)
{
    return diffusion_condensed_assembler_vector<Mesh>(msh, hdi);
}


template<typename Mesh>
auto
make_hho_stokes(const Mesh& msh, const typename Mesh::cell_type& cl,
                const hho_degree_info& hdi, const bool& use_sym_grad)
{
    if(use_sym_grad)
        return make_hho_vector_symmetric_laplacian(msh, cl, hdi);
    else
        return make_hho_vector_laplacian(msh, cl, hdi);
}


template<typename Mesh>
auto
make_hlow_stokes(const Mesh& msh, const typename Mesh::cell_type& cl,
                const hho_degree_info& hdi, const bool& use_sym_grad)
{
    if(use_sym_grad)
        return make_hho_sym_gradrec_matrix(msh, cl, hdi);
    else
        return make_hho_gradrec_matrix(msh, cl, hdi);
}


template<typename Mesh>
class stokes_assembler
{
    using T = typename Mesh::coordinate_type;
    typedef disk::mechanics::BoundaryConditions<Mesh>    boundary_type;

    std::vector<size_t>                 compress_table;
    std::vector<size_t>                 expand_table;

    boundary_type                       m_bnd;
    hho_degree_info                     di;
    std::vector< Triplet<T> >           triplets;

    size_t      num_all_faces, num_dirichlet_faces, num_other_faces;
    size_t      cbs_A, cbs_B, fbs_A;
    size_t      system_size;

    class assembly_index
    {
        size_t  idx;
        bool    assem;

    public:
        assembly_index(size_t i, bool as)
            : idx(i), assem(as)
        {}

        operator size_t() const
        {
            if (!assem)
                throw std::logic_error("Invalid assembly_index");

            return idx;
        }

        bool assemble() const
        {
            return assem;
        }

        friend std::ostream& operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };

public:
  typedef Matrix<T, Dynamic, Dynamic> matrix_type;
  typedef Matrix<T, Dynamic, 1>       vector_type;

  SparseMatrix<T> LHS;
  vector_type     RHS;

  stokes_assembler(const Mesh& msh, const hho_degree_info& hdi, const boundary_type& bnd) : di(hdi), m_bnd(bnd)
  {
      auto is_dirichlet = [&](const typename Mesh::face& fc) -> bool {
          auto fc_id = msh.lookup(fc);
          return bnd.is_dirichlet_face(fc_id);
      };

      num_all_faces       = msh.faces_size();
      num_dirichlet_faces = std::count_if(msh.faces_begin(), msh.faces_end(), is_dirichlet);
      num_other_faces     = num_all_faces - num_dirichlet_faces;

      compress_table.resize(num_all_faces);
      expand_table.resize(num_other_faces);

      size_t compressed_offset = 0;
      for (size_t i = 0; i < num_all_faces; i++)
      {
          auto fc = *std::next(msh.faces_begin(), i);
          if (!is_dirichlet(fc))
          {
              compress_table.at(i)               = compressed_offset;
              expand_table.at(compressed_offset) = i;
              compressed_offset++;
          }
      }

      cbs_A = vector_basis_size(di.cell_degree(), Mesh::dimension, Mesh::dimension);
      fbs_A = vector_basis_size(di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
      cbs_B = scalar_basis_size(di.face_degree(), Mesh::dimension);

      system_size = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * msh.cells_size() + 1;

      LHS = SparseMatrix<T>(system_size, system_size);
      RHS = vector_type::Zero(system_size);

      // testing Boundary module of Nicolas
      auto   num_face_dofs = vector_basis_size(di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
      size_t face_dofs     = 0;
      for (size_t face_id = 0; face_id < msh.faces_size(); face_id++)
          face_dofs += num_face_dofs - m_bnd.dirichlet_imposed_dofs(face_id, di.face_degree());

      assert(face_dofs == fbs_A * num_other_faces);
    }

#if 0
    void dump_tables() const
    {
        std::cout << "Compress table: " << std::endl;
        for (size_t i = 0; i < compress_table.size(); i++)
            std::cout << i << " -> " << compress_table.at(i) << std::endl;
    }
#endif
    void
    initialize()
    {
        LHS = SparseMatrix<T>( system_size, system_size );
        RHS = vector_type::Zero( system_size );
        return;
    }

    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const matrix_type& lhs_A,
             const matrix_type& lhs_B,
             const vector_type& rhs)
    {
        auto fcs = faces(msh, cl);

        std::vector<assembly_index> asm_map;
        asm_map.reserve(cbs_A + fcs.size()*fbs_A);

        auto cell_offset        = priv::offset(msh, cl);
        auto cell_LHS_offset    = cell_offset * cbs_A;

        auto B_offset = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * cell_offset;

        for (size_t i = 0; i < cbs_A; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );

        vector_type dirichlet_data = vector_type::Zero(cbs_A + fcs.size()*fbs_A);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto face_offset = priv::offset(msh, fc);
            auto face_LHS_offset = cbs_A * msh.cells_size() + compress_table.at(face_offset)*fbs_A;

            auto fc_id = msh.lookup(fc);
            bool dirichlet = m_bnd.is_dirichlet_face(fc_id);

            for (size_t i = 0; i < fbs_A; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );

            if (dirichlet)
            {
                auto fb = make_vector_monomial_basis(msh, fc, di.face_degree());
                auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
                if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

                const auto face_id                  = eid.second;

                auto dirichlet_fun  = m_bnd.dirichlet_boundary_func(face_id);

                matrix_type mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                vector_type rhs = make_rhs(msh, fc, fb, dirichlet_fun, di.face_degree());
                dirichlet_data.block(cbs_A + face_i*fbs_A, 0, fbs_A, 1) = mass.ldlt().solve(rhs);
            }
        }

        assert( asm_map.size() == lhs_A.rows() && asm_map.size() == lhs_A.cols() );

        for (size_t i = 0; i < lhs_A.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs_A.cols(); j++)
            {
                if ( asm_map[j].assemble() )
                    triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs_A(i,j)) );
                else
                    RHS(asm_map[i]) -= lhs_A(i,j) * dirichlet_data(j);
            }
        }

        for (size_t i = 0; i < lhs_B.rows(); i++)
        {
            for (size_t j = 0; j < lhs_B.cols(); j++)
            {
                auto global_i = B_offset + i;
                auto global_j = asm_map[j];
                if ( asm_map[j].assemble() )
                {
                    triplets.push_back( Triplet<T>(global_i, global_j, lhs_B(i,j)) );
                    triplets.push_back( Triplet<T>(global_j, global_i, lhs_B(i,j)) );
                }
                else
                    RHS(global_i) -= lhs_B(i,j)*dirichlet_data(j);
            }
        }

        auto scalar_cell_basis = make_scalar_monomial_basis(msh, cl, di.face_degree());
        auto qps = integrate(msh, cl, di.face_degree());
        vector_type mult = vector_type::Zero( scalar_cell_basis.size() );
        for (auto& qp : qps)
        {
            auto phi = scalar_cell_basis.eval_functions(qp.point());
            mult += qp.weight() * phi;
        }
        auto mult_offset = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * msh.cells_size();

        for (size_t i = 0; i < mult.rows(); i++)
        {
            triplets.push_back( Triplet<T>(B_offset+i, mult_offset, mult(i)) );
            triplets.push_back( Triplet<T>(mult_offset, B_offset+i, mult(i)) );
        }

        RHS.block(cell_LHS_offset, 0, cbs_A, 1) += rhs.block(0, 0, cbs_A, 1);

    } // assemble()

    void
    assemble_alg(const Mesh& msh, const typename Mesh::cell_type& cl,
             const matrix_type& lhs_A,
             const matrix_type& lhs_B,
             const vector_type& rhs)
    {
        auto fcs = faces(msh, cl);

        std::vector<assembly_index> asm_map;
        asm_map.reserve(cbs_A + fcs.size()*fbs_A);

        auto cell_offset        = priv::offset(msh, cl);
        auto cell_LHS_offset    = cell_offset * cbs_A;

        auto B_offset = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * cell_offset;

        for (size_t i = 0; i < cbs_A; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );

        vector_type dirichlet_data = vector_type::Zero(cbs_A + fcs.size()*fbs_A);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto face_offset = priv::offset(msh, fc);
            auto face_LHS_offset = cbs_A * msh.cells_size() + compress_table.at(face_offset)*fbs_A;

            auto fc_id = msh.lookup(fc);
            bool dirichlet = m_bnd.is_dirichlet_face(fc_id);

            for (size_t i = 0; i < fbs_A; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );

            if (dirichlet)
            {
                auto fb = make_vector_monomial_basis(msh, fc, di.face_degree());
                auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
                if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

                const auto face_id                  = eid.second;

                auto dirichlet_fun  = m_bnd.dirichlet_boundary_func(face_id);

                matrix_type mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                vector_type rhs_bnd = make_rhs(msh, fc, fb, dirichlet_fun, di.face_degree());
                dirichlet_data.block(cbs_A + face_i*fbs_A, 0, fbs_A, 1) = mass.ldlt().solve(rhs_bnd);
            }
        }

        assert( asm_map.size() == lhs_A.rows() && asm_map.size() == lhs_A.cols() );
        assert( asm_map.size() == rhs.rows() );

        for (size_t i = 0; i < lhs_A.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            RHS(asm_map[i]) += rhs(i);

            for (size_t j = 0; j < lhs_A.cols(); j++)
            {
                if ( asm_map[j].assemble() )
                    triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs_A(i,j)) );
                else
                    RHS(asm_map[i]) -= lhs_A(i,j) * dirichlet_data(j);
            }
        }

        for (size_t i = 0; i < lhs_B.rows(); i++)
        {
            for (size_t j = 0; j < lhs_B.cols(); j++)
            {
                auto global_i = B_offset + i;
                auto global_j = asm_map[j];
                if ( asm_map[j].assemble() )
                {
                    triplets.push_back( Triplet<T>(global_i, global_j, lhs_B(i,j)) );
                    triplets.push_back( Triplet<T>(global_j, global_i, lhs_B(i,j)) );
                }
                else
                    RHS(global_i) -= lhs_B(i,j)*dirichlet_data(j);
            }
        }

        auto scalar_cell_basis = make_scalar_monomial_basis(msh, cl, di.face_degree());
        auto qps = integrate(msh, cl, di.face_degree());
        vector_type mult = vector_type::Zero( scalar_cell_basis.size() );
        for (auto& qp : qps)
        {
            auto phi = scalar_cell_basis.eval_functions(qp.point());
            mult += qp.weight() * phi;
        }
        auto mult_offset = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * msh.cells_size();

        for (size_t i = 0; i < mult.rows(); i++)
        {
            triplets.push_back( Triplet<T>(B_offset+i, mult_offset, mult(i)) );
            triplets.push_back( Triplet<T>(mult_offset, B_offset+i, mult(i)) );
        }

        //RHS.block(cell_LHS_offset, 0, cbs_A, 1) += rhs.block(0, 0, cbs_A, 1);

    } // assemble_alg()



    /*
    void
    impose_neumann_boundary_conditions(const Mesh& msh, const boundary_type& bnd)
    {
        const auto face_degree   = di.face_degree();
        const auto num_face_dofs = vector_basis_size(face_degree, Mesh::dimension - 1, Mesh:dimension);

        if (bnd.nb_faces_neumann() > 0)
        {
            for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
            {
                const auto bfc = *itor;

                const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
                if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

                const auto face_id = eid.second;

                if (bnd.is_neumann_face(face_id))
                {
                    const size_t      face_offset = face_compress_map.at(face_id);
                    auto              fb = make_vector_monomial_basis(msh, bfc, face_degree);
                    const vector_type neumann =
                    make_rhs(msh, bfc, fb, bnd.neumann_boundary_func(face_id));

                    assert(neumann.size() == num_face_dofs);

                    if (bnd.is_dirichlet_face(face_id))
                    {
                        switch (bnd.dirichlet_boundary_type(face_id))
                        {
                            case disk::mechanics::DIRICHLET:
                                throw std::invalid_argument("You tried to impose"
                                "both Dirichlet and Neumann conditions on the same face");
                                break;
                            default:
                                throw std::logic_error("Unknown Dirichlet Conditions");
                            break;
                        }
                    }
                    else
                        RHS.segment(face_offset, num_face_dofs) += neumann;
                }
            }
        }
    }
    */
    void finalize(void)
    {
        LHS.setFromTriplets( triplets.begin(), triplets.end() );
        triplets.clear();
    }
    size_t num_assembled_faces() const
    {
        return num_other_faces;
    }

    size_t global_system_size() const
    {
        return system_size;
    }

    vector_type
    take_velocity(  const Mesh& msh, const typename Mesh::cell_type& cl,
                    const vector_type& sol) const
    {
        auto num_faces = howmany_faces(msh, cl);
        auto dim = Mesh::dimension;
        auto cell_ofs = priv::offset(msh, cl);

        vector_type svel(cbs_A + num_faces * fbs_A );
        svel.block(0, 0, cbs_A, 1) = sol.block(cell_ofs * cbs_A, 0, cbs_A, 1);
        auto fcs = faces(msh, cl);
        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto fc = fcs[i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
            const auto face_id                  = eid.second;

            if (m_bnd.is_dirichlet_face( face_id))
            {
                auto fb = make_vector_monomial_basis(msh, fc, di.face_degree());
                matrix_type mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                auto velocity = m_bnd.dirichlet_boundary_func(face_id);
                vector_type rhs = make_rhs(msh, fc, fb, velocity, di.face_degree());
                svel.block(cbs_A + i * fbs_A, 0, fbs_A, 1) = mass.ldlt().solve(rhs);
            }
            else
            {
                auto face_ofs = priv::offset(msh, fc);
                auto global_ofs = cbs_A * msh.cells_size() + compress_table.at(face_ofs)*fbs_A;
                svel.block(cbs_A + i*fbs_A, 0, fbs_A, 1) = sol.block(global_ofs, 0, fbs_A, 1);
            }
        }
        return svel;
    }

    vector_type
    take_pressure( const Mesh& msh, const typename Mesh::cell_type& cl,
                  const vector_type& sol) const
    {
        auto cell_ofs = priv::offset(msh, cl);
        auto pres_ofs = cbs_A * msh.cells_size() + fbs_A * num_other_faces
                                                            + cbs_B * cell_ofs;

        vector_type spres = sol.block(pres_ofs, 0, cbs_B, 1);
        return spres;
    }

    auto
    global_face_offset(const Mesh& msh, const typename Mesh::face_type& fc )
    {
        auto cbs_A = vector_basis_size(di.cell_degree(), Mesh::dimension, Mesh::dimension);
        auto fbs_A = vector_basis_size(di.face_degree(), Mesh::dimension-1, Mesh::dimension);
        auto cbs_B = scalar_basis_size(di.face_degree(), Mesh::dimension);

        auto face_offset = priv::offset(msh, fc);
        return cbs_A * msh.cells_size() + compress_table.at(face_offset)*fbs_A;
    }

    auto
    global_face_offset(const Mesh& msh, const typename Mesh::face_type& fc ) const
    {
        auto cbs_A = vector_basis_size(di.cell_degree(), Mesh::dimension, Mesh::dimension);
        auto fbs_A = vector_basis_size(di.face_degree(), Mesh::dimension-1, Mesh::dimension);
        auto cbs_B = scalar_basis_size(di.face_degree(), Mesh::dimension);

        auto face_offset = priv::offset(msh, fc);
        return cbs_A * msh.cells_size() + compress_table.at(face_offset)*fbs_A;
    }
};

template<typename Mesh>
auto make_stokes_assembler(const Mesh& msh, hho_degree_info hdi,
                            const  disk::mechanics::BoundaryConditions<Mesh>& bnd)
{
    return stokes_assembler<Mesh>(msh, hdi, bnd);
}

template<typename Mesh>
class stokes_assembler_alg
{
    using T = typename Mesh::coordinate_type;
    typedef disk::mechanics::BoundaryConditions<Mesh>    boundary_type;

    std::vector<size_t>                 compress_table;
    std::vector<size_t>                 expand_table;

    boundary_type                       m_bnd;
    hho_degree_info                     di;
    std::vector< Triplet<T> >           triplets;
    Matrix<T, Dynamic, 1>               RHS_DIRICHLET;

    size_t      num_all_faces, num_dirichlet_faces, num_other_faces;
    size_t      cbs_A, cbs_B, fbs_A;
    size_t      system_size;

    class assembly_index
    {
        size_t  idx;
        bool    assem;

    public:
        assembly_index(size_t i, bool as)
            : idx(i), assem(as)
        {}

        operator size_t() const
        {
            if (!assem)
                throw std::logic_error("Invalid assembly_index");

            return idx;
        }

        bool assemble() const
        {
            return assem;
        }

        friend std::ostream& operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };


public:
  typedef Matrix<T, Dynamic, Dynamic> matrix_type;
  typedef Matrix<T, Dynamic, 1>       vector_type;

  SparseMatrix<T>                     LHS;
  vector_type                         RHS;

  stokes_assembler_alg(const Mesh& msh, const hho_degree_info& hdi, const boundary_type& bnd) : di(hdi), m_bnd(bnd)
  {
      auto is_dirichlet = [&](const typename Mesh::face& fc) -> bool {
          auto fc_id = msh.lookup(fc);
          return bnd.is_dirichlet_face(fc_id);
      };

      num_all_faces       = msh.faces_size();
      num_dirichlet_faces = std::count_if(msh.faces_begin(), msh.faces_end(), is_dirichlet);
      num_other_faces     = num_all_faces - num_dirichlet_faces;

      compress_table.resize(num_all_faces);
      expand_table.resize(num_other_faces);

      size_t compressed_offset = 0;
      for (size_t i = 0; i < num_all_faces; i++)
      {
          auto fc = *std::next(msh.faces_begin(), i);
          if (!is_dirichlet(fc))
          {
              compress_table.at(i)               = compressed_offset;
              expand_table.at(compressed_offset) = i;
              compressed_offset++;
          }
      }

      cbs_A = vector_basis_size(di.cell_degree(), Mesh::dimension, Mesh::dimension);
      fbs_A = vector_basis_size(di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
      cbs_B = scalar_basis_size(di.face_degree(), Mesh::dimension);

      system_size = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * msh.cells_size() + 1;

      LHS           = SparseMatrix<T>(system_size, system_size);
      RHS           = vector_type::Zero(system_size);
      RHS_DIRICHLET = vector_type::Zero(system_size);

      // testing Boundary module of Nicolas
      auto   num_face_dofs = vector_basis_size(di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
      size_t face_dofs     = 0;
      for (size_t face_id = 0; face_id < msh.faces_size(); face_id++)
          face_dofs += num_face_dofs - m_bnd.dirichlet_imposed_dofs(face_id, di.face_degree());

      assert(face_dofs == fbs_A * num_other_faces);
    }

    void
    initialize_lhs()
    {
        LHS = SparseMatrix<T>( system_size, system_size );
        RHS_DIRICHLET = vector_type::Zero( system_size );
        return;
    }

    void
    initialize_rhs()
    {
        RHS = vector_type::Zero( system_size );
        return;
    }

    void
    assemble_lhs(const Mesh& msh, const typename Mesh::cell_type& cl,
                const matrix_type& lhs_A,
                const matrix_type& lhs_B)
    {
        auto fcs = faces(msh, cl);
        std::vector<assembly_index>         asm_map;
        asm_map.reserve(cbs_A + fcs.size()*fbs_A);

        auto cell_offset        = priv::offset(msh, cl);
        auto cell_LHS_offset    = cell_offset * cbs_A;

        auto B_offset = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * cell_offset;

        for (size_t i = 0; i < cbs_A; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );

        vector_type dirichlet_data = vector_type::Zero(cbs_A + fcs.size()*fbs_A);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto face_offset = priv::offset(msh, fc);
            auto face_LHS_offset = cbs_A * msh.cells_size() + compress_table.at(face_offset)*fbs_A;

            auto fc_id = msh.lookup(fc);
            bool dirichlet = m_bnd.is_dirichlet_face(fc_id);

            for (size_t i = 0; i < fbs_A; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );

            if (dirichlet)
            {
                auto fb = make_vector_monomial_basis(msh, fc, di.face_degree());
                auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
                if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

                const auto face_id                  = eid.second;

                auto dirichlet_fun  = m_bnd.dirichlet_boundary_func(face_id);

                matrix_type mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                vector_type rhs_bnd = make_rhs(msh, fc, fb, dirichlet_fun, di.face_degree());
                dirichlet_data.block(cbs_A + face_i*fbs_A, 0, fbs_A, 1) = mass.ldlt().solve(rhs_bnd);
            }
        }

        assert( asm_map.size() == lhs_A.rows() && asm_map.size() == lhs_A.cols() );

        for (size_t i = 0; i < lhs_A.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;


            for (size_t j = 0; j < lhs_A.cols(); j++)
            {
                if ( asm_map[j].assemble() )
                    triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs_A(i,j)) );
                else
                    RHS_DIRICHLET(asm_map[i]) -= lhs_A(i,j) * dirichlet_data(j);
            }
        }

        for (size_t i = 0; i < lhs_B.rows(); i++)
        {
            for (size_t j = 0; j < lhs_B.cols(); j++)
            {
                auto global_i = B_offset + i;
                auto global_j = asm_map[j];
                if ( asm_map[j].assemble() )
                {
                    triplets.push_back( Triplet<T>(global_i, global_j, lhs_B(i,j)) );
                    triplets.push_back( Triplet<T>(global_j, global_i, lhs_B(i,j)) );
                }
                else
                    RHS_DIRICHLET(global_i) -= lhs_B(i,j)*dirichlet_data(j);
            }
        }

        auto scalar_cell_basis = make_scalar_monomial_basis(msh, cl, di.face_degree());
        auto qps = integrate(msh, cl, di.face_degree());
        vector_type mult = vector_type::Zero( scalar_cell_basis.size() );
        for (auto& qp : qps)
        {
            auto phi = scalar_cell_basis.eval_functions(qp.point());
            mult += qp.weight() * phi;
        }
        auto mult_offset = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * msh.cells_size();

        for (size_t i = 0; i < mult.rows(); i++)
        {
            triplets.push_back( Triplet<T>(B_offset+i, mult_offset, mult(i)) );
            triplets.push_back( Triplet<T>(mult_offset, B_offset+i, mult(i)) );
        }

        //RHS.block(cell_LHS_offset, 0, cbs_A, 1) += rhs.block(0, 0, cbs_A, 1);

    } // assemble_alg()
    void
    assemble_rhs(const Mesh& msh, const typename Mesh::cell_type& cl,
             const vector_type& rhs)
    {
        auto fcs = faces(msh, cl);
        std::vector<assembly_index>         asm_map;
        asm_map.reserve(cbs_A + fcs.size()*fbs_A);

        auto cell_offset        = priv::offset(msh, cl);
        auto cell_LHS_offset    = cell_offset * cbs_A;

        auto B_offset = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * cell_offset;

        for (size_t i = 0; i < cbs_A; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto face_offset = priv::offset(msh, fc);
            auto face_LHS_offset = cbs_A * msh.cells_size() + compress_table.at(face_offset)*fbs_A;

            auto fc_id = msh.lookup(fc);
            bool dirichlet = m_bnd.is_dirichlet_face(fc_id);

            for (size_t i = 0; i < fbs_A; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
        }

        assert( asm_map.size() == rhs.rows() );

        for (size_t i = 0; i < rhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            RHS(asm_map[i]) += rhs(i);
        }
    } // assemble_alg()

    void finalize_lhs(void)
    {
        LHS.setFromTriplets( triplets.begin(), triplets.end() );
        triplets.clear();
    }

    void finalize_rhs(void)
    {
        RHS += RHS_DIRICHLET;
    }

    size_t global_system_size() const
    {
        return system_size;
    }

    vector_type
    take_velocity(  const Mesh& msh, const typename Mesh::cell_type& cl,
                    const vector_type& sol) const
    {
        auto num_faces = howmany_faces(msh, cl);
        auto dim = Mesh::dimension;
        auto cell_ofs = priv::offset(msh, cl);

        vector_type svel(cbs_A + num_faces * fbs_A );
        svel.block(0, 0, cbs_A, 1) = sol.block(cell_ofs * cbs_A, 0, cbs_A, 1);
        auto fcs = faces(msh, cl);
        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto fc = fcs[i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

            const auto face_id                  = eid.second;

            if (m_bnd.is_dirichlet_face(face_id) )
            {
                auto fb = make_vector_monomial_basis(msh, fc, di.face_degree());
                matrix_type mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                auto velocity = m_bnd.dirichlet_boundary_func(face_id);
                vector_type rhs = make_rhs(msh, fc, fb, velocity, di.face_degree());
                svel.block(cbs_A + i * fbs_A, 0, fbs_A, 1) = mass.ldlt().solve(rhs);
            }
            else
            {
                auto face_ofs = priv::offset(msh, fc);
                auto global_ofs = cbs_A * msh.cells_size() + compress_table.at(face_ofs)*fbs_A;
                svel.block(cbs_A + i*fbs_A, 0, fbs_A, 1) = sol.block(global_ofs, 0, fbs_A, 1);
            }
        }
        return svel;
    }

    vector_type
    take_pressure( const Mesh& msh, const typename Mesh::cell_type& cl,
                  const vector_type& sol) const
    {
        auto cell_ofs = priv::offset(msh, cl);
        auto pres_ofs = cbs_A * msh.cells_size() + fbs_A * num_other_faces
                                                            + cbs_B * cell_ofs;

        vector_type spres = sol.block(pres_ofs, 0, cbs_B, 1);
        return spres;
    }
};

template<typename Mesh>
auto make_stokes_assembler_alg(const Mesh& msh, hho_degree_info hdi,
                            const  disk::mechanics::BoundaryConditions<Mesh>& bnd)
{
    return stokes_assembler_alg<Mesh>(msh, hdi, bnd);
}

template<typename Mesh>
class assembler_mechanics
{
   typedef disk::mechanics::BoundaryConditions<Mesh>    bnd_type;
   typedef Mesh                                         mesh_type;
   typedef typename mesh_type::coordinate_type              scalar_type;
   typedef typename mesh_type::cell                     cell_type;
   typedef typename mesh_type::face                     face_type;

   typedef dynamic_matrix<scalar_type>  matrix_type;
   typedef dynamic_vector<scalar_type>  vector_type;
   typedef sparse_matrix<scalar_type>   sparse_type;
   typedef triplet<scalar_type>         triplet_type;

   const static size_t dimension = mesh_type::dimension;

   std::vector<triplet_type> m_triplets;
   size_t                    m_num_unknowns;
   std::vector<size_t>       face_compress_map;
   hho_degree_info           m_hdi;

 public:

   sparse_type         LHS;
   vector_type         RHS;

   assembler_mechanics(){}

   assembler_mechanics(const mesh_type&       msh,
                       const hho_degree_info& hdi,
                       const bnd_type&        bnd) :
     m_hdi(hdi)
   {
      const auto num_face_dofs = vector_basis_size(m_hdi.face_degree(), dimension-1, dimension);

      face_compress_map.resize(msh.faces_size());

      size_t total_dofs = 0;
      for (size_t face_id = 0; face_id < msh.faces_size(); face_id++) {

         face_compress_map.at(face_id) = total_dofs;
         const auto free_dofs = num_face_dofs - bnd.dirichlet_imposed_dofs(face_id, m_hdi.face_degree());
         total_dofs += free_dofs;
      }
      m_num_unknowns = total_dofs;
      LHS            = sparse_type(m_num_unknowns, m_num_unknowns);
      RHS            = vector_type::Zero(m_num_unknowns);
   }

    // don't forget to reset RHS at each Newton iteration
   void
   setZeroRhs()
   {
       RHS.setZero();
   }

   template<typename LocalContrib>
   void
   assemble(const mesh_type& msh, const cell_type& cl, const bnd_type& bnd, const LocalContrib& lc, int di = 0)
   {
      const size_t      face_degree   = m_hdi.face_degree();
      const auto        num_face_dofs = vector_basis_size(face_degree, dimension - 1, dimension);
      const scalar_type zero          = 0;

      const auto          fcs = faces(msh, cl);
      std::vector<size_t> l2g(fcs.size() * num_face_dofs);
      std::vector<bool> l2l(fcs.size() * num_face_dofs, true);
      vector_type rhs_bc = vector_type::Zero(fcs.size() * num_face_dofs);

      for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
         const auto fc = fcs[face_i];

         auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const auto face_id                  = eid.second;
         const bool fc_is_dirichlet_boundary = bnd.is_dirichlet_face(face_id);
         const auto face_offset              = face_compress_map.at(face_id);
         const auto pos                      = face_i * num_face_dofs;

         if (!fc_is_dirichlet_boundary) {
            for (size_t i = 0; i < num_face_dofs; i++) {
               l2g.at(pos + i) = face_offset + i;
            }
         } else {
            size_t ind_sol = 0;

            vector_type proj_bcf =
              project_function(msh, fc, face_degree, bnd.dirichlet_boundary_func(face_id), di);

            bool ind_ok = false;
            for (size_t face_j = 0; face_j < fcs.size(); face_j++) {
               const auto fcj  = fcs[face_j];
               auto       eidj = find_element_id(msh.faces_begin(), msh.faces_end(), fcj);
               if (!eidj.first) throw std::invalid_argument("This is a bug: face not found");

               const auto face_idj                  = eidj.second;
               const bool fcj_is_dirichlet_boundary = bnd.is_dirichlet_face(face_idj);

               matrix_type mat_Fj =
                 lc.first.block(face_j * num_face_dofs, pos, num_face_dofs, num_face_dofs);

               switch (bnd.dirichlet_boundary_type(face_id)) {
                  case disk::mechanics::DIRICHLET: {
                     if (!ind_ok) {
                        for (size_t i = 0; i < num_face_dofs; i++) {
                           l2g.at(pos + i) = 0xDEADBEEF;
                           l2l.at(pos + i) = false;
                        }
                        ind_ok = true;
                     }
                     break;
                  }
                  case disk::mechanics::CLAMPED: {
                     proj_bcf.setZero();
                     mat_Fj.setZero();
                     if (!ind_ok) {
                        for (size_t i = 0; i < num_face_dofs; i++) {
                           l2g.at(pos + i) = 0xDEADBEEF;
                           l2l.at(pos + i) = false;
                        }
                        ind_ok = true;
                     }
                     break;
                  }
                  case disk::mechanics::DX: {
                     for (size_t i = 0; i < num_face_dofs; i += dimension) {
                        mat_Fj.col(i + 1).setZero();
                        proj_bcf(i + 1) = zero;
                        if (dimension == 3) {
                           mat_Fj.col(i + 2).setZero();
                           proj_bcf(i + 2) = zero;
                        }
                        if (!ind_ok) {
                           l2g.at(pos + i)     = 0xDEADBEEF;
                           l2l.at(pos + i)     = false;
                           l2g.at(pos + i + 1) = face_offset + ind_sol++;
                           if (dimension == 3) {
                              l2g.at(pos + i + 2) = face_offset + ind_sol++;
                           }
                        }
                     }
                     ind_ok = true;
                     break;
                  }
                  case disk::mechanics::DY: {
                     for (size_t i = 0; i < num_face_dofs; i += dimension) {
                        mat_Fj.col(i).setZero();
                        proj_bcf(i) = zero;
                        if (dimension == 3) {
                           mat_Fj.col(i + 2).setZero();
                           proj_bcf(i + 2) = zero;
                        }
                        if (!ind_ok) {
                           l2g.at(pos + i)     = face_offset + ind_sol++;
                           l2g.at(pos + i + 1) = 0xDEADBEEF;
                           l2l.at(pos + i + 1) = false;
                           if (dimension == 3) {
                              l2g.at(pos + i + 2) = face_offset + ind_sol++;
                           }
                        }
                     }
                     ind_ok = true;
                     break;
                  }
                  case disk::mechanics::DZ: {
                     if (dimension != 3) throw std::invalid_argument("You are not in 3D");
                     for (size_t i = 0; i < num_face_dofs; i += dimension) {
                        mat_Fj.col(i).setZero();
                        proj_bcf(i) = zero;
                        mat_Fj.col(i + 1).setZero();
                        proj_bcf(i + 1) = zero;
                        if (!ind_ok) {
                           l2g.at(pos + i)     = face_offset + ind_sol++;
                           l2g.at(pos + i + 1) = face_offset + ind_sol++;
                           l2g.at(pos + i + 2) = 0xDEADBEEF;
                           l2l.at(pos + i + 2) = false;
                        }
                     }
                     ind_ok = true;
                     break;
                  }
                  case disk::mechanics::DXDY: {
                     for (size_t i = 0; i < num_face_dofs; i += dimension) {
                        if (dimension == 3) {
                           mat_Fj.col(i + 2).setZero();
                           proj_bcf(i + 2) = zero;
                        }
                        if (!ind_ok) {
                           l2g.at(pos + i)     = 0xDEADBEEF;
                           l2g.at(pos + i + 1) = 0xDEADBEEF;
                           l2l.at(pos + i)     = false;
                           l2l.at(pos + i + 1) = false;
                           if (dimension == 3) {
                              l2g.at(pos + i + 2) = face_offset + ind_sol++;
                           }
                        }
                     }
                     ind_ok = true;
                     break;
                  }
                  case disk::mechanics::DXDZ: {
                     if (dimension != 3) throw std::invalid_argument("You are not in 3D");
                     for (size_t i = 0; i < num_face_dofs; i += dimension) {
                        mat_Fj.col(i + 1).setZero();
                        proj_bcf(i + 1) = zero;
                        if (!ind_ok) {
                           l2g.at(pos + i)     = 0xDEADBEEF;
                           l2l.at(pos + i)     = false;
                           l2g.at(pos + i + 1) = face_offset + ind_sol++;
                           l2g.at(pos + i + 2) = 0xDEADBEEF;
                           l2l.at(pos + i + 2) = false;
                        }
                     }
                     ind_ok = true;
                     break;
                  }
                  case disk::mechanics::DYDZ: {
                     if (dimension != 3) throw std::invalid_argument("You are not in 3D");
                     for (size_t i = 0; i < num_face_dofs; i += dimension) {
                        mat_Fj.col(i).setZero();
                        proj_bcf(i) = zero;
                        if (!ind_ok) {
                           l2g.at(pos + i)     = face_offset + ind_sol++;
                           l2g.at(pos + i + 1) = 0xDEADBEEF;
                           l2g.at(pos + i + 2) = 0xDEADBEEF;
                           l2l.at(pos + i + 1) = false;
                           l2l.at(pos + i + 2) = false;
                        }
                     }
                     ind_ok = true;
                     break;
                  }
                  default: {
                     throw std::logic_error("Unknown Dirichlet Conditions (assembler)");
                     break;
                  }
               }

               rhs_bc.segment(face_j * num_face_dofs, num_face_dofs) += mat_Fj * proj_bcf;
            }
         }
      }
      assert(lc.first.rows() == lc.first.cols());
      assert(lc.first.rows() == lc.second.size());
      assert(lc.second.size() == l2g.size());
      assert(lc.second.size() == rhs_bc.size());

#ifdef FILL_COLMAJOR
      for (size_t j = 0; j < lc.first.cols(); j++) {
         if (l2l[j]){
            for (size_t i = 0; i < lc.first.rows(); i++) {
                if (l2l[i]){
                    m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));
                }
            }
            RHS(l2g.at(j)) += lc.second(j) - rhs_bc(j);
         }
      }
#else
      for (size_t i = 0; i < lc.first.rows(); i++) {
         if (l2l[i]){
            for (size_t j = 0; j < lc.first.cols(); j++) {
                if (l2l[j]){
                    m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));
                }
            }
            RHS(l2g.at(i)) += lc.second(i) - rhs_bc(i);
         }
      }
#endif
   }

   vector_type
   expand_solution(const mesh_type& msh, const bnd_type& bnd, const vector_type& solution, int di = 0)
   {
      assert(solution.size() == m_num_unknowns);
      const auto face_degree   = m_hdi.face_degree();
      const auto num_face_dofs = vector_basis_size(face_degree, dimension-1, dimension);

      vector_type ret = vector_type::Zero(num_face_dofs * msh.faces_size());

      for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++) {
         const auto bfc = *itor;
         const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const auto face_id         = eid.second;
         const auto face_offset     = face_id * num_face_dofs;
         const auto compress_offset = face_compress_map.at(face_id);

         if (bnd.is_dirichlet_face(face_id)) {
            size_t sol_ind = 0;

            const vector_type proj_bcf =
              project_function(msh, bfc, face_degree, bnd.dirichlet_boundary_func(face_id), di);

            assert(proj_bcf.size() == num_face_dofs);

            switch (bnd.dirichlet_boundary_type(face_id)) {
               case disk::mechanics::DIRICHLET: {
                  ret.segment(face_offset, num_face_dofs) = proj_bcf;
                  break;
               }
               case disk::mechanics::CLAMPED: {
                  ret.segment(face_offset, num_face_dofs).setZero();
                  break;
               }
               case disk::mechanics::DX: {

                  for (size_t i = 0; i < num_face_dofs; i += dimension) {
                     ret(face_offset + i)     = proj_bcf(i);
                     ret(face_offset + i + 1) = solution(compress_offset + sol_ind++);
                     if (dimension == 3) {
                        ret(face_offset + i + 2) = solution(compress_offset + sol_ind++);
                     }
                  }
                  break;
               }
               case disk::mechanics::DY: {
                  for (size_t i = 0; i < num_face_dofs; i += dimension) {
                     ret(face_offset + i)     = solution(compress_offset + sol_ind++);
                     ret(face_offset + i + 1) = proj_bcf(i + 1);
                     if (dimension == 3) {
                        ret(face_offset + i + 2) = solution(compress_offset + sol_ind++);
                     }
                  }
                  break;
               }
               case disk::mechanics::DZ: {
                  if (dimension != 3) throw std::invalid_argument("You are not in 3D");
                  for (size_t i = 0; i < num_face_dofs; i += dimension) {
                     ret(face_offset + i)     = solution(compress_offset + sol_ind++);
                     ret(face_offset + i + 1) = solution(compress_offset + sol_ind++);
                     ret(face_offset + i + 2) = proj_bcf(i + 2);
                  }
                  break;
               }
               case disk::mechanics::DXDY: {
                  for (size_t i = 0; i < num_face_dofs; i += dimension) {
                     ret(face_offset + i)     = proj_bcf(i);
                     ret(face_offset + i + 1) = proj_bcf(i + 1);
                     if (dimension == 3) {
                        ret(face_offset + i + 2) = solution(compress_offset + sol_ind++);
                     }
                  }
                  break;
               }
               case disk::mechanics::DXDZ: {
                  if (dimension != 3) throw std::invalid_argument("You are not in 3D");
                  for (size_t i = 0; i < num_face_dofs; i += dimension) {
                     ret(face_offset + i)     = proj_bcf(i);
                     ret(face_offset + i + 1) = solution(compress_offset + sol_ind++);
                     ret(face_offset + i + 2) = proj_bcf(i + 2);
                  }
                  break;
               }
               case disk::mechanics::DYDZ: {
                  if (dimension != 3) throw std::invalid_argument("You are not in 3D");
                  for (size_t i = 0; i < num_face_dofs; i += dimension) {
                     ret(face_offset + i)     = solution(compress_offset + sol_ind++);
                     ret(face_offset + i + 1) = proj_bcf(i + 1);
                     ret(face_offset + i + 2) = proj_bcf(i + 2);
                  }
                  break;
               }
               default: {
                  throw std::logic_error("Unknown Dirichlet Conditions (assembler)");
                  break;
               }
            }
         } else {
            ret.segment(face_offset, num_face_dofs) =
              solution.segment(compress_offset, num_face_dofs);
         }
      }

      return ret;
   }

   template<typename LocalContrib>
   void
   assemble_nl(const mesh_type&                msh,
               const cell_type&                cl,
               const bnd_type&                 bnd,
               const LocalContrib&             lc,
               const std::vector<vector_type>& sol_F,
               int di = 0)
   {
      assert(sol_F.size() == msh.faces_size());
      const size_t      face_degree   = m_hdi.face_degree();
      const auto        num_face_dofs = vector_basis_size(face_degree, dimension - 1, dimension);
      const scalar_type zero          = 0;

      const auto          fcs = faces(msh, cl);
      std::vector<size_t> l2g(fcs.size() * num_face_dofs);
      std::vector<bool>   l2l(fcs.size() * num_face_dofs, true);
      vector_type         rhs_bc = vector_type::Zero(fcs.size() * num_face_dofs);

      for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
         const auto fc = fcs[face_i];

         auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const auto face_id                  = eid.second;
         const bool fc_is_dirichlet_boundary = bnd.is_dirichlet_face(face_id);
         const auto face_offset              = face_compress_map.at(face_id);
         const auto pos                      = face_i * num_face_dofs;

         if (!fc_is_dirichlet_boundary) {
            for (size_t i = 0; i < num_face_dofs; i++) {
               l2g.at(pos + i) = face_offset + i;
            }
         } else {
            size_t ind_sol = 0;

            const vector_type proj_bcf =
              project_function(msh, fc, face_degree, bnd.dirichlet_boundary_func(face_id), di);
            assert(proj_bcf.size() == sol_F[face_id].size());

            vector_type incr   = proj_bcf - sol_F[face_id];
            bool        ind_ok = false;
            for (size_t face_j = 0; face_j < fcs.size(); face_j++) {
               const auto fcj  = fcs[face_j];
               auto       eidj = find_element_id(msh.faces_begin(), msh.faces_end(), fcj);
               if (!eidj.first) throw std::invalid_argument("This is a bug: face not found");

               const auto face_idj                  = eidj.second;
               const bool fcj_is_dirichlet_boundary = bnd.is_dirichlet_face(face_idj);

               matrix_type mat_Fj =
                 lc.first.block(face_j * num_face_dofs, pos, num_face_dofs, num_face_dofs);

               switch (bnd.dirichlet_boundary_type(face_id)) {
                  case disk::mechanics::DIRICHLET: {
                     if (!ind_ok) {
                        for (size_t i = 0; i < num_face_dofs; i++) {
                           l2g.at(pos + i) = 0xDEADBEEF;
                           l2l.at(pos + i) = false;
                        }
                        ind_ok = true;
                     }
                     break;
                  }
                  case disk::mechanics::CLAMPED: {
                     incr = -sol_F[face_id];
                     if (!ind_ok) {
                        for (size_t i = 0; i < num_face_dofs; i++) {
                           l2g.at(pos + i) = 0xDEADBEEF;
                           l2l.at(pos + i) = false;
                        }
                        ind_ok = true;
                     }
                     break;
                  }
                  case disk::mechanics::DX: {
                     for (size_t i = 0; i < num_face_dofs; i += dimension) {
                        mat_Fj.col(i + 1).setZero();
                        incr(i + 1) = zero;
                        if (dimension == 3) {
                           mat_Fj.col(i + 2).setZero();
                           incr(i + 2) = zero;
                        }
                        if (!ind_ok) {
                           l2g.at(pos + i)     = 0xDEADBEEF;
                           l2l.at(pos + i)     = false;
                           l2g.at(pos + i + 1) = face_offset + ind_sol++;
                           if (dimension == 3) {
                              l2g.at(pos + i + 2) = face_offset + ind_sol++;
                           }
                        }
                     }
                     ind_ok = true;
                     break;
                  }
                  case disk::mechanics::DY: {
                     for (size_t i = 0; i < num_face_dofs; i += dimension) {
                        mat_Fj.col(i).setZero();
                        incr(i) = zero;
                        if (dimension == 3) {
                           mat_Fj.col(i + 2).setZero();
                           incr(i + 2) = zero;
                        }
                        if (!ind_ok) {
                           l2g.at(pos + i)     = face_offset + ind_sol++;
                           l2g.at(pos + i + 1) = 0xDEADBEEF;
                           l2l.at(pos + i + 1) = false;
                           if (dimension == 3) {
                              l2g.at(pos + i + 2) = face_offset + ind_sol++;
                           }
                        }
                     }
                     ind_ok = true;
                     break;
                  }
                  case disk::mechanics::DZ: {
                     if (dimension != 3) throw std::invalid_argument("You are not in 3D");
                     for (size_t i = 0; i < num_face_dofs; i += dimension) {
                        mat_Fj.col(i).setZero();
                        incr(i) = zero;
                        mat_Fj.col(i + 1).setZero();
                        incr(i + 1) = zero;
                        if (!ind_ok) {
                           l2g.at(pos + i)     = face_offset + ind_sol++;
                           l2g.at(pos + i + 1) = face_offset + ind_sol++;
                           l2g.at(pos + i + 2) = 0xDEADBEEF;
                           l2l.at(pos + i + 2) = false;
                        }
                     }
                     ind_ok = true;
                     break;
                  }
                  case disk::mechanics::DXDY: {
                     for (size_t i = 0; i < num_face_dofs; i += dimension) {
                        if (dimension == 3) {
                           mat_Fj.col(i + 2).setZero();
                           incr(i + 2) = zero;
                        }
                        if (!ind_ok) {
                           l2g.at(pos + i)     = 0xDEADBEEF;
                           l2g.at(pos + i + 1) = 0xDEADBEEF;
                           l2l.at(pos + i)     = false;
                           l2l.at(pos + i + 1) = false;
                           if (dimension == 3) {
                              l2g.at(pos + i + 2) = face_offset + ind_sol++;
                           }
                        }
                     }
                     ind_ok = true;
                     break;
                  }
                  case disk::mechanics::DXDZ: {
                     if (dimension != 3) throw std::invalid_argument("You are not in 3D");
                     for (size_t i = 0; i < num_face_dofs; i += dimension) {
                        mat_Fj.col(i + 1).setZero();
                        incr(i + 1) = zero;
                        if (!ind_ok) {
                           l2g.at(pos + i)     = 0xDEADBEEF;
                           l2l.at(pos + i)     = false;
                           l2g.at(pos + i + 1) = face_offset + ind_sol++;
                           l2g.at(pos + i + 2) = 0xDEADBEEF;
                           l2l.at(pos + i + 2) = false;
                        }
                     }
                     ind_ok = true;
                     break;
                  }
                  case disk::mechanics::DYDZ: {
                     if (dimension != 3) throw std::invalid_argument("You are not in 3D");
                     for (size_t i = 0; i < num_face_dofs; i += dimension) {
                        mat_Fj.col(i).setZero();
                        incr(i) = zero;
                        if (!ind_ok) {
                           l2g.at(pos + i)     = face_offset + ind_sol++;
                           l2g.at(pos + i + 1) = 0xDEADBEEF;
                           l2g.at(pos + i + 2) = 0xDEADBEEF;
                           l2l.at(pos + i + 1) = false;
                           l2l.at(pos + i + 2) = false;
                        }
                     }
                     ind_ok = true;
                     break;
                  }
                  default: {
                     throw std::logic_error("Unknown Dirichlet Conditions");
                     break;
                  }
               }

               rhs_bc.segment(face_j * num_face_dofs, num_face_dofs) += mat_Fj * incr;
            }
         }
      }
      assert(lc.first.rows() == lc.first.cols());
      assert(lc.first.rows() == lc.second.size());
      assert(lc.second.size() == l2g.size());
      assert(lc.second.size() == rhs_bc.size());

#ifdef FILL_COLMAJOR
      for (size_t j = 0; j < lc.first.cols(); j++) {
         if (l2l[j]){
            for (size_t i = 0; i < lc.first.rows(); i++) {
                if (l2l[i]){
                    m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));
                }
            }
            RHS(l2g.at(j)) += lc.second(j) - rhs_bc(j);
         }
      }
#else
      for (size_t i = 0; i < lc.first.rows(); i++) {
         if (l2l[i]){
            for (size_t j = 0; j < lc.first.cols(); j++) {
                if (l2l[j]){
                    m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));
                }
            }
            RHS(l2g.at(i)) += lc.second(i) - rhs_bc(i);
         }
      }
#endif
   }

   vector_type
   expand_solution_nl(const mesh_type&                msh,
                      const bnd_type&                 bnd,
                      const vector_type&              solution,
                      const std::vector<vector_type>& sol_F,
                      int di = 0)
   {
      assert(solution.size() == m_num_unknowns);
      assert(sol_F.size() == msh.faces_size());
      const auto face_degree   = m_hdi.face_degree();
      const auto num_face_dofs = vector_basis_size(face_degree, dimension - 1, dimension);

      vector_type ret = vector_type::Zero(num_face_dofs * msh.faces_size());

      for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++) {
         const auto bfc = *itor;
         const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
         if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

         const auto face_id         = eid.second;
         const auto face_offset     = face_id * num_face_dofs;
         const auto compress_offset = face_compress_map.at(face_id);

         if (bnd.is_dirichlet_face(face_id)) {
            size_t sol_ind = 0;

            const vector_type proj_bcf =
              project_function(msh, bfc, face_degree, bnd.dirichlet_boundary_func(face_id), di);
            vector_type incr = proj_bcf - sol_F[face_id];
            assert(proj_bcf.size() == num_face_dofs);

            switch (bnd.dirichlet_boundary_type(face_id)) {
               case disk::mechanics::DIRICHLET: {
                  ret.segment(face_offset, num_face_dofs) = incr;
                  break;
               }
               case disk::mechanics::CLAMPED: {
                  incr                                    = -sol_F[face_id];
                  ret.segment(face_offset, num_face_dofs) = incr;
                  break;
               }
               case disk::mechanics::DX: {

                  for (size_t i = 0; i < num_face_dofs; i += dimension) {
                     ret(face_offset + i)     = incr(i);
                     ret(face_offset + i + 1) = solution(compress_offset + sol_ind++);
                     if (dimension == 3) {
                        ret(face_offset + i + 2) = solution(compress_offset + sol_ind++);
                     }
                  }
                  break;
               }
               case disk::mechanics::DY: {
                  for (size_t i = 0; i < num_face_dofs; i += dimension) {
                     ret(face_offset + i)     = solution(compress_offset + sol_ind++);
                     ret(face_offset + i + 1) = incr(i + 1);
                     if (dimension == 3) {
                        ret(face_offset + i + 2) = solution(compress_offset + sol_ind++);
                     }
                  }
                  break;
               }
               case disk::mechanics::DZ: {
                  if (dimension != 3) throw std::invalid_argument("You are not in 3D");
                  for (size_t i = 0; i < num_face_dofs; i += dimension) {
                     ret(face_offset + i)     = solution(compress_offset + sol_ind++);
                     ret(face_offset + i + 1) = solution(compress_offset + sol_ind++);
                     ret(face_offset + i + 2) = incr(i + 2);
                  }
                  break;
               }
               case disk::mechanics::DXDY: {
                  for (size_t i = 0; i < num_face_dofs; i += dimension) {
                     ret(face_offset + i)     = incr(i);
                     ret(face_offset + i + 1) = incr(i + 1);
                     if (dimension == 3) {
                        ret(face_offset + i + 2) = solution(compress_offset + sol_ind++);
                     }
                  }
                  break;
               }
               case disk::mechanics::DXDZ: {
                  if (dimension != 3) throw std::invalid_argument("You are not in 3D");
                  for (size_t i = 0; i < num_face_dofs; i += dimension) {
                     ret(face_offset + i)     = incr(i);
                     ret(face_offset + i + 1) = solution(compress_offset + sol_ind++);
                     ret(face_offset + i + 2) = incr(i + 2);
                  }
                  break;
               }
               case disk::mechanics::DYDZ: {
                  if (dimension != 3) throw std::invalid_argument("You are not in 3D");
                  for (size_t i = 0; i < num_face_dofs; i += dimension) {
                     ret(face_offset + i)     = solution(compress_offset + sol_ind++);
                     ret(face_offset + i + 1) = incr(i + 1);
                     ret(face_offset + i + 2) = incr(i + 2);
                  }
                  break;
               }
               default: {
                  throw std::logic_error("Unknown Dirichlet Conditions");
                  break;
               }
            }
         } else {
            ret.segment(face_offset, num_face_dofs) =
              solution.segment(compress_offset, num_face_dofs);
         }
      }

      return ret;
   }

   void
   impose_neumann_boundary_conditions(const mesh_type& msh, const bnd_type& bnd)
   {
      const auto face_degree   = m_hdi.face_degree();
      const auto num_face_dofs = vector_basis_size(face_degree, dimension - 1, dimension);

      if (bnd.nb_faces_neumann() > 0) {
         for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++) {
            const auto bfc = *itor;

            const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");

            const auto face_id = eid.second;

            if (bnd.is_neumann_face(face_id)) {
               const size_t      face_offset = face_compress_map.at(face_id);
               auto              fb = make_vector_monomial_basis(msh, bfc, face_degree);
               const vector_type neumann =
                 make_rhs(msh, bfc, fb, bnd.neumann_boundary_func(face_id));

               assert(neumann.size() == num_face_dofs);

               if (bnd.is_dirichlet_face(face_id)) {
                  switch (bnd.dirichlet_boundary_type(face_id)) {
                     case disk::mechanics::DIRICHLET: {
                        throw std::invalid_argument("You tried to impose both Dirichlet and "
                                                    "Neumann conditions on the same face");
                        break;
                     }
                     case disk::mechanics::CLAMPED: {
                        throw std::invalid_argument("You tried to impose both Dirichlet and "
                                                    "Neumann conditions on the same face");
                        break;
                     }
                     case disk::mechanics::DX: {
                        for (size_t i = 0; i < num_face_dofs; i += dimension) {
                           RHS(face_offset + i + 1) += neumann(i + 1);
                           if (dimension == 3) {
                              RHS(face_offset + i + 2) += neumann(i + 2);
                           }
                        }
                        break;
                     }
                     case disk::mechanics::DY: {
                        for (size_t i = 0; i < num_face_dofs; i += dimension) {
                           RHS(face_offset + i) = neumann(i);
                           if (dimension == 3) {
                              RHS(face_offset + i + 2) += neumann(i + 2);
                           }
                        }

                        break;
                     }
                     case disk::mechanics::DZ: {
                        if (dimension != 3) throw std::invalid_argument("You are not in 3D");
                        for (size_t i = 0; i < num_face_dofs; i += dimension) {
                           RHS(face_offset + i) += neumann(i);
                           RHS(face_offset + i + 1) += neumann(i + 1);
                        }
                        break;
                     }
                     case disk::mechanics::DXDY: {
                        for (size_t i = 0; i < num_face_dofs; i += dimension) {
                           if (dimension == 3) {
                              RHS(face_offset + i + 2) += neumann(i + 2);
                           }
                        }
                        break;
                     }
                     case disk::mechanics::DXDZ: {
                        if (dimension != 3) throw std::invalid_argument("You are not in 3D");
                        for (size_t i = 0; i < num_face_dofs; i += dimension) {
                           RHS(face_offset + i + 1) += neumann(i + 1);
                        }
                        break;
                     }
                     case disk::mechanics::DYDZ: {
                        if (dimension != 3) throw std::invalid_argument("You are not in 3D");
                        for (size_t i = 0; i < num_face_dofs; i += dimension) {
                           RHS(face_offset + i) += neumann(i);
                        }
                        break;
                     }
                     default: {
                        throw std::logic_error("Unknown Dirichlet Conditions");
                        break;
                     }
                  }
               } else {
                  RHS.segment(face_offset, num_face_dofs) += neumann;
               }
            }
         }
      }
   }

   void
   finalize()
   {
      LHS.setFromTriplets(m_triplets.begin(), m_triplets.end());
      m_triplets.clear();
   }
};

template<typename Mesh>
auto
make_mechanics_assembler( const Mesh&           msh,
                          const hho_degree_info hdi,
                          const disk::mechanics::BoundaryConditions<Mesh>& bnd)
{
   return assembler_mechanics<Mesh>(msh, hdi, bnd);
}

template<typename Mesh>
class static_condensation_vector
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell            cell_type;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    const static size_t dimension = mesh_type::dimension;

  public:
    matrix_type AL;
    vector_type bL;
    static_condensation_vector() {}

    std::pair<matrix_type, vector_type>
    compute(const mesh_type&      msh,
            const cell_type&      cl,
            const matrix_type&    local_mat,
            const vector_type&    cell_rhs,
            const hho_degree_info hdi)
    {
        const size_t num_cell_dofs = vector_basis_size(hdi.cell_degree(), dimension, dimension);
        const size_t num_face_dofs = vector_basis_size(hdi.face_degree(), dimension - 1, dimension);

        const auto   fcs       = faces(msh, cl);
        const size_t num_faces = fcs.size();

        const size_t cell_size = num_cell_dofs;
        const size_t face_size = num_faces * num_face_dofs;

        assert(cell_size == cell_rhs.rows() && "wrong rhs dimension");
        assert((cell_size + face_size) == local_mat.rows() && "wrong lhs rows dimension");
        assert((cell_size + face_size) == local_mat.cols() && "wrong lhs cols dimension");

        const matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
        const matrix_type K_TF = local_mat.topRightCorner(cell_size, face_size);
        const matrix_type K_FT = local_mat.bottomLeftCorner(face_size, cell_size);
        const matrix_type K_FF = local_mat.bottomRightCorner(face_size, face_size);

        assert(K_TT.cols() == cell_size && "wrong K_TT dimension");
        assert(K_TT.cols() + K_TF.cols() == local_mat.cols());
        assert(K_TT.rows() + K_FT.rows() == local_mat.rows());
        assert(K_TF.rows() + K_FF.rows() == local_mat.rows());
        assert(K_FT.cols() + K_FF.cols() == local_mat.cols());

        const auto K_TT_ldlt = K_TT.ldlt();
        AL                   = K_TT_ldlt.solve(K_TF);
        if (K_TT_ldlt.info() != Eigen::Success)
        {
            throw std::invalid_argument("static_condenstion: the matrix is not symmetrix positive definite");
        }
        bL = K_TT_ldlt.solve(cell_rhs);
        if (K_TT_ldlt.info() != Eigen::Success)
        {
            throw std::invalid_argument("static_condenstion: the matrix is not symmetrix positive definite");
        }

        const matrix_type AC = K_FF - K_FT * AL;
        const vector_type bC = /* no projection on faces, eqn. 26*/ -K_FT * bL;

        return std::make_pair(AC, bC);
    }

    std::pair<matrix_type, vector_type>
    compute_rhsfull(const mesh_type&      msh,
                    const cell_type&      cl,
                    const matrix_type&    local_mat,
                    const vector_type&    rhs,
                    const hho_degree_info hdi)
    {
        const size_t num_cell_dofs = vector_basis_size(hdi.cell_degree(), dimension, dimension);
        const size_t num_face_dofs = vector_basis_size(hdi.face_degree(), dimension - 1, dimension);

        const auto   fcs       = faces(msh, cl);
        const size_t num_faces = fcs.size();

        const size_t cell_size        = num_cell_dofs;
        const size_t total_faces_size = num_faces * num_face_dofs;

        const vector_type cell_rhs  = rhs.head(cell_size);
        const vector_type faces_rhs = rhs.tail(total_faces_size);

        assert(cell_size + total_faces_size == rhs.rows() && "wrong  rhs dimension");
        assert(cell_size == cell_rhs.rows() && "wrong cell rhs dimension");
        assert(total_faces_size == faces_rhs.rows() && "wrong face rhs dimension");
        assert((cell_size + total_faces_size) == local_mat.rows() && "wrong lhs rows dimension");
        assert((cell_size + total_faces_size) == local_mat.cols() && "wrong lhs cols dimension");

        const matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
        const matrix_type K_TF = local_mat.topRightCorner(cell_size, total_faces_size);
        const matrix_type K_FT = local_mat.bottomLeftCorner(total_faces_size, cell_size);
        const matrix_type K_FF = local_mat.bottomRightCorner(total_faces_size, total_faces_size);

        assert(K_TT.cols() == cell_size && "wrong K_TT dimension");
        assert(K_TT.cols() + K_TF.cols() == local_mat.cols());
        assert(K_TT.rows() + K_FT.rows() == local_mat.rows());
        assert(K_TF.rows() + K_FF.rows() == local_mat.rows());
        assert(K_FT.cols() + K_FF.cols() == local_mat.cols());

        const auto K_TT_ldlt = K_TT.ldlt();
        AL                   = K_TT_ldlt.solve(K_TF);

        if (K_TT_ldlt.info() != Eigen::Success)
        {
            const auto K_TT_lu = K_TT.fullPivLu();

            if(K_TT_lu.isInvertible())
            {
                AL = K_TT_lu.solve(K_TF);
                bL = K_TT_lu.solve(cell_rhs);
            }
            else
            {
                throw std::invalid_argument("static_condenstion: the matrix is not inversible");
            }
        }
        else
        {
            bL = K_TT_ldlt.solve(cell_rhs);
            if (K_TT_ldlt.info() != Eigen::Success)
            {
                throw std::invalid_argument("static_condenstion: the matrix is not symmetrix positive definite");
            }
        }

        const matrix_type AC = K_FF - K_FT * AL;
        const vector_type bC = faces_rhs - K_FT * bL;

        return std::make_pair(AC, bC);
    }
};

// define some optimization
namespace priv {

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
compute_lhs_vector(const Mesh&                msh,
                   const typename Mesh::cell& cl,
                   const hho_degree_info& hdi,
                   const dynamic_matrix<typename Mesh::coordinate_type>& lhs_scalar)
{
    typedef typename Mesh::coordinate_type scalar_type;

    const int dimension     = Mesh::dimension;
    const auto num_cell_dofs = vector_basis_size(hdi.cell_degree(), dimension, dimension);
    const auto num_face_dofs = vector_basis_size(hdi.face_degree(), dimension - 1, dimension);

    const auto   fcs     = faces(msh, cl);
    const auto num_faces = fcs.size();

    const auto total_dofs = num_cell_dofs + num_faces * num_face_dofs;

    dynamic_matrix<scalar_type> lhs = dynamic_matrix<scalar_type>::Zero(total_dofs, total_dofs);

    const auto scal_cell_dofs = scalar_basis_size(hdi.cell_degree(), dimension);
    const auto scal_face_dofs = scalar_basis_size(hdi.face_degree(), dimension - 1);
    const auto scal_total_dofs = scal_cell_dofs + num_faces * scal_face_dofs;

    assert(lhs_scalar.rows() == scal_total_dofs);
    assert(lhs_scalar.cols() == scal_total_dofs);

    int row, col;
#ifdef FILL_COLMAJOR
    for (int j = 0; j < scal_total_dofs; j++)
    {
        col = j * dimension;
        for (int i = 0; i < scal_total_dofs; i++)
        {
            row = i * dimension;
            for (int k = 0; k < dimension; k++)
            {
                lhs(row + k, col + k) = lhs_scalar(i, j);
            }
        }
    }
#else
    for (int i = 0; i < scal_total_dofs; i++)
    {
        row = i * dimension;
        for (int j = 0; j < scal_total_dofs; j++)
        {
            col = j * dimension;
            for (int k = 0; k < dimension; k++)
            {
                lhs(row + k, col + k) = lhs_scalar(i, j);
            }
        }
    }
#endif

    return lhs;
}

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
compute_grad_vector(const Mesh&                                       msh,
                    const typename Mesh::cell&                        cl,
                    const hho_degree_info&                            hdi,
                    const dynamic_matrix<typename Mesh::coordinate_type>& grad_scalar)
{
    typedef typename Mesh::coordinate_type scalar_type;

    const int  dimension     = Mesh::dimension;
    const auto rbs           = vector_basis_size(hdi.reconstruction_degree(), dimension, dimension);
    const auto num_cell_dofs = vector_basis_size(hdi.cell_degree(), dimension, dimension);
    const auto num_face_dofs = vector_basis_size(hdi.face_degree(), dimension - 1, dimension);

    const auto fcs       = faces(msh, cl);
    const auto num_faces = fcs.size();

    const auto total_dofs = num_cell_dofs + num_faces * num_face_dofs;

    dynamic_matrix<scalar_type> grad = dynamic_matrix<scalar_type>::Zero(rbs - dimension, total_dofs);

    const auto scal_rbs        = scalar_basis_size(hdi.reconstruction_degree(), dimension);
    const auto scal_cell_dofs  = scalar_basis_size(hdi.cell_degree(), dimension);
    const auto scal_face_dofs  = scalar_basis_size(hdi.face_degree(), dimension - 1);
    const auto scal_total_dofs = scal_cell_dofs + num_faces * scal_face_dofs;

    assert(grad_scalar.rows() == scal_rbs - 1);
    assert(grad_scalar.cols() == scal_total_dofs);

    int row, col;
#ifdef FILL_COLMAJOR
    for (int j = 0; j < scal_total_dofs; j++)
    {
        col = j * dimension;
        for (int i = 0; i < scal_total_dofs; i++)
        {
            row = i * dimension;
            for (int k = 0; k < dimension; k++)
            {
                grad(row + k, col + k) = grad_scalar(i, j);
            }
        }
    }
#else
    for (int i = 0; i < scal_rbs-1; i++)
    {
        row = i * dimension;
        for (int j = 0; j < scal_total_dofs; j++)
        {
            col = j * dimension;
            for (int k = 0; k < dimension; k++)
            {
                grad(row + k, col + k) = grad_scalar(i, j);
            }
        }
    }
#endif

    return grad;
}

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
compute_grad_matrix(const Mesh&                                           msh,
                    const typename Mesh::cell&                            cl,
                    const hho_degree_info&                                hdi,
                    const dynamic_matrix<typename Mesh::coordinate_type>& grad_vector)
{
    typedef typename Mesh::coordinate_type scalar_type;

    const int  dimension     = Mesh::dimension;
    const auto gbs           = matrix_basis_size(hdi.grad_degree(), dimension, dimension);
    const auto num_cell_dofs = vector_basis_size(hdi.cell_degree(), dimension, dimension);
    const auto num_face_dofs = vector_basis_size(hdi.face_degree(), dimension - 1, dimension);

    const auto fcs       = faces(msh, cl);
    const auto num_faces = fcs.size();

    const auto total_dofs = num_cell_dofs + num_faces * num_face_dofs;

    dynamic_matrix<scalar_type> grad = dynamic_matrix<scalar_type>::Zero(gbs, total_dofs);

    const auto vec_gbs         = vector_basis_size(hdi.grad_degree(), dimension, dimension);
    const auto scal_cell_dofs  = scalar_basis_size(hdi.cell_degree(), dimension);
    const auto scal_face_dofs  = scalar_basis_size(hdi.face_degree(), dimension - 1);
    const auto scal_total_dofs = scal_cell_dofs + num_faces * scal_face_dofs;

    assert(grad_vector.rows() == vec_gbs);
    assert(grad_vector.cols() == scal_total_dofs);

    int row, col;
#ifdef FILL_COLMAJOR
    for (int j = 0; j < scal_total_dofs; j++)
    {
        col = j * dimension;
        for (int i = 0; i < vec_gbs; i++)
        {
            row = i * dimension;
            for (int k = 0; k < dimension; k++)
            {
                grad(row + k, col + k) = grad_vector(i, j);
            }
        }
    }
#else
    for (int i = 0; i < vec_gbs; i++)
    {
        row = i * dimension;
        for (int j = 0; j < scal_total_dofs; j++)
        {
            col = j * dimension;
            for (int k = 0; k < dimension; k++)
            {
                grad(row + k, col + k) = grad_vector(i, j);
            }
        }
    }
#endif

    return grad;
}
} // end priv

template<typename Mesh>
std::pair<Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
          Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>>
make_hho_vector_laplacian(const Mesh&                     msh,
                          const typename Mesh::cell_type& cl,
                          const hho_degree_info&          hdi,
                          const std::pair < Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
                          Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic> >& hho_scalar_laplacian)
{
    const auto oper = priv::compute_grad_vector(msh, cl, hdi, hho_scalar_laplacian.first);
    const auto data = priv::compute_lhs_vector(msh, cl, hdi, hho_scalar_laplacian.second);

    return std::make_pair(oper, data);
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hdg_vector_stabilization(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    const auto hdg_scalar_stab = make_hdg_scalar_stabilization(msh, cl, di);

    return priv::compute_lhs_vector(msh, cl, di, hdg_scalar_stab);
}

// doesn't work for symmetric gradient
template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hho_vector_stabilization_optim(const Mesh& msh, const typename Mesh::cell_type& cl,
                              const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>& reconstruction_scalar,
                              const hho_degree_info& hdi)
{
    const auto hho_scalar_stab = make_hho_scalar_stabilization(msh, cl, reconstruction_scalar, hdi);

    return priv::compute_lhs_vector(msh, cl, hdi, hho_scalar_stab);
}

#if 1
template<typename Mesh>
std::pair<Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
          Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>>
make_hho_gradrec_matrix(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& hdi)

{
    const auto hho_gradrec_vector = make_hho_gradrec_vector(msh, cl, hdi);
    const auto lhs                = priv::compute_lhs_vector(msh, cl, hdi, hho_gradrec_vector.second);
    const auto oper               = priv::compute_grad_matrix(msh, cl, hdi, hho_gradrec_vector.first);

    return std::make_pair(oper, lhs);
}

#else
template<typename Mesh>
std::pair<Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
          Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>>
make_hho_gradrec_matrix(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1> vector_type;

    const size_t N = Mesh::dimension;

    const auto graddeg = di.grad_degree();
    const auto celdeg = di.cell_degree();
    const auto facdeg = di.face_degree();

    const auto gb = make_matrix_monomial_basis(msh, cl, graddeg);
    const auto cb = make_vector_monomial_basis(msh, cl, celdeg);

    const auto gbs = matrix_basis_size(graddeg, Mesh::dimension, Mesh::dimension);
    const auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
    const auto fbs = vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension);

    const auto num_faces = howmany_faces(msh, cl);

    matrix_type gr_lhs = matrix_type::Zero(gbs, gbs);
    matrix_type gr_rhs = matrix_type::Zero(gbs, cbs + num_faces * fbs);
    const size_t dim2 = N * N;

    // this is very costly to build it
    const auto qps = integrate(msh, cl, 2 * graddeg);
    for (auto& qp : qps)
    {
        const auto gphi = gb.eval_functions(qp.point());

        for (size_t j = 0; j < gbs; j += dim2)
        {
            const auto qp_gphi_j = priv::inner_product(qp.weight(), gphi[j]);
            for (size_t i = j; i < gbs; i += dim2)
            {
                gr_lhs(i, j) += priv::inner_product(gphi[i], qp_gphi_j);
            }
        }
    }

    // upper part
    for (size_t j = 0; j < gbs; j++)
    {
        for (size_t i = 0; i < j; i++)
        {
            gr_lhs(i, j) = gr_lhs(j, i);
        }
    }

    // copy of each cols
    for (size_t j = 0; j < gbs; j += dim2)
    {
        for (size_t col = 1; col < dim2; col++)
        {
            gr_lhs.block(col, j + col, gbs - dim2 + 1, 1) = gr_lhs.block(0, j, gbs - dim2 + 1, 1);
        }
    }

    // compute rhs
    if (celdeg > 0)
    {
        const auto qpc = integrate(msh, cl, graddeg + celdeg - 1);
        for (auto& qp : qpc)
        {
            const auto gphi = gb.eval_functions(qp.point());
            const auto dphi = cb.eval_gradients(qp.point());

            for (size_t j = 0; j < cbs; j++)
            {
                const auto qp_dphi_j = priv::inner_product(qp.weight(), dphi[j]);
                for (size_t i = 0; i < gbs; i++)
                {
                    gr_rhs(i, j) += priv::inner_product(gphi[i], qp_dphi_j);
                }
            }
        } // end qp
    }

    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n = normal(msh, cl, fc);
        const auto fb = make_vector_monomial_basis(msh, fc, facdeg);

        const auto qps_f = integrate(msh, fc, graddeg + std::max(celdeg, facdeg));
        for (auto& qp : qps_f)
        {
            const auto cphi = cb.eval_functions(qp.point());
            const auto gphi = gb.eval_functions(qp.point());
            const auto fphi = fb.eval_functions(qp.point());

            const auto gphi_n = priv::outer_product(gphi, n);
            gr_rhs.block(0, cbs + i * fbs, gbs, fbs) += qp.weight() * priv::outer_product(fphi, gphi_n);
            gr_rhs.block(0, 0, gbs, cbs) -= qp.weight() * priv::outer_product(cphi, gphi_n);
        }
    }

    matrix_type oper = gr_lhs.ldlt().solve(gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

#endif

} // end revolution
