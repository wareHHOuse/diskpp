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
#include "revolution/bases"
#include "revolution/quadratures"
#include "mechanics/BoundaryConditions.hpp"

using namespace Eigen;

namespace revolution
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
    const size_t DIM = Mesh::dimension;

    auto recdeg = di.reconstruction_degree();
    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    auto cb = make_scalar_monomial_basis(msh, cl, recdeg);

    auto rbs = scalar_basis_size(recdeg, Mesh::dimension);
    auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);

    auto num_faces = howmany_faces(msh, cl);

    Matrix<T, Dynamic, Dynamic> stiff = Matrix<T, Dynamic, Dynamic>::Zero(rbs, rbs);
    Matrix<T, Dynamic, Dynamic> gr_lhs = Matrix<T, Dynamic, Dynamic>::Zero(rbs-1, rbs-1);
    Matrix<T, Dynamic, Dynamic> gr_rhs = Matrix<T, Dynamic, Dynamic>::Zero(rbs-1, cbs + num_faces*fbs);

    auto qps = integrate(msh, cl, 2*recdeg);
    for (auto& qp : qps)
    {
        auto dphi = cb.eval_gradients(qp.point());
        stiff += qp.weight() * dphi * dphi.transpose();
    }

    gr_lhs = stiff.block(1, 1, rbs-1, rbs-1);
    gr_rhs.block(0, 0, rbs-1, cbs) = stiff.block(1, 0, rbs-1, cbs);

    auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        auto fc = fcs[i];
        auto n = normal(msh, cl, fc);
        auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

        auto qps_f = integrate(msh, fc, 2*facdeg);
        for (auto& qp : qps_f)
        {
            Matrix<T, Dynamic, 1> c_phi_tmp = cb.eval_functions(qp.point());
            Matrix<T, Dynamic, 1> c_phi = c_phi_tmp.head(cbs);
            Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
            Matrix<T, Dynamic, DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
            Matrix<T, Dynamic, 1> f_phi = fb.eval_functions(qp.point());
            gr_rhs.block(0, cbs+i*fbs, rbs-1, fbs) += qp.weight() * (c_dphi * n) * f_phi.transpose();
            gr_rhs.block(0, 0, rbs-1, cbs) -= qp.weight() * (c_dphi * n) * c_phi.transpose();
        }
    }

    Matrix<T, Dynamic, Dynamic> oper = gr_lhs.llt().solve(gr_rhs);
    Matrix<T, Dynamic, Dynamic> data = gr_rhs.transpose() * oper;

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

    auto recdeg = di.reconstruction_degree();
    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    auto cb = make_vector_monomial_basis(msh, cl, recdeg);

    auto rbs = vector_basis_size(recdeg, Mesh::dimension, Mesh::dimension);
    auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
    auto fbs = vector_basis_size(facdeg, Mesh::dimension-1, Mesh::dimension);

    auto num_faces = howmany_faces(msh, cl);

    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, N, N>             gradient_type;
    typedef Matrix<T, Dynamic, N>       function_type;

    matrix_type stiff  = matrix_type::Zero(rbs, rbs);
    matrix_type gr_lhs = matrix_type::Zero(rbs-N, rbs-N);
    matrix_type gr_rhs = matrix_type::Zero(rbs-N, cbs + num_faces*fbs);

    auto qps = integrate(msh, cl, 2*recdeg);
    for (auto& qp : qps)
    {
        auto dphi = cb.eval_gradients(qp.point());
        stiff += qp.weight() * priv::outer_product(dphi, dphi);
    }

    gr_lhs = stiff.block(N, N, rbs-N, rbs-N);
    gr_rhs.block(0, 0, rbs-N, cbs) = stiff.block(N, 0, rbs-N, cbs);

    auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        auto fc = fcs[i];
        auto n  = normal(msh, cl, fc);
        auto fb = make_vector_monomial_basis(msh, fc, facdeg);

        auto qps_f = integrate(msh, fc, 2*facdeg);
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

    Matrix<T, Dynamic, Dynamic> oper = gr_lhs.llt().solve(gr_rhs);
    Matrix<T, Dynamic, Dynamic> data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}
//#if 0
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

    auto recdeg = di.reconstruction_degree();
    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    auto cb = make_vector_monomial_basis(msh, cl, recdeg);

    auto rbs = vector_basis_size(recdeg, N, N);
    auto cbs = vector_basis_size(celdeg, N, N);
    auto fbs = vector_basis_size(facdeg, N-1, N);

    auto num_faces = howmany_faces(msh, cl);

    typedef Matrix<T, Dynamic, 1> vector_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, N, N>             gradient_type;
    typedef Matrix<T, Dynamic, N>       function_type;

    size_t rbs_ho = rbs - N;
    size_t num_total_dofs = cbs + num_faces*fbs;
    size_t nb_lag         = priv::nb_lag(N);

    matrix_type stiff          = matrix_type::Zero(rbs, rbs);
    matrix_type gr_lhs = matrix_type::Zero( rbs_ho + nb_lag, rbs_ho + nb_lag);
    matrix_type gr_rhs = matrix_type::Zero( rbs_ho + nb_lag, num_total_dofs);

    auto qps = integrate(msh, cl, 2*recdeg);
    for (auto& qp : qps)
    {
        auto dphi = cb.eval_sgradients(qp.point());
        stiff += qp.weight() * priv::outer_product(dphi, dphi);
    }

    gr_lhs.block(0, 0, rbs_ho, rbs_ho) = stiff.block(N, N, rbs_ho, rbs_ho);
    gr_rhs.block(0, 0, rbs_ho, cbs)    = stiff.block(N, 0, rbs_ho, cbs);

    auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        auto fc = fcs[i];
        auto n = normal(msh, cl, fc);
        auto fb = make_vector_monomial_basis(msh, fc, facdeg);

        auto qps_f = integrate(msh, fc, 2*facdeg);
        for (auto& qp : qps_f)
        {
            function_type   c_phi_tmp = cb.eval_functions(qp.point());
            function_type   c_phi     = c_phi_tmp.block(0, 0, cbs, N);

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

    vector_type rot = vector_type::Zero(rbs);
    for (auto& qp : qps)
    {
        auto rphi = cb.eval_curls(qp.point());
        rot += qp.weight() * rphi;
    }

    gr_lhs.block(0, rbs_ho, rbs_ho, nb_lag ) = rot.tail(rbs_ho);
    gr_lhs.block(rbs_ho, 0, nb_lag, rbs_ho ) = rot.tail(rbs_ho).transpose();

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
make_hho_gradrec_matrix(const Mesh&                     msh,
                          const typename Mesh::cell_type& cl,
                          const hho_degree_info&          di)
{
   using T        = typename Mesh::coordinate_type;
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

   typedef Matrix<T, Dynamic, Dynamic> matrix_type;

   matrix_type gr_lhs = matrix_type::Zero(gbs, gbs);
   matrix_type gr_rhs = matrix_type::Zero(gbs, cbs + num_faces * fbs);

   const size_t dim2 = N * N;

    // this is very costly to build it
   const auto qps = integrate(msh, cl, 2 * graddeg);
   for (auto& qp : qps) {
      const auto gphi = gb.eval_functions(qp.point());

      for (size_t j = 0; j < gbs; j += dim2) {
         const auto qp_gphi_j = priv::inner_product(qp.weight(), gphi[j]);
         for (size_t i = j; i < gbs; i += dim2) {
            gr_lhs(i, j) += priv::inner_product(gphi[i], qp_gphi_j);
         }
      }
   }

   // upper part
   for (size_t j = 0; j < gbs; j++) {
      for (size_t i = 0; i < j; i++) {
         gr_lhs(i, j) = gr_lhs(j, i);
      }
   }

   // copy of each cols
   for (size_t j = 0; j < gbs; j += dim2) {
      for (size_t col = 1; col < dim2; col++) {
         gr_lhs.block(col, j + col, gbs - dim2 + 1, 1) = gr_lhs.block(0, j, gbs - dim2 + 1, 1);
      }
   }

    // compute rhs
   const auto qpc = integrate(msh, cl, std::max(int(graddeg + celdeg) -1,0));
   for (auto& qp : qpc) {
      const auto gphi = gb.eval_functions(qp.point());
      const auto dphi = cb.eval_gradients(qp.point());

      for (size_t j = 0; j < cbs; j++) {
         const auto qp_dphi_j = priv::inner_product(qp.weight(), dphi[j]);
         for (size_t i = 0; i < gbs; i++) {
            gr_rhs(i, j) += priv::inner_product(gphi[i], qp_dphi_j);
         }
      }
   } // end qp

   const auto fcs = faces(msh, cl);
   for (size_t i = 0; i < fcs.size(); i++) {
      const auto fc = fcs[i];
      const auto n  = normal(msh, cl, fc);
      const auto fb = make_vector_monomial_basis(msh, fc, facdeg);

      const auto qps_f = integrate(msh, fc, graddeg + std::max(celdeg, facdeg));
      for (auto& qp : qps_f) {
         const auto cphi = cb.eval_functions(qp.point());
         const auto gphi = gb.eval_functions(qp.point());
         const auto fphi = fb.eval_functions(qp.point());

         const auto gphi_n = priv::outer_product(gphi, n);
         gr_rhs.block(0, cbs + i * fbs, gbs , fbs) +=
           qp.weight() * priv::outer_product(fphi, gphi_n);
         gr_rhs.block(0, 0, gbs, cbs) -= qp.weight() * priv::outer_product(cphi, gphi_n);
      }
   }

   Matrix<T, Dynamic, Dynamic> oper = gr_lhs.llt().solve(gr_rhs);
   Matrix<T, Dynamic, Dynamic> data = gr_rhs.transpose() * oper;

   return std::make_pair(oper, data);
}

template<typename Mesh>
std::pair<Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
          Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>>
make_hho_sym_gradrec_matrix(const Mesh&                     msh,
                        const typename Mesh::cell_type& cl,
                        const hho_degree_info&          di)
{
   using T        = typename Mesh::coordinate_type;
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

   typedef Matrix<T, Dynamic, Dynamic> matrix_type;

   matrix_type gr_lhs = matrix_type::Zero(gbs, gbs);
   matrix_type gr_rhs = matrix_type::Zero(gbs, cbs + num_faces * fbs);

   const size_t dim2 = N * N;

   // this is very costly to build it
   const auto qps = integrate(msh, cl, 2 * graddeg);

   size_t dec = 0;
   if (N == 3) {
      dec = 6;
   } else if (N == 2) {
      dec = 3;
   } else
      std::logic_error("Expected 3 >= dim > 1");

   for (auto& qp : qps) {
      const auto gphi = gb.eval_functions(qp.point());

      for (size_t j = 0; j < gbs; j++) {
         const auto qp_gphi_j = priv::inner_product(qp.weight(), gphi[j]);
         for (size_t i = j; i < gbs; i += dec) {
            gr_lhs(i, j) += priv::inner_product(gphi[i], qp_gphi_j);
         }
      }
   }

   // upper part
   for (size_t j = 0; j < gbs; j++) {
      for (size_t i = 0; i < j; i++) {
         gr_lhs(i, j) = gr_lhs(j, i);
      }
   }

   // compute rhs
   const auto qpc = integrate(msh, cl, std::max(int(graddeg + celdeg) - 1, 0));
   for (auto& qp : qpc) {
      const auto gphi = gb.eval_functions(qp.point());
      const auto dphi = cb.eval_sgradients(qp.point());

      for (size_t j = 0; j < cbs; j++) {
         const auto qp_dphi_j = priv::inner_product(qp.weight(), dphi[j]);
         for (size_t i = 0; i < gbs; i++) {
            gr_rhs(i, j) += priv::inner_product(gphi[i], qp_dphi_j);
         }
      }
   } // end qp

   const auto fcs = faces(msh, cl);
   for (size_t i = 0; i < fcs.size(); i++) {
      const auto fc = fcs[i];
      const auto n  = normal(msh, cl, fc);
      const auto fb = make_vector_monomial_basis(msh, fc, facdeg);

      const auto qps_f = integrate(msh, fc, graddeg + std::max(celdeg, facdeg));
      for (auto& qp : qps_f) {
         const auto cphi = cb.eval_functions(qp.point());
         const auto gphi = gb.eval_functions(qp.point());
         const auto fphi = fb.eval_functions(qp.point());

         const auto gphi_n = priv::outer_product(gphi, n);
         gr_rhs.block(0, cbs + i * fbs, gbs, fbs) +=
           qp.weight() * priv::outer_product(fphi, gphi_n);
         gr_rhs.block(0, 0, gbs, cbs) -= qp.weight() * priv::outer_product(cphi, gphi_n);
      }
   }

   Matrix<T, Dynamic, Dynamic> oper = gr_lhs.llt().solve(gr_rhs);
   Matrix<T, Dynamic, Dynamic> data = gr_rhs.transpose() * oper;

   return std::make_pair(oper, data);
}

template<typename Mesh>
std::pair<   Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>  >
make_hho_divergence_reconstruction(const Mesh& msh, const typename Mesh::cell_type& cl,
                                   const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;

    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    auto cbas_v = make_vector_monomial_basis(msh, cl, celdeg);
    auto cbas_s = make_scalar_monomial_basis(msh, cl, facdeg);

    auto rbs = scalar_basis_size(facdeg, Mesh::dimension);
    auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
    auto fbs = vector_basis_size(facdeg, Mesh::dimension-1, Mesh::dimension);

    auto num_faces = howmany_faces(msh, cl);

    Matrix<T, Dynamic, Dynamic> dr_lhs = Matrix<T, Dynamic, Dynamic>::Zero(rbs, rbs);
    Matrix<T, Dynamic, Dynamic> dr_rhs = Matrix<T, Dynamic, Dynamic>::Zero(rbs, cbs + num_faces*fbs);

    auto qps = integrate(msh, cl, 2*facdeg);
    for (auto& qp : qps)
    {
        Matrix<T, Dynamic, 1> s_phi  = cbas_s.eval_functions(qp.point());
        auto s_dphi = cbas_s.eval_gradients(qp.point());
        auto v_phi  = cbas_v.eval_functions(qp.point());

        dr_lhs += qp.weight() * priv::outer_product(s_phi, s_phi);
        dr_rhs.block(0, 0, rbs, cbs) -= qp.weight() * priv::outer_product(v_phi, s_dphi);
    }

    auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        auto fc = fcs[i];
        auto n = normal(msh, cl, fc);
        auto fbas_v = make_vector_monomial_basis(msh, fc, facdeg);

        auto qps_f = integrate(msh, fc, 2*facdeg);
        for (auto& qp : qps_f)
        {
            Matrix<T, Dynamic, 1> s_phi = cbas_s.eval_functions(qp.point());
            auto f_phi = fbas_v.eval_functions(qp.point());

            Matrix<T, Dynamic, 2> s_phi_n = (s_phi * n.transpose());//priv::outer_product(s_phi, n);
            dr_rhs.block(0, cbs + i*fbs, rbs, fbs) +=
                    qp.weight() * priv::outer_product(f_phi, s_phi_n);

        }
    }

    Matrix<T, Dynamic, Dynamic> oper = dr_lhs.llt().solve(dr_rhs);
    Matrix<T, Dynamic, Dynamic> data = dr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hdg_scalar_stabilization(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;

    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);

    auto num_faces = howmany_faces(msh, cl);

    Matrix<T, Dynamic, Dynamic> data = Matrix<T, Dynamic, Dynamic>::Zero(cbs+num_faces*fbs, cbs+num_faces*fbs);
    Matrix<T, Dynamic, Dynamic> If = Matrix<T, Dynamic, Dynamic>::Identity(fbs, fbs);

    auto cb = make_scalar_monomial_basis(msh, cl, celdeg);
    auto fcs = faces(msh, cl);

    for (size_t i = 0; i < num_faces; i++)
    {
        auto fc = fcs[i];
        auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

        Matrix<T, Dynamic, Dynamic> oper = Matrix<T, Dynamic, Dynamic>::Zero(fbs, cbs+num_faces*fbs);
        Matrix<T, Dynamic, Dynamic> tr = Matrix<T, Dynamic, Dynamic>::Zero(fbs, cbs+num_faces*fbs);
        Matrix<T, Dynamic, Dynamic> mass = Matrix<T, Dynamic, Dynamic>::Zero(fbs, fbs);
        Matrix<T, Dynamic, Dynamic> trace = Matrix<T, Dynamic, Dynamic>::Zero(fbs, cbs);

        oper.block(0, cbs+i*fbs, fbs, fbs) = -If;

        auto qps = integrate(msh, fc, 2*facdeg);
        for (auto& qp : qps)
        {
            auto c_phi = cb.eval_functions(qp.point());
            auto f_phi = fb.eval_functions(qp.point());

            mass += qp.weight() * f_phi * f_phi.transpose();
            trace += qp.weight() * f_phi * c_phi.transpose();
        }

        tr.block(0, cbs+i*fbs, fbs, fbs) = -mass;
        tr.block(0, 0, fbs, cbs) = trace;
        auto h = measure(msh, fc);
        oper.block(0, 0, fbs, cbs) = mass.llt().solve(trace);
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

    auto K_TT_ldlt = K_TT.llt();
    matrix_type AL = K_TT_ldlt.solve(K_TF);
    vector_type bL = K_TT_ldlt.solve(cell_rhs);

    matrix_type AC = K_FF - K_FT * AL;
    vector_type bC = /* no projection on faces, eqn. 26*/ - K_FT * bL;

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

    vector_type solT = K_TT.llt().solve(cell_rhs - K_TF*solF);

    ret.head(cell_size)         = solT;
    ret.tail(all_faces_size)    = solF;

    return ret;
}



template<typename Mesh>
[[deprecated("Please use 'make_hho_scalar_stabilization()'")]]
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hho_fancy_stabilization(const Mesh& msh, const typename Mesh::cell_type& cl,
                             const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic> reconstruction,
                             const hho_degree_info& di)
{
    return make_hho_scalar_stabilization(msh, cl, reconstruction, di);
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hho_scalar_stabilization(const Mesh& msh, const typename Mesh::cell_type& cl,
                              const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic> reconstruction,
                              const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;

    auto recdeg = di.reconstruction_degree();
    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    auto rbs = scalar_basis_size(recdeg, Mesh::dimension);
    auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);

    auto cb = make_scalar_monomial_basis(msh, cl, recdeg);

    Matrix<T, Dynamic, Dynamic> mass_mat = Matrix<T, Dynamic, Dynamic>::Zero(rbs, rbs);
    auto cell_quadpoints = integrate(msh, cl, 2*recdeg);
    for (auto& qp : cell_quadpoints)
    {
        auto c_phi = cb.eval_functions(qp.point());
        mass_mat += qp.weight() * c_phi * c_phi.transpose();
    }

    // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

    //Step 1: compute \pi_T^k p_T^k v (third term).
    Matrix<T, Dynamic, Dynamic> M1 = mass_mat.block(0, 0, cbs, cbs);
    Matrix<T, Dynamic, Dynamic> M2 = mass_mat.block(0, 1, cbs, rbs-1);
    Matrix<T, Dynamic, Dynamic> proj1 = -M1.llt().solve(M2*reconstruction);

    //Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
    Matrix<T, Dynamic, Dynamic> I_T = Matrix<T, Dynamic, Dynamic>::Identity(cbs, cbs);
    proj1.block(0, 0, cbs, cbs) += I_T;

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    Matrix<T, Dynamic, Dynamic> data = Matrix<T, Dynamic, Dynamic>::Zero(cbs+num_faces*fbs, cbs+num_faces*fbs);

    auto h = diameter(msh, cl);

    // Step 3: project on faces (eqn. 21)
    for (size_t face_i = 0; face_i < num_faces; face_i++)
    {
        auto fc = fcs[face_i];
        auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

        Matrix<T, Dynamic, Dynamic> face_mass_matrix    = Matrix<T, Dynamic, Dynamic>::Zero(fbs, fbs);
        Matrix<T, Dynamic, Dynamic> face_trace_matrix   = Matrix<T, Dynamic, Dynamic>::Zero(fbs, rbs);

        auto face_quadpoints = integrate(msh, fc, 2*recdeg);
        for (auto& qp : face_quadpoints)
        {
            auto f_phi = fb.eval_functions(qp.point());
            auto c_phi = cb.eval_functions(qp.point());
            Matrix<T, Dynamic, 1> q_f_phi = qp.weight() * f_phi;
            face_mass_matrix += q_f_phi * f_phi.transpose();
            face_trace_matrix += q_f_phi * c_phi.transpose();
        }

        LLT<Matrix<T, Dynamic, Dynamic>> piKF;
        piKF.compute(face_mass_matrix);

        // Step 3a: \pi_F^k( v_F - p_T^k v )
        Matrix<T, Dynamic, Dynamic> MR1 = face_trace_matrix.block(0, 1, fbs, rbs-1);

        Matrix<T, Dynamic, Dynamic> proj2 = piKF.solve(MR1*reconstruction);
        Matrix<T, Dynamic, Dynamic> I_F = Matrix<T, Dynamic, Dynamic>::Identity(fbs, fbs);
        proj2.block(0, cbs+face_i*fbs, fbs, fbs) -= I_F;

        // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
        Matrix<T, Dynamic, Dynamic> MR2 = face_trace_matrix.block(0, 0, fbs, cbs);
        Matrix<T, Dynamic, Dynamic> proj3 = piKF.solve(MR2*proj1);
        Matrix<T, Dynamic, Dynamic> BRF = proj2 + proj3;

        data += BRF.transpose() * face_mass_matrix * BRF / h;
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

    auto recdeg = di.reconstruction_degree();
    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    auto rbs = scalar_basis_size(recdeg, Mesh::dimension);
    auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);

    auto cb = make_scalar_monomial_basis(msh, cl, recdeg);

    Matrix<T, Dynamic, Dynamic> mass_mat = Matrix<T, Dynamic, Dynamic>::Zero(rbs, rbs);
    auto cell_quadpoints = integrate(msh, cl, 2*recdeg);
    for (auto& qp : cell_quadpoints)
    {
        auto c_phi = cb.eval_functions(qp.point());
        mass_mat += qp.weight() * c_phi * c_phi.transpose();
    }

    // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

    //Step 1: compute \pi_T^k p_T^k v (third term).
    Matrix<T, Dynamic, Dynamic> M1 = mass_mat.block(0, 0, cbs, cbs);
    Matrix<T, Dynamic, Dynamic> M2 = mass_mat.block(0, 1, cbs, rbs-1);
    Matrix<T, Dynamic, Dynamic> proj1 = -M1.llt().solve(M2*reconstruction);

    //Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
    Matrix<T, Dynamic, Dynamic> I_T = Matrix<T, Dynamic, Dynamic>::Identity(cbs, cbs);
    proj1.block(0, 0, cbs, cbs) += I_T;

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    Matrix<T, Dynamic, Dynamic> data = Matrix<T, Dynamic, Dynamic>::Zero(cbs+num_faces*fbs, cbs+num_faces*fbs);

    auto h = diameter(msh, cl);

    // Step 3: project on faces (eqn. 21)
    for (size_t face_i = 0; face_i < num_faces; face_i++)
    {
        auto fc = fcs[face_i];
        auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

        Matrix<T, Dynamic, Dynamic> face_mass_matrix    = Matrix<T, Dynamic, Dynamic>::Zero(fbs, fbs);
        Matrix<T, Dynamic, Dynamic> face_trace_matrix   = Matrix<T, Dynamic, Dynamic>::Zero(fbs, rbs);

        auto face_quadpoints = integrate(msh, fc, 2*recdeg);
        for (auto& qp : face_quadpoints)
        {
            auto f_phi = fb.eval_functions(qp.point());
            auto c_phi = cb.eval_functions(qp.point());
            Matrix<T, Dynamic, 1> q_f_phi = qp.weight() * f_phi;
            face_mass_matrix += q_f_phi * f_phi.transpose();
            face_trace_matrix += q_f_phi * c_phi.transpose();
        }

        LLT<Matrix<T, Dynamic, Dynamic>> piKF;
        piKF.compute(face_mass_matrix);

        // Step 3a: \pi_F^k( v_F - p_T^k v )
        Matrix<T, Dynamic, Dynamic> MR1 = face_trace_matrix.block(0, 1, fbs, rbs-1);
        Matrix<T, Dynamic, Dynamic> A = MR1 * reconstruction;

        Matrix<T, Dynamic, Dynamic> proj2 = piKF.solve(A);
        Matrix<T, Dynamic, Dynamic> I_F = Matrix<T, Dynamic, Dynamic>::Identity(fbs, fbs);
        proj2.block(0, cbs+face_i*fbs, fbs, fbs) -= I_F;
        A.block(0, cbs+face_i*fbs, fbs, fbs) -= face_mass_matrix;

        // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
        Matrix<T, Dynamic, Dynamic> MR2 = face_trace_matrix.block(0, 0, fbs, cbs);
        Matrix<T, Dynamic, Dynamic> B = MR2 * proj1;
        Matrix<T, Dynamic, Dynamic> proj3 = piKF.solve(B);
        Matrix<T, Dynamic, Dynamic> BRF = proj2 + proj3;

        data += BRF.transpose() * (A+B) / h;
    }

    return data;
}



template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hho_fancy_stabilization_vector(const Mesh& msh, const typename Mesh::cell_type& cl,
                      const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic> reconstruction,
                      const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;

    auto recdeg = di.reconstruction_degree();
    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    auto rbs = vector_basis_size(recdeg, Mesh::dimension, Mesh::dimension);
    auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
    auto fbs = vector_basis_size(facdeg, Mesh::dimension-1, Mesh::dimension);

    size_t N = Mesh::dimension;

    auto cb = make_vector_monomial_basis(msh, cl, recdeg);

    Matrix<T, Dynamic, Dynamic> mass_mat = Matrix<T, Dynamic, Dynamic>::Zero(rbs, rbs);
    auto cell_quadpoints = integrate(msh, cl, 2*recdeg);
    for (auto& qp : cell_quadpoints)
    {
        auto c_phi = cb.eval_functions(qp.point());
        mass_mat += qp.weight() * priv::outer_product(c_phi, c_phi);
    }

    // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

    //Step 1: compute \pi_T^k p_T^k v (third term).
    Matrix<T, Dynamic, Dynamic> M1 = mass_mat.block(0, 0, cbs, cbs);
    Matrix<T, Dynamic, Dynamic> M2 = mass_mat.block(0, N, cbs, rbs-N);
    Matrix<T, Dynamic, Dynamic> proj1 = -M1.llt().solve(M2*reconstruction);

    //Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
    Matrix<T, Dynamic, Dynamic> I_T = Matrix<T, Dynamic, Dynamic>::Identity(cbs, cbs);
    proj1.block(0, 0, cbs, cbs) += I_T;

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    Matrix<T, Dynamic, Dynamic> data = Matrix<T, Dynamic, Dynamic>::Zero(cbs+num_faces*fbs, cbs+num_faces*fbs);

    // Step 3: project on faces (eqn. 21)
    for (size_t face_i = 0; face_i < num_faces; face_i++)
    {
        auto h = diameter(msh, fcs[face_i]);
        auto fc = fcs[face_i];
        auto fb = make_vector_monomial_basis(msh, fc, facdeg);

        Matrix<T, Dynamic, Dynamic> face_mass_matrix    = Matrix<T, Dynamic, Dynamic>::Zero(fbs, fbs);
        Matrix<T, Dynamic, Dynamic> face_trace_matrix   = Matrix<T, Dynamic, Dynamic>::Zero(fbs, rbs);

        auto face_quadpoints = integrate(msh, fc, 2*recdeg);
        for (auto& qp : face_quadpoints)
        {
            auto f_phi = fb.eval_functions(qp.point());
            auto c_phi = cb.eval_functions(qp.point());
            Matrix<T, Dynamic, Mesh::dimension> q_f_phi = qp.weight() * f_phi;
            face_mass_matrix += priv::outer_product(f_phi, q_f_phi);
            face_trace_matrix += priv::outer_product(c_phi, q_f_phi);
        }

        LLT<Matrix<T, Dynamic, Dynamic>> piKF;
        piKF.compute(face_mass_matrix);

        // Step 3a: \pi_F^k( v_F - p_T^k v )
        Matrix<T, Dynamic, Dynamic> MR1 = face_trace_matrix.block(0, N, fbs, rbs-N);

        Matrix<T, Dynamic, Dynamic> proj2 = piKF.solve(MR1*reconstruction);
        Matrix<T, Dynamic, Dynamic> I_F = Matrix<T, Dynamic, Dynamic>::Identity(fbs, fbs);
        proj2.block(0, cbs+face_i*fbs, fbs, fbs) -= I_F;

        // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
        Matrix<T, Dynamic, Dynamic> MR2 = face_trace_matrix.block(0, 0, fbs, cbs);
        Matrix<T, Dynamic, Dynamic> proj3 = piKF.solve(MR2*proj1);
        Matrix<T, Dynamic, Dynamic> BRF = proj2 + proj3;

        data += BRF.transpose() * face_mass_matrix * BRF / h;
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

    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto num_faces = howmany_faces(msh, cl);

    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs+num_faces*fbs);

    auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
    Matrix<T, Dynamic, Dynamic> cell_mm = make_mass_matrix(msh, cl, cb, di);
    Matrix<T, Dynamic, 1> cell_rhs = make_rhs(msh, cl, cb, f, di);
    ret.block(0, 0, cbs, 1) = cell_mm.llt().solve(cell_rhs);

    auto fcs = faces(msh, cl);
    for (size_t i = 0; i < num_faces; i++)
    {
        auto fc = fcs[i];
        auto fb = make_scalar_monomial_basis(msh, fc, hdi.face_degree());
        Matrix<T, Dynamic, Dynamic> face_mm = make_mass_matrix(msh, fc, fb, di);
        Matrix<T, Dynamic, 1> face_rhs = make_rhs(msh, fc, fb, f, di);
        ret.block(cbs+i*fbs, 0, fbs, 1) = face_mm.llt().solve(face_rhs);
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

    auto cbs = vector_basis_size(hdi.cell_degree(), Mesh::dimension, Mesh::dimension);
    auto fbs = vector_basis_size(hdi.face_degree(), Mesh::dimension-1, Mesh::dimension);
    auto num_faces = howmany_faces(msh, cl);

    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs+num_faces*fbs);

    auto cb = make_vector_monomial_basis(msh, cl, hdi.cell_degree());
    Matrix<T, Dynamic, Dynamic> cell_mm = make_mass_matrix(msh, cl, cb, di);
    Matrix<T, Dynamic, 1> cell_rhs = make_rhs(msh, cl, cb, f, di);
    ret.block(0, 0, cbs, 1) = cell_mm.llt().solve(cell_rhs);

    auto fcs = faces(msh, cl);
    for (size_t i = 0; i < num_faces; i++)
    {
        auto fc = fcs[i];
        auto fb = make_vector_monomial_basis(msh, fc, hdi.face_degree());
        Matrix<T, Dynamic, Dynamic> face_mm = make_mass_matrix(msh, fc, fb, di);
        Matrix<T, Dynamic, 1> face_rhs = make_rhs(msh, fc, fb, f, di);
        ret.block(cbs+i*fbs, 0, fbs, 1) = face_mm.llt().solve(face_rhs);
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

    auto fbs       = vector_basis_size(hdi.face_degree(), Mesh::dimension - 1, Mesh::dimension);

    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(fbs);

    auto                        fb       = make_vector_monomial_basis(msh, fc, hdi.face_degree());
    Matrix<T, Dynamic, Dynamic> face_mm  = make_mass_matrix(msh, fc, fb, di);
    Matrix<T, Dynamic, 1>       face_rhs = make_rhs(msh, fc, fb, f, di);
    ret.block(0, 0, fbs, 1)  = face_mm.llt().solve(face_rhs);

   return ret;
}

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

    SparseMatrix<T>         LHS;
    Matrix<T, Dynamic, 1>   RHS;

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
        RHS = Matrix<T, Dynamic, 1>::Zero( system_size );
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
             const Matrix<T, Dynamic, Dynamic>& lhs,
             const Matrix<T, Dynamic, 1>& rhs,
             const Function& dirichlet_bf)
    {
        auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
            return msh.is_boundary(fc);
        };

        auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);

        auto fcs = faces(msh, cl);

        std::vector<assembly_index> asm_map;
        asm_map.reserve(fcs.size()*fbs);

        Matrix<T, Dynamic, 1> dirichlet_data = Matrix<T, Dynamic, 1>::Zero(fcs.size()*fbs);

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

                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_bf, di.face_degree());
                dirichlet_data.block(face_i*fbs, 0, fbs, 1) = mass.llt().solve(rhs);
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
    Matrix<T, Dynamic, 1>
    take_local_data(const Mesh& msh, const typename Mesh::cell_type& cl,
    const Matrix<T, Dynamic, 1>& solution, const Function& dirichlet_bf)
    {
        auto facdeg = di.face_degree();
        auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);
        auto fcs = faces(msh, cl);

        auto num_faces = fcs.size();

        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(num_faces*fbs);

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

                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_bf, di.face_degree());
                ret.block(face_i*fbs, 0, fbs, 1) = mass.llt().solve(rhs);
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
class stokes_assembler
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

    SparseMatrix<T>         LHS;
    Matrix<T, Dynamic, 1>   RHS;

    stokes_assembler(const Mesh& msh, hho_degree_info hdi)
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

        auto cbs_A = vector_basis_size(hdi.cell_degree(), Mesh::dimension, Mesh::dimension);
        auto fbs_A = vector_basis_size(hdi.face_degree(), Mesh::dimension-1, Mesh::dimension);
        auto cbs_B = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

        auto system_size = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * msh.cells_size() + 1;

        LHS = SparseMatrix<T>( system_size, system_size );
        RHS = Matrix<T, Dynamic, 1>::Zero( system_size );
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
             const Matrix<T, Dynamic, Dynamic>& lhs_A,
             const Matrix<T, Dynamic, Dynamic>& lhs_B,
             const Matrix<T, Dynamic, 1>& rhs,
             const Function& dirichlet_bf)
    {
        auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
            return msh.is_boundary(fc);
        };

        auto cbs_A = vector_basis_size(di.cell_degree(), Mesh::dimension, Mesh::dimension);
        auto fbs_A = vector_basis_size(di.face_degree(), Mesh::dimension-1, Mesh::dimension);
        auto cbs_B = scalar_basis_size(di.cell_degree(), Mesh::dimension);

        auto fcs = faces(msh, cl);

        std::vector<assembly_index> asm_map;
        asm_map.reserve(cbs_A + fcs.size()*fbs_A);

        auto cell_offset        = priv::offset(msh, cl);
        auto cell_LHS_offset    = cell_offset * cbs_A;

        auto B_offset = cbs_A * msh.cells_size() + fbs_A * num_other_faces + cbs_B * cell_offset;

        for (size_t i = 0; i < cbs_A; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );

        Matrix<T, Dynamic, 1> dirichlet_data = Matrix<T, Dynamic, 1>::Zero(cbs_A + fcs.size()*fbs_A);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto face_offset = priv::offset(msh, fc);
            auto face_LHS_offset = cbs_A * msh.cells_size() + compress_table.at(face_offset)*fbs_A;

            bool dirichlet = is_dirichlet(fc);

            for (size_t i = 0; i < fbs_A; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );

            if (dirichlet)
            {
                auto fb = make_vector_monomial_basis(msh, fc, di.face_degree());

                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_bf, di.face_degree());
                dirichlet_data.block(cbs_A + face_i*fbs_A, 0, fbs_A, 1) = mass.llt().solve(rhs);
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

        auto scalar_cell_basis = make_scalar_monomial_basis(msh, cl, di.cell_degree());
        auto qps = integrate(msh, cl, di.cell_degree());
        Matrix<T, Dynamic, 1> mult = Matrix<T, Dynamic, 1>::Zero( scalar_cell_basis.size() );
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

#if 0
    template<typename Function>
    Matrix<T, Dynamic, 1>
    take_local_data(const Mesh& msh, const typename Mesh::cell_type& cl,
    const Matrix<T, Dynamic, 1>& solution, const Function& dirichlet_bf)
    {
        auto celdeg = di.cell_degree();
        auto facdeg = di.face_degree();

        auto cbs = cell_basis<Mesh,T>::size(celdeg);
        auto fbs = face_basis<Mesh,T>::size(facdeg);

        auto cell_offset        = offset(msh, cl);
        auto cell_SOL_offset    = cell_offset * cbs;

        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs + 4*fbs);
        ret.block(0, 0, cbs, 1) = solution.block(cell_SOL_offset, 0, cbs, 1);

        auto fcs = faces(msh, cl);
        for (size_t face_i = 0; face_i < 4; face_i++)
        {
            auto fc = fcs[face_i];

            bool dirichlet = fc.is_boundary && fc.bndtype == boundary::DIRICHLET;

            if (dirichlet)
            {
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, facdeg);
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, facdeg, dirichlet_bf);
                ret.block(cbs+face_i*fbs, 0, fbs, 1) = mass.llt().solve(rhs);
            }
            else
            {
                auto face_offset = offset(msh, fc);
                auto face_SOL_offset = cbs * msh.cells.size() + compress_table.at(face_offset)*fbs;
                ret.block(cbs+face_i*fbs, 0, fbs, 1) = solution.block(face_SOL_offset, 0, fbs, 1);
            }
        }

        return ret;
    }
#endif

    void finalize(void)
    {
        LHS.setFromTriplets( triplets.begin(), triplets.end() );
        triplets.clear();
    }
    size_t num_assembled_faces() const
    {
        return num_other_faces;
    }

};



template<typename Mesh>
auto make_stokes_assembler(const Mesh& msh, hho_degree_info hdi)
{
    return stokes_assembler<Mesh>(msh, hdi);
}

template<typename Mesh>
class assembler_mechanics
{
   typedef disk::mechanics::BoundaryConditions<Mesh>    bnd_type;
   typedef Mesh                                         mesh_type;
   typedef typename mesh_type::scalar_type              scalar_type;
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

   template<typename LocalContrib>
   void
   assemble(const mesh_type& msh, const cell_type& cl, const bnd_type& bnd, const LocalContrib& lc)
   {
      const size_t      face_degree   = m_hdi.face_degree();
      const auto        num_face_dofs = vector_basis_size(face_degree, dimension - 1, dimension);
      const scalar_type zero          = 0;

      const auto          fcs = faces(msh, cl);
      std::vector<size_t> l2g(fcs.size() * num_face_dofs);
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

            vector_type proj_bcf =
              project_function(msh, fc, face_degree, bnd.dirichlet_boundary_func(face_id));

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
                           l2g.at(pos + i + 1) = face_offset + ind_sol++;
                           l2g.at(pos + i + 2) = 0xDEADBEEF;
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
         if (l2g[j] == 0xDEADBEEF) continue;

         for (size_t i = 0; i < lc.first.rows(); i++) {
            if (l2g[i] == 0xDEADBEEF) continue;

            m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));
         }
         RHS(l2g.at(j)) += lc.second(j) - rhs_bc(j);
      }
#else
      for (size_t i = 0; i < lc.first.rows(); i++) {
         if (l2g[i] == 0xDEADBEEF) continue;

         for (size_t j = 0; j < lc.first.cols(); j++) {
            if (l2g[j] == 0xDEADBEEF) continue;

            m_triplets.push_back(triplet_type(l2g.at(i), l2g.at(j), lc.first(i, j)));
         }
         RHS(l2g.at(i)) += lc.second(i) - rhs_bc(i);
      }
#endif
   }

   vector_type
   expand_solution(const mesh_type& msh, const bnd_type& bnd, const vector_type& solution)
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
              project_function(msh, bfc, face_degree, bnd.dirichlet_boundary_func(face_id));

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
   typedef Mesh                            mesh_type;
   typedef typename mesh_type::scalar_type scalar_type;
   typedef typename mesh_type::cell        cell_type;

   typedef dynamic_matrix<scalar_type> matrix_type;
   typedef dynamic_vector<scalar_type> vector_type;

   const static size_t dimension = mesh_type::dimension;

 public:
   matrix_type AL;
   vector_type bL;
   static_condensation_vector() {}

   std::pair<matrix_type, vector_type>
   compute(const mesh_type&   msh,
           const cell_type&   cl,
           const matrix_type& local_mat,
           const vector_type& cell_rhs,
           const hho_degree_info hdi)
   {
      const size_t num_cell_dofs = vector_basis_size(hdi.cell_degree(), dimension, dimension);
      const size_t num_face_dofs = vector_basis_size(hdi.face_degree(), dimension - 1, dimension);

      const auto   fcs       = faces(msh, cl);
      const size_t num_faces = fcs.size();

      disk::dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

      const size_t cell_size = dsr.cell_range().size();
      const size_t face_size = dsr.all_faces_range().size();

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

      const auto K_TT_ldlt = K_TT.llt();
      AL                   = K_TT_ldlt.solve(K_TF);
      bL                   = K_TT_ldlt.solve(cell_rhs);

      const matrix_type AC = K_FF - K_FT * AL;
      const vector_type bC = /* no projection on faces, eqn. 26*/ -K_FT * bL;

      return std::make_pair(AC, bC);
   }
};

} // revolution
