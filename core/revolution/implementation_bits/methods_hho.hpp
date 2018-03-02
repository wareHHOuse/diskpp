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

        std::cout << cell_deg << " " << face_deg << " " << reconstruction_deg << std::endl;
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

       std::cout << cell_deg << " " << face_deg << " " << reconstruction_deg << " "
                 << grad_deg << std::endl;
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
};

template<typename Mesh>
std::pair<   Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>  >
make_hho_scalar_laplacian(const Mesh& msh, const typename Mesh::cell_type& cl,
                          const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;

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
            Matrix<T, Dynamic, 2> c_dphi_tmp = cb.eval_gradients(qp.point());
            Matrix<T, Dynamic, 2> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, 2);
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
    matrix_type gr_lhs = matrix_type::Zero(rbs-2, rbs-2);
    matrix_type gr_rhs = matrix_type::Zero(rbs-2, cbs + num_faces*fbs);

    auto qps = integrate(msh, cl, 2*recdeg);
    for (auto& qp : qps)
    {
        auto dphi = cb.eval_gradients(qp.point());
        stiff += qp.weight() * priv::outer_product(dphi, dphi);
    }

    gr_lhs = stiff.block(2, 2, rbs-2, rbs-2);
    gr_rhs.block(0, 0, rbs-2, cbs) = stiff.block(2, 0, rbs-2, cbs);

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
            function_type   c_phi = c_phi_tmp.block(0, 0, cbs, 2);

            std::vector<gradient_type> c_dphi_tmp = cb.eval_gradients(qp.point());

            auto begin_iter = std::next(c_dphi_tmp.begin(), 2);
            std::vector<gradient_type> c_dphi;
            c_dphi.resize(rbs - 2);
            assert( std::distance(begin_iter, c_dphi_tmp.end()) == c_dphi.size() );
            std::copy(begin_iter, c_dphi_tmp.end(), c_dphi.begin());

            function_type   f_phi = fb.eval_functions(qp.point());

            Matrix<T, Dynamic, N> c_dphi_n = priv::outer_product(c_dphi, n);
            gr_rhs.block(0, cbs + i*fbs, rbs-2, fbs) +=
                    qp.weight() * priv::outer_product(f_phi, c_dphi_n);
            gr_rhs.block(0, 0, rbs-2, cbs) -=
                    qp.weight() * priv::outer_product(c_phi, c_dphi_n);
        }
    }

    Matrix<T, Dynamic, Dynamic> oper = gr_lhs.llt().solve(gr_rhs);
    Matrix<T, Dynamic, Dynamic> data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}
//#if 0
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

    auto rbs = vector_basis_size(recdeg, Mesh::dimension, Mesh::dimension);
    auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
    auto fbs = vector_basis_size(facdeg, Mesh::dimension-1, Mesh::dimension);

    auto num_faces = howmany_faces(msh, cl);

    typedef Matrix<T, Dynamic, 1> vector_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, N, N>             gradient_type;
    typedef Matrix<T, Dynamic, N>       function_type;

    size_t rbs_ho = rbs - Mesh::dimension;
    size_t num_total_dofs = cbs + num_faces*fbs;
    matrix_type stiff  = matrix_type::Zero(rbs, rbs);
    matrix_type gr_lhs = matrix_type::Zero( rbs_ho + 1, rbs_ho + 1);
    matrix_type gr_rhs = matrix_type::Zero( rbs_ho + 1, num_total_dofs);

    auto qps = integrate(msh, cl, 2*recdeg);
    for (auto& qp : qps)
    {
        auto dphi = cb.eval_sgradients(qp.point());
        stiff += qp.weight() * priv::outer_product(dphi, dphi);
    }

    gr_lhs.block(0, 0, rbs_ho, rbs_ho) = stiff.block(2, 2, rbs_ho, rbs_ho);
    gr_rhs.block(0, 0, rbs_ho, cbs) = stiff.block(2, 0, rbs_ho, cbs);

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
            function_type   c_phi = c_phi_tmp.block(0, 0, cbs, 2);

            std::vector<gradient_type> c_dphi_tmp = cb.eval_sgradients(qp.point());

            auto begin_iter = std::next(c_dphi_tmp.begin(), 2);
            std::vector<gradient_type> c_dphi(rbs_ho);
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

    gr_lhs.block(0, rbs_ho, rbs_ho, 1 ) = rot.tail(rbs_ho);
    gr_lhs.block(rbs_ho, 0, 1, rbs_ho ) = rot.tail(rbs_ho).transpose();

    matrix_type sol  = gr_lhs.ldlt().solve(gr_rhs);
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

    auto recdeg = di.cell_degree();
    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    auto cbas_v = make_vector_monomial_basis(msh, cl, celdeg);
    auto cbas_s = make_scalar_monomial_basis(msh, cl, recdeg);

    auto rbs = scalar_basis_size(recdeg, Mesh::dimension);
    auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
    auto fbs = vector_basis_size(facdeg, Mesh::dimension-1, Mesh::dimension);

    auto num_faces = howmany_faces(msh, cl);

    Matrix<T, Dynamic, Dynamic> dr_lhs = Matrix<T, Dynamic, Dynamic>::Zero(rbs, rbs);
    Matrix<T, Dynamic, Dynamic> dr_rhs = Matrix<T, Dynamic, Dynamic>::Zero(rbs, cbs + num_faces*fbs);

    Matrix<T, Dynamic, 1> avgs = compute_averages(msh, cl, cbas_s);

    auto qps = integrate(msh, cl, 2*recdeg);
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

    Matrix<T, Dynamic, Dynamic> oper = dr_rhs;//dr_lhs.llt().solve(dr_rhs);
    Matrix<T, Dynamic, Dynamic> data = dr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hho_naive_stabilization(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
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
    auto h = measure(msh, cl);

    for (size_t i = 0; i < num_faces; i++)
    {
        auto fc = fcs[i];
        auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

        Matrix<T, Dynamic, Dynamic> oper = Matrix<T, Dynamic, Dynamic>::Zero(fbs, cbs+num_faces*fbs);
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

        oper.block(0, 0, fbs, cbs) = mass.llt().solve(trace);

        data += oper.transpose() * mass * oper * (1./h);
    }

    return data;
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hho_fancy_stabilization(const Mesh& msh, const typename Mesh::cell_type& cl,
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
    auto cell_quadpoints = integrate(msh, cl, 2*recdeg+2);
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

    // Step 3: project on faces (eqn. 21)
    for (size_t face_i = 0; face_i < num_faces; face_i++)
    {
        auto h = diameter(msh, /*fcs[face_i]*/cl);
        auto fc = fcs[face_i];
        auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

        Matrix<T, Dynamic, Dynamic> face_mass_matrix    = Matrix<T, Dynamic, Dynamic>::Zero(fbs, fbs);
        Matrix<T, Dynamic, Dynamic> face_trace_matrix   = Matrix<T, Dynamic, Dynamic>::Zero(fbs, rbs);

        auto face_quadpoints = integrate(msh, fc, 2*recdeg+2);
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
    Matrix<T, Dynamic, Dynamic> M2 = mass_mat.block(0, 2, cbs, rbs-2);
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
        auto h = diameter(msh, /*fcs[face_i]*/cl);
        auto fc = fcs[face_i];
        auto fb = make_vector_monomial_basis(msh, fc, facdeg);

        Matrix<T, Dynamic, Dynamic> face_mass_matrix    = Matrix<T, Dynamic, Dynamic>::Zero(fbs, fbs);
        Matrix<T, Dynamic, Dynamic> face_trace_matrix   = Matrix<T, Dynamic, Dynamic>::Zero(fbs, rbs);

        auto face_quadpoints = integrate(msh, fc, 2*facdeg);
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
        Matrix<T, Dynamic, Dynamic> MR1 = face_trace_matrix.block(0, 2, fbs, rbs-2);

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
using SRT = typename Mesh::coordinate_type;

template<typename Mesh>
using VRT = Matrix<typename Mesh::coordinate_type, Mesh::dimension, 1>;

template<typename Mesh>
using PT = typename Mesh::point_type;


template<typename Mesh>
using scalar_rhs_function = std::function<SRT<Mesh>(PT<Mesh>)>;

template<typename Mesh>
using vector_rhs_function = std::function<VRT<Mesh>(PT<Mesh>)>;

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
    return stokes_assembler_temp<Mesh>(msh, hdi);
}



} // revolution
