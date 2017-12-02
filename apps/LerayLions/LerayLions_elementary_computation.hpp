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
#include "hho/hho_bq.hpp"
#include "bases/bases_utils.hpp"
#include "BehaviorLaws/behaviorlaws.hpp"
#include "BehaviorLaws/maths_tensor.hpp"
#include "timecounter.h"

//#define USE_BLAS

namespace LerayLions {


   template<typename BQData>
   class LerayLions
   {
      typedef typename BQData::mesh_type          mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;

      typedef dynamic_matrix<scalar_type>         matrix_type;
      typedef dynamic_vector<scalar_type>         vector_type;

      const BQData&                               m_bqd;


   public:
      matrix_type     K_int;
      vector_type     RTF;
      double time_law;

      LerayLions(const BQData& bqd) : m_bqd(bqd)
      {}

      template<typename Function>
      void
      compute(const mesh_type& msh, const cell_type& cl, const Function& load, const matrix_type& GT,
              const vector_type& uTF, const scalar_type leray_param)
      {
         const size_t DIM= msh.dimension;
         const size_t cell_degree = m_bqd.cell_degree();
         const size_t face_degree = m_bqd.face_degree();
         const size_t grad_degree = m_bqd.grad_degree();
         const size_t cell_basis_size = (m_bqd.cell_basis.range(0, cell_degree)).size();
         const size_t grad_basis_size = (m_bqd.grad_basis.range(0, grad_degree)).size();
         const size_t face_basis_size = (m_bqd.face_basis.range(0, face_degree)).size();


         const size_t cpk = binomial(cell_degree + DIM, cell_degree);
         const size_t fpk = binomial(face_degree  + DIM -1, face_degree);
         const size_t gpk = DIM * binomial(grad_degree + DIM, grad_degree );

         assert(cell_basis_size == cpk);
         assert(grad_basis_size == gpk);
         assert(face_basis_size == fpk);

         time_law = 0.0;
         timecounter tc;

         auto fcs = faces(msh, cl);
         const size_t num_faces = fcs.size();
         const size_t num_total_dofs = cell_basis_size + num_faces * face_basis_size;


         matrix_type AT = matrix_type::Zero(grad_basis_size, grad_basis_size);
         vector_type aT = vector_type::Zero(grad_basis_size);

         RTF = vector_type::Zero(num_total_dofs);

         const vector_type GT_uTF = GT * uTF;

         pLaplaceLaw<scalar_type>  law(leray_param);

         auto grad_quadpoints = m_bqd.grad_quadrature.integrate(msh, cl);

         for (auto& qp : grad_quadpoints)
         {
            auto c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);
            auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());

            // Compute local gradient and norm
            auto GT_iqn = disk::compute_gradient_matrix_pt(GT_uTF, gphi);

            //Compute bahavior
            tc.tic();
            auto tensor_behavior = law.compute_whole(GT_iqn);
            tc.toc();
            time_law += tc.to_double();

            for(std::size_t i = 0; i < grad_basis_size; i++) {
               auto Agphi_i = disk::mm_prod(tensor_behavior.second, gphi[i]);
               for(std::size_t j = i; j < grad_basis_size; j++) {
                  //compute (Gkt v, A(u) : Gkt du)
                  AT(i,j) += qp.weight() * disk::mm_prod(Agphi_i, gphi[j]);
               }
            }

               // compute (PK1(u), G^k_T v)_T
            for(std::size_t i = 0; i < grad_basis_size; i += DIM){
               size_t row = i;
               for(size_t k = 0; k < DIM; k++){
                  aT(row) += qp.weight() * tensor_behavior.first(k) * gphi[row](k);
                  row++;
               }
            }

            //compute (f,v)_T
            for(std::size_t i = 0; i < cell_basis_size; i++) {
               RTF(i) += qp.weight() * disk::mm_prod(load(qp.point()) , c_phi[i]);
            }
         }

         //lower part AT
         for(std::size_t i = 1; i < grad_basis_size; i++)
            for(std::size_t j = 0; j < i; j++)
               AT(i,j) = AT(j,i);

         K_int = GT.transpose() * AT * GT;
         RTF -= GT.transpose() * aT;

         assert(K_int.rows() == num_total_dofs);
         assert(K_int.cols() == num_total_dofs);
         assert(RTF.rows() == num_total_dofs);
      }
};

}//end namespace LerayLions
