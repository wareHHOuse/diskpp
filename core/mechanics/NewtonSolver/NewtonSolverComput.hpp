/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework.
 * M. Abbas, A. Ern, N. Pignet.
 * International Journal of Numerical Methods in Engineering (2019)
 * 120(3), 303-327
 * DOI: 10.1002/nme.6137
 */

#pragma once

#include <cassert>

#include "bases/bases.hpp"
#include "common/eigen.hpp"
#include "mechanics/behaviors/laws/behaviorlaws.hpp"
#include "mechanics/deformation_tensors.hpp"
#include "methods/hho"
#include "quadratures/quadratures.hpp"

#include "timecounter.h"

namespace disk
{

namespace mechanics
{

template<typename MeshType>
class mechanical_computation
{
    typedef MeshType                            mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell            cell_type;

    typedef NewtonSolverParameter<scalar_type>    param_type;
    typedef vector_boundary_conditions<mesh_type> bnd_type;
    typedef Behavior<mesh_type>                   behavior_type;

    typedef static_matrix<scalar_type, mesh_type::dimension, mesh_type::dimension> static_matrix_type;
    typedef static_tensor<scalar_type, mesh_type::dimension>                       static_tensor_type;

    const static int dimension = mesh_type::dimension;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    bool two_dim;

    /**
     * @brief Compute A : gphi for finite deformation
     *
     * @param A fourth-order tangent modulus
     * @param gphi set of basis function for gradient reconstruction
     * @return eigen_compatible_stdvector<static_matrix_type> A : gphi
     */
    eigen_compatible_stdvector<static_matrix_type>
    compute_A_gphi(const static_tensor_type& A, const eigen_compatible_stdvector<static_matrix_type>& gphi) const
    {
        const auto grad_basis_size = gphi.size();
        const auto DIM2            = dimension * dimension;

        disk::eigen_compatible_stdvector<static_matrix_type> Aphi;
        Aphi.reserve(grad_basis_size);

        // poly classique
        for (size_t i = 0; i < grad_basis_size; i += DIM2)
        {
            size_t row = i;
            for (size_t k = 0; k < dimension; k++)
            { // depend de l'ordre des bases
                for (size_t l = 0; l < dimension; l++)
                { // depend de l'ordre des bases
                    Aphi.push_back(disk::tm_prod(A, gphi[row], l, k));
                    row++;
                }
            }
        }

        return Aphi;
    }

    /**
     * @brief Compute the stress tensor and the tangent modulus at the given qudature point
     *  - in small deformation Cauchy and Cep
     *  - in finite deformation PK1 and A
     * @tparam QPType Type of law at the quadrature point
     * @tparam LawData Type for material data
     * @param qp A behavior quadrature point where we compute the constitutive law
     * @param RkT_iqn Reconstruction operator evaluated at the quadrature point (symmetric gradient)
     * in small deformation and gradient in finite deformation)
     * @param small_def small deformation yes or no
     * @return std::pair<static_matrix_type, static_tensor_type> stress tensor and the tangent modulus
     */

    std::pair<static_matrix_type, static_tensor_type>
    compute_behavior(behavior_type&            behavior,
                     const size_t&             cell_id,
                     const size_t&             qp_id,
                     const static_matrix_type& RkT_iqn,
                     const bool                small_def) const
    {
        if (small_def)
        {
            return behavior.compute_whole(cell_id, qp_id, RkT_iqn, false);
        }
        else
        {
            const auto FT_iqn = convertGtoF(RkT_iqn);
            return behavior.compute_whole(cell_id, qp_id, FT_iqn, false);
        }
    }

    void
    compute_rigidity(const static_tensor_type&                             Cep,
                     const eigen_compatible_stdvector<static_matrix_type>& gphi,
                     const bool&                                           small_def,
                     const scalar_type                                     weight,
                     const size_t                                          grad_dim_dofs,
                     const size_t                                          grad_basis_size,
                     matrix_type&                                          AT) const
    {
        // std::cout << "module : " << Cep << std::endl;
        if (small_def)
        {
            // upper part
            for (size_t j = 0; j < grad_basis_size; j++)
            {
                const static_matrix_type Cgphi_j = weight * tm_prod(Cep, gphi[j]);
                // std::cout << j << std::endl;
                // std::cout << gphi[j] << std::endl;
                // std::cout << Cgphi_j << std::endl;
                for (size_t i = 0; i <= j; i += grad_dim_dofs)
                {
                    // compute (Ekt v, C(u) : Ekt du)
                    if (two_dim)
                    {
                        AT(i, j) += Cgphi_j(0, 0) * gphi[i](0, 0);
                        AT(i + 1, j) += 2 * Cgphi_j(0, 1) * gphi[i + 1](0, 1);
                        AT(i + 2, j) += Cgphi_j(1, 1) * gphi[i + 2](1, 1);
                    }
                    else
                    {
                        AT(i, j) += Cgphi_j(0, 0) * gphi[i](0, 0);
                        AT(i + 1, j) += 2 * Cgphi_j(0, 1) * gphi[i + 1](0, 1);
                        AT(i + 2, j) += Cgphi_j(1, 1) * gphi[i + 2](1, 1);
                        AT(i + 3, j) += 2 * Cgphi_j(0, 2) * gphi[i + 3](0, 2);
                        AT(i + 4, j) += 2 * Cgphi_j(1, 2) * gphi[i + 4](1, 2);
                        AT(i + 5, j) += Cgphi_j(2, 2) * gphi[i + 5](2, 2);
                    }
                    // AT(i, j) += disk::mm_prod(gphi[i], Agphi_j);
                }
            }
        }
        else
        {
            // lower part
            const auto qp_A_gphi = compute_A_gphi(weight * Cep, gphi);

            for (size_t j = 0; j < grad_basis_size; j += grad_dim_dofs)
            {
                size_t col = j;
                for (size_t k = 0; k < dimension; k++)
                { // depend de l'ordre des bases
                    for (size_t l = 0; l < dimension; l++)
                    { // depend de l'ordre des bases
                        for (size_t i = col; i < grad_basis_size; i++)
                        {
                            AT(i, col) += qp_A_gphi[i](l, k) * gphi[col](l, k);
                        }
                        col++;
                    }
                }
            }
        }
    }

    void
    symmetrized_rigidity_matrix(const size_t grad_basis_size, matrix_type& AT, const bool small_def) const
    {
        if (small_def)
        {
            // lower part AT
            for (size_t j = 0; j < grad_basis_size; j++)
                for (size_t i = j; i < grad_basis_size; i++)
                    AT(i, j) = AT(j, i);
        }
        else
        {
            // upper part AT
            for (size_t i = 0; i < grad_basis_size; i++)
                for (size_t j = i; j < grad_basis_size; j++)
                    AT(i, j) = AT(j, i);
        }
    }

    template<typename Function>
    void
    compute_external_forces(const mesh_type& msh,
                            const cell_type& cl,
                            const Function&  load,
                            const size_t     cell_degree,
                            vector_type&     RTF)
    {
        // compute (f,v)_T
        const auto cb       = make_vector_monomial_basis(msh, cl, cell_degree);
        RTF.head(cb.size()) = make_rhs(msh, cl, cb, load);
    }

    void
    compute_internal_forces(const static_matrix_type&                             stress,
                            const eigen_compatible_stdvector<static_matrix_type>& gphi,
                            const bool&                                           small_def,
                            const scalar_type                                     weight,
                            const size_t                                          grad_dim_dofs,
                            const size_t                                          grad_basis_size,
                            vector_type&                                          aT) const
    {
        //   std::cout << "stress" << std::endl;
        //   std::cout << stress << std::endl;
        const static_matrix_type stress_qp = weight * stress;

        if (small_def)
        {
            // compute (sigma(u), E^k_T v)_T
            for (size_t i = 0; i < grad_basis_size; i += grad_dim_dofs)
            {
                if (two_dim)
                {
                    aT(i) += stress_qp(0, 0) * gphi[i](0, 0);
                    aT(i + 1) += 2 * stress_qp(0, 1) * gphi[i + 1](0, 1);
                    aT(i + 2) += stress_qp(1, 1) * gphi[i + 2](1, 1);
                }
                else
                {
                    aT(i) += stress_qp(0, 0) * gphi[i](0, 0);
                    aT(i + 1) += 2 * stress_qp(0, 1) * gphi[i + 1](0, 1);
                    aT(i + 2) += stress_qp(1, 1) * gphi[i + 2](1, 1);
                    aT(i + 3) += 2 * stress_qp(0, 2) * gphi[i + 3](0, 2);
                    aT(i + 4) += 2 * stress_qp(1, 2) * gphi[i + 4](1, 2);
                    aT(i + 5) += stress_qp(2, 2) * gphi[i + 5](2, 2);
                }
            }
        }
        else
        {
            // compute (PK1(u), G^k_T v)_T
            for (size_t i = 0; i < grad_basis_size; i += grad_dim_dofs)
            {
                size_t row = i;
                for (size_t k = 0; k < dimension; k++)
                { // depend de l'ordre des bases
                    for (size_t l = 0; l < dimension; l++)
                    { // depend de l'ordre des bases
                        // compute (PK1(u), G^k_T v)_T
                        aT(row) += stress_qp(l, k) * gphi[row](l, k);
                        row++;
                    }
                }
            }
        }
    }

    void
    compute_contact_terms(void) const
    {
        // to fill in the future
    }

    size_t
    num_grad_dim_dofs(const bool& small_def) const
    {
        if (small_def)
        {
            if (two_dim)
                return 3;
            else
                return 6;
        }
        else
        {
            if (two_dim)
                return 4;
            else
                return 9;
        }
    }

  public:
    matrix_type K_int;
    vector_type RTF;
    vector_type F_int;
    double      time_law;
    double      time_contact;

    mechanical_computation(void)
    {
        if (dimension == 2)
            two_dim = true;
        else if (dimension == 3)
            two_dim = false;
        else
            assert(false);
    }

    template<typename Function>
    void
    compute(const mesh_type&                 msh,
            const cell_type&                 cl,
            const bnd_type&                  bnd,
            const param_type&                rp,
            const MeshDegreeInfo<mesh_type>& degree_infos,
            const Function&                  load,
            const matrix_type&               RkT,
            const vector_type&               uTF,
            behavior_type&                   behavior,
            const bool                       small_def)
    {
        time_law     = 0.0;
        time_contact = 0.0;
        timecounter tc;

        const auto cell_infos  = degree_infos.cellDegreeInfo(msh, cl);
        const auto faces_infos = cell_infos.facesDegreeInfo();

        const auto cell_degree = cell_infos.cell_degree();
        const auto grad_degree = cell_infos.grad_degree();

        const auto cell_basis_size = vector_basis_size(cell_degree, dimension, dimension);

        size_t gb_size = 0;
        if (small_def)
        {
            gb_size = sym_matrix_basis_size(grad_degree, dimension, dimension);
        }
        else
        {
            gb_size = matrix_basis_size(grad_degree, dimension, dimension);
        }

        const auto grad_basis_size = gb_size;

        const auto num_faces_dofs = vector_faces_dofs(msh, faces_infos);
        const auto num_total_dofs = cell_basis_size + num_faces_dofs;
        const auto grad_dim_dofs  = num_grad_dim_dofs(small_def);

        matrix_type AT = matrix_type::Zero(grad_basis_size, grad_basis_size);
        vector_type aT = vector_type::Zero(grad_basis_size);

        RTF   = vector_type::Zero(num_total_dofs);
        F_int = vector_type::Zero(num_total_dofs);

        assert(RkT.cols() == uTF.rows());
        assert(RkT.rows() == grad_basis_size);

        //   std::cout << "sol" << std::endl;
        //   std::cout << uTF.transpose() << std::endl;

        const vector_type RkT_uTF = RkT * uTF;

        //  std::cout << "RkT: " << RkT.norm() << std::endl;
        //  std::cout << "RkT_Utf: " << RkT_uTF.transpose() << std::endl;

        const auto gb  = make_matrix_monomial_basis(msh, cl, grad_degree);
        const auto gbs = make_sym_matrix_monomial_basis(msh, cl, grad_degree);

        eigen_compatible_stdvector<static_matrix_type> gphi;

        const auto cell_id = msh.lookup(cl);
        const auto nb_qp   = behavior.numberOfQP(cell_id);
        for (int i_qp = 0; i_qp < nb_qp; i_qp++)
        {
            //  std::cout << "qp: " << qp.point() << std::endl;
            // Compute gradient basis function
            const auto qp = behavior.quadrature_point(cell_id, i_qp);
            if (small_def)
            {
                gphi = gbs.eval_functions(qp.point());
            }
            else
            {
                gphi = gb.eval_functions(qp.point());
            }

            assert(gphi.size() == grad_basis_size);

            // Compute local gradient and norm
            const auto RkT_iqn = eval(RkT_uTF, gphi);

            // std::cout << "RkT_iqn" << std::endl;
            // std::cout << RkT_iqn << std::endl;

            // Compute behavior
            tc.tic();
            const auto [stress, Cep] = compute_behavior(behavior, cell_id, i_qp, RkT_iqn, small_def);
            // std::cout << "stress" << std::endl;
            // std::cout << stress << std::endl;
            // std::cout << "Cep" << std::endl;
            // std::cout << Cep << std::endl;
            tc.toc();
            time_law += tc.to_double();

            // Compute rigidity
            this->compute_rigidity(Cep, gphi, small_def, qp.weight(), grad_dim_dofs, grad_basis_size, AT);
            // Compute internal force
            this->compute_internal_forces(stress, gphi, small_def, qp.weight(), grad_dim_dofs, grad_basis_size, aT);
            // Compute contact terms
            tc.tic();
            this->compute_contact_terms();
            tc.toc();
            time_contact += tc.to_double();
        }

        // compute external forces
        this->compute_external_forces(msh, cl, load, cell_degree, RTF);

        // Symmetrize rigidity matrix
        this->symmetrized_rigidity_matrix(grad_basis_size, AT, small_def);

        // std::cout << "AT: " << AT.norm() << std::endl;
        // std::cout << AT << std::endl;
        // std::cout << "aT: " << aT.norm() << std::endl;

        K_int = RkT.transpose() * AT * RkT;
        F_int = RkT.transpose() * aT;
        RTF -= F_int;

        // std::cout << "K: " << K_int.norm() << std::endl;
        // std::cout << K_int << std::endl;
        // std::cout << "F: " << F_int.norm() << std::endl;
        // std::cout << "RTF: " << RTF.norm() << std::endl;

        assert(K_int.rows() == num_total_dofs);
        assert(K_int.cols() == num_total_dofs);
        assert(RTF.rows() == num_total_dofs);
    }
};
}

} // end namespace diskpp