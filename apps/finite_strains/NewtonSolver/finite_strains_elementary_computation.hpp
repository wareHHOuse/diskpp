/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
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

#include <cassert>

#include "common/eigen.hpp"
#include "mechanics/behaviors/laws/behaviorlaws.hpp"
#include "mechanics/deformation_tensors.hpp"
#include "revolution/bases"
#include "revolution/methods/hho"
#include "revolution/quadratures"

#include "timecounter.h"

namespace NLE
{

template<typename T>
struct MaterialParameters
{
    T lambda;
    T mu;
    T H;
    T K;
    T sigma_y0;

    T
    converttomu(const T E, const T nu)
    {
        return E / (2 * (1 + nu));
    }

    T
    converttolambda(const T E, const T nu)
    {
        return E * nu / ((1 + nu) * (1 - 2 * nu));
    }

    T
    converttoH(const T E, const T ET, const T K)
    {
        return E * ET / (E - ET) - 1.5 * K;
    }
};

template<typename MeshType>
class finite_strains
{
    typedef MeshType                             mesh_type;
    typedef typename mesh_type::scalar_type      scalar_type;
    typedef typename mesh_type::cell             cell_type;
    typedef typename revolution::hho_degree_info hdi_type;

    const static int dimension = mesh_type::dimension;

    typedef static_matrix<scalar_type, dimension, dimension> gvt;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    const mesh_type& m_msh;
    const hdi_type&  m_hdi;

    bool two_dim;

    int
    num_dofs_dim() const
    {
        if (two_dim)
            return 4;
        else
            return 9;
    }

    template<int DIM>
    eigen_compatible_stdvector<static_matrix<scalar_type, DIM, DIM>>
    compute_A_gphi(const static_tensor<scalar_type, DIM>&                                  tens,
                   const eigen_compatible_stdvector<static_matrix<scalar_type, DIM, DIM>>& gphi) const
    {
        const int grad_basis_size = gphi.size();
        const int DIM2            = DIM * DIM;

        eigen_compatible_stdvector<static_matrix<scalar_type, DIM, DIM>> Aphi;
        Aphi.reserve(grad_basis_size);

        // poly classique
        for (int i = 0; i < grad_basis_size; i += DIM2)
        {
            int row = i;
            for (int k = 0; k < DIM; k++)
            { // depend de l'ordre des bases
                for (int l = 0; l < DIM; l++)
                { // depend de l'ordre des bases
                    Aphi.push_back(disk::tm_prod(tens, gphi[row], l, k));
                    row++;
                }
            }
        }

        return Aphi;
    }

  public:
    matrix_type K_int;
    vector_type RTF;
    vector_type F_int;
    double      time_law;

    finite_strains(const mesh_type& msh, const hdi_type& hdi) : m_msh(msh), m_hdi(hdi)
    {
        if (dimension == 2)
            two_dim = true;
        else if (dimension == 3)
            two_dim = false;
        else
            assert(false);
    }

    template<typename Function, typename LawCell, typename LawData>
    void
    compute(const cell_type&   cl,
            const Function&    load,
            const matrix_type& GT,
            const vector_type& uTF,
            LawCell&           law,
            const LawData&     material_data,
            bool               elatic_modulus)
    {
        const auto cell_degree = m_hdi.cell_degree();
        const auto grad_degree = m_hdi.grad_degree();
        const auto face_degree = m_hdi.face_degree();

        const auto cell_basis_size = revolution::vector_basis_size(cell_degree, dimension, dimension);
        const auto grad_basis_size = revolution::matrix_basis_size(grad_degree, dimension, dimension);
        const auto face_basis_size = revolution::vector_basis_size(face_degree, dimension - 1, dimension);

        time_law = 0.0;
        timecounter tc;

        const auto fcs            = faces(m_msh, cl);
        const auto num_faces      = fcs.size();
        const auto num_total_dofs = cell_basis_size + num_faces * face_basis_size;
        const auto dim_dofs       = num_dofs_dim();

        matrix_type AT = matrix_type::Zero(grad_basis_size, grad_basis_size);
        vector_type aT = vector_type::Zero(grad_basis_size);

        RTF   = vector_type::Zero(num_total_dofs);
        F_int = vector_type::Zero(num_total_dofs);

        assert(GT.cols() == uTF.rows());
        assert(GT.rows() == grad_basis_size);

        // std::cout << "sol" << std::endl;
        // std::cout << uTF.transpose() << std::endl;

        const vector_type GT_uTF = GT * uTF;

        // std::cout << "ET: " << GsT.norm() << std::endl;
        // std::cout << GsT << std::endl;

        auto& law_quadpoints = law.getQPs();

        auto gb = revolution::make_matrix_monomial_basis(m_msh, cl, grad_degree);

        //std::cout << "nb: " << law_quadpoints.size() << std::endl;
        for (auto& qp : law_quadpoints)
        {
            //std::cout << "qp: " << qp.point() << std::endl;
            const auto gphi = gb.eval_functions(qp.point());

            assert(gphi.size() == grad_basis_size);

            // Compute local gradient and norm
            //std::cout << "GT_utf: " << GsT_uTF << std::endl;
            const auto GT_iqn = revolution::eval(GT_uTF, gphi);
            const gvt  incr_F = GT_iqn - qp.getTotalStrainPrev();
            // std::cout << "Em" << std::endl;
            // std::cout << qp.getTotalStrainPrev() << std::endl;
            // std::cout << "dE" << std::endl;
            // std::cout << incr_strain << std::endl;

            // Compute bahavior
            tc.tic();
            const auto tensor_behavior = qp.compute_whole(incr_F, material_data, !elatic_modulus);
            tc.toc();
            time_law += tc.to_double();

            //std::cout << "module " << tensor_behavior.second << std::endl;
            const auto qp_A_gphi = compute_A_gphi(qp.weight() * tensor_behavior.second, gphi);

            for (int j = 0; j < grad_basis_size; j += dim_dofs)
            {
                int col = j;
                for (int k = 0; k < dimension; k++)
                { // depend de l'ordre des bases
                    for (int l = 0; l < dimension; l++)
                    { // depend de l'ordre des bases
                        for (int i = col; i < grad_basis_size; i++)
                        {
                            AT(i, col) += qp_A_gphi[i](l, k) * gphi[col](l, k);
                        }
                        col++;
                    }
                }
            }

            // compute (PK1(u), G^k_T v)_T
            const auto stress_qp = revolution::priv::inner_product(qp.weight(), tensor_behavior.first);
            // std::cout << "stress" << std::endl;
            // std::cout << tensor_behavior.first << std::endl;

            for (int i = 0; i < grad_basis_size; i += dim_dofs)
            {
                int row = i;
                for (int k = 0; k < dimension; k++)
                { // depend de l'ordre des bases
                    for (int l = 0; l < dimension; l++)
                    { // depend de l'ordre des bases
                        // compute (PK1(u), G^k_T v)_T
                        aT(row) += stress_qp(l, k) * gphi[row](l, k);
                        row++;
                    }
                }
            }
        }

        // compute (f,v)_T
        auto cb                         = revolution::make_vector_monomial_basis(m_msh, cl, cell_degree);
        RTF.segment(0, cell_basis_size) = make_rhs(m_msh, cl, cb, load);

        // lower part AT
        for (int i = 0; i < grad_basis_size; i++)
            for (int j = i; j < grad_basis_size; j++)
                AT(i, j) = AT(j, i);

        // std::cout << "AT: " << AT.norm() << std::endl;
        // std::cout << AT << std::endl;
        // std::cout << "at: " << aT.norm() << std::endl;

        K_int = GT.transpose() * AT * GT;
        F_int = GT.transpose() * aT;
        RTF -= F_int;

        // std::cout << "K: " << K_int.norm() << std::endl;
        // std::cout << K_int << std::endl;
        // std::cout << "F: " << F_int.norm() << std::endl;

        assert(K_int.rows() == num_total_dofs);
        assert(K_int.cols() == num_total_dofs);
        assert(RTF.rows() == num_total_dofs);
    }
};

} // end namespace NLE