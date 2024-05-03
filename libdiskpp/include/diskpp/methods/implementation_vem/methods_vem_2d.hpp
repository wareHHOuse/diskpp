/* This file is copyright of the following authors:
 * Karol Cascavita (C) 2024         karol.cascavita@polito.it
 *
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 */

#pragma once

#include "diskpp/bases/bases.hpp"

namespace disk
{

namespace vem_2d
{

template<typename T>
using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T>
using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;


/**
 * @brief Compute stiffness matrix \f$ (\nabla m_i,\nabla m_j) i,j = 1,...,dim( \mathcal{P}^{k}(T)) \f$ 
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl  cell T
 * @param degree  vem degree k 
 * @return dynamic_matrix<typename Mesh::coordinate_type> stiffness matrix in \f$ \mathcal{P}^{k}(T)\f$ 
 */

template<typename Mesh>
matrix_type<typename Mesh::coordinate_type>
stiffness_matrix_polk(const Mesh&    msh,
         const typename Mesh::cell_type&  cl,
         const size_t   degree)
{
    //Matrix G in papers
    using T = typename Mesh::coordinate_type;
    const auto num_faces = howmany_faces(msh, cl);

    const auto cbs = scalar_basis_size(degree, Mesh::dimension);
    matrix_type<T> stiffness = make_stiffness_matrix(msh, cl, cbs);
    

    const auto cb  = make_scalar_monomial_basis(msh, cl, degree);
    const auto qps = integrate(msh, cl, 2 * (degree - 1));
    for (auto& qp : qps)
    {
        const auto phi    = cb.eval_functions(qp.point());
        const auto qp_phi = priv::inner_product(qp.weight(), phi);
        stiffness.row(0) += priv::outer_product(qp_phi, phi);
    }
    stiffness.row(0) /= measure(msh, cl); 

    return stiffness;
}



template<typename Mesh>
matrix_type<typename Mesh::coordinate_type>

matrix_BF(const Mesh&    msh,
         const typename Mesh::cell_type&  cl,
         const size_t   degree)
{
    using T = typename Mesh::coordinate_type;

    const auto num_faces = howmany_faces(msh, cl);
    const auto cbs = scalar_basis_size(degree, Mesh::dimension);

    matrix_type<T> BF = matrix_type<T>::Zero(cbs, num_faces * degree);

    return BF;
}


template<typename Mesh>
matrix_type<typename Mesh::coordinate_type>
matrix_BT(const Mesh&    msh,
         const typename Mesh::cell_type&  cl,
         const size_t   degree)
{
    using T = typename Mesh::coordinate_type;

    const auto cbs = scalar_basis_size(degree, Mesh::dimension);
    const auto lbs = scalar_basis_size(degree-2, Mesh::dimension);

    const auto num_faces = howmany_faces(msh, cl);
    const auto num_dofs = num_faces * degree + lbs; 

    const auto cb  = make_scalar_monomial_basis(msh, cl, degree);
    const auto lcb = make_scalar_monomial_basis(msh, cl, degree-2);

    matrix_type<T> BT = matrix_type<T>::Zero(cbs, lbs);

    size_t ipol = 2;

    for (size_t k = 2; k <= degree; k++)
    {
        for (size_t i = 0; i <= k; i++) 
        {
            ipol++;

            size_t powx = k-i;
            size_t powy = i;

            size_t dxx_row = k-2;
            size_t dyy_row = dxx_row;
            size_t dxx_col = i;
            size_t dyy_col = i-2;

            if(dxx_row >= i && dxx_col >= 0)
            {
                size_t ipolxx = 0.5 * (dxx_row + 1) * dxx_row + dxx_col;
                size_t coeff_xx = powx * (powx - 1);
                BT(ipol, ipolxx) = coeff_xx;
            }

            if(dyy_row >= 0 && dyy_col >= 0)
            {
                size_t ipolyy = 0.5 * (dyy_row + 1) * dyy_row + dyy_col;          
                size_t coeff_yy = powy * (powy - 1);
                BT(ipol, ipolyy) = coeff_yy;    
            }
        }
    }

    BT *= measure(msh, cl); 
    return BT;
}


} // end namespace vem 2d
}// end namespace diskpp