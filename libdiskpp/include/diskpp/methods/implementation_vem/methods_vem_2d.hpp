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

#include <Eigen/Dense>
#include "diskpp/bases/bases.hpp"
#include "diskpp/quadratures/bits/quad_raw_gauss_lobatto.hpp"
#include <cassert>

namespace disk
{

/* In DiSk++ 2D faces (= edges) are stored in lexicographically ordered
 * pairs (n0, n1). Node n0 is always the one with the smallest number, n1
 * the one with the highest number. This is true both globally (inside the
 * Mesh data structure) and locally (the data returned by `faces(msh, cl)`).
 * In addition, for compatibility with the DGA, only generic_mesh stores
 * point ids clockwise.
 * To simplify the implementation of 2D VEM, the following functions reorder
 * the faces so that:
 *  1) they appear CCW in the list returned by faces_ccw(msh, cl)
 *  2) they have a flag attached that indicates whether they should be flipped
 */
template<typename T>
std::array< std::pair< typename simplicial_mesh<T,2>::face_type, bool >, 3>
faces_ccw(const simplicial_mesh<T,2>& msh, const typename simplicial_mesh<T,2>::cell_type& cl)
{
    using face_type = typename simplicial_mesh<T,2>::face_type;
    auto fcs = faces(msh, cl);
    auto ptids = cl.point_ids();
    std::array< std::pair<face_type, bool>, 3> reorder;

    reorder[0] = { fcs[0], false };
    reorder[1] = { fcs[1], false };
    reorder[2] = { fcs[2], true };

    return reorder;
}

template<typename T>
std::array< std::pair< typename cartesian_mesh<T,2>::face_type, bool >, 4>
faces_ccw(const cartesian_mesh<T,2>& msh, const typename cartesian_mesh<T,2>::cell_type& cl)
{
    using face_type = typename cartesian_mesh<T,2>::face_type;
    auto fcs = faces(msh, cl);
    auto ptids = cl.point_ids();
    std::array< std::pair<face_type, bool>, 4> reorder;

    reorder[0] = { fcs[0], false };
    reorder[1] = { fcs[2], false };
    reorder[2] = { fcs[3], true };
    reorder[3] = { fcs[1], true };

    return reorder;
}


template<typename T>
std::vector< std::pair< typename generic_mesh<T,2>::face_type, bool > >
faces_ccw(const generic_mesh<T,2>& msh, const typename generic_mesh<T,2>::cell_type& cl)
{
    /* This could be potentially expensive on polygons with
     * many faces because it does n^2 comparisons. It could
     * be improved taking into account the fact that faces()
     * return value is ordered lexicographically. We'll address
     * this when/if really needed.
     */
    using face_type = typename generic_mesh<T,2>::face_type;
    auto fcs = faces(msh, cl);
    auto ptids = cl.point_ids();
    std::vector< std::pair<face_type, bool> > reorder;
    for (size_t i = 0; i < ptids.size(); i++)
    {
        bool flipped = false;
        auto n0 = ptids[i];
        auto n1 = ptids[(i+1)%ptids.size()];
        if (n0 > n1) {
            std::swap(n0, n1);
            flipped = true;
        }
        
        for (size_t j = 0; j < fcs.size(); j++) {
            auto fc_ptids = fcs[j].point_ids();
            assert(fc_ptids.size() == 2);
            assert(fc_ptids[0] < fc_ptids[1]);
            if (fc_ptids[0] == n0 and fc_ptids[1] == n1)
                reorder.push_back({fcs[j], flipped});
        }
    }

    assert(fcs.size() == reorder.size());
    return reorder;
}

template<disk::mesh_2D Mesh>
class dof_mapper {

    std::vector<std::vector<size_t>>    l2g_byface;
    std::vector<std::vector<size_t>>    l2g_bycell;

public:
    dof_mapper(const Mesh& msh, size_t k) {
        recompute(msh, k);
    }
    
    void recompute(const Mesh& msh, size_t k) {
        if (k < 1)
            throw std::invalid_argument("dof_mapper: degree must be > 0");

        size_t pts_base = msh.points_size();

        /* For each edge in the mesh (edges are ordered lexicographically)
         * compute the offset of all its edge-based dofs */
        for (auto& fc : faces(msh)) {
            auto ptids = fc.point_ids();
            std::vector<size_t> l2g(k+1);
            l2g[0] = ptids[0];
            for (size_t i = 1; i < k; i++)
                l2g[i] = pts_base++;
            l2g[k] = ptids[1];
            l2g_byface.push_back( std::move(l2g) );
        }

        for (auto& cl : msh)
        {
            auto fcs_ccw = faces_ccw(msh, cl);

            /* OK, now we build the whole element l2g */
            std::vector<size_t> cl_l2g;
            for (size_t i = 0; i < fcs_ccw.size(); i++) {
                auto [fc, flip] = fcs_ccw[i];
                auto& l2g = l2g_byface[ offset(msh, fc) ];
                /* The last node of this segment is the first of
                 * the next, so we discard it. */
                if (flip) {
                    cl_l2g.insert(cl_l2g.end(), l2g.rbegin(), l2g.rend()-1);
                }
                else {
                    cl_l2g.insert(cl_l2g.end(), l2g.begin(), l2g.end()-1);
                }
            }
            l2g_bycell.push_back( std::move(cl_l2g) );
        }
    }

    auto face_to_global(size_t face_num) {
        assert(face_num < l2g_byface.size());
        return l2g_byface[face_num];
    }

    auto cell_to_global(size_t cell_num) {
        assert(cell_num < l2g_bycell.size());
        return l2g_bycell[cell_num];
    }
};


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

template<disk::mesh_2D Mesh>
matrix_type<typename Mesh::coordinate_type>
matrix_G(const Mesh&    msh,
         const typename Mesh::cell_type&  cl,
         const size_t   degree)
{
    using T = typename Mesh::coordinate_type;
    const auto num_faces = howmany_faces(msh, cl);

    const auto cbs = scalar_basis_size(degree, Mesh::dimension);
    const auto cb  = make_scalar_monomial_basis(msh, cl, degree);
    matrix_type<T> stiffness = make_stiffness_matrix(msh, cl, cb);
    
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



template<disk::mesh_2D Mesh>
matrix_type<typename Mesh::coordinate_type>
matrix_BF(const Mesh&    msh,
         const typename Mesh::cell_type&  cl,
         const size_t   degree)
{
    assert(degree > 0);

    using T = typename Mesh::coordinate_type;

    const auto num_faces = howmany_faces(msh, cl);
    const auto cb  = make_scalar_monomial_basis(msh, cl, degree);
    const auto cbs = cb.size(); 
    const auto num_dofs_bnd = num_faces * degree;

    matrix_type<T> BF = matrix_type<T>::Zero(cbs, num_dofs_bnd);

    auto fcs_ccw = faces_ccw(msh, cl);

/*
    std::cout << "Original : " << std::endl;
    for(const auto&  fc: fcs)
    {
        std::cout << " * face(" <<offset(msh,fc)<<") : [";
        
        for (auto pid : fc.point_ids())        
            std::cout << pid << ", ";
        std::cout <<"]" <<std::endl;
    }

    std::cout << "CCW : "<< fcs_ccw.size() << std::endl;

    for(const auto&  [fc, flip] : fcs_ccw)
    {
        std::cout << " * face(" <<offset(msh,fc)<< "," << flip <<") : [";
        
        for (auto pid : fc.point_ids())        
            std::cout << pid << ", ";
        std::cout <<"]" <<std::endl;
    }

    std::cout << " ---------------------------------"<<std::endl;
*/

    auto iface = 0;
    for(const auto [fc, flip] : fcs_ccw)
    {
        const auto n   = normal(msh, cl, fc);
        const auto pts = points(msh, fc);

        size_t idx0(0), idx1(1); 
        if(flip)
            std::swap(idx0,idx1);

        auto qps = disk::quadrature::gauss_lobatto(2 * degree -1, pts[idx0], pts[idx1]);

        auto qcount = iface * degree; 
        //std::cout << " * face(" <<offset(msh,fc)<<") :  with qcount = " << qcount << std::endl;

        for (auto& qp : qps)
        {       
            const auto dphi   = cb.eval_gradients(qp.point());
            const auto dphi_n = qp.weight() * (dphi * n);

            size_t index_dof = qcount%(num_dofs_bnd); 
            BF.block(0, index_dof, cbs, 1) += dphi_n;

            //std::cout << "   - index_dof = " << index_dof << std::endl;

            qcount++;
        }  
        iface++;
    }

    return BF;
}


template<disk::mesh_2D Mesh>
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

    int ipol = 2;

    auto h = diameter(msh, cl); 
    auto h_powk = h;

    for (int k = 2; k <= degree; k++)
    {
        h_powk *= h;

        for (int i = 0; i <= k; i++) 
        {
            ipol++;

            int powx = k-i;
            int powy = i;

            int coeff_xx = powx * (powx - 1);
            int coeff_yy = powy * (powy - 1);

            if(coeff_xx > 0)
            {
                int dxx_row = k-2;
                int dxx_col = i;

                int ipolxx = 0.5 * (dxx_row + 1) * dxx_row + dxx_col;
                BT(ipol, ipolxx) = -coeff_xx / h_powk ;
            }

            if(coeff_yy > 0) 
            {
                int dyy_row = k-2;
                int dyy_col = i-2;

                int ipolyy = 0.5 * (dyy_row + 1) * dyy_row + dyy_col;          
                BT(ipol, ipolyy) = -coeff_yy / h_powk;    
            }
        }
    }

    BT *= measure(msh, cl); 

    auto P0v = 1;
    BT(0, 0) = P0v; 

    return BT;
}


template<disk::mesh_2D Mesh>
matrix_type<typename Mesh::coordinate_type>
matrix_B(const Mesh&    msh,
         const typename Mesh::cell_type&  cl,
         const size_t   degree)
{
    using T = typename Mesh::coordinate_type;

    const auto cbs = scalar_basis_size(degree, Mesh::dimension);
    const auto lbs = scalar_basis_size(degree-2, Mesh::dimension);

    const auto num_faces = howmany_faces(msh, cl);
    const auto num_dofs = num_faces * degree + lbs; 

    matrix_type<T> B = matrix_type<T>::Zero(cbs, num_dofs);

    B.block(0, 0, cbs, num_faces * degree)  = matrix_BF(msh, cl, degree);
    B.block(0,num_faces * degree, cbs, lbs) = matrix_BT(msh, cl, degree);

    return B;
}


template<disk::mesh_2D Mesh>
matrix_type<typename Mesh::coordinate_type>
matrix_D(const Mesh&    msh,
         const typename Mesh::cell_type&  cl,
         const size_t   degree)
{
    using T = typename Mesh::coordinate_type;

    const auto cbs = scalar_basis_size(degree, Mesh::dimension);
    const auto lbs = scalar_basis_size(degree-2, Mesh::dimension);

    const auto cb  = make_scalar_monomial_basis(msh, cl, degree);
    const auto lb  = make_scalar_monomial_basis(msh, cl, degree-2);

    const auto num_faces = howmany_faces(msh, cl);
    const auto num_dofs = num_faces * degree + lbs; 

    matrix_type<T> D = matrix_type<T>::Zero(num_dofs, cbs);

    auto count = 0;
    for(const auto& [fc, flip] : faces_ccw(msh, cl))
    {
        const auto pts = points(msh, fc);

        size_t idx0(0), idx1(1); 
        if(flip)
            std::swap(idx0,idx1);

        auto face_qps = disk::quadrature::gauss_lobatto(2 * degree -1, pts[idx0], pts[idx1]);

        for (size_t iq = 0; iq < face_qps.size() - 1;  iq++)
        {  
            auto qp = face_qps[iq];     
            const auto phi   = cb.eval_functions(qp.point());

            D.block(count, 0, 1, cbs) = phi.transpose();
            count++;
        }  
    }

    auto area = measure(msh, cl);
    const auto cell_qps = integrate(msh, cl, 2*degree-2);

    auto offset_cell_dofs = num_faces * degree;

    for (const auto& qp : cell_qps)
    {
        const auto phi     = cb.eval_functions(qp.point());
        const auto phi_dof = lb.eval_functions(qp.point());

        D.block(offset_cell_dofs, 0, lbs, cbs)  += qp.weight() * phi_dof * phi.transpose();
    }
    D.block(offset_cell_dofs, 0, lbs, cbs) /= area; 

    return D;
}



} // end namespace vem 2d
}// end namespace diskpp