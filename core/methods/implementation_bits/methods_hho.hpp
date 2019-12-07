 /*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018, 2019
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018, 2019   matteo.cicuttin@enpc.fr
 * Karol Cascavita (C) 2018                     karol.cascavita@enpc.fr
 * Nicolas Pignet  (C) 2018, 2019               nicolas.pignet@enpc.fr
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
#include "boundary_conditions/boundary_conditions.hpp"

using namespace Eigen;

namespace disk
{

/**
 * @brief this class contains the different polynomial degrees which can used with HHO methods
 *
 */
class hho_degree_info
{
    /**
     * @brief
     *  cell_deg : polynomial degree used for the cells
     *  face_deg : polynomial degree used for the faces
     *  grad_deg : polynomial degree used for the gradient
     *
     */
    size_t  cell_deg, face_deg, grad_deg;

public:
    /**
     * @brief Construct a new hho degree info object
     * The default polynomial degree is 1 for the face and cell degrees
     *
     */   
    hho_degree_info()
        : cell_deg(1), face_deg(1),  grad_deg(1)
    {}

    /**
     * @brief Construct a new hho degree info object
     *
     * @param degree polynomial degree used for the cells, faces and gradient
     */
    explicit hho_degree_info(size_t degree)
        : cell_deg(degree), face_deg(degree),  grad_deg(degree)
    {}

    /**
     * @brief Construct a new hho degree info object
     *
     * @param cd polynomial degree used for the cells
     * @param fd polynomial degree used for the faces
     *
     * @brief Note that
     * fd >= 0 and cd >=0
     * fd - 1 <= cd <= fd +1
     */
    hho_degree_info(size_t cd, size_t fd)
    {
        bool c1 = fd > 0  && (cd == fd-1 || cd == fd || cd == fd+1);
        bool c2 = fd == 0 && (cd == fd || cd == fd+1);
        if ( c1 || c2 )
        {
            cell_deg            = cd;
            face_deg            = fd;
            grad_deg            = fd;
        }
        else
        {
            std::cout << "Invalid cell degree. Reverting to equal-order" << std::endl;
            cell_deg            = fd;
            face_deg            = fd;
            grad_deg            = fd;
        }
    }

    /**
     * @brief Construct a new hho degree info object
     *
     * @param cd polynomial degree used for the cells
     * @param fd polynomial degree used for the faces
     * @param gd polynomial degree used for the gradient
     *
     * @brief Note that
     * fd >= 0, cd >=0 and gd >=0
     * fd - 1 <= cd <= fd +1
     * gd >= fd
     */
    hho_degree_info(size_t cd, size_t fd, size_t gd)
    {
       bool c1 = fd > 0 && (cd == fd - 1 || cd == fd || cd == fd + 1);
       bool c2 = fd == 0 && (cd == fd || cd == fd + 1);
       bool c3 = gd >= fd;
       if (c1 || c2 || c3) {
          cell_deg           = cd;
          face_deg           = fd;
          grad_deg           = gd;
       } else {
          std::cout << "Invalid cell degree. Reverting to equal-order" << std::endl;
          cell_deg           = fd;
          face_deg           = fd;
          grad_deg           = fd;
       }
    }

    /**
     * @brief Return the polynomial degree used for the cells
     *
     * @return size_t polynomial degree used for the cells
     */
    size_t cell_degree() const
    {
        return cell_deg;
    }

    /**
     * @brief Return the polynomial degree used for the faces
     *
     * @return size_t polynomial degree used for the faces
     */
    size_t face_degree() const
    {
        return face_deg;
    }

    /**
     * @brief Return the polynomial degree used for the reconstruction operator
     *
     * @return size_t polynomial degree used for the reconstruction operator
     */
    size_t reconstruction_degree() const
    {
        return face_deg+1;
    }

    /**
     * @brief Return the polynomial degree used for the gradient reconstruction
     *
     * @return size_t polynomial degree used for the gradient reconstruction
     */
    size_t
    grad_degree() const
    {
       return grad_deg;
    }

    /**
     * @brief Print the different degree
     *
     */
    void
    info_degree() const
    {
       std::cout << "cell degree: " << cell_deg << std::endl;
       std::cout << "face degree: " << face_deg << std::endl;
       std::cout << "gradient degree: " << grad_deg << std::endl;
       std::cout << "reconstruction degree: " << face_deg+1 << std::endl;
    }
};

/**
 * @brief Compute the interpolation operator \f$ \hat{I}(u) := (\Pi^{cd}_T(u), \Pi^{fd}_{\partial T}(u)) \f$
 * of a function \f$ u \in H^1(T;\mathbb{R}) \f$ where \f$ \Pi^{cd}_T \f$, respc. \f$ \Pi^{fd}_{\partial T} \f$ are the L2-orthogonal projector on the cell and faces
 *
 * @tparam Mesh type of the mesh
 * @param msh  mesh
 * @param cl  cell
 * @param hdi  hho degree informations
 * @param u  function  \f$ u \in H^1(T;\mathbb{R}) \f$ to project
 * @param di additonal quadratue degree
 * @return dynamic_vector<typename Mesh::coordinate_type> coefficient of the the interpolation operator \f$ \hat{I}(u) \f$
 */
template<typename Mesh>
dynamic_vector<typename Mesh::coordinate_type>
project_function(const Mesh&                      msh,
                 const typename Mesh::cell_type&  cl,
                 const hho_degree_info&           hdi,
                 const scalar_rhs_function<Mesh>& u,
                 size_t                           di = 0)
{
    typedef dynamic_vector<typename Mesh::coordinate_type> vector_type;

    const auto cbs       = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    const auto fbs       = scalar_basis_size(hdi.face_degree(), Mesh::dimension - 1);
    const auto num_faces = howmany_faces(msh, cl);

    vector_type ret = vector_type::Zero(cbs + num_faces * fbs);

    ret.block(0, 0, cbs, 1) = project_function(msh, cl, hdi.cell_degree(), u, di);

    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < num_faces; i++)
    {
        ret.segment(cbs + i * fbs, fbs) = project_function(msh, fcs[i], hdi.face_degree(), u, di);
    }

    return ret;
}

/**
 * @brief Compute the interpolation operator \f$ \hat{I}(u) := (\Pi^{cd}_T(u), \Pi^{fd}_{\partial T}(u)) \f$
 * of a function \f$ u \in H^1(T;\mathbb{R}^d) \f$ where \f$ \Pi^{cd}_T \f$, respc. \f$ \Pi^{fd}_{\partial T} \f$ are the L2-orthogonal projector on the cell and faces
 *
 * @tparam Mesh type of the mesh
 * @param msh  mesh
 * @param cl  cell
 * @param hdi  hho degree informations
 * @param u  function  \f$ u \in H^1(T;\mathbb{R}^d) \f$ to project
 * @param di additonal quadratue degree
 * @return dynamic_vector<typename Mesh::coordinate_type> coefficient of the the interpolation operator \f$ \hat{I}(u) \f$
 */
template<typename Mesh>
dynamic_vector<typename Mesh::coordinate_type>
project_function(const Mesh&                      msh,
                 const typename Mesh::cell_type&  cl,
                 const hho_degree_info&           hdi,
                 const vector_rhs_function<Mesh>& u,
                 size_t                           di = 0)
{
    typedef dynamic_vector<typename Mesh::coordinate_type> vector_type;

    const auto cbs       = vector_basis_size(hdi.cell_degree(), Mesh::dimension, Mesh::dimension);
    const auto fbs       = vector_basis_size(hdi.face_degree(), Mesh::dimension - 1, Mesh::dimension);
    const auto num_faces = howmany_faces(msh, cl);

    vector_type ret = vector_type::Zero(cbs + num_faces * fbs);

    ret.block(0, 0, cbs, 1) = project_function(msh, cl, hdi.cell_degree(), u, di);

    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < num_faces; i++)
    {
        ret.segment(cbs + i * fbs, fbs) = project_function(msh, fcs[i], hdi.face_degree(), u, di);
    }

    return ret;
}


/**
 * @brief Compute the term \f$ (\nabla R^{k+1}_T(\hat{u}), \nabla R^{k+1}_T(\hat{v}) )_T   \f$  where \f$ R^{k+1}_T \f$ is the high-order
 * reconstruction operator which if solution of
 *
 * \f$ (\nabla R^{k+1}_T(\hat{u}), \nabla w )_T = (\nabla u_T, \nabla w )_T + (u_{\partial T} - u_T, \nabla w n_T)_{\partial T}, \quad \forall w \f$
 *
 *  The term \f$ (\nabla R^{k+1}_T(\hat{u}), \nabla R^{k+1}_T(\hat{v}) )_T   \f$ is the HHO conterpart of \f$ (\nabla u, \nabla v )_T   \f$
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl  cell
 * @param di  hho degree information
 * @return std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
dynamic_matrix<typename Mesh::coordinate_type>  > the first term is \f$ \nabla R^{k+1}_T \f$ and the second term is \f$ (\nabla R^{k+1}_T(\hat{u}), \nabla R^{k+1}_T(\hat{v}) )_T   \f$
 */
template<typename Mesh>
std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
             dynamic_matrix<typename Mesh::coordinate_type>  >
make_scalar_hho_laplacian(const Mesh& msh, const typename Mesh::cell_type& cl,
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

    const matrix_type stiff  = make_stiffness_matrix(msh, cl, cb);
    matrix_type gr_lhs = matrix_type::Zero(rbs-1, rbs-1);
    matrix_type gr_rhs = matrix_type::Zero(rbs-1, cbs + num_faces*fbs);

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

namespace priv {
template<typename Mesh, typename GradBasis>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>,
          dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_gradrec_impl(const Mesh&                     msh,
                             const typename Mesh::cell_type& cl,
                             const hho_degree_info&          di,
                             const GradBasis&                gb)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    const auto celdeg  = di.cell_degree();
    const auto facdeg  = di.face_degree();
    const auto graddeg = gb.degree();

    const auto cb = make_scalar_monomial_basis(msh, cl, celdeg);

    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    const auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);
    const auto gbs = gb.size();

    const auto num_faces = howmany_faces(msh, cl);

    const matrix_type gr_lhs = make_mass_matrix(msh, cl, gb);
    matrix_type       gr_rhs = matrix_type::Zero(gbs, cbs + num_faces * fbs);

    // (vT, div(tau))_T
    if (graddeg > 0)
    {
        const auto qps = integrate(msh, cl, celdeg + graddeg - 1);
        for (auto& qp : qps)
        {
            const auto c_phi = cb.eval_functions(qp.point());
            const auto g_dphi  = gb.eval_divergences(qp.point());
            const vector_type qp_g_dphi = qp.weight() * g_dphi;

            gr_rhs.block(0, 0, gbs, cbs) -= priv::outer_product(qp_g_dphi, c_phi);
        }
    }

    // (vF, tau.n)_F
    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = normal(msh, cl, fc);
        const auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

        const auto qps_f = integrate(msh, fc, graddeg + facdeg);
        for (auto& qp : qps_f)
        {
            const vector_type f_phi      = fb.eval_functions(qp.point());
            const auto        g_phi      = gb.eval_functions(qp.point());
            const vector_type qp_g_phi_n = g_phi * (qp.weight() * n);

            gr_rhs.block(0, cbs + i * fbs, gbs, fbs) += priv::outer_product(qp_g_phi_n, f_phi);
        }
    }

    matrix_type oper = gr_lhs.ldlt().solve(gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}
}


/**
 * @brief Compute the term \f$ (G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$  where \f$ G^{k}_T \f$ is the hho gradient reconstruction which if solution of
 *
 * \f$ (G^{k}_T(\hat{u})(\hat{u}), \tau )_T = (\nabla u_T, \tau )_T + (u_{\partial T} - u_T, \tau n_T)_{\partial T}, \quad \forall w \f$
 *
 *  The term \f$ (G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$ is the HHO conterpart of \f$ (\nabla u, \nabla v )_T   \f$
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl  cell
 * @param di  hho degree information
 * @return std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
dynamic_matrix<typename Mesh::coordinate_type>  > the first term is \f$ G^{k}_T \f$ and the second term is \f$ (G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$
 */
template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>,
          dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_gradrec(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    const auto graddeg = di.grad_degree();
    const auto gb      = make_vector_monomial_basis(msh, cl, graddeg);

    return priv::make_vector_hho_gradrec_impl(msh, cl, di, gb);
}

/**
 * @brief Compute the term \f$ (K G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$  where \f$ G^{k}_T \f$ is the hho gradient reconstruction which if solution of
 *
 * \f$ (G^{k}_T(\hat{u})(\hat{u}), \tau )_T = (\nabla u_T, \tau )_T + (u_{\partial T} - u_T, \tau n_T)_{\partial T}, \quad \forall w \f$
 *
 *  The term \f$ (K G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$ is the HHO conterpart of \f$ (K \nabla u, \nabla v )_T   \f$
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl  cell
 * @param di  hho degree information
 * @param mfield is a material fiel which can depend on the coordiantes (mfield is symmetric positive definite matrix)
 * @return std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
dynamic_matrix<typename Mesh::coordinate_type>  > the first term is \f$ G^{k}_T \f$ and the second term is \f$ (K G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$
 */
template<typename Mesh, typename TensorField>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>,
          dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_gradrec(const Mesh&                     msh,
                        const typename Mesh::cell_type& cl,
                        const hho_degree_info&          di,
                        const TensorField&              mfield)
{
    const auto gradrec = make_vector_hho_gradrec(msh, cl, di);

    const auto graddeg = di.grad_degree();
    const auto gb      = make_vector_monomial_basis(msh, cl, graddeg);

    const auto mass = make_mass_matrix(msh, cl, gb, mfield);
    const auto LHS  = gradrec.first.transpose() * mass * gradrec.first;

    return std::make_pair(gradrec.first, LHS);
}

/**
 * @brief Compute the term \f$ (G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$  where \f$ G^{k}_T \f$ is the hho gradient reconstruction which if solution of
 *
 * \f$ (G^{k}_T(\hat{u})(\hat{u}), \tau )_T = (\nabla u_T, \tau )_T + (u_{\partial T} - u_T, \tau n_T)_{\partial T}, \quad \forall w \f$
 *
 *  The term \f$ (G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$ is the HHO conterpart of \f$ (K \nabla u, \nabla v )_T   \f$
 *
 *  Here the reconstructed gradient \f$ G^{k}_T \f$ is reconstructed in the Raviart-Thomas space (works only for simplical cell)
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl  simplicial cell (triangle or tetrahedron)
 * @param di  hho degree information
 * @return std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
dynamic_matrix<typename Mesh::coordinate_type>  > the first term is \f$ G^{k}_T \f$ and the second term is \f$ ( G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$
 */

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>,
          dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_gradrec_RT(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    const auto graddeg = di.grad_degree();
    const auto gb      = make_vector_monomial_basis_RT(msh, cl, graddeg);

    return make_vector_hho_gradrec_impl(msh, cl, di, gb);
}

/**
 * @brief Compute the term \f$ (K G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$  where \f$ G^{k}_T \f$ is the hho gradient reconstruction which if solution of
 *
 * \f$ (G^{k}_T(\hat{u})(\hat{u}), \tau )_T = (\nabla u_T, \tau )_T + (u_{\partial T} - u_T, \tau n_T)_{\partial T}, \quad \forall w \f$
 *
 *  The term \f$ (K G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$ is the HHO conterpart of \f$ (K \nabla u, \nabla v )_T   \f$
 *
 *  Here the reconstructed gradient \f$ G^{k}_T \f$ is reconstructed in the Raviart-Thomas space (works only for simplical cell)
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl  simplicial cell (triangle or tetrahedron)
 * @param di  hho degree information
 * @param mfield is a material fiel which can depend on the coordiantes (mfield is symmetric positive definite matrix)
 * @return std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
dynamic_matrix<typename Mesh::coordinate_type>  > the first term is \f$ G^{k}_T \f$ and the second term is \f$ (K G^{k}_T(\hat{u}), G^{k}_T(\hat{v}) )_T   \f$
 */
template<typename Mesh, typename TensorField>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>,
          dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_gradrec_RT(const Mesh&                     msh,
                           const typename Mesh::cell_type& cl,
                           const hho_degree_info&          di,
                           const TensorField&              mfield)
{
    const auto gradrec_RT = make_vector_hho_gradrec_RT(msh, cl, di);

    const auto graddeg = di.grad_degree();
    const auto gb      = make_vector_monomial_basis_RT(msh, cl, graddeg);

    const auto mass = make_mass_matrix(msh, cl, gb, mfield);
    const auto LHS  = gradrec_RT.first.transpose() * mass * gradrec_RT.first;

    return std::make_pair(gradrec_RT.first, LHS);
}

namespace priv
{
size_t
nb_lag(const size_t dim)
{
    size_t lag = 1;
    if (dim == 3)
        lag = 3;
    return lag;
}

template<typename Mesh, typename BasisCell>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>,
          dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_symmetric_laplacian_impl_facezero(const Mesh&                     msh,
                                                        const typename Mesh::cell_type& cl,
                                                        const BasisCell&                cb)
{
    using T        = typename Mesh::coordinate_type;
    const size_t N = Mesh::dimension;

    const auto recdeg = 1;
    const auto celdeg = cb.degree();
    const auto facdeg = 0;

    const auto rb = make_vector_monomial_basis(msh, cl, recdeg);

    const auto rbs = vector_basis_size(recdeg, N, N);
    const auto cbs = cb.size();
    const auto fbs = vector_basis_size(facdeg, N - 1, N);

    const auto num_faces = howmany_faces(msh, cl);

    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, N, N>             gradient_type;
    typedef Matrix<T, Dynamic, N>       function_type;

    const size_t rbs_ho         = rbs - N;
    const size_t num_total_dofs = cbs + num_faces * fbs;
    const size_t nb_lag         = priv::nb_lag(N);

    matrix_type stiff  = matrix_type::Zero(rbs, rbs);
    matrix_type gr_lhs = matrix_type::Zero(rbs_ho + nb_lag, rbs_ho + nb_lag);
    matrix_type gr_rhs = matrix_type::Zero(rbs_ho + nb_lag, num_total_dofs);

    auto qps = integrate(msh, cl, 0);
    for (auto& qp : qps)
    {
        const auto dphi = rb.eval_sgradients(qp.point());
        const auto qp_dphi = priv::inner_product(qp.weight(), dphi);
        stiff += priv::outer_product(qp_dphi, dphi);
    }

    gr_lhs.block(0, 0, rbs_ho, rbs_ho) = stiff.block(N, N, rbs_ho, rbs_ho);

    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = normal(msh, cl, fc);
        const auto fb = make_vector_monomial_basis(msh, fc, facdeg);

        const auto qps_f = integrate(msh, fc, 0);
        for (auto& qp : qps_f)
        {
            eigen_compatible_stdvector<gradient_type> r_dphi_tmp = rb.eval_sgradients(qp.point());
            auto                                      begin_iter = std::next(r_dphi_tmp.begin(), N);

            eigen_compatible_stdvector<gradient_type> r_dphi(rbs_ho);
            std::copy(begin_iter, r_dphi_tmp.end(), r_dphi.begin());

            const function_type f_phi    = fb.eval_functions(qp.point());
            const function_type qp_r_dphi_n = priv::inner_product(r_dphi, priv::inner_product(qp.weight(), n));
            gr_rhs.block(0, cbs + i * fbs, rbs_ho, fbs) += priv::outer_product(qp_r_dphi_n, f_phi);
        }
    }

    qps = integrate(msh, cl, recdeg);

    matrix_type rot = matrix_type::Zero(rbs, nb_lag);
    for (auto& qp : qps)
    {
        const auto rphi = rb.eval_curls(qp.point());
        rot += qp.weight() * rphi;
    }
    gr_lhs.block(0, rbs_ho, rbs_ho, nb_lag) += rot.bottomLeftCorner(rbs_ho, nb_lag);
    gr_lhs.block(rbs_ho, 0, nb_lag, rbs_ho) += rot.bottomLeftCorner(rbs_ho, nb_lag).transpose();

    // use LU solver because lhs is only symmetric and positive
    matrix_type sol  = gr_lhs.lu().solve(gr_rhs);
    matrix_type oper = sol.block(0, 0, rbs_ho, num_total_dofs);
    matrix_type gr   = gr_rhs.block(0, 0, rbs_ho, num_total_dofs);
    matrix_type data = gr.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh, typename BasisCell>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>,
          dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_symmetric_laplacian_impl(const Mesh& msh, const typename Mesh::cell_type& cl, const BasisCell& cb, const hho_degree_info& di)
{
    using T        = typename Mesh::coordinate_type;
    const size_t N = Mesh::dimension;

    const auto recdeg = di.reconstruction_degree();
    const auto celdeg = cb.degree();
    const auto facdeg = di.face_degree();

    if (facdeg == 0)
    {
        return make_vector_hho_symmetric_laplacian_impl_facezero(msh, cl, cb);
    }

    const auto rb = make_vector_monomial_basis(msh, cl, recdeg);

    const auto rbs = vector_basis_size(recdeg, N, N);
    const auto cbs = cb.size();
    const auto fbs = vector_basis_size(facdeg, N - 1, N);

    const auto num_faces = howmany_faces(msh, cl);

    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, N, N>             gradient_type;
    typedef Matrix<T, Dynamic, N>       function_type;

    const size_t rbs_ho         = rbs - N;
    const size_t num_total_dofs = cbs + num_faces * fbs;
    const size_t nb_lag         = priv::nb_lag(N);

    matrix_type stiff  = matrix_type::Zero(rbs, rbs);
    matrix_type gr_lhs = matrix_type::Zero(rbs_ho + nb_lag, rbs_ho + nb_lag);
    matrix_type gr_rhs = matrix_type::Zero(rbs_ho + nb_lag, num_total_dofs);

    const auto qps = integrate(msh, cl, 2 * (recdeg - 1));
    for (auto& qp : qps)
    {
        const auto dphi = rb.eval_sgradients(qp.point());
        const auto qp_dphi = priv::inner_product(qp.weight(), dphi);
        stiff += priv::outer_product(qp_dphi, dphi);
    }

    gr_lhs.block(0, 0, rbs_ho, rbs_ho) = stiff.block(N, N, rbs_ho, rbs_ho);
    gr_rhs.block(0, 0, rbs_ho, cbs) = stiff.block(N, 0, rbs_ho, cbs);

    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = normal(msh, cl, fc);
        const auto fb = make_vector_monomial_basis(msh, fc, facdeg);

        const auto qps_f = integrate(msh, fc, std::max(facdeg, celdeg) + recdeg - 1);
        for (auto& qp : qps_f)
        {
            eigen_compatible_stdvector<gradient_type> r_dphi_tmp = rb.eval_sgradients(qp.point());

            auto  begin_iter = std::next(r_dphi_tmp.begin(), N);
            eigen_compatible_stdvector<gradient_type> r_dphi(rbs_ho);
            std::copy(begin_iter, r_dphi_tmp.end(), r_dphi.begin());

            const function_type c_phi    = cb.eval_functions(qp.point());
            const function_type f_phi    = fb.eval_functions(qp.point());
            const function_type qp_r_dphi_n = qp.weight() * priv::inner_product(r_dphi, n);
            gr_rhs.block(0, cbs + i * fbs, rbs_ho, fbs) += priv::outer_product(qp_r_dphi_n, f_phi);
            gr_rhs.block(0, 0, rbs_ho, cbs) -= priv::outer_product(qp_r_dphi_n, c_phi);
        }
    }

    const auto qps_2 = integrate(msh, cl, recdeg);

    matrix_type rot = matrix_type::Zero(rbs, nb_lag);
    for (auto& qp : qps_2)
    {
        const auto rphi = rb.eval_curls(qp.point());
        rot += qp.weight() * rphi;
    }
    gr_lhs.block(0, rbs_ho, rbs_ho, nb_lag) += rot.bottomLeftCorner(rbs_ho, nb_lag);
    gr_lhs.block(rbs_ho, 0, nb_lag, rbs_ho) += rot.bottomLeftCorner(rbs_ho, nb_lag).transpose();

    // use LU solver because lhs is only symmetric and positive
    matrix_type sol  = gr_lhs.lu().solve(gr_rhs);
    matrix_type oper = sol.block(0, 0, rbs_ho, num_total_dofs);
    matrix_type gr   = gr_rhs.block(0, 0, rbs_ho, num_total_dofs);
    matrix_type data = gr.transpose() * oper;

    return std::make_pair(oper, data);
}
}


template<typename Mesh>
std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
             dynamic_matrix<typename Mesh::coordinate_type>  >
make_vector_hho_symmetric_laplacian(const Mesh& msh,
                          const typename Mesh::cell_type& cl,
                          const hho_degree_info& di)
{
    const auto cb = make_vector_monomial_basis(msh, cl, di.cell_degree());

    return priv::make_vector_hho_symmetric_laplacian_impl(msh, cl, cb, di);
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>,
          dynamic_matrix<typename Mesh::coordinate_type>>
make_matrix_symmetric_gradrec(const Mesh&                     msh,
                              const typename Mesh::cell_type& cl,
                              const hho_degree_info&          di)
{
   using T        = typename Mesh::coordinate_type;
   typedef Matrix<T, Dynamic, Dynamic> matrix_type;

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
    if (celdeg > 0)
    {
        const auto qpc = integrate(msh, cl, graddeg + celdeg - 1);
        for (auto& qp : qpc)
        {
            const auto gphi    = gb.eval_functions(qp.point());
            const auto dphi    = cb.eval_sgradients(qp.point());
            const auto qp_dphi = priv::inner_product(qp.weight(), dphi);

            gr_rhs.block(0, 0, gbs, cbs) += priv::outer_product(gphi, qp_dphi);

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

            const auto qp_gphi_n = priv::inner_product(gphi, priv::inner_product(qp.weight(), n));
            gr_rhs.block(0, cbs + i * fbs, gbs, fbs) += priv::outer_product(qp_gphi_n, fphi);
            gr_rhs.block(0, 0, gbs, cbs) -= priv::outer_product(qp_gphi_n, cphi);
        }
    }

    matrix_type oper = gr_lhs.ldlt().solve(gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>,
          dynamic_matrix<typename Mesh::coordinate_type>>
make_matrix_symmetric_gradrec_RT(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const size_t N = Mesh::dimension;

    const auto graddeg = di.grad_degree();
    const auto celdeg  = di.cell_degree();
    const auto facdeg  = di.face_degree();

    const auto gb = make_sym_matrix_monomial_basis_RT(msh, cl, graddeg);
    const auto cb = make_vector_monomial_basis(msh, cl, celdeg);

    const auto gbs = sym_matrix_basis_size_RT(graddeg, Mesh::dimension, Mesh::dimension);
    const auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
    const auto fbs = vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension);

    const auto num_faces = howmany_faces(msh, cl);

    const matrix_type gr_lhs = make_mass_matrix(msh, cl, gb);
    matrix_type gr_rhs = matrix_type::Zero(gbs, cbs + num_faces * fbs);

    // compute rhs
    if (celdeg > 0)
    {
        const auto qpc = integrate(msh, cl, graddeg + celdeg - 1);
        for (auto& qp : qpc)
        {
            const auto gphi = gb.eval_functions(qp.point());
            const auto dphi = cb.eval_sgradients(qp.point());
            const auto qp_dphi = priv::inner_product(qp.weight(), dphi);

            gr_rhs.block(0, 0, gbs, cbs) += priv::outer_product(gphi, qp_dphi);

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

            const auto qp_gphi_n = priv::inner_product(gphi, priv::inner_product(qp.weight(), n));
            gr_rhs.block(0, cbs + i * fbs, gbs, fbs) += priv::outer_product(qp_gphi_n, fphi);
            gr_rhs.block(0, 0, gbs, cbs) -= priv::outer_product(qp_gphi_n, cphi);
        }
    }

    matrix_type oper = gr_lhs.ldlt().solve(gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
std::pair<   dynamic_matrix<typename Mesh::coordinate_type>,
             dynamic_matrix<typename Mesh::coordinate_type>  >
make_hho_divergence_reconstruction(const Mesh& msh, const typename Mesh::cell_type& cl,
                                   const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    auto cbas_s = make_scalar_monomial_basis(msh, cl, di.face_degree());

    const auto dr_lhs = make_mass_matrix(msh, cl, cbas_s);
    const auto dr_rhs = make_hho_divergence_reconstruction_rhs(msh, cl, di);

    matrix_type oper = dr_lhs.ldlt().solve(dr_rhs);
    matrix_type data = dr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_hho_divergence_reconstruction_rhs(const Mesh& msh, const typename Mesh::cell_type& cl,
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

            dr_rhs.block(0, 0, rbs, cbs) -= qp.weight() * priv::outer_product(s_dphi, v_phi);
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

            const auto qp_f_phi_n = priv::inner_product(f_phi, priv::inner_product(qp.weight(), n));
            dr_rhs.block(0, cbs + i * fbs, rbs, fbs) += priv::outer_product(s_phi, qp_f_phi_n);
        }
    }

    return dr_rhs;
}

/**
 * @brief compute the stabilization term \f$ \sum_{F \in F_T} 1/h_F(u_F-\Pi^k_F(u_T), v_F-\Pi^k_F(v_T))_F \f$
 * for scalar HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param di hho degree information
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_hdg_stabilization(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

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
        matrix_type mass  = make_mass_matrix(msh, fc, fb);
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

            trace += priv::outer_product(priv::inner_product(qp.weight(), f_phi), c_phi);
        }

        tr.block(0, cbs + i * fbs, fbs, fbs) = -mass;
        tr.block(0, 0, fbs, cbs) = trace;

        oper.block(0, 0, fbs, cbs) = mass.ldlt().solve(trace);
        data += oper.transpose() * tr * (1./h);
    }

    return data;
}

/**
 * @brief compute the stabilization term \f$ \sum_{F \in F_T} 1/h_F(u_F-u_T, v_F-v_T)_F \f$
 * for scalar HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param di hho degree information
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_dg_stabilization(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto celdeg = di.cell_degree();
    const auto facdeg = di.face_degree();

    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    const auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);

    const auto num_faces  = howmany_faces(msh, cl);
    const auto total_dofs = cbs + num_faces * fbs;

    matrix_type data = matrix_type::Zero(total_dofs, total_dofs);

    auto cb = make_scalar_monomial_basis(msh, cl, celdeg);

    const auto fcs = faces(msh, cl);

    for (size_t i = 0; i < num_faces; i++)
    {
        const auto fc = fcs[i];
        const auto hf = diameter(msh, fc);
        auto       fb = make_scalar_monomial_basis(msh, fc, facdeg);

        matrix_type mass_F = make_mass_matrix(msh, fc, fb);
        matrix_type mass_T = make_mass_matrix(msh, fc, cb);
        matrix_type trace  = matrix_type::Zero(fbs, cbs);

        const auto qps = integrate(msh, fc, facdeg + celdeg);
        for (auto& qp : qps)
        {
            const auto c_phi = cb.eval_functions(qp.point());
            const auto f_phi = fb.eval_functions(qp.point());

            trace += priv::outer_product(priv::inner_product(qp.weight(), f_phi), c_phi);
        }

        const auto offset = cbs + i * fbs;

        data.block(0, 0, cbs, cbs) += mass_T / hf;
        data.block(offset, offset, fbs, fbs) += mass_F / hf;
        data.block(0, offset, cbs, fbs) -= trace.transpose() / hf;
        data.block(offset, 0, fbs, cbs) -= trace / hf;
    }

    return data;
}

/**
 * @brief compute the stabilization term \f$\sum_{F \in F_T} 1/h_F(u_F - u_T + \Pi^k_T R^{k+1}_T(\hat{u}_T) - R^{k+1}_T(\hat{u}_T), v_F - v_T + \Pi^k_T R^{k+1}_T(\hat{v}_T) - R^{k+1}_T(\hat{v}_T))_F \f$
 * for scalar HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param reconstruction reconstruction operator \f$ R^{k+1}_T \f$
 * @param di hho degree information
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_scalar_hho_stabilization(const Mesh&                                                     msh,
                              const typename Mesh::cell_type&                                 cl,
                              const dynamic_matrix<typename Mesh::coordinate_type>& reconstruction,
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

    const auto cb = make_scalar_monomial_basis(msh, cl, recdeg);

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

        matrix_type face_mass_matrix  = make_mass_matrix(msh, fc, fb);
        matrix_type face_trace_matrix = matrix_type::Zero(fbs, rbs);

        const auto face_quadpoints = integrate(msh, fc, recdeg + facdeg);
        for (auto& qp : face_quadpoints)
        {
            const auto        f_phi   = fb.eval_functions(qp.point());
            const auto        c_phi   = cb.eval_functions(qp.point());
            face_trace_matrix += priv::outer_product(priv::inner_product(qp.weight(), f_phi), c_phi);
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

/**
 * @brief compute the stabilization term \f$\sum_{F \in F_T} 1/h_F(u_F - u_T + \Pi^k_T R^{k+1}_T(\hat{u}_T) - R^{k+1}_T(\hat{u}_T), v_F - v_T + \Pi^k_T R^{k+1}_T(\hat{v}_T) - R^{k+1}_T(\hat{v}_T))_F \f$
 * for vectorial HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param reconstruction reconstruction operator \f$ R^{k+1}_T \f$
 * @param di hho degree information
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_vector_hho_stabilization(const Mesh&                                                     msh,
                              const typename Mesh::cell_type&                                 cl,
                              const dynamic_matrix<typename Mesh::coordinate_type>& reconstruction,
                              const hho_degree_info&                                          di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto recdeg = di.reconstruction_degree();
    const auto celdeg = di.cell_degree();
    const auto facdeg = di.face_degree();

    const auto rbs = vector_basis_size(recdeg, Mesh::dimension, Mesh::dimension);
    const auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
    const auto fbs = vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension);

    const size_t N = Mesh::dimension;

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

        matrix_type face_mass_matrix  = make_mass_matrix(msh, fc, fb);
        matrix_type face_trace_matrix = matrix_type::Zero(fbs, rbs);

        const auto face_quadpoints = integrate(msh, fc, recdeg + facdeg);
        for (auto& qp : face_quadpoints)
        {
            const auto        f_phi   = fb.eval_functions(qp.point());
            const auto        c_phi   = cb.eval_functions(qp.point());
            face_trace_matrix += priv::outer_product(priv::inner_product(qp.weight(), f_phi), c_phi);
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

namespace priv
{
// static condensation
template<typename Mesh, typename T>
auto
static_condensation_impl(const Mesh&                                                      msh,
                         const typename Mesh::cell_type&                                  cl,
                         const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                         const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs,
                         const size_t                                                     num_cell_dofs,
                         const size_t                                                     num_face_dofs)
{
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    const auto fcs            = faces(msh, cl);
    const auto num_faces      = fcs.size();
    const auto num_faces_dofs = num_faces * num_face_dofs;
    const auto num_total_dofs = num_cell_dofs + num_faces_dofs;

    assert(lhs.rows() == lhs.cols());
    assert(lhs.cols() == num_total_dofs);

    if ((rhs.size() != num_cell_dofs) && (rhs.size() != num_total_dofs))
    {
        throw std::invalid_argument("static condensation: incorrect size of the rhs");
    }

    const matrix_type K_TT = lhs.topLeftCorner(num_cell_dofs, num_cell_dofs);
    const matrix_type K_TF = lhs.topRightCorner(num_cell_dofs, num_faces_dofs);
    const matrix_type K_FT = lhs.bottomLeftCorner(num_faces_dofs, num_cell_dofs);
    const matrix_type K_FF = lhs.bottomRightCorner(num_faces_dofs, num_faces_dofs);

    assert(K_TT.cols() == num_cell_dofs);
    assert(K_TT.cols() + K_TF.cols() == lhs.cols());
    assert(K_TT.rows() + K_FT.rows() == lhs.rows());
    assert(K_TF.rows() + K_FF.rows() == lhs.rows());
    assert(K_FT.cols() + K_FF.cols() == lhs.cols());

    const vector_type cell_rhs  = rhs.head(num_cell_dofs);
    vector_type       faces_rhs = vector_type::Zero(num_faces_dofs);

    if (rhs.size() == num_total_dofs)
    {
        faces_rhs = rhs.tail(num_faces_dofs);
    }

    const auto K_TT_ldlt = K_TT.ldlt();
    if (K_TT_ldlt.info() != Eigen::Success)
    {
        throw std::invalid_argument("static condensation: K_TT is not positive definite");
    }

    const matrix_type AL = K_TT_ldlt.solve(K_TF);
    const vector_type bL = K_TT_ldlt.solve(cell_rhs);

    const matrix_type AC = K_FF - K_FT * AL;
    const vector_type bC = faces_rhs - K_FT * bL;

    return std::make_tuple(std::make_pair(AC, bC), AL, bL);
}

// static decondensation for primal scalar problem
template<typename Mesh, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
static_decondensation_impl(const Mesh&                                                      msh,
                           const typename Mesh::cell_type&                                  cl,
                           const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                           const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs,
                           const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              solF,
                           const size_t                                                     num_cell_dofs,
                           const size_t                                                     num_face_dofs)
{
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    const auto fcs            = faces(msh, cl);
    const auto num_faces      = fcs.size();
    const auto num_faces_dofs = num_faces * num_face_dofs;
    const auto num_total_dofs = num_cell_dofs + num_faces_dofs;

    assert(lhs.rows() == lhs.cols());
    assert(lhs.cols() == num_total_dofs);

    if ((rhs.size() < num_cell_dofs))
    {
        throw std::invalid_argument("static condensation: incorrect size of the rhs");
    }

    const matrix_type K_TT = lhs.topLeftCorner(num_cell_dofs, num_cell_dofs);
    const matrix_type K_TF = lhs.topRightCorner(num_cell_dofs, num_faces_dofs);

    const vector_type solT = K_TT.ldlt().solve(rhs.head(num_cell_dofs) - K_TF * solF);

    vector_type ret          = vector_type::Zero(num_total_dofs);
    ret.head(num_cell_dofs)  = solT;
    ret.tail(num_faces_dofs) = solF;

    return ret;
}
} // namespace priv

// static condensation for primal scalar problem like diffusion
template<typename Mesh, typename T>
auto
make_scalar_static_condensation_withMatrix(const Mesh&                                                 msh,
                                      const typename Mesh::cell_type&                                  cl,
                                      const hho_degree_info&                                           hdi,
                                      const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                                      const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs)
{
    const auto num_cell_dofs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    const auto num_face_dofs = scalar_basis_size(hdi.face_degree(), Mesh::dimension - 1);

    return priv::static_condensation_impl(msh, cl, lhs, rhs, num_cell_dofs, num_face_dofs);
}

template<typename Mesh, typename T>
auto
make_scalar_static_condensation(const Mesh&                                                 msh,
                           const typename Mesh::cell_type&                                  cl,
                           const hho_degree_info&                                           hdi,
                           const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                           const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs)
{
    return std::get<0>(make_scalar_static_condensation_withMatrix(msh, cl, hdi, lhs, rhs));
}

// static decondensation for primal scalar problem
template<typename Mesh, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
make_scalar_static_decondensation(const Mesh&                                                 msh,
                             const typename Mesh::cell_type&                                  cl,
                             const hho_degree_info                                            hdi,
                             const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                             const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs,
                             const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              solF)
{
    const auto num_cell_dofs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    const auto num_face_dofs = scalar_basis_size(hdi.face_degree(), Mesh::dimension - 1);

    return static_decondensation_impl(msh, cl, lhs, rhs, solF, num_cell_dofs, num_face_dofs);
}

// static decondensation for primal scalar problem
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
make_scalar_static_decondensation_withMatrix(const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& AL,
                                             const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              bL,
                                             const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              solF)
{
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    vector_type ret          = vector_type::Zero(bL.size() + solF.size());
    ret.head(bL.size())      = bL - AL * solF;
    ret.tail(solF.size())    = solF;

    return ret;
}

// static condensation for primal vectorial problem like elasticity
template<typename Mesh, typename T>
auto
make_vector_static_condensation_withMatrix(const Mesh&                                                 msh,
                                      const typename Mesh::cell_type&                                  cl,
                                      const hho_degree_info&                                           hdi,
                                      const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                                      const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs)
{
    const auto num_cell_dofs = vector_basis_size(hdi.cell_degree(), Mesh::dimension, Mesh::dimension);
    const auto num_face_dofs = vector_basis_size(hdi.face_degree(), Mesh::dimension - 1, Mesh::dimension);

    return priv::static_condensation_impl(msh, cl, lhs, rhs, num_cell_dofs, num_face_dofs);
}

template<typename Mesh, typename T>
auto
make_vector_static_condensation(const Mesh&                                                 msh,
                           const typename Mesh::cell_type&                                  cl,
                           const hho_degree_info&                                           hdi,
                           const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                           const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs)
{
    return std::get<0>(make_vector_static_condensation_withMatrix(msh, cl, hdi, lhs, rhs));
}

// static decondensation for primal vector problem
template<typename Mesh, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
make_vector_static_decondensation(const Mesh&                                                 msh,
                             const typename Mesh::cell_type&                                  cl,
                             const hho_degree_info                                            hdi,
                             const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                             const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs,
                             const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              solF)
{
    const auto num_cell_dofs = vector_basis_size(hdi.cell_degree(), Mesh::dimension, Mesh::dimension);
    const auto num_face_dofs = vector_basis_size(hdi.face_degree(), Mesh::dimension - 1, Mesh::dimension);

    return static_decondensation_impl(msh, cl, lhs, rhs, solF, num_cell_dofs, num_face_dofs);
}

// static decondensation for primal vector problem
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
make_vector_static_decondensation_withMatrix(const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& AL,
                                             const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              bL,
                                             const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              solF)
{
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    vector_type ret          = vector_type::Zero(bL.size() + solF.size());
    ret.head(bL.size())      = bL - AL * solF;
    ret.tail(solF.size())    = solF;

    return ret;
}


template<typename Mesh>
auto
make_hho_stokes(const Mesh& msh, const typename Mesh::cell_type& cl,
                const hho_degree_info& hdi, const bool& use_sym_grad)
{
    if(use_sym_grad)
        return make_vector_hho_symmetric_laplacian(msh, cl, hdi);
    else
        return make_vector_hho_laplacian(msh, cl, hdi);
}

template<typename Mesh>
auto
make_hlow_stokes(const Mesh& msh, const typename Mesh::cell_type& cl,
                const hho_degree_info& hdi, const bool& use_sym_grad)
{
    if(use_sym_grad)
        return make_matrix_symmetric_gradrec(msh, cl, hdi);
    else
        return make_marix_hho_gradrec(msh, cl, hdi);
}

// define some optimization
namespace priv
{



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
std::pair<dynamic_matrix<typename Mesh::coordinate_type>,
          dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_laplacian(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& hdi)
{
    const auto hho_scalar_laplacian = make_scalar_hho_laplacian(msh, cl, hdi);

    return make_vector_hho_laplacian(msh, cl, hdi, hho_scalar_laplacian);
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>,
          dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_laplacian(const Mesh&                     msh,
                          const typename Mesh::cell_type& cl,
                          const hho_degree_info&          hdi,
                          const std::pair < dynamic_matrix<typename Mesh::coordinate_type>,
                          dynamic_matrix<typename Mesh::coordinate_type> >& hho_scalar_laplacian)
{
    const auto oper = priv::compute_grad_vector(msh, cl, hdi, hho_scalar_laplacian.first);
    const auto data = priv::compute_lhs_vector(msh, cl, hdi, hho_scalar_laplacian.second);

    return std::make_pair(oper, data);
}

/**
 * @brief compute the stabilization term \f$ \sum_{F \in F_T} 1/h_F(u_F-\Pi^k_F(u_T), v_F-\Pi^k_F(v_T))_F \f$
 * for vectorial HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param di hho degree information
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_vector_hdg_stabilization(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    const auto hdg_scalar_stab = make_scalar_hdg_stabilization(msh, cl, di);

    return priv::compute_lhs_vector(msh, cl, di, hdg_scalar_stab);
}

/**
 * @brief compute the stabilization term \f$ \sum_{F \in F_T} 1/h_F(u_F-u_T, v_F-v_T)_F \f$
 * for vectorial HHO unknowns
 *
 * @tparam Mesh type of the mesh
 * @param msh mesh
 * @param cl cell
 * @param di hho degree information
 * @return dynamic_matrix<typename Mesh::coordinate_type> return the stabilization term
 */
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_vector_dg_stabilization(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    const auto dg_scalar_stab = make_scalar_dg_stabilization(msh, cl, di);

    return priv::compute_lhs_vector(msh, cl, di, dg_scalar_stab);
}

// doesn't work for symmetric gradient
template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_vector_hho_stabilization_optim(const Mesh& msh, const typename Mesh::cell_type& cl,
                              const dynamic_matrix<typename Mesh::coordinate_type>& reconstruction_scalar,
                              const hho_degree_info& hdi)
{
    const auto hho_scalar_stab = make_scalar_hho_stabilization(msh, cl, reconstruction_scalar, hdi);

    return priv::compute_lhs_vector(msh, cl, hdi, hho_scalar_stab);
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>,
          dynamic_matrix<typename Mesh::coordinate_type>>
make_marix_hho_gradrec(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& hdi)

{
    const auto hho_gradrec_vector = make_vector_hho_gradrec(msh, cl, hdi);
    const auto lhs                = priv::compute_lhs_vector(msh, cl, hdi, hho_gradrec_vector.second);
    const auto oper               = priv::compute_grad_matrix(msh, cl, hdi, hho_gradrec_vector.first);

    return std::make_pair(oper, lhs);
}


} // end revolution
