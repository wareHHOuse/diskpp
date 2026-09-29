#pragma once

#include "diskpp/methods/hho"
#include "diskpp/bases/bases.hpp"
#include "diskpp/boundary_conditions/boundary_conditions.hpp"
#include "diskpp/adaptivity/adaptivity.hpp"

namespace disk
{


/* ==========================================================================
 * Primitives
 * ========================================================================*/

/* Traction operator at a quadrature point: row j = sigma(phi_j) . n
 *
 *   sigma(v) = 2 mu E(v) + lambda D(v) I,    D(v) = trace(E(v))
 * with E the discrete symmetric gradient (grad_op = gradrec.first).       */
template<typename Mesh, typename Basis>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
eval_traction_operator(const Mesh&                                                     msh,
                       const Basis&                                                    gb,
                       const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>& grad_op,
                       const typename Mesh::point_type&                                pt,
                       const static_vector<typename Mesh::coordinate_type, Mesh::dimension>& n,
                       const typename Mesh::coordinate_type                            mu,
                       const typename Mesh::coordinate_type                            lam)
{
    using T           = typename Mesh::coordinate_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;

    const size_t N   = Mesh::dimension;
    const size_t gbs = grad_op.rows();

    auto grad_basis = gb.eval_functions(pt);

    matrix_type traction_on_grad_basis(gbs, N);
    for (size_t i = 0; i < gbs; i++)
    {
        const auto& E_i = grad_basis[i];
        traction_on_grad_basis.row(i) =
            (2.0*mu*(E_i*n) + lam*E_i.trace()*n).transpose();
    }

    return grad_op.transpose() * traction_on_grad_basis;   // (ndofs x N)
}

/* Trace operator at a quadrature point on face F: row j = phi^F_j(pt).
 * Face-based variant     */
template<typename Mesh, typename Basis>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
eval_trace_operator(const Mesh&                      msh,
                    const Basis&                     fb,
                    const typename Mesh::point_type& pt,
                    const size_t                     num_total_dofs,
                    const size_t                     face_offset)
{
    using T           = typename Mesh::coordinate_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;

    const size_t N   = Mesh::dimension;
    const size_t fbs = fb.size();

    matrix_type trace_op = matrix_type::Zero(num_total_dofs, N);

    trace_op.block(face_offset, 0, fbs, N) = fb.eval_functions(pt);

    return trace_op;
}

/* Trace operator at a quadrature point on face F using the CELL basis:
 * row j = phi^T_j(pt), written into the CELL block (offset 0).
 *
 * This is the cell-version (Cascavita-Chouly-Ern) counterpart of
 * eval_trace_operator.  In the Nitsche penalty P_n(u) = sigma_n(u) - gamma u_n
 * the displacement u_n is taken here as the trace of the CELL unknown u_T|_F
 * instead of the face unknown u_F.  Requires the cell degree to exceed the
 * face degree (mixed order l = k+1), which spatial_operator already sets via
 * di(m, k+1, k).  The face-version helper eval_trace_operator is kept intact
 * so the two versions can be compared by swapping the call site.            */
template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
eval_cell_trace_operator(const Mesh&                      msh,
                         const typename Mesh::cell_type&  cl,
                         const size_t                     cell_degree,
                         const typename Mesh::point_type& pt,
                         const size_t                     num_total_dofs)
{
    using T           = typename Mesh::coordinate_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;

    const size_t N   = Mesh::dimension;
    auto         cb  = make_vector_monomial_basis(msh, cl, cell_degree);
    const size_t cbs = cb.size();

    matrix_type trace_op = matrix_type::Zero(num_total_dofs, N);

    trace_op.block(0, 0, cbs, N) = cb.eval_functions(pt);   // CELL block

    return trace_op;
}

/* ==========================================================================
 * Nitsche contact terms
 *
 * Sign convention: terms are returned WITHOUT their leading sign.
 * ========================================================================*/

/* Bilinear Nitsche consistency:
 *   (theta/gamma_n) \int sigma_n(u) sigma_n(v)
 * + (theta/gamma_t) \int sigma_t(u).sigma_t(v)  */
template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_vector_hho_nitsche(const Mesh&                     msh,
                        const typename Mesh::cell_type& cl,
                        const MeshDegreeInfo<Mesh>&     degree_infos,
                        const Matrix<typename Mesh::coordinate_type,
                                     Dynamic, Dynamic>& grad_op,
                        const typename Mesh::coordinate_type gamma_n_0,
                        const typename Mesh::coordinate_type gamma_t_0,
                        const typename Mesh::coordinate_type theta,
                        const typename Mesh::coordinate_type mu,
                        const typename Mesh::coordinate_type lam,
                        const vector_boundary_conditions<Mesh>& bnd)
{
    using T           = typename Mesh::coordinate_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t N = Mesh::dimension;

    const auto cell_infos = degree_infos.cellDegreeInfo(msh, cl);
    const auto grad_deg   = cell_infos.grad_degree();
    const auto quad_deg   = 2 * cell_infos.reconstruction_degree();
    const auto cbs        = vector_basis_size(cell_infos.cell_degree(), N, N);
    const auto fcs        = faces(msh, cl);
    const auto num_total_dofs = cbs + vector_faces_dofs(msh, cl, degree_infos);

    auto gb = make_sym_matrix_monomial_basis(msh, cl, grad_deg);

    matrix_type nitsche_matrix =
        matrix_type::Zero(num_total_dofs, num_total_dofs);

    for (auto& fc : fcs)
    {
        if (!bnd.is_contact_face(fc))
            continue;

        const auto n   = normal(msh, cl, fc);
        const auto qps = integrate(msh, fc, quad_deg);

        const T h_F       = diameter(msh, fc);
        const T gamma_n_h = gamma_n_0 / h_F;
        const T gamma_t_h = gamma_t_0 / h_F;

        for (auto& qp : qps)
        {
            matrix_type traction_op = eval_traction_operator(msh, gb, grad_op,
                                                             qp.point(), n, mu, lam);

            vector_type sigma_n_op = traction_op * n;                 // (ndofs)
            matrix_type sigma_t_op = traction_op
                                     - sigma_n_op * n.transpose();    // (ndofs x N)

            nitsche_matrix += qp.weight() * theta *
                ( (1.0/gamma_n_h) * sigma_n_op * sigma_n_op.transpose()
                + (1.0/gamma_t_h) * sigma_t_op * sigma_t_op.transpose() );
        }
    }

    return nitsche_matrix;    // theta and the penalties are already folded in
}

/* Normal contact term, evaluated at u_prev:
 *   (1/gamma_n) \int [P_n(u_prev)]_-  P_n^theta(v)
 * with  P_n(u)       = sigma_n(u)       - gamma_n u_dT.n
 *       P_n^theta(v) = theta sigma_n(v) - gamma_n v_dT.n
 *
 * bilateral = true enforces u_n = 0: the negative part is dropped, the
 * term is always active and becomes linear (standard Nitsche-Dirichlet on
 * the normal component).                                                  */
template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
make_vector_hho_contact_rhs(const Mesh&                     msh,
                            const typename Mesh::cell_type& cl,
                            const MeshDegreeInfo<Mesh>&     degree_infos,
                            const Matrix<typename Mesh::coordinate_type,
                                         Dynamic, Dynamic>& grad_op,
                            const typename Mesh::coordinate_type gamma_n_0,
                            const typename Mesh::coordinate_type theta,
                            const typename Mesh::coordinate_type mu,
                            const typename Mesh::coordinate_type lam,
                            const vector_boundary_conditions<Mesh>& bnd,
                            const Matrix<typename Mesh::coordinate_type,
                                         Dynamic, 1>& u_prev,
                            const bool bilateral = false)
{
    using T           = typename Mesh::coordinate_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t N = Mesh::dimension;

    const auto cell_infos = degree_infos.cellDegreeInfo(msh, cl);
    const auto grad_deg   = cell_infos.grad_degree();
    const auto quad_deg   = 2 * cell_infos.reconstruction_degree();
    const auto cbs        = vector_basis_size(cell_infos.cell_degree(), N, N);
    const auto fcs        = faces(msh, cl);
    const auto num_total_dofs = cbs + vector_faces_dofs(msh, cl, degree_infos);

    auto gb = make_sym_matrix_monomial_basis(msh, cl, grad_deg);

    vector_type contact_rhs = vector_type::Zero(num_total_dofs);

    size_t face_offset = cbs;
    for (auto& fc : fcs)
    {
        const auto   fdi = degree_infos.degreeInfo(msh, fc);
        const size_t nfd = vector_face_dofs(msh, fdi);

        if (bnd.is_contact_face(fc))
        {

            const auto n   = normal(msh, cl, fc);
            const auto qps = integrate(msh, fc, quad_deg);

            const T gamma_n_h = gamma_n_0 / diameter(msh, fc);

            for (auto& qp : qps)
            {
                matrix_type traction_op = eval_traction_operator(msh, gb, grad_op,
                                                                 qp.point(), n,
                                                                 mu, lam);
                matrix_type trace_op    = eval_cell_trace_operator(msh, cl, cell_infos.cell_degree(),
                                                              qp.point(), num_total_dofs);

                vector_type sigma_n_op = traction_op * n;
                vector_type trace_n_op = trace_op * n;

                vector_type P_u_op = sigma_n_op       - gamma_n_h * trace_n_op;
                vector_type P_v_op = theta*sigma_n_op - gamma_n_h * trace_n_op;

                const T P_u_value = P_u_op.dot(u_prev);
                const T P_u_eff   = bilateral ? P_u_value
                                              : std::min(P_u_value, T(0));

                contact_rhs += (qp.weight() / gamma_n_h) * P_u_eff * P_v_op;
            }
        }

        face_offset += nfd;
    }

    return contact_rhs;
}

/* Tresca friction term, evaluated at u_prev with a GIVEN constant
 * threshold s:
 *   (1/gamma_t) \int [P_t(u_prev)]_s . P_t^theta(v)
 * with  P_t(u) = sigma_t(u) - gamma_t u_dT,t   (tangential, vector)
 * and   [q]_s = q          if |q| <= s   (stick)
 *             = s q / |q|  otherwise      (slip)                          */
template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
make_vector_hho_tresca_rhs(const Mesh&                     msh,
                           const typename Mesh::cell_type& cl,
                           const MeshDegreeInfo<Mesh>&     degree_infos,
                           const Matrix<typename Mesh::coordinate_type,
                                        Dynamic, Dynamic>& grad_op,
                           const typename Mesh::coordinate_type gamma_t_0,
                           const typename Mesh::coordinate_type theta,
                           const typename Mesh::coordinate_type mu,
                           const typename Mesh::coordinate_type lam,
                           const typename Mesh::coordinate_type s_threshold,
                           const vector_boundary_conditions<Mesh>& bnd,
                           const Matrix<typename Mesh::coordinate_type,
                                        Dynamic, 1>& u_prev)
{
    using T           = typename Mesh::coordinate_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;
    using tensor_vec  = static_vector<T, Mesh::dimension>;

    const size_t N = Mesh::dimension;

    const auto cell_infos = degree_infos.cellDegreeInfo(msh, cl);
    const auto grad_deg   = cell_infos.grad_degree();
    const auto quad_deg   = 2 * cell_infos.reconstruction_degree();
    const auto cbs        = vector_basis_size(cell_infos.cell_degree(), N, N);
    const auto fcs        = faces(msh, cl);
    const auto num_total_dofs = cbs + vector_faces_dofs(msh, cl, degree_infos);

    auto gb = make_sym_matrix_monomial_basis(msh, cl, grad_deg);

    vector_type tresca_rhs = vector_type::Zero(num_total_dofs);

    size_t face_offset = cbs;
    for (auto& fc : fcs)
    {
        const auto   fdi = degree_infos.degreeInfo(msh, fc);
        const size_t nfd = vector_face_dofs(msh, fdi);

        if (bnd.is_contact_face(fc))
        {

            const auto n   = normal(msh, cl, fc);
            const auto qps = integrate(msh, fc, quad_deg);

            const T gamma_t_h = gamma_t_0 / diameter(msh, fc);

            for (auto& qp : qps)
            {
                matrix_type traction_op = eval_traction_operator(msh, gb, grad_op,
                                                                 qp.point(), n,
                                                                 mu, lam);
                matrix_type trace_op    = eval_cell_trace_operator(msh, cl, cell_infos.cell_degree(),
                                                              qp.point(), num_total_dofs);

                matrix_type sigma_t_op = traction_op - (traction_op * n) * n.transpose();
                matrix_type trace_t_op = trace_op    - (trace_op * n)    * n.transpose();

                matrix_type P_t_u_op = sigma_t_op       - gamma_t_h * trace_t_op;
                matrix_type P_t_v_op = theta*sigma_t_op - gamma_t_h * trace_t_op;

                tensor_vec P_t_u_value = P_t_u_op.transpose() * u_prev;

                const T P_t_u_norm = P_t_u_value.norm();

                tensor_vec P_t_u_proj;
                if (P_t_u_norm <= s_threshold)
                    P_t_u_proj = P_t_u_value;                                // stick
                else
                    P_t_u_proj = (s_threshold / P_t_u_norm) * P_t_u_value;   // slip

                tresca_rhs += (qp.weight() / gamma_t_h) * (P_t_v_op * P_t_u_proj);
            }
        }

        face_offset += nfd;
    }

    return tresca_rhs;
}

/* Normal contact Jacobian, linearised at u_lin:
 *   (1/gamma_n) \int H(P_n(u)) P_n(du) P_n^theta(v)
 * with H(x) = 1 if x <= 0 (active contact), 0 otherwise.
 * bilateral = true: always active.                                        */
template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_vector_hho_contact_jacobian(const Mesh&                     msh,
                                 const typename Mesh::cell_type& cl,
                                 const MeshDegreeInfo<Mesh>&     degree_infos,
                                 const Matrix<typename Mesh::coordinate_type,
                                              Dynamic, Dynamic>& grad_op,
                                 const typename Mesh::coordinate_type gamma_n_0,
                                 const typename Mesh::coordinate_type theta,
                                 const typename Mesh::coordinate_type mu,
                                 const typename Mesh::coordinate_type lam,
                                 const vector_boundary_conditions<Mesh>& bnd,
                                 const Matrix<typename Mesh::coordinate_type,
                                              Dynamic, 1>& u_lin,
                                 const bool bilateral = false)
{
    using T           = typename Mesh::coordinate_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t N = Mesh::dimension;

    const auto cell_infos = degree_infos.cellDegreeInfo(msh, cl);
    const auto grad_deg   = cell_infos.grad_degree();
    const auto quad_deg   = 2 * cell_infos.reconstruction_degree();
    const auto cbs        = vector_basis_size(cell_infos.cell_degree(), N, N);
    const auto fcs        = faces(msh, cl);
    const auto num_total_dofs = cbs + vector_faces_dofs(msh, cl, degree_infos);

    auto gb = make_sym_matrix_monomial_basis(msh, cl, grad_deg);

    matrix_type contact_jacobian =
        matrix_type::Zero(num_total_dofs, num_total_dofs);

    size_t face_offset = cbs;
    for (auto& fc : fcs)
    {
        const auto   fdi = degree_infos.degreeInfo(msh, fc);
        const size_t nfd = vector_face_dofs(msh, fdi);

        if (bnd.is_contact_face(fc))
        {

            const auto n   = normal(msh, cl, fc);
            const auto qps = integrate(msh, fc, quad_deg);

            const T gamma_n_h = gamma_n_0 / diameter(msh, fc);

            for (auto& qp : qps)
            {
                matrix_type traction_op = eval_traction_operator(msh, gb, grad_op,
                                                                 qp.point(), n,
                                                                 mu, lam);
                matrix_type trace_op    = eval_cell_trace_operator(msh, cl, cell_infos.cell_degree(),
                                                              qp.point(), num_total_dofs);

                vector_type sigma_n_op = traction_op * n;
                vector_type trace_n_op = trace_op * n;

                vector_type P_u_op = sigma_n_op       - gamma_n_h * trace_n_op;
                vector_type P_v_op = theta*sigma_n_op - gamma_n_h * trace_n_op;

                const T P_u_value = P_u_op.dot(u_lin);

                if (bilateral || P_u_value <= T(0))
                    contact_jacobian += (qp.weight() / gamma_n_h)
                                        * (P_v_op * P_u_op.transpose());
            }
        }

        face_offset += nfd;
    }

    return contact_jacobian;
}

/* Tresca friction Jacobian, linearised at u_lin, GIVEN threshold s:
 *   (1/gamma_t) \int [ Ds(P_t(u)) P_t(du) ] . P_t^theta(v)
 *   Ds(q) = I                                if |q| <= s   (stick)
 *         = s/|q| (I - (q x q)/|q|^2)        otherwise     (slip)
 *  */
template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_vector_hho_tresca_jacobian(const Mesh&                     msh,
                                const typename Mesh::cell_type& cl,
                                const MeshDegreeInfo<Mesh>&     degree_infos,
                                const Matrix<typename Mesh::coordinate_type,
                                             Dynamic, Dynamic>& grad_op,
                                const typename Mesh::coordinate_type gamma_t_0,
                                const typename Mesh::coordinate_type theta,
                                const typename Mesh::coordinate_type mu,
                                const typename Mesh::coordinate_type lam,
                                const typename Mesh::coordinate_type s_threshold,
                                const vector_boundary_conditions<Mesh>& bnd,
                                const Matrix<typename Mesh::coordinate_type,
                                             Dynamic, 1>& u_lin)
{
    using T           = typename Mesh::coordinate_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using tensor_type = static_matrix<T, Mesh::dimension, Mesh::dimension>;
    using tensor_vec  = static_vector<T, Mesh::dimension>;

    const size_t N = Mesh::dimension;

    const auto cell_infos = degree_infos.cellDegreeInfo(msh, cl);
    const auto grad_deg   = cell_infos.grad_degree();
    const auto quad_deg   = 2 * cell_infos.reconstruction_degree();
    const auto cbs        = vector_basis_size(cell_infos.cell_degree(), N, N);
    const auto fcs        = faces(msh, cl);
    const auto num_total_dofs = cbs + vector_faces_dofs(msh, cl, degree_infos);

    const tensor_type Id = tensor_type::Identity();

    auto gb = make_sym_matrix_monomial_basis(msh, cl, grad_deg);

    matrix_type tresca_jacobian =
        matrix_type::Zero(num_total_dofs, num_total_dofs);

    size_t face_offset = cbs;
    for (auto& fc : fcs)
    {
        const auto   fdi = degree_infos.degreeInfo(msh, fc);
        const size_t nfd = vector_face_dofs(msh, fdi);

        if (bnd.is_contact_face(fc))
        {

            const auto n   = normal(msh, cl, fc);
            const auto qps = integrate(msh, fc, quad_deg);

            const T gamma_t_h = gamma_t_0 / diameter(msh, fc);

            for (auto& qp : qps)
            {
                matrix_type traction_op = eval_traction_operator(msh, gb, grad_op,
                                                                 qp.point(), n,
                                                                 mu, lam);
                matrix_type trace_op    = eval_cell_trace_operator(msh, cl, cell_infos.cell_degree(),
                                                              qp.point(), num_total_dofs);

                matrix_type sigma_t_op = traction_op - (traction_op * n) * n.transpose();
                matrix_type trace_t_op = trace_op    - (trace_op * n)    * n.transpose();

                matrix_type P_t_u_op = sigma_t_op       - gamma_t_h * trace_t_op;
                matrix_type P_t_v_op = theta*sigma_t_op - gamma_t_h * trace_t_op;

                tensor_vec P_t_u_value = P_t_u_op.transpose() * u_lin;

                const T P_t_u_norm = P_t_u_value.norm();

                tensor_type Ds;
                if (s_threshold <= T(0))
                    Ds.setZero();          // [q]_0 = 0 identically: derivative 0
                else if (P_t_u_norm <= s_threshold)
                    Ds = Id;                                                 // stick
                else
                    Ds = (s_threshold / P_t_u_norm)
                         * (Id - (P_t_u_value * P_t_u_value.transpose())
                                 / P_t_u_value.squaredNorm());               // slip

                tresca_jacobian += (qp.weight() / gamma_t_h)
                                   * (P_t_v_op * Ds * P_t_u_op.transpose());
            }
        }

        face_offset += nfd;
    }

    return tresca_jacobian;
}

/* ==========================================================================
 * Coulomb friction
 *
 *     s(u) = - F [ P_n(u) ]_-  ,   P_n(u) = sigma_n(u) - gamma_n u_dT.n
 *
 * Since [.]_- <= 0, s >= 0 automatically.  The threshold is FROZEN at the
 * outer Picard iterate u_picard
 *
 *     Picard:  u_picard <- u ;  solve Tresca with s(u_picard)  ->  u
 * ========================================================================*/

/* Coulomb friction term: Tresca rhs with s = -F [P_n(u_picard)]_-        */
template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
make_vector_hho_coulomb_rhs(const Mesh&                     msh,
                            const typename Mesh::cell_type& cl,
                            const MeshDegreeInfo<Mesh>&     degree_infos,
                            const Matrix<typename Mesh::coordinate_type,
                                         Dynamic, Dynamic>& grad_op,
                            const typename Mesh::coordinate_type gamma_n_0,
                            const typename Mesh::coordinate_type gamma_t_0,
                            const typename Mesh::coordinate_type theta,
                            const typename Mesh::coordinate_type mu,
                            const typename Mesh::coordinate_type lam,
                            const typename Mesh::coordinate_type friction_coeff,
                            const vector_boundary_conditions<Mesh>& bnd,
                            const Matrix<typename Mesh::coordinate_type,
                                         Dynamic, 1>& u_prev,
                            const Matrix<typename Mesh::coordinate_type,
                                         Dynamic, 1>& u_picard)
{
    using T           = typename Mesh::coordinate_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;
    using tensor_vec  = static_vector<T, Mesh::dimension>;

    const size_t N = Mesh::dimension;

    const auto cell_infos = degree_infos.cellDegreeInfo(msh, cl);
    const auto grad_deg   = cell_infos.grad_degree();
    const auto quad_deg   = 2 * cell_infos.reconstruction_degree();
    const auto cbs        = vector_basis_size(cell_infos.cell_degree(), N, N);
    const auto fcs        = faces(msh, cl);
    const auto num_total_dofs = cbs + vector_faces_dofs(msh, cl, degree_infos);

    auto gb = make_sym_matrix_monomial_basis(msh, cl, grad_deg);

    vector_type coulomb_rhs = vector_type::Zero(num_total_dofs);

    size_t face_offset = cbs;
    for (auto& fc : fcs)
    {
        const auto   fdi = degree_infos.degreeInfo(msh, fc);
        const size_t nfd = vector_face_dofs(msh, fdi);

        if (bnd.is_contact_face(fc))
        {

            const auto n   = normal(msh, cl, fc);
            const auto qps = integrate(msh, fc, quad_deg);

            const T h_F       = diameter(msh, fc);
            const T gamma_n_h = gamma_n_0 / h_F;
            const T gamma_t_h = gamma_t_0 / h_F;

            for (auto& qp : qps)
            {
                matrix_type traction_op = eval_traction_operator(msh, gb, grad_op,
                                                                 qp.point(), n,
                                                                 mu, lam);
                matrix_type trace_op    = eval_cell_trace_operator(msh, cl, cell_infos.cell_degree(),
                                                              qp.point(), num_total_dofs);

                // ---- Coulomb threshold, frozen at the Picard iterate ----
                vector_type sigma_n_op = traction_op * n;
                vector_type trace_n_op = trace_op * n;
                vector_type P_n_op = sigma_n_op - gamma_n_h * trace_n_op;

                const T P_n_picard  = P_n_op.dot(u_picard);
                const T s_threshold = -friction_coeff
                                      * std::min(P_n_picard, T(0));   // >= 0

                // ---- from here on: identical to Tresca ----
                matrix_type sigma_t_op = traction_op - sigma_n_op * n.transpose();
                matrix_type trace_t_op = trace_op    - trace_n_op * n.transpose();

                matrix_type P_t_u_op = sigma_t_op       - gamma_t_h * trace_t_op;
                matrix_type P_t_v_op = theta*sigma_t_op - gamma_t_h * trace_t_op;

                tensor_vec P_t_u_value = P_t_u_op.transpose() * u_prev;

                const T P_t_u_norm = P_t_u_value.norm();

                tensor_vec P_t_u_proj;
                if (P_t_u_norm <= s_threshold)
                    P_t_u_proj = P_t_u_value;                                // stick
                else
                    P_t_u_proj = (s_threshold / P_t_u_norm) * P_t_u_value;   // slip

                coulomb_rhs += (qp.weight() / gamma_t_h) * (P_t_v_op * P_t_u_proj);
            }
        }

        face_offset += nfd;
    }

    return coulomb_rhs;
}

/* Coulomb friction Jacobian: Tresca Jacobian with the frozen threshold.
 * Exact for the inner (Tresca) problem*/
template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_vector_hho_coulomb_jacobian(const Mesh&                     msh,
                                 const typename Mesh::cell_type& cl,
                                 const MeshDegreeInfo<Mesh>&     degree_infos,
                                 const Matrix<typename Mesh::coordinate_type,
                                              Dynamic, Dynamic>& grad_op,
                                 const typename Mesh::coordinate_type gamma_n_0,
                                 const typename Mesh::coordinate_type gamma_t_0,
                                 const typename Mesh::coordinate_type theta,
                                 const typename Mesh::coordinate_type mu,
                                 const typename Mesh::coordinate_type lam,
                                 const typename Mesh::coordinate_type friction_coeff,
                                 const vector_boundary_conditions<Mesh>& bnd,
                                 const Matrix<typename Mesh::coordinate_type,
                                              Dynamic, 1>& u_lin,
                                 const Matrix<typename Mesh::coordinate_type,
                                              Dynamic, 1>& u_picard)
{
    using T           = typename Mesh::coordinate_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;
    using tensor_type = static_matrix<T, Mesh::dimension, Mesh::dimension>;
    using tensor_vec  = static_vector<T, Mesh::dimension>;

    const size_t N = Mesh::dimension;

    const auto cell_infos = degree_infos.cellDegreeInfo(msh, cl);
    const auto grad_deg   = cell_infos.grad_degree();
    const auto quad_deg   = 2 * cell_infos.reconstruction_degree();
    const auto cbs        = vector_basis_size(cell_infos.cell_degree(), N, N);
    const auto fcs        = faces(msh, cl);
    const auto num_total_dofs = cbs + vector_faces_dofs(msh, cl, degree_infos);

    const tensor_type Id = tensor_type::Identity();

    auto gb = make_sym_matrix_monomial_basis(msh, cl, grad_deg);

    matrix_type coulomb_jacobian =
        matrix_type::Zero(num_total_dofs, num_total_dofs);

    size_t face_offset = cbs;
    for (auto& fc : fcs)
    {
        const auto   fdi = degree_infos.degreeInfo(msh, fc);
        const size_t nfd = vector_face_dofs(msh, fdi);

        if (bnd.is_contact_face(fc))
        {

            const auto n   = normal(msh, cl, fc);
            const auto qps = integrate(msh, fc, quad_deg);

            const T h_F       = diameter(msh, fc);
            const T gamma_n_h = gamma_n_0 / h_F;
            const T gamma_t_h = gamma_t_0 / h_F;

            for (auto& qp : qps)
            {
                matrix_type traction_op = eval_traction_operator(msh, gb, grad_op,
                                                                 qp.point(), n,
                                                                 mu, lam);
                matrix_type trace_op    = eval_cell_trace_operator(msh, cl, cell_infos.cell_degree(),
                                                              qp.point(), num_total_dofs);

                // ---- Coulomb threshold, frozen at the Picard iterate ----
                vector_type sigma_n_op = traction_op * n;
                vector_type trace_n_op = trace_op * n;
                vector_type P_n_op = sigma_n_op - gamma_n_h * trace_n_op;

                const T P_n_picard  = P_n_op.dot(u_picard);
                const T s_threshold = friction_coeff
                                      * std::max(-1.0*P_n_picard, T(0));   // >= 0

                // ---- from here on: identical to Tresca ----
                matrix_type sigma_t_op = traction_op - sigma_n_op * n.transpose();
                matrix_type trace_t_op = trace_op    - trace_n_op * n.transpose();

                matrix_type P_t_u_op = sigma_t_op       - gamma_t_h * trace_t_op;
                matrix_type P_t_v_op = theta*sigma_t_op - gamma_t_h * trace_t_op;

                tensor_vec P_t_u_value = P_t_u_op.transpose() * u_lin;

                const T P_t_u_norm = P_t_u_value.norm();

                tensor_type Ds;
                if (s_threshold <= T(0))
                    Ds.setZero();          // [q]_0 = 0 identically: derivative 0
                else if (P_t_u_norm <= s_threshold)
                    Ds = Id;                                                 // stick
                else
                    Ds = (s_threshold / P_t_u_norm)
                         * (Id - (P_t_u_value * P_t_u_value.transpose())
                                 / P_t_u_value.squaredNorm());               // slip

                coulomb_jacobian += (qp.weight() / gamma_t_h)
                                    * (P_t_v_op * Ds * P_t_u_op.transpose());
            }
        }

        face_offset += nfd;
    }

    return coulomb_jacobian;
}


/* ==========================================================================
 * Weak symmetry (sliding) boundary:  u.n = 0,  sigma_t = 0.
 *
 * The normal component is imposed a la Nitsche, the tangential component is free
 *
 *   consistency (goes into A with the same minus sign as the contact one):
 *       (theta/gamma_n) \int_\GammaS sigma_n(u) sigma_n(v)
 *   residual:
 *       (1/gamma_n) \int_\GammaS P_n(u_prev) P_n^theta(v)      (always active)
 *   jacobian:
 *       (1/gamma_n) \int_\GammaS P_n(du) P_n^theta(v)          (linear: constant)
 * ========================================================================*/

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_vector_hho_symmetry_nitsche(const Mesh&                     msh,
                                 const typename Mesh::cell_type& cl,
                                 const MeshDegreeInfo<Mesh>&     degree_infos,
                                 const Matrix<typename Mesh::coordinate_type,
                                              Dynamic, Dynamic>& grad_op,
                                 const typename Mesh::coordinate_type gamma_n_0,
                                 const typename Mesh::coordinate_type theta,
                                 const typename Mesh::coordinate_type mu,
                                 const typename Mesh::coordinate_type lam,
                                 const vector_boundary_conditions<Mesh>& bnd)
{
    using T           = typename Mesh::coordinate_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t N = Mesh::dimension;

    const auto cell_infos = degree_infos.cellDegreeInfo(msh, cl);
    const auto grad_deg   = cell_infos.grad_degree();
    const auto quad_deg   = 2 * cell_infos.reconstruction_degree();
    const auto cbs        = vector_basis_size(cell_infos.cell_degree(), N, N);
    const auto num_total_dofs = cbs + vector_faces_dofs(msh, cl, degree_infos);

    auto gb = make_sym_matrix_monomial_basis(msh, cl, grad_deg);

    matrix_type nitsche_matrix =
        matrix_type::Zero(num_total_dofs, num_total_dofs);

    for (auto& fc : faces(msh, cl))
    {
        if (!bnd.is_contact_face(fc))
            continue;

        const auto n   = normal(msh, cl, fc);
        const auto qps = integrate(msh, fc, quad_deg);

        const T gamma_n_h = gamma_n_0 / diameter(msh, fc);

        for (auto& qp : qps)
        {
            matrix_type traction_op = eval_traction_operator(msh, gb, grad_op,
                                                             qp.point(), n, mu, lam);
            vector_type sigma_n_op = traction_op * n;

            nitsche_matrix += qp.weight() * (theta / gamma_n_h)
                              * sigma_n_op * sigma_n_op.transpose();
        }
    }

    return nitsche_matrix;
}

/* Residual of the weak symmetry condition, always active (no [.]_-):
 * this is the bilateral normal term specialised to its own name.          */
template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
make_vector_hho_symmetry_rhs(const Mesh&                     msh,
                             const typename Mesh::cell_type& cl,
                             const MeshDegreeInfo<Mesh>&     degree_infos,
                             const Matrix<typename Mesh::coordinate_type,
                                          Dynamic, Dynamic>& grad_op,
                             const typename Mesh::coordinate_type gamma_n_0,
                             const typename Mesh::coordinate_type theta,
                             const typename Mesh::coordinate_type mu,
                             const typename Mesh::coordinate_type lam,
                             const vector_boundary_conditions<Mesh>& bnd,
                             const Matrix<typename Mesh::coordinate_type,
                                          Dynamic, 1>& u_prev)
{
    return make_vector_hho_contact_rhs(msh, cl, degree_infos, grad_op,
                                       gamma_n_0, theta, mu, lam, bnd, u_prev,
                                       /* bilateral = */ true);
}

/* Jacobian of the weak symmetry condition: linear, hence u-independent.   */
template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_vector_hho_symmetry_jacobian(const Mesh&                     msh,
                                  const typename Mesh::cell_type& cl,
                                  const MeshDegreeInfo<Mesh>&     degree_infos,
                                  const Matrix<typename Mesh::coordinate_type,
                                               Dynamic, Dynamic>& grad_op,
                                  const typename Mesh::coordinate_type gamma_n_0,
                                  const typename Mesh::coordinate_type theta,
                                  const typename Mesh::coordinate_type mu,
                                  const typename Mesh::coordinate_type lam,
                                  const vector_boundary_conditions<Mesh>& bnd)
{
    const Matrix<typename Mesh::coordinate_type, Dynamic, 1> u_dummy =
        Matrix<typename Mesh::coordinate_type, Dynamic, 1>::Zero(grad_op.cols());

    return make_vector_hho_contact_jacobian(msh, cl, degree_infos, grad_op,
                                            gamma_n_0, theta, mu, lam, bnd, u_dummy,
                                            /* bilateral = */ true);
}

} // namespace disk