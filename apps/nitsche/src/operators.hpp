#pragma once
#include "asm.hpp"

enum hho_mode {
    standard,
    nitsche
};

template<typename Mesh>
using DM = disk::dynamic_matrix<typename Mesh::coordinate_type>;

template<typename Mesh>
using DV = disk::dynamic_vector<typename Mesh::coordinate_type>;

template<typename Mesh>
auto
hho_mixedhigh_symlapl(const Mesh& msh,
    const typename Mesh::cell_type& cl, size_t degree,
    hho_mode mode, const std::vector<bc>& bcs)
{
    const static size_t DIM = Mesh::dimension;
    static_assert(DIM==2 or DIM==3, "Symmetric laplacian: only DIM = 2 or DIM = 3");

    using scalar_type = typename Mesh::coordinate_type;
    /* Reconstruction space basis */
    auto rb = disk::make_vector_monomial_basis(msh, cl, degree+1);
    auto rbs = rb.size();

    auto cb = disk::make_vector_monomial_basis(msh, cl, degree+1);
    auto cbs = cb.size();

    auto fcs = faces(msh, cl);
    auto fbs = disk::vector_basis_size(degree, DIM-1, DIM);
    auto n_allfacedofs = fcs.size() * fbs;

    /* The local problem has rbs rows, minus the DIM dofs
     * for constants in each direction, plus the rigid
     * body constraint (1 row in 2D, 3 rows in 3D) */
    constexpr auto curldofs = ((DIM == 2) ? 1 : 3);
    auto lprows = rbs - DIM + curldofs;

    /* Stiffness */
    disk::dynamic_matrix<scalar_type> K =
        disk::dynamic_matrix<scalar_type>::Zero(rbs, rbs);

    /* Local problem LHS */
    disk::dynamic_matrix<scalar_type> LHS =
        disk::dynamic_matrix<scalar_type>::Zero(lprows, lprows);
    
    /* Local problem RHS */
    disk::dynamic_matrix<scalar_type> RHS =
        disk::dynamic_matrix<scalar_type>::Zero(lprows, cbs + n_allfacedofs);

    /* (symgrad, symgrad) */
    auto qps = disk::integrate(msh, cl, 2*degree);
    for (auto& qp : qps) {
        auto sym_dphi = rb.eval_sgradients(qp.point());
        K += qp.weight() * disk::priv::outer_product(sym_dphi, sym_dphi);
    }
    LHS.block(0, 0, rbs-DIM, rbs-DIM) = K.bottomRightCorner(rbs-DIM, rbs-DIM);
    RHS.block(0, 0, rbs-DIM, cbs) = K.block(DIM, 0, rbs-DIM, cbs);

    /* Rigid body constraints (LHS), lagrange multiplier then eqn */
    auto qps_lm = disk::integrate(msh, cl, degree);
    for (auto& qp : qps_lm) {
        auto rphi = rb.eval_curls(qp.point());
        LHS.block(0, LHS.cols()-curldofs, rbs-DIM, curldofs) =
            rphi.bottomRows(rbs-DIM);
    }
    LHS.block(LHS.cols()-curldofs, 0, curldofs, rbs-DIM) =
        LHS.block(0, LHS.cols()-curldofs, rbs-DIM, curldofs).transpose();
    
    /* Now the face stuff */
    for (size_t fcnum = 0; fcnum < fcs.size(); fcnum++) {
        auto fcofs = cbs+fbs*fcnum;
        const auto& fc = fcs[fcnum];
        auto bi = msh.boundary_info(fc);
        auto fcid = offset(msh, fc);

        /* Dirichlet boundary conditions are always applied HHO-style,
         * so we don't need to add additional face-based terms to the
         * bilinear form and we do not have user-specified penalization
         * parameter. Neumann boundary conditions are applied either
         * HHO-style or Nitsche-style, depending on the user choice.
         * Don't forget boundary HHO terms on the internal interfaces.*/
        bool std_hho = (mode == hho_mode::standard);
        bool nitsche_not_neumann = (mode == hho_mode::nitsche) and (bcs[fcid] != bc::neumann);
        if (std_hho or nitsche_not_neumann) {
            auto fb = disk::make_vector_monomial_basis(msh, fc, degree);
            auto n = normal(msh, cl, fc);
            auto fqps = disk::integrate(msh, fc, 2*degree+1);
            for (auto& qp : fqps) {
                auto c_phi = cb.eval_functions(qp.point());
                auto f_phi = fb.eval_functions(qp.point());
                auto r_dphi = rb.eval_sgradients(qp.point());
                auto r_dphi_n = disk::priv::inner_product(r_dphi, (qp.weight()*n).eval());
                RHS.block(0, fcofs, rbs-DIM, fbs) +=
                    disk::priv::outer_product(r_dphi_n, f_phi).bottomRows(rbs-DIM);
                RHS.block(0, 0, rbs-DIM, cbs) -=
                    disk::priv::outer_product(r_dphi_n, c_phi).block(DIM, 0, rbs-DIM, cbs);
            }
        } else {
            /* Nothing to do */
        }
    }

    Eigen::FullPivLU<disk::dynamic_matrix<scalar_type>> fact(LHS);
    disk::dynamic_matrix<scalar_type> sol = fact.solve(RHS);
    //if (fact.info() != Eigen::Success) {
    //    std::cout << "Symgrad: cannot factorize local matrix" << std::endl;
    //}
    disk::dynamic_matrix<scalar_type> oper = sol.topRows(rbs-DIM);
    disk::dynamic_matrix<scalar_type> data = RHS.topRows(rbs-DIM).transpose()*oper;
    return std::pair(oper, data);
}

template<typename Mesh>
auto
hho_mixedhigh_divrec(const Mesh& msh,
    const typename Mesh::cell_type& cl, size_t degree,
    hho_mode mode, const std::vector<bc>& bcs)
{
    const static size_t DIM = Mesh::dimension;

    using scalar_type = typename Mesh::coordinate_type;
    /* Reconstruction space basis */
    auto rb = disk::make_scalar_monomial_basis(msh, cl, degree+1);
    auto rbs = rb.size();

    auto cb = disk::make_vector_monomial_basis(msh, cl, degree+1);
    auto cbs = cb.size();

    auto fcs = faces(msh, cl);
    auto fbs = disk::vector_basis_size(degree, DIM-1, DIM);
    auto n_allfacedofs = fcs.size() * fbs;

    disk::dynamic_matrix<scalar_type> LHS =
        disk::dynamic_matrix<scalar_type>::Zero(rbs, rbs);

    disk::dynamic_matrix<scalar_type> RHS =
        disk::dynamic_matrix<scalar_type>::Zero(rbs, cbs + n_allfacedofs);
    
    auto qps = disk::integrate(msh, cl, 2*degree+2);
    for (auto& qp : qps) {
        auto scalphi = rb.eval_functions(qp.point());
        LHS += qp.weight() * scalphi * scalphi.transpose();

        auto dscalphi = rb.eval_gradients(qp.point());
        auto vecphi = cb.eval_functions(qp.point());

        RHS.leftCols(cbs) -= qp.weight() * dscalphi * vecphi.transpose();
    }

    for (size_t fcnum = 0; fcnum < fcs.size(); fcnum++) {
        auto fcofs = cbs+fbs*fcnum;
        const auto& fc = fcs[fcnum];
        auto fcid = offset(msh, fc);
        auto bi = msh.boundary_info(fc);
        /* Dirichlet boundary conditions are always applied HHO-style,
         * so we don't need to add additional face-based terms to the
         * bilinear form and we do not have user-specified penalization
         * parameter. Neumann boundary conditions are applied either
         * HHO-style or Nitsche-style, depending on the user choice.
         * Don't forget boundary HHO terms on the internal interfaces.*/
         bool std_hho = (mode == hho_mode::standard);
         bool nitsche_not_neumann = (mode == hho_mode::nitsche) and (bcs[fcid] != bc::neumann);
         if (std_hho or nitsche_not_neumann) {   
            auto fb = disk::make_vector_monomial_basis(msh, fc, degree);
            auto n = normal(msh, cl, fc);
            auto fqps = disk::integrate(msh, fc, 2*degree+1);
            for (auto& qp : fqps) {
                auto r_phi = rb.eval_functions(qp.point());
                auto f_phi = fb.eval_functions(qp.point());
                
                RHS.block(0, fcofs, rbs, fbs) +=
                    qp.weight()*r_phi*(f_phi*n).transpose();
            }
        }
        else {
            /* Divergence reconstruction needs to be adjusted
             * on Nitsche-Neumann boundaries */
            auto n = normal(msh, cl, fc);
            auto fqps = disk::integrate(msh, fc, 2*degree+2);
            for (auto& qp : fqps) {
                auto r_phi = rb.eval_functions(qp.point());
                auto c_phi = cb.eval_functions(qp.point());
                
                RHS.block(0, 0, rbs, cbs) +=
                    qp.weight()*r_phi*(c_phi*n).transpose();
            }
        }
    }

    Eigen::LLT<disk::dynamic_matrix<scalar_type>> fact(LHS);
    if (fact.info() != Eigen::Success) {
        std::cout << "Divrec: cannot factorize local matrix" << std::endl;
    }
    disk::dynamic_matrix<scalar_type> oper = fact.solve(RHS);
    disk::dynamic_matrix<scalar_type> data = RHS.transpose() * oper; 

    return std::pair(oper, data);
}

template<typename Mesh>
disk::dynamic_matrix<typename Mesh::coordinate_type>
vstab(const Mesh& msh,
    const typename Mesh::cell_type& cl, size_t degree,
    hho_mode mode, const std::vector<bc>& bcs)
{
    /* Nitsche-HHO as implemented here is mixed-order (k+1 on cells
     * and k on faces). We use a standard Lehrenfeld-Schoeberl
     * stabilization. We need to stabilize only on the internal
     * interfaces, not on the domain boundary. */
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto cb = disk::make_vector_monomial_basis(msh, cl, degree+1);
    const auto cbs = cb.size();

    const auto fcs = faces(msh, cl);
    const auto fbs = disk::vector_basis_size(degree, Mesh::dimension-1, Mesh::dimension);
    const auto num_faces_dofs = fbs*fcs.size();
    const auto total_dofs     = cbs + num_faces_dofs;

    matrix_type data = matrix_type::Zero(total_dofs, total_dofs);

    T hT = diameter(msh, cl);
    T stabparam = 1.0/hT;

    for (size_t i = 0; i < fcs.size(); i++) {
        size_t fcofs = cbs+i*fbs;
        const auto fc = fcs[i];
        auto fcid = offset(msh, fc);
        stabparam = 1./diameter(msh, fc);

        auto bi = msh.boundary_info(fc);
        if ( (mode == hho_mode::nitsche) and 
             bi.is_boundary() and (bcs[fcid] == bc::neumann) ) {
            /* Nitsche-HHO boundary faces don't need stabilization
             * contributions, as there are no face unknowns there.
             * All other faces do (remember that we do Dirichlet
             * HHO-style)
             */
            continue;
        }

        /* Compute standard L-S stabilization otherwise. */
        const auto fb  = make_vector_monomial_basis(msh, fc, degree);

        const matrix_type If    = matrix_type::Identity(fbs, fbs);
        matrix_type       oper  = matrix_type::Zero(fbs, total_dofs);
        matrix_type       mass  = matrix_type::Zero(fbs, fbs);
        matrix_type       trace = matrix_type::Zero(fbs, cbs);

        const auto qps = integrate(msh, fc, 2*degree+1);
        for (auto& qp : qps)
        {
            const auto c_phi = cb.eval_functions(qp.point());
            const auto f_phi = fb.eval_functions(qp.point());

            assert(c_phi.rows() == cbs);
            assert(f_phi.rows() == fbs);
            assert(c_phi.cols() == f_phi.cols());
            mass += (qp.weight() * f_phi) * f_phi.transpose();
            trace += (qp.weight() * f_phi) * c_phi.transpose();
        }

        oper.block(0, 0, fbs, cbs) = mass.ldlt().solve(trace);
        oper.block(0, fcofs, fbs, fbs) = -If;
        data += oper.transpose() * mass * oper * stabparam;
    }

    return data;
}