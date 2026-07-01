/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2016-2026
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 *
 * This file is copyright of the following authors:
 * Guillaume Delay (C) 2020-2021         guillaume.delay@sorbonne-universite.fr
 * Sorbonne Universite
 * Laboratoire Jacques-Louis Lions (LJLL)
 *
 * Matteo Cicuttin (C) 2016-2026
 */

#include "Lagrange_basis.hpp"
#include "diskpp/methods/hho"

//////////////////////////////////////////////////////////////////////////////////
////////////////////////////   ASSEMBLY ROUTINES   ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

template<typename Mesh>
std::pair<Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
          Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>>
make_vector_hho_gradrec_Lag(const Mesh&                     msh,
                            const typename Mesh::cell_type& cl,
                            const hho_degree_info&          di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    const auto celdeg  = di.cell_degree();
    const auto facdeg  = di.face_degree();
    const auto graddeg = di.grad_degree();

    const auto cb = make_scalar_Lagrange_basis(msh, cl, celdeg);
    const auto gb = make_vector_Lagrange_basis(msh, cl, graddeg);

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
        const auto fb = make_scalar_Lagrange_basis(msh, fc, facdeg);

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



// we compute the stabilisation 1/h_F(uF-pi^k_F(uT), vF-pi^k_F(vT))_F
template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_scalar_hdg_stabilization_Lag(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
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

    auto cb = make_scalar_Lagrange_basis(msh, cl, celdeg);
    const auto fcs = faces(msh, cl);

    for (size_t i = 0; i < num_faces; i++)
    {
        const auto fc = fcs[i];
        const auto h  = diameter(msh, fc);
        auto fb = make_scalar_Lagrange_basis(msh, fc, facdeg);

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

