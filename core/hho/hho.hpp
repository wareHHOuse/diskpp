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
#include "bases/bases_utils.hpp"
#include "bases/bases_ranges.hpp"
#include "timecounter.h"

namespace disk {

/* Build a matrix starting from the two specified precomputed bases. */

template<typename PrecompBasisA, typename PrecompBasisB>
class local_matrix
{
    typedef PrecompBasisA                                   precomputed_basis_a_type;
    typedef PrecompBasisB                                   precomputed_basis_b_type;
    typedef typename precomputed_basis_a_type::value_type   a_prebasis_value_type;
    typedef typename precomputed_basis_b_type::value_type   b_prebasis_value_type;

    static_assert(std::is_same<typename precomputed_basis_a_type::scalar_type,
                               typename precomputed_basis_b_type::scalar_type>::value,
                  "The two bases do not have the same scalar type");

    typedef typename precomputed_basis_a_type::scalar_type  scalar_type;
    typedef dynamic_matrix<scalar_type>                     matrix_type;

    std::shared_ptr<matrix_type>                            m_local_matrix;

    matrix_type ret;

    template<typename Function>
    void
    compute_impl(const precomputed_basis_a_type& pba,
                 const dof_range& a_range,
                 const precomputed_basis_b_type& pbb,
                 const dof_range& b_range,
                 const Function& mm_prod_lambda)
    {
        /* pba and pbb must be evaluated at the same points! */
        assert(pba.size_eval_points() == pbb.size_eval_points());
        size_t eval_points = pba.size_eval_points();

        auto qpts_itor = pba.eval_points_begin();

        auto a_bs = pba.size();
        auto a_begin = pba.data_begin();
        auto a_end = pba.data_begin();
        std::advance(a_begin, a_range.min());
        std::advance(a_end, a_range.max());

        auto b_bs = pbb.size();
        auto b_begin = pbb.data_begin();
        auto b_end = pbb.data_begin();
        std::advance(b_begin, b_range.min());
        std::advance(b_end, b_range.max());

        (*m_local_matrix) = dynamic_matrix<scalar_type>::Zero(b_range.size(), a_range.size());

        for (size_t qpi = 0; qpi < eval_points; qpi++)
        {
            auto weight = (*qpts_itor++).weight();

            size_t i = 0;
            for (auto itor = b_begin; itor != b_end; itor++, i++)
            {
                size_t j = 0;
                for (auto jtor = a_begin; jtor != a_end; jtor++, j++)
                {
                    assert(i < b_range.size());
                    assert(j < a_range.size());
                    (*m_local_matrix)(i,j) += weight * mm_prod_lambda( *jtor, *itor );
                }
            }

            std::advance(a_begin, a_bs);
            std::advance(a_end,   a_bs);
            std::advance(b_begin, b_bs);
            std::advance(b_end,   b_bs);
        }
    }



public:
    local_matrix()
        //: m_local_matrix(nullptr)
    {
        m_local_matrix = std::make_shared<matrix_type>();
    }

    void
    compute(const precomputed_basis_a_type& pba,
            const precomputed_basis_b_type& pbb)
    {
        auto tp = [](const a_prebasis_value_type& a,
                     const b_prebasis_value_type& b) -> auto {
            return mm_prod(a, b);
        };

        compute_impl(pba, pba.full_range(), pbb, pbb.full_range(), tp);
    }

    void
    compute(const precomputed_basis_a_type& pba,
            const dof_range& a_range,
            const precomputed_basis_b_type& pbb,
            const dof_range& b_range)
    {
        auto tp = [](const a_prebasis_value_type& a,
                     const b_prebasis_value_type& b) -> auto {
            return mm_prod(a, b);
        };

        compute_impl(pba, a_range, pbb, b_range, tp);
    }

    template<typename NormalVec>
    void
    compute(const precomputed_basis_a_type& pba,
            const precomputed_basis_b_type& pbb,
            const NormalVec& n)
    {
        auto tp = [&](const a_prebasis_value_type& a,
                      const b_prebasis_value_type& b) -> auto {
            auto t = mm_prod(b, n);
            return mm_prod(a, t);
        };

        compute_impl(pba, pba.full_range(), pbb, pbb.full_range(), tp);
    }

    template<typename NormalVec>
    void
    compute(const precomputed_basis_a_type& pba,
            const dof_range& a_range,
            const precomputed_basis_b_type& pbb,
            const dof_range& b_range,
            const NormalVec& n)
    {
        auto tp = [&](const a_prebasis_value_type& a,
                      const b_prebasis_value_type& b) -> auto {
            auto t = mm_prod(b, n);
            return mm_prod(a, t);
        };

        compute_impl(pba, a_range, pbb, b_range, tp);
    }

    matrix_type
    get(void) const
    {
        assert(m_local_matrix);
        return *m_local_matrix;
    }
};

/* Build a RHS starting from the specified precomputed basis */
template<typename PrecompBasis>
class local_rhs
{
    typedef PrecompBasis                                    precomputed_basis_type;
    typedef typename precomputed_basis_type::value_type     prebasis_value_type;
    typedef typename precomputed_basis_type::scalar_type    scalar_type;
    typedef dynamic_vector<scalar_type>                     vector_type;

    std::shared_ptr<vector_type>                            m_local_rhs;

public:
    local_rhs()
    {
        m_local_rhs = std::make_shared<vector_type>();
    }

    template<typename Function>
    void
    compute(const precomputed_basis_type& pb, const Function& f)
    {
        m_local_rhs->resize(pb.size());

        auto qpts_itor = pb.eval_points_begin();
        prebasis_value_type function_val = (*qpts_itor).weight() * f((*qpts_itor).point());

        auto begin = pb.data_begin();
        auto end = pb.data_begin();
        std::advance(end, pb.size());

        size_t vi = 0;
        for (auto itor = begin; itor != end; itor++, vi++)
            (*m_local_rhs)(vi) = mm_prod(function_val, *itor);


        for (size_t i = 1; i < pb.size_eval_points(); i++)
        {
            std::advance(qpts_itor, 1);
            std::advance(begin, pb.size());
            std::advance(end, pb.size());
            function_val = (*qpts_itor).weight() * f((*qpts_itor).point());
            vi = 0;
            for (auto itor = begin; itor != end; itor++, vi++)
                (*m_local_rhs)(vi) += mm_prod(function_val, *itor);
        }
    }

    vector_type
    get(void) const
    {
        assert(m_local_rhs);
        return *m_local_rhs;
    }
};

/* Project the function f on the cell pertaining to LocalData.
 */
template<typename LocalData, typename Function>
dynamic_vector<typename LocalData::scalar_type>
project(LocalData& ld, size_t degree, const Function& f)
{
    typedef typename LocalData::scalar_type     scalar_type;

    auto proj_stuff     = ld.get_projector_stuff();
    auto mm_wrapper     = std::get<0>(proj_stuff);
    auto cell_basis     = std::get<1>(proj_stuff);

    local_rhs<decltype(cell_basis)> rhs_wrapper;
    rhs_wrapper.compute(cell_basis, f);

    auto dof_range      = cell_basis.range(0, degree);

    dynamic_matrix<scalar_type> mass = take(mm_wrapper.get(), dof_range, dof_range);
    dynamic_vector<scalar_type> rhs = take(rhs_wrapper.get(), dof_range);

    return mass.llt().solve(rhs);
}

template<typename LocalData, typename Function>
dynamic_vector<typename LocalData::scalar_type>
project(LocalData& ld, const typename LocalData::cell_type& cl, const Function& f)
{
    typedef typename LocalData::scalar_type     scalar_type;
    typedef typename LocalData::mesh_type       mesh_type;

    mesh_type msh           = ld.get_mesh();
    auto degree             = ld.get_degree();

    auto proj_stuff         = ld.get_projector_stuff();
    auto cell_mm_wrapper    = std::get<0>(proj_stuff);
    auto cell_basis         = std::get<1>(proj_stuff);
    auto face_mm_wrappers   = std::get<2>(proj_stuff);
    auto face_bases         = std::get<3>(proj_stuff);

    typename LocalData::cell_rhs_type       cell_rhs_wrapper;
    cell_rhs_wrapper.compute(cell_basis, f);

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    dofspace_ranges dsr(ld);

    dynamic_vector<scalar_type> projection;
    projection.resize(dsr.total_size());

    auto cell_range = dsr.cell_range();
    dynamic_matrix<scalar_type> mass = take(cell_mm_wrapper.get(), cell_range, cell_range);
    dynamic_vector<scalar_type> rhs = take(cell_rhs_wrapper.get(), cell_range);

    projection.head(cell_range.size()) = mass.ldlt().solve(rhs);

    for (size_t i = 0; i < num_faces; i++)
    {
        typename LocalData::face_rhs_type   face_rhs_wrapper;
        face_rhs_wrapper.compute(face_bases[i], f);

        dynamic_matrix<scalar_type> f_mass = face_mm_wrappers[i].get();
        dynamic_vector<scalar_type> f_rhs = face_rhs_wrapper.get();

        auto face_range = dsr.face_range(i);

        projection.head(face_range.max()).tail(face_range.size()) =
            f_mass.llt().solve(f_rhs);
    }

    return projection;
}

template<typename CellBasisType, typename CellQuadType, typename Mesh,
         typename Function>
dynamic_vector<typename Mesh::scalar_type>
compute_rhs(const Mesh& msh, const typename Mesh::cell& cl,
            const Function& f, size_t degree)
{
    typedef dynamic_vector<typename Mesh::scalar_type> vector_type;

    auto cell_basis     = CellBasisType(degree);
    auto cell_quad      = CellQuadType(2*degree);

    vector_type ret = vector_type::Zero(cell_basis.size());

    auto cell_quadpoints = cell_quad.integrate(msh, cl);
    for (auto& qp : cell_quadpoints)
    {
        auto phi = cell_basis.eval_functions(msh, cl, qp.point());
        auto fval = f(qp.point());
        for (size_t i = 0; i < cell_basis.size(); i++)
            ret(i) += qp.weight() * mm_prod(fval, phi[i]);
    }

    return ret;
}

template<typename Mesh, typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType>
class projector_nopre
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    typedef CellBasisType                       cell_basis_type;
    typedef CellQuadType                        cell_quadrature_type;
    typedef FaceBasisType                       face_basis_type;
    typedef FaceQuadType                        face_quadrature_type;

    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;

    cell_basis_type                             cell_basis;
    cell_quadrature_type                        cell_quadrature;

    face_basis_type                             face_basis;
    face_quadrature_type                        face_quadrature;

    size_t                                      m_degree;

public:

    projector_nopre()
        : m_degree(1)
    {
        cell_basis          = cell_basis_type(m_degree);
        cell_quadrature     = cell_quadrature_type(2*m_degree);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }

    projector_nopre(size_t degree)
        : m_degree(degree)
    {
        cell_basis          = cell_basis_type(m_degree);
        cell_quadrature     = cell_quadrature_type(2*m_degree);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }

    matrix_type cell_mm;

    template<typename Function>
    vector_type
    compute_cell(const mesh_type& msh, const cell_type& cl, const Function& f)
    {
        matrix_type mm = matrix_type::Zero(cell_basis.size(), cell_basis.size());
        vector_type rhs = vector_type::Zero(cell_basis.size());

        auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
        for (auto& qp : cell_quadpoints)
        {
            auto phi = cell_basis.eval_functions(msh, cl, qp.point());
            auto fval = f(qp.point());
            for (size_t i = 0; i < phi.size(); i++)
            {
                for (size_t j = 0; j < phi.size(); j++)
                    mm(i,j) += qp.weight() * mm_prod(phi[i], phi[j]);

                rhs(i) += qp.weight() * mm_prod(fval, phi[i]);
            }
        }

        cell_mm = mm;
        return mm.llt().solve(rhs);
    }

    template<typename Function>
    vector_type
    compute_whole(const mesh_type& msh, const cell_type& cl, const Function& f)
    {
        auto fcs = faces(msh, cl);
        vector_type ret(cell_basis.size() + fcs.size()*face_basis.size());

        ret.block(0, 0, cell_basis.size(), 1) = compute_cell(msh, cl);

        size_t face_offset = cell_basis.size();
        for (auto& fc : fcs)
        {
            matrix_type mm = matrix_type::Zero(face_basis.size(), face_basis.size());
            vector_type rhs = vector_type::Zero(face_basis.size());

            auto face_quadpoints = face_quadrature.integrate(msh, cl);
            for (auto& qp : face_quadpoints)
            {
                auto phi = face_basis.eval_functions(msh, cl, qp.point());
                auto fval = f(qp.point());
                for (size_t i = 0; i < phi.size(); i++)
                {
                    for (size_t j = 0; j < phi.size(); j++)
                        mm(i,j) += qp.weight() * mm_prod(phi[i], phi[j]);

                    rhs(i) += qp.weight() * mm_prod(fval, phi[i]);
                }
            }

            ret.block(face_offset, 0, face_basis.size(), 1) = mm.llt().solve(rhs);
            face_offset += face_basis.size();
        }

        return ret;
    }
};

template<typename LocalData>
class gradient_reconstruction
{
    typedef typename LocalData::cpbf_type       cpbf_type;
    typedef typename LocalData::cpbg_type       cpbg_type;
    typedef typename LocalData::fpbf_type       fpbf_type;

    typedef typename LocalData::scalar_type     scalar_type;
    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;

    typedef typename LocalData::mesh_type       mesh_type;

public:
    matrix_type     oper;
    matrix_type     data;

    gradient_reconstruction()
    {}

    void
    compute(LocalData& ld)
    {
        auto msh            = ld.get_mesh();
        auto cl             = ld.get_cell();
        auto degree         = ld.get_degree();

        auto gr_stuff       = ld.get_gr_stuff();
        auto gr_matrix      = std::get<0>(gr_stuff);    //stiffness matrix
        auto cell_basis     = std::get<1>(gr_stuff);    //cell basis pre-eval'd on cell
        auto face_bases     = std::get<2>(gr_stuff);    //face basis pre-eval'd on faces
        auto cell_basis_f   = std::get<3>(gr_stuff);    //cell basis pre-eval'd on faces
        auto cell_grads_f   = std::get<4>(gr_stuff);    //cell grads pre-eval'd on faces

        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();
        assert(face_bases.size() > 0 && face_bases.size() == num_faces);
        assert(cell_basis_f.size() == num_faces);
        assert(cell_grads_f.size() == num_faces);

        /* LHS: take basis functions derivatives from degree 1 to K+1 */
        auto MG_rowcol_range = cell_basis.range(1, degree+1);
        matrix_type MG = take(gr_matrix.get(), MG_rowcol_range, MG_rowcol_range);

        /* RHS, volumetric part. */
        auto BG_row_range = cell_basis.range(1, degree+1);
        auto BG_col_range = cell_basis.range(0, degree);

        dofspace_ranges dsr(ld);

        matrix_type BG = matrix_type::Zero(BG_row_range.size(), dsr.total_size());

        BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
                        take(gr_matrix.get(), BG_row_range, BG_col_range);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto current_face_range = dsr.face_range(face_i);
            auto fc = fcs[face_i];
            auto n = normal(msh, cl, fc);

            local_matrix<cpbf_type, cpbg_type> lm1;
            lm1.compute(cell_basis_f[face_i], BG_col_range,
                        cell_grads_f[face_i], BG_row_range, n);

            BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) -= lm1.get();

            local_matrix<fpbf_type, cpbg_type> lm2;
            lm2.compute(face_bases[face_i], face_bases[face_i].full_range(),
                        cell_grads_f[face_i], BG_row_range, n);

            BG.block(0, current_face_range.min(),
                     BG_row_range.size(), current_face_range.size()) += lm2.get();

        }

        oper  = MG.llt().solve(BG);    // GT
        data  = BG.transpose() * oper;  // A
    }
};

template<typename Mesh, typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType>
class gradient_reconstruction_nopre
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    typedef CellBasisType                       cell_basis_type;
    typedef CellQuadType                        cell_quadrature_type;
    typedef FaceBasisType                       face_basis_type;
    typedef FaceQuadType                        face_quadrature_type;

    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;

    cell_basis_type                             cell_basis;
    cell_quadrature_type                        cell_quadrature;

    face_basis_type                             face_basis;
    face_quadrature_type                        face_quadrature;

    size_t                                      m_degree;

public:
    matrix_type     oper;
    matrix_type     data;

    gradient_reconstruction_nopre()
        : m_degree(1)
    {
        cell_basis          = cell_basis_type(m_degree+1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }

    gradient_reconstruction_nopre(size_t degree)
        : m_degree(degree)
    {
        cell_basis          = cell_basis_type(m_degree+1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }

    void compute(const mesh_type& msh, const cell_type& cl)
    {
        matrix_type stiff_mat = matrix_type::Zero(cell_basis.size(), cell_basis.size());

        auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
        for (auto& qp : cell_quadpoints)
        {
            auto c_dphi = cell_basis.eval_gradients(msh, cl, qp.point());
            for (size_t i = 0; i < cell_basis.size(); i++)
                for (size_t j = 0; j < cell_basis.size(); j++)
                    stiff_mat(i,j) += qp.weight() * mm_prod(c_dphi[i], c_dphi[j]);
        }

        /* LHS: take basis functions derivatives from degree 1 to K+1 */
        auto MG_rowcol_range = cell_basis.range(1, m_degree+1);
        matrix_type MG = take(stiff_mat, MG_rowcol_range, MG_rowcol_range);

        /* RHS, volumetric part. */
        auto BG_row_range = cell_basis.range(1, m_degree+1);
        auto BG_col_range = cell_basis.range(0, m_degree);

        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        auto num_cell_dofs = cell_basis.range(0, m_degree).size();
        auto num_face_dofs = face_basis.size();

        dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

        matrix_type BG = matrix_type::Zero(BG_row_range.size(), dsr.total_size());

        BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
                                take(stiff_mat, BG_row_range, BG_col_range);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto current_face_range = dsr.face_range(face_i);
            auto fc = fcs[face_i];
            auto n = normal(msh, cl, fc);
            auto face_quadpoints = face_quadrature.integrate(msh, fc);

            //matrix_type m1 = matrix_type::Zeros(BG_row_range.size(), BG_col_range.size());
            for (auto& qp : face_quadpoints)
            {
                auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());
                auto c_dphi = cell_basis.eval_gradients(msh, cl, qp.point());
                for (size_t i = BG_row_range.min(), ii = 0; i < BG_row_range.max(); i++, ii++)
                {
                    for (size_t j = BG_col_range.min(), jj = 0; j < BG_col_range.max(); j++, jj++)
                    {
                        auto p1 = mm_prod(c_dphi[i], n);
                        auto p2 = mm_prod(c_phi[j], p1);
                        BG(ii,jj) -= qp.weight() * p2;
                    }
                }

                auto f_phi = face_basis.eval_functions(msh, fc, qp.point());

                for (size_t i = BG_row_range.min(), ii = 0; i < BG_row_range.max(); i++, ii++)
                {
                    for (size_t j = 0, jj = current_face_range.min();
                                j < current_face_range.size();
                                j++, jj++)
                    {
                        auto p1 = mm_prod(c_dphi[i], n);
                        auto p2 = mm_prod(f_phi[j], p1);
                        BG(ii,jj) += qp.weight() * p2;
                    }
                }
            }
        }

        oper  = MG.llt().solve(BG);    // GT
        data  = BG.transpose() * oper;  // A
    }
};

template<typename Mesh, typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType,
                        typename DivCellBasisType, typename DivCellQuadType>
class divergence_reconstruction_nopre
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    typedef CellBasisType                       cell_basis_type;
    typedef CellQuadType                        cell_quadrature_type;
    typedef FaceBasisType                       face_basis_type;
    typedef FaceQuadType                        face_quadrature_type;
    typedef DivCellBasisType                    div_cell_basis_type;
    typedef DivCellQuadType                     div_cell_quadrature_type;

    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;

    cell_basis_type                             cell_basis;
    cell_quadrature_type                        cell_quadrature;

    face_basis_type                             face_basis;
    face_quadrature_type                        face_quadrature;

    div_cell_basis_type                         div_cell_basis;
    div_cell_quadrature_type                    div_cell_quadrature;

    size_t                                      m_degree;

public:
    matrix_type     oper;
    matrix_type     data;

    divergence_reconstruction_nopre()
        : m_degree(1)
    {
        cell_basis          = cell_basis_type(m_degree+1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
        div_cell_basis      = div_cell_basis_type(m_degree);
        div_cell_quadrature = div_cell_quadrature_type(2*m_degree);
    }

    divergence_reconstruction_nopre(size_t degree)
        : m_degree(degree)
    {
        cell_basis          = cell_basis_type(m_degree+1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
        div_cell_basis      = div_cell_basis_type(m_degree);
        div_cell_quadrature = div_cell_quadrature_type(2*m_degree);
    }

    void compute(const mesh_type& msh, const cell_type& cl)
    {
        auto dcbs = div_cell_basis.size();
        matrix_type MD = matrix_type::Zero(dcbs, dcbs);

        auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
        for (auto& qp : cell_quadpoints)
        {
            auto phi_d = div_cell_basis.eval_functions(msh, cl, qp.point());
            for (size_t i = 0; i < dcbs; i++)
                for (size_t j = 0; j < dcbs; j++)
                    MD(i,j) += qp.weight() * mm_prod(phi_d[i], phi_d[j]);
        }

        /* RHS, volumetric part. */
        auto fcs = faces(msh, cl);
        auto num_cell_dofs = cell_basis.range(0, m_degree).size();
        auto num_face_dofs = face_basis.size();
        auto num_faces = fcs.size();

        dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

        matrix_type BD = matrix_type::Zero(dcbs, dsr.total_size());
        for (auto& qp : cell_quadpoints)
        {
            auto phi = cell_basis.eval_functions(msh, cl, qp.point());
            auto dphi_d = div_cell_basis.eval_gradients(msh, cl, qp.point());

            for (size_t i = 0; i < dphi_d.size(); i++)
                for (size_t j = 0; j < num_cell_dofs; j++)
                    BD(i,j) -= qp.weight() * mm_prod(dphi_d[i], phi[j]);
        }

        size_t face_offset = num_cell_dofs;

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];
            auto n = normal(msh, cl, fc);
            auto face_quadpoints = face_quadrature.integrate(msh, fc);

            for (auto& qp : face_quadpoints)
            {
                auto phi_d = div_cell_basis.eval_functions(msh, cl, qp.point());
                auto phi = face_basis.eval_functions(msh, fc, qp.point());
                for (size_t i = 0; i < phi_d.size(); i++)
                {
                    for (size_t j = 0; j < face_basis.size(); j++)
                    {
                        auto p1 = mm_prod(phi[j], n);
                        scalar_type p2 = mm_prod(p1, phi_d[i]);
                        BD(i,face_offset+j) += qp.weight() * p2;
                    }
                }
            }

            face_offset += face_basis.size();
        }

        oper = MD.partialPivLu().solve(BD);
        data = BD.transpose() * oper;
    }
};

template<typename LocalData>
class diffusion_like_stabilization
{
    typedef typename LocalData::scalar_type     scalar_type;
    typedef dynamic_matrix<scalar_type>         matrix_type;

public:
    dynamic_matrix<scalar_type> data;

    diffusion_like_stabilization()
    {}

    void
    compute(LocalData& ld, const matrix_type& GT)
    {
        auto msh    = ld.get_mesh();
        auto cl     = ld.get_cell();
        auto degree = ld.get_degree();

        auto stab_stuff         = ld.get_stab_stuff();
        auto stab_matrix        = std::get<0>(stab_stuff);  //cell mass matrix
        auto cell_basis         = std::get<1>(stab_stuff);  //cell basis pre-eval'd on cell
        auto face_bases         = std::get<2>(stab_stuff);  //face basis pre-eval'd on faces
        auto stab_f_matrices    = std::get<3>(stab_stuff);  //face mass matrices
        auto stab_tr_matrices   = std::get<4>(stab_stuff);  //face trace matrices

        auto zero_range         = cell_basis.range(0, degree);
        auto one_range          = cell_basis.range(1, degree+1);

        // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

        //Step 1: compute \pi_T^k p_T^k v (third term).
        matrix_type M1 = take(stab_matrix.get(), zero_range, zero_range);
        matrix_type M2 = take(stab_matrix.get(), zero_range, one_range);
        matrix_type proj1 = -M1.llt().solve(M2*GT);

        //Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
        matrix_type I_T = matrix_type::Identity(zero_range.size(), zero_range.size());
        proj1.block(0, 0, zero_range.size(), zero_range.size()) += I_T;

        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();
        dofspace_ranges dsr(ld);

        data = matrix_type::Zero(dsr.total_size(), dsr.total_size());

        // Step 3: project on faces (eqn. 21)
        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto current_face_range = dsr.face_range(face_i);
            auto fbs = current_face_range.size();
            auto h = diameter(msh, /*fcs[face_i]*/cl);

            Eigen::LLT<matrix_type> piKF;
            piKF.compute(stab_f_matrices[face_i].get());

            // Step 3a: \pi_F^k( v_F - p_T^k v )
            auto face_range = current_face_range.remove_offset();
            matrix_type MR1 = take(stab_tr_matrices[face_i].get(), face_range, one_range);
            matrix_type proj2 = piKF.solve(MR1*GT);
            matrix_type I_F = matrix_type::Identity(fbs, fbs);
            auto block_offset = current_face_range.min();
            proj2.block(0, block_offset, fbs, fbs) -= I_F;

            // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
            matrix_type MR2 = take(stab_tr_matrices[face_i].get(), face_range, zero_range);
            matrix_type proj3 = piKF.solve(MR2*proj1);

            matrix_type BRF = proj2 + proj3;

            data += BRF.transpose() * stab_f_matrices[face_i].get() * BRF / h;
        }
    }
};

template<typename Mesh, typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType>
class diffusion_like_stabilization_nopre
{
    typedef Mesh                            mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;

    typedef CellBasisType                   cell_basis_type;
    typedef CellQuadType                    cell_quadrature_type;
    typedef FaceBasisType                   face_basis_type;
    typedef FaceQuadType                    face_quadrature_type;

    typedef dynamic_matrix<scalar_type>     matrix_type;
    typedef dynamic_vector<scalar_type>     vector_type;

    cell_basis_type                         cell_basis;
    cell_quadrature_type                    cell_quadrature;

    face_basis_type                         face_basis;
    face_quadrature_type                    face_quadrature;

    size_t                                  m_degree;

public:
    matrix_type     data;

    diffusion_like_stabilization_nopre()
        : m_degree(1)
    {
        cell_basis          = cell_basis_type(m_degree+1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }

    diffusion_like_stabilization_nopre(size_t degree)
        : m_degree(degree)
    {
        cell_basis          = cell_basis_type(m_degree+1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }

    void compute(const mesh_type& msh, const cell_type& cl, const matrix_type& gradrec_oper)
    {
        matrix_type mass_mat = matrix_type::Zero(cell_basis.size(), cell_basis.size());

        auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
        for (auto& qp : cell_quadpoints)
        {
            auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());
            for (size_t i = 0; i < cell_basis.size(); i++)
                for (size_t j = 0; j < cell_basis.size(); j++)
                    mass_mat(i,j) += qp.weight() * mm_prod(c_phi[i], c_phi[j]);
        }

        auto zero_range         = cell_basis.range(0, m_degree);
        auto one_range          = cell_basis.range(1, m_degree+1);

        // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

        //Step 1: compute \pi_T^k p_T^k v (third term).
        matrix_type M1 = take(mass_mat, zero_range, zero_range);
        matrix_type M2 = take(mass_mat, zero_range, one_range);
        matrix_type proj1 = -M1.llt().solve(M2*gradrec_oper);

        //Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
        matrix_type I_T = matrix_type::Identity(zero_range.size(), zero_range.size());
        proj1.block(0, 0, zero_range.size(), zero_range.size()) += I_T;

        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        auto num_cell_dofs = cell_basis.range(0, m_degree).size();
        auto num_face_dofs = face_basis.size();

        dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

        data = matrix_type::Zero(dsr.total_size(), dsr.total_size());

        // Step 3: project on faces (eqn. 21)
        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto current_face_range = dsr.face_range(face_i);
            auto fbs = current_face_range.size();
            auto h = diameter(msh, /*fcs[face_i]*/cl);
            auto fc = fcs[face_i];

            matrix_type face_mass_matrix    = matrix_type::Zero(fbs, fbs);
            matrix_type face_trace_matrix   = matrix_type::Zero(fbs, cell_basis.size());

            auto face_quadpoints = face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
                auto f_phi = face_basis.eval_functions(msh, fc, qp.point());
                auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());

                for (size_t i = 0; i < face_basis.size(); i++)
                    for (size_t j = 0; j < face_basis.size(); j++)
                        face_mass_matrix(i,j) += qp.weight() * mm_prod(f_phi[i], f_phi[j]);

                for (size_t i = 0; i < face_basis.size(); i++)
                    for (size_t j = 0; j < cell_basis.size(); j++)
                        face_trace_matrix(i,j) += qp.weight() * mm_prod(f_phi[i], c_phi[j]);
            }

            Eigen::LLT<matrix_type> piKF;
            piKF.compute(face_mass_matrix);

            // Step 3a: \pi_F^k( v_F - p_T^k v )
            auto face_range = current_face_range.remove_offset();
            matrix_type MR1 = take(face_trace_matrix, face_range, one_range);
            matrix_type proj2 = piKF.solve(MR1*gradrec_oper);
            matrix_type I_F = matrix_type::Identity(fbs, fbs);
            auto block_offset = current_face_range.min();
            proj2.block(0, block_offset, fbs, fbs) -= I_F;

            // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
            matrix_type MR2 = take(face_trace_matrix, face_range, zero_range);
            matrix_type proj3 = piKF.solve(MR2*proj1);

            matrix_type BRF = proj2 + proj3;

            data += BRF.transpose() * face_mass_matrix * BRF / h;
        }
    }
};

template<typename LocalData>
class diffusion_like_static_condensation
{
    typedef typename LocalData::scalar_type     scalar_type;
    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;

public:
    diffusion_like_static_condensation()
    {}

    auto
    compute(LocalData& ld, const matrix_type& local_mat)
    {
        auto msh    = ld.get_mesh();
        //auto cl     = ld.get_cell();
        auto rhs    = ld.get_cell_rhs();

        dofspace_ranges dsr(ld);

        size_t cell_size = dsr.cell_range().size();
        size_t face_size = dsr.all_faces_range().size();
        assert(cell_size == ld.num_cell_dofs());

        matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
        matrix_type K_TF = local_mat.topRightCorner(cell_size, face_size);
        matrix_type K_FT = local_mat.bottomLeftCorner(face_size, cell_size);
        matrix_type K_FF = local_mat.bottomRightCorner(face_size, face_size);

        assert(K_TT.cols() == cell_size);
        assert(K_TT.cols() + K_TF.cols() == local_mat.cols());
        assert(K_TT.rows() + K_FT.rows() == local_mat.rows());
        assert(K_TF.rows() + K_FF.rows() == local_mat.rows());
        assert(K_FT.cols() + K_FF.cols() == local_mat.cols());

        auto K_TT_ldlt = K_TT.llt();
        matrix_type AL = K_TT_ldlt.solve(K_TF);
        vector_type bL = K_TT_ldlt.solve(rhs.get().block(0, 0, cell_size, 1));

        matrix_type AC = K_FF - K_FT * AL;
        vector_type bC = - K_FT * bL;

        return std::make_pair(AC, bC);
    }

    vector_type
    recover(LocalData& ld, const matrix_type& local_mat, const vector_type& solF)
    {
        auto msh    = ld.get_mesh();
        //auto cl     = ld.get_cell();
        auto rhs    = ld.get_cell_rhs();

        dofspace_ranges dsr(ld);

        size_t cell_size        = dsr.cell_range().size();
        size_t all_faces_size   = dsr.all_faces_range().size();

        vector_type ret( dsr.total_size() );

        matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
        matrix_type K_TF = local_mat.topRightCorner(cell_size, all_faces_size);

        vector_type solT = K_TT.llt().solve(rhs.get().block(0, 0, cell_size, 1) - K_TF*solF);

        ret.head(cell_size)         = solT;
        ret.tail(all_faces_size)    = solF;

        return ret;
    }

};

template<typename Mesh, typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType>
class diffusion_like_static_condensation_nopre
{
    typedef Mesh                            mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;

    typedef CellBasisType                   cell_basis_type;
    typedef CellQuadType                    cell_quadrature_type;
    typedef FaceBasisType                   face_basis_type;
    typedef FaceQuadType                    face_quadrature_type;

    typedef dynamic_matrix<scalar_type>     matrix_type;
    typedef dynamic_vector<scalar_type>     vector_type;

    cell_basis_type                         cell_basis;
    cell_quadrature_type                    cell_quadrature;

    face_basis_type                         face_basis;
    face_quadrature_type                    face_quadrature;

    size_t                                  m_degree;

public:
    diffusion_like_static_condensation_nopre()
        : m_degree(1)
    {
        cell_basis          = cell_basis_type(m_degree);
        cell_quadrature     = cell_quadrature_type(2*m_degree);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }

    diffusion_like_static_condensation_nopre(size_t degree)
        : m_degree(degree)
    {
        cell_basis          = cell_basis_type(m_degree);
        cell_quadrature     = cell_quadrature_type(2*m_degree);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }

    auto
    compute(const mesh_type& msh, const cell_type& cl,
            const matrix_type& local_mat,
            const vector_type& cell_rhs)
    {
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        auto num_cell_dofs = cell_basis.size();
        auto num_face_dofs = face_basis.size();

        dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

        size_t cell_size = dsr.cell_range().size();
        size_t face_size = dsr.all_faces_range().size();

        matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
        matrix_type K_TF = local_mat.topRightCorner(cell_size, face_size);
        matrix_type K_FT = local_mat.bottomLeftCorner(face_size, cell_size);
        matrix_type K_FF = local_mat.bottomRightCorner(face_size, face_size);

        assert(K_TT.cols() == cell_size);
        assert(K_TT.cols() + K_TF.cols() == local_mat.cols());
        assert(K_TT.rows() + K_FT.rows() == local_mat.rows());
        assert(K_TF.rows() + K_FF.rows() == local_mat.rows());
        assert(K_FT.cols() + K_FF.cols() == local_mat.cols());

        auto K_TT_ldlt = K_TT.llt();
        matrix_type AL = K_TT_ldlt.solve(K_TF);
        vector_type bL = K_TT_ldlt.solve(cell_rhs);

        matrix_type AC = K_FF - K_FT * AL;
        vector_type bC = /* no projection on faces, eqn. 26*/ - K_FT * bL;

        return std::make_pair(AC, bC);
    }

    vector_type
    recover(const mesh_type& msh, const cell_type& cl,
            const matrix_type& local_mat,
            const vector_type& cell_rhs,
            const vector_type& solF)
    {
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        auto num_cell_dofs = cell_basis.size();
        auto num_face_dofs = face_basis.size();

        dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

        size_t cell_size        = dsr.cell_range().size();
        size_t all_faces_size   = dsr.all_faces_range().size();

        vector_type ret( dsr.total_size() );

        matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
        matrix_type K_TF = local_mat.topRightCorner(cell_size, all_faces_size);

        vector_type solT = K_TT.llt().solve(cell_rhs - K_TF*solF);

        ret.head(cell_size)         = solT;
        ret.tail(all_faces_size)    = solF;

        return ret;
    }

};

template<typename LocalData>
class assembler
{
    typedef typename LocalData::scalar_type     scalar_type;
    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;

    typedef Eigen::SparseMatrix<scalar_type>    sparse_matrix_type;
    typedef Eigen::Triplet<scalar_type>         triplet_type;

    std::vector<triplet_type>                   m_triplets;
    size_t                                      m_num_unknowns;

public:

    sparse_matrix_type      matrix;
    vector_type             rhs;

    assembler()                 = delete;
    assembler(const assembler&) = delete;
    assembler(assembler&&)      = delete;

    assembler(LocalData& ld)
    {
        auto msh        = ld.get_mesh();

        m_num_unknowns = ld.num_face_dofs() * (msh.faces_size() + msh.boundary_faces_size());
        matrix = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
        rhs = vector_type::Zero(m_num_unknowns);
    }

    template<typename LocalContrib>
    void
    assemble(LocalData& ld, const LocalContrib& lc)
    {
        auto msh        = ld.get_mesh();
        auto cl         = ld.get_cell();

        auto fcs = faces(msh, cl);
        std::vector<size_t> l2g(fcs.size() * ld.num_face_dofs());
        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            auto face_offset = face_id * ld.num_face_dofs();

            auto pos = face_i * ld.num_face_dofs();

            for (size_t i = 0; i < ld.num_face_dofs(); i++)
                l2g[pos+i] = face_offset+i;
        }

        assert(lc.first.rows() == lc.first.cols());
        assert(lc.first.rows() == lc.second.size());
        assert(lc.second.size() == l2g.size());

        for (size_t i = 0; i < lc.first.rows(); i++)
        {
            for (size_t j = 0; j < lc.first.cols(); j++)
                m_triplets.push_back( triplet_type( l2g.at(i), l2g.at(j), lc.first(i,j) ) );

            rhs(l2g.at(i)) += lc.second(i);
        }
    }

    template<typename Function>
    void
    impose_boundary_conditions(LocalData& ld, const Function& bc)
    {
        auto degree     = ld.get_degree();
        auto msh        = ld.get_mesh();

        typename LocalData::face_quadrature_type    fq(2*degree);
        typename LocalData::face_basis_type         fb(degree);

        size_t fbs = fb.size();
        size_t face_i = 0;
        for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
        {
            auto bfc = *itor;

            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            auto face_offset = face_id * fbs;
            auto face_offset_lagrange = (msh.faces_size() + face_i) * fbs;

            auto fqd = fq.integrate(msh, bfc);

            matrix_type MFF     = matrix_type::Zero(fbs, fbs);
            vector_type rhs_f   = vector_type::Zero(fbs);

            for (auto& qp : fqd)
            {
                auto f_phi = fb.eval_functions(msh, bfc, qp.point());

                for (size_t i = 0; i < f_phi.size(); i++)
                {
                    for (size_t j = 0; j < f_phi.size(); j++)
                        MFF(i,j) += qp.weight() * f_phi[i] * f_phi[j];

                    rhs_f(i) += qp.weight() * f_phi[i] * bc(qp.point());
                }
            }


            for (size_t i = 0; i < MFF.rows(); i++)
            {
                for (size_t j = 0; j < MFF.cols(); j++)
                {
                    m_triplets.push_back( triplet_type(face_offset + i, face_offset_lagrange + j, MFF(i,j)) );
                    m_triplets.push_back( triplet_type(face_offset_lagrange + j, face_offset + i, MFF(i,j)) );
                }
                rhs(face_offset_lagrange+i) = rhs_f(i);
            }

            face_i++;
        }
    }

    void
    finalize()
    {
        matrix.setFromTriplets(m_triplets.begin(), m_triplets.end());
        m_triplets.clear();
    }
};


template<typename Mesh, //typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType>
class assembler_nopre
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    //typedef CellBasisType                       cell_basis_type;
    //typedef CellQuadType                        cell_quadrature_type;
    typedef FaceBasisType                       face_basis_type;
    typedef FaceQuadType                        face_quadrature_type;

    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;

    //cell_basis_type                             cell_basis;
    //cell_quadrature_type                        cell_quadrature;

    face_basis_type                             face_basis;
    face_quadrature_type                        face_quadrature;

    size_t                                      m_degree;

    typedef Eigen::SparseMatrix<scalar_type>    sparse_matrix_type;
    typedef Eigen::Triplet<scalar_type>         triplet_type;

    std::vector<triplet_type>                   m_triplets;
    size_t                                      m_num_unknowns;

public:

    sparse_matrix_type      matrix;
    vector_type             rhs;

    assembler_nopre()                 = delete;
    assembler_nopre(const assembler_nopre&) = delete;
    assembler_nopre(assembler_nopre&&)      = delete;

    assembler_nopre(const mesh_type& msh, size_t degree)
        : m_degree(degree)
    {
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);

        m_num_unknowns = face_basis.size() * (msh.faces_size() + msh.boundary_faces_size());
        matrix = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
        rhs = vector_type::Zero(m_num_unknowns);
    }

    template<typename LocalContrib>
    void
    assemble(const mesh_type& msh, const cell_type& cl, const LocalContrib& lc)
    {
        auto fcs = faces(msh, cl);
        std::vector<size_t> l2g(fcs.size() * face_basis.size());
        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            auto face_offset = face_id * face_basis.size();

            auto pos = face_i * face_basis.size();

            for (size_t i = 0; i < face_basis.size(); i++)
                l2g[pos+i] = face_offset+i;
        }

        assert(lc.first.rows() == lc.first.cols());
        assert(lc.first.rows() == lc.second.size());
        assert(lc.second.size() == l2g.size());

        //std::cout << lc.second.size() << " " << l2g.size() << std::endl;

        for (size_t i = 0; i < lc.first.rows(); i++)
        {
            for (size_t j = 0; j < lc.first.cols(); j++)
                m_triplets.push_back( triplet_type( l2g.at(i), l2g.at(j), lc.first(i,j) ) );

            rhs(l2g.at(i)) += lc.second(i);
        }
    }

    template<typename Function>
    void
    impose_boundary_conditions(mesh_type& msh, const Function& bc)
    {
        size_t fbs = face_basis.size();
        size_t face_i = 0;
        for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
        {
            auto bfc = *itor;

            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            auto face_offset = face_id * fbs;
            auto face_offset_lagrange = (msh.faces_size() + face_i) * fbs;

            auto fqd = face_quadrature.integrate(msh, bfc);

            matrix_type MFF     = matrix_type::Zero(fbs, fbs);
            vector_type rhs_f   = vector_type::Zero(fbs);

            for (auto& qp : fqd)
            {
                auto f_phi = face_basis.eval_functions(msh, bfc, qp.point());

                for (size_t i = 0; i < f_phi.size(); i++)
                {
                    for (size_t j = 0; j < f_phi.size(); j++)
                        MFF(i,j) += qp.weight() * mm_prod(f_phi[i], f_phi[j]);

                    rhs_f(i) += qp.weight() * mm_prod(f_phi[i], bc(qp.point()));
                }
            }


            for (size_t i = 0; i < MFF.rows(); i++)
            {
                for (size_t j = 0; j < MFF.cols(); j++)
                {
                    m_triplets.push_back( triplet_type(face_offset + i, face_offset_lagrange + j, MFF(i,j)) );
                    m_triplets.push_back( triplet_type(face_offset_lagrange + j, face_offset + i, MFF(i,j)) );
                }
                rhs(face_offset_lagrange+i) = rhs_f(i);
            }

            face_i++;
        }
    }

    void
    finalize()
    {
        matrix.setFromTriplets(m_triplets.begin(), m_triplets.end());
        m_triplets.clear();
    }
};



template<typename LocalData, typename GR>
auto
high_order_reconstruction(LocalData& ld, const GR& gr,
                          const dynamic_vector<double>& v)
{
    typedef typename LocalData::scalar_type     scalar_type;
    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;

    auto msh        = ld.get_mesh();
    auto cl         = ld.get_cell();
    auto degree     = ld.get_degree();
    auto mm         = ld.get_cell_mass_matrix();
    auto cell_basis = ld.get_cell_precomputed_basis_integrated_on_cell();

    // Use eqn. (22) to do high order reconstruction.
    auto cbs = cell_basis.size();
    auto zero_range = cell_basis.range(0, degree);
    auto one_range = cell_basis.range(1, degree+1);

    vector_type P = vector_type::Zero(cbs);
    vector_type vT = v.head(zero_range.size());

    vector_type grad = gr.oper * v;
    P.tail(one_range.size()) = grad;

    matrix_type M1 = take(mm.get(), zero_range, zero_range);//cell_mass_matrix.block(0, 0, basis_k_size, basis_k_size);
    matrix_type M2 = take(mm.get(), zero_range, one_range);//cell_mass_matrix.block(0, 0, basis_k_size, cbs);
    matrix_type R = vT - M1.ldlt().solve(M2*P);

    P.head(zero_range.size()) += R;

    return P;
}

/* Compute the L2 error between the function f and the function approximated
 * by the specified dofs.
 *
 * LocalData interface must have:
 *  - get_cell_mass_matrix()
 *  - get_cell_precomputed_basis_integrated_on_cell()
 */
template<typename LocalData, typename Function>
typename LocalData::scalar_type
compute_L2_error(LocalData& ld, size_t degree, const Function& f,
                 const dynamic_vector<typename LocalData::scalar_type>& dofs)
{
    typedef typename LocalData::scalar_type     scalar_type;

    auto mm_wrapper     = ld.get_cell_mass_matrix();
    auto cell_basis     = ld.get_cell_precomputed_basis_integrated_on_cell();

    auto dof_range      = cell_basis.range(0, degree);

    dynamic_matrix<scalar_type> mass = take(mm_wrapper.get(), dof_range, dof_range);
    dynamic_vector<scalar_type> proj = project(ld, degree, f);

    assert( proj.size() == dof_range.size() );

    dynamic_vector<scalar_type> diff = proj - dofs.head(dof_range.size());

    return diff.dot(mass*diff);
}

} // namespace disk
