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

#include <iomanip>
#include <map>
#include "common/eigen.hpp"
#include "quadratures/quadratures.hpp"
#include "bases/bases.hpp"
#include "geometry/geometry.hpp"

#include "hho/hho.hpp"

namespace disk {




template<typename Mesh, typename CellQuadrature, typename CellBasis,
         typename FaceQuadrature, typename FaceBasis>
class diffusion_local_data
{
    size_t                                                  m_degree;

public:
    /* MANDATORY PUBLIC TYPES BEGIN */
    typedef Mesh                                            mesh_type;
    typedef typename mesh_type::cell                        cell_type;
    typedef typename mesh_type::face                        face_type;
    typedef typename mesh_type::point_type                  point_type;
    typedef typename mesh_type::scalar_type                 scalar_type;
    //typedef typename dynamic_matrix<scalar_type>            matrix_type;
    //typedef typename dynamic_vector<scalar_type>            vector_type;

    typedef CellQuadrature                                  cell_quadrature_type;
    typedef CellBasis                                       cell_basis_type;

    typedef FaceQuadrature                                  face_quadrature_type;
    typedef FaceBasis                                       face_basis_type;

    typedef typename cell_basis_type::function_value_type   solution_value_type;
    typedef typename cell_basis_type::gradient_value_type   solution_gradient_type;
    /* MANDATORY PUBLIC TYPES END */

    typedef precomputed_basis_functions<mesh_type, cell_basis_type>     cpbf_type;
    typedef precomputed_basis_gradients<mesh_type, cell_basis_type>     cpbg_type;
    typedef precomputed_basis_functions<mesh_type, face_basis_type>     fpbf_type;

    /* MANDATORY PUBLIC TYPES BEGIN */
    typedef local_matrix<cpbf_type, cpbf_type>              cell_mass_matrix_type;
    typedef local_matrix<cpbg_type, cpbg_type>              cell_stif_matrix_type;
    typedef local_rhs<cpbf_type>                            cell_rhs_type;
    typedef local_matrix<fpbf_type, fpbf_type>              face_mass_matrix_type;
    typedef local_matrix<cpbf_type, fpbf_type>              face_trace_matrix_type;
    typedef local_rhs<fpbf_type>                            face_rhs_type;
    /* MANDATORY PUBLIC TYPES END */

private:
    mesh_type                               m_msh;
    cpbf_type                               cpbf_cell;
    cpbg_type                               cpbg_cell;

    cell_quadrature_type                    cell_quadrature;
    face_quadrature_type                    face_quadrature;

    std::shared_ptr<cell_basis_type>        cell_basis;
    std::shared_ptr<face_basis_type>        face_basis;

    ///std::shared_ptr<matrix_type>            cell_mass_matrix_ptr;
    cell_mass_matrix_type                   cell_mass_matrix;
    cell_stif_matrix_type                   cell_stif_matrix;
    cell_rhs_type                           cell_rhs;

    std::vector<fpbf_type>                  fpbf_faces;
    std::vector<face_mass_matrix_type>      face_mass_matrices;
    std::vector<face_trace_matrix_type>     face_trace_matrices;

    std::vector<cpbf_type>                  cpbf_faces;
    std::vector<cpbg_type>                  cpbg_faces;

    cell_type                               m_current_cell;

    std::map<std::string, double>           timings;


public:
    diffusion_local_data(const mesh_type& msh, size_t degree)
        : m_degree(degree), m_msh(msh)
    {
        /* Init objects. This can be slow, is called only once per thread */
        cell_basis          = std::make_shared<cell_basis_type>(m_degree+1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        cpbf_cell           = cpbf_type(cell_basis);
        cpbg_cell           = cpbg_type(cell_basis);

        face_basis          = std::make_shared<face_basis_type>(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }

    template<typename Function>
    void
    compute(const cell_type& cl, const Function& f)
    {
        //timecounter tc;
        /* Compute stuff on the cell. This should run as fast as possible */
        m_current_cell = cl;

        /* Precompute basis functions on the cell */
        //tc.tic();
        auto cell_quad_points = cell_quadrature.integrate(m_msh, cl);
        //tc.toc();
        //timings["int"] += tc.to_double();

        //tc.tic();
        cpbf_cell.compute(m_msh, cl, cell_quad_points); // functions
        //tc.toc();
        //timings["pf"] += tc.to_double();

        //tc.tic();
        cpbg_cell.compute(m_msh, cl, cell_quad_points); // gradients
        //tc.toc();
        //timings["pg"] += tc.to_double();

        //tc.tic();
        cell_mass_matrix.compute(cpbf_cell, cpbf_cell); // build mass matrix
        //tc.toc();
        //timings["mm"] += tc.to_double();

        //tc.tic();
        cell_stif_matrix.compute(cpbg_cell, cpbg_cell); // build stif matrix
        //tc.toc();
        //timings["sm"] += tc.to_double();

        //tc.tic();
        cell_rhs.compute(cpbf_cell, f);                 // build rhs
        //tc.toc();
        //timings["rhs"] += tc.to_double();


        /* Now precompute basis functions on the faces, then build matrices */
        //tc.tic();
        auto fcs = faces(m_msh, cl);
        auto num_faces = fcs.size();

        bool need_resize = false;

        if ( fpbf_faces.size() != num_faces )
        {
            fpbf_faces.resize(fcs.size());
            face_mass_matrices.resize(fcs.size());
            cpbf_faces.resize(fcs.size());
            cpbg_faces.resize(fcs.size());
            face_trace_matrices.resize(fcs.size());
            need_resize = true;
        }

        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto fc = fcs[i];

            //tc.tic();
            auto face_quad_points = face_quadrature.integrate(m_msh, fc);
            //tc.toc();
            //timings["f_int"] += tc.to_double();

            //tc.tic();
            if (need_resize)
            {
                fpbf_faces[i] = fpbf_type(face_basis);
                cpbf_faces[i] = cpbf_type(cell_basis);
                cpbg_faces[i] = cpbg_type(cell_basis);
            }
            //tc.toc();
            //timings["f_alloc"] += tc.to_double();

            //tc.tic();
            fpbf_faces[i].compute(m_msh, fc, face_quad_points);
            cpbf_faces[i].compute(m_msh, cl, face_quad_points);
            cpbg_faces[i].compute(m_msh, cl, face_quad_points);
            //tc.toc();
            //timings["f_precomp"] += tc.to_double();

            //tc.tic();
            face_mass_matrices[i].compute(fpbf_faces[i], fpbf_faces[i]);
            face_trace_matrices[i].compute(cpbf_faces[i], fpbf_faces[i]);
            //tc.toc();
            //timings["f_mats"] += tc.to_double();
        }
        //tc.toc();
        //timings["facestuff"] += tc.to_double();
    }

    void dumptimings(void) const {
        for (auto& t : timings)
            std::cout << t.first << ": " << t.second << std::endl;
    }

    // return the mesh
    auto get_mesh() const {
        return m_msh;
    }

    auto get_cell() const {
        return m_current_cell;
    }

    size_t num_cell_dofs() const {
        return cell_basis->range(0, m_degree).size();
    }

    size_t num_face_dofs() const {
        return face_basis->range(0, m_degree).size();
    }

    // return the degree of
    auto get_degree() const {
        return m_degree;
    }

    // should be explained?
    auto get_cell_mass_matrix() {
        return cell_mass_matrix;
    }

    auto get_cell_precomputed_basis_integrated_on_cell() {
        return cpbf_cell;
    }

    auto get_cell_rhs() {
        return cell_rhs;
    }

    auto get_projector_stuff() {
        return std::make_tuple(cell_mass_matrix,
                               cpbf_cell,
                               face_mass_matrices,
                               fpbf_faces);
    }

    auto get_gr_stuff() {
        return std::make_tuple(cell_stif_matrix,
                               cpbg_cell,
                               fpbf_faces,
                               cpbf_faces,
                               cpbg_faces);
    }

    auto get_stab_stuff() {
        return std::make_tuple(cell_mass_matrix,
                               cpbf_cell,
                               fpbf_faces,
                               face_mass_matrices,
                               face_trace_matrices);
    }

    auto get_statcond_stuff() {
        return std::make_tuple(cpbf_cell,
                               fpbf_faces,
                               cell_rhs);
    }

};









template<typename Mesh, typename CellQuadrature, typename CellBasis,
         typename FaceQuadrature, typename FaceBasis>
class diffusion_template
{
    size_t                                                  m_degree;

    typedef Mesh                                            mesh_type;
    typedef typename mesh_type::cell                        cell_type;
    typedef typename mesh_type::face                        face_type;
    typedef typename mesh_type::point_type                  point_type;
    typedef typename mesh_type::scalar_type                 scalar_type;

    typedef CellQuadrature                                  cell_quadrature_type;
    typedef  std::vector<typename cell_quadrature_type::quadpoint_type>    cell_quaddata_type;
    typedef CellBasis                                       cell_basis_type;

    typedef FaceQuadrature                                  face_quadrature_type;
    typedef std::vector<typename face_quadrature_type::quadpoint_type>    face_quaddata_type;
    typedef FaceBasis                                       face_basis_type;

    typedef dynamic_matrix<scalar_type>                     matrix_type;
    typedef dynamic_vector<scalar_type>                     vector_type;

    cell_basis_type                                         cell_basis;
    cell_quadrature_type                                    cell_quadrature;
    cell_quaddata_type                                      cell_quadrature_points;
    matrix_type                                             cell_mass_matrix;
    matrix_type                                             cell_stiffness_matrix;
    vector_type                                             cell_rhs;

    face_basis_type                                         face_basis;
    face_quadrature_type                                    face_quadrature;

    std::vector<face_quaddata_type>                         face_quadrature_points_vec;
    std::vector<matrix_type>                                face_mass_matrix_vec;
    std::vector<matrix_type>                                face_trace_matrix_vec;
    std::vector<vector_type>                                face_rhs_vec;

    matrix_type GT, A, S;

    void
    make_bases_and_quadratures()
    {
        cell_quadrature  = cell_quadrature_type(2*(m_degree+1));
        cell_basis       = cell_basis_type(m_degree+1);

        face_quadrature  = face_quadrature_type(2*m_degree);
        face_basis       = face_basis_type(m_degree);
    }

    void
    prepare_common_data(const mesh_type& msh, const cell_type& cl)
    {
        cell_quadrature_points = cell_quadrature.integrate(msh, cl);

        auto fcs = faces(msh, cl);
        size_t num_faces = fcs.size();

        face_quadrature_points_vec.resize(num_faces);
        face_mass_matrix_vec.resize(num_faces);
        face_trace_matrix_vec.resize(num_faces);
        face_rhs_vec.resize(num_faces);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];
            face_quadrature_points_vec[face_i] = face_quadrature.integrate(msh, fc);
        }

    }

    template<typename Function>
    void
    build_cell(const mesh_type& msh, const cell_type& cl, const Function& load_function)
    {
        size_t cbs            = cell_basis.size();

        cell_mass_matrix      = matrix_type::Zero(cbs, cbs);
        cell_stiffness_matrix = matrix_type::Zero(cbs, cbs);
        cell_rhs              = vector_type::Zero(cbs);

        for (auto& qp : cell_quadrature_points)
        {
            auto phi  = cell_basis.eval_functions(msh, cl, qp.point());
            auto dphi = cell_basis.eval_gradients(msh, cl, qp.point());

            assert(cbs = phi.size());
            assert(phi.size() == dphi.size());

            for (size_t i = 0; i < cbs; i++)
                for (size_t j = 0; j < cbs; j++)
                    cell_mass_matrix(i,j) += qp.weight() * mm_prod(phi[i], phi[j]);

            for (size_t i = 0; i < cbs; i++)
                for (size_t j = 0; j < cbs; j++)
                    cell_stiffness_matrix(i,j) += qp.weight() * mm_prod(dphi[i], dphi[j]);;

            for (size_t i = 0; i < cbs; i++)
                cell_rhs(i) += qp.weight() * mm_prod(phi[i], load_function(qp.point()));
        }

        return;
    }

    template<typename Function>
    void
    build_faces(const mesh_type& msh, const cell_type& cl, const Function& load_function)
    {
        size_t cbs       = cell_basis.size();
        size_t fbs       = face_basis.size();

        auto fcs = faces(msh, cl);
        size_t num_faces = fcs.size();

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];

            face_mass_matrix_vec[face_i] = matrix_type::Zero(fbs, fbs);
            face_trace_matrix_vec[face_i] = matrix_type::Zero(fbs, cbs);
            face_rhs_vec[face_i] = vector_type::Zero(fbs);

            for (auto& qp : face_quadrature_points_vec[face_i])
            {
                auto f_phi = face_basis.eval_functions(msh, fc, qp.point());
                auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());

                for (size_t i = 0; i < fbs; i++)
                {
                    for (size_t j = 0; j < fbs; j++)
                    {
                        scalar_type mass_matrix_contrib = mm_prod(f_phi[i], f_phi[j]);
                        face_mass_matrix_vec[face_i](i,j) += qp.weight() * mass_matrix_contrib;
                    } //for j

                    for (size_t j = 0; j < c_phi.size(); j++)
                        face_trace_matrix_vec[face_i](i,j) += qp.weight() * mm_prod(f_phi[i], c_phi[j]);

                    face_rhs_vec[face_i](i) += qp.weight() * mm_prod(f_phi[i], load_function(qp.point()));
                } // for i
            } // for qp
        } // for face_i
    }

    void
    build_gradient_reconstruction(const mesh_type& msh, const cell_type& cl)
    {
        //auto cbs = cell_basis.size();
        auto fbs = face_basis.size();
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        /* LHS: take basis functions derivatives from degree 1 to K+1 */
        auto MG_rowcol_range = cell_basis.range(1, m_degree+1);
        matrix_type MG = take(cell_stiffness_matrix, MG_rowcol_range, MG_rowcol_range);

        /* RHS, volumetric part. */
        auto BG_row_range = cell_basis.range(1, m_degree+1);
        auto BG_col_range = cell_basis.range(0, m_degree);

        dofspace_ranges dsr(BG_col_range.size(), fbs, num_faces);

        matrix_type BG = matrix_type::Zero(BG_row_range.size(), dsr.total_size());


        BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
            take(cell_stiffness_matrix, BG_row_range, BG_col_range);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto current_face_range = dsr.face_range(face_i);

            auto fc = fcs[face_i];
            auto n = normal(msh, cl, fc);
            auto face_quadpoints = face_quadrature_points_vec[face_i];

            for (auto& qp : face_quadpoints)
            {
                auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());
                auto c_dphi = cell_basis.eval_gradients(msh, cl, qp.point());

                for (size_t i = 0; i < BG_row_range.size(); i++)
                {
                    for (size_t j = 0; j < BG_col_range.size(); j++)
                    {
                        auto p1 = mm_prod(c_dphi.at(i+BG_row_range.min()), n);
                        scalar_type p2 = mm_prod(c_phi.at(j+BG_col_range.min()), p1);
                        BG(i,j) -= qp.weight() * p2;
                    }
                }

                auto f_phi = face_basis.eval_functions(msh, fc, qp.point());

                for (size_t i = 0; i < BG_row_range.size(); i++)
                {
                    for (size_t j = 0; j < current_face_range.size(); j++)
                    {
                        auto p1 = mm_prod(c_dphi.at(i+BG_row_range.min()), n);
                        scalar_type p2 = mm_prod(f_phi.at(j), p1);
                        BG(i,j+current_face_range.min()) += qp.weight() * p2;
                    }
                }
            }
        }

        GT = MG.ldlt().solve(BG);
        A = BG.transpose() * GT;
    }

    void
    build_stabilization(const mesh_type& msh, const cell_type& cl)
    {
        //size of basis of degree K == degree_index(K+1)
        auto basis_k_size = cell_basis.degree_index(m_degree+1);
        auto deg1_ofs = cell_basis.degree_index(1);
        auto cbs = cell_basis.size();
        auto fbs = face_basis.size();
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        S = matrix_type::Zero(basis_k_size+num_faces*fbs, basis_k_size+num_faces*fbs);

        // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

        //Step 1: compute \pi_T^k p_T^k v (third term).
        matrix_type M1 = cell_mass_matrix.block(0, 0, basis_k_size, basis_k_size);
        matrix_type M2 = cell_mass_matrix.block(0, deg1_ofs, basis_k_size, cbs-deg1_ofs);

        matrix_type proj1 = -M1.ldlt().solve(M2*GT);

        //Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
        matrix_type I_T = matrix_type::Identity(basis_k_size, basis_k_size);
        proj1.block(0, 0, basis_k_size, basis_k_size) += I_T;

        // Step 3: project on faces (eqn. 21)
        for (size_t i = 0; i < num_faces; i++)
        {
            auto h = diameter(msh, fcs[i]);
            Eigen::LDLT<matrix_type> piKF;
            piKF.compute(face_mass_matrix_vec[i]);
            // Step 3a: \pi_F^k( v_F - p_T^k v )
            matrix_type MR1 = face_trace_matrix_vec[i].block(0, deg1_ofs, fbs, cbs-deg1_ofs);
            matrix_type proj2 = piKF.solve(MR1*GT);
            matrix_type I_F = matrix_type::Identity(fbs, fbs);
            auto block_offset = basis_k_size + i * fbs;
            proj2.block(0, block_offset, fbs, fbs) -= I_F;

            // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
            matrix_type MR2 = face_trace_matrix_vec[i].topLeftCorner(fbs, basis_k_size);
            matrix_type proj3 = piKF.solve(MR2*proj1);

            matrix_type BRF = proj2 + proj3;

            matrix_type ss = BRF.transpose() * face_mass_matrix_vec[i] * BRF / h;

            S += ss;
        }
    }

    std::pair<matrix_type, vector_type>
    static_condensation(const mesh_type& msh, const cell_type& cl)
    {
        size_t num_faces = faces(msh, cl).size();

        matrix_type local_mat = A + S;

        auto basis_k_size = cell_basis.degree_index(m_degree+1);
        auto fbs = face_basis.size();

        matrix_type K_TT = local_mat.topLeftCorner(basis_k_size, basis_k_size);
        matrix_type K_TF = local_mat.topRightCorner(basis_k_size, num_faces*fbs);
        matrix_type K_FT = local_mat.bottomLeftCorner(num_faces*fbs, basis_k_size);
        matrix_type K_FF = local_mat.bottomRightCorner(num_faces*fbs, num_faces*fbs);

        assert(K_TT.cols() == basis_k_size);
        assert(K_TT.cols() + K_TF.cols() == local_mat.cols());
        assert(K_TT.rows() + K_FT.rows() == local_mat.rows());
        assert(K_TF.rows() + K_FF.rows() == local_mat.rows());
        assert(K_FT.cols() + K_FF.cols() == local_mat.cols());

        auto K_TT_ldlt = K_TT.ldlt();
        matrix_type AL = K_TT_ldlt.solve(K_TF);
        vector_type bL = K_TT_ldlt.solve(cell_rhs.block(0, 0, basis_k_size, 1));

        matrix_type AC = K_FF - K_FT * AL;
        vector_type bC = - K_FT * bL;

        return std::make_pair(AC, bC);
    }

    template<typename Function>
    vector_type
    project(const mesh_type& msh, const cell_type& cl, const Function& f)
    {
        auto cbrange = cell_basis.range(0, m_degree);
        auto fbrange = face_basis.range(0, m_degree);
        size_t cell_basis_k_size = cbrange.size();
        size_t face_basis_k_size = fbrange.size();

        std::cout << cell_basis_k_size << std::endl;

        matrix_type cell_mm = cell_mass_matrix.block(0, 0, cell_basis_k_size, cell_basis_k_size);
        vector_type cell_r = cell_rhs.block(0, 0, cell_basis_k_size, 1);
        vector_type vT = cell_mm.ldlt().solve(cell_r);

        std::vector<vector_type> vF;
        size_t num_faces = face_mass_matrix_vec.size();
        vF.resize( num_faces );
        for (size_t i = 0; i < num_faces; i++)
        {
            matrix_type face_mm = face_mass_matrix_vec[i].block(0, 0, face_basis_k_size, face_basis_k_size);
            vector_type face_r = face_rhs_vec[i].block(0, 0, face_basis_k_size, 1);
            vF[i] = face_mm.ldlt().solve(face_r);
        }

        vector_type ret;
        ret.resize(cell_basis_k_size + num_faces*face_basis_k_size);
        ret.block(0, 0, cell_basis_k_size, 1) = vT;
        for (size_t i = 0; i < num_faces; i++)
        {
            size_t start    = cell_basis_k_size + i*face_basis_k_size;
            size_t stop     = cell_basis_k_size + (i+1)*face_basis_k_size;
            ret.block(start, 0, stop-start, 1) = vF[i];

        }

        return ret;
    }

public:

    diffusion_template()
        : m_degree(1)
    {
        make_bases_and_quadratures();
    }

    diffusion_template(size_t degree)
        : m_degree(degree)
    {
        make_bases_and_quadratures();
    }

    template<typename Function>
    scalar_type
    compute_cell_error(const mesh_type& msh, const cell_type& cl, const vector_type& xT, const Function& sol)
    {
        auto cbs = cell_basis.size();

        vector_type eT = vector_type::Zero(cbs);

        for (auto& qp : cell_quadrature_points)
        {
            auto phi  = cell_basis.eval_functions(msh, cl, qp.point());

            for (size_t i = 0; i < cbs; i++)
                eT(i) += qp.weight() * mm_prod(phi[i], sol(qp.point()));
        }

        vector_type dT = cell_mass_matrix.ldlt().solve(eT);

        vector_type diff = dT - xT;

        scalar_type e = diff.dot(cell_mass_matrix * diff);

        return e;
    }

    bool
    test_operators(const mesh_type& msh, const cell_type& cl)
    {
        auto f = [](const point_type& pt) -> scalar_type {
            return pt.x() * pt.y();
        };

        auto gradf = [](const point_type& pt) -> static_vector<scalar_type, mesh_type::dimension> {
            static_vector<scalar_type, mesh_type::dimension> grad = static_vector<scalar_type, mesh_type::dimension>::Zero();
            grad(0) = pt.y();
            grad(1) = pt.x();
            //grad(2) = 0.;
            return grad;
        };

        prepare_common_data(msh, cl);
        build_cell(msh, cl, f);
        build_faces(msh, cl, f);
        build_gradient_reconstruction(msh, cl);

        size_t cell_basis_k_size = cell_basis.degree_index(m_degree+1);

        vector_type v = project(msh, cl, f);

        vector_type rgrad = GT*v;

        auto tps = make_test_points(msh, cl);

        for (auto& tp : tps)
        {
            scalar_type expfun = f(tp);
            vector_type expgrad = gradf(tp);
            scalar_type compfun = 0.;
            static_vector<scalar_type, mesh_type::dimension> compgrad = static_vector<scalar_type, mesh_type::dimension>::Zero();

            auto phi = cell_basis.eval_functions(msh, cl, tp);
            for (size_t i = 0; i < cell_basis_k_size; i++)
                compfun += v(i) * phi[i];

            auto dphi = cell_basis.eval_gradients(msh, cl, tp);
            for (size_t i = 0; i < dphi.size()-1; i++)
                compgrad += rgrad(i) * dphi[i+1];
                std::cout << "comp" << std::endl;
            std::cout << std::setprecision(4);

            std::cout << expfun - compfun << std::endl;

            std::cout << (expgrad - compgrad).transpose() << std::endl;
            std::cout << expgrad.transpose() << std::endl;
            std::cout << compgrad.transpose() << std::endl;
        }

        return true;
    }

    template<typename Function>
    std::pair<matrix_type, vector_type>
    build_local_contrib(const mesh_type& msh, const cell_type& cl, const Function& load_function)
    {
        prepare_common_data(msh, cl);
        build_cell(msh, cl, load_function);
        build_faces(msh, cl, load_function);
        build_gradient_reconstruction(msh, cl);
        build_stabilization(msh, cl);
        return static_condensation(msh, cl);
    }

    template<typename Function>
    vector_type
    recover_full_solution(const mesh_type& msh, const cell_type& cl,
                          const vector_type& solF, const Function& load_function)
    {
        prepare_common_data(msh, cl);
        build_cell(msh, cl, load_function);
        build_faces(msh, cl, load_function);
        build_gradient_reconstruction(msh, cl);
        build_stabilization(msh, cl);
        matrix_type local_mat = A + S;

        size_t num_faces = faces(msh, cl).size();
        auto basis_k_size = cell_basis.degree_index(m_degree+1);
        auto fbs = face_basis.size();

        vector_type ret = vector_type::Zero(basis_k_size + num_faces*fbs);

        matrix_type K_TT = local_mat.topLeftCorner(basis_k_size, basis_k_size);
        matrix_type K_TF = local_mat.topRightCorner(basis_k_size, num_faces*fbs);

        vector_type solT = K_TT.ldlt().solve(cell_rhs.block(0, 0, basis_k_size, 1) - K_TF*solF);

        ret.head(basis_k_size) = solT;
        ret.tail(num_faces*fbs) = solF;

        return ret;
    }

    vector_type
    high_order_reconstruction(const mesh_type& msh, const cell_type& cl,
                              const vector_type& v)
    {
        // Use eqn. (22) to do high order reconstruction.
        auto cbs = cell_basis.size();
        auto basis_k_size = cell_basis.degree_index(m_degree+1);
        auto deg1_ofs = cell_basis.degree_index(1);

        vector_type P = vector_type::Zero(cbs);
        vector_type vT = v.head(basis_k_size);

        vector_type grad = GT*v;
        P.tail(cbs-deg1_ofs) = grad;

        matrix_type M1 = cell_mass_matrix.block(0, 0, basis_k_size, basis_k_size);
        matrix_type M2 = cell_mass_matrix.block(0, 0, basis_k_size, cbs);
        matrix_type R = vT - M1.ldlt().solve(M2*P);

        P.head(basis_k_size) += R;

        return P;
    }


    vector_type
    high_order_reconstruction_2(const mesh_type& msh, const cell_type& cl,
                                const vector_type& v)
    {
        // Use eqn. (22) to do high order reconstruction.
        auto cbs = cell_basis.size();
        auto basis_k_size = cell_basis.degree_index(m_degree+1);
        auto deg1_ofs = cell_basis.degree_index(1);

        vector_type P = vector_type::Zero(cbs);
        vector_type vT = v.head(basis_k_size);

        vector_type grad = GT*v;
        P.tail(cbs-deg1_ofs) = grad;
        P(0) = vT(0);

        return P;
    }

    scalar_type
    evaluate_at_point(const mesh_type& msh, const cell_type& cl,
                      const vector_type& data, const point_type& pt)
    {
        auto phi = cell_basis.eval_functions(msh, cl, pt);
        assert (data.size() == phi.size());

        scalar_type acc = 0.;
        for (size_t i = 0; i < phi.size(); i++)
            acc += data(i) * phi[i];

        return acc;
    }

    size_t cell_basis_size() const
    {
        return cell_basis.size();
    }

    size_t face_basis_size() const
    {
        return face_basis.size();
    }

};

} // namespace disk
