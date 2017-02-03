/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016, 2017 - matteo.cicuttin@enpc.fr
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

#include "hho.hpp"

namespace disk {

template<typename Mesh>
class multiscale_local_problem
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;
    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;
    typedef Eigen::SparseMatrix<scalar_type>    sparse_matrix_type;
    typedef Eigen::Triplet<scalar_type>         triplet_type;

    typedef basis_quadrature_data<mesh_type,
                                  scaled_monomial_scalar_basis,
                                  quadrature>   bqdata_type;

    size_t                        m_degree;
    size_t                        m_refinement_levels;

public:
    sparse_matrix_type                          matrix;
    matrix_type                                 rhs;

    multiscale_local_problem()
    {
        m_degree                = 1;
        m_refinement_levels     = 2;
    }

    multiscale_local_problem(size_t degree, size_t refinement_levels)
    {
        m_degree                = degree;
        m_refinement_levels     = refinement_levels;
    }

    void assemble(const mesh_type& coarse_msh,
                  const cell_type& coarse_cl)
    {
        submesher<Mesh>                                 submesher;
        bqdata_type                                     bqd(m_degree, m_degree);
        gradient_reconstruction_bq<bqdata_type>         gradrec(bqd);
        diffusion_like_stabilization_bq<bqdata_type>    stab(bqd);

        auto msh = submesher.generate_mesh(coarse_msh, coarse_cl, m_refinement_levels);

        for (auto& cl : msh)
        {
            std::cout << cl << std::endl;
            auto pts = points(msh, cl);
            for (auto& pt : pts)
                std::cout << " " << pt << std::endl;

            auto fcs = faces(msh, cl);
            for (auto& fc : fcs)
            {
                std::cout << fc << std::endl;
                auto pts = points(msh, fc);
                for (auto& pt : pts)
                    std::cout << " " << pt << std::endl;
            }
        }

        auto basis = make_scaled_monomial_scalar_basis(coarse_msh, m_degree-1, m_degree);
        auto coarse_cell_basis = basis.first;
        auto coarse_face_basis = basis.second;

        auto num_cell_faces = howmany_faces(coarse_msh, coarse_cl);

        auto num_cell_dofs = howmany_dofs(bqd.cell_basis);
        auto num_face_dofs = howmany_dofs(bqd.face_basis);

        auto num_cell_funcs = howmany_dofs(coarse_cell_basis);
        auto num_face_funcs = howmany_dofs(coarse_face_basis);

        std::cout << "Fine scale cell dofs: " << num_cell_dofs << std::endl;
        std::cout << "Fine scale face dofs: " << num_face_dofs << std::endl;
        std::cout << "Coarse scale cell funcs: " << num_cell_funcs << std::endl;
        std::cout << "Coarse scale face funcs: " << num_face_funcs << std::endl;

        auto matrix_face_offset = num_cell_dofs * msh.cells_size();
        auto matrix_mult_offset = matrix_face_offset + num_face_dofs * msh.faces_size();
        auto system_size = num_cell_dofs * msh.cells_size() +
                           num_face_dofs * msh.faces_size() +
                           num_face_funcs * num_cell_faces;

        std::cout << "System size is " << system_size << std::endl;
        std::cout << "  " << num_cell_dofs << " " << msh.cells_size() << std::endl;
        std::cout << "  " << num_face_dofs << " " << msh.faces_size() << std::endl;
        std::cout << "  " << num_face_funcs << " " << num_cell_faces << std::endl;

        matrix = sparse_matrix_type(system_size, system_size);

        auto rhs_size = num_cell_funcs + num_face_funcs * num_cell_faces;
        rhs = matrix_type::Zero(system_size, rhs_size);

        std::vector<triplet_type> triplets;

        /* Assemble standard HHO part */

        size_t cell_idx = 0;
        for (auto& cl : msh)
        {
            std::vector<size_t> l2g(num_cell_dofs + num_cell_faces * num_face_dofs);

            /* Build DOF offset table: cell */
            for (size_t i = 0; i < num_cell_dofs; i++)
                l2g[i] = cell_idx * num_cell_dofs + i;

            /* Build DOF offset table: faces */
            auto fcs = faces(msh, cl);
            for (size_t i = 0; i < fcs.size(); i++)
            {
                auto fc = fcs[i];
                auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
                if (!eid.first)
                    throw std::invalid_argument("This is a bug: face not found");

                auto face_id = eid.second;
                /* global offset of current face */
                auto face_offset = matrix_face_offset + face_id * num_face_dofs;

                /* offset in the DOF table */
                auto dt_ofs = num_cell_dofs + i * num_face_dofs;

                for (size_t j = 0; j < num_face_dofs; j++)
                    l2g[dt_ofs+j] = face_offset+j;
            }

            /* Compute HHO element contribution */
            gradrec.compute(msh, cl);
            stab.compute(msh, cl, gradrec.oper);
            dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
            assert(loc.rows() == l2g.size());
            assert(loc.cols() == l2g.size());

            /* Assemble into the matrix */
            for (size_t i = 0; i < l2g.size(); i++)
                for (size_t j = 0; j < l2g.size(); j++)
                    triplets.push_back( triplet_type(l2g[i], l2g[j], loc(i,j)) );

            /* Now compute the multiscale-specific stuff */
            matrix_type cell_rhs;
            cell_rhs = matrix_type::Zero(num_cell_dofs, num_cell_funcs);

            auto cell_quadpoints = bqd.cell_quadrature.integrate(msh, cl);
            for (auto& qp : cell_quadpoints)
            {
                auto phi = bqd.cell_basis.eval_functions(msh, cl, qp.point());
                auto phi_coarse = coarse_cell_basis.eval_functions(coarse_msh, coarse_cl, qp.point());
                cell_rhs += qp.weight() * phi * phi_coarse.transpose();
            }

            auto cell_offset = cell_idx * num_cell_dofs;
            rhs.block(cell_offset, 0, num_cell_dofs, num_cell_funcs) = cell_rhs;

            cell_idx++;
        }

        auto bid_list = msh.boundary_id_list();

        for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
        {
            auto fc = *itor;
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            /* row */
            auto face_offset = matrix_face_offset + face_id * num_face_dofs;

            size_t coarse_id = msh.boundary_id(fc);
            size_t fine_id;
            for (size_t i = 0; i < bid_list.size(); i++)
                if (bid_list[i] == coarse_id)
                {
                    fine_id = i;
                    break;
                }

            assert (fine_id < 3);

            auto coarse_fc = *(coarse_msh.faces_begin() + coarse_id);

            matrix_type face_matrix, ff;
            face_matrix = matrix_type::Zero(num_face_dofs, num_face_funcs);
            ff = matrix_type::Zero(num_face_funcs, num_face_funcs);

            auto face_quadpoints = bqd.face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
                auto phi = bqd.face_basis.eval_functions(msh, fc, qp.point());
                auto phi_coarse = coarse_face_basis.eval_functions(coarse_msh, coarse_fc, qp.point());
                face_matrix += qp.weight() * phi * phi_coarse.transpose();
            }

            auto coarse_face_quadpoints = bqd.face_quadrature.integrate(coarse_msh, coarse_fc);
            for (auto& qp : coarse_face_quadpoints)
            {
                auto phi_coarse = coarse_face_basis.eval_functions(coarse_msh, coarse_fc, qp.point());
                ff += qp.weight() * phi_coarse * phi_coarse.transpose();
            }

            /* col */
            auto mult_offset = matrix_mult_offset + fine_id * num_face_funcs;
            auto rhs_offset = num_cell_funcs + fine_id*num_face_funcs;

            for (size_t j = 0; j < face_matrix.rows(); j++)
            {
                for (size_t k = 0; k < face_matrix.cols(); k++)
                {
                    size_t row = face_offset + j;
                    size_t col = mult_offset + k;

                    triplets.push_back( triplet_type(row, col, -face_matrix(j,k)) );
                    triplets.push_back( triplet_type(col, row, face_matrix(j,k)) );

                }
            }

            rhs.block(mult_offset, rhs_offset, ff.rows(), ff.cols()) = ff;

        }


        std::ofstream ofsr("rhs.dat");

        for (size_t i = 0; i < rhs.rows(); i++)
        {
            for (size_t j = 0; j < rhs.cols(); j++)
                ofsr << rhs(i,j) << " ";
            ofsr << std::endl;
        }

        ofsr.close();


        matrix.setFromTriplets(triplets.begin(), triplets.end());


        std::stringstream ssw;
        ssw << "matrix.mat";
        std::ofstream ofs(ssw.str());

        for (int k=0; k<matrix.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(matrix,k); it; ++it)
                ofs << it.row() << " " << it.col() << " " << it.value() << std::endl;

        ofs.close();

        triplets.clear();



#ifdef HAVE_INTEL_MKL
        Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
#else
        Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
#endif

        solver.analyzePattern(matrix);
        solver.factorize(matrix);
        matrix_type X = solver.solve(rhs);

        auto eid = find_element_id(coarse_msh.cells_begin(), coarse_msh.cells_end(), coarse_cl);
        if (!eid.first)
            throw std::invalid_argument("This is a bug: cell not found");

        auto cell_id = eid.second;

        for (size_t funcnum = 0; funcnum < num_cell_funcs + 3*num_face_funcs; funcnum++)
        {
            std::stringstream ss;
            ss << "mshho_k_" << m_degree << "_rl_" << m_refinement_levels << "_fn_" << funcnum << ".dat";

            std::ofstream ofs_sol(ss.str());
            cell_idx = 0;
            for (auto& cl : msh)
            {
                auto cell_sol = X.block(num_cell_dofs * cell_idx, funcnum, num_cell_dofs, 1);
                //auto qps = bqd.cell_quadrature.integrate(msh, cl);
                auto qps = make_test_points(msh, cl, 15);
                for (auto& qp : qps)
                {
                    auto phi = bqd.cell_basis.eval_functions(msh, cl, qp);

                    scalar_type pot = 0.0;
                    for (size_t i = 0; i < bqd.cell_basis.range(0, m_degree).size(); i++)
                        pot += phi[i] * cell_sol(i, 0);

                    auto tp = qp;
                    for (size_t i = 0; i < mesh_type::dimension; i++)
                        ofs_sol << tp[i] << " ";
                    ofs_sol << pot << std::endl;
                }

                cell_idx++;
            }

            ofs_sol.close();
        }
    }
};


} //namespace disk
