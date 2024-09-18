/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
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

/* Solve the linear elasticity problem with HHO method
 * A hybrid high-order locking-free method for linear elasticity on general meshes
 * Di Pietro, Daniele A. and Ern, Alexandre.
 * Comput. Methods Appl. Mech. Engrg. (2015)
 *  DOI: 10.1016/j.cma.2014.09.009
 */

#include <iostream>
#include <sstream>

#include "diskpp/methods/hho"
#include "diskpp/solvers/solver.hpp"
#include "diskpp/output/gmshConvertMesh.hpp"
#include "diskpp/output/gmshDisk.hpp"
#include "diskpp/output/postMesh.hpp"
#include "diskpp/common/timecounter.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

struct assembly_info
{
    size_t linear_system_size;
    double time_gradrec, time_statcond, time_stab, time_assembly, time_divrec;
};

struct solver_info
{
    double time_solver;
};

struct postprocess_info
{
    double time_postprocess;
};

struct ElasticityParameters
{
    double lambda;
    double mu;
};

template<typename Mesh>
class linear_elasticity_solver
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    typedef disk::vector_boundary_conditions<mesh_type> bnd_type;

    typedef disk::dynamic_matrix<scalar_type> matrix_dynamic;
    typedef disk::dynamic_vector<scalar_type> vector_dynamic;

    typedef disk::assembler_mechanics<mesh_type> assembler_type;

    size_t m_cell_degree, m_face_degree;

    const static size_t dimension = mesh_type::dimension;

    const bnd_type&                      m_bnd;
    const mesh_type&                     m_msh;
    typename disk::hho_degree_info m_hdi;
    assembler_type                       m_assembler;

    std::vector<vector_dynamic> m_solution_data;
    vector_dynamic              m_system_solution;

    std::vector<vector_dynamic> m_bL;
    std::vector<matrix_dynamic> m_AL;

    bool m_verbose;

    ElasticityParameters m_elas_parameters;

  public:
    linear_elasticity_solver(const mesh_type&           msh,
                             const bnd_type&            bnd,
                             const ElasticityParameters data,
                             int                        degree,
                             int                        l = 0) :
      m_msh(msh),
      m_verbose(false), m_bnd(bnd)
    {
        if (l < -1 or l > 1)
        {
            std::cout << "'l' should be -1, 0 or 1. Reverting to 0." << std::endl;
            l = 0;
        }

        if (degree == 0 && l == -1)
        {
            std::cout << "'l' should be 0 or 1. Reverting to 0." << std::endl;
            l = 0;
        }

        m_cell_degree = degree + l;
        m_face_degree = degree;

        m_elas_parameters.mu     = data.mu;
        m_elas_parameters.lambda = data.lambda;

        m_hdi       = disk::hho_degree_info(m_cell_degree, m_face_degree);
        m_assembler = disk::make_mechanics_assembler(m_msh, m_hdi, m_bnd);
        m_AL.clear();
        m_AL.reserve(m_msh.cells_size());

        m_bL.clear();
        m_bL.reserve(m_msh.cells_size());

        if (m_verbose)
            m_hdi.info_degree();
    }

    void
    changeElasticityParameters(const ElasticityParameters data)
    {
        m_elas_parameters.mu     = data.mu;
        m_elas_parameters.lambda = data.lambda;
    }

    bool
    verbose(void) const
    {
        return m_verbose;
    }
    void
    verbose(bool v)
    {
        m_verbose = v;
    }

    size_t
    getDofs()
    {
        return m_msh.faces_size() * disk::vector_basis_size(m_hdi.face_degree(), dimension - 1, dimension);
    }

    template<typename LoadFunction>
    assembly_info
    assemble(const LoadFunction& lf)
    {
        assembly_info ai;
        bzero(&ai, sizeof(ai));

        timecounter tc;

        for (auto& cl : m_msh)
        {
            tc.tic();
            const auto sgr = disk::make_vector_hho_symmetric_laplacian(m_msh, cl, m_hdi);
            const auto sg  = disk::make_matrix_hho_symmetric_gradrec(m_msh, cl, m_hdi);
            tc.toc();
            ai.time_gradrec += tc.elapsed();

            tc.tic();
            const auto dr = make_hho_divergence_reconstruction(m_msh, cl, m_hdi);
            tc.toc();
            ai.time_divrec += tc.elapsed();

            tc.tic();
            matrix_dynamic stab;
            if (m_hdi.cell_degree() == (m_hdi.face_degree() + 1))
            {
                stab = make_vector_hdg_stabilization(m_msh, cl, m_hdi);
            }
            else
            {
                stab = make_vector_hho_stabilization(m_msh, cl, sgr.first, m_hdi);
            }
            tc.toc();
            ai.time_stab += tc.elapsed();

            tc.tic();
            auto                 cb       = disk::make_vector_monomial_basis(m_msh, cl, m_hdi.cell_degree());
            const auto           cell_rhs = make_rhs(m_msh, cl, cb, lf, 1);
            const matrix_dynamic loc =
              2.0 * m_elas_parameters.mu * (sg.second + stab) + m_elas_parameters.lambda * dr.second;
            const auto scnp = make_vector_static_condensation_withMatrix(m_msh, cl, m_hdi, loc, cell_rhs);

            m_AL.push_back(std::get<1>(scnp));
            m_bL.push_back(std::get<2>(scnp));
            tc.toc();
            ai.time_statcond += tc.elapsed();

            m_assembler.assemble(m_msh, cl, m_bnd, std::get<0>(scnp), 2);
        }

        m_assembler.impose_neumann_boundary_conditions(m_msh, m_bnd);
        m_assembler.finalize();

        ai.linear_system_size = m_assembler.LHS.rows();
        ai.time_assembly      = ai.time_gradrec + ai.time_divrec + ai.time_stab + ai.time_statcond;
        return ai;
    }

    solver_info
    solve(void)
    {
        solver_info si;

        size_t systsz = m_assembler.LHS.rows();
        size_t nnz    = m_assembler.LHS.nonZeros();

        if (verbose())
        {
            std::cout << "Starting linear solver..." << std::endl;
            std::cout << " * Solving for " << systsz << " unknowns." << std::endl;
            std::cout << " * Matrix fill: " << 100.0 * double(nnz) / (systsz * systsz) << "%" << std::endl;
        }

        timecounter tc;

        tc.tic();
        m_system_solution = vector_dynamic::Zero(systsz);

        disk::solvers::pardiso_params<scalar_type> pparams;
        mkl_pardiso(pparams, m_assembler.LHS, m_assembler.RHS, m_system_solution);

        tc.toc();
        si.time_solver = tc.elapsed();

        return si;
    }

    template<typename LoadFunction>
    postprocess_info
    postprocess(const LoadFunction& lf)
    {
        const auto   fbs = disk::vector_basis_size(m_hdi.face_degree(), dimension - 1, dimension);
        const auto   cbs = disk::vector_basis_size(m_hdi.cell_degree(), dimension, dimension);

        postprocess_info pi;

        m_solution_data.reserve(m_msh.cells_size());

        const auto solF = m_assembler.expand_solution(m_msh, m_bnd, m_system_solution, 2);

        timecounter tc;
        tc.tic();
        size_t cell_i = 0;
        for (auto& cl : m_msh)
        {
            const auto fcs        = faces(m_msh, cl);
            const auto num_faces  = fcs.size();
            const auto total_dofs = cbs + num_faces * fbs;

            vector_dynamic xFs = vector_dynamic::Zero(num_faces * fbs);
            vector_dynamic x   = vector_dynamic::Zero(total_dofs);

            for (size_t face_i = 0; face_i < num_faces; face_i++)
            {
                const auto fc  = fcs[face_i];
                const auto eid = find_element_id(m_msh.faces_begin(), m_msh.faces_end(), fc);
                if (!eid.first)
                    throw std::invalid_argument("This is a bug: face not found");

                const auto face_id             = eid.second;
                xFs.segment(face_i * fbs, fbs) = solF.segment(face_id * fbs, fbs);
            }

            const vector_dynamic xT          = m_bL.at(cell_i) - m_AL.at(cell_i) * xFs;
            x.segment(0, cbs)                = xT;
            x.segment(cbs, total_dofs - cbs) = xFs;
            m_solution_data.push_back(x);

            cell_i++;
        }
        tc.toc();
        pi.time_postprocess = tc.elapsed();

        return pi;
    }

    template<typename AnalyticalSolution>
    scalar_type
    compute_l2_displacement_error(const AnalyticalSolution& as)
    {
        scalar_type err_dof = 0;

        const size_t cbs      = disk::vector_basis_size(m_hdi.cell_degree(), dimension, dimension);
        const int    diff_deg = m_hdi.face_degree() - m_hdi.cell_degree();
        const int    di       = std::max(diff_deg, 1);

        size_t cell_i = 0;

        for (auto& cl : m_msh)
        {
            const auto x = m_solution_data.at(cell_i++);

            const vector_dynamic true_dof = disk::project_function(m_msh, cl, m_hdi.cell_degree(), as, di);

            auto                 cb   = disk::make_vector_monomial_basis(m_msh, cl, m_hdi.cell_degree());
            const matrix_dynamic mass = disk::make_mass_matrix(m_msh, cl, cb);

            const vector_dynamic comp_dof = x.head(cbs);
            const vector_dynamic diff_dof = (true_dof - comp_dof);
            assert(comp_dof.size() == true_dof.size());
            err_dof += diff_dof.dot(mass * diff_dof);
        }

        return sqrt(err_dof);
    }

    template<typename AnalyticalSolution>
    scalar_type
    compute_l2_stress_error(const AnalyticalSolution& stress) const
    {
        const auto face_degree = m_hdi.face_degree();
        const auto cell_degree = m_hdi.cell_degree();
        const auto rec_degree  = m_hdi.reconstruction_degree();

        size_t      cell_i = 0;
        scalar_type error_stress(0.0);

        for (auto& cl : m_msh)
        {
            const auto           x   = m_solution_data.at(cell_i++);
            const auto           sgr = make_vector_hho_symmetric_laplacian(m_msh, cl, m_hdi);
            const vector_dynamic GTu = sgr.first * x;

            const auto           dr   = make_hho_divergence_reconstruction(m_msh, cl, m_hdi);
            const vector_dynamic divu = dr.first * x;

            auto cbas_v = disk::make_vector_monomial_basis(m_msh, cl, rec_degree);
            auto cbas_s = disk::make_scalar_monomial_basis(m_msh, cl, face_degree);

            auto qps = disk::integrate(m_msh, cl, 2 * rec_degree);
            for (auto& qp : qps)
            {
                const auto gphi   = cbas_v.eval_sgradients(qp.point());
                const auto GT_iqn = disk::eval(GTu, gphi, dimension);

                const auto divphi   = cbas_s.eval_functions(qp.point());
                const auto divu_iqn = disk::eval(divu, divphi);

                const auto sigma =
                  2.0 * m_elas_parameters.mu * GT_iqn +
                  m_elas_parameters.lambda * divu_iqn * disk::static_matrix<scalar_type, dimension, dimension>::Identity();

                const auto stress_diff = (stress(qp.point()) - sigma).eval();

                error_stress += qp.weight() * stress_diff.squaredNorm();
            }
        }

        return sqrt(error_stress);
    }

    void
    compute_continuous_displacement(const std::string& filename) const
    {
        const size_t cell_degree = m_hdi.cell_degree();
        const size_t cbs         = disk::vector_basis_size(cell_degree, dimension, dimension);

        // compute mesh for post-processing
        disk::PostMesh<mesh_type> post_mesh = disk::PostMesh<mesh_type>(m_msh);
        gmsh::Gmesh               gmsh      = disk::convertMesh(post_mesh);
        auto                      storage   = post_mesh.mesh().backend_storage();

        const disk::static_vector<scalar_type, dimension> vzero = disk::static_vector<scalar_type, dimension>::Zero();

        const size_t nb_nodes(gmsh.getNumberofNodes());

        // first(number of data at this node), second(cumulated value)
        std::vector<std::pair<size_t, disk::static_vector<scalar_type, dimension>>> value(nb_nodes, std::make_pair(0, vzero));

        size_t cell_i = 0;
        for (auto& cl : m_msh)
        {
            vector_dynamic x          = m_solution_data.at(cell_i).head(cbs);
            const auto     cell_nodes = post_mesh.nodes_cell(cell_i);
            auto           cbas       = disk::make_vector_monomial_basis(m_msh, cl, cell_degree);

            // Loop on the nodes of the cell
            for (auto& point_id : cell_nodes)
            {
                const auto pt = storage->points[point_id];

                const auto phi = cbas.eval_functions(pt);
                assert(phi.rows() == cbs);
                const auto depl = disk::eval(x, phi);

                // Add displacement at node
                value[point_id].first++;
                value[point_id].second += depl;
            }
            cell_i++;
        }

        std::vector<gmsh::Data>    data;    // create data
        std::vector<gmsh::SubData> subdata; // create subdata
        data.reserve(nb_nodes);             // data has a size of nb_node

        // Compute the average value and save it
        for (size_t i_node = 0; i_node < value.size(); i_node++)
        {
            const disk::static_vector<scalar_type, dimension> depl_avr = value[i_node].second / double(value[i_node].first);

            const gmsh::Data tmp_data(i_node + 1, disk::convertToVectorGmsh(depl_avr));
            data.push_back(tmp_data);
        }

        // Create and init a nodedata view
        gmsh::NodeData nodedata(3, 0.0, "depl_node_cont", data, subdata);
        // Save the view
        nodedata.saveNodeData(filename, gmsh);
    }
};
