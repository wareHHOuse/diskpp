/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
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

// NewtonRaphson_solver

#pragma once

#include <iostream>
#include <sstream>
#include <vector>

#include "../Informations.hpp"
#include "../Parameters.hpp"

#include "diskpp/boundary_conditions/boundary_conditions.hpp"
#include "newton_step.hpp"

#include "diskpp/common/timecounter.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

namespace NLE
{

template<typename MeshType>
class NewtonRaphson_solver_finite_strains
{
    typedef MeshType                            mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef ParamRun<scalar_type>               param_type;
    typedef typename disk::hho_degree_info      hdi_type;

    typedef disk::dynamic_matrix<scalar_type> matrix_type;
    typedef disk::dynamic_vector<scalar_type> vector_type;

    typedef disk::vector_boundary_conditions<mesh_type> bnd_type;

    const mesh_type& m_msh;
    const hdi_type& m_hdi;
    const bnd_type& m_bnd;
    const param_type& m_rp;

    std::vector<vector_type> m_solution_data;
    std::vector<vector_type> m_solution_cells, m_solution_faces;

    bool m_verbose;
    bool m_convergence;

  public:
    NewtonRaphson_solver_finite_strains(const mesh_type&  msh,
                                        const hdi_type&   hdi,
                                        const bnd_type&   bnd,
                                        const param_type& rp) :
      m_msh(msh),
      m_hdi(hdi), m_bnd(bnd), m_rp(rp), m_verbose(rp.m_verbose), m_convergence(false)
    {
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

    void
    initialize(const std::vector<vector_type>& initial_solution_cells,
               const std::vector<vector_type>& initial_solution_faces,
               const std::vector<vector_type>& initial_solution)
    {
        m_solution_cells.clear();
        m_solution_cells = initial_solution_cells;
        assert(m_msh.cells_size() == m_solution_cells.size());

        m_solution_faces.clear();
        m_solution_faces = initial_solution_faces;
        assert(m_msh.faces_size() == m_solution_faces.size());

        m_solution_data.clear();
        m_solution_data = initial_solution;
        assert(m_msh.cells_size() == m_solution_data.size());
    }

    template<typename LoadIncrement, typename Law>
    NewtonSolverInfo
    compute(const LoadIncrement&               lf,
            const std::vector<matrix_type>& gradient_precomputed,
            const std::vector<matrix_type>& stab_precomputed,
            Law&                               law)
    {
        NewtonSolverInfo ni;
        timecounter      tc;
        tc.tic();

        // initialise the NewtonRaphson_step
        NewtonRaphson_step_finite_strains<mesh_type> newton_step(m_msh, m_hdi, m_bnd, m_rp);

        newton_step.initialize(m_solution_cells, m_solution_faces, m_solution_data);

        m_convergence           = false;
        bool elastic_prediction = false;

        for (int iter = 0; iter < m_rp.m_iter_max; iter++)
        {
            // assemble lhs and rhs
            AssemblyInfo assembly_info;
            try
            {
                assembly_info =
                  newton_step.assemble(lf, gradient_precomputed, stab_precomputed, law, elastic_prediction);
            }
            catch (const std::invalid_argument& ia)
            {
                std::cerr << "Invalid argument: " << ia.what() << std::endl;
                m_convergence = false;
                tc.toc();
                ni.m_time_newton = tc.elapsed();
                return ni;
            }

            elastic_prediction = false;

            ni.updateAssemblyInfo(assembly_info);
            // test convergence
            m_convergence = newton_step.test_convergence(iter);
            if (m_convergence)
            {
                newton_step.save_solutions(m_solution_cells, m_solution_faces, m_solution_data);
                tc.toc();
                ni.m_time_newton = tc.elapsed();
                return ni;
            }

            // solve the global system
            SolveInfo solve_info = newton_step.solve();
            ni.updateSolveInfo(solve_info);
            // update unknowns
            ni.m_assembly_info.m_time_postpro += newton_step.postprocess();

            ni.m_iter++;
        }

        tc.toc();
        ni.m_time_newton = tc.elapsed();
        return ni;
    }

    bool
    test_convergence() const
    {
        return m_convergence;
    }

    void
    save_solutions(std::vector<vector_type>& solution_cells,
                   std::vector<vector_type>& solution_faces,
                   std::vector<vector_type>& solution)
    {
        solution_cells.clear();
        solution_cells = m_solution_cells;
        assert(m_solution_cells.size() == solution_cells.size());

        solution_faces.clear();
        solution_faces = m_solution_faces;
        assert(m_solution_faces.size() == solution_faces.size());

        solution.clear();
        solution = m_solution_data;
        assert(m_solution_data.size() == solution.size());
    }
};

} // end NLE
