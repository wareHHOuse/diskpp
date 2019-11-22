/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework.
 * M. Abbas, A. Ern, N. Pignet.
 * International Journal of Numerical Methods in Engineering (2019)
 * 120(3), 303-327
 * DOI: 10.1002/nme.6137
 */

// NewtonRaphson_step

#pragma once

#include <iostream>
#include <sstream>
#include <vector>

#include "NewtonIteration.hpp"
#include "NewtonSolverInformations.hpp"
#include "NewtonSolverParameters.hpp"
#include "mechanics/behaviors/laws/behaviorlaws.hpp"

#include "adaptivity/adaptivity.hpp"
#include "boundary_conditions/boundary_conditions.hpp"
#include "methods/hho"

#include "timecounter.h"

namespace disk
{

namespace mechanics
{

/**
 * @brief Newton-Raphson step for nonlinear solid mechanics
 *
 *  Specialized for HHO methods
 *
 *  Options :  - small and finite deformations
 *             - plasticity, hyperelasticity (various laws)
 *
 * @tparam MeshType type of the mesh
 */
template<typename MeshType>
class NewtonStep
{
    typedef MeshType                            mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    typedef NewtonSolverParameter<scalar_type>    param_type;
    typedef vector_boundary_conditions<mesh_type> bnd_type;
    typedef hho_degree_info                       hdi_type;
    typedef Behavior<mesh_type>                   behavior_type;

    const mesh_type&  m_msh;
    const hdi_type&   m_hdi;
    const bnd_type&   m_bnd;
    const param_type& m_rp;

    std::vector<vector_type> m_solution, m_solution_faces;

    bool m_verbose;
    bool m_convergence;

  public:
    NewtonStep(const mesh_type& msh, const hdi_type& hdi, const bnd_type& bnd, const param_type& rp) :
      m_msh(msh), m_hdi(hdi), m_bnd(bnd), m_rp(rp), m_verbose(rp.m_verbose), m_convergence(false)
    {
        m_solution.clear();
        m_solution_faces.clear();
    }

    /**
     * @brief return a boolean to know if the verbosity mode is activated
     *
     */
    bool
    verbose(void) const
    {
        return m_verbose;
    }

    /**
     * @brief Set the verbosity mode
     *
     * @param v boolean to activate or desactivate the verbosity mode
     */
    void
    verbose(bool v)
    {
        m_verbose = v;
    }

    /**
     * @brief Initialize the initial guess of the Newton's step with a given guess
     *
     *  @param initial_solution solution \f$ u_T, u_{\partial T} \f$ for each cell
     *  @param initial_solution_faces solution \f$ u_{F} \f$ for each face
     */

    void
    initialize(const std::vector<vector_type>& initial_solution, const std::vector<vector_type>& initial_solution_faces)
    {
        m_solution_faces.clear();
        m_solution_faces = initial_solution_faces;
        assert(m_msh.faces_size() == m_solution_faces.size());

        m_solution.clear();
        m_solution = initial_solution;
        assert(m_msh.cells_size() == m_solution.size());
    }

    /**
     * @brief Compute the Newton's step until convergence or stopped criterion
     *
     * @tparam LoadIncrement Type of the loading function
     * @param lf loading function
     * @param gradient_precomputed contains the precomputed gradient for HHO methods (can be empty)
     * @param stab_precomputed contains the precomputed stabilization operators for HHO methods (can be empty)
     * @return NewtonSolverInfo Informations about the Newton's step during the computation
     */
    template<typename LoadIncrement>
    NewtonSolverInfo
    compute(const LoadIncrement&            lf,
            const std::vector<matrix_type>& gradient_precomputed,
            const std::vector<matrix_type>& stab_precomputed,
            const MeshDegreeInfo<mesh_type>&    degree_infos,
            behavior_type & behavior)
    {
        NewtonSolverInfo ni;
        timecounter      tc;
        tc.tic();

        // initialise the NewtonRaphson iteration
        NewtonIteration<mesh_type> newton_iter(m_msh, m_hdi, m_bnd, m_rp, degree_infos);

        newton_iter.initialize(m_solution, m_solution_faces);

        m_convergence = false;

        for (size_t iter = 0; iter < m_rp.m_iter_max; iter++)
        {
            // assemble lhs and rhs
            AssemblyInfo assembly_info;
            try
            {
                assembly_info = newton_iter.assemble(lf, gradient_precomputed, stab_precomputed, degree_infos, behavior);
            }
            catch (const std::invalid_argument& ia)
            {
                std::cerr << "Invalid argument: " << ia.what() << std::endl;
                m_convergence = false;
                tc.toc();
                ni.m_time_newton = tc.to_double();
                return ni;
            }

            ni.updateAssemblyInfo(assembly_info);
            // test convergence
            m_convergence = newton_iter.convergence(iter);
            if (m_convergence)
            {
                newton_iter.save_solutions(m_solution, m_solution_faces);
                tc.toc();
                ni.m_time_newton = tc.to_double();
                return ni;
            }

            // solve the global system
            SolveInfo solve_info = newton_iter.solve();
            ni.updateSolveInfo(solve_info);
            // update unknowns
            ni.m_assembly_info.m_time_postpro += newton_iter.postprocess(degree_infos);

            ni.m_iter++;
        }

        tc.toc();
        ni.m_time_newton = tc.to_double();
        return ni;
    }

    /**
     * @brief Test convergence of the Newton's iteration
     *
     * @return true if the norm of the residual is lower that a given criterion
     * @return false else
     */
    bool
    convergence() const
    {
        return m_convergence;
    }

    /**
     * @brief Save solution of the Newton's step
     *
     */
    void
    save_solutions(std::vector<vector_type>& solution, std::vector<vector_type>& solution_faces)
    {
        solution_faces.clear();
        solution_faces = m_solution_faces;
        assert(m_solution_faces.size() == solution_faces.size());

        solution.clear();
        solution = m_solution;
        assert(m_solution.size() == solution.size());
    }
};
}

} // end disk
