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
#include "TimeManager.hpp"

#include "diskpp/mechanics/behaviors/laws/behaviorlaws.hpp"
#include "diskpp/adaptivity/adaptivity.hpp"
#include "diskpp/boundary_conditions/boundary_conditions.hpp"
#include "diskpp/methods/hho"
#include "diskpp/common/timecounter.hpp"

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
    typedef Behavior<mesh_type>                   behavior_type;

    std::vector<vector_type> m_displ, m_displ_faces;
    std::vector<vector_type> m_velocity, m_acce;

    bool m_verbose;
    bool m_convergence;

  public:
    NewtonStep(const param_type& rp) : m_verbose(rp.m_verbose), m_convergence(false)
    {
        m_displ.clear();
        m_displ_faces.clear();
        m_velocity.clear();
        m_acce.clear();
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
     *  @param initial_displ solution \f$ u_T, u_{\partial T} \f$ for each cell
     *  @param initial_displ_faces solution \f$ u_{F} \f$ for each face
     */

    void
    initialize(const std::vector<vector_type>& initial_displ,
               const std::vector<vector_type>& initial_displ_faces,
               const std::vector<vector_type>& initial_velocity,
               const std::vector<vector_type>& initial_acce)
    {
        m_displ_faces.clear();
        m_displ_faces = initial_displ_faces;

        m_displ.clear();
        m_displ = initial_displ;

        m_velocity.clear();
        m_velocity = initial_velocity;

        m_acce.clear();
        m_acce = initial_acce;
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
    compute(const mesh_type&                 msh,
            const bnd_type&                  bnd,
            const param_type&                rp,
            const MeshDegreeInfo<mesh_type>& degree_infos,
            const LoadIncrement&             lf,
            const TimeStep<scalar_type>&     current_step,
            const std::vector<matrix_type>&  gradient_precomputed,
            const std::vector<matrix_type>&  stab_precomputed,
            behavior_type&                   behavior,
            StabCoeffManager<scalar_type>&   stab_manager)
    {
        NewtonSolverInfo ni;
        timecounter      tc;
        tc.tic();

        // initialise the NewtonRaphson iteration
        NewtonIteration<mesh_type> newton_iter(msh, bnd, rp, degree_infos, current_step);

        newton_iter.initialize(msh, rp, m_displ, m_displ_faces, m_velocity, m_acce);

        m_convergence = false;

        for (size_t iter = 0; iter < rp.m_iter_max; iter++)
        {
            // assemble lhs and rhs
            AssemblyInfo assembly_info;
            try
            {
                assembly_info = newton_iter.assemble(msh,
                                                     bnd,
                                                     rp,
                                                     degree_infos,
                                                     lf,
                                                     gradient_precomputed,
                                                     stab_precomputed,
                                                     behavior,
                                                     stab_manager);
            }
            catch (const std::invalid_argument& ia)
            {
                std::cerr << "Invalid argument: " << ia.what() << std::endl;
                m_convergence = false;
                tc.toc();
                ni.m_time_newton = tc.elapsed();
                return ni;
            }

            ni.updateAssemblyInfo(assembly_info);
            // test convergence
            try
            {
                m_convergence = newton_iter.convergence(rp, iter);
            }
            catch (const std::runtime_error& ia)
            {
                std::cerr << "Runtime error: " << ia.what() << std::endl;
                m_convergence = false;
                tc.toc();
                ni.m_time_newton = tc.elapsed();
                return ni;
            }

            if (m_convergence)
            {
                newton_iter.save_solutions(m_displ, m_displ_faces, m_velocity, m_acce);
                tc.toc();
                ni.m_time_newton = tc.elapsed();
                return ni;
            }

            // solve the global system
            SolveInfo solve_info = newton_iter.solve();
            ni.updateSolveInfo(solve_info);
            // update unknowns
            ni.m_assembly_info.m_time_postpro += newton_iter.postprocess(msh, bnd, rp, degree_infos);

            ni.m_iter++;
        }

        tc.toc();
        ni.m_time_newton = tc.elapsed();
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
    save_solutions(std::vector<vector_type>& displ,
                   std::vector<vector_type>& displ_faces,
                   std::vector<vector_type>& velocity,
                   std::vector<vector_type>& acce)
    {
        displ_faces.clear();
        displ_faces = m_displ_faces;
        assert(m_displ_faces.size() == displ_faces.size());

        displ.clear();
        displ = m_displ;
        assert(m_displ.size() == displ.size());

        velocity.clear();
        velocity = m_velocity;
        assert(m_velocity.size() == velocity.size());

        acce.clear();
        acce = m_acce;
        assert(m_acce.size() == acce.size());
    }
};
}

} // end disk
