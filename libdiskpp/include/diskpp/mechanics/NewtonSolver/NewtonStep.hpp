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

#include "diskpp/adaptivity/adaptivity.hpp"
#include "diskpp/boundary_conditions/boundary_conditions.hpp"
#include "diskpp/common/timecounter.hpp"
#include "diskpp/mechanics/NewtonSolver/NewtonIteration.hpp"
#include "diskpp/mechanics/NewtonSolver/NewtonSolverInformations.hpp"
#include "diskpp/mechanics/NewtonSolver/NewtonSolverParameters.hpp"
#include "diskpp/mechanics/NewtonSolver/TimeManager.hpp"
#include "diskpp/mechanics/behaviors/laws/behaviorlaws.hpp"
#include "diskpp/methods/hho"

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

    bool m_verbose;
    bool m_convergence;

  public:
    NewtonStep(const param_type& rp) : m_verbose(rp.m_verbose), m_convergence(false)
    {
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
     * @brief Compute the Newton's step until convergence or stopped criterion
     *
     * @tparam LoadIncrement Type of the loading function
     * @param lf loading function
     * @param gradient_precomputed contains the precomputed gradient for HHO methods (can be empty)
     * @param stab_precomputed contains the precomputed stabilization operators for HHO methods (can be empty)
     * @return NewtonSolverInfo Informations about the Newton's step during the computation
     */
    template <typename LoadIncrement>
    NewtonSolverInfo
    compute(const mesh_type &msh,
            const bnd_type &bnd,
            const param_type &rp,
            const MeshDegreeInfo<mesh_type> &degree_infos,
            const LoadIncrement &lf,
            const TimeStep<scalar_type> &current_step,
            const std::vector<matrix_type> &gradient_precomputed,
            const std::vector<matrix_type> &stab_precomputed,
            behavior_type &behavior,
            StabCoeffManager<scalar_type> &stab_manager,
            MultiTimeField<scalar_type> &fields)
    {
        NewtonSolverInfo ni;
        timecounter      tc;
        tc.tic();

        // initialise the NewtonRaphson iteration
        NewtonIteration<mesh_type> newton_iter(msh, bnd, rp, degree_infos, current_step);

        newton_iter.initialize(msh, fields);

        m_convergence = false;

        for (size_t iter = 0; iter < rp.m_iter_max; iter++)
        {
            // assemble lhs and rhs
            AssemblyInfo assembly_info;
            try
            {
                assembly_info = newton_iter.assemble(
                    msh, bnd, rp, degree_infos, lf, gradient_precomputed, stab_precomputed, behavior, stab_manager, fields);
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
                tc.toc();
                ni.m_time_newton = tc.elapsed();
                return ni;
            }

            // solve the global system
            SolveInfo solve_info = newton_iter.solve();
            ni.updateSolveInfo(solve_info);
            // update unknowns
            ni.m_assembly_info.m_time_postpro += newton_iter.postprocess(msh, bnd, rp, degree_infos, fields);

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
};
}

} // end disk
