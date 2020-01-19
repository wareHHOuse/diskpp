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

// NewtonRaphson_solver

#pragma once

#include <iostream>
#include <sstream>
#include <vector>

#include "NewtonSolverInformations.hpp"
#include "NewtonSolverParameters.hpp"
#include "NewtonStep.hpp"
#include "TimeManager.hpp"
#include "mechanics/behaviors/laws/behaviorlaws.hpp"
#include "mechanics/contact/ContactManager.hpp"

#include "adaptivity/adaptivity.hpp"
#include "boundary_conditions/boundary_conditions.hpp"
#include "methods/hho"

#include "timecounter.h"

namespace disk
{

namespace mechanics
{

/**
 * @brief Newton-Raphson solver for nonlinear solid mechanics
 *
 *  Specialized for HHO methods
 *
 *  Options :  - small and finite deformations
 *             - plasticity, hyperelasticity (various laws)
 *
 * @tparam Mesh type of the mesh
 */
template<typename Mesh>
class NewtonSolver
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    typedef NewtonSolverParameter<scalar_type>    param_type;
    typedef vector_boundary_conditions<mesh_type> bnd_type;
    typedef Behavior<mesh_type>                   behavior_type;

    bnd_type                  m_bnd;
    const mesh_type&          m_msh;
    param_type                m_rp;
    behavior_type             m_behavior;
    MeshDegreeInfo<mesh_type> m_degree_infos;
    ContactManager<mesh_type> m_contact_manager;

    std::vector<vector_type> m_solution, m_solution_faces, m_solution_mult;
    std::vector<matrix_type> m_gradient_precomputed, m_stab_precomputed;

    bool m_verbose, m_convergence;

    void
    init_degree(const size_t cell_degree, const size_t face_degree, const size_t grad_degree)
    {
        m_degree_infos = MeshDegreeInfo<mesh_type>(m_msh, cell_degree, face_degree, grad_degree);

        if (m_bnd.nb_faces_contact() > 0)
        {
            for (auto itor = m_msh.boundary_faces_begin(); itor != m_msh.boundary_faces_end(); itor++)
            {
                const auto bfc     = *itor;
                const auto face_id = m_msh.lookup(bfc);

                if (m_bnd.is_contact_face(face_id))
                {
                    switch (m_bnd.contact_boundary_type(face_id))
                    {
                        case SIGNORINI_FACE:
                        {
                            m_degree_infos.degree(m_msh, bfc, face_degree + 1);
                            break;
                        }
                        default: { throw std::invalid_argument("Invalid contact type");
                        }
                    }
                }
            }
        }
    }

    // Initializa data structures
    void
    init(void)
    {
        const auto dimension = mesh_type::dimension;

        size_t total_dof = 0;

        m_solution.clear();
        m_solution_faces.clear();

        m_solution.reserve(m_msh.cells_size());
        m_solution_faces.reserve(m_msh.faces_size());

        for (auto& cl : m_msh)
        {
            const auto di            = m_degree_infos.degreeInfo(m_msh, cl);
            const auto num_cell_dofs = vector_basis_size(di.degree(), dimension, dimension);
            total_dof += num_cell_dofs;

            const auto fcs    = faces(m_msh, cl);
            const auto fcs_di = m_degree_infos.degreeInfo(m_msh, fcs);

            size_t num_faces_dofs = 0;
            for (auto& fc : fcs_di)
            {
                if (fc.hasUnknowns())
                {
                    num_faces_dofs += vector_basis_size(fc.degree(), dimension - 1, dimension);
                }
            }

            m_solution.push_back(vector_type::Zero(num_cell_dofs + num_faces_dofs));
        }

        for (auto itor = m_msh.faces_begin(); itor != m_msh.faces_end(); itor++)
        {
            const auto fc = *itor;
            const auto di = m_degree_infos.degreeInfo(m_msh, fc);

            size_t num_face_dofs = 0;
            if (di.hasUnknowns())
            {
                num_face_dofs = vector_basis_size(di.degree(), dimension - 1, dimension);
            }
            total_dof += num_face_dofs;

            m_solution_faces.push_back(vector_type::Zero(num_face_dofs));
        }

        if (m_verbose)
        {
            std::cout << "** Numbers of cells: " << m_msh.cells_size() << std::endl;
            std::cout << "** Numbers of faces: " << m_msh.faces_size()
                      << "  ( boundary faces: " << m_msh.boundary_faces_size() << " )" << std::endl;
            std::cout << "** Numbers of dofs: " << std::endl;
            std::cout << "   ** Before static condensation: " << total_dof << std::endl;
            std::cout << "   ** After static condensation: " << this->numberOfDofs() << std::endl;
            std::cout << " " << std::endl;
        }
    }

    /**
     * @brief Precompute the gradient reconstruction and stabilization operator for HHO methods
     *
     */
    void
    pre_computation(void)
    {
        m_gradient_precomputed.clear();
        m_gradient_precomputed.reserve(m_msh.cells_size());

        m_stab_precomputed.clear();
        m_stab_precomputed.reserve(m_msh.cells_size());

        for (auto& cl : m_msh)
        {
            //std::cout << m_degree_infos.cellDegreeInfo(m_msh, cl) << std::endl;
            /////// Gradient Reconstruction /////////
            if (m_behavior.getDeformation() == SMALL_DEF)
            {
                const auto sgr = make_matrix_symmetric_gradrec(m_msh, cl, m_degree_infos);
                m_gradient_precomputed.push_back(sgr.first);
            }
            else
            {
                const auto gr = make_marix_hho_gradrec(m_msh, cl, m_degree_infos);
                m_gradient_precomputed.push_back(gr.first);
            }

            if (m_rp.m_stab)
            {
                switch (m_rp.m_stab_type)
                {
                    case HHO:
                    {
                        if (m_behavior.getDeformation() == SMALL_DEF)
                        {
                            const auto recons = make_vector_hho_symmetric_laplacian(m_msh, cl, m_degree_infos);
                            m_stab_precomputed.push_back(
                              make_vector_hho_stabilization(m_msh, cl, recons.first, m_degree_infos));
                        }
                        else
                        {
                            const auto recons_scalar = make_scalar_hho_laplacian(m_msh, cl, m_degree_infos);
                            m_stab_precomputed.push_back(
                              make_vector_hho_stabilization_optim(m_msh, cl, recons_scalar.first, m_degree_infos));
                        }
                        break;
                    }
                    case HDG:
                    {
                        m_stab_precomputed.push_back(make_vector_hdg_stabilization(m_msh, cl, m_degree_infos));
                        break;
                    }
                    case DG:
                    {
                        m_stab_precomputed.push_back(make_vector_dg_stabilization(m_msh, cl, m_degree_infos));
                        break;
                    }
                    case NO: { break;
                    }
                    default: throw std::invalid_argument("Unknown stabilization");
                }
            }
        }
    }

  public:
    NewtonSolver(const mesh_type& msh, const bnd_type& bnd, const param_type& rp) :
      m_msh(msh), m_verbose(rp.m_verbose), m_convergence(false), m_rp(rp), m_bnd(bnd)
    {
        if (m_verbose)
        {
            std::cout << "------------------------------------------------------------------------"
                         "----------------------"
                      << std::endl;
            std::cout
              << "|********************** Nonlinear Newton's Solver for solid mechanics ***********************|"
              << std::endl;
            std::cout << "------------------------------------------------------------------------"
                         "----------------------"
                      << std::endl;
        }
        int face_degree = rp.m_face_degree;
        if (rp.m_face_degree < 0)
        {
            std::cout << "'face_degree' should be > 0. Reverting to 1." << std::endl;
            face_degree = 1;
        }

        m_rp.m_face_degree = face_degree;

        int cell_degree = rp.m_cell_degree;
        if ((face_degree - 1 > cell_degree) or (cell_degree > face_degree + 1))
        {
            std::cout << "'cell_degree' should be 'face_degree + 1' =>"
                      << "'cell_degree' => 'face_degree -1'. Reverting to 'face_degree'." << std::endl;
            cell_degree = face_degree;
        }

        m_rp.m_cell_degree = cell_degree;

        int grad_degree = rp.m_grad_degree;
        if (grad_degree < face_degree)
        {
            std::cout << "'grad_degree' should be >= 'face_degree'. Reverting to 'face_degree'." << std::endl;
            grad_degree = face_degree;
        }

        if (m_verbose)
        {
            m_rp.infos();
        }

        // Initialization
        if (m_verbose)
            std::cout << "Initialization ..." << std::endl;
        this->init_degree(cell_degree, face_degree, grad_degree);
        this->init();

        // Precomputation
        if (m_rp.m_precomputation)
        {
            timecounter t1;
            t1.tic();
            this->pre_computation();
            t1.toc();
            if (m_verbose)
                std::cout << "Precomputation: " << t1.to_double() << " sec" << std::endl;
        }
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
     * @brief Initialize the inital guess with a given function
     *
     * @param func given function
     */
    void
    initial_guess(const vector_rhs_function<mesh_type> func)
    {
        size_t cell_i = 0;

        for (auto& cl : m_msh)
        {
            m_solution.at(cell_i++) = project_function(m_msh, cl, m_degree_infos, func, 2);
        }

        for (auto itor = m_msh.faces_begin(); itor != m_msh.faces_end(); itor++)
        {
            const auto bfc     = *itor;
            const auto face_id = m_msh.lookup(bfc);
            const auto fdi     = m_degree_infos.degreeInfo(m_msh, bfc);

            if (m_bnd.contact_boundary_type(face_id) == SIGNORINI_FACE)
            {
                const auto proj_bcf = project_function(m_msh, bfc, fdi.degree() + 1, func, 2);
                assert(m_solution_faces[face_id].size() == proj_bcf.size());

                m_solution_faces[face_id] = proj_bcf;
            }
            else if (m_bnd.contact_boundary_type(face_id) == SIGNORINI_CELL)
            {
                assert(m_solution_faces[face_id].size() == 0);
            }
            else
            {
                const auto proj_bcf = project_function(m_msh, bfc, fdi.degree(), func, 2);
                assert(m_solution_faces[face_id].size() == proj_bcf.size());

                m_solution_faces[face_id] = proj_bcf;
            }
        }
    }

    /**
     * @brief Add a behavior for materials
     *
     * @param deformation Type of deformation
     * @param law Type of Law
     */
    void
    addBehavior(const size_t deformation, const size_t law)
    {
        m_behavior = behavior_type(m_msh, 2 * m_rp.m_grad_degree, deformation, law);

        if (m_verbose)
        {
            std::cout << "Add behavior ..." << std::endl;
            std::cout << "** Number of integration points: " << m_behavior.numberOfQP() << std::endl;
        }
    }

    /**
     * @brief Add a behavior for materials (by copy)
     *
     * @param behavior Given behavior
     */
    void
    addBehavior(const behavior_type& behavior)
    {
        m_behavior = behavior;
        if (m_verbose)
        {
            std::cout << "Add behavior ..." << std::endl;
            std::cout << "** Number of integration points: " << m_behavior.numberOfQP() << std::endl;
        }
    }

    /**
     * @brief Add material properties for the behavior
     *
     * @param material_data material properties
     */
    void
    addMaterialData(const MaterialData<scalar_type>& material_data)
    {
        m_behavior.addMaterialData(material_data);

        if (m_verbose)
        {
            std::cout << "Add material ..." << std::endl;
            m_behavior.getMaterialData().print();
        }
    }

    template<typename LoadFunction>
    SolverInfo
    compute(const LoadFunction& lf)
    {
        SolverInfo  si;
        timecounter ttot;
        ttot.tic();

        // list of time step
        ListOfTimeStep<scalar_type> list_time_step(m_rp.m_time_step);

        // time of saving
        bool time_saving = false;
        if (m_rp.m_n_time_save > 0)
        {
            time_saving = true;
        }

        // Newton step
        NewtonStep<mesh_type> newton_step(m_rp);
        newton_step.initialize(m_solution, m_solution_faces, m_solution_mult);

        // Loop on time step
        while (!list_time_step.empty())
        {
            const auto current_step = list_time_step.getCurrentTimeStep();
            const auto current_time = current_step.end_time();

            if (m_verbose)
            {
                list_time_step.printCurrentTimeStep();
            }

            auto rlf = [&lf, &current_time ](const point<scalar_type, mesh_type::dimension>& p) -> auto
            {
                return priv::inner_product(current_time, lf(p));
            };

            m_bnd.multiplyAllFunctionsByAFactor(current_time);

            //  Newton correction
            NewtonSolverInfo newton_info = newton_step.compute(m_msh, m_bnd, m_rp, m_degree_infos,
              rlf, m_gradient_precomputed, m_stab_precomputed, m_behavior, m_contact_manager);
            si.updateInfo(newton_info);

            if (m_verbose)
            {
                newton_info.printInfo();
            }

            // Test convergence
            m_convergence = newton_step.convergence();

            if (!m_convergence)
            {
                if (current_step.level() + 1 > m_rp.m_sublevel)
                {
                    std::cout << "***********************************************************" << std::endl;
                    std::cout << "***** PROBLEM OF CONVERGENCE: We stop the calcul here *****" << std::endl;
                    std::cout << "***********************************************************" << std::endl;
                    break;
                }
                else
                {
                    if (m_verbose)
                    {
                        std::cout << "***********************************************************" << std::endl;
                        std::cout << "*****     NO CONVERGENCE: We split the time step     ******" << std::endl;
                        std::cout << "***********************************************************" << std::endl;
                    }

                    list_time_step.splitCurrentTimeStep();
                }
            }
            else
            {
                list_time_step.removeCurrentTimeStep();
                m_behavior.update();

                if (time_saving)
                {
                    if (m_rp.m_time_save.front() < current_time + 1E-5)
                    {
                        newton_step.save_solutions(m_solution, m_solution_faces, m_solution_mult);

                        std::cout << "** Save results" << std::endl;
                        std::string name =
                          "result" + std::to_string(mesh_type::dimension) + "D_t" + std::to_string(current_time) + "_";

                        m_rp.m_time_save.pop_front();
                        if (m_rp.m_time_save.empty())
                            time_saving = false;
                    }
                }
            }
        }

        // save solutions
        newton_step.save_solutions(m_solution, m_solution_faces, m_solution_mult);
        si.m_time_step = list_time_step.numberOfTimeStep();

        ttot.toc();
        si.m_time_solver = ttot.to_double();

        return si;
    }

    bool
    convergence() const
    {
        return m_convergence;
    }

    size_t
    numberOfDofs()
    {
        const auto dimension      = mesh_type::dimension;
        size_t     num_faces_dofs = 0;
        for (auto itor = m_msh.faces_begin(); itor != m_msh.faces_end(); itor++)
        {
            const auto fc = *itor;
            const auto di = m_degree_infos.degreeInfo(m_msh, fc);

            if (di.hasUnknowns())
            {
                num_faces_dofs += vector_basis_size(di.degree(), dimension - 1, dimension);
            }
        }
        return num_faces_dofs;
    }

    void
    printSolutionCell() const
    {
        size_t cell_i = 0;
        std::cout << "Solution at the cells:" << std::endl;
        for (auto& cl : m_msh)
        {
            const auto di            = m_degree_infos.degreeInfo(m_msh, cl);
            const auto num_cell_dofs = vector_basis_size(di.cell_degree(), mesh_type::dimension, mesh_type::dimension);
            std::cout << "cell " << cell_i << ": " << std::endl;
            std::cout << m_solution.at(cell_i++).head(num_cell_dofs).transpose() << std::endl;
        }
    }

    // compute l2 error
    template<typename AnalyticalSolution>
    scalar_type
    compute_l2_displacement_error(const AnalyticalSolution& as)
    {
        scalar_type err_dof = 0;

        size_t cell_i = 0;

        for (auto& cl : m_msh)
        {
            const auto cdi           = m_degree_infos.degreeInfo(m_msh, cl);
            const auto num_cell_dofs = vector_basis_size(cdi.degree(), mesh_type::dimension, mesh_type::dimension);
            const vector_type comp_dof = m_solution.at(cell_i++).head(num_cell_dofs);
            const vector_type true_dof = project_function(m_msh, cl, cdi.degree(), as, 2);

            const auto        cb   = make_vector_monomial_basis(m_msh, cl, cdi.degree());
            const matrix_type mass = make_mass_matrix(m_msh, cl, cb);

            const vector_type diff_dof = (true_dof - comp_dof);
            assert(comp_dof.size() == true_dof.size());
            err_dof += diff_dof.dot(mass * diff_dof);
        }

        return sqrt(err_dof);
    }

    // compute l2 error
    template<typename AnalyticalSolution>
    scalar_type
    compute_H1_error(const AnalyticalSolution& as)
    {
        scalar_type err_dof = 0;

        size_t cell_i = 0;

        matrix_type grad;
        matrix_type stab;

        for (auto& cl : m_msh)
        {
            // std::cout << m_degree_infos.cellDegreeInfo(m_msh, cl) << std::endl;
            /////// Gradient Reconstruction /////////
            if (m_behavior.getDeformation() == SMALL_DEF)
            {
                grad = make_matrix_symmetric_gradrec(m_msh, cl, m_degree_infos).second;
            }
            else
            {
                grad = make_marix_hho_gradrec(m_msh, cl, m_degree_infos).second;
            }

            if (m_rp.m_stab)
            {
                switch (m_rp.m_stab_type)
                {
                    case HHO:
                    {
                        if (m_behavior.getDeformation() == SMALL_DEF)
                        {
                            const auto recons = make_vector_hho_symmetric_laplacian(m_msh, cl, m_degree_infos);
                            stab = make_vector_hho_stabilization(m_msh, cl, recons.first, m_degree_infos);
                        }
                        else
                        {
                            const auto recons_scalar = make_scalar_hho_laplacian(m_msh, cl, m_degree_infos);
                            stab = make_vector_hho_stabilization_optim(m_msh, cl, recons_scalar.first, m_degree_infos);
                        }
                        break;
                    }
                    case HDG:
                    {
                        stab = make_vector_hdg_stabilization(m_msh, cl, m_degree_infos);
                        break;
                    }
                    case DG:
                    {
                        stab = make_vector_dg_stabilization(m_msh, cl, m_degree_infos);
                        break;
                    }
                    case NO: { break;
                        stab.setZero();
                    }
                    default: throw std::invalid_argument("Unknown stabilization");
                }
            }

            const auto Ah = grad + stab;

            const vector_type comp_dof = m_solution.at(cell_i);
            const vector_type true_dof = project_function(m_msh, cl, m_degree_infos, as, 2);

            const vector_type diff_dof = (true_dof - comp_dof);
            assert(comp_dof.size() == true_dof.size());
            err_dof += diff_dof.dot(Ah * diff_dof);

            cell_i++;
        }

        return sqrt(err_dof);
    }
};
}

} // end disk
