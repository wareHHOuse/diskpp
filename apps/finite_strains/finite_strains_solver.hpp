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

#include <iostream>
#include <sstream>

#include <list>
#include <vector>

#include "revolution/bases"
#include "revolution/methods/hho"
#include "revolution/quadratures"

#include "Informations.hpp"
#include "NewtonSolver/newton_solver.hpp"
#include "Parameters.hpp"
#include "mechanics/BoundaryConditions.hpp"
#include "mechanics/behaviors/laws/behaviorlaws.hpp"

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"
#include "output/postMesh.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

struct time_step
{
    double time;
    size_t level;
};

template<typename Mesh>
class finite_strains_solver
{
    typedef Mesh                                 mesh_type;
    typedef typename mesh_type::scalar_type      scalar_type;
    typedef ParamRun<scalar_type>                param_type;
    typedef NLE::MaterialParameters<scalar_type> data_type;

    typedef dynamic_matrix<scalar_type> matrix_dynamic;
    typedef dynamic_vector<scalar_type> vector_dynamic;

    typedef disk::mechanics::BoundaryConditions<mesh_type>        bnd_type;
    typedef disk::LinearLaw<mesh_type> law_type;

    typename revolution::hho_degree_info m_hdi;
    bnd_type                             m_bnd;
    const mesh_type&                     m_msh;

    disk::PostMesh<mesh_type> post_mesh;

    std::vector<vector_dynamic> m_solution_data;
    std::vector<vector_dynamic> m_solution_cells, m_solution_faces;
    std::vector<matrix_dynamic> m_gradient_precomputed, m_stab_precomputed;

    law_type m_law;

    bool m_verbose, m_convergence;

    param_type m_rp;

    const static size_t dimension = mesh_type::dimension;

    int total_dof_depl_static;

    void
    init()
    {
        m_solution_data.clear();
        m_solution_cells.clear();
        m_solution_faces.clear();

        m_solution_data.reserve(m_msh.cells_size());
        m_solution_cells.reserve(m_msh.cells_size());
        m_solution_faces.reserve(m_msh.faces_size());

        const auto num_cell_dofs = revolution::vector_basis_size(m_hdi.cell_degree(), dimension, dimension);
        const auto num_face_dofs = revolution::vector_basis_size(m_hdi.face_degree(), dimension - 1, dimension);
        const auto total_dof     = m_msh.cells_size() * num_cell_dofs + m_msh.faces_size() * num_face_dofs;

        for (auto& cl : m_msh)
        {
            const auto fcs       = faces(m_msh, cl);
            const auto num_faces = fcs.size();
            m_solution_data.push_back(vector_dynamic::Zero(num_cell_dofs + num_faces * num_face_dofs));
            m_solution_cells.push_back(vector_dynamic::Zero(num_cell_dofs));
        }

        for (int i = 0; i < m_msh.faces_size(); i++)
        {
            m_solution_faces.push_back(vector_dynamic::Zero(num_face_dofs));
        }

        if (m_verbose)
        {
            std::cout << "** Numbers of cells: " << m_msh.cells_size() << std::endl;
            std::cout << "** Numbers of faces: " << m_msh.faces_size()
                      << " ( boundary faces: " << m_msh.boundary_faces_size() << " )" << std::endl;
            std::cout << "** Numbers of dofs: " << total_dof << std::endl;
            std::cout << "** After static condensation: " << std::endl;
            std::cout << "** Numbers of dofs: " << m_msh.faces_size() * num_face_dofs << std::endl;
            std::cout << "** Number of integration points: " << m_law.getNumberOfQP() << std::endl;
            std::cout << " " << std::endl;

            m_law.getMaterialData().print();
        }

        total_dof_depl_static = m_msh.faces_size() * num_face_dofs;
    }

    void
    pre_computation()
    {
        m_gradient_precomputed.clear();
        m_gradient_precomputed.reserve(m_msh.cells_size());

        m_stab_precomputed.clear();
        m_stab_precomputed.reserve(m_msh.cells_size());

        for (auto& cl : m_msh)
        {
            /////// Gradient Reconstruction /////////
            const auto gr = make_hho_gradrec_matrix(m_msh, cl, m_hdi);
            m_gradient_precomputed.push_back(gr.first);

            if (m_rp.m_stab)
            {
                switch (m_rp.m_stab_type)
                {
                    case HHO:
                    {
                        const auto recons_scalar = make_hho_scalar_laplacian(m_msh, cl, m_hdi);
                        const auto stab_HHO = make_hho_vector_stabilization_optim(m_msh, cl, recons_scalar.first, m_hdi);
                        m_stab_precomputed.push_back(stab_HHO);
                        break;
                    }
                    case HDG:
                    {
                        const auto stab_HDG = make_hdg_vector_stabilization(m_msh, cl, m_hdi);
                        m_stab_precomputed.push_back(stab_HDG);
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
    finite_strains_solver(const mesh_type&  msh,
                          const bnd_type&   bnd,
                          const param_type& rp,
                          const data_type&  material_data) :
      m_msh(msh),
      m_verbose(rp.m_verbose), m_convergence(false), m_rp(rp), m_bnd(bnd)
    {
        int face_degree = rp.m_face_degree;
        if (rp.m_face_degree < 0)
        {
            std::cout << "'face_degree' should be > 0. Reverting to 1." << std::endl;
            face_degree = 1;
        }

        m_rp.m_face_degree = face_degree;

        int cell_degree = rp.m_cell_degree;
        // if (face_degree - 1 > cell_degree or cell_degree > face_degree + 1) {
        //    std::cout << "'cell_degree' should be 'face_degree + 1' => 'cell_degree' => 'face_degree
        //    "
        //                 "-1'. Reverting to 'face_degree'."
        //              << std::endl;
        //    cell_degree = face_degree;
        // }

        m_rp.m_cell_degree = cell_degree;

        int grad_degree = rp.m_grad_degree;
        if (grad_degree < face_degree)
        {
            std::cout << "'grad_degree' should be >= 'face_degree'. Reverting to 'face_degree'." << std::endl;
            grad_degree = face_degree;
        }

        m_rp.m_grad_degree = grad_degree;

        m_hdi = revolution::hho_degree_info(cell_degree, face_degree, grad_degree);

        m_law = law_type(m_msh, 2 * m_hdi.grad_degree());
        m_law.addMaterialData(
          material_data.lambda, material_data.mu);

        if (m_verbose)
        {
            m_hdi.info_degree();
            m_rp.infos();
        }

        init();
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

    template<typename LoadFunction>
    SolverInfo
    compute(const LoadFunction& lf)
    {

        SolverInfo  si;
        timecounter ttot;
        ttot.tic();

        // time step
        std::list<time_step> list_step;

        scalar_type time1 = 0.0;
        scalar_type time2 = 0.0;
        for (int n = 0; n < m_rp.m_time_step.size(); n++)
        {
            auto time_info            = m_rp.m_time_step[n];
            time2                     = time_info.first;
            const scalar_type delta_t = (time2 - time1) / time_info.second;
            for (int i = 0; i < time_info.second; i++)
            {
                time_step step;
                step.time  = time1 + (i + 1) * delta_t;
                step.level = 1;
                list_step.push_back(step);
            }
            time1 = time2;
        }

        int current_step = 0;
        int total_step   = list_step.size();

        scalar_type old_time = 0.0;

        // time of saving
        bool time_saving(false);
        if (m_rp.m_n_time_save > 0)
        {
            time_saving = true;
        }

        // Precomputation
        if (m_rp.m_precomputation)
        {
            timecounter t1;
            t1.tic();
            this->pre_computation();
            t1.toc();
            if (m_verbose)
                std::cout << "-Precomputation: " << t1.to_double() << " sec" << std::endl;
        }

        // Newton solver
        NLE::NewtonRaphson_solver_finite_strains<mesh_type> newton_solver(m_msh, m_hdi, m_bnd, m_rp);

        newton_solver.initialize(m_solution_cells, m_solution_faces, m_solution_data);

        // loading
        while (!list_step.empty())
        {
            current_step += 1;
            time_step         step         = list_step.front();
            const scalar_type current_time = step.time;

            if (m_verbose)
            {
                std::cout << "------------------------------------------------------------------------"
                             "----------------------"
                          << std::endl;
                std::cout << "************************ Time : " << current_time << " sec (step: " << current_step << "/"
                          << total_step << ", sublevel: " << step.level << " ) *************************|" << std::endl;
            }

            auto rlf = [&lf, &current_time ](const point<scalar_type, dimension>& p) -> auto
            {
                return revolution::priv::inner_product(current_time, lf(p));
            };

            m_bnd.multiplyAllFunctionsOfAFactor(current_time);

            // correction
            NewtonSolverInfo newton_info =
              newton_solver.compute(rlf, m_gradient_precomputed, m_stab_precomputed, m_law);
            si.updateInfo(newton_info);

            if (m_verbose)
            {
                std::cout << "** Time in this step " << newton_info.m_time_newton << " sec" << std::endl;
                std::cout << "**** Assembly time: " << newton_info.m_assembly_info.m_time_assembly << " sec"
                          << std::endl;
                std::cout << "****** Gradient reconstruction: " << newton_info.m_assembly_info.m_time_gradrec << " sec"
                          << std::endl;
                std::cout << "****** Stabilisation: " << newton_info.m_assembly_info.m_time_stab << " sec" << std::endl;
                std::cout << "****** Mechanical computation: " << newton_info.m_assembly_info.m_time_elem << " sec"
                          << std::endl;
                std::cout << "       *** Behavior computation: " << newton_info.m_assembly_info.m_time_law << " sec"
                          << std::endl;
                std::cout << "****** Static condensation: " << newton_info.m_assembly_info.m_time_statcond << " sec"
                          << std::endl;
                std::cout << "****** Postprocess time: " << newton_info.m_assembly_info.m_time_postpro << " sec"
                          << std::endl;
                std::cout << "**** Solver time: " << newton_info.m_solve_info.m_time_solve << " sec" << std::endl;
            }

            m_convergence = newton_solver.test_convergence();

            if (!m_convergence)
            {
                if (step.level > m_rp.m_sublevel)
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
                    total_step += 1;
                    current_step -= 1;
                    time_step new_step;
                    new_step.time  = old_time + (current_time - old_time) / 2.0;
                    new_step.level = step.level + 1;
                    list_step.push_front(new_step);
                }
            }
            else
            {
                old_time = current_time;
                list_step.pop_front();
                newton_solver.save_solutions(m_solution_cells, m_solution_faces, m_solution_data);
                m_law.update();

                if (time_saving)
                {
                    if (m_rp.m_time_save.front() < old_time + 1E-5)
                    {
                        m_rp.m_time_save.pop_front();
                        if (m_rp.m_time_save.empty())
                            time_saving = false;
                    }
                }
            }
        }

        si.m_time_step = total_step;

        ttot.toc();
        si.m_time_solver = ttot.to_double();

        // compute mesh for post-processing
        post_mesh = disk::PostMesh<mesh_type>(m_msh);

        return si;
    }

    bool
    test_convergence() const
    {
        return m_convergence;
    }

    int
    getDofs()
    {
        return total_dof_depl_static;
    }

    void printSolutionCell() const
    {
        int cell_i = 0;
        std::cout << "Solution at the cells:" << std::endl;
        for (auto& cl : m_msh)
        {
            std::cout << "cell " << cell_i << ": " << std::endl;
            std::cout << m_solution_cells.at(cell_i++) << std::endl;
        }
    }

    // compute l2 error
    template<typename AnalyticalSolution>
    scalar_type
    compute_l2_displacement_error(const AnalyticalSolution& as)
    {
        scalar_type err_dof = 0;

        const auto cbs      = revolution::vector_basis_size(m_hdi.cell_degree(), dimension, dimension);
        const int  diff_deg = m_hdi.face_degree() - m_hdi.cell_degree();
        const int  di       = std::max(diff_deg, 0);

        int cell_i = 0;

        for (auto& cl : m_msh)
        {
            const auto x = m_solution_cells.at(cell_i++);

            const vector_dynamic true_dof = revolution::project_function(m_msh, cl, m_hdi.cell_degree(), as, di);

            auto                 cb   = revolution::make_vector_monomial_basis(m_msh, cl, m_hdi.cell_degree());
            const matrix_dynamic mass = revolution::make_mass_matrix(m_msh, cl, cb);

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
        const auto face_degree   = m_hdi.face_degree();
        const auto material_data = m_law.getMaterialData();

        int         cell_i = 0;
        scalar_type error_stress(0.0);

        for (auto& cl : m_msh)
        {
            const auto law_quadpoints = m_law.getCellIVs(cell_i).getIVs();

            // Loop on nodes
            for (auto& qp : law_quadpoints)
            {
                const auto stress_comp = qp.compute_stress(material_data);
                const auto stress_diff = (stress(qp.point()) - stress_comp).eval();

                error_stress += qp.weight() * revolution::priv::inner_product(stress_diff, stress_diff);
            }
        }

        return sqrt(error_stress);
    }

    // post-processing
    scalar_type
    energy_mechanic() const
    {
        const auto  material_data = m_law.getMaterialData();
        scalar_type energy        = 0.;

        int cell_i = 0;
        for (auto& cl : m_msh)
        {
            const auto law_quadpoints = m_law.getCellIVs(cell_i).getIVs();

            // Loop on nodes
            for (auto& qp : law_quadpoints)
            {
                const auto                stress = qp.compute_stress(material_data);
                const std::vector<double> tens   = disk::convertToVectorGmsh(stress);

                energy += qp.weight() * revolution::priv::inner_product(stress, stress);
            }
            cell_i++;
        }

        return sqrt(energy);
    }

    void
    compute_discontinuous_displacement(const std::string& filename) const
    {
        const auto cell_degree = m_hdi.cell_degree();

        gmsh::Gmesh gmsh(dimension);
        auto        storage = m_msh.backend_storage();

        std::vector<gmsh::Data>          data;    // create data (not used)
        const std::vector<gmsh::SubData> subdata; // create subdata to save soution at gauss point

        int cell_i   = 0;
        int nb_nodes = 0;
        for (auto& cl : m_msh)
        {
            auto                    cb         = revolution::make_vector_monomial_basis(m_msh, cl, cell_degree);
            const vector_dynamic    x          = m_solution_cells.at(cell_i++);
            const auto              cell_nodes = disk::cell_nodes(m_msh, cl);
            std::vector<gmsh::Node> new_nodes;

            // loop on the nodes of the cell
            for (int i = 0; i < cell_nodes.size(); i++)
            {
                nb_nodes++;
                const auto point_ids = cell_nodes[i];
                const auto pt        = storage->points[point_ids];

                const auto phi = cb.eval_functions(pt);

                const auto depl = revolution::eval(x, phi);

                const std::vector<double>   deplv = disk::convertToVectorGmsh(depl);
                const std::array<double, 3> coor  = disk::init_coor(pt);

                // Add a node
                const gmsh::Node tmp_node(coor, nb_nodes, 0);
                new_nodes.push_back(tmp_node);
                gmsh.addNode(tmp_node);

                const gmsh::Data datatmp(nb_nodes, deplv);
                data.push_back(datatmp);
            }
            // Add new element
            disk::add_element(gmsh, new_nodes);
        }

        // Create and init a nodedata view
        gmsh::NodeData nodedata(3, 0.0, "depl_node_disc", data, subdata);

        // Save the view
        nodedata.saveNodeData(filename, gmsh);
    }

    void
    compute_continuous_displacement(const std::string& filename) const
    {
        const auto cell_degree = m_hdi.cell_degree();

        gmsh::Gmesh gmsh    = disk::convertMesh(post_mesh);
        auto        storage = post_mesh.mesh().backend_storage();

        const static_vector<scalar_type, dimension> vzero = static_vector<scalar_type, dimension>::Zero();

        const size_t nb_nodes(gmsh.getNumberofNodes());

        // first(number of data at this node), second(cumulated value)
        std::vector<std::pair<size_t, static_vector<scalar_type, dimension>>> value(nb_nodes, std::make_pair(0, vzero));

        int cell_i = 0;
        for (auto& cl : m_msh)
        {
            auto                 cb         = revolution::make_vector_monomial_basis(m_msh, cl, cell_degree);
            const vector_dynamic x          = m_solution_cells.at(cell_i);
            const auto           cell_nodes = post_mesh.nodes_cell(cell_i);

            // Loop on the nodes of the cell
            for (auto& point_id : cell_nodes)
            {
                const auto pt = storage->points[point_id];

                const auto phi  = cb.eval_functions(pt);
                const auto depl = revolution::eval(x, phi);

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
        for (int i_node = 0; i_node < value.size(); i_node++)
        {
            const static_vector<scalar_type, dimension> depl_avr = value[i_node].second / double(value[i_node].first);

            const gmsh::Data tmp_data(i_node + 1, disk::convertToVectorGmsh(depl_avr));
            data.push_back(tmp_data);
        }

        // Create and init a nodedata view
        gmsh::NodeData nodedata(3, 0.0, "depl_node_cont", data, subdata);
        // Save the view
        nodedata.saveNodeData(filename, gmsh);
    }

    void
    compute_stress_GP(const std::string& filename) const
    {
        gmsh::Gmesh gmsh = disk::convertMesh(post_mesh);

        std::vector<gmsh::Data>    data;    // create data (not used)
        std::vector<gmsh::SubData> subdata; // create subdata to save soution at gauss point
        size_t                     nb_nodes(gmsh.getNumberofNodes());

        const auto material_data = m_law.getMaterialData();

        int cell_i = 0;
        for (auto& cl : m_msh)
        {

            const auto law_quadpoints = m_law.getCellIVs(cell_i).getIVs();

            // Loop on nodes
            for (auto& qp : law_quadpoints)
            {
                const auto                stress = qp.compute_stress(material_data);
                const std::vector<double> tens   = disk::convertToVectorGmsh(stress);

                // Add GP
                // Create a node at gauss point
                nb_nodes++;
                const gmsh::Node    new_node = disk::convertPoint(qp.point(), nb_nodes);
                const gmsh::SubData sdata(tens, new_node);
                subdata.push_back(sdata); // add subdata
            }
            cell_i++;
        }

        // Save
        gmsh::NodeData nodedata(9, 0.0, "Stress_GP", data, subdata); // create and init a nodedata view

        nodedata.saveNodeData(filename, gmsh); // save the view
    }

    void
    compute_discontinuous_stress(const std::string& filename) const
    {
        gmsh::Gmesh gmsh(dimension);
        auto        storage = m_msh.backend_storage();

        std::vector<gmsh::Data>          data;    // create data (not used)
        const std::vector<gmsh::SubData> subdata; // create subdata to save soution at gauss point

        const auto material_data = m_law.getMaterialData();

        int    cell_i   = 0;
        size_t nb_nodes = 0;
        for (auto& cl : m_msh)
        {
            auto                    gb = revolution::make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());
            const auto              law_cell     = m_law.getCellIVs(cell_i);
            const vector_dynamic    stress_coeff = law_cell.projectStressOnCell(m_msh, cl, m_hdi, material_data);
            const auto              cell_nodes   = disk::cell_nodes(m_msh, cl);
            std::vector<gmsh::Node> new_nodes;

            // loop on the nodes of the cell
            for (int i = 0; i < cell_nodes.size(); i++)
            {
                nb_nodes++;
                const auto point_ids = cell_nodes[i];
                const auto pt        = storage->points[point_ids];

                const auto gphi   = gb.eval_functions(pt);
                const auto stress = revolution::eval(stress_coeff, gphi);

                const std::vector<double>   tens = disk::convertToVectorGmsh(stress);
                const std::array<double, 3> coor = disk::init_coor(pt);

                // Add a node
                const gmsh::Node tmp_node(coor, nb_nodes, 0);
                new_nodes.push_back(tmp_node);
                gmsh.addNode(tmp_node);

                const gmsh::Data datatmp(nb_nodes, tens);
                data.push_back(datatmp);
            }

            // Add new element
            disk::add_element(gmsh, new_nodes);
            cell_i++;
        }

        // Create and init a nodedata view
        gmsh::NodeData nodedata(9, 0.0, "stress_node_disc", data, subdata);

        // Save the view
        nodedata.saveNodeData(filename, gmsh);
    }

    void
    compute_continuous_stress(const std::string& filename) const
    {
        gmsh::Gmesh gmsh    = disk::convertMesh(post_mesh);
        auto        storage = post_mesh.mesh().backend_storage();

        const static_matrix<scalar_type, dimension, dimension> vzero =
          static_matrix<scalar_type, dimension, dimension>::Zero();

        const size_t nb_nodes(gmsh.getNumberofNodes());

        // first(number of data at this node), second(cumulated value)
        std::vector<std::pair<size_t, static_matrix<scalar_type, dimension, dimension>>> value(
          nb_nodes, std::make_pair(0, vzero));

        const auto material_data = m_law.getMaterialData();

        int cell_i = 0;

        for (auto& cl : m_msh)
        {
            auto                 gb       = revolution::make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());
            const auto           law_cell = m_law.getCellIVs(cell_i);
            const vector_dynamic stress_coeff = law_cell.projectStressOnCell(m_msh, cl, m_hdi, material_data);
            const auto           cell_nodes   = post_mesh.nodes_cell(cell_i);

            // Loop on the nodes of the cell
            for (auto& point_id : cell_nodes)
            {
                const auto pt = storage->points[point_id];

                const auto gphi   = gb.eval_functions(pt);
                const auto stress = revolution::eval(stress_coeff, gphi);

                value[point_id].first++;
                value[point_id].second += stress;
            }
            cell_i++;
        }

        std::vector<gmsh::Data>    data;    // create data
        std::vector<gmsh::SubData> subdata; // create subdata
        data.reserve(nb_nodes);             // data has a size of nb_node

        // Compute the average value and save it
        for (int i_node = 0; i_node < value.size(); i_node++)
        {
            const static_matrix<scalar_type, dimension, dimension> stress_avr =
              value[i_node].second / double(value[i_node].first);

            const gmsh::Data tmp_data(i_node + 1, disk::convertToVectorGmsh(stress_avr));
            data.push_back(tmp_data);
        }

        // Create and init a nodedata view
        gmsh::NodeData nodedata(9, 0.0, "stress_node_cont", data, subdata);
        // Save the view
        nodedata.saveNodeData(filename, gmsh);
    }

    void
    compute_equivalent_plastic_strain_GP(const std::string& filename) const
    {
        gmsh::Gmesh gmsh = disk::convertMesh(post_mesh);

        std::vector<gmsh::Data>    data;    // create data (not used)
        std::vector<gmsh::SubData> subdata; // create subdata to save soution at gauss point
        size_t                     nb_nodes(gmsh.getNumberofNodes());

        const auto material_data = m_law.getMaterialData();

        int cell_i = 0;
        for (auto& cl : m_msh)
        {

            const auto law_quadpoints = m_law.getCellIVs(cell_i).getIVs();

            // Loop on nodes
            for (auto& qp : law_quadpoints)
            {
                const auto                p   = qp.getAccumulatedPlasticStrain();
                const std::vector<double> p_s = disk::convertToVectorGmsh(p);

                // Add GP
                // Create a node at gauss point
                nb_nodes++;
                const gmsh::Node    new_node = disk::convertPoint(qp.point(), nb_nodes);
                const gmsh::SubData sdata(p_s, new_node);
                subdata.push_back(sdata); // add subdata
            }
            cell_i++;
        }

        // Save
        gmsh::NodeData nodedata(1, 0.0, "p_GP", data, subdata); // create and init a nodedata view

        nodedata.saveNodeData(filename, gmsh); // save the view
    }

    void
    compute_discontinuous_equivalent_plastic_strain(const std::string& filename) const
    {
        gmsh::Gmesh gmsh(dimension);
        auto        storage = m_msh.backend_storage();

        const auto grad_degree = m_hdi.grad_degree();

        std::vector<gmsh::Data>          data;    // create data (not used)
        const std::vector<gmsh::SubData> subdata; // create subdata to save soution at gauss point

        int    cell_i   = 0;
        size_t nb_nodes = 0;
        for (auto& cl : m_msh)
        {
            auto       pb       = revolution::make_scalar_monomial_basis(m_msh, cl, grad_degree);
            const auto law_cell = m_law.getCellIVs(cell_i);

            const vector_dynamic    p_coeff    = law_cell.projectPOnCell(m_msh, cl, m_hdi);
            const auto              cell_nodes = disk::cell_nodes(m_msh, cl);
            std::vector<gmsh::Node> new_nodes;

            // loop on the nodes of the cell
            for (int i = 0; i < cell_nodes.size(); i++)
            {
                nb_nodes++;
                const auto point_ids = cell_nodes[i];
                const auto pt        = storage->points[point_ids];

                const auto  pphi = pb.eval_functions(pt);
                scalar_type p    = revolution::eval(p_coeff, pphi);

                if (p <= scalar_type(0))
                    p = 0;

                const std::vector<double>   tens(1, p);
                const std::array<double, 3> coor = disk::init_coor(pt);

                // Add a node
                const gmsh::Node tmp_node(coor, nb_nodes, 0);
                new_nodes.push_back(tmp_node);
                gmsh.addNode(tmp_node);

                const gmsh::Data datatmp(nb_nodes, tens);
                data.push_back(datatmp);
            }

            // Add new element
            disk::add_element(gmsh, new_nodes);
            cell_i++;
        }

        // Create and init a nodedata view
        gmsh::NodeData nodedata(1, 0.0, "p_node_disc", data, subdata);

        // Save the view
        nodedata.saveNodeData(filename, gmsh);
    }

    void
    compute_continuous_equivalent_plastic_strain(const std::string& filename) const
    {
        const auto grad_degree = m_hdi.grad_degree();

        gmsh::Gmesh gmsh    = disk::convertMesh(post_mesh);
        auto        storage = post_mesh.mesh().backend_storage();

        const size_t nb_nodes(gmsh.getNumberofNodes());

        // first(number of data at this node), second(cumulated value)
        std::vector<std::pair<size_t, scalar_type>> value(nb_nodes, std::make_pair(0, 0.0));

        int cell_i = 0;
        for (auto& cl : m_msh)
        {
            auto       pb       = revolution::make_scalar_monomial_basis(m_msh, cl, grad_degree);
            const auto law_cell = m_law.getCellIVs(cell_i);

            const dynamic_vector<scalar_type> p_coeff    = law_cell.projectPOnCell(m_msh, cl, m_hdi);
            const auto                        cell_nodes = post_mesh.nodes_cell(cell_i);

            // Loop on the nodes of the cell
            for (auto& point_id : cell_nodes)
            {
                const auto pt = storage->points[point_id];

                const auto  pphi = pb.eval_functions(pt);
                scalar_type p    = revolution::eval(p_coeff, pphi);

                if (p <= scalar_type{0})
                    p = 0;

                // Add plastic at node
                value[point_id].first++;
                value[point_id].second += p;
            }
            cell_i++;
        }

        std::vector<gmsh::Data>    data;    // create data
        std::vector<gmsh::SubData> subdata; // create subdata
        data.reserve(nb_nodes);             // data has a size of nb_node

        // Compute the average value and save it
        for (int i_node = 0; i_node < value.size(); i_node++)
        {
            const scalar_type p_avr = value[i_node].second / double(value[i_node].first);

            const gmsh::Data tmp_data(i_node + 1, disk::convertToVectorGmsh(p_avr));
            data.push_back(tmp_data);
        }

        // Create and init a nodedata view
        gmsh::NodeData nodedata(1, 0.0, "p_node_const", data, subdata);

        // Save the view
        nodedata.saveNodeData(filename, gmsh);
    }

    void
    compute_is_plastic_continuous(const std::string& filename) const
    {
        const auto grad_degree = m_hdi.grad_degree();

        gmsh::Gmesh gmsh    = disk::convertMesh(post_mesh);
        auto        storage = post_mesh.mesh().backend_storage();

        const size_t nb_nodes(gmsh.getNumberofNodes());

        // first(number of data at this node), second(cumulated value)
        std::vector<std::pair<size_t, scalar_type>> value(nb_nodes, std::make_pair(0, 0.0));

        int cell_i = 0;
        for (auto& cl : m_msh)
        {
            auto       pb       = revolution::make_scalar_monomial_basis(m_msh, cl, grad_degree);
            const auto law_cell = m_law.getCellIVs(cell_i);

            const vector_dynamic state_coeff = law_cell.projectStateOnCell(m_msh, cl, m_hdi);
            const auto           cell_nodes  = post_mesh.nodes_cell(cell_i);

            // Loop on the nodes of the cell
            for (auto& point_id : cell_nodes)
            {
                const auto pt = storage->points[point_id];

                const auto  pphi = pb.eval_functions(pt);
                scalar_type p    = revolution::eval(state_coeff, pphi);

                if (p <= scalar_type(0))
                {
                    p = 0;
                }
                else if (p >= scalar_type(1))
                {
                    p = 1;
                }

                // Add plastic at node
                value[point_id].first++;
                value[point_id].second += p;
            }
            cell_i++;
        }

        std::vector<gmsh::Data>    data;    // create data
        std::vector<gmsh::SubData> subdata; // create subdata
        data.reserve(nb_nodes);             // data has a size of nb_node

        // Compute the average value and save it
        for (int i_node = 0; i_node < value.size(); i_node++)
        {
            const scalar_type p_avr = value[i_node].second / double(value[i_node].first);

            const gmsh::Data tmp_data(i_node + 1, disk::convertToVectorGmsh(p_avr));
            data.push_back(tmp_data);
        }

        // Create and init a nodedata view
        gmsh::NodeData nodedata(1, 0.0, "state_node_const", data, subdata);

        // Save the view
        nodedata.saveNodeData(filename, gmsh);
    }

    void
    compute_is_plastic_GP(const std::string& filename) const
    {
        gmsh::Gmesh gmsh = disk::convertMesh(post_mesh);

        std::vector<gmsh::Data>    data;    // create data (not used)
        std::vector<gmsh::SubData> subdata; // create subdata to save soution at gauss point
        size_t                     nb_nodes(gmsh.getNumberofNodes());

        const auto material_data = m_law.getMaterialData();

        int cell_i = 0;
        for (auto& cl : m_msh)
        {

            const auto law_quadpoints = m_law.getCellIVs(cell_i).getIVs();

            // Loop on nodes
            for (auto& qp : law_quadpoints)
            {
                scalar_type p = 0;
                if (qp.is_plastic())
                    p = 1;

                const std::vector<double> p_s = disk::convertToVectorGmsh(p);

                // Add GP
                // Create a node at gauss point
                nb_nodes++;
                const gmsh::Node    new_node = disk::convertPoint(qp.point(), nb_nodes);
                const gmsh::SubData sdata(p_s, new_node);
                subdata.push_back(sdata); // add subdata
            }
            cell_i++;
        }

        // Save
        gmsh::NodeData nodedata(1, 0.0, "state_GP", data, subdata); // create and init a nodedata view

        nodedata.saveNodeData(filename, gmsh); // save the view
    }

    void
    compute_continuous_deformed(const std::string& filename) const
    {
        const auto cell_degree = m_hdi.cell_degree();

        gmsh::Gmesh gmsh(dimension);
        auto        storage = m_msh.backend_storage();

        const static_vector<scalar_type, dimension> vzero = static_vector<scalar_type, dimension>::Zero();
        const size_t                                nb_nodes(m_msh.points_size());

        // first(number of data at this node), second(cumulated value)
        std::vector<std::pair<size_t, static_vector<scalar_type, dimension>>> value(nb_nodes, std::make_pair(0, vzero));

        int cell_i = 0;
        for (auto& cl : m_msh)
        {
            auto                 cb         = revolution::make_vector_monomial_basis(m_msh, cl, cell_degree);
            const vector_dynamic x          = m_solution_cells.at(cell_i++);
            const auto           cell_nodes = disk::cell_nodes(m_msh, cl);

            // Loop on the nodes of the cell
            for (int i = 0; i < cell_nodes.size(); i++)
            {
                const auto point_ids = cell_nodes[i];
                const auto pt        = storage->points[point_ids];

                const auto phi  = cb.eval_functions(pt);
                const auto depl = revolution::eval(x, phi);

                // Add displacement at node
                value[point_ids].first++;
                value[point_ids].second += depl;
            }
        }

        // New coordinate
        int i_node = 0;
        for (auto itor = m_msh.points_begin(); itor != m_msh.points_end(); itor++)
        {
            const auto            pt   = *itor;
            std::array<double, 3> coor = disk::init_coor(pt);

            const static_vector<scalar_type, dimension> depl_avr = value[i_node].second / double(value[i_node].first);

            for (int j = 0; j < dimension; j++)
                coor[j] += depl_avr(j);

            i_node++;
            const gmsh::Node tmp_node(coor, i_node, 0);
            gmsh.addNode(tmp_node);
        }
        const auto Nodes = gmsh.getNodes();

        // Add new elements
        for (auto& cl : m_msh)
        {
            const auto              cell_nodes = disk::cell_nodes(m_msh, cl);
            std::vector<gmsh::Node> new_nodes;

            // Loop on nodes of the cell
            for (int i = 0; i < cell_nodes.size(); i++)
            {
                const auto point_ids = cell_nodes[i];

                new_nodes.push_back(Nodes[point_ids]);
            }
            // Add new element
            disk::add_element(gmsh, new_nodes);
        }
        // Save mesh
        gmsh.writeGmesh(filename, 2);
    }

    void
    compute_discontinuous_deformed(const std::string& filename) const
    {
        const auto cell_degree = m_hdi.cell_degree();

        gmsh::Gmesh gmsh(dimension);
        auto        storage = m_msh.backend_storage();

        int    cell_i   = 0;
        size_t nb_nodes = 0;
        for (auto& cl : m_msh)
        {
            auto                    cb         = revolution::make_vector_monomial_basis(m_msh, cl, cell_degree);
            const vector_dynamic    x          = m_solution_cells.at(cell_i++);
            const auto              cell_nodes = disk::cell_nodes(m_msh, cl);
            std::vector<gmsh::Node> new_nodes;

            // Loop on nodes of the cell
            for (int i = 0; i < cell_nodes.size(); i++)
            {
                nb_nodes++;
                const auto point_ids = cell_nodes[i];
                const auto pt        = storage->points[point_ids];

                const auto phi  = cb.eval_functions(pt);
                const auto depl = revolution::eval(x, phi);

                std::array<double, 3> coor = disk::init_coor(pt);
                // Compute new coordinates
                for (int j = 0; j < dimension; j++)
                    coor[j] += depl(j);

                // Save node
                const gmsh::Node tmp_node(coor, nb_nodes, 0);
                new_nodes.push_back(tmp_node);
                gmsh.addNode(tmp_node);
            }
            // Add new element
            disk::add_element(gmsh, new_nodes);
        }
        // Save mesh
        gmsh.writeGmesh(filename, 2);
    }
};
