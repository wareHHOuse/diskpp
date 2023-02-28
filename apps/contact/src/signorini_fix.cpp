/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Karol Cascavita (C) 2018                     karol.cascavita@enpc.fr
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
 #include <iomanip>
 #include <regex>
 #include <unistd.h>

 #include "diskpp/methods/hho"
 #include "diskpp/loaders/loader.hpp"
 #include "diskpp/output/silo.hpp"
 #include "diskpp/solvers/solver.hpp"
 #include "common.hpp"

template<typename Mesh, typename T, typename Function>
auto
fix_point_solver_faces(const Mesh& msh,
                const disk::scalar_boundary_conditions<Mesh>& bnd,
                const algorithm_parameters<T>& ap,
                const Function& rhs_fun)
{
    hho_degree_info hdi(ap.degree); // allow also (degree + 1, degree)

    std::cout << " * cell degree :"<< hdi.cell_degree() << std::endl;
    std::cout << " * face degree :"<< hdi.face_degree() << std::endl;

    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    auto is_contact_vector = make_is_contact_vector(msh, bnd);
    auto max_iter = 10000;
    auto tol = 1.e-6;

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs * msh.cells_size() + 2 * fbs * msh.faces_size()
                                        - fbs * msh.boundary_faces_size() ;
    disk::dynamic_vector<T> full_sol_old = disk::dynamic_vector<T>::Zero(num_full_dofs);

    auto offset_vector = full_offset(msh, hdi);

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        int cl_count = 0;
        auto assembler = make_diffusion_assembler(msh, hdi, bnd);

        for (auto& cl : msh)
        {
            auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
            auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);

            matrix_type Ah  = gr.second + stab;
            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());

            vector_type Bnegative  = vector_type::Zero(num_total_dofs);
            matrix_type Anitsche   = matrix_type::Zero(num_total_dofs, num_total_dofs);

            if (is_contact_vector.at(cl_count) == 1)
            {
                auto cell_ofs = offset_vector.at(cl_count);
                vector_type u_full = full_sol_old.block(cell_ofs, 0, num_total_dofs, 1);

            	Anitsche  = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                Bnegative = make_hho_negative_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);
            }

            matrix_type A = Ah - Anitsche;
            vector_type b = - Bnegative;
            b.block(0, 0, cbs, 1) += Lh;

            auto sc = make_scalar_static_condensation(msh, cl, hdi, A, b);
            assembler.assemble(msh, cl, bnd, sc.first, sc.second);
            cl_count++;
        }

        assembler.impose_neumann_boundary_conditions(msh, bnd);
        assembler.finalize();

        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        //std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
        //std::cout << "Mesh faces: " << msh.faces_size() << std::endl;
        //std::cout << "Dofs: " << systsz << std::endl;

        disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(systsz);

        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = true;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

        T error = 0.0 ;
        disk::dynamic_vector<T> full_sol = disk::dynamic_vector<T>::Zero(num_full_dofs);

        cl_count = 0;

        for (auto& cl : msh)
        {
            auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
            auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);

            matrix_type Ah  = gr.second + stab;
            vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());

            vector_type Bnegative  = vector_type::Zero(num_total_dofs);
            matrix_type Anitsche   = matrix_type::Zero(num_total_dofs, num_total_dofs);

            auto cell_ofs = offset_vector.at(cl_count);
            vector_type u_full_old = full_sol_old.block(cell_ofs, 0, num_total_dofs, 1);

            if (is_contact_vector.at(cl_count))
            {
                Anitsche  = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                Bnegative = make_hho_negative_faces(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full_old);
            }

            matrix_type A = Ah - Anitsche;
            vector_type cell_rhs = Lh - Bnegative.block(0, 0, cbs, 1);

            vector_type u_faces  = assembler.take_local_data(msh, cl, bnd, sol);
            vector_type u_full   = make_scalar_static_decondensation(msh, cl, hdi, A, cell_rhs, u_faces);

            full_sol.block(cell_ofs, 0, num_total_dofs, 1) = u_full;
            vector_type diff = u_full - u_full_old;

            error += diff.dot(A * diff);

            auto bar = barycenter(msh, cl);
            cl_count++;
        }

        if(iter % 100 == 0)
            std::cout << "  "<< iter << "  "<< std::sqrt(error)<< std::endl;

        if( std::sqrt(error)  < tol)
        {
            std::ofstream efs("solution_fhho_faces.dat");

            if(!efs.is_open())
                std::cout<< "Error opening file"<<std::endl;

            auto cl_count = 0;
            for(auto& cl : msh)
            {
                auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
                auto cell_ofs = offset_vector.at(cl_count++);
                vector_type u_bar = full_sol.block(cell_ofs, 0, 1, 1);
                auto bar = barycenter(msh, cl);
                efs << bar.x() << " " << bar.y() <<" "<< u_bar(0) << std::endl;
            }

            efs.close();

            std::cout << "  "<< iter << "  "<< std::sqrt(error)<< std::endl;
            return 0;
        }

        full_sol_old = full_sol;
    }
    return 1;
}

template<typename Mesh, typename T, typename Function>
auto
fix_point_solver_cells(const Mesh& msh,
                const disk::scalar_boundary_conditions<Mesh>& bnd,
                const algorithm_parameters<T>& ap,
                const Function& rhs_fun)
{
    hho_degree_info hdi(ap.degree + 1, ap.degree); // Not allow (degree + 1, degree)

    std::cout << " * cell degree :"<< hdi.cell_degree() << std::endl;
    std::cout << " * face degree :"<< hdi.face_degree() << std::endl;

    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    auto is_contact_vector = make_is_contact_vector(msh, bnd);

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension - 1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                        - fbs*msh.boundary_faces_size() ;
    disk::dynamic_vector<T> full_sol_old = disk::dynamic_vector<T>::Zero(num_full_dofs);

    auto offset_vector = full_offset(msh, hdi);
    auto max_iter = 10000;
    auto tol = 1.e-8;

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        int cl_count = 0;
        auto assembler = make_diffusion_assembler(msh, hdi, bnd);

        for (auto& cl : msh)
        {
            auto cb  = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());

            if (is_contact_vector.at(cl_count) == 1)
            {
                auto cell_ofs = offset_vector.at(cl_count);
                auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
                vector_type u_full = full_sol_old.block(cell_ofs, 0, num_total_dofs, 1);

                auto gr   = make_hho_contact_scalar_laplacian(msh, cl, hdi, bnd);
                auto stab = make_scalar_hdg_stabilization(msh, cl, hdi);
                //auto gr    = make_scalar_hho_laplacian(msh, cl, hdi);
                //auto stab  = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);

                matrix_type Ah = gr.second + stab;
                vector_type Lh = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());

            	matrix_type Anitsche  = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                vector_type Bnegative = make_hho_negative(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full);

                matrix_type A = Ah - Anitsche;
                vector_type b = - Bnegative;
                b.block(0, 0, cbs, 1) += Lh;

                const auto sc = make_scalar_static_condensation(msh, cl, hdi, A, b);
                assembler.assemble(msh, cl, bnd, sc.first, sc.second);
            }
            else
            {
                auto gr   = make_scalar_hho_laplacian(msh, cl, hdi);
                auto stab = make_scalar_hdg_stabilization(msh, cl, hdi);
                //auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);

                vector_type Lh = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
                matrix_type Ah = gr.second + stab;

                const auto sc = make_scalar_static_condensation(msh, cl, hdi, Ah, Lh);
                assembler.assemble(msh, cl, bnd, sc.first, sc.second);
            }
            cl_count++;
        }

        assembler.impose_neumann_boundary_conditions(msh, bnd);
        assembler.finalize();

        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        //std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
        //std::cout << "Mesh faces: " << msh.faces_size() << std::endl;
        //std::cout << "Dofs: " << systsz << std::endl;

        disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(systsz);

        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = true;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

        T error = 0.0 ;
        cl_count = 0;
        disk::dynamic_vector<T> full_sol = disk::dynamic_vector<T>::Zero(num_full_dofs);

        for (auto& cl : msh)
        {
            auto cb  = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto cell_ofs = offset_vector.at(cl_count);
            auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;

            vector_type u_full_old = full_sol_old.block(cell_ofs, 0, num_total_dofs, 1);

            if (is_contact_vector.at(cl_count))
            {
                auto gr   = make_hho_contact_scalar_laplacian(msh, cl, hdi, bnd);
                auto stab = make_scalar_hdg_stabilization(msh, cl, hdi);

                //auto gr   = make_scalar_hho_laplacian(msh, cl, hdi);
                //auto stab = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);

                matrix_type Ah = gr.second + stab;
                vector_type Lh = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());

                matrix_type Anitsche  = make_hho_nitsche(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd );
                vector_type Bnegative = make_hho_negative(msh, cl, hdi, gr.first, ap.gamma_0, ap.theta, bnd, u_full_old);

                matrix_type A = Ah - Anitsche;
                vector_type cell_rhs = Lh - Bnegative.block(0, 0, cbs, 1);

                vector_type u_faces  = assembler.take_local_data(msh, cl, bnd, sol);
                vector_type u_full   = make_scalar_static_decondensation(msh, cl, hdi, A, cell_rhs, u_faces);

                vector_type diff = u_full - u_full_old;
                error += diff.dot(A * diff);
                full_sol.block(cell_ofs, 0, num_total_dofs, 1) = u_full;
            }
            else
            {
                auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
                auto stab   = make_scalar_hdg_stabilization(msh, cl, hdi);
                //auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);

                matrix_type Ah  = gr.second + stab;
                vector_type Lh  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());

                vector_type u_faces = assembler.take_local_data(msh, cl, bnd, sol);
                vector_type u_full  = make_scalar_static_decondensation(msh, cl, hdi, Ah, Lh, u_faces);

                vector_type diff = u_full - u_full_old;
                error += diff.dot(Ah * diff);

                full_sol.block(cell_ofs, 0, num_total_dofs, 1) = u_full;
            }
            cl_count++;
        }

        if(iter % 100 == 0)
            std::cout << "  "<< iter << "  "<< std::sqrt(error)<< std::endl;

        if( std::sqrt(error)  < tol)
        {
            std::ofstream efs("solution_fhho_cells.dat");

            if(!efs.is_open())
                std::cout<< "Error opening file"<<std::endl;

            auto cl_count = 0;
            for(auto& cl : msh)
            {
                auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
                auto cell_ofs = offset_vector.at(cl_count++);
                vector_type u_bar = full_sol.block(cell_ofs, 0, 1, 1);
                auto bar = barycenter(msh, cl);
                efs << bar.x() << " " << bar.y() <<" "<< u_bar(0) << std::endl;
            }

            efs.close();

            std::cout << "  "<< iter << "  "<< std::sqrt(error)<< std::endl;

            return 0;
        }

        full_sol_old = full_sol;
    }
    return 1;
}


template<typename Mesh, typename T>
void
run_signorini(  const Mesh& msh, const algorithm_parameters<T>& ap)
{
    typedef typename Mesh::point_type  point_type;

    auto zero_fun = [](const point_type& p) -> T {
        return 0;
    };

    auto rhs_fun = [](const point_type& p) -> T {
        return - 2.* M_PI *  std::sin(2. * M_PI * p.x());
    };

    disk::scalar_boundary_conditions<Mesh>     bnd(msh);

    bnd.addDirichletBC(disk::DIRICHLET,1,zero_fun); //TOP
    bnd.addNeumannBC(disk::NEUMANN,  2,  zero_fun); //
    bnd.addNeumannBC(disk::NEUMANN,  4,  zero_fun); //
    bnd.addContactBC(disk::SIGNORINI,3); //BOTTOM

    switch(ap.solver)
    {
        case EVAL_ON_FACES:
            fix_point_solver_faces( msh, bnd, ap, rhs_fun);
            break;
        case EVAL_IN_CELLS:
            fix_point_solver_cells( msh, bnd, ap, rhs_fun);
            break;
        default:
            throw std::invalid_argument("Invalid solver");
    }

    return;
}

int main(int argc, char **argv)
{
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    char    *filename       = nullptr;
    using T = double;

    int ch;
    algorithm_parameters<T> ap;
    int degree = 1;

    while ( (ch = getopt(argc, argv, "k:g:npzfc")) != -1 )
    {
        switch(ch)
        {
            case 'k':
                degree = atoi(optarg);
                if (degree < 0)
                {
                    std::cout << "Degree must be positive. Falling back to 1." << std::endl;
                    degree = 1;
                }
                ap.degree = degree;
                break;
            case 'g':
                ap.gamma_0 = atof(optarg);
                if (ap.gamma_0 <= 0)
                {
                    std::cout << "gamma_0 must be > 0. Falling back to 0.1" << std::endl;
                    ap.gamma_0 = 0.1;
                }
                break;
            case 'n':
                ap.theta = -1.;
                break;
            case 'p':
                ap.theta = 1.;
                break;
            case 'z':
                ap.theta = 0.;
                break;
            case 'f':
                ap.solver = EVAL_ON_FACES;
                break;
            case 'c':
                ap.solver = EVAL_IN_CELLS;
                break;
            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }
    std::cout << ap << std::endl;

    argc -= optind;
    argv += optind;

    filename = argv[0];

    /* Netgen 2d*/
    if (std::regex_match(filename, std::regex(".*\\.mesh2d$") ))
    {
        typedef disk::simplicial_mesh<T, 2>  mesh_type;
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        mesh_type msh;
        disk::netgen_mesh_loader<T, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        run_signorini(msh, ap);
    }
    #if 0
    /* Medit 2d*/
    if (std::regex_match(filename, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        typedef disk::generic_mesh<T, 2>  mesh_type;
        mesh_type msh = disk::load_medit_2d_mesh<T>(filename);
        run_signorini(msh, ap);
    }
    #endif

    return 0;
}
