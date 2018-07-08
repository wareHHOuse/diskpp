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

 #include "revolution/bases"
 #include "revolution/quadratures"
 #include "revolution/methods/hho"

 #include "core/loaders/loader.hpp"

 #include "output/silo.hpp"
 #include "common.hpp"

enum method_type
{
    POSITIVE,
    NEGATIVE,
    ZERO
};

template<typename Mesh>
class newton_solver
{

    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    matrix_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>                 vector_type;

    using T = typename Mesh::scalar_type;
    size_t cbs, fcs;
    hho_degree_info hdi;

    template<typename Function>
    auto
    solve(const mesh_type&  msh, const Function& rhs_fun)
    {
        auto cl_count = 0;

        auto assembler = make_diffusion_assembler2(msh, hdi, bnd);
        for (auto& cl : msh)
        {
            auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
            auto rhs_force  = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            matrix_type A   = gr.second + stab;
            vector_type rhs = vector_type::Zero(num_total_dofs);
            rhs.block(0, 0, cbs, 1) = rhs_force;

            auto cell_ofs = offset_vector.at(cl_count);
            vector_type u_full_old = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            if (is_contact_vector.at(cl_count))
            {
                //matrix An(du, v)
            	A -= make_nitsche(msh, cl, hdi, gr.first, gamma_0, theta, bnd );

                //rhs: int [Pn(u_old,v)]_Pn[v]
                rhs -= make_negative(msh, cl, hdi, gr.first, gamma_0, theta, bnd, u_full_old);
            }
            //rhs: An(u_old,v)
            rhs -= A * u_full_old;

            if (is_contact_vector.at(cl_count))
               A -=  make_heaviside(msh, cl, hdi, gr.first, gamma_0, theta, bnd, u_full_old);

            auto sc = diffusion_static_condensation_compute_full(msh, cl, hdi, A, rhs);
            assembler.assemble(msh, cl, sc.first, sc.second);
            cl_count++;
        }

        assembler.impose_neumann_boundary_conditions(msh, bnd);
        assembler.finalize();

        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        //std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
        //std::cout << "Mesh faces: " << msh.faces_size() << std::endl;
        //std::cout << "Dofs: " << systsz << std::endl;

        dynamic_vector<T> dsol = dynamic_vector<T>::Zero(systsz);

        disk::solvers::pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = true;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, dsol);
        return dsol;
    }


    template<typename Function>
    auto
    recover(const mesh_type& msh, const dynamic_vector<T>& dsol, const Function& rhs_fun)
    {
        T error = 0.0 ;
        std::ofstream ofs("sol_"+ tostr(iter) +".dat");
        if(!ofs.is_open())
            std::cout<< "Error opening file"<<std::endl;

        auto cl_count = 0;
        for (auto& cl : msh)
        {
            auto num_total_dofs = cbs + howmany_faces(msh,cl) * fbs;
            auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto gr     = make_hho_scalar_laplacian(msh, cl, hdi);
            auto stab   = make_hho_scalar_stabilization(msh, cl, gr.first, hdi);
            auto rhs_force = make_rhs(msh, cl, cb, rhs_fun, hdi.cell_degree());
            matrix_type A   = gr.second + stab;

            vector_type rhs = vector_type::Zero(num_total_dofs);
            rhs.block(0, 0, cbs, 1) = rhs_force;
            auto cell_ofs = offset_vector.at(cl_count);
            vector_type u_full_old = full_sol.block(cell_ofs, 0, num_total_dofs, 1);

            if (is_contact_vector.at(cl_count))
            {
                //matrix An(du, v)
            	A -= make_nitsche(msh, cl, hdi, gr.first, gamma_0, theta, bnd );

                //rhs: int [Pn(u_old,v)]_Pn[v]
                rhs -= make_negative(msh, cl, hdi, gr.first, gamma_0, theta, bnd, u_full_old);
            }
            //rhs: An(u_old,v)
            rhs -= A * u_full_old;

            //matrix Heaviside term: Dont use this before rhs -= A * u_full;
            if (is_contact_vector.at(cl_count))
               A -=  make_heaviside(msh, cl, hdi, gr.first, gamma_0, theta, bnd, u_full_old);

            vector_type cell_rhs = rhs.block(0, 0, cbs, 1);
            vector_type du_faces = assembler.take_local_data(msh, cl, dsol);
            vector_type du_full  =
                diffusion_static_condensation_recover(msh, cl, hdi, A, cell_rhs, du_faces);
            vector_type  diff = du_full;
            error += diff.dot(A * diff);

            full_sol.block(cell_ofs, 0, num_total_dofs ,1) += du_full;
            auto bar = barycenter(msh, cl);
            ofs << bar.x() << " " << bar.y() <<" "<< du_full(0) << std::endl;
            cl_count++;
        }
        ofs.close();

        std::cout << "  "<< iter << "  "<< std::sqrt(error)<< std::endl;
        if( std::sqrt(error)  < tol)
            return 0;
    }
public:
    auto
    newton_solver(const Mesh& msh, const hho_degree_info& hdi,
                    const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
                    const algorithm_parameters<T>& ap,
                    const Function& rhs_fun,
                    dynamic_vector<T>& full_sol)
    {
        std::cout << "Solver : Generalized Newton" << std::endl;
        std::cout << "Theta  : "<< theta << std::endl;
        std::cout << "Gamma0 : "<< gamma_0 << std::endl;

        auto is_contact_vector = make_is_contact_vector(msh, bnd);

        fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
        cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

        auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                        - fbs*msh.boundary_faces_size() ;


        offset_vector = full_offset(msh, hdi);

        return;
    }

    run_hho_newton()
    {
        auto max_iter = 1000;
        auto tol = 1.e-6;

        for(size_t iter = 0; iter < max_iter; iter++)
        {
            solve();
            recover();
        }
        return;
    }
}


template<typename Mesh, typename T>
void
run_signorini(  const Mesh& msh,
                algorithm_parameters<T>& ap,
                const method_type& tt,
                const solver_type& ss)
{
    typedef typename mesh_type::point_type  point_type;

    const size_t degree = 0 ;
    hho_degree_info hdi(degree, degree);

    auto zero_fun = [](const point_type& p) -> T {
        return 0;
    }

    auto rhs_fun = [](const point_type& p) -> T {
        return - 2.* M_PI *  std::sin(2. * M_PI * pt.x());
    }

    typedef disk::mechanics::BoundaryConditionsScalar<Mesh> boundary_type;
    boundary_type  m_bnd(msh);

    m_bnd.addDirichletBC(disk::mechanics::DIRICHLET,1,zero_fun); //TOP
    m_bnd.addNeumannBC(disk::mechanics::NEUMANN,3,zero_fun); //
    m_bnd.addNeumannBC(disk::mechanics::NEUMANN,4,zero_fun); //
    m_bnd.addContactBC(disk::mechanics::SIGNORINI,2); //BOTTOM

    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    auto num_full_dofs = cbs*msh.cells_size() + 2 * fbs*msh.faces_size()
                                        - fbs*msh.boundary_faces_size() ;

    dynamic_vector<T> full_sol = dynamic_vector<T>::Zero(num_full_dofs);

    newton_solver ns(msh, hdi, m_bnd, ap, rhs_fun);
    ns.run_hho_newton();

    //full_sol = fix_point_solver( msh, hdi, m_bnd, ap, rhs_fun);
    //newton_solver( msh, hdi, m_bnd, ap, rhs_fun, full_sol);

    return;
}

int main(int argc, char **argv)
{
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    char    *filename       = nullptr;
    using T = double;

    int ch;
    algorithm_parameters<T> ap;
    method_type tt = POSITIVE;
    solver_type ss = FIX_POINT;

    while ( (ch = getopt(argc, argv, "g:npzfwm")) != -1 )
    {
        switch(ch)
        {
            case 'g':
                std::cout << "choosing gamma" << std::endl;
                ap.gamma_0 = atof(optarg);
                if (ap.gamma_0 <= 0)
                {
                    std::cout << "gamma_0 must be >0. Falling back to 0.1" << std::endl;
                    ap.gamma_0 = 0.1;
                }
                break;
            case 'n':
                ap.theta = -1.;
                std::cout << "theta negative chosen" << std::endl;
                break;
            case 'p':
                ap.theta = 1.;
                std::cout << "theta positive chosen" << std::endl;
                break;
            case 'z':
                ap.theta = 0.;
                std::cout << "theta zero chosen" << std::endl;
                break;
            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    filename = argv[0];

    /* Medit 2d*/
    if (std::regex_match(filename, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        typedef disk::generic_mesh<T, 2>  mesh_type;
        mesh_type msh = disk::load_medit_2d_mesh<T>(filename);
        run_signorini(msh, ap, tt, ss);
    }

    return 0;
}
