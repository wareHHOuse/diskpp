/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
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

#include "bases/bases.hpp"
#include "quadratures/quadratures.hpp"
#include "methods/hho"

#include "core/loaders/loader.hpp"

#include "solvers/solver.hpp"

#include "output/silo.hpp"


template<typename Mesh>
class augmented_lagrangian_diffusion
{
    typedef Mesh mesh_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;
    typedef typename mesh_type::scalar_type scalar_type;

    typedef disk::mechanics::BoundaryConditions<mesh_type> boundary_type;

    using point_type = typename mesh_type::point_type;

    typedef dynamic_vector<scalar_type>       vector_type;
    typedef dynamic_matrix<scalar_type>       matrix_type;

    typedef Matrix<scalar_type, 2, 1>         tensor_type;
    typedef Matrix<scalar_type, 2, 1>         vector2d_type;
    typedef std::function<scalar_type   (const point_type &)> scalar_funtion_type;

    scalar_funtion_type     sol_fun, rhs_fun;
    vector_type             sol, sol_old;
    vector_type             multiplier;
    vector_type             auxiliar;
    Matrix<scalar_type, Dynamic, 1> rhs_all;
    typename disk::hho_degree_info di;
    scalar_type             viscosity;
    scalar_type             alpha;
    size_t                  cbs, fbs, sbs;

public:
    scalar_type             convergence;

    augmented_lagrangian_diffusion(const Mesh& msh,
                            const typename disk::hho_degree_info & hdi):
                            di(hdi)
    {
        rhs_fun = [](const point_type& pt) -> scalar_type
        {
            auto sin_px = std::sin(M_PI * pt.x());
            auto sin_py = std::sin(M_PI * pt.y());
            return 2.0 * M_PI * M_PI * sin_px * sin_py;
        };

        sol_fun = [](const point_type& pt) -> scalar_type
        {
            auto sin_px = std::sin(M_PI * pt.x());
            auto sin_py = std::sin(M_PI * pt.y());
            return sin_px * sin_py;
        };

        viscosity = 1.;
        alpha = 0.1;
        convergence = 0.;
        auto dim =  Mesh::dimension;

        cbs = disk::scalar_basis_size(di.cell_degree(), dim);
        fbs = disk::scalar_basis_size(di.face_degree(), dim - 1);
        sbs = disk::vector_basis_size(di.cell_degree(), dim, dim);
    };

    auto
    initialize(const mesh_type& msh)
    {
        size_t full_dofs = 0;

        for(auto & cl : msh)
            full_dofs += cbs + fbs * howmany_faces(msh, cl);

        sol = vector_type::Zero(full_dofs);
        sol_old = vector_type::Zero(full_dofs);

        rhs_all = Matrix<scalar_type, Dynamic, 1>::Zero(msh.cells_size() * cbs);

        multiplier = vector_type::Zero(msh.cells_size() * sbs);
        auxiliar  = vector_type::Zero(msh.cells_size() * sbs);

        return;
    }

    auto
    compute_errors( const mesh_type& msh)
    {
        auto dim =  Mesh::dimension;

        scalar_type error(0), error_vel(0), error_pres(0);

        for (auto& cl : msh)
        {
        	auto bar = barycenter(msh, cl);
        	Matrix<scalar_type, Dynamic, 1> p = project_function(msh, cl, di, sol_fun);
        	auto cell_ofs = disk::priv::offset(msh, cl);
        	Matrix<scalar_type, Dynamic, 1> s = sol.block(cell_ofs * cbs, 0, cbs, 1);
        	Matrix<scalar_type, Dynamic, 1> diff = s - p.head(cbs);
        	auto cb = disk::make_scalar_monomial_basis(msh, cl, di.cell_degree());
        	Matrix<scalar_type, Dynamic, Dynamic> mm = disk::make_mass_matrix(msh, cl, cb);
        	error += diff.dot(mm*diff);
        	//ofs << bar.x() << " " << bar.y() << " " << s(0) << " " << s(1) << std::endl;

            //energy error

            auto num_total_dofs = cbs + fbs * howmany_faces(msh, cl);
            Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> svel =
                sol.block(cell_ofs * num_total_dofs, 0, num_total_dofs, 1);

            Matrix<scalar_type, Dynamic, 1> diff_vel = svel - p;

            auto gr   = make_hho_scalar_laplacian(msh, cl, di);
            auto stab = make_hho_vector_stabilization(msh, cl, gr.first, di);
            auto G    = make_hho_gradrec_vector(msh, cl, di);

            error_vel += diff_vel.dot( viscosity * (G.second + stab)*diff_vel);
        }

        //ofs.close();
        return std::make_pair(std::sqrt(error), std::sqrt(error_pres));
    }

    void
    save_auxiliar(  const mesh_type& msh,
                        const cell_type& cl)
    {
        auto gamma = compute_auxiliar(msh, cl);
        auto cell_ofs  = disk::priv::offset(msh, cl);
        auxiliar.block(cell_ofs * sbs, 0, sbs, 1) = gamma;
        return;
    }

    Matrix<scalar_type, Dynamic, 1>
    compute_auxiliar(   const mesh_type& msh,
                        const cell_type& cl)
    {
        auto sb = disk::make_vector_monomial_basis(msh, cl, di.cell_degree());

        auto cell_ofs  = disk::priv::offset(msh, cl);
        auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;

        Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> u_TF =
            sol_old.block( cell_ofs * num_total_dofs, 0, num_total_dofs, 1);

        auto value = 1./ (viscosity + alpha);
        auto G     = make_hho_gradrec_vector(msh, cl, di);

        Matrix<scalar_type, Dynamic, 1> Gu = G.first * u_TF;
        Matrix<scalar_type, Dynamic, 1> stress = multiplier.block(cell_ofs * sbs, 0, sbs, 1);
        Matrix<scalar_type, Dynamic, 1> gamma = value * (stress +  alpha * Gu);

        return gamma;
    }

    auto
    update_multiplier(const mesh_type& msh)
    {
        convergence = 0.;
        auto dim = Mesh::dimension;

        for(auto cl: msh)
        {
            auto cell_ofs  = disk::priv::offset(msh, cl);
            auto num_total_dofs = cbs + howmany_faces(msh, cl) * fbs;

            Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> u_TF =
                sol.block(cell_ofs * num_total_dofs, 0, num_total_dofs, 1);

            auto                            G         = make_hho_gradrec_vector(msh, cl, di);
            Matrix<scalar_type, Dynamic, 1> Gu        = G.first * u_TF;
            Matrix<scalar_type, Dynamic, 1> gamma     = compute_auxiliar(msh, cl);
            Matrix<scalar_type, Dynamic, 1> gamma_old = auxiliar.block(cell_ofs * sbs, 0, sbs, 1);

            Matrix<scalar_type, Dynamic, 1> diff_stress  = alpha * (Gu - gamma);
            Matrix<scalar_type, Dynamic, 1> diff_gamma   = alpha * (gamma - gamma_old);

            auto sb =  disk::make_vector_monomial_basis(msh, cl, di.cell_degree());
            Matrix<scalar_type, Dynamic, Dynamic> mass = disk::make_mass_matrix(msh, cl, sb);

            convergence += diff_stress.dot(mass * diff_stress) + diff_gamma.dot(mass * diff_gamma);

            multiplier.block(cell_ofs * sbs, 0, sbs, 1) +=  alpha * diff_stress;
        }

        return;
    }

    Matrix<scalar_type, Dynamic, 1>
    make_rhs_alg(   const mesh_type& msh,
                    const cell_type& cl)
    {
        auto G         = make_hho_gradrec_vector(msh, cl, di);
        auto cb = disk::make_scalar_monomial_basis(msh, cl, di.cell_degree());
        auto sb = disk::make_vector_monomial_basis(msh, cl, di.cell_degree());
        auto cell_ofs  = disk::priv::offset(msh, cl);
        auto num_faces = howmany_faces(msh, cl);

        auto stress = multiplier.block( sbs * cell_ofs,  0, sbs, 1);
        auto gamma  = compute_auxiliar( msh,  cl);
        vector_type str_agam = stress - alpha * gamma;

        Matrix<scalar_type, Dynamic, Dynamic> mm = disk::make_mass_matrix(msh, cl, sb);

        Matrix<scalar_type, Dynamic, 1> rhs =
                    Matrix<scalar_type, Dynamic, 1>::Zero(cbs + fbs * num_faces);

        //(f, v_T)
        rhs.block( 0, 0, cbs, 1) = make_rhs(msh, cl, cb, rhs_fun);

        //(stress - factor * alpha * gamma, Gv)
        rhs -=  G.first.transpose() * mm * str_agam;

        rhs_all.block( cbs * cell_ofs, 0, cbs, 1) = rhs.block( 0, 0, cbs, 1);
        return rhs;
    }

    auto
    run_diffusion_like(const mesh_type& msh)
    {
        sol_old = sol;

        auto assembler = disk::make_diffusion_assembler(msh, di);
        for (auto& cl : msh)
        {
            auto cb = disk::make_scalar_monomial_basis(msh, cl, di.cell_degree());
            auto G    = make_hho_gradrec_vector(msh, cl, di);
            auto gr   = make_hho_scalar_laplacian(msh, cl, di);
            auto stab = make_hho_scalar_stabilization(msh, cl, gr.first, di);

            Matrix<scalar_type, Dynamic, 1> local_rhs =
                Matrix<scalar_type, Dynamic, 1>::Zero(cbs + fbs * howmany_faces(msh,cl));
            local_rhs = make_rhs_alg(msh, cl);

            Matrix<scalar_type, Dynamic, Dynamic> A = alpha * G.second + viscosity * stab;
            auto sc = diffusion_static_condensation_compute_alg(msh, cl, di, A, local_rhs);
            assembler.assemble(msh, cl, sc.first, sc.second, sol_fun);

            save_auxiliar(msh, cl);
        }

        assembler.finalize();

        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        //std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
        //std::cout << "Mesh faces: " << msh.faces_size() << std::endl;
        //std::cout << "Dofs: " << systsz << std::endl;

        dynamic_vector<scalar_type> sol_faces = dynamic_vector<scalar_type>::Zero(systsz);

        disk::solvers::pardiso_params<scalar_type> pparams;
        pparams.report_factorization_Mflops = false;
        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol_faces);

        scalar_type error = 0.0;

        //std::ofstream ofs("sol.dat");

        for (auto& cl : msh)
        {
            auto cb = disk::make_scalar_monomial_basis(msh, cl, di.cell_degree());
            auto G    = make_hho_gradrec_vector(msh, cl, di);
            auto gr   = make_hho_scalar_laplacian(msh, cl, di);
            auto stab = make_hho_scalar_stabilization(msh, cl, gr.first, di);

            auto cell_ofs  = disk::priv::offset(msh, cl);
            auto num_total_dofs = cbs + fbs * howmany_faces(msh, cl);

            Matrix<scalar_type, Dynamic, 1> cell_rhs = rhs_all.block( cbs*cell_ofs, 0, cbs, 1);

            Matrix<scalar_type, Dynamic, Dynamic> A = alpha*G.second + viscosity*stab;

            Matrix<scalar_type, Dynamic, 1> locsol =
                assembler.take_local_data(msh, cl, sol_faces, sol_fun);

            Matrix<scalar_type, Dynamic, 1> fullsol =
                diffusion_static_condensation_recover(msh, cl, di,  A, cell_rhs, locsol);

            sol.block(cell_ofs * num_total_dofs, 0, num_total_dofs, 1) = fullsol;

            Matrix<scalar_type, Dynamic, 1> realsol = project_function(msh, cl, di, sol_fun);

            auto diff = realsol - fullsol;

            Matrix<scalar_type, Dynamic, Dynamic> B = viscosity*G.second + viscosity*stab;

            error += diff.dot( B *diff);
        }

        return std::sqrt(error);
    }
    void
    quiver(  const mesh_type& msh, std::ofstream & ofs)
    {

        for (auto& cl: msh)
        {
            auto cell_ofs = disk::priv::offset(msh, cl);
            vector_type stress = multiplier.block(cell_ofs * sbs, 0, sbs, 1);
            auto gamma  = compute_auxiliar( msh,  cl);

            auto sb = disk::make_vector_monomial_basis(msh, cl, di.cell_degree());
            auto qps = disk::integrate(msh, cl, 2 * di.face_degree());
            for (auto& qp : qps)
            {
                auto phi = sb.eval_functions(qp.point());
                vector_type  st = stress.transpose() * phi;
                vector_type  gm = gamma.transpose() * phi;

                ofs << qp.point().x() << "  "<< qp.point().y()<< "   ";
                ofs << st(0) << "  "<< st(1) << "    ";
                ofs << gm(0) << "  "<< gm(1) << std::endl;
            }
        }
    }
};

template< typename T>
std::string
tostr(const T & var)
{
    std::ostringstream  ostr;
    ostr << var;
    return ostr.str();
}

template<typename Mesh>
auto
run_alg_diffusion(const Mesh& msh, size_t degree,
                    std::ofstream & ofs, std::ofstream & qfs)
{
    using T = typename Mesh::coordinate_type;
    T tolerance = 10.e-9, Ninf = 10.e+10;
    size_t max_iters = 10000;

    typename disk::hho_degree_info hdi(degree+1, degree);
    augmented_lagrangian_diffusion<Mesh> ald(msh, hdi);

    ald.initialize(msh);

    T error = 0.;
    size_t i = 0;

    for( i = 0; i < max_iters; i++)
    {
        error = ald.run_diffusion_like(msh);
        ald.update_multiplier(msh);
        auto convergence = std::sqrt(ald.convergence);
        ofs << i << " "<< error << "   "<<  convergence << std::endl;

        assert(convergence < Ninf);
        if(convergence < tolerance)
            break;
    }

    ald.quiver(msh, qfs);
    if(i == max_iters-1)
        std::cout << "Max_iters reached" << std::endl;
    //auto error = ald.compute_errors(msh, assembler);

    return error;
}

void convergence_test_typ1(void)
{
    using T = double;
    bool use_sym_grad = false;
    std::vector<std::string> meshfiles;
    /*
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_1.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_2.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_3.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_4.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_5.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_6.typ1");
    */

    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_1.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_2.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_3.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_4.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_5.typ1");

    /*
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_1.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_2.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_3.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_4.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_5.typ1");
    */
    std::cout << "                   velocity H1-error";
    std::cout << "    -     pressure L2-error "<< std::endl;

    for (size_t k = 0; k < 3; k++)
    {
        std::cout << "DEGREE " << k << std::endl;

        std::ofstream ofs("errors_k" + tostr(k) + ".data");
        std::ofstream qfs("quiver_k" + tostr(k) + ".data");

        if (!ofs.is_open())
            std::cout << "Error opening errors "<<std::endl;

        if (!qfs.is_open())
            std::cout << "Error opening quiver "<<std::endl;

        std::vector<T> mesh_hs;
        std::vector<T> errors;

        //for (size_t i = 0; i < meshfiles.size(); i++)
        for (size_t i = 0; i < 1; i++)
        {
            typedef disk::generic_mesh<T, 2>  mesh_type;

            std::cout << " Mesh : "<< i << std::endl;
            mesh_type msh;
            disk::fvca5_mesh_loader<T, 2> loader;
            if (!loader.read_mesh(meshfiles.at(i)))
            {
                std::cout << "Problem loading mesh." << std::endl;
                continue;
            }
            loader.populate_mesh(msh);
            auto error = run_alg_diffusion(msh, k, ofs, qfs);
            mesh_hs.push_back( disk::average_diameter(msh) );
            errors.push_back(error);

            ofs << " " << std::endl;
            ofs << " " << std::endl;
        }

        ofs.close();
        qfs.close();

        for (size_t i = 0; i < mesh_hs.size(); i++)
        {
            if (i == 0)
            {
                std::cout << "    ";
                std::cout << std::scientific << std::setprecision(4) << mesh_hs.at(i) << "    ";
                std::cout << std::scientific << std::setprecision(4) << errors.at(i) << std::endl;
            }
            else
            {
                auto rate = std::log( errors.at(i)/errors.at(i-1) ) /
                            std::log( mesh_hs.at(i)/mesh_hs.at(i-1) );
                std::cout << "    ";
                std::cout << std::scientific  << std::setprecision(4) << mesh_hs.at(i) << "    ";
                std::cout << std::scientific  << std::setprecision(4) << errors.at(i)  << "    ";
                std::cout << std::fixed<< std::setprecision(2) << rate <<  std::endl;
            }
        }
    }
}

int main(void)
{
    convergence_test_typ1();
}

#if 0

int main(int argc, char **argv)
{
    using RealType = double;

    char    *filename       = nullptr;
    int ch;
    size_t degree = 1;

    while ( (ch = getopt(argc, argv, "k:")) != -1 )
    {
        switch(ch)
        {
            case 'k':
                degree = atoi(optarg);
                break;

            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    if (argc != 1)
    {
        std::cout << "Please specify a 2D mesh" << std::endl;

        return 0;
    }

    filename = argv[0];

    if (std::regex_match(filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;

        typedef disk::generic_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::fvca5_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        run_stokes(msh, degree);
    }

    if (std::regex_match(filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;

        typedef disk::simplicial_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::netgen_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        run_stokes(msh, degree);
    }

/*
    if (std::regex_match(filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: Cartesian 2D" << std::endl;

        typedef disk::cartesian_mesh<RealType, 2>   mesh_type;

        mesh_type msh;
        disk::cartesian_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        test_bases(msh);
    }
    */

    return 0;
}
#endif
