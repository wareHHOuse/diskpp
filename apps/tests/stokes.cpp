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
 * Karol Cascavita (C) 2018                     klcascavitam@unal.edu.co
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

#include "revolution/bases"
#include "revolution/quadratures"
#include "revolution/methods/hho"

#include "core/loaders/loader.hpp"

#include "solvers/solver.hpp"


template<typename Mesh>
typename Mesh::scalar_type
run_stokes(const Mesh& msh, size_t degree)
{

    typedef Mesh mesh_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;
    typedef typename mesh_type::scalar_type scalar_type;

    typedef dynamic_matrix<scalar_type>     matrix_type;
    typedef disk::mechanics::BoundaryConditions<mesh_type> boundary_type;

    using point_type = typename mesh_type::point_type;

    auto rhs_fun = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
        Matrix<scalar_type, 2, 1> ret;

        scalar_type x1 = p.x();
        scalar_type x2 = std::pow(p.x(), 2.);
        scalar_type x3 = std::pow(p.x(), 3.);
        scalar_type x4 = std::pow(p.x(), 4.);
        scalar_type y1 = p.y();
        scalar_type y2 = std::pow(p.y(), 2.);
        scalar_type y3 = std::pow(p.y(), 3.);
        scalar_type y4 = std::pow(p.y(), 4.);

        ret(0) = -(12.* x2 - 12.* x1 + 2.) * ( 4. * y3 - 6. * y2 + 2.* y1 )
                 -(x4 - 2. * x3 + x2 ) * (24. * y1 - 12.) + 5.* x4;

        ret(1) = +(12.* y2 - 12.* y1 + 2.) * ( 4. * x3 - 6. * x2 + 2.* x1 )
                 +(y4 - 2. * y3 + y2 ) * (24. * x1 - 12.) +  5.* y4;

        return ret;
    };

    auto sol_fun = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
        Matrix<scalar_type, 2, 1> ret;

        scalar_type x1 = p.x();
        scalar_type x2 = std::pow(p.x(), 2.);
        scalar_type x3 = std::pow(p.x(), 3.);
        scalar_type x4 = std::pow(p.x(), 4.);
        scalar_type y1 = p.y();
        scalar_type y2 = std::pow(p.y(), 2.);
        scalar_type y3 = std::pow(p.y(), 3.);
        scalar_type y4 = std::pow(p.y(), 4.);

        ret(0) =  (x4 - 2. * x3 + x2)  * ( 4. * y3 - 6. * y2 + 2.* y1 );
        ret(1) = -(y4 - 2. * y3 + y2 ) * ( 4. * x3 - 6. * x2 + 2.* x1 );

        return ret;
    };

    typename revolution::hho_degree_info hdi(degree);
    boundary_type bnd(msh);
    bnd.addDirichletEverywhere(sol_fun);

    auto assembler = revolution::make_stokes_assembler(msh, hdi, bnd);

    for (auto cl : msh)
    {
        auto gr = make_hho_vector_laplacian(msh, cl, hdi);

        Matrix<scalar_type, Dynamic, Dynamic> stab;
        stab = make_hho_fancy_stabilization_vector(msh, cl, gr.first, hdi);

        auto dr = make_hho_divergence_reconstruction_stokes_rhs(msh, cl, hdi);

        auto face_basis = revolution::make_vector_monomial_basis(msh, cl, hdi.face_degree());

        auto rhs = make_rhs(msh, cl, face_basis, rhs_fun);

        assembler.assemble(msh, cl, (gr.second + stab), -dr, rhs);
    }

    assembler.finalize();

    //dump_sparse_matrix(assembler.LHS, "stokes.txt");

    size_t systsz = assembler.LHS.rows();
    size_t nnz = assembler.LHS.nonZeros();

    dynamic_vector<scalar_type> sol = dynamic_vector<scalar_type>::Zero(systsz);

    disk::solvers::pardiso_params<scalar_type> pparams;
    mkl_pardiso_ldlt(pparams, assembler.LHS, assembler.RHS, sol);

    //std::ofstream ofs("velocity.dat");

    scalar_type error = 0.0;

    for (auto& cl : msh)
    {
    	auto bar = barycenter(msh, cl);

    	Matrix<scalar_type, Dynamic, 1> p = project_function(msh, cl, hdi, sol_fun);

    	auto cbs = revolution::vector_basis_size(degree, Mesh::dimension, Mesh::dimension);

    	auto cell_ofs = revolution::priv::offset(msh, cl);
    	Matrix<scalar_type, Dynamic, 1> s = sol.block(cell_ofs * cbs, 0, cbs, 1);

    	Matrix<scalar_type, Dynamic, 1> diff = s - p.head(cbs);

    	auto cb = revolution::make_vector_monomial_basis(msh, cl, hdi.cell_degree());

    	Matrix<scalar_type, Dynamic, Dynamic> mm = revolution::make_mass_matrix(msh, cl, cb);

    	error += diff.dot(mm*diff);

    	//ofs << bar.x() << " " << bar.y() << " " << s(0) << " " << s(1) << std::endl;
    }

    //ofs.close();

    return std::sqrt(error);

}

void convergence_test_typ1(void)
{
    using T = double;

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

    for (size_t k = 0; k < 5; k++)
    {
        std::cout << "DEGREE " << k << std::endl;

        std::vector<T> mesh_hs;
        std::vector<T> l2_velocity_errors;

        for (size_t i = 0; i < meshfiles.size(); i++)
        {
            typedef disk::generic_mesh<T, 2>  mesh_type;

            mesh_type msh;
            disk::fvca5_mesh_loader<T, 2> loader;
            if (!loader.read_mesh(meshfiles.at(i)))
            {
                std::cout << "Problem loading mesh." << std::endl;
                continue;
            }
            loader.populate_mesh(msh);

            auto error = run_stokes(msh, k);
            mesh_hs.push_back( disk::mesh_h(msh) );
            l2_velocity_errors.push_back(error);
        }

        for (size_t i = 0; i < mesh_hs.size(); i++)
        {
            if (i == 0)
            {
                std::cout << "    ";
                std::cout << std::scientific << std::setprecision(5) << mesh_hs.at(i) << "    ";
                std::cout << std::scientific << std::setprecision(5) << l2_velocity_errors.at(i);
                std::cout << "     -- " << std::endl;
            }
            else
            {
                auto rate = std::log( l2_velocity_errors.at(i)/l2_velocity_errors.at(i-1) ) /
                            std::log( mesh_hs.at(i)/mesh_hs.at(i-1) );
                std::cout << "    ";
                std::cout << std::scientific << std::setprecision(5) << mesh_hs.at(i) << "    ";
                std::cout << std::scientific << std::setprecision(5) << l2_velocity_errors.at(i) << "    ";
                std::cout << std::defaultfloat << std::setprecision(3) << rate << std::endl;
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
