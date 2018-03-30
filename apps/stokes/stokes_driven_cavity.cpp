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

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"
#include "output/postMesh.hpp"

#include "output/silo.hpp"

template<typename Mesh>
void
compute_discontinuous_velocity(const Mesh& msh,
                        const dynamic_vector< typename Mesh::scalar_type>& sol,
                        const typename revolution::hho_degree_info& hdi,
                        const std::string& filename)
{
    typedef Mesh mesh_type;
    typedef typename Mesh::scalar_type scalar_type;
    const size_t cell_degree   = hdi.cell_degree();
    const size_t cbs = revolution::vector_basis_size(cell_degree,
                                    Mesh::dimension, Mesh::dimension);
    // compute mesh for post-processing
    disk::PostMesh<mesh_type> post_mesh = disk::PostMesh<mesh_type>(msh);
    gmsh::Gmesh gmsh    = disk::convertMesh(post_mesh);
    auto        storage = post_mesh.mesh().backend_storage();

    const static_vector<scalar_type, Mesh::dimension> vzero =
                            static_vector<scalar_type, Mesh::dimension>::Zero();

    const size_t nb_nodes(gmsh.getNumberofNodes());

    // first(number of data at this node), second(cumulated value)
    std::vector<std::pair<size_t, static_vector<scalar_type, Mesh::dimension>>>
                                    value(nb_nodes, std::make_pair(0, vzero));

    size_t cell_i = 0;

    for (auto& cl : msh)
    {
        auto cell_ofs = revolution::priv::offset(msh, cl);
        Matrix<scalar_type, Dynamic, 1> x = sol.block(cell_ofs * cbs, 0, cbs, 1);

        const auto cell_nodes = post_mesh.nodes_cell(cell_i);
        auto  cbas = revolution::make_vector_monomial_basis(msh, cl, cell_degree);

       // Loop on the nodes of the cell
        for (auto& point_id : cell_nodes)
        {
            const auto pt = storage->points[point_id];

            const auto phi = cbas.eval_functions(pt);
            assert(phi.rows() == cbs);
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
    for (size_t i_node = 0; i_node < value.size(); i_node++)
    {
        const static_vector<scalar_type, Mesh::dimension> depl_avr =
                            value[i_node].second / double(value[i_node].first);

        const gmsh::Data tmp_data(i_node + 1, disk::convertToVectorGmsh(depl_avr));
        data.push_back(tmp_data);
    }

    // Create and init a nodedata view
    gmsh::NodeData nodedata(3, 0.0, "depl_node_cont", data, subdata);
    // Save the view
    nodedata.saveNodeData(filename, gmsh);

    return;
}

template<typename Mesh>
auto
run_stokes(const Mesh& msh, size_t degree, bool use_sym_grad = true)
{
    typedef Mesh mesh_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;
    typedef typename mesh_type::scalar_type scalar_type;

    typedef dynamic_matrix<scalar_type>     matrix_type;
    typedef disk::mechanics::BoundaryConditions<mesh_type> boundary_type;

    using point_type = typename mesh_type::point_type;

    auto rhs_fun = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
        return Matrix<scalar_type, 2, 1>::Zero();
    };

    auto wall = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
        return Matrix<scalar_type, 2, 1>::Zero();
    };
    auto movingWall = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
        return Matrix<scalar_type, 2, 1>{1,0};
    };

    typename revolution::hho_degree_info hdi(degree, degree);
    boundary_type bnd(msh);

    bnd.addDirichletBC(0, 1, movingWall);
    bnd.addDirichletBC(0, 2, wall);
    bnd.addDirichletBC(0, 3, wall);
    bnd.addDirichletBC(0, 4, wall);

    auto assembler = revolution::make_stokes_assembler(msh, hdi, bnd);

    scalar_type factor = (use_sym_grad)? 2. : 1.;

    for (auto cl : msh)
    {
        auto gr = revolution::make_hho_stokes(msh, cl, hdi, use_sym_grad);
        Matrix<scalar_type, Dynamic, Dynamic> stab;
        stab = make_hho_fancy_stabilization_vector(msh, cl, gr.first, hdi);
        auto dr = make_hho_divergence_reconstruction_stokes_rhs(msh, cl, hdi);
        auto cb = revolution::make_vector_monomial_basis(msh, cl, hdi.cell_degree());
        auto rhs = make_rhs(msh, cl, cb, rhs_fun);
        assembler.assemble(msh, cl, factor * (gr.second + stab), -dr, rhs);
    }

    assembler.finalize();

    //dump_sparse_matrix(assembler.LHS, "stokes.txt");

    size_t systsz = assembler.LHS.rows();
    size_t nnz = assembler.LHS.nonZeros();

    dynamic_vector<scalar_type> sol = dynamic_vector<scalar_type>::Zero(systsz);

    disk::solvers::pardiso_params<scalar_type> pparams;
    mkl_pardiso_ldlt(pparams, assembler.LHS, assembler.RHS, sol);

    std::ofstream ofs("velocity.dat");
    for (auto& cl : msh)
    {
        auto cbs   = revolution::vector_basis_size(degree, Mesh::dimension, Mesh::dimension);
        auto fbs   = revolution::vector_basis_size(degree, Mesh::dimension - 1, Mesh::dimension);
        auto cbs_B = revolution::scalar_basis_size(degree, Mesh::dimension);

        auto cell_ofs = revolution::priv::offset(msh, cl);
        auto num_other_faces = assembler.num_assembled_faces();
        auto offset_B = cbs * msh.cells_size() + fbs * num_other_faces + cbs_B * cell_ofs;

        Matrix<scalar_type, Dynamic, 1> s = sol.block(cell_ofs * cbs, 0, cbs, 1);
        Matrix<scalar_type, Dynamic, 1> press = sol.block( offset_B, 0, cbs_B, 1);
        auto bar = barycenter(msh, cl);
        ofs << bar.x() << " " << bar.y() << " " << s(0) << " " << s(1)<< " ";
        ofs << press(0) << std::endl;
        #if 0
        auto fcs = faces(msh, cl);
        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto fc = fcs[i];
            if (msh.is_boundary(fc))
            {
                auto face_offset = assembler.global_face_offset(msh, fc);
                Matrix<scalar_type, Dynamic, 1> s = sol.block(face_offset, 0, fbs, 1);

                auto bar = barycenter(msh, fc);
                ofs << bar.x() << " " << bar.y() << " " << s(0) << " " << s(1)<< " ";
                std::cout << bar.x() << " " << bar.y() << " " << s(0) << " " << s(1)<< " " << std::endl;
                ofs << press(0) << std::endl;
            }
        }
        #endif
    }

    ofs.close();

    compute_discontinuous_displacement( msh, sol, hdi, "depl2d.msh");

    return;

}

//#if 0
int main(int argc, char **argv)
{
    using RealType = double;
    bool use_sym_grad = false;

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

        run_stokes(msh, degree, use_sym_grad);
        std::cout << "fini" << std::endl;
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

        run_stokes(msh, degree, use_sym_grad);
    }

    return 0;
}
//#endif
