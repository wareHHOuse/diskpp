/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2024
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */
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

#include "diskpp/methods/dga"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/output/silo.hpp"

#include "mumps.hpp"

#define MU_0 (4e-7*M_PI)
#define EPS_0 (8.8541878188e-12)

int main(int argc, char **argv)
{
    using CoordT = double;
    using SolT = std::complex<double>;

    /*
    if (argc != 2)
    {
        std::cout << "Please specify file name." << std::endl;
        return 1;
    }

    char *meshfile = argv[1];

    typedef disk::simplicial_mesh<CoordT, 3>  mesh_type;

    mesh_type msh;
    disk::netgen_mesh_loader<CoordT, 3> loader;

    if (!loader.read_mesh(meshfile))
    {
        std::cout << "Problem loading mesh." << std::endl;
        return 1;
    }
    loader.populate_mesh(msh);
    */
    
    using mesh_type = disk::simplicial_mesh<CoordT, 3>;

    mesh_type msh;
    auto mesher = make_simple_mesher(msh);
    for (int i = 0; i < 3; i++)
        mesher.refine();

    typedef Eigen::SparseMatrix<SolT>  sparse_matrix_type;
    typedef Eigen::Triplet<SolT>       triplet_type;

    auto num_all_edges = edges(msh).size();

    std::vector<bool> source_edges( num_all_edges );
    std::vector<bool> dirichlet_edges( num_all_edges );
    for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
    {
        auto fc = *itor;
        if (msh.boundary_id(fc) == 1 or msh.boundary_id(fc) == 4)
            continue;

        auto edgs = edges(msh, fc);
        dirichlet_edges.at( offset(msh, edgs[0]) ) = true;
        dirichlet_edges.at( offset(msh, edgs[1]) ) = true;
        dirichlet_edges.at( offset(msh, edgs[2]) ) = true;

        if (msh.boundary_id(fc) == 0) {
            source_edges.at( offset(msh, edgs[0]) ) = true;
            source_edges.at( offset(msh, edgs[1]) ) = true;
            source_edges.at( offset(msh, edgs[2]) ) = true;
        }
    }

    std::vector<size_t> compress_map, expand_map;
    compress_map.resize( num_all_edges );
    size_t system_size = std::count_if(dirichlet_edges.begin(), dirichlet_edges.end(), [](bool d) -> bool {return !d;});
    expand_map.resize( system_size );

    auto nnum = 0;
    for (size_t i = 0; i < num_all_edges; i++)
    {
        if ( dirichlet_edges.at(i) )
            continue;

        expand_map.at(nnum) = i;
        compress_map.at(i) = nnum++;
    }

    sparse_matrix_type          gA(system_size, system_size);
    disk::dynamic_vector<SolT>  gb(system_size), gx(system_size);
    gb = disk::dynamic_vector<SolT>::Zero(system_size);

    disk::vec3<SolT> src = {0.0, 1.0, 0.0};

    std::vector<triplet_type>       triplets;
    for (auto& cl : msh)
    {
        double omega = 2*M_PI*300e6;
        auto C = curl_matrix(msh, cl);
        auto Me = edge_matrix(msh, cl, 1.0*EPS_0);
        auto Mf = face_matrix(msh, cl, 1.0/MU_0);
        Eigen::Matrix<SolT,6,6> maxop = C.transpose() * Mf * C - omega*omega*Me;

        auto pev = primal_edge_vectors(msh, cl);

        auto eids = edge_ids(msh, cl);

        for (size_t i = 0; i < maxop.rows(); i++)
        {
            if ( dirichlet_edges.at(eids[i]) )
                continue;

            for (size_t j = 0; j < maxop.cols(); j++)
            {
                if ( dirichlet_edges.at(eids[j]) ) {
                    if (source_edges.at(eids[j])) {
                        gb(compress_map.at(eids[i])) = -maxop(i,j)*src.dot(pev[j]);
                    } else continue;
                }

                triplets.push_back( triplet_type(compress_map.at(eids[i]),
                                                 compress_map.at(eids[j]),
                                                 maxop(i,j)) );
            }

            //auto rhs_val = rhs_fun(pts[i]);
            //gb(compress_map.at(eids[i])) += rhs_val * vol * 0.25;
        }
    }

    gA.setFromTriplets(triplets.begin(), triplets.end());

    std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
    std::cout << "Dofs: " << gA.rows() << std::endl;

    std::cout << "Running MUMPS" << std::endl;
    gx = mumps_lu(gA, gb);
    std::cout << gx.norm() << std::endl;
    disk::dynamic_vector<SolT> sol = disk::dynamic_vector<SolT>::Zero( num_all_edges );

    for (size_t i = 0; i < gx.size(); i++)
        sol(expand_map.at(i)) = gx(i);

    disk::silo_database silo;
    silo.create("maxwell_dga.silo");
    silo.add_mesh(msh, "mesh");
    //silo.add_variable("mesh", "u", sol, disk::nodal_variable_t);

    for (auto& cl : msh)
    {
        auto eids = edge_ids(msh, cl);
        Eigen::Matrix<SolT, 6, 1> locsol;
        for (int i = 0; i < 6; i++)
            locsol(i) = sol( eids[i] );
    }

    return 0;
}
