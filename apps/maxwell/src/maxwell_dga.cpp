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
#include "diskpp/solvers/solver.hpp"

#include "mumps.hpp"

#define MU_0 (4e-7*M_PI)
#define EPS_0 (8.8541878188e-12)

int main(int argc, char **argv)
{
    using CoordT = double;
    using SolT = std::complex<double>;

    using mesh_type = disk::simplicial_mesh<CoordT, 3>;

    mesh_type msh;
    auto mesher = make_simple_mesher(msh);
    for (int i = 0; i < 4; i++)
        mesher.refine();
    
    typedef Eigen::SparseMatrix<SolT>  sparse_matrix_type;
    typedef Eigen::Triplet<SolT>       triplet_type;

    auto num_all_edges = edges(msh).size();

    std::vector<bool> source_edges( num_all_edges );
    std::vector<bool> dirichlet_edges( num_all_edges );
    for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
    {
        auto fc = *itor;
        if (msh.boundary_id(fc) == 2 or msh.boundary_id(fc) == 5)
            continue;

        auto edgids = edge_ids(msh, fc);
        assert(edgids.size() == 3);
        dirichlet_edges.at( edgids[0] ) = true;
        dirichlet_edges.at( edgids[1] ) = true;
        dirichlet_edges.at( edgids[2] ) = true;

        if (msh.boundary_id(fc) == 0) {
            source_edges.at( edgids[0] ) = true;
            source_edges.at( edgids[1] ) = true;
            source_edges.at( edgids[2] ) = true;
        }
    }

    std::vector<size_t> compress_map, expand_map;
    compress_map.resize( num_all_edges );
    size_t system_size = std::count_if(dirichlet_edges.begin(), dirichlet_edges.end(), [](bool d) -> bool {return !d;});
    expand_map.resize( system_size );

    for (size_t i = 0, nnum = 0; i < num_all_edges; i++)
    {
        if ( dirichlet_edges.at(i) )
            continue;

        expand_map.at(nnum) = i;
        compress_map.at(i) = nnum++;
    }

    sparse_matrix_type          gA(system_size, system_size);
    disk::dynamic_vector<SolT>  gb(system_size), gx(system_size);
    gb = disk::dynamic_vector<SolT>::Zero(system_size);

    double omega = 2*M_PI*300e6;
    disk::vec3<SolT> src = {0.0, 1.0, 0.0};

    std::vector<triplet_type>       triplets;
    for (auto& cl : msh)
    {
        auto C = curl_matrix(msh, cl);
        auto Me = edge_matrix(msh, cl, 1.0*EPS_0);
        auto Mf = face_matrix(msh, cl, 1.0/MU_0);
        Eigen::Matrix<SolT,6,6> maxop = C.transpose() * Mf * C - omega*omega*Me;

        auto pev = primal_edge_vectors(msh, cl);
        auto dav = dual_area_vectors(msh, cl);

        auto edgs = edges(msh, cl);
        auto eids = edge_ids(msh, cl);

        for (size_t i = 0; i < maxop.rows(); i++)
        {
            if ( dirichlet_edges.at(eids[i]) )
                continue;

            for (size_t j = 0; j < maxop.cols(); j++)
            {
                if ( dirichlet_edges.at(eids[j]) )
                {
                    if (source_edges.at(eids[j]))
                    {
                        auto gi = compress_map.at(eids[i]);
                        gb(gi) += -maxop(i,j)*src.dot(pev[j]);
                    }
                    
                    continue;
                }

                triplets.push_back( triplet_type(compress_map.at(eids[i]),
                                                 compress_map.at(eids[j]),
                                                 maxop(i,j)) );
            }

            //auto bar = barycenter(msh, cl);
            //gb(compress_map.at(eids[i])) += (M_PI*M_PI/MU_0 - omega*omega*EPS_0)*src.dot(dav[i])*std::sin(M_PI*bar.x())*std::sin(M_PI*bar.z());
        }
    }

    gA.setFromTriplets(triplets.begin(), triplets.end());

    std::cout << "Mesh elements: " << msh.cells_size() << std::endl;
    std::cout << "Dofs: " << gA.rows() << std::endl;

    std::cout << "Running MUMPS" << std::endl;
    gx = mumps_lu(gA, gb);
    
    //Eigen::SparseLU<Eigen::SparseMatrix<SolT>> solver;
    //                solver.compute(gA);
    //                if(solver.info() != Eigen::Success) {
    //                    std::cout << "SparseLU failed" << std::endl;
    //                }
    //                gx = solver.solve(gb);
    
    disk::dynamic_vector<SolT> sol = disk::dynamic_vector<SolT>::Zero( num_all_edges );

    for (size_t i = 0; i < gx.size(); i++)
        sol(expand_map.at(i)) = gx(i);

    std::vector<CoordT> ex, ey, ez, emag;

    for (auto& cl : msh)
    {
        auto eids = edge_ids(msh, cl);
        Eigen::Matrix<SolT, 6, 1> locsol = Eigen::Matrix<SolT, 6, 1>::Zero();
        for (int i = 0; i < 6; i++)
            locsol(i) = sol( eids[i] );

        //auto pev = primal_edge_vectors(msh, cl);
        //for (int i = 0; i < 6; i++)
        //    locsol(i) = src.dot(pev[i]);

        auto pav = primal_area_vectors(msh, cl);
        auto vf = 12*volume_signed(msh, cl);

        disk::vec3<SolT> locfield;
        locfield =  (locsol(0)*pav[1] - locsol(1)*pav[2] + locsol(2)*pav[3])/vf;
        locfield += (locsol(0)*pav[0] - locsol(3)*pav[2] + locsol(4)*pav[3])/vf;
        locfield += (locsol(1)*pav[0] - locsol(3)*pav[1] + locsol(5)*pav[3])/vf;
        locfield += (locsol(2)*pav[0] - locsol(4)*pav[1] + locsol(5)*pav[2])/vf;

        ex.push_back( real(locfield(0)) );
        ey.push_back( real(locfield(1)) );
        ez.push_back( real(locfield(2)) );
        emag.push_back( sqrt(real(locfield.dot(locfield))) );
    }

    SolT ee = 0.0;
    SolT me = 0.0;
    for (auto& cl : msh)
    {
        auto eids = edge_ids(msh, cl);
        auto Me = edge_matrix(msh, cl, 1.0*EPS_0);
        Eigen::Matrix<SolT, 6, 1> Us;
        for (int i = 0; i < 6; i++)
            Us(i) = sol( eids[i] );
        ee += Us.dot(Me*Us);

        auto Mf = face_matrix(msh, cl, 1.0/MU_0);
        auto C = curl_matrix(msh, cl);
        Eigen::Matrix<SolT, 4, 1> Phis = C*Us/omega;
        me += Phis.dot(Mf*Phis);
    }

    std::cout << "Electric energy: " << ee << std::endl;
    std::cout << "Magnetic energy: " << me << std::endl;


    disk::silo_database silo;
    silo.create("maxwell_dga.silo");
    silo.add_mesh(msh, "mesh");
    silo.add_variable("mesh", "ex", ex, disk::zonal_variable_t);
    silo.add_variable("mesh", "ey", ey, disk::zonal_variable_t);
    silo.add_variable("mesh", "ez", ez, disk::zonal_variable_t);
    silo.add_variable("mesh", "emag", emag, disk::zonal_variable_t);

    return 0;
}
