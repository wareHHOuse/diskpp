/*
 *       /\        Matteo Cicuttin (C) 2016-2021
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Guillaume Delay  (C) 2021         guillaume.delay@sorbonne-universite.fr
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

/*
 * Attention : penser a changer le triangle des ponts pour le symbole de Sorbonne ...
 * + mettre mes coordonnees
 */

/* this file deals with the heat equation with a discontinuous Galerkin method in time */
#include <iostream>
#include <regex>
#include <sstream>
#include <iomanip>

#include <unistd.h>

#include "loaders/loader.hpp"
#include "methods/hho"
#include "solvers/solver.hpp"
#include "output/silo.hpp"
#include "timecounter.h"
#include "colormanip.h"
#include "bases/bases.hpp"

using namespace disk;
using namespace Eigen;
using namespace std;


////////////////// time dG Assembler /////////////////
template<typename Mesh>
class heat_dG_assembler
{
    using T = typename Mesh::coordinate_type;

    std::vector<size_t>     compress_table;
    std::vector<size_t>     expand_table;
    hho_degree_info         di;
    size_t                  time_degree;
    std::vector<Triplet<T>> triplets, triplets_MAT_RHS;

    size_t num_all_faces, num_dirichlet_faces, num_other_faces, num_cells, space_system_size, system_size;

    class assembly_index
    {
        size_t idx;
        bool   assem;

      public:
        assembly_index(size_t i, bool as) : idx(i), assem(as) {}

        operator size_t() const
        {
            if (!assem)
                throw std::logic_error("Invalid assembly_index");

            return idx;
        }

        bool
        assemble() const
        {
            return assem;
        }

        friend std::ostream&
        operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };

  public:
    typedef dynamic_matrix<T> matrix_type;
    typedef dynamic_vector<T> vector_type;

    SparseMatrix<T> LHS, MAT_RHS;
    vector_type     RHS, RHS_F;

    heat_dG_assembler(const Mesh& msh, hho_degree_info hdi, size_t t_degree)
	: di(hdi), time_degree(t_degree)
    {
	auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool { return msh.is_boundary(fc); };

        num_all_faces       = msh.faces_size();
        num_dirichlet_faces = std::count_if(msh.faces_begin(), msh.faces_end(), is_dirichlet);
        num_other_faces     = num_all_faces - num_dirichlet_faces;
        num_cells = msh.cells_size();

        compress_table.resize(num_all_faces);
        expand_table.resize(num_other_faces);

        size_t compressed_offset = 0;
        for (size_t i = 0; i < num_all_faces; i++)
        {
            const auto fc = *std::next(msh.faces_begin(), i);
            if (!is_dirichlet(fc))
            {
                compress_table.at(i)               = compressed_offset;
                expand_table.at(compressed_offset) = i;
                compressed_offset++;
            }
        }

        const auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension - 1);
        const auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
        space_system_size    = fbs * num_other_faces + cbs * msh.cells_size();
        system_size = space_system_size * (time_degree + 1);

        LHS = SparseMatrix<T>(system_size, system_size);
        MAT_RHS = SparseMatrix<T>(system_size, system_size);
        RHS = vector_type::Zero(system_size);
        RHS_F = vector_type::Zero(system_size);
    }

    // here the Dirichlet data are not taken into account
    // the rhs function is not time-dependent
    template<typename Function, typename Function2>
    void
    assemble(const Mesh&                     msh,
             const typename Mesh::cell_type& cl,
             const matrix_type&              lhs,
             const vector_type&              rhs,
             const matrix_type&              mat_rhs,
             const Function&                 dirichlet_bf,
             const Function2&                init_fun)
    {
	auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool { return msh.is_boundary(fc); };

        const auto fbs    = scalar_basis_size(di.face_degree(), Mesh::dimension - 1);
        const auto cbs    = scalar_basis_size(di.cell_degree(), Mesh::dimension);
        const auto fcs    = faces(msh, cl);
        const auto fcs_id = faces_id(msh, cl);

        std::vector<assembly_index> asm_map;
        size_t loc_size = ( fcs.size() * fbs + cbs ) * (time_degree + 1);
        asm_map.reserve(loc_size);

        auto cell_offset = offset(msh, cl);
        for(size_t i = 0; i < cbs*(time_degree+1); i++)
            asm_map.push_back(assembly_index(cell_offset * cbs * (time_degree+1) + i, true));

        vector_type dirichlet_data = vector_type::Zero(fcs.size() * fbs);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            const auto fc              = fcs[face_i];
            const auto face_offset     = fcs_id[face_i]; // offset(msh, fc);
            const auto face_LHS_offset = num_cells * cbs * (time_degree+1)
                + compress_table.at(face_offset) * fbs * (time_degree+1);

            const bool dirichlet = is_dirichlet(fc);

            for (size_t i = 0; i < fbs*(time_degree+1); i++)
                asm_map.push_back(assembly_index(face_LHS_offset + i, !dirichlet));

            // Dirichlet data are not taken into account for the moment
            // if (dirichlet)
            // {
            //     dirichlet_data.block(face_i * fbs, 0, fbs, 1) =
            // 	    project_function(msh, fc, di.face_degree(), dirichlet_bf, di.face_degree());
            // }
        }

	
	// compute initial datum contribution to RHS
        vector_type u0 = vector_type::Zero(loc_size);
        u0.block(0,0,cbs,1) = project_function(msh, cl, di.cell_degree(), init_fun, di.cell_degree());
        auto rhs_modif = mat_rhs * u0 + rhs;

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if (asm_map[j].assemble())
                {
                    triplets.push_back(Triplet<T>(asm_map[i], asm_map[j], lhs(i, j)));
                    triplets_MAT_RHS.push_back(Triplet<T>(asm_map[i], asm_map[j], mat_rhs(i, j)));
		}
                // Dirichlet not taken into account
                // else
                //     RHS(asm_map[i]) -= lhs(i, j) * dirichlet_data(j);
            }

            RHS[ asm_map[i] ] += rhs_modif(i);
            RHS_F[ asm_map[i] ] += rhs(i);
        }
    } // assemble()

    template<typename Function>
    vector_type
    take_local_solution(const Mesh&                     msh,
                        const typename Mesh::cell_type& cl,
                        const vector_type&              solution,
                        const Function&                 dirichlet_bf)
    {
        const auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension - 1);
        const auto cbs = scalar_basis_size(di.cell_degree(), Mesh::dimension);
        const auto fcs = faces(msh, cl);
        const auto num_faces = fcs.size();

        vector_type ret = vector_type::Zero( (num_faces * fbs + cbs) * (time_degree + 1) );

        auto cell_offset = offset(msh, cl);

        // cell components
        ret.block(0, 0, (time_degree+1)*cbs, 1)
            = solution.block(cell_offset * cbs * (time_degree+1), 0, (time_degree+1)*cbs, 1);


        // face components
        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            const auto fc = fcs[face_i];

            auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool { return msh.is_boundary(fc); };

            const bool dirichlet = is_dirichlet(fc);

            // Dirichlet data not taken into account
            // if (dirichlet)
            // {
            // 	for(int l=0; l <= time_degree; l++)
            // 	    ret.block(cbs * (time_degree+1) + (face_i * time_degree + l) * fbs, 0, fbs, 1) =
            // 		project_function(msh, fc, di.face_degree(), dirichlet_bf, di.face_degree());
            // }

            if(!dirichlet)
            {
                const auto face_offset     = offset(msh, fc);
                const auto face_SOL_offset = num_cells * cbs * (time_degree+1)
                    + compress_table.at(face_offset) * fbs * (time_degree+1);

		
                ret.block(cbs * (time_degree+1) + face_i * (time_degree+1) * fbs, 0, (time_degree+1) * fbs, 1)
                    = solution.block(face_SOL_offset, 0, fbs * (time_degree+1), 1);
            }
        }

        return ret;
    }

    void
    finalize(void)
    {
        LHS.setFromTriplets(triplets.begin(), triplets.end());
        MAT_RHS.setFromTriplets(triplets_MAT_RHS.begin(), triplets_MAT_RHS.end());
        triplets.clear();
        triplets_MAT_RHS.clear();
    }

    size_t
    num_assembled_faces() const
    {
        return num_other_faces;
    }

};


template<typename Mesh>
auto
make_heat_dG_assembler(const Mesh& msh, const hho_degree_info& hdi, size_t time_degree)
{
    return heat_dG_assembler<Mesh>(msh, hdi, time_degree);
}

//////////////////////////////////////////////////////


template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
void
export_to_silo(const Mesh<T, 2, Storage>& msh,
               const Matrix<T, Dynamic, 1>& data, int cycle = -1)
{
    disk::silo_database silo;

    if (cycle == -1)
        silo.create("heat.silo");
    else
    {
        std::stringstream ss;
        ss << "out_" << cycle << ".silo";
        silo.create(ss.str());
    }

    silo.add_mesh(msh, "mesh");

    silo.add_variable("mesh", "sol", data, disk::zonal_variable_t );
    silo.close();
}

///////////////////////////////////////

template<typename Mesh>
bool
unsteady_laplacian_solver(const Mesh& msh, size_t degree, size_t time_steps, size_t time_degree)
{
    typedef typename Mesh::coordinate_type  scalar_type;
    typedef typename Mesh::point_type       point_type;
    using T = scalar_type;

    hho_degree_info hdi(degree, degree, degree+1);

    auto num_cells = msh.cells_size();
    auto num_faces = msh.faces_size();

    auto assembler = make_heat_dG_assembler(msh, hdi, time_degree);

    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);

    timecounter tc;

    int K = 1; // index for the exact solution
    scalar_type final_time = 2.;

#if false
    auto rhs_fun = [K](const point_type& pt) -> auto { return 2*M_PI*M_PI*K*K*std::sin(M_PI*K*pt.x()) * std::sin(M_PI*K*pt.y()); };
    auto bcs_fun = [](const point_type& pt) -> auto { return 0.0; };
    auto solution= [K](const double t, const point_type& pt) -> auto { return std::sin(M_PI*K*pt.x()) * std::sin(M_PI*K*pt.y()); };
    auto init_fun = [K](const point_type& pt) -> scalar_type {
			return std::sin(M_PI*K*pt.x()) * std::sin(M_PI*K*pt.y());
		    };
#else
    auto rhs_fun = [](const point_type& pt) -> auto { return 0.0; };
    auto bcs_fun = [](const point_type& pt) -> auto { return 0.0; };
    auto solution= [K](const double t, const point_type& pt) -> auto { return std::exp(-2*M_PI*M_PI*K*K*t) * std::sin(M_PI*K*pt.x()) * std::sin(M_PI*K*pt.y()); };
    auto init_fun = [K](const point_type& pt) -> scalar_type {
			return std::sin(M_PI*K*pt.x()) * std::sin(M_PI*K*pt.y());
		    };
#endif


    cout << "time_steps = " << time_steps << endl;
    disk::generic_mesh<T, 1>  time_mesh;
    disk::uniform_mesh_loader<T, 1> time_loader(0, final_time, time_steps);
    time_loader.populate_mesh(time_mesh);

    auto time_cell = *time_mesh.cells_begin();
    auto next_time_cell = *(++time_mesh.cells_begin());

    auto time_cb = make_scalar_monomial_basis(time_mesh, time_cell, time_degree);
    auto time_cb_next = make_scalar_monomial_basis(time_mesh, next_time_cell, time_degree);
    auto time_mass = make_mass_matrix(time_mesh, time_cell, time_cb);

    Matrix<T, Dynamic, Dynamic> time_deriv = Matrix<T, Dynamic, Dynamic>::Zero(time_degree+1, time_degree+1);
    auto qps_t = integrate(time_mesh, time_cell, 2*time_degree);
    for (auto& qp : qps_t)
    {
        auto phi   = time_cb.eval_functions( qp.point() );
        auto phi_t = time_cb.eval_gradients( qp.point() );
        time_deriv += qp.weight() * phi * phi_t.transpose();
    }
    
    //////
    Matrix<T, Dynamic, Dynamic> time_loc = Matrix<T, Dynamic, Dynamic>::Zero(time_degree+1, time_degree+1);
    auto t_fcs = faces(time_mesh, time_cell);
    auto qps_f_t = integrate(time_mesh, t_fcs[0], 2*time_degree);
    for (auto& qp : qps_f_t)
    {
        auto phi   = time_cb.eval_functions( qp.point() );
        time_loc += qp.weight() * phi * phi.transpose();
    }

    //////
    Matrix<T, Dynamic, Dynamic> time_loc_bis = Matrix<T, Dynamic, Dynamic>::Zero(time_degree+1, time_degree+1);
    auto qps_f_t_bis = integrate(time_mesh, t_fcs[1], 2*time_degree);
    for (auto& qp : qps_f_t_bis)
    {
        auto phi_prev   = time_cb.eval_functions( qp.point() );
        auto phi_next   = time_cb_next.eval_functions( qp.point() );
        time_loc_bis += qp.weight() * phi_next * phi_prev.transpose();
    }

    tc.tic();

    for (auto& cl : msh)
    {
        auto fcs    = faces(msh, cl);
        auto num_faces = fcs.size();
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
        auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);
        auto mass   = make_mass_matrix(msh, cl, cb);
        Matrix<scalar_type, Dynamic, Dynamic> ah = gr.second + stab;

        Matrix<scalar_type, Dynamic, Dynamic> lhs = Matrix<scalar_type, Dynamic, Dynamic>::Zero(ah.cols()*(time_degree+1), ah.cols()*(time_degree+1));

        // diffusion term
        // fill lhs with the tensor product between time and space
        // cell - cell
        for(size_t l1 = 0; l1 <= time_degree; l1++)
            for(size_t l2 = 0; l2 <= time_degree; l2++)
                lhs.block(l1*cbs, l2*cbs, cbs, cbs) = ah.block(0, 0, cbs, cbs) * time_mass(l1,l2);
	// cell - face
        for(size_t l1 = 0; l1 <= time_degree; l1++)
            for(size_t l2 = 0; l2 <= time_degree; l2++)
                for(size_t face_j = 0; face_j < num_faces; face_j++)
                    lhs.block(l1*cbs, cbs*(time_degree+1) + face_j*(time_degree+1)*fbs + l2*fbs,
                              cbs, fbs) = ah.block(0, cbs+face_j*fbs, cbs, fbs) * time_mass(l1,l2);
        // face - cell
        for(size_t l1 = 0; l1 <= time_degree; l1++)
            for(size_t l2 = 0; l2 <= time_degree; l2++)
                for(size_t face_i = 0; face_i < num_faces; face_i++)
                    lhs.block(cbs*(time_degree+1) + face_i*(time_degree+1)*fbs + l1*fbs, l2*cbs,
                              fbs, cbs) = ah.block(cbs + face_i * fbs, 0, fbs, cbs) * time_mass(l1,l2);
        // face - face
        for(size_t l1 = 0; l1 <= time_degree; l1++)
            for(size_t l2 = 0; l2 <= time_degree; l2++)
                for(size_t face_i = 0; face_i < num_faces; face_i++)
		    for(size_t face_j = 0; face_j < num_faces; face_j++)
                        lhs.block(cbs*(time_degree+1) + face_i*(time_degree+1)*fbs + l1*fbs,
                                  cbs*(time_degree+1) + face_j*(time_degree+1)*fbs + l2*fbs,
                                  fbs, fbs) = ah.block(cbs + face_i*fbs, cbs + face_j*fbs, fbs, fbs*num_faces) * time_mass(l1,l2);

        // derivative term (cell - cell only)
        for(size_t l1 = 0; l1 <= time_degree; l1++)
            for(size_t l2 = 0; l2 <= time_degree; l2++)
                lhs.block(l1*cbs, l2*cbs, cbs, cbs) += mass.block(0, 0, cbs, cbs) * time_deriv(l1,l2);

        // jump term (cell - cell only)
        for(size_t l1 = 0; l1 <= time_degree; l1++)
            for(size_t l2 = 0; l2 <= time_degree; l2++)
                lhs.block(l1*cbs, l2*cbs, cbs, cbs) += mass.block(0, 0, cbs, cbs) * time_loc(l1,l2);

        Matrix<scalar_type, Dynamic, Dynamic> mat_rhs = Matrix<scalar_type, Dynamic, Dynamic>::Zero(ah.cols()*(time_degree+1), ah.cols()*(time_degree+1));
        // jump term (cell - cell only)
        for(size_t l1 = 0; l1 <= time_degree; l1++)
            for(size_t l2 = 0; l2 <= time_degree; l2++)
                mat_rhs.block(l1*cbs, l2*cbs, cbs, cbs) += mass.block(0, 0, cbs, cbs) * time_loc_bis(l1,l2);

        // the rhs is computed for constant functions
        Matrix<scalar_type, Dynamic, 1> rhs = Matrix<scalar_type, Dynamic, 1>::Zero(lhs.cols());
        auto space_rhs = make_rhs(msh, cl, cb, rhs_fun, 1);
        for(size_t l1 = 0; l1 <= time_degree; l1++)
            rhs.block(l1*cbs, 0, cbs, 1) = space_rhs * time_mass(l1,0);


        assembler.assemble(msh, cl, lhs, rhs, mat_rhs, bcs_fun, init_fun);
    }
    cout << "end assembly loop" << endl;

    assembler.finalize();

    auto LHS = assembler.LHS;
    auto MAT_RHS = assembler.MAT_RHS;
    auto RHS = assembler.RHS;
    auto RHS_F = assembler.RHS_F;

    cout << "LHS.rows() = " << LHS.rows() << endl;

    tc.toc();
    std::cout << " Assembly time: " << tc << std::endl;    


    Matrix<scalar_type, Dynamic, 1> u;
    scalar_type t = 0.0;
    scalar_type dt = final_time/time_steps;
    scalar_type L2H1_error = 0.;

    size_t freq_exp = 300;

    // time loop
    for(size_t step_i = 0; step_i<time_steps; step_i++)
    {
        if(step_i % freq_exp == 0)
            std::cout << "Step " << step_i << std::endl;
        disk::solvers::pardiso_params<scalar_type> pparams;
        mkl_pardiso(pparams, LHS, RHS, u);

        Matrix<scalar_type, Dynamic, 1> sol_silo = Matrix<scalar_type, Dynamic, 1>::Zero(msh.cells_size());

        if(step_i % freq_exp == 0)
        {
            for (size_t i = 0; i < msh.cells_size(); i++)
                sol_silo(i) = u(i*cbs*(time_degree+1));
            export_to_silo( msh, sol_silo, step_i );
        }

        auto t_cell = *(time_mesh.cells_begin()+step_i);
        auto t_cb = make_scalar_monomial_basis(time_mesh, t_cell, time_degree);
        const auto qpst = integrate(time_mesh, t_cell , 2*time_degree);

        // update RHS with the previous solution
        RHS = RHS_F + MAT_RHS * u;
        size_t cell_i = 0;
        for (auto& cl : msh)
        {
            auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());

            /* compute the L2 projection in space and time */
            Matrix<scalar_type, Dynamic, 1> rhs_proj
                = Matrix<scalar_type, Dynamic, 1>::Zero( cbs * (time_degree+1) );
            const auto qps = integrate(msh, cl, 2*hdi.cell_degree());
            for(auto& qpt : qpst)
            {
                auto t_phi   = t_cb.eval_functions( qpt.point() );
                T time_point = qpt.point().x();

                for(auto& qp : qps)
                {
                    auto x_phi   = cb.eval_functions( qp.point() );
                    for(size_t l1 = 0; l1 <= time_degree; l1++)
                    {
                        rhs_proj.block(l1*cbs, 0, cbs, 1) += qp.weight() * qpt.weight() * t_phi[l1] * solution( time_point, qp.point() ) * x_phi;
                    }
                }
            }

            // use the mass matrix to compute the coordinates of the projection
            auto cell_mass   = make_mass_matrix(msh, cl, cb);
            Matrix<scalar_type, Dynamic, Dynamic> mass_matrix
                = Matrix<scalar_type, Dynamic, Dynamic>::Zero(cbs*(time_degree+1) , cbs*(time_degree+1));
            for(size_t l1 = 0; l1 <= time_degree; l1++)
            {
                for(size_t l2 = 0; l2 <= time_degree; l2++)
                {
                    mass_matrix.block(l1*cbs, l2*cbs, cbs, cbs) = time_mass(l1,l2) * cell_mass;
                }
            }
            Matrix<scalar_type, Dynamic, 1> proj
                = Matrix<scalar_type, Dynamic, 1>::Zero( cbs * (time_degree+1) );
            LLT< Matrix<scalar_type, Dynamic, Dynamic> > mat_llt;
            mat_llt.compute(mass_matrix);
            proj = mat_llt.solve(rhs_proj);


            /* compute L2-H1-error of the current time step */
            Matrix<scalar_type, Dynamic, 1> diff
                = u.block(cell_i*cbs*(time_degree+1), 0, cbs*(time_degree+1), 1) - proj;

            Matrix<scalar_type, Dynamic, Dynamic> grad_matrix
                = Matrix<scalar_type, Dynamic, Dynamic>::Zero(cbs , cbs);
            for(auto& qp : qps)
            {
                const auto g_phi = cb.eval_gradients( qp.point() );
                grad_matrix += qp.weight() * g_phi * g_phi.transpose();
            }


            for(size_t l1 = 0; l1 <= time_degree; l1++)
                for(size_t l2 = 0; l2 <= time_degree; l2++)
                    for(size_t I = 0; I < cbs; I++)
                        for(size_t J = 0; J < cbs; J++)
                        {
                            L2H1_error += diff(l1 * cbs + I) * diff(l2 * cbs + J) * time_mass(l1,l2) * grad_matrix(I,J);
                        }

            cell_i++;
        }
        t += dt;
    } // time loop

    cout << "final time = " << t << endl;
    std::cout << "L2-H1-error = " << std::sqrt(L2H1_error) << std::endl;

    return true;
}

/* run main with :
   ./heat_dG -m ../../../diskpp/meshes/2D_quads/diskpp/testmesh-16-16.quad -k 1 -N 8 -l 0
*/


int main(int argc, char **argv)
{
    using T = double;
    disk::cartesian_mesh<T, 2> msh;

    size_t      degree = 1;
    size_t      time_degree = 0;
    size_t      N = 8;
    char *      mesh_filename = nullptr;
    int ch;
    while ( (ch = getopt(argc, argv, "k:m:N:l:")) != -1 )
    {
	switch(ch)
        {
            case 'k':
                degree = std::stoi(optarg);
                break;

            case 'm':
                mesh_filename = optarg;
                break;

            case 'N':
                N = std::stoi(optarg);
                break;

            case 'l':
                time_degree = std::stoi(optarg);
                break;

            default:
                std::cout << "Invalid option" << std::endl;
                return 1;
	}
    }

    msh = load_cartesian_2d_mesh<T>(mesh_filename);

    std::cout << "Mesh loaded ..." << std::endl;
    unsteady_laplacian_solver(msh, degree, N, time_degree);
}

