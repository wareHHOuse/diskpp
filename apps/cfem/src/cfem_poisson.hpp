/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2026
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
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

#pragma once

#include "diskpp/quadratures/quadratures.hpp"

namespace disk::cfem::poisson {

template<typename T>
struct solver_state
{
    using mesh_type = disk::simplicial_mesh<T,2>;
    using spmat = Eigen::SparseMatrix<T>;
    using dvec = disk::dynamic_vector<T>;

    mesh_type                               msh;
    std::vector<std::optional<T>>           dirichlet_values;
    std::vector<std::optional<size_t>>      compress_map;
    std::vector<size_t>                     expand_map;

    spmat           K;
    dvec            u, f;
};



template<typename T>
void init(solver_state<T>& state)
{
    std::vector<std::pair<int, T>>  dirichlet_vals {
        {1, 0.0},
        {2, 0.0},
        {3, 0.0},
        {4, 0.0},
    };

    state.dirichlet_values.resize(
        state.msh.points_size()
    );

    for (auto& fc : faces(state.msh)) {
        auto bi = state.msh.boundary_info(fc);
        if ( bi.is_boundary() ) {
            auto ptids = fc.point_ids();
            for (const auto& [tag, value] : dirichlet_vals) {
                if (tag == bi.tag()) {
                    state.dirichlet_values[ptids[0]] = value;
                    state.dirichlet_values[ptids[1]] = value;
                }
            }
        }
    }

    size_t system_size = 0;
    for (const auto& dv : state.dirichlet_values) {
        if (not dv) {
            system_size++;
        }
    }

    state.compress_map.resize( state.msh.points_size() );
    state.expand_map.resize( system_size );

    size_t nnum = 0;
    for (size_t i = 0; i < state.msh.points_size(); i++) {
        if ( state.dirichlet_values[i] ) {
            continue;
        }

        state.expand_map[nnum] = i;
        state.compress_map[i] = nnum++;
    }

    state.K = typename solver_state<T>::spmat(system_size, system_size);
    state.u = disk::dynamic_vector<T>::Zero(system_size);
    state.f = disk::dynamic_vector<T>::Zero(system_size);
}



template<typename T>
void assemble(solver_state<T>& state)
{
    using triplet_type = Eigen::Triplet<T>;
    std::vector<triplet_type>       triplets;

    auto f = [](const typename disk::simplicial_mesh<T, 2>::point_type& pt) -> auto {
        return 2.0 * M_PI * M_PI *std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
    };

    for (const auto& cl : state.msh)
    {
        disk::static_matrix<T, 2, 2> kappa = disk::static_matrix<T, 2, 2>::Zero();

        kappa(0,0) = 1.0;
        kappa(1,1) = 1.0;

        auto loc_K = disk::cfem::stiffness_matrix(state.msh, cl, kappa);
        auto loc_f = disk::cfem::make_rhs(state.msh, cl, f);

        auto ptids = cl.point_ids();

        for (size_t i = 0; i < loc_K.rows(); i++)
        {
            if ( state.dirichlet_values[ptids[i]] ) {
                continue;
            }

            auto ci = state.compress_map[ptids[i]].value();

            for (size_t j = 0; j < loc_K.cols(); j++)
            {
                if ( state.dirichlet_values[ptids[j]] ) {
                    auto value = *state.dirichlet_values[ptids[j]];
                    state.f(ci) -= value * loc_K(i,j);
                    continue;
                }
                auto cj = state.compress_map[ptids[j]].value();
                triplets.push_back( triplet_type(ci, cj, loc_K(i,j)) );
            }

            state.f(ci) += loc_f(i);
        }
    }

    state.K.setFromTriplets(triplets.begin(), triplets.end());
}



template<typename T>
void solve(solver_state<T>& state)
{
    Eigen::SparseLU<typename solver_state<T>::spmat> solver(state.K);
    state.u = solver.solve(state.f);
}

template<typename T>
disk::dynamic_vector<T>
expand_solution(const solver_state<T>& state)
{
    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(
        state.msh.points_size()
    );

    for (size_t i = 0; i < state.u.size(); i++) {
        sol[ state.expand_map[i] ] = state.u(i);
    }

    for (size_t i = 0; i < state.dirichlet_values.size(); i++) {
        if ( state.dirichlet_values[i] ) {
            sol[i] += *state.dirichlet_values[i];
        }
    }

    return sol;
}

template<typename T>
void picard(solver_state<T>& state)
{
    using triplet_type = Eigen::Triplet<T>;

    /* Manufactured solution */
    auto f = [](const typename disk::simplicial_mesh<T, 2>::point_type& pt) -> auto {
        auto cpix = std::cos(M_PI*pt.x());
        auto cpiy = std::cos(M_PI*pt.y());
        auto spix = std::sin(M_PI*pt.x());
        auto spiy = std::sin(M_PI*pt.y());
        auto u = spix * spiy;
        auto modgradu2 = M_PI*M_PI*( cpix*cpix*spiy*spiy + spix*spix*cpiy*cpiy );
        return 2*M_PI*M_PI*u*(1+u*u) - 2*u*modgradu2;
    };

    /* Iterate */
    for (size_t pi = 0; pi < 20; pi++) {
        /* Clear state */
        state.K.setZero();
        state.f.setZero();
        std::vector<triplet_type> triplets;

        /* Full solution including dirichlet dofs */
        auto u_full = expand_solution(state);

        for (const auto& cl : state.msh)
        {
            static_matrix<T, 3, 3> loc_K = static_matrix<T, 3, 3>::Zero();
            /* evaluate ( kappa(u[n-1])grad(u[n]), grad(v) )_T = (f,v)_T */
            auto ptids = cl.point_ids();
            auto qps = disk::integrate(state.msh, cl, 2);
            for (const auto& qp : qps) {
                auto phi = disk::cfem::eval_basis(state.msh, cl, qp.point());
                auto dphi = disk::cfem::eval_basis_grad(state.msh, cl);
                /* evaluate u in the quadrature point */
                T u_pt = 0.0;
                for (size_t i = 0; i < phi.size(); i++) {
                    u_pt += u_full[ptids[i]] * phi[i];
                }
                /* integrate */
                loc_K += qp.weight() * (1 + u_pt*u_pt)*(dphi * dphi.transpose());
            }
            auto loc_f = disk::cfem::make_rhs(state.msh, cl, f);

            /* Assemble */
            for (size_t i = 0; i < loc_K.rows(); i++)
            {
                if ( state.dirichlet_values[ptids[i]] ) {
                    continue;
                }

                auto ci = state.compress_map[ptids[i]].value();

                for (size_t j = 0; j < loc_K.cols(); j++)
                {
                    if ( state.dirichlet_values[ptids[j]] ) {
                        auto value = *state.dirichlet_values[ptids[j]];
                        state.f(ci) -= value * loc_K(i,j);
                        continue;
                    }
                    auto cj = state.compress_map[ptids[j]].value();
                    triplets.push_back( triplet_type(ci, cj, loc_K(i,j)) );
                }

                state.f(ci) += loc_f(i);
            }
        }

        /* build matrix */
        state.K.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();

        /* solve */
        Eigen::SparseLU<typename solver_state<T>::spmat> solver(state.K);
        disk::dynamic_vector<T> u_picard = solver.solve(state.f);
        
        T omega = 1;

        disk::dynamic_vector<T> u_next = (1.0 - omega) * state.u + omega * u_picard;

        /* check if converged */
        auto diffnorm = (u_next - state.u).norm();
        std::cout << "Picard iteration " << pi << ", norm: " << diffnorm << "\n";
        state.u = u_next;
        if (diffnorm < 1e-6) {
            break;
        }
    } // for
}

template<typename T>
void postpro(solver_state<T>& state)
{
    auto sol = expand_solution(state);
    
    disk::silo_database silo_db;
    silo_db.create("test2.silo");
    silo_db.add_mesh(state.msh, "mesh");
    silo_db.add_variable("mesh", "u", sol, disk::nodal_variable_t);
}

} // namespace disk::cfem::poisson