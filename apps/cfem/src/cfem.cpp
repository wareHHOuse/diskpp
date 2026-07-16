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

#include <iostream>
#include <regex>
#include <optional>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/loaders/loader_gmsh.hpp"

#include "diskpp/cfem/cfem.hpp"
#include "diskpp/output/silo.hpp"
#include "diskpp/mesh/meshgen.hpp"


template<typename T>
struct cfem_solver_state
{
    using mesh_type = disk::simplicial_mesh<T,2>;
    using spmat = Eigen::SparseMatrix<T>;
    using dvec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    mesh_type                               msh;
    std::vector<std::optional<T>>           dirichlet_values;
    std::vector<std::optional<size_t>>      compress_map;
    std::vector<size_t>                     expand_map;

    spmat           K;
    dvec            u, f;
};

template<typename T>
struct cfem_solver_state_transient
{
    using mesh_type = disk::simplicial_mesh<T,2>;
    using spmat = Eigen::SparseMatrix<T>;
    using dvec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    mesh_type                               msh;
    std::vector<std::optional<T>>           dirichlet_values;

    Eigen::SparseLU<spmat> lu_M;

    spmat           K, M;
    dvec            u, f;
    T               dt;
};

template<typename T>
void init_maps_steady(cfem_solver_state<T>& state)
{
    std::vector<std::pair<int, T>>  dirichlet_vals {
        {11, 1.0},
        {10, 0.0},
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

    state.K = typename cfem_solver_state<T>::spmat(system_size, system_size);
    state.u = cfem_solver_state<T>::dvec::Zero(system_size);
    state.f = cfem_solver_state<T>::dvec::Zero(system_size);
}

template<typename T>
void init(cfem_solver_state_transient<T>& state)
{
    std::vector<std::pair<int, T>>  dirichlet_vals {
        {13, 1.0},
        {14, 1.0},
        {15, 1.0},
        {11, 0.0},
        {10, 0.0}
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

    size_t system_size = state.msh.points_size();

    state.K = typename cfem_solver_state<T>::spmat(system_size, system_size);
    state.M = typename cfem_solver_state<T>::spmat(system_size, system_size);
    state.u = cfem_solver_state<T>::dvec::Zero(system_size);
    state.f = cfem_solver_state<T>::dvec::Zero(system_size);
}

template<typename T>
void assemble_steady(cfem_solver_state<T>& state)
{
    using triplet_type = Eigen::Triplet<T>;
    std::vector<triplet_type>       triplets;

    auto f = [](const typename disk::simplicial_mesh<T, 2>::point_type& pt) -> auto {
        return 0.0;
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
void assemble(cfem_solver_state_transient<T>& state)
{
    using triplet_type = Eigen::Triplet<T>;
    std::vector<triplet_type>       triplets_K;
    std::vector<triplet_type>       triplets_M;

    auto f = [](const typename disk::simplicial_mesh<T, 2>::point_type& pt) -> auto {
        return 0.0;
    };

    for (const auto& cl : state.msh)
    {
        disk::static_matrix<T, 2, 2> kappa = disk::static_matrix<T, 2, 2>::Zero();

        auto bar = barycenter(state.msh, cl);

        if (bar.x() < 0.9) {
            kappa(0,0) = 0.1;
            kappa(1,1) = 0.1;
        } else {
            kappa(0,0) = 0.1;
            kappa(1,1) = 0.001;
        }

        auto loc_K = disk::cfem::stiffness_matrix(state.msh, cl, kappa);
        auto loc_M = disk::cfem::mass_matrix(state.msh, cl);
        auto loc_f = disk::cfem::make_rhs(state.msh, cl, f);

        auto ptids = cl.point_ids();

        for (size_t i = 0; i < loc_K.rows(); i++)
        {
            auto ci = ptids[i];
            for (size_t j = 0; j < loc_K.cols(); j++)
            {
                auto cj = ptids[j];
                triplets_K.push_back( triplet_type(ci, cj, loc_K(i,j)) );
                triplets_M.push_back( triplet_type(ci, cj, loc_M(i,j)) );
            }

            state.f(ci) += loc_f(i);
        }
    }

    state.K.setFromTriplets(triplets_K.begin(), triplets_K.end());
    state.M.setFromTriplets(triplets_M.begin(), triplets_M.end());
    state.lu_M.compute(state.M);
}

template<typename T>
void solve_steady(cfem_solver_state<T>& state)
{
    Eigen::SparseLU<typename cfem_solver_state<T>::spmat> solver(state.K);
    state.u = solver.solve(state.f);
}

template<typename T>
typename cfem_solver_state_transient<T>::dvec
timestep(cfem_solver_state_transient<T>& state,
    typename cfem_solver_state_transient<T>::dvec& u_curr)
{
    typename cfem_solver_state_transient<T>::dvec u_next;

    u_next = u_curr - state.dt * state.lu_M.solve(state.K * u_curr);

    for (size_t i = 0; i < state.dirichlet_values.size(); i++) {
        if (state.dirichlet_values[i]) {
            u_next[i] = *state.dirichlet_values[i];
        }
    }

    return u_next;
}

template<typename T>
void postpro_steady(cfem_solver_state<T>& state)
{
    typename cfem_solver_state<T>::dvec sol = cfem_solver_state<T>::dvec::Zero(
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

    disk::silo_database silo_db;
    silo_db.create("test2.silo");
    silo_db.add_mesh(state.msh, "mesh");
    silo_db.add_variable("mesh", "u", sol, disk::nodal_variable_t);
}




int main(int argc, char **argv)
{
    using RealType = double;

    if (argc < 2) {
        return 1;
    }

    using mesh_type = disk::simplicial_mesh<RealType, 2>;

    disk::gmsh_geometry_loader<mesh_type> loader;
    loader.read_mesh(argv[1]);


    cfem_solver_state_transient<RealType> state;
    state.dt = 1e-3;

    loader.populate_mesh(state.msh);

    std::cout << "init & asm\n"; 

    init(state);
    assemble(state);

    using dvec = typename cfem_solver_state_transient<RealType>::dvec;
    dvec u_curr = dvec::Zero(state.msh.points_size());
    dvec u_next = dvec::Zero(state.msh.points_size());

    for (size_t i = 0; i < state.dirichlet_values.size(); i++) {
        if (state.dirichlet_values[i]) {
            u_curr[i] = *state.dirichlet_values[i];
        }
    }


    for (size_t ts = 0; ts < 10000; ts++) {
        if (ts % 10 == 0) {
            std::cout << "ts " << ts << "\n";
            disk::silo_database silo_db;
            std::string fname = "ts_" + std::to_string(ts) + ".silo";
            silo_db.create(fname);
            silo_db.add_mesh(state.msh, "mesh");
            silo_db.add_variable("mesh", "u", u_curr, disk::nodal_variable_t);
        }

        u_next = timestep(state, u_curr);

        u_curr = u_next;
    }


    return 0;
}
