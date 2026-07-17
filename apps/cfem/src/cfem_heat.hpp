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

namespace disk::cfem::heat {

template<typename T>
struct solver_state
{
    using mesh_type = disk::simplicial_mesh<T,2>;
    using spmat = Eigen::SparseMatrix<T>;
    using dvec = disk::dynamic_vector<T>;

    mesh_type                               msh;
    std::vector<std::optional<T>>           dirichlet_values;

    Eigen::SparseLU<spmat> lu_M;

    spmat           K, M;
    dvec            f;
    T               dt;
};

template<typename T>
void init(solver_state<T>& state)
{
    std::vector<std::pair<int, T>>  dirichlet_vals {
        //{13, 1.0},
        //{14, 1.0},
        {15, 1.0},
        //{11, 0.0},
        //{10, 0.0}
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

    state.K = typename solver_state<T>::spmat(system_size, system_size);
    state.M = typename solver_state<T>::spmat(system_size, system_size);
    state.f = dynamic_vector<T>::Zero(system_size);
}

template<typename T>
void assemble(solver_state<T>& state)
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
            kappa(1,1) = 0.1;
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
dynamic_vector<T>
timestep(solver_state<T>& state, dynamic_vector<T>& u_curr)
{
    dynamic_vector<T> u_next;

    u_next = u_curr - state.dt * state.lu_M.solve(state.K * u_curr);

    for (size_t i = 0; i < state.dirichlet_values.size(); i++) {
        if (state.dirichlet_values[i]) {
            u_next[i] = *state.dirichlet_values[i];
        }
    }

    return u_next;
}

} // namespace diskpp::cfem::heat