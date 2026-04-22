/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2026
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#pragma once

#include <iostream>

namespace disk::solvers {

enum class feast_inner_solver {
    eigen_sparselu,
    mumps
};

enum class feast_status {
    success,
    did_not_converge,
    invalid_input,
    inner_solver_problem,
    subspace_too_small
};

inline std::ostream&
operator<<(std::ostream& os, const feast_status& fs)
{
    switch (fs) {
        case feast_status::success:
            os << "success";
            break;
        
        case feast_status::did_not_converge:
            os << "did_not_converge";
            break;

        case feast_status::invalid_input:
            os << "invalid_input";
            break;

        case feast_status::inner_solver_problem:
            os << "inner_solver_problem";
            break;

        case feast_status::subspace_too_small:
            os << "subspace_too_small";
            break;
    }

    return os;
}

template<typename T>
struct feast_eigensolver_params
{
    bool    verbose;
    int     tolerance;
    T       min_eigval;
    T       max_eigval;
    int     subspace_size;
    int     eigvals_found;
    int     feast_info;
    size_t  max_iter;
    feast_inner_solver fis = feast_inner_solver::eigen_sparselu;
};

static double quadrature_xs[] = {
    0.183434642495649, -0.183434642495649,
    0.525532409916328, -0.525532409916328,
    0.796666477413626, -0.796666477413626,
    0.960289856497536, -0.960289856497536
};

static double quadrature_ws[] = {
    0.362683783378361, 0.362683783378361,
    0.313706645877887, 0.313706645877887,
    0.222381034453374, 0.222381034453374,
    0.101228536290376, 0.101228536290376
};

} // namespace disk::solvers