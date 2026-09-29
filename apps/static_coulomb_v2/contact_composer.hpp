/*
 * contact_composer.hpp -- ContactTerms: a composable list of boundary terms.
 */

#pragma once

#include <vector>
#include <string>
#include <stdexcept>

#include "contact_terms.hpp"
#include "spatial_operator.hpp"

namespace hho_contact {

template<typename Mesh>
class ContactTerms
{
public:
    using T           = typename Mesh::coordinate_type;
    using matrix_type = disk::dynamic_matrix<T>;
    using vector_type = disk::dynamic_vector<T>;
    using cell_type   = typename Mesh::cell_type;
    using bc_type     = disk::vector_boundary_conditions<Mesh>;

private:
    enum class kind { contact_unilateral, contact_symmetry_condition, tresca, coulomb, symmetry };

    struct term
    {
        kind           knd;
        const bc_type* bnd;      // the view this term integrates over
        T gn0, gt0, theta;
        T coeff;                 // tresca: s, coulomb: F, otherwise unused
    };

    const spatial_operator<Mesh>& op;
    std::vector<term>             terms_;
    bool                          picard_ = false;

public:
    explicit ContactTerms(const spatial_operator<Mesh>& op_) : op(op_) {}

    // ---- composition ----------------------------------------------------
    void add_contact(const bc_type& bnd, T gn0, T gt0, T theta,
                     bool symmetry_condition = false)
    {
        terms_.push_back({symmetry_condition ? kind::contact_symmetry_condition
                                    : kind::contact_unilateral,
                          &bnd, gn0, gt0, theta, T(0)});
    }

    void add_tresca(const bc_type& bnd, T gt0, T theta, T s)
    {
        terms_.push_back({kind::tresca, &bnd, T(0), gt0, theta, s});
    }

    void add_coulomb(const bc_type& bnd, T gn0, T gt0, T theta, T F)
    {
        terms_.push_back({kind::coulomb, &bnd, gn0, gt0, theta, F});
        picard_ = true;
    }

    void add_symmetry(const bc_type& bnd, T gn0, T theta)
    {
        terms_.push_back({kind::symmetry, &bnd, gn0, T(0), theta, T(0)});
    }

    bool needs_picard() const { return picard_; }

    /* Nitsche consistency, to be SUBTRACTED from K. */
    matrix_type consistency(const cell_type& cl, size_t ci) const
    {
        const auto&  G = op.G[ci];
        const size_t n = op.total_dofs(cl);
        matrix_type  C = matrix_type::Zero(n, n);

        for (const auto& t : terms_)
        {
            switch (t.knd)
            {
            case kind::contact_unilateral:
            case kind::contact_symmetry_condition:
                C += disk::make_vector_hho_nitsche(op.msh, cl, op.di, G,
                        t.gn0, t.gt0, t.theta, op.mu, op.lam, *t.bnd,
                        /* include_normal = */ true,
                        /* include_tangential = */ false);
                break;
            case kind::tresca:
            case kind::coulomb:
                C += disk::make_vector_hho_nitsche(op.msh, cl, op.di, G,
                        t.gn0, t.gt0, t.theta, op.mu, op.lam, *t.bnd,
                        /* include_normal = */ false,
                        /* include_tangential = */ true);
                break;
            case kind::symmetry:
                C += disk::make_vector_hho_symmetry_nitsche(op.msh, cl, op.di, G,
                        t.gn0, t.theta, op.mu, op.lam, *t.bnd);
                break;
            default: break;
            }
        }
        return C;
    }

    /* Constant (u-independent) penalty jacobians, to be ADDED to K.
     * Their residual contribution J.u is then produced by A.u - L.        */
    matrix_type linear_jacobian(const cell_type& cl, size_t ci) const
    {
        const auto&  G = op.G[ci];
        const size_t n = op.total_dofs(cl);
        matrix_type  J = matrix_type::Zero(n, n);

        for (const auto& t : terms_)
        {
            switch (t.knd)
            {
            case kind::contact_symmetry_condition:
            {
                const vector_type u0 = vector_type::Zero(G.cols());
                J += disk::make_vector_hho_contact_jacobian(op.msh, cl, op.di, G,
                        t.gn0, t.theta, op.mu, op.lam, *t.bnd, u0,
                        /* symmetry_condition = */ true,
                        /* symmetry is Dirichlet-type */
                        hho_contact::trace_variant::face);
                break;
            }
            case kind::symmetry:
                J += disk::make_vector_hho_symmetry_jacobian(op.msh, cl, op.di, G,
                        t.gn0, t.theta, op.mu, op.lam, *t.bnd,
                        hho_contact::trace_variant::face);
                break;
            default: break;
            }
        }
        return J;
    }

    // ---- per-Newton-iteration pieces ------------------------------------

    /* Accumulate u-dependent (J, R) in place.  p is the frozen Picard
     * iterate                                              */
    void accumulate_nonlinear(const cell_type& cl, size_t ci,
                              const vector_type& u, const vector_type& p,
                              matrix_type& J, vector_type& R) const
    {
        const auto& G = op.G[ci];

        for (const auto& t : terms_)
        {
            switch (t.knd)
            {
            case kind::contact_unilateral:
                J += disk::make_vector_hho_contact_jacobian(op.msh, cl, op.di, G,
                        t.gn0, t.theta, op.mu, op.lam, *t.bnd, u,
                        /* symmetry_condition = */ false, op.variant);
                R += disk::make_vector_hho_contact_rhs(op.msh, cl, op.di, G,
                        t.gn0, t.theta, op.mu, op.lam, *t.bnd, u,
                        /* symmetry_condition = */ false, op.variant);
                break;
            case kind::tresca:
                J += disk::make_vector_hho_tresca_jacobian(op.msh, cl, op.di, G,
                        t.gt0, t.theta, op.mu, op.lam, t.coeff, *t.bnd, u,
                        op.variant);
                R += disk::make_vector_hho_tresca_rhs(op.msh, cl, op.di, G,
                        t.gt0, t.theta, op.mu, op.lam, t.coeff, *t.bnd, u,
                        op.variant);
                break;
            case kind::coulomb:
                J += disk::make_vector_hho_coulomb_jacobian(op.msh, cl, op.di, G,
                        t.gn0, t.gt0, t.theta, op.mu, op.lam, t.coeff, *t.bnd, u, p,
                        op.variant);
                R += disk::make_vector_hho_coulomb_rhs(op.msh, cl, op.di, G,
                        t.gn0, t.gt0, t.theta, op.mu, op.lam, t.coeff, *t.bnd, u, p,
                        op.variant);
                break;
            default: break;      // linear terms already live in A
            }
        }
    }

    /* True if any term still needs per-iteration work: lets the solver
     * skip the switch entirely for a purely linear problem.               */
    bool has_nonlinear() const
    {
        for (const auto& t : terms_)
            if (t.knd == kind::contact_unilateral ||
                t.knd == kind::tresca ||
                t.knd == kind::coulomb)
                return true;
        return false;
    }
};

} // namespace hho_contact