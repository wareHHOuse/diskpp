/*
 * contact_composer.hpp -- ContactTerms: a composable list of boundary terms.
 *
 * Every term the drivers used lives here, wrapped from contact_terms.hpp,
 * and is classified by WHERE its work belongs:
 *
 *   term        consistency (into A, once)   linear jac (into A, once)   nonlinear (per Newton it.)
 *   ---------   --------------------------   -------------------------   --------------------------
 *   contact,    full n+t Nitsche             --                          [.]_- jacobian + rhs
 *   unilateral
 *   contact,    full n+t Nitsche             bilateral normal penalty    --
 *   bilateral                                (constant: no [.]_-)
 *   tresca      --                           --                          ball projection jac + rhs
 *   coulomb     --                           --                          jac + rhs, needs frozen p
 *   symmetry    normal-only Nitsche          bilateral normal penalty    --
 *
 * The classification is the efficiency point: bilateral contact and
 * symmetry are LINEAR (their jacobian is u-independent -- see
 * make_vector_hho_symmetry_jacobian, which calls the contact jacobian with
 * u = 0), and their residual is exactly J.u.  Folding J into the cached
 * system matrix A makes R = A.u - L include them for free, and removes
 * their re-assembly from every Newton iteration.  In the Bostan-Han driver
 * this leaves the Tresca ball projection as the ONLY per-iteration term.
 *
 * A Coulomb term flips needs_picard(): NewtonSolver reads it to decide
 * whether to wrap the Newton loop in the outer fixed-point iteration.
 * Symmetry is NOT tied to Coulomb: any experiment adds it -- or not --
 * with its own boundary view.
 */

#pragma once

#include <vector>

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
    enum class kind { contact_unilateral, contact_bilateral, tresca, coulomb, symmetry };

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
                     bool bilateral = false)
    {
        terms_.push_back({bilateral ? kind::contact_bilateral
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
            case kind::contact_bilateral:
                C += disk::make_vector_hho_nitsche(op.msh, cl, op.di, G,
                        t.gn0, t.gt0, t.theta, op.mu, op.lam, *t.bnd);
                break;
            case kind::symmetry:
                C += disk::make_vector_hho_symmetry_nitsche(op.msh, cl, op.di, G,
                        t.gn0, t.theta, op.mu, op.lam, *t.bnd);
                break;
            default: break;      // tresca / coulomb carry no consistency of their own
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
            case kind::contact_bilateral:
            {
                const vector_type u0 = vector_type::Zero(G.cols());
                J += disk::make_vector_hho_contact_jacobian(op.msh, cl, op.di, G,
                        t.gn0, t.theta, op.mu, op.lam, *t.bnd, u0,
                        /* bilateral = */ true);
                break;
            }
            case kind::symmetry:
                J += disk::make_vector_hho_symmetry_jacobian(op.msh, cl, op.di, G,
                        t.gn0, t.theta, op.mu, op.lam, *t.bnd);
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
                        t.gn0, t.theta, op.mu, op.lam, *t.bnd, u);
                R += disk::make_vector_hho_contact_rhs(op.msh, cl, op.di, G,
                        t.gn0, t.theta, op.mu, op.lam, *t.bnd, u);
                break;
            case kind::tresca:
                J += disk::make_vector_hho_tresca_jacobian(op.msh, cl, op.di, G,
                        t.gt0, t.theta, op.mu, op.lam, t.coeff, *t.bnd, u);
                R += disk::make_vector_hho_tresca_rhs(op.msh, cl, op.di, G,
                        t.gt0, t.theta, op.mu, op.lam, t.coeff, *t.bnd, u);
                break;
            case kind::coulomb:
                J += disk::make_vector_hho_coulomb_jacobian(op.msh, cl, op.di, G,
                        t.gn0, t.gt0, t.theta, op.mu, op.lam, t.coeff, *t.bnd, u, p);
                R += disk::make_vector_hho_coulomb_rhs(op.msh, cl, op.di, G,
                        t.gn0, t.gt0, t.theta, op.mu, op.lam, t.coeff, *t.bnd, u, p);
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
