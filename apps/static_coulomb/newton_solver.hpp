#pragma once

#include <vector>
#include <utility>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "diskpp/methods/hho"
#include "diskpp/solvers/direct_solvers.hpp"

#include "spatial_operator.hpp"
#include "contact_composer.hpp"
#include "hho_params.hpp"

namespace hho_contact {

template<typename T>
struct hho_state
{
    std::vector<disk::dynamic_vector<T>> cells;
    disk::dynamic_vector<T>              faces;
};

/* Per-solve iteration record, filled by solve().  For a Newton-only
 * problem picard_iters = 1 and newton_per_picard has a single entry.     */
struct solve_stats
{
    bool                converged   = false;
    size_t              picard_iters = 0;   // outer fixed-point passes
    size_t              newton_total = 0;   // sum over all passes
    std::vector<size_t> newton_per_picard;  // one entry per Picard pass
};

template<typename Mesh>
class NewtonSolver
{
public:
    using T           = typename Mesh::coordinate_type;
    using matrix_type = disk::dynamic_matrix<T>;
    using vector_type = disk::dynamic_vector<T>;
    using cell_type   = typename Mesh::cell_type;
    using bc_type     = disk::vector_boundary_conditions<Mesh>;

    using assembler_type =
        decltype(disk::make_vector_primal_hho_assembler(
                     std::declval<const Mesh&>(),
                     std::declval<const disk::MeshDegreeInfo<Mesh>&>(),
                     std::declval<const bc_type&>()));

private:
    const spatial_operator<Mesh>& op;
    const ContactTerms<Mesh>&     ct;
    const bc_type&                bnd;        // state gather (loads, faces)
    const bc_type&                bnd_incr;   // increment assembly / scatter

    assembler_type assembler;

    std::vector<matrix_type> A;    // K - consistency + linear penalties   (cached)
    std::vector<vector_type> L;    // body + Neumann loads                 (cached)
    std::vector<matrix_type> AL;   // per-iteration condensation storage
    std::vector<vector_type> bL;

public:
    template<typename BodyF>
    NewtonSolver(const spatial_operator<Mesh>& op_,
                 const ContactTerms<Mesh>&     ct_,
                 const bc_type& bnd_state, const bc_type& bnd_increment,
                 BodyF&& body)
        : op(op_), ct(ct_), bnd(bnd_state), bnd_incr(bnd_increment),
          assembler(disk::make_vector_primal_hho_assembler(op_.msh, op_.di,
                                                           bnd_state))
    {
        const size_t nc = op.num_cells();
        A.resize(nc); L.resize(nc); AL.resize(nc); bL.resize(nc);

        size_t ci = 0;
        for (auto& cl : op.msh)
        {
            // ---- constant system matrix: every linear term leaves the loop
            A[ci] = op.K[ci]
                    - ct.consistency(cl, ci)
                    + ct.linear_jacobian(cl, ci);

            // ---- loads
            auto cb = disk::make_vector_monomial_basis(op.msh, cl, op.cell_deg);
            L[ci]   = vector_type::Zero(op.total_dofs(cl));
            L[ci].head(op.cbs) = disk::make_rhs(op.msh, cl, cb, body);

            size_t ofs = op.cbs;
            for (auto& fc : disk::faces(op.msh, cl))
            {
                const auto   fdi = op.di.degreeInfo(op.msh, fc);
                const size_t nfd = disk::vector_face_dofs(op.msh, fdi);
                if (bnd.is_neumann_face(fc))
                {
                    auto fb = disk::make_vector_monomial_basis(op.msh, fc,
                                                               fdi.degree());
                    L[ci].segment(ofs, nfd) += disk::make_rhs(op.msh, fc, fb,
                        bnd.neumann_boundary_func(op.msh.lookup(fc)));
                }
                ofs += nfd;
            }
            ci++;
        }
    }

    // ---- state ----------------------------------------------------------
    hho_state<T> make_state() const
    {
        hho_state<T> st;
        st.cells.assign(op.num_cells(), vector_type::Zero(op.cbs));
        st.faces = vector_type::Zero(assembler.LHS.rows());
        return st;
    }

    vector_type gather(const cell_type& cl, size_t ci,
                       const hho_state<T>& st) const
    {
        vector_type u = vector_type::Zero(A[ci].rows());
        u.head(op.cbs) = st.cells[ci];
        vector_type uF = assembler.take_local_solution(op.msh, cl, bnd, st.faces);
        u.tail(uF.size()) = uF;
        return u;
    }

    // ---- driver entry point ---------------------------------------------
    bool solve(hho_state<T>& st, const solver_opts& o)
    {
        solve_stats stats;
        return solve(st, o, stats);
    }

    bool solve(hho_state<T>& st, const solver_opts& o, solve_stats& stats)
    {
        stats = solve_stats{};

        if (!ct.needs_picard())
        {
            size_t nit = 0;
            const bool ok = newton_loop(st, nullptr, o, nit, o.verbose);
            stats.converged    = ok;
            stats.picard_iters = 1;
            stats.newton_total = nit;
            stats.newton_per_picard = { nit };
            if (o.verbose)
            {
                if (!ok) std::cout << "  *** NEWTON DID NOT CONVERGE ***\n";
                std::cout << "  Newton iterations: " << nit << "\n";
            }
            return ok;
        }

        // ---- Picard on the frozen Coulomb threshold ----------------------
        hho_state<T> frozen = st;
        size_t total_newton = 0;
        bool   picard_ok    = false;

        for (size_t pit = 0; pit < o.max_picard; pit++)
        {
            size_t nit = 0;
            const bool newton_ok = newton_loop(st, &frozen, o, nit, false);
            total_newton += nit;
            stats.newton_per_picard.push_back(nit);

            T diff2 = (st.faces - frozen.faces).squaredNorm();
            for (size_t ci = 0; ci < st.cells.size(); ci++)
                diff2 += (st.cells[ci] - frozen.cells[ci]).squaredNorm();
            const T pdiff = std::sqrt(diff2);

            if (o.verbose)
                std::cout << "  picard " << std::setw(2) << pit
                          << "   newton " << std::setw(2) << nit
                          << (newton_ok ? " " : "*")
                          << "   |u - u_picard| = " << std::scientific
                          << std::setprecision(3) << pdiff << "\n";

            frozen = st;

            if (pdiff < o.tol_picard) { picard_ok = true; break; }
        }

        stats.converged    = picard_ok;
        stats.picard_iters = stats.newton_per_picard.size();
        stats.newton_total = total_newton;

        if (o.verbose)
        {
            if (!picard_ok) std::cout << "  *** PICARD DID NOT CONVERGE ***\n";
            std::cout << "  total inner Newton iterations: "
                      << total_newton << "\n";
        }
        return picard_ok;
    }

private:
    /* One condensed Newton step; returns (|du_faces|, |u_faces|).        */
    std::pair<T,T> step(hho_state<T>& st, const hho_state<T>* frozen)
    {
        const bool nonlin = ct.has_nonlinear();

        assembler.initialize();

        size_t ci = 0;
        for (auto& cl : op.msh)
        {
            vector_type u_loc = gather(cl, ci, st);

            matrix_type J = A[ci];
            vector_type R = A[ci]*u_loc - L[ci];

            if (nonlin)
            {
                // Coulomb reads p; everyone else ignores it, so avoid a
                // second gather unless an outer iteration is running.
                if (frozen)
                {
                    vector_type p_loc = gather(cl, ci, *frozen);
                    ct.accumulate_nonlinear(cl, ci, u_loc, p_loc, J, R);
                }
                else
                    ct.accumulate_nonlinear(cl, ci, u_loc, u_loc, J, R);
            }

            vector_type negR = -R;
            auto scT = disk::make_vector_static_condensation_withMatrix(
                           op.msh, cl, op.di, J, negR);
            auto sc  = std::get<0>(scT);
            AL[ci]   = std::get<1>(scT);
            bL[ci]   = std::get<2>(scT);

            assembler.assemble(op.msh, cl, bnd_incr, sc.first, sc.second);
            ci++;
        }
        assembler.finalize();

        vector_type du_faces = vector_type::Zero(assembler.LHS.rows());
#ifdef HAVE_PARDISO
        disk::solvers::sparse_lu(assembler.LHS, assembler.RHS, du_faces,
                                 disk::solvers::direct_solver::pardiso);
#else
        disk::solvers::sparse_lu(assembler.LHS, assembler.RHS, du_faces,
                                 disk::solvers::direct_solver::sparselu);
#endif

        ci = 0;
        for (auto& cl : op.msh)
        {
            vector_type duF = assembler.take_local_solution(op.msh, cl,
                                                             bnd_incr, du_faces);
            st.cells[ci] += bL[ci] - AL[ci]*duF;
            ci++;
        }
        st.faces += du_faces;

        return { du_faces.norm(), st.faces.norm() };
    }

    bool newton_loop(hho_state<T>& st, const hho_state<T>* frozen,
                     const solver_opts& o, size_t& nit, bool trace)
    {
        for (nit = 0; nit < o.max_newton; nit++)
        {
            auto [du, u] = step(st, frozen);

            if (trace)
                std::cout << "  newton " << std::setw(2) << nit
                          << "   |du| = " << std::scientific
                          << std::setprecision(3) << du << "\n";

            // absolute OR relative: whichever is looser at the current scale
            if (du < o.tol_newton || du < o.tol_newton_rel*(u + 1e-14))
            { nit++; return true; }
        }
        return false;
    }
};

} // namespace hho_contact