#pragma once

#include <vector>
#include <utility>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <string>

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

/* Per-solve iteration record, filled by solve()    */
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

    // ---- assemble the local DOF solution (uT,uF) ------------------------
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
        const auto t_start = std::chrono::steady_clock::now();

        if (!ct.needs_picard())
        {
            size_t nit = 0;
            const bool ok = newton_loop(st, nullptr, o, nit);
            stats.converged    = ok;
            stats.picard_iters = 1;
            stats.newton_total = nit;
            stats.newton_per_picard = { nit };

            report_summary(stats, ok, T(0), o, t_start, /* picard = */ false);
            return ok;
        }

        // ---- Picard on the frozen Coulomb threshold ----------------------
        hho_state<T> frozen = st;
        size_t total_newton = 0;
        size_t newton_failed = 0;
        bool   picard_ok    = false;
        T      last_pdiff   = T(0);

        for (size_t pit = 0; pit < o.max_picard; pit++)
        {
            std::cout << "\n*** Picard pass " << pit << " ***\n";
            size_t nit = 0;
            const bool newton_ok = newton_loop(st, &frozen, o, nit);
            total_newton += nit;
            stats.newton_per_picard.push_back(nit);

            T diff2 = (st.faces - frozen.faces).squaredNorm();
            T norm2 = st.faces.squaredNorm();
            for (size_t ci = 0; ci < st.cells.size(); ci++)
            {
                diff2 += (st.cells[ci] - frozen.cells[ci]).squaredNorm();
                norm2 += st.cells[ci].squaredNorm();
            }
            const T abs_err = std::sqrt(diff2);
            const T rel_err = abs_err / (std::sqrt(norm2) + 1e-14);

            if (!newton_ok) newton_failed++;
            last_pdiff = abs_err;

            {
                std::ios::fmtflags pf(std::cout.flags());
                const auto pprec = std::cout.precision();
                std::cout << "Picard " << pit
                          << ": |u - u_picard| = " << std::scientific
                          << std::setprecision(5) << abs_err
                          << "   relative = " << rel_err
                          << "   (inner Newton: " << nit << " iters, "
                          << (newton_ok ? "converged" : "NOT converged") << ")\n";
                std::cout.flags(pf);
                std::cout.precision(pprec);
            }

            frozen = st;

            if (abs_err < o.tol_picard && newton_ok) { picard_ok = true; break; }
        }

        stats.converged    = picard_ok;
        stats.picard_iters = stats.newton_per_picard.size();
        stats.newton_total = total_newton;

        report_summary(stats, picard_ok, last_pdiff, o, t_start,
                       /* picard = */ true, newton_failed);
        return picard_ok;
    }

private:
    /* Closing block, in the style of
     * libdiskpp mechanics/NewtonSolver :: SolverInfo::printInfo().       */
    void report_summary(const solve_stats& stats, bool ok, T last_pdiff,
                        const solver_opts& o,
                        std::chrono::steady_clock::time_point t_start,
                        bool picard, size_t newton_failed = 0) const
    {
        const double secs = std::chrono::duration<double>(
                                std::chrono::steady_clock::now() - t_start).count();

        std::ios::fmtflags f(std::cout.flags());
        const auto prec = std::cout.precision();

        std::cout << "\n";
        std::cout << "------------------------------------------------------- \n";
        std::cout << "Summaring: \n";

        if (picard)
        {
            std::cout << "Total Newton's iterations: " << stats.newton_total
                      << " in " << stats.picard_iters << " Picard passes\n";
            std::cout << "Picard: " << (ok ? "converged" : "NOT converged")
                      << "   |u - u_picard| = " << std::scientific
                      << std::setprecision(5) << last_pdiff
                      << "   (tol " << o.tol_picard << ")\n";
            std::cout << "Newton: "
                      << (stats.picard_iters - newton_failed) << " of "
                      << stats.picard_iters << " passes converged";
            if (newton_failed > 0)
                std::cout << "   <-- " << newton_failed
                          << " pass(es) hit max_newton = " << o.max_newton;
            std::cout << "\n";
        }
        else
        {
            std::cout << "Total Newton's iterations: " << stats.newton_total
                      << "\n";
            std::cout << "Newton: " << (ok ? "converged" : "NOT converged")
                      << "   (tol " << std::scientific << std::setprecision(5)
                      << o.tol_newton << ", max_newton = " << o.max_newton
                      << ")\n";
        }

        std::cout << "Total time to solve the problem: " << std::fixed
                  << std::setprecision(3) << secs << " sec\n";
        std::cout << "------------------------------------------------------- \n";
        std::cout << " \n";

        std::cout.flags(f);
        std::cout.precision(prec);
    }

public:

private:
    /* What one Newton step reports back: the increment it applied, the
     * state norm it applied it to, and the residual norm it was formed
     * from (i.e. the residual AT the incoming iterate).                  */
    struct step_info { T du; T u; T residual; };

    /* One condensed Newton step. */
    step_info step(hho_state<T>& st, const hho_state<T>* frozen)
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

        /* Residual of the system actually being solved: the assembled rhs
         * is the CONDENSED -R, with the Dirichlet rows already eliminated.
         * The raw per-cell R is not usable here -- it still carries the
         * Dirichlet reactions, which are large and never vanish.         */
        const T res_norm = assembler.RHS.norm();

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

        return { du_faces.norm(), st.faces.norm(), res_norm };
    }

    // ---- reporting, after libdiskpp mechanics/NewtonSolver -------------
    static void table_rule()
    {
        std::cout << "------------------------------------------------------"
                     "---------------------\n";
    }

    static void table_header()
    {
        table_rule();
        std::cout << "| Iteration | Norme l2 incr | Relative incr |  "
                     "Residual l2  | Relative res  |\n";
        table_rule();
    }

    static void table_row(size_t iter, T incr, T rel_incr, T res, T rel_res)
    {
        std::ios::fmtflags f(std::cout.flags());
        const auto prec = std::cout.precision();
        std::cout.precision(5);
        std::cout.setf(std::iostream::scientific, std::iostream::floatfield);

        std::string s_iter = "   " + std::to_string(iter) + "            ";
        s_iter.resize(9);

        std::cout << "| " << s_iter << " |   " << incr
                  << " |   " << rel_incr
                  << " |   " << res
                  << " |   " << rel_res << " |\n";
        std::cout.flags(f);
        std::cout.precision(prec);
    }

    bool newton_loop(hho_state<T>& st, const hho_state<T>* frozen,
                     const solver_opts& o, size_t& nit)
    {
        T res0 = 0;
        table_header();

        for (nit = 0; nit < o.max_newton; nit++)
        {
            const auto si = step(st, frozen);

            const T abs_err = si.du;
            const T rel_err = si.du / (si.u + 1e-14);

            if (nit == 0) res0 = si.residual;
            const T rel_res = (res0 > 1e-30) ? si.residual / res0 : T(1);

            table_row(nit, abs_err, rel_err, si.residual, rel_res);

            // absolute OR relative: whichever is looser at the current scale
            if (abs_err < o.tol_newton || rel_err < o.tol_newton_rel)
            { nit++; table_rule(); return true; }
        }
        table_rule();
        return false;
    }
};

} // end namespace hho_contact
