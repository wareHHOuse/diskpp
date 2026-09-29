/*
 * experiment.hpp : where a benchmark is set up.
 *
 * Experiment owns Params and the mesh, and wires the pipeline in run():
 *
 *     build_mesh -> setup_bc -> spatial_operator -> ContactTerms
 *       -> NewtonSolver -> solve -> ContactReport
 */

#pragma once

#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <filesystem>
#include <map>
#include <iostream>

#include "hho_params.hpp"
#include "spatial_operator.hpp"
#include "contact_composer.hpp"
#include "newton_solver.hpp"
#include "contact_report.hpp"

namespace hho_contact {

/* Diagnostic: centroid of each native boundary id.  y ~ 0 is bottom,
 * x ~ 0 is left, and so on.  NOTE: adjust the boundary-id accessor if the
 * diskpp revision at hand names it differently.                           */
template<typename Mesh>
void report_native_tags(const Mesh& msh)
{
    using T = typename Mesh::coordinate_type;
    constexpr size_t DIM = Mesh::dimension;

    std::map<size_t, std::pair<std::array<T,3>, size_t>> acc;

    for (auto itor = msh.boundary_faces_begin();
         itor != msh.boundary_faces_end(); itor++)
    {
        const auto& fc = *itor;
        const auto  bi = msh.boundary_id(fc);          // native tag
        const auto  b  = disk::barycenter(msh, fc);

        auto& [sum, cnt] = acc[bi];
        sum[0] += b.x();
        sum[1] += b.y();
        if constexpr (DIM == 3) sum[2] += b.z();
        cnt++;
    }

    std::cout << "native boundary tags:\n";
    for (auto& [id, sc] : acc)
    {
        const auto& [sum, cnt] = sc;
        std::cout << "  tag " << id << " (" << cnt << " faces)  centroid ( "
                  << sum[0]/cnt << " , " << sum[1]/cnt;
        if constexpr (DIM == 3) std::cout << " , " << sum[2]/cnt;
        std::cout << " )\n";
    }
}

template<typename Mesh>
class Experiment
{
public:
    using T       = typename Mesh::coordinate_type;
    using bc_type = disk::vector_boundary_conditions<Mesh>;

    Params prm;
    Mesh   msh;

    virtual ~Experiment() = default;

    // ---- the four hooks a benchmark implements --------------------------
    virtual std::string name() const = 0;

    virtual void build_mesh() = 0;

    /* Fill the boundary views and return the wiring the solver needs.
     * Views are OWNED by the concrete experiment (members), because their
     * number varies: Bostan-Han uses 2, SINUM 6.2 uses 4.                */
    struct bc_wiring
    {
        const bc_type* state;     // assembler + state gather (real loads)
        const bc_type* incr;      // increment assembly       (zero loads)
        const bc_type* contact;   // which faces the report treats as contact
    };
    virtual bc_wiring setup_bc() = 0;

    virtual void setup_terms(ContactTerms<Mesh>& ct) = 0;

    /* Filename stem encoding the case: <name>_L<level>_k<k>_g<gamma0>.
     * gamma0 is formatted compactly (e.g. 1e+03 -> g1e3).                 */
    std::string case_tag() const
    {
        std::ostringstream g;
        g << std::scientific << std::setprecision(0) << prm.gamma_0; // e.g. 1e+03
        std::string gs = g.str();
        gs.erase(std::remove(gs.begin(), gs.end(), '+'), gs.end());   // 1e03
        // strip leading zero in exponent: 1e03 -> 1e3
        auto epos = gs.find('e');
        if (epos != std::string::npos)
            while (epos + 1 < gs.size() && gs[epos+1] == '0' && epos + 2 < gs.size())
                gs.erase(epos + 1, 1);
        return name() + "_L" + std::to_string(prm.mesh_level)
                      + "_k" + std::to_string(prm.k)
                      + "_g" + gs;
    }

    /* Output path stem: paraview-data/<name>/<case_tag>, creating the
     * per-example directory if needed.  All writers append their own
     * extension (.vtu, _contact.vtu, .csv).                              */
    std::string output_stem() const
    {
        namespace fs = std::filesystem;
        const fs::path dir = fs::path("paraview-data") / name();
        std::error_code ec;
        fs::create_directories(dir, ec);   // no-op if it already exists
        if (ec)
            std::cout << "  warning: could not create " << dir
                      << " (" << ec.message() << "), writing to cwd\n";
        return (dir / case_tag()).string();
    }

    /* Model-specific output; default dispatches on the term list.        */
    virtual void report(const ContactReport<Mesh>& rep,
                        const hho_state<T>& st,
                        bool coulomb_mode)
    {
        const std::string stem = output_stem();
        rep.write_bulk(stem, st);
        if (coulomb_mode)
        {
            rep.write_coulomb(stem, st, prm.friction);
            rep.write_csv(stem + ".csv", st, prm.friction, /*coulomb=*/true);
        }
        else
        {
            rep.write_tresca(stem, st, prm.friction);
            rep.write_csv(stem + ".csv", st, prm.friction, /*coulomb=*/false);
        }
    }

    virtual Eigen::Matrix<T, Mesh::dimension, 1>
    body_force(const typename Mesh::point_type&) const
    {
        Eigen::Matrix<T, Mesh::dimension, 1> f; f.setZero(); return f;
    }

    // ---- common pipeline ------------------------------------------------
    int run()
    {
        build_mesh();
        report_native_tags(msh);

        prm.derive_lame();
        const auto wiring = setup_bc();

        std::cout << name() << "\n"
                  << "  " << msh.cells_size() << " cells   k=" << prm.k
                  << "  friction=" << prm.friction
                  << "  theta=" << prm.theta
                  << "  gamma_0=" << prm.gamma_0 << "\n"
                  << "  mu = " << prm.mu << "   lambda = " << prm.lam << "\n";

        spatial_operator<Mesh> op(msh, prm.k, prm.mu, prm.lam);

        ContactTerms<Mesh> ct(op);
        setup_terms(ct);

        auto body = [this](const auto& pt) { return this->body_force(pt); };
        NewtonSolver<Mesh> solver(op, ct, *wiring.state, *wiring.incr, body);

        auto st = solver.make_state();
        solve_stats stats;
        const bool ok = solver.solve(st, prm.opts, stats);

        // one machine-parseable line for gamma_0 / convergence sweeps
        std::cout << "SUMMARY"
                  << " name=" << name()
                  << " level=" << prm.mesh_level
                  << " k=" << prm.k
                  << " theta=" << prm.theta
                  << " gamma_0=" << prm.gamma_0
                  << " friction=" << prm.friction
                  << " converged=" << (stats.converged ? 1 : 0)
                  << " picard=" << stats.picard_iters
                  << " newton_total=" << stats.newton_total
                  << " newton_per_picard=";
        for (size_t i = 0; i < stats.newton_per_picard.size(); i++)
            std::cout << stats.newton_per_picard[i]
                      << (i + 1 < stats.newton_per_picard.size() ? "," : "");
        std::cout << "\n";

        ContactReport<Mesh> rep(solver, op, *wiring.contact,
                                prm.gamma_0, prm.gamma_0);
        report(rep, st, ct.needs_picard());

        return ok ? 0 : 1;
    }
};

} // namespace hho_contact  