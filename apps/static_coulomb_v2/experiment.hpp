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

template<typename Mesh>
class Experiment
{
public:
    using T       = typename Mesh::coordinate_type;
    using bc_type = disk::vector_boundary_conditions<Mesh>;

    Params prm;
    Mesh   msh;

    virtual ~Experiment() = default;

    size_t signorini_tag() const
    {
        return prm.variant == trace_variant::cell ? disk::SIGNORINI_CELL
                                                  : disk::SIGNORINI_FACE;
    }

    // ---- the four hooks a benchmark implements --------------------------
    virtual std::string name() const = 0;

    virtual void build_mesh() = 0;

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
                      + "_g" + gs
                      + (prm.variant == trace_variant::cell ? "_cell" : "_face");
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

        prm.derive_lame();
        const auto wiring = setup_bc();

        std::cout << name()
                  << "  (" << msh.cells_size() << " cells, k=" << prm.k
                  << ", "
                  << (prm.variant == trace_variant::cell ? "cell" : "face")
                  << ")\n";

        /* The contact view already knows which faces are Signorini; in cell
         * mode the operator uses it to strip their unknowns.               */
        auto is_contact = [bc = wiring.contact]
                          (const typename Mesh::face_type& fc)
                          { return bc->is_contact_face(fc); };

        spatial_operator<Mesh> op(msh, prm.k, prm.mu, prm.lam,
                                  /* with_mass = */ false,
                                  prm.variant, is_contact);

        ContactTerms<Mesh> ct(op);
        setup_terms(ct);

        auto body = [this](const auto& pt) { return this->body_force(pt); };
        NewtonSolver<Mesh> solver(op, ct, *wiring.state, *wiring.incr, body);

        auto st = solver.make_state();
        solve_stats stats;
        const bool ok = solver.solve(st, prm.opts, stats);


        std::ios::fmtflags sf(std::cout.flags());
        const auto sprec = std::cout.precision();
        std::cout << std::scientific << std::setprecision(3);

        std::cout << "SUMMARY"
                  << " name=" << name()
                  << " level=" << prm.mesh_level
                  << " k=" << prm.k
                  << " theta=" << prm.theta
                  << " gamma_0=" << prm.gamma_0
                  << " gamma_n_0=" << prm.gn0()
                  << " gamma_t_0=" << prm.gt0()
                  << " friction=" << prm.friction
                  << " converged=" << (stats.converged ? 1 : 0)
                  << " picard=" << stats.picard_iters
                  << " newton_total=" << stats.newton_total
                  << " newton_per_picard=";
        for (size_t i = 0; i < stats.newton_per_picard.size(); i++)
            std::cout << stats.newton_per_picard[i]
                      << (i + 1 < stats.newton_per_picard.size() ? "," : "");
        std::cout << "\n";

        std::cout.flags(sf);
        std::cout.precision(sprec);

        ContactReport<Mesh> rep(solver, op, *wiring.contact,
                                prm.gn0(), prm.gt0());
        report(rep, st, ct.needs_picard());

        return ok ? 0 : 1;
    }
};

} // namespace hho_contact