/*
 * SINUM2, Example 6.2 -- unilateral contact + Coulomb friction, with a WEAK
 * symmetry edge.  Refactored driver.
 *
 * The composition below is the whole difference against Bostan-Han:
 *   - contact is unilateral (nonlinear, [.]_- active)
 *   - a Coulomb term is added  -> ContactTerms flips needs_picard() and the
 *     solver wraps Newton in the outer fixed-point automatically
 *   - the symmetry term is added EXPLICITLY with its own view: it belongs
 *     to this benchmark's geometry, not to the Coulomb model.  Being
 *     linear, it is folded into the cached A and costs nothing per
 *     iteration.
 *
 * Reference zones on Gamma_C: separation x < 0.26, slip 0.26 < x < 0.47,
 * stick x > 0.47.
 */

#include "diskpp/loaders/loader.hpp"
#include "diskpp/mesh/meshgen.hpp"

#include <memory>

#include "hho_params.hpp"
#include "experiment.hpp"

using namespace hho_contact;

using mesh_type = disk::generic_mesh<double, 2>;

struct Sinum62 : Experiment<mesh_type>
{
    // native tags, confirmed from report_native_tags(): 0=bottom 1=right 2=top 3=left
    enum : size_t { BOTTOM = 0, RIGHT = 1, TOP = 2, LEFT = 3 };

    // views are built LAZILY in setup_bc(), i.e. AFTER build_mesh(): a
    // vector_boundary_conditions bound to an empty mesh throws on addBC.
    std::unique_ptr<bc_type> bnd;      // state: Gamma_C AND Gamma_S carry face unknowns
    std::unique_ptr<bc_type> bnd_incr; // same, zero loads
    std::unique_ptr<bc_type> bnd_C;    // only Gamma_C counts as contact
    std::unique_ptr<bc_type> bnd_S;    // only Gamma_S counts as contact

    std::string name() const override { return "sinum62"; }

    void build_mesh() override
    {
        if (!prm.mesh_path.empty())
        {
            // TODO: adjust to the poly2d loader entry point in use
            disk::load_mesh_poly(prm.mesh_path.c_str(), msh);
        }
        else
        {
            auto mesher = disk::make_fvca5_hex_mesher(msh);
            mesher.make_level(prm.mesh_level);
        }
    }

    bc_wiring setup_bc() override
    {
        auto zero = [](const auto&) {
            Eigen::Matrix<double,2,1> z; z.setZero(); return z; };
        // inward load on the right half of the top edge: outward n = (0,1)
        auto load_top = [](const auto& pt) {
            Eigen::Matrix<double,2,1> g; g.setZero();
            if (pt.x() > 0.5) g(1) = -1.0;
            return g; };
        // inward load on the upper half of the left edge: outward n = (-1,0)
        auto load_left = [](const auto& pt) {
            Eigen::Matrix<double,2,1> g; g.setZero();
            if (pt.y() > 0.5) g(0) = 1.0;
            return g; };

        bnd      = std::make_unique<bc_type>(msh);
        bnd_incr = std::make_unique<bc_type>(msh);
        bnd_C    = std::make_unique<bc_type>(msh);
        bnd_S    = std::make_unique<bc_type>(msh);

        bnd->addNeumannBC(disk::NEUMANN,        TOP,    load_top);
        bnd->addContactBC(disk::SIGNORINI_FACE, RIGHT);            // Gamma_S
        bnd->addContactBC(disk::SIGNORINI_FACE, BOTTOM);           // Gamma_C
        bnd->addNeumannBC(disk::NEUMANN,        LEFT,   load_left);

        bnd_incr->addNeumannBC(disk::NEUMANN,        TOP,    zero);
        bnd_incr->addContactBC(disk::SIGNORINI_FACE, RIGHT);
        bnd_incr->addContactBC(disk::SIGNORINI_FACE, BOTTOM);
        bnd_incr->addNeumannBC(disk::NEUMANN,        LEFT,   zero);

        bnd_C->addNeumannBC(disk::NEUMANN,        TOP,    zero);
        bnd_C->addNeumannBC(disk::NEUMANN,        RIGHT,  zero);
        bnd_C->addContactBC(disk::SIGNORINI_FACE, BOTTOM);
        bnd_C->addNeumannBC(disk::NEUMANN,        LEFT,   zero);

        bnd_S->addNeumannBC(disk::NEUMANN,        TOP,    zero);
        bnd_S->addContactBC(disk::SIGNORINI_FACE, RIGHT);
        bnd_S->addNeumannBC(disk::NEUMANN,        BOTTOM, zero);
        bnd_S->addNeumannBC(disk::NEUMANN,        LEFT,   zero);

        return { bnd.get(), bnd_incr.get(), bnd_C.get() };
    }

    void setup_terms(ContactTerms<mesh_type>& ct) override
    {
        const double g0 = prm.gamma_0, th = prm.theta;
        ct.add_contact (*bnd_C, g0, g0, th);                 // unilateral
        ct.add_coulomb (*bnd_C, g0, g0, th, prm.friction);   // friction = F -> Picard
        ct.add_symmetry(*bnd_S, g0, th);                     // this benchmark only
    }
};

int main(int argc, char** argv)
{
    Sinum62 exp;
    exp.prm = (argc > 1) ? Params::load(argv[1]) : Params{};

    // benchmark defaults when no .dat overrides them
    if (exp.prm.E == 0.0 && exp.prm.mu == 0.0)
    { exp.prm.E = 1.0e4; exp.prm.nu = 0.2; }
    if (exp.prm.friction == 0.0) exp.prm.friction = 0.5;

    return exp.run();
}