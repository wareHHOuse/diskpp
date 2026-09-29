/*
 * SINUM2, Example 6.2 -- unilateral contact + Coulomb friction, with a WEAK
 * symmetry edge. 
 *
 * Reference zones on Gamma_C: separation x < 0.26, slip 0.26 < x < 0.47,
 * stick x > 0.47.
 */

#include "diskpp/loaders/loader.hpp"
#include "diskpp/mesh/meshgen.hpp"

#include "hho_params.hpp"
#include "experiment.hpp"

using namespace hho_contact;

using mesh_type = disk::generic_mesh<double, 2>;

struct Sinum62 : Experiment<mesh_type>
{
    // native tags -- VERIFY with report_native_tags() output
    enum : size_t { BOTTOM = 0, LEFT = 1, TOP = 2, RIGHT = 3 };

    bc_type bnd      { msh };   // assembler/state: Gamma_C AND Gamma_S carry face unknowns
    bc_type bnd_incr { msh };   // same, zero loads
    bc_type bnd_C    { msh };   // only Gamma_C counts as contact
    bc_type bnd_S    { msh };   // only Gamma_S counts as contact

    Sinum62() : bnd(msh), bnd_incr(msh), bnd_C(msh), bnd_S(msh) {}

    std::string name() const override { return "sinum62"; }

    void build_mesh() override
    {
        if (!prm.mesh_path.empty())
        {
            // TODO: adjust to the poly2d loader entry point in use
            disk::load_mesh_poly2d<double>(prm.mesh_path.c_str(), msh);
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

        bnd.addNeumannBC(disk::NEUMANN,        TOP,    load_top);
        bnd.addContactBC(disk::SIGNORINI_FACE, RIGHT);            // Gamma_S
        bnd.addContactBC(disk::SIGNORINI_FACE, BOTTOM);           // Gamma_C
        bnd.addNeumannBC(disk::NEUMANN,        LEFT,   load_left);

        bnd_incr.addNeumannBC(disk::NEUMANN,        TOP,    zero);
        bnd_incr.addContactBC(disk::SIGNORINI_FACE, RIGHT);
        bnd_incr.addContactBC(disk::SIGNORINI_FACE, BOTTOM);
        bnd_incr.addNeumannBC(disk::NEUMANN,        LEFT,   zero);

        bnd_C.addNeumannBC(disk::NEUMANN,        TOP,    zero);
        bnd_C.addNeumannBC(disk::NEUMANN,        RIGHT,  zero);
        bnd_C.addContactBC(disk::SIGNORINI_FACE, BOTTOM);
        bnd_C.addNeumannBC(disk::NEUMANN,        LEFT,   zero);

        bnd_S.addNeumannBC(disk::NEUMANN,        TOP,    zero);
        bnd_S.addContactBC(disk::SIGNORINI_FACE, RIGHT);
        bnd_S.addNeumannBC(disk::NEUMANN,        BOTTOM, zero);
        bnd_S.addNeumannBC(disk::NEUMANN,        LEFT,   zero);

        return { &bnd, &bnd_incr, &bnd_C };
    }

    void setup_terms(ContactTerms<mesh_type>& ct) override
    {
        const double g0 = prm.gamma_0, th = prm.theta;
        ct.add_contact (bnd_C, g0, g0, th);                 // unilateral
        ct.add_coulomb (bnd_C, g0, g0, th, prm.friction);   // friction = F -> Picard
        ct.add_symmetry(bnd_S, g0, th);                     // this benchmark only
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
