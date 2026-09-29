/*
 * SINUM2, Example 6.2 -- unilateral contact + Coulomb friction, with a WEAK
 * symmetry edge. 
 *
 * Reference zones on Gamma_C: separation x < 0.26, slip 0.26 < x < 0.47,
 * stick x > 0.47.
 */

#include <memory>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/mesh/meshgen.hpp"

#include "hho_params.hpp"
#include "experiment.hpp"

using namespace hho_contact;

using mesh_type = disk::generic_mesh<double, 2>;

struct Sinum62 : Experiment<mesh_type>
{
    // native tags
    //   tag 0 (0.515625, 0) BOTTOM   tag 1 (1, 0.5) RIGHT
    //   tag 2 (0.515625, 1) TOP      tag 3 (0, 0.5) LEFT
    enum : size_t { BOTTOM = 0, RIGHT = 1, TOP = 2, LEFT = 3 };
    
    std::unique_ptr<bc_type> bnd;       // assembler/state: Gamma_C AND Gamma_S carry face unknowns
    std::unique_ptr<bc_type> bnd_incr;  // same, zero loads
    std::unique_ptr<bc_type> bnd_C;     // only Gamma_C counts as contact
    std::unique_ptr<bc_type> bnd_S;     // only Gamma_S counts as contact

    std::string name() const override { return "sinum62"; }

    void build_mesh() override
    {
        if (!prm.mesh_path.empty())
        {
            // TODO: adjust to the poly2d loader entry point in use
            disk::load_mesh_poly<double>(prm.mesh_path.c_str(), msh);
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

        auto fric = [F = prm.friction](const auto&) { return F; };
        auto no_fric = [](const auto&){ return 0.0; };
        auto no_gap = [](const auto&, const auto&) { return 0.0; };   

        bnd      = std::make_unique<bc_type>(msh);
        bnd_incr = std::make_unique<bc_type>(msh);
        bnd_C    = std::make_unique<bc_type>(msh);
        bnd_S    = std::make_unique<bc_type>(msh);

        bnd->addNeumannBC(disk::NEUMANN,        TOP,    load_top);
        bnd->addContactBC(disk::SIGNORINI_FACE, RIGHT, no_fric, no_gap);            // Gamma_S
        bnd->addContactBC(signorini_tag(),     BOTTOM, fric, no_gap);           // Gamma_C
        bnd->addNeumannBC(disk::NEUMANN,        LEFT,   load_left);

        bnd_incr->addNeumannBC(disk::NEUMANN,        TOP,    zero);
        bnd_incr->addContactBC(disk::SIGNORINI_FACE, RIGHT, no_fric, no_gap);
        bnd_incr->addContactBC(signorini_tag(),     BOTTOM, fric, no_gap);
        bnd_incr->addNeumannBC(disk::NEUMANN,        LEFT,   zero);

        bnd_C->addNeumannBC(disk::NEUMANN,        TOP,    zero);
        bnd_C->addNeumannBC(disk::NEUMANN,        RIGHT,  zero);
        bnd_C->addContactBC(signorini_tag(),     BOTTOM, fric, no_gap);
        bnd_C->addNeumannBC(disk::NEUMANN,        LEFT,   zero);

        bnd_S->addNeumannBC(disk::NEUMANN,        TOP,    zero);
        bnd_S->addContactBC(disk::SIGNORINI_FACE, RIGHT, no_fric, no_gap);
        bnd_S->addNeumannBC(disk::NEUMANN,        BOTTOM, zero);
        bnd_S->addNeumannBC(disk::NEUMANN,        LEFT,   zero);

        return { bnd.get(), bnd_incr.get(), bnd_C.get() };
    }

    void setup_terms(ContactTerms<mesh_type>& ct) override
    {
        const double gn0 = prm.gn0(), gt0 = prm.gt0(), th = prm.theta;
        const double gs0 = prm.gs0();   // weak-symmetry edge penalty
        ct.add_contact (*bnd_C, gn0, gt0, th);                 // unilateral
        ct.add_coulomb (*bnd_C, gn0, gt0, th, prm.friction);   // friction = F -> Picard
        ct.add_symmetry(*bnd_S, gs0, th);                     // this benchmark only
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
