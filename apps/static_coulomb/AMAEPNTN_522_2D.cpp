/*
 * AMAEPNTN Test 5.2.2 in 2D  (small-strain short-sliding analogue)
 */

#include "diskpp/loaders/loader.hpp"
#include "diskpp/mesh/meshgen.hpp"

#include "hho_params.hpp"
#include "experiment.hpp"

#include <memory>

using namespace hho_contact;

using mesh_type = disk::generic_mesh<double, 2>;

struct small_strain_short_sliding : Experiment<mesh_type>
{
    enum : size_t { BOTTOM = 0, LEFT = 1, TOP = 2, RIGHT = 3 };

    // Built in setup_bc(), AFTER build_mesh() has populated msh.
    std::unique_ptr<bc_type> bnd;        // state: real loads
    std::unique_ptr<bc_type> bnd_incr;   // increment: zero loads
    std::unique_ptr<bc_type> bnd_C;      // only Gamma_C counts as contact

    std::string name() const override { return "small_strain_short_sliding"; }

    void build_mesh() override
    {
        if (!prm.mesh_path.empty())
        {
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
        bnd      = std::make_unique<bc_type>(msh);   // msh already loaded here
        bnd_incr = std::make_unique<bc_type>(msh);
        bnd_C    = std::make_unique<bc_type>(msh);

        auto zero = [](const auto&) {
            Eigen::Matrix<double,2,1> z; z.setZero(); return z; };

        // press-then-slide on the top edge (small strain: u_t~0.03, u_n~0.01)
        auto depl_top = [](const auto&) {
            Eigen::Matrix<double,2,1> g;
            g(0) =  0.03;    // tangential slide (+x)
            g(1) = -0.01;    // normal press toward rigid plane at BOTTOM
            return g; };

        bnd->addNeumannBC(disk::NEUMANN, RIGHT, zero);
        bnd->addNeumannBC(disk::NEUMANN, LEFT,  zero);
        bnd->addDirichletBC(disk::DIRICHLET, TOP, depl_top);
        bnd->addContactBC(disk::SIGNORINI_FACE, BOTTOM);

        bnd_incr->addNeumannBC(disk::NEUMANN, RIGHT, zero);
        bnd_incr->addNeumannBC(disk::NEUMANN, LEFT,  zero);
        bnd_incr->addDirichletBC(disk::DIRICHLET, TOP, zero);
        bnd_incr->addContactBC(disk::SIGNORINI_FACE, BOTTOM);

        bnd_C->addNeumannBC(disk::NEUMANN, RIGHT, zero);
        bnd_C->addNeumannBC(disk::NEUMANN, LEFT,  zero);
        bnd_C->addDirichletBC(disk::DIRICHLET, TOP, zero);
        bnd_C->addContactBC(disk::SIGNORINI_FACE, BOTTOM);

        return { bnd.get(), bnd_incr.get(), bnd_C.get() };
    }

    void setup_terms(ContactTerms<mesh_type>& ct) override
    {
        const double g0 = prm.gamma_0, th = prm.theta;
        ct.add_contact (*bnd_C, g0, g0, th);                 // unilateral
        ct.add_coulomb (*bnd_C, g0, g0, th, prm.friction);   // friction = F -> Picard
    }
};

int main(int argc, char** argv)
{
    small_strain_short_sliding exp;
    exp.prm = (argc > 1) ? Params::load(argv[1]) : Params{};

    if (exp.prm.E == 0.0 && exp.prm.mu == 0.0)
    { exp.prm.E = 1.0e8; exp.prm.nu = 0.3; }
    if (exp.prm.friction == 0.0) exp.prm.friction = 0.1;

    return exp.run();
}