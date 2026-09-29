/*
 * Rahmi et al. 2023, §4.2.2  --  frictional full-sliding contact (2D plane strain)
 *
 *   Foundation  = (0,48) x (0,30) mm,  E=100 MPa, nu=0.3
 *   Rigid flat punch = 2a=6 mm, centred at x=24 on the top edge -> Gamma_C
 *
 * Faithful mapping (fixed punch, driven base):
 *   contact : top strip x in [21,27], y=30   Signorini + Coulomb vs fixed punch
 *   base    : y=0                             Dirichlet (u_x=slide, u_y=+indent)
 *   free    : sides + top outside the strip   homogeneous Neumann
 *
 * Single solve targets the full-sliding steady state (global slip), so no
 * incremental Coulomb is needed. Validate p(x) on the strip against eq.(46).
 *
 * Mesh: gmsh -> msh2poly2d (rahmi_422.poly2d). Boundary-id mapping from the
 * converter:  base -> 0,  free -> 1,  contact -> 2.
 */

#include <memory>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/mesh/meshgen.hpp"

#include "hho_params.hpp"
#include "experiment.hpp"

using namespace hho_contact;

using mesh_type = disk::generic_mesh<double, 2>;

struct rahmi_422 : Experiment<mesh_type>
{
    enum : size_t { BASE = 0, FREE = 1, CONTACT = 2 };

    static constexpr double indent = 0.05;   // mm, base pushed up into the punch
    static constexpr double slide  = 0.10;   // mm, tangential drive (sweep for full slip)

    std::unique_ptr<bc_type> bnd, bnd_incr, bnd_C;

    std::string name() const override { return "rahmi_422"; }

    void build_mesh() override
    {
        if (prm.mesh_path.empty())
        {
            std::cerr << "rahmi_422: needs a .poly2d mesh (set 'mesh' in the .dat)\n";
            std::exit(1);
        }
        disk::load_mesh_poly<double>(prm.mesh_path.c_str(), msh);
    }

    bc_wiring setup_bc() override
    {
        auto zero = [](const auto&) {
            Eigen::Matrix<double,2,1> z; z.setZero(); return z; };
        // base drive: press strip into the fixed punch (+y) and slide (+x)
        auto depl_base = [](const auto&) {
            Eigen::Matrix<double,2,1> g;
            g << slide, indent; return g; };

        bnd      = std::make_unique<bc_type>(msh);
        bnd_incr = std::make_unique<bc_type>(msh);
        bnd_C    = std::make_unique<bc_type>(msh);

        bnd->addNeumannBC  (disk::NEUMANN,        FREE,    zero);
        bnd->addContactBC  (disk::SIGNORINI_FACE, CONTACT);
        bnd->addDirichletBC(disk::DIRICHLET,      BASE,    depl_base);

        bnd_incr->addNeumannBC  (disk::NEUMANN,        FREE,    zero);
        bnd_incr->addContactBC  (disk::SIGNORINI_FACE, CONTACT);
        bnd_incr->addDirichletBC(disk::DIRICHLET,      BASE,    zero);

        bnd_C->addNeumannBC  (disk::NEUMANN,        FREE,    zero);
        bnd_C->addContactBC  (disk::SIGNORINI_FACE, CONTACT);
        bnd_C->addDirichletBC(disk::DIRICHLET,      BASE,    zero);

        return { bnd.get(), bnd_incr.get(), bnd_C.get() };
    }

    void setup_terms(ContactTerms<mesh_type>& ct) override
    {
        const double g0 = prm.gamma_0, th = prm.theta;
        ct.add_contact(*bnd_C, g0, g0, th);
        ct.add_coulomb(*bnd_C, g0, g0, th, prm.friction);   // F = 0.3 -> Picard
    }
};

int main(int argc, char** argv)
{
    rahmi_422 exp;
    exp.prm = (argc > 1) ? Params::load(argv[1]) : Params{};

    if (exp.prm.E == 0.0 && exp.prm.mu == 0.0) { exp.prm.E = 100.0; exp.prm.nu = 0.3; }
    if (exp.prm.friction == 0.0) exp.prm.friction = 0.3;

    return exp.run();
}