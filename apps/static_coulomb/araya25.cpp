/*
 * Araya & Chouly [20] -- three-dimensional Bostan-Han with Coulomb friction.
 *
 *   Omega   = (0,8) x (0,4) x (0,4)
 *   Gamma_C = (0,8) x {0}   x (0,4)   unilateral contact + Coulomb (F = 0.8)
 *   Gamma_D = (0,8) x {4}   x (0,4)   Dirichlet u = 0
 *   Gamma_N1= {0}   x (0,4) x (0,4)   Neumann  F = (400, 0, 0)
 *   Gamma_N2= rest of dOmega          homogeneous Neumann
 *
 *   E = 1e3,  nu = 0.3,  f = 0,  theta = -1,  k = 1.   No symmetry edge.
 *
 * Gamma_N2 is defined BY EXCLUSION ("the rest of the boundary"), which the
 * native mesh markers do not express cleanly.  So this driver uses
 * tag-by-plane (mesh_tagging.hpp): paint the whole boundary as free Neumann
 * (marker 2), then override the three named faces.  Order matters.
 *
 *   marker 1 : Gamma_D   (y = 4)   Dirichlet
 *   marker 2 : Gamma_N2  (rest)    free Neumann
 *   marker 3 : Gamma_C   (y = 0)   contact + Coulomb
 *   marker 4 : Gamma_N1  (x = 0)   loaded Neumann
 *
 * Mesh: unit tetrahedral mesh, refined `level` times, scaled to {8,4,4}.
 */

#include <memory>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/mesh/meshgen.hpp"

#include "hho_params.hpp"
#include "experiment.hpp"
#include "mesh_tagging.hpp"

using namespace hho_contact;

using mesh_type = disk::tetrahedral_mesh<double>;

struct Araya25 : Experiment<mesh_type>
{
    enum : size_t { DIRICH_Y4 = 1, FREE_N2 = 2, CONTACT_Y0 = 3, LOAD_X0 = 4 };

    std::unique_ptr<bc_type> bnd;
    std::unique_ptr<bc_type> bnd_incr;
    std::unique_ptr<bc_type> bnd_C;

    std::string name() const override { return "araya25"; }

    void build_mesh() override
    {
        auto mesher = disk::make_simple_mesher(msh);
        for (size_t i = 0; i < prm.mesh_level; i++)
            mesher.refine();
        disk::scale_mesh(msh, {8.0, 4.0, 4.0});

        // Gamma_N2 by exclusion: paint all boundary free, then override.
        disk::tag_boundary(msh, FREE_N2, [](const auto&) { return true; });
        disk::tag_boundary_plane(msh, CONTACT_Y0, 1, 0.0);   // y = 0  contact
        disk::tag_boundary_plane(msh, LOAD_X0,    0, 0.0);   // x = 0  loaded
        disk::tag_boundary_plane(msh, DIRICH_Y4,  1, 4.0);   // y = 4  Dirichlet
    }

    bc_wiring setup_bc() override
    {
        auto zero   = [](const auto&) {
            Eigen::Matrix<double,3,1> z; z.setZero(); return z; };
        auto load_x0 = [](const auto&) {
            Eigen::Matrix<double,3,1> g; g << 400.0, 0.0, 0.0; return g; };

        bnd      = std::make_unique<bc_type>(msh);
        bnd_incr = std::make_unique<bc_type>(msh);
        bnd_C    = std::make_unique<bc_type>(msh);

        
        bnd->addNeumannBC  (disk::NEUMANN,        FREE_N2,    zero);
        bnd->addContactBC  (disk::SIGNORINI_FACE, CONTACT_Y0);
        bnd->addNeumannBC  (disk::NEUMANN,        LOAD_X0,    load_x0);
        bnd->addDirichletBC(disk::DIRICHLET,      DIRICH_Y4,  zero);

        
        bnd_incr->addNeumannBC  (disk::NEUMANN,        FREE_N2,    zero);
        bnd_incr->addContactBC  (disk::SIGNORINI_FACE, CONTACT_Y0);
        bnd_incr->addNeumannBC  (disk::NEUMANN,        LOAD_X0,    zero);
        bnd_incr->addDirichletBC(disk::DIRICHLET,      DIRICH_Y4,  zero);

        
        bnd_C->addNeumannBC  (disk::NEUMANN,        FREE_N2,    zero);
        bnd_C->addContactBC  (disk::SIGNORINI_FACE, CONTACT_Y0);
        bnd_C->addNeumannBC  (disk::NEUMANN,        LOAD_X0,    zero);
        bnd_C->addDirichletBC(disk::DIRICHLET,      DIRICH_Y4,  zero);

        return { bnd.get(), bnd_incr.get(), bnd_C.get() };
    }

    void setup_terms(ContactTerms<mesh_type>& ct) override
    {
        const double g0 = prm.gamma_0, th = prm.theta;
        ct.add_contact(*bnd_C, g0, g0, th);                 // unilateral, 3D
        ct.add_coulomb(*bnd_C, g0, g0, th, prm.friction);   // F = 0.8 -> Picard
    }
};

int main(int argc, char** argv)
{
    Araya25 exp;
    exp.prm = (argc > 1) ? Params::load(argv[1]) : Params{};

    if (exp.prm.E == 0.0 && exp.prm.mu == 0.0) { exp.prm.E = 1.0e3; exp.prm.nu = 0.3; }
    if (exp.prm.friction == 0.0)  exp.prm.friction = 0.8;

    return exp.run();
}