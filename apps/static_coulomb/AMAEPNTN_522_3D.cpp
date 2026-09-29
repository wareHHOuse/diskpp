/*
 * AMAEPNTN Test 5.2.2 in 3D  (small-strain short-sliding analogue)
 *
 *   Omega   = (0,4) x (0,1) x (0,4)     flat slab (low H to limit tipping)
 *   Gamma_C = (0,4) x {0}   x (0,4)     unilateral contact + Coulomb
 *   Gamma_T = (0,4) x {1}   x (0,4)     prescribed press-then-slide
 *   Gamma_F = rest of dOmega            free Neumann
 *
 *   depl_top = (u_t, -u_n, 0) : slide in +x, press toward rigid plane at y=0.
 *   Small strain / short sliding: u_t ~ 0.03, u_n ~ 0.01 on this geometry.
 *
 * Gamma_F is "the rest", so tag-by-plane (mesh_tagging.hpp): paint all
 * boundary free, then override contact (y=0) and top (y=H). Top last so it
 * wins on shared edges with the lateral faces.
 *
 *   marker 1 : Gamma_T   (y = 1)   Dirichlet, prescribed
 *   marker 2 : Gamma_F   (rest)    free Neumann
 *   marker 3 : Gamma_C   (y = 0)   contact + Coulomb
 */

#include <memory>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/mesh/meshgen.hpp"

#include "hho_params.hpp"
#include "experiment.hpp"
#include "mesh_tagging.hpp"

using namespace hho_contact;

using mesh_type = disk::tetrahedral_mesh<double>;

struct small_strain_short_sliding_3d : Experiment<mesh_type>
{
    enum : size_t { TOP_YH = 1, FREE_N = 2, CONTACT_Y0 = 3 };

    // press-then-slide amplitudes (small strain: keep well below H)
    static constexpr double u_t =  0.08;   // tangential slide (+x)
    static constexpr double u_n =  0.01;   // normal press (-y, toward y=0)

    std::unique_ptr<bc_type> bnd;
    std::unique_ptr<bc_type> bnd_incr;
    std::unique_ptr<bc_type> bnd_C;

    std::string name() const override { return "small_strain_short_sliding_3d"; }

    void build_mesh() override
    {
        auto mesher = disk::make_simple_mesher(msh);
        for (size_t i = 0; i < prm.mesh_level; i++)
            mesher.refine();
        disk::scale_mesh(msh, {1.0, 1.0, 1.0});   // flat slab, aspect 4:1:4

        // Gamma_F by exclusion: paint all boundary free, then override.
        disk::tag_boundary(msh, FREE_N, [](const auto&) { return true; });
        disk::tag_boundary_plane(msh, CONTACT_Y0, 1, 0.0);   // y = 0  contact
        disk::tag_boundary_plane(msh, TOP_YH,     1, 1.0);   // y = 1  prescribed
    }

    bc_wiring setup_bc() override
    {
        auto zero = [](const auto&) {
            Eigen::Matrix<double,3,1> z; z.setZero(); return z; };
        auto depl_top = [](const auto&) {
            Eigen::Matrix<double,3,1> g;
            g << u_t, -u_n, 0.0; return g; };

        bnd      = std::make_unique<bc_type>(msh);
        bnd_incr = std::make_unique<bc_type>(msh);
        bnd_C    = std::make_unique<bc_type>(msh);

        bnd->addNeumannBC  (disk::NEUMANN,        FREE_N,     zero);
        bnd->addContactBC  (disk::SIGNORINI_FACE, CONTACT_Y0);
        bnd->addDirichletBC(disk::DIRICHLET,      TOP_YH,     depl_top);

        bnd_incr->addNeumannBC  (disk::NEUMANN,        FREE_N,     zero);
        bnd_incr->addContactBC  (disk::SIGNORINI_FACE, CONTACT_Y0);
        bnd_incr->addDirichletBC(disk::DIRICHLET,      TOP_YH,     zero);

        bnd_C->addNeumannBC  (disk::NEUMANN,        FREE_N,     zero);
        bnd_C->addContactBC  (disk::SIGNORINI_FACE, CONTACT_Y0);
        bnd_C->addDirichletBC(disk::DIRICHLET,      TOP_YH,     zero);

        return { bnd.get(), bnd_incr.get(), bnd_C.get() };
    }

    void setup_terms(ContactTerms<mesh_type>& ct) override
    {
        const double g0 = prm.gamma_0, th = prm.theta;
        ct.add_contact(*bnd_C, g0, g0, th);                 // unilateral, 3D
        ct.add_coulomb(*bnd_C, g0, g0, th, prm.friction);   // F -> Picard
    }
};

int main(int argc, char** argv)
{
    small_strain_short_sliding_3d exp;
    exp.prm = (argc > 1) ? Params::load(argv[1]) : Params{};

    if (exp.prm.E == 0.0 && exp.prm.mu == 0.0) { exp.prm.E = 1.0e4; exp.prm.nu = 0.2; }
    if (exp.prm.friction == 0.0) exp.prm.friction = 0.1;

    return exp.run();
}