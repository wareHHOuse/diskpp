/*
 * Clamped square sliding through a wall 
 *
 *   Omega   = (0,1)^2
 *   Gamma_D = {0} x (0,1)          Dirichlet u = 0   (clamped, left)
 *   Gamma_C = {1} x (0,1)          contact + Coulomb (F = 0.2)  (right, the wall)
 *   Gamma_N = (0,1) x ({0}u{1})    homogeneous Neumann (bottom, top)
 *
 *   lambda = 576923, mu = 384615   (E = 1e6, nu = 0.3)
 *   body force  f = (0, -76518)    (gravity),  no external tractions
 *   Coulomb F = 0.2,  theta = -1
 *   Reference: separation point at y ~ 0.68 on Gamma_C. 
 */

#include <memory>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/mesh/meshgen.hpp"

#include "hho_params.hpp"
#include "experiment.hpp"

using namespace hho_contact;

using mesh_type = disk::generic_mesh<double, 2>;

struct WallSlide : Experiment<mesh_type>
{
    // native tags -- CONFIRM with report_native_tags():
    //   x=0 -> left/Dirichlet, x=1 -> right/contact, y=0 bottom, y=1 top
    enum : size_t { BOTTOM = 0, RIGHT = 1, TOP = 2, LEFT = 3 };

    std::unique_ptr<bc_type> bnd;
    std::unique_ptr<bc_type> bnd_incr;
    std::unique_ptr<bc_type> bnd_C;

    std::string name() const override { return "wall_slide"; }

    void build_mesh() override
    {
        if (!prm.mesh_path.empty())
            disk::load_mesh_poly(prm.mesh_path.c_str(), msh);
        else
        {
            auto mesher = disk::make_fvca5_hex_mesher(msh);
            mesher.make_level(prm.mesh_level);
        }
    }

    // gravity: the whole novelty of this benchmark
    Eigen::Matrix<double,2,1>
    body_force(const typename mesh_type::point_type&) const override
    {
        Eigen::Matrix<double,2,1> f;
        f << 0.0, -76518.0;
        return f;
    }

    bc_wiring setup_bc() override
    {
        auto zero = [](const auto&) {
            Eigen::Matrix<double,2,1> z; z.setZero(); return z; };

        bnd      = std::make_unique<bc_type>(msh);
        bnd_incr = std::make_unique<bc_type>(msh);
        bnd_C    = std::make_unique<bc_type>(msh);

        bnd->addDirichletBC(disk::DIRICHLET,      LEFT,   zero);   // clamped
        bnd->addContactBC  (disk::SIGNORINI_FACE, RIGHT);          // the wall
        bnd->addNeumannBC  (disk::NEUMANN,        BOTTOM, zero);
        bnd->addNeumannBC  (disk::NEUMANN,        TOP,    zero);

        bnd_incr->addDirichletBC(disk::DIRICHLET,      LEFT,   zero);
        bnd_incr->addContactBC  (disk::SIGNORINI_FACE, RIGHT);
        bnd_incr->addNeumannBC  (disk::NEUMANN,        BOTTOM, zero);
        bnd_incr->addNeumannBC  (disk::NEUMANN,        TOP,    zero);

        bnd_C->addDirichletBC(disk::DIRICHLET,      LEFT,   zero);
        bnd_C->addContactBC  (disk::SIGNORINI_FACE, RIGHT);
        bnd_C->addNeumannBC  (disk::NEUMANN,        BOTTOM, zero);
        bnd_C->addNeumannBC  (disk::NEUMANN,        TOP,    zero);

        return { bnd.get(), bnd_incr.get(), bnd_C.get() };
    }

    void setup_terms(ContactTerms<mesh_type>& ct) override
    {
        const double g0 = prm.gamma_0, th = prm.theta;
        ct.add_contact(*bnd_C, g0, g0, th);                 // unilateral
        ct.add_coulomb(*bnd_C, g0, g0, th, prm.friction);   // F = 0.2 -> Picard
    }
};

int main(int argc, char** argv)
{
    WallSlide exp;
    exp.prm = (argc > 1) ? Params::load(argv[1]) : Params{};

    if (exp.prm.E == 0.0 && exp.prm.mu == 0.0) { exp.prm.E = 1.0e6; exp.prm.nu = 0.3; }
    if (exp.prm.friction == 0.0) exp.prm.friction = 0.2;

    return exp.run();
}