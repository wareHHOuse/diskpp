#define SOL_ALL_SAFETIES_ON 1
#include <sol/sol.hpp>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/loaders/loader_utils.hpp"
#include "diskpp/mesh/meshgen.hpp"

#include "diskpp/bases/bases.hpp"
#include "diskpp/methods/hho_slapl.hpp"

static lua_State *gL = nullptr;
struct config {
    size_t      degree = 0;
    size_t      vertices = 4;
};

namespace priv {

sol::table
init(sol::this_state L) {
    sol::state_view lua(L);
    sol::table module = lua.create_table();
    
    module.new_usertype<config>("config",
        "degree", &config::degree,
        "vertices", &config::vertices
    );

    using T = double;
    using point_type = disk::point<T,2>;
    module.new_usertype<point_type>("point",
        sol::constructors<point_type(),
            point_type(const T&, const T&)
            >(),
        "x", sol::resolve<T(void) const>(&point_type::template x<T>),
        "y", sol::resolve<T(void) const>(&point_type::template y<T>)
    );

    return module;
}
};

template<typename Mesh>
auto
hho_vlapl(const Mesh& msh,
    const typename Mesh::cell_type& cl, disk::hho::slapl::degree_info di)
{
    const static size_t DIM = Mesh::dimension;
    static_assert(DIM==2 or DIM==3, "Laplacian: only DIM = 2 or DIM = 3");

    using scalar_type = typename Mesh::coordinate_type;
    /* Reconstruction space basis */
    auto rb = disk::make_vector_monomial_basis(msh, cl, di.face);
    auto rbs = rb.size();

    auto cb = disk::make_scalar_monomial_basis(msh, cl, di.cell);
    auto cbs = cb.size();

    auto fcs = faces(msh, cl);
    auto fbs = disk::scalar_basis_size(di.face, DIM-1);
    auto n_allfacedofs = fcs.size() * fbs;

    /* Mass */
    disk::dynamic_matrix<scalar_type> M =
        disk::dynamic_matrix<scalar_type>::Zero(rbs, rbs);
    
    /* Local problem RHS */
    disk::dynamic_matrix<scalar_type> RHS =
        disk::dynamic_matrix<scalar_type>::Zero(rbs, cbs + n_allfacedofs);

    /* mass */
    auto qps = disk::integrate(msh, cl, di.cell+di.face);
    for (auto& qp : qps) {
        auto rphi = rb.eval_functions(qp.point());
        M += qp.weight() * rphi * rphi.transpose();

        auto cphi = cb.eval_functions(qp.point());
        auto div_rphi = rb.eval_divergences(qp.point());
        RHS.block(0, 0, rbs, cbs) -= qp.weight() * div_rphi * cphi.transpose();
    }
    
    /* Now the face stuff */
    for (size_t fcnum = 0; fcnum < fcs.size(); fcnum++) {
        auto fcofs = cbs+fbs*fcnum;
        const auto& fc = fcs[fcnum];

        auto fb = disk::make_scalar_monomial_basis(msh, fc, di.face);
        auto n = normal(msh, cl, fc);
        auto fqps = disk::integrate(msh, fc, 2*di.face+1);
        for (auto& qp : fqps) {
            auto f_phi = fb.eval_functions(qp.point());
            auto r_phi = rb.eval_functions(qp.point());
            auto r_dphi_n = disk::priv::inner_product(r_phi, (qp.weight()*n).eval());
            RHS.block(0, fcofs, rbs, fbs) +=
                disk::priv::outer_product(r_dphi_n, f_phi);
        }
    }

    Eigen::FullPivLU<disk::dynamic_matrix<scalar_type>> fact(M);

    std::cout << "**********************\n";
    std::cout << RHS << std::endl;
    std::cout << "**********************\n";

    disk::dynamic_matrix<scalar_type> oper = RHS;// fact.solve(RHS);
    return oper;
}

void run(const config& cfg)
{
    auto lua = sol::state_view(gL);

    using T = double;
    using mesh_type = disk::generic_mesh<T,2>;
    using point_type = mesh_type::point_type;

    mesh_type msh;
    double scale = 1.0;
    disk::make_single_element_mesh(msh, scale, cfg.vertices);

    using namespace disk::basis;
    using namespace disk::hho::slapl;

    //degree_info di(cfg.degree+1, cfg.degree);
    degree_info di(cfg.degree);

    auto fun = lua["pfun"];
    if (not fun.valid()) {
        std::cout << "pfun() not defined on the Lua side\n";
        return;
    }

    auto v = [&](const point_type& pt) -> T {
        return fun(pt);
    };

    for (auto& cl : msh)
    {
        auto [R, A] = local_operator(msh, cl, di);
        disk::dynamic_vector<T> Iv = local_reduction(msh, cl, di, v);

        std::cout << "Iv: " << Iv.transpose() << std::endl;

        std::cout << "**********************\n";
        std::cout << R << std::endl;
        std::cout << "**********************\n";

        std::cout << "Gradrec S: " << (R*Iv).transpose() << std::endl;
        
        
        auto G = hho_vlapl(msh, cl, di);
        
        std::cout << "Gradrec V: " << (G*Iv).transpose() << std::endl;
    }
}

extern "C" int luaopen_rkernel(lua_State* L) {
    gL = L;
    sol::state_view lua(L);

    lua["run"] = &run;
    auto r = priv::init(L);
    r.push();

    return 1;
}