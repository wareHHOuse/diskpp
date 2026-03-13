
#define SOL_ALL_SAFETIES_ON 1
#include <sol/sol.hpp>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/loaders/loader_utils.hpp"
#include "diskpp/mesh/meshgen.hpp"

#include "diskpp/methods/hho"

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
    const typename Mesh::cell_type& cl, disk::hho_degree_info di)
{
    const static size_t DIM = Mesh::dimension;
    static_assert(DIM==2 or DIM==3, "Laplacian: only DIM = 2 or DIM = 3");

    using scalar_type = typename Mesh::coordinate_type;
    /* Reconstruction space basis */
    auto rb = disk::make_vector_monomial_basis(msh, cl, di.face_degree());
    auto rbs = rb.size();

    auto cb = disk::make_scalar_monomial_basis(msh, cl, di.cell_degree());
    auto cbs = cb.size();

    auto fcs = faces(msh, cl);
    auto fbs = disk::scalar_basis_size(di.face_degree(), DIM-1);
    auto n_allfacedofs = fcs.size() * fbs;

    /* Mass */
    disk::dynamic_matrix<scalar_type> M =
        disk::dynamic_matrix<scalar_type>::Zero(rbs, rbs);
    
    /* Local problem RHS */
    disk::dynamic_matrix<scalar_type> RHS =
        disk::dynamic_matrix<scalar_type>::Zero(rbs, cbs + n_allfacedofs);

    /* mass */
    auto qps = disk::integrate(msh, cl, di.cell_degree()+di.face_degree());
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

        auto fb = disk::make_scalar_monomial_basis(msh, fc, di.face_degree());
        auto n = normal(msh, cl, fc);
        auto fqps = disk::integrate(msh, fc, 2*di.face_degree()+1);
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

template<typename Mesh>
disk::dynamic_matrix<typename Mesh::coordinate_type>
corr(const Mesh& msh,
    const typename Mesh::cell_type& cl, disk::hho_degree_info di)
{
    const static size_t DIM = Mesh::dimension;
    static_assert(DIM==2 or DIM==3, "Laplacian: only DIM = 2 or DIM = 3");

    using scalar_type = typename Mesh::coordinate_type;
    /* Reconstruction space basis */
    auto rb = disk::make_scalar_harmonic_top_basis(msh, cl, di.reconstruction_degree());
    rb.maximum_polynomial_degree(di.cell_degree()+2);
    auto rbs = rb.size();

    auto cb = disk::make_scalar_monomial_basis(msh, cl, di.cell_degree());
    auto cbs = cb.size();

    auto fcs = faces(msh, cl);
    auto fbs = disk::scalar_basis_size(di.face_degree(), DIM-1);
    auto n_allfacedofs = fcs.size() * fbs;

    using dv = disk::dynamic_vector<scalar_type>;
    using dm = disk::dynamic_matrix<scalar_type>;

    dm oper = dm::Zero(rbs, cbs + n_allfacedofs);

    //auto qps = disk::integrate(msh, cl, rb.degree() + cb.degree());
    //for (auto& qp : qps) {
    //    auto r_dphi = rb.eval_gradients(qp.point());
    //    auto c_dphi = cb.eval_gradients(qp.point());
    //    oper.block(0,0,rbs,cbs) += qp.weight() * r_dphi * c_dphi.transpose();
    //}
    
    /* Now the face stuff */
    for (size_t fcnum = 0; fcnum < fcs.size(); fcnum++) {
        auto fcofs = cbs+fbs*fcnum;
        const auto& fc = fcs[fcnum];

        dm FF = dm::Zero(fbs, fbs);
        dm TF = dm::Zero(fbs, cbs);
        dm FR = dm::Zero(rbs, fbs);

        auto fb = disk::make_scalar_monomial_basis(msh, fc, di.face_degree());
        auto n = normal(msh, cl, fc);
        auto ideg = rb.degree() + std::max(fb.degree(), cb.degree());
        auto fqps = disk::integrate(msh, fc, ideg );
        for (auto& qp : fqps) {
            auto f_phi = fb.eval_functions(qp.point());
            FF += qp.weight() * f_phi * f_phi.transpose();

            auto c_phi = cb.eval_functions(qp.point());
            TF += qp.weight() * f_phi * c_phi.transpose();

            auto r_phi = rb.eval_gradients(qp.point());
            dv r_dphi_n = r_phi*n;
            
            FR += qp.weight() * r_dphi_n * f_phi.transpose();

            oper.block(0, 0, rbs, cbs) += qp.weight() *  r_dphi_n * c_phi.transpose();
            //oper.block(0, fcofs, rbs, fbs) += qp.weight() *  r_dphi_n * f_phi.transpose();
        }

        //oper.block(0, 0, rbs, cbs) -= FR * FF.ldlt().solve(TF);
    }

    std::cout << "**********************\n";
    std::cout << oper/*.leftCols(cbs)*/ << std::endl;
    std::cout << "**********************\n";

    return oper;
}

template<disk::mesh_2D Mesh>
void
adjust_stabfree_recdeg(const Mesh& msh, const typename Mesh::cell_type& cl,
    disk::hho_degree_info& hdi)
{
    size_t cd = hdi.cell_degree();
    size_t fd = hdi.face_degree();
    bool is_mixed_high = (hdi.cell_degree() > hdi.face_degree());
    size_t n = faces(msh, cl).size();   
    size_t rpd = cd+2;

    /* HHO space dofs */
    size_t from = ((cd+2)*(cd+1))/2 + n*(fd+1);
    /* Reconstruction dofs */
    size_t to = ((rpd+2)*(rpd+1))/2;

    if (from <= to) {
        hdi.reconstruction_degree(rpd);
    }
    else {
        /* Every harmonic degree provides 2 additional dofs, therefore
         * we need an increment that it is sufficient to accomodate
         * (from-to) dofs => ((from - to) + (2-1))/2 */
        size_t incr = (from - to + 1)/2;
        hdi.reconstruction_degree(rpd+incr);
    }

    std::cout << "F: " << hdi.face_degree() << ", C: " << hdi.cell_degree();
    std::cout << ", R: " << hdi.reconstruction_degree() << std::endl;
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

    auto tr = lua["transform"];
    if (tr.valid()) {
        msh.transform(tr);
    }

    disk::hho_degree_info hdi(cfg.degree+1, cfg.degree);

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
        auto [R, A] = disk::make_scalar_hho_laplacian(msh, cl, hdi);
        disk::dynamic_vector<T> Iv = project_function(msh, cl, hdi, v);
        if ( lua["kill_interior"].get_or(false) ) {
            auto k = hdi.cell_degree();
            auto cbs = (k+2)*(k+1)/2;
            Iv.segment(0,cbs) = disk::dynamic_vector<T>::Zero(cbs);
        }

        std::cout << "Iv: " << Iv.transpose() << std::endl;

        std::cout << "**********************\n";
        std::cout << R << std::endl;
        std::cout << "**********************\n";

        std::cout << "Gradrec S: " << (R*Iv).transpose() << std::endl;
        
        
        auto G = hho_vlapl(msh, cl, hdi);
        
        std::cout << "Gradrec V: " << (G*Iv).transpose() << std::endl;
    
        adjust_stabfree_recdeg(msh, cl, hdi);

        auto C = corr(msh, cl, hdi);
        std::cout << "Corr V: " << (C*Iv).transpose() << std::endl;
        std::cout << "norm corr:    " << (C*Iv).norm() << std::endl;
        std::cout << "norm gradrec: " << (R*Iv).norm() << std::endl;


        auto oper = disk::make_shl_face_proj_harmonic(msh, cl, hdi);
        std::cout << "**********************\n";
        std::cout << oper.first << std::endl;
        std::cout << "**********************\n";
        std::cout << (oper.first*Iv).transpose() << std::endl;

        auto oper2 = disk::make_shl_harmonic(msh, cl, hdi);
        std::cout << "**********************\n";
        std::cout << oper2.first << std::endl;
        std::cout << "**********************\n";
        std::cout << (oper2.first*Iv).transpose() << std::endl;
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
