#include "diskpp/quadratures/quadratures.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/bases/bases_new.hpp"
#include "diskpp/bases/bases_operations.hpp"
#include "diskpp/methods/hho_slapl.hpp"



int main(void)
{
    using T = double;

    using mesh_type = disk::simplicial_mesh<T,2>;
    mesh_type msh;
    using point_type = typename mesh_type::point_type;
    auto mesher = disk::make_simple_mesher(msh);
    mesher.refine();
    mesher.refine();

    using namespace disk::basis;
    using namespace disk::hho::slapl;

    degree_info di(1);

    auto f = [](const point_type& pt) {
        return std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
    };

    for (auto& cl : msh)
    {
        auto [R, A] = local_operator(msh, cl, di);
        //std::cout << R << std::endl;
        auto S = local_stabilization(msh, cl, di, R);
        //std::cout << S << std::endl;
    
        disk::dynamic_vector<T> redf = local_reduction(msh, cl, di, f);
        disk::dynamic_vector<T> recf = R*redf;
        std::cout << recf.transpose() << std::endl;

        auto phiR = typename hho_space<mesh_type>::reco_basis_type(msh, cl, di.reco);
        disk::dynamic_vector<T> proj = L2_project(msh, cl, f, phiR);
        std::cout << proj.transpose().tail(phiR.size()-1) << std::endl;
    }

    return 0;
}