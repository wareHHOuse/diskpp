/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2025
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <iostream>

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/bases/bases.hpp"
#include "diskpp/methods/hho"
#include "diskpp/methods/implementation_hho/curl.hpp"
#include "diskpp/methods/hho_slapl.hpp"
#include "diskpp/methods/hho_assemblers.hpp"
#include "mumps.hpp"
#include "diskpp/common/timecounter.hpp"
#include "diskpp/output/silo.hpp"

#include "operators.hpp"

int main(void) {

    using T = double;
    //using mesh_type = disk::simplicial_mesh<T,2>;
    using mesh_type = disk::cartesian_mesh<T,2>;
    using point_type = typename mesh_type::point_type;

    mesh_type msh;
    auto mesher = make_simple_mesher(msh);
    mesher.refine();
    mesher.refine();
    mesher.refine();

    auto fun = [](const point_type& pt) {
        disk::static_vector<T,2> ret;
        ret(0) = std::sin(M_PI*pt.x());
        ret(1) = std::cos(M_PI*pt.x());
        return ret;
    };

    auto div_fun = [](const point_type& pt) {
        return M_PI*(std::cos(M_PI*pt.x()) + std::sin(M_PI*pt.y()));
    };

    size_t degree = 1;

    disk::hho_degree_info hdi(degree+1, degree);

    T L2err;
    for (auto& cl : msh) {
        disk::dynamic_vector<T> fun_hho =
            disk::project_function(msh, cl, hdi, fun);
        auto [DR, Adr] = hho_mixedhigh_divrec(msh, cl, degree);

        disk::dynamic_vector<T> div_hho = DR*fun_hho;

        auto cb = disk::make_scalar_monomial_basis(msh, cl, degree+1);
        auto div_proj = disk::project_function(msh, cl, cb, div_fun);

        disk::dynamic_matrix<T> mass = disk::make_mass_matrix(msh, cl, cb);
        disk::dynamic_vector<T> diff = div_proj - div_hho;

        L2err += diff.dot(mass*diff); 
    }

    std::cout << std::sqrt(L2err) << std::endl;

}