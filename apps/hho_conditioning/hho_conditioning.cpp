/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
 * Nicolas Pignet  (C) 2018                     nicolas.pignet@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */

#include <iostream>
#include <sstream>
#include <iomanip>
#include <regex>

#include <Eigen/Eigenvalues>

#include "loaders/loader.hpp"
#include "methods/hho"
#include "output/silo.hpp"
#include "solvers/solver.hpp"

#include "contrib/sol2/sol.hpp"
#include "contrib/timecounter.h"

template<typename Mesh>
bool
estimate_element_cond(sol::state& lua, const Mesh& msh)
{
    typedef Mesh                                       mesh_type;
    typedef typename mesh_type::coordinate_type            scalar_type;
    typedef typename mesh_type::cell                   cell_type;
    typedef typename mesh_type::face                   face_type;

    size_t  degree  = lua["config"]["degree"].get_or(1);

    disk::hho_degree_info hdi(degree);

    auto lua_rhs_fun = lua["right_hand_side"];
    if ( !lua_rhs_fun.valid() )
    {
        std::cout << "[config] right_hand_side() not defined" << std::endl;
        return false;
    }

    auto rhs_fun = [&](const typename mesh_type::point_type& pt) -> scalar_type { return lua_rhs_fun(pt.x(), pt.y()); };

    for (auto& cl : msh)
    {
        const auto cb   = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        const auto gr   = make_scalar_hho_laplacian(msh, cl, hdi);
        const auto stab = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);
        const auto rhs  = make_rhs(msh, cl, cb, rhs_fun);
        const auto A    = gr.second + stab;

        auto cbs = cb.size();
        auto cloc = A.block(0,0,cbs,cbs);

        Eigen::JacobiSVD<dynamic_matrix<scalar_type>> svd(cloc);
        auto sigma_max = svd.singularValues()(0);
        auto sigma_min = svd.singularValues()(svd.singularValues().size()-1);
        auto cond =  sigma_max / sigma_min;
        std::cout << "Condition number: " << cond << std::endl;
        //std::cout << "Sigma max: " << sigma_max << std::endl;
        //std::cout << "Sigma min: " << sigma_min << std::endl;
    }

    return true;
}

int main(int argc, char **argv)
{
    using scalar_type = double;

    if (argc != 2)
    {
        std::cout << "Please specify configuration file" << std::endl;
        return 1;
    }

    sol::state lua;
    lua.open_libraries(sol::lib::math);
    lua["config"] = lua.create_table();

    auto r = lua.do_file(argv[1]);

    if ( !r.valid() )
    {
        std::cout << "Problems opening configuration file" << std::endl;
        return 1;
    }

    std::string     input_mesh  = lua["config"]["input_mesh"];


    if (std::regex_match(input_mesh, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;

        typedef disk::generic_mesh<scalar_type, 2>  mesh_type;

        mesh_type msh;
        disk::fvca5_mesh_loader<scalar_type, 2> loader;

        if (!loader.read_mesh(input_mesh))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }

        loader.populate_mesh(msh);

        estimate_element_cond(lua, msh);
    }
}
