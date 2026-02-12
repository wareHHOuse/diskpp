/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2026
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <cstddef>
#include <iostream>
#include <regex>
#include <set>
#include <string>
#include <filesystem>

#include "diskpp/common/util.h"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/loaders/loader_gmsh.hpp"
#include "diskpp/loaders/loader_utils.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/output/silo.hpp"
#include "diskpp/methods/hho_slapl.hpp"
#include "diskpp/methods/hho"
#include "diskpp/methods/dg"
#include "rasdd.hpp"
#include "common.hpp"
#include "solvers.hpp"
#include "diskpp_git_revision.h"


int main(void)
{
    using T = double;

    //disk::generic_mesh<T,1> msh;
    //disk::uniform_mesh_loader<T,1> loader(0,1,100);
    //loader.populate_mesh(msh);

    /*
    disk::generic_mesh<T, 2> msh;
    auto mesher = disk::make_fvca5_hex_mesher(msh);
    mesher.make_level(2);
    */

    disk::cartesian_mesh<T, 2> msh;
    auto mesher = disk::make_simple_mesher(msh);
    mesher.refine();
    mesher.refine();
    mesher.refine();
    mesher.refine();

    msh.transform(
        [](const disk::point<T,2>& pt) {
            auto newx = std::pow(pt.x(), 2);
            auto newy = pt.y()+std::sin(2.0*M_PI*pt.y())*0.1;
            //auto newx = pt.x() + 0.25*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
            //auto newy = pt.y();
            return disk::point<T,2>( newx, newy );
        }
    );
    
    disk::silo_database db;
    db.create("test_1d.silo");
    db.add_mesh(msh, "srcmesh");
    db.add_mesh(msh, "dstmesh");
    hho_diffusion_solver(msh, 4, db);
    dg_diffusion_solver(msh, 5, 10.0, db);

    Eigen::Matrix<T, Eigen::Dynamic, 2> axes_a = Eigen::Matrix<T, Eigen::Dynamic, 2>::Zero( msh.cells_size(), 2 );
    Eigen::Matrix<T, Eigen::Dynamic, 2> axes_b = Eigen::Matrix<T, Eigen::Dynamic, 2>::Zero( msh.cells_size(), 2 );
    for (size_t i = 0; i < msh.cells_size(); i++) {
        const auto& cl = msh.cell_at(i);
        auto axes = disk::scaled_inertia_axes(msh, cl);
        axes_a.row(i) = axes.col(0).transpose();
        axes_b.row(i) = axes.col(1).transpose();
    }

    db.add_variable("srcmesh", "axes_a", axes_a, disk::zonal_variable_t);
    db.add_variable("srcmesh", "axes_b", axes_b, disk::zonal_variable_t);


    return 0;
}