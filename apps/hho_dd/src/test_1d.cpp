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


int main(int argc, char **argv)
{
    using T = double;

    size_t  k_hho = 0;
    size_t  k_dg = 1;
    size_t  n_refs = 1;
    double  Cpen = 10;

    int opt;
    while ((opt = getopt(argc, argv, "e:h:k:r:")) != -1) {
        int         arg_i;
        double      arg_d;
        
        switch (opt) {

        case 'e':
            arg_d = std::stod(optarg);
            if (arg_d > 0.0) Cpen = arg_d;
            break;

        case 'h':
            arg_i = std::stoi(optarg);
            if (arg_i >= 0) k_hho = arg_i;
            break;

        case 'k':    
            arg_i = std::stoi(optarg);
            if (arg_i >= 1) k_dg = arg_i;        
            break;

        case 'r':    
            arg_i = std::stoi(optarg);
            if (arg_i > 0) n_refs = arg_i;        
            break;

        default:
            std::cout << "Invalid option" << std::endl;
            return 1;
        }
    }





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
    for (size_t i = 0; i < n_refs; i++) {
        mesher.refine();
    }

    /*
    msh.transform(
        [](const disk::point<T,2>& pt) {
            auto newx = std::pow(pt.x(), 2);
            auto newy = pt.y()+std::sin(2.0*M_PI*pt.y())*0.1;
            //auto newx = pt.x() + 0.25*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
            //auto newy = pt.y();
            return disk::point<T,2>( newx, newy );
        }
    );
    */

    disk::silo_database db;
    db.create("test_1d.silo");
    db.add_mesh(msh, "srcmesh");
    db.add_mesh(msh, "dstmesh");
    disk::hho_diffusion_solver hho_solver(msh, k_hho);
    auto f = disk::make_rhs_function(msh);
    hho_solver.solve( f );

    dg_diffusion_solver(msh, k_dg, Cpen, db);

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