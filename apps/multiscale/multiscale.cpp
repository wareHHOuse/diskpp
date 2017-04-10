/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016, 2017 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#include <iostream>
#include <fstream>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>
#include <fstream>

#include "common/eigen.hpp"

#include <map>

#include "colormanip.h"

#include "../../config.h"

#include "loaders/loader.hpp"
#include "hho/hho_multiscale.hpp"


#include "contrib/sol2/sol.hpp"


template<typename T>
std::ostream&
operator<<(std::ostream& os, const std::vector<T>& vec)
{
    for (auto& elem : vec)
        os << elem << " ";

    return os;
}


int main2(int argc, char **argv)
{
    using RealType = double;

    typedef Eigen::SparseMatrix<RealType>            sparse_matrix_type;

    if ( !std::regex_match(argv[1], std::regex(".*\\.mesh2d$") ) )
    {
        std::cout << "Mesh format must be Netgen 2D" << std::endl;
        return 1;
    }

    auto msh = disk::load_netgen_2d_mesh<RealType>(argv[1]);

    sol::state lua;
    lua.do_file("params.lua");

    size_t k_outer = lua["k_outer"];
    size_t k_inner = lua["k_inner"];
    size_t rl = lua["rl"];



    std::cout << k_outer << " " << k_inner << " " << rl << std::endl;

    auto basis = make_scaled_monomial_scalar_basis(msh, k_outer-1, k_outer);
    auto ccb = basis.first;
    auto cfb = basis.second;

    disk::multiscale_local_problem<decltype(msh)> mlp(k_inner, rl);

    auto num_cells          = msh.cells_size();
    auto num_faces          = msh.faces_size();
    auto num_cell_dofs      = ccb.size();
    auto num_face_dofs      = cfb.size();
    auto num_boundary_faces = msh.boundary_faces_size();
    auto num_internal_faces = msh.internal_faces_size();

    auto system_size = num_cells*num_cell_dofs + num_internal_faces*num_face_dofs;

    sparse_matrix_type              A(system_size, system_size);
    dynamic_vector<RealType>        b(system_size), x(system_size);

    b = dynamic_vector<RealType>::Zero(system_size);

    typedef Eigen::Triplet<RealType>         triplet_type;
    std::vector<triplet_type>                   triplets;

    std::vector<size_t> face_compress_map( num_faces );
    std::vector<size_t> face_expand_map( num_internal_faces );
    size_t fn = 0, fi = 0;
    for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++, fi++)
    {
        if( msh.is_boundary(*itor) )
            continue;

        face_compress_map.at(fi) = fn;
        face_expand_map.at(fn) = fi;
        fn++;
    }

    std::ofstream mat_ofs("sysmat.dat");
    std::ofstream rhs_ofs("sysrhs.dat");

    size_t elemnum = 0;
    for (auto& cl : msh)
    {
        std::cout << "Assembly: " << elemnum << "/" << msh.cells_size() << std::endl;

        disk::multiscale_local_basis<decltype(msh)> mlb(msh, cl, ccb, cfb, k_inner, rl);
        disk::gradient_reconstruction_multiscale<decltype(msh), decltype(ccb), decltype(cfb)> gradrec(msh, cl, mlb, ccb, cfb);
        
        auto proj = disk::make_projector(msh, cl, ccb, cfb);

        auto f = [](const typename decltype(msh)::point_type& pt) ->
            typename decltype(msh)::point_type::value_type
        {
            //return 2 * M_PI * M_PI * sin(pt.x() * M_PI) * sin(pt.y() * M_PI);
            return sin(pt.x())*sin(pt.y());
        };

        dynamic_vector<RealType> dofs = make_rhs(msh, cl, ccb, cfb, f).head(num_cell_dofs);

        auto fcs = faces(msh, cl);
        std::sort(fcs.begin(), fcs.end());
        std::vector<size_t> l2g(num_cell_dofs + fcs.size() * num_face_dofs);

        for (size_t cell_i = 0; cell_i < num_cell_dofs; cell_i++)
            l2g[cell_i] = elemnum * num_cell_dofs + cell_i;

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];

            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;
            auto face_offset = num_cells * num_cell_dofs + face_compress_map.at(face_id) * num_face_dofs;
            auto pos = face_i * num_face_dofs;

            for (size_t i = 0; i < num_face_dofs; i++)
            {
                if ( msh.is_boundary(fc) )
                    l2g[num_cell_dofs+pos+i] = 0xDEADBEEF;
                else
                    l2g[num_cell_dofs+pos+i] = face_offset+i;
            }
        }

        size_t dsz = gradrec.data.rows();
        for (size_t i = 0; i < dsz; i++)
            for (size_t j = 0; j < dsz; j++)
            {
                if (l2g[i] == 0xDEADBEEF || l2g[j] == 0xDEADBEEF)
                    continue;
                triplets.push_back( triplet_type(l2g[i], l2g[j], gradrec.data(i,j)) );
                mat_ofs << l2g[i]+1 << " " << l2g[j]+1 << " " << gradrec.data(i,j) << std::endl;
            }   
        b.block(elemnum * num_cell_dofs, 0, num_cell_dofs, 1) = dofs;

        elemnum++;

    }

    mat_ofs.close();
    rhs_ofs.close();

    A.setFromTriplets(triplets.begin(), triplets.end());

#ifdef HAVE_INTEL_MKL
    Eigen::PardisoLU<Eigen::SparseMatrix<RealType>>  solver;
#else
    Eigen::SparseLU<Eigen::SparseMatrix<RealType>>   solver;
#endif

    solver.analyzePattern(A);
    solver.factorize(A);
    x = solver.solve(b);

    std::cout << x.transpose() << std::endl;

    std::ofstream sol_ofs("solution.dat");
    std::ofstream sol_ofs_ms("solution_ms.dat");

    //std::cout << " *** POSTPROCESSING ***" << std::endl;
    elemnum = 0;
    for (auto& cl : msh)
    {
        std::cout << cl << std::endl;
        auto fcs = faces(msh, cl);
        std::sort(fcs.begin(), fcs.end());

        dynamic_vector<RealType> dofs =
            dynamic_vector<RealType>::Zero(num_cell_dofs + fcs.size()*num_face_dofs);

        dofs.head(num_cell_dofs) = x.block(elemnum * num_cell_dofs, 0, num_cell_dofs, 1);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            std::cout << " - " << fc << " ";
            if (msh.is_boundary(fc))
            {
                std::cout << "skipped" << std::endl;
                continue;
            }
            std::cout << std::endl;

            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;
            auto face_offset = num_cells * num_cell_dofs + face_compress_map.at(face_id) * num_face_dofs;
            auto pos = face_i * num_face_dofs;

            dofs.block(num_cell_dofs + pos, 0, num_face_dofs, 1) =
                x.block(face_offset, 0, num_face_dofs, 1);
        }

        std::cout << dofs.transpose() << std::endl;

        auto tps = make_test_points(msh, cl, 5);
        for(auto& tp : tps)
        {
            auto phi = ccb.eval_functions(msh, cl, tp);
            auto val = dofs.head(num_cell_dofs).dot(phi);
            
            sol_ofs << tp.x() << " " << tp.y() << " " << val << std::endl;
        }
//#if 0
        disk::multiscale_local_basis<decltype(msh)> mlb(msh, cl, ccb, cfb, k_inner, rl);
        disk::gradient_reconstruction_multiscale<decltype(msh), decltype(ccb), decltype(cfb)> gradrec(msh, cl, mlb, ccb, cfb);
        dynamic_vector<RealType> R = gradrec.oper * dofs;

        auto ms_inner_mesh = mlb.inner_mesh();
        for (auto& icl : ms_inner_mesh)
        {
            auto itps = make_test_points(ms_inner_mesh, icl, 5);
            for(auto& tp : itps)
            {
                auto msphi = mlb.eval_functions(icl, tp);
                auto msval = R.dot(msphi);
            
                sol_ofs_ms << tp.x() << " " << tp.y() << " " << msval << std::endl;
            }
        }
//#endif
        elemnum++;
    }

    sol_ofs.close();
    sol_ofs_ms.close();

}
















int main(int argc, char **argv)
{
    using RealType = double;

    if ( !std::regex_match(argv[1], std::regex(".*\\.mesh2d$") ) )
    {
        std::cout << "Mesh format must be Netgen 2D" << std::endl;
        return 1;
    }

    auto msh = disk::load_netgen_2d_mesh<RealType>(argv[1]);

    sol::state lua;
    lua.do_file("params.lua");

    size_t k_outer = lua["k_outer"];
    size_t k_inner = lua["k_inner"];
    size_t rl = lua["rl"];



    std::cout << k_outer << " " << k_inner << " " << rl << std::endl;

    auto basis = make_scaled_monomial_scalar_basis(msh, k_outer-1, k_outer);
    auto ccb = basis.first;
    auto cfb = basis.second;

    disk::multiscale_local_problem<decltype(msh)> mlp(k_inner, rl);


    std::stringstream ss;
    ss << "testgr_" << k_outer << "_" << k_inner << "_" << rl;

    std::ofstream ofs(ss.str());

    size_t howmany = lua["howmany"].get_or(9999999);
    for (auto& cl : msh)
    {
        //std::cout << "Basis" << std::endl;
        disk::multiscale_local_basis<decltype(msh)> mlb(msh, cl, ccb, cfb, k_inner, rl);
        
        //std::cout << "GR" << std::endl;
        disk::gradient_reconstruction_multiscale<decltype(msh), decltype(ccb), decltype(cfb)> gradrec(msh, cl, mlb, ccb, cfb);

        auto proj = disk::make_projector(msh, cl, ccb, cfb);

        auto f = [&](const typename decltype(msh)::point_type& pt) -> typename decltype(msh)::point_type::value_type
        {
            return lua["process_point"](pt.x(), pt.y());
        };

        auto dofs = proj.project(f);

        //auto bar = barycenter(msh, cl);

        //bar = bar+bar/4;

        //std::cout << dofs.transpose() << std::endl;
        //std::cout << f(bar) << std::endl;
        //std::cout << disk::evaluate_dofs(msh, cl, ccb, dofs, bar) << std::endl;

        dofs(0) = lua["dofs0"].get_or(dofs(0));
        dofs(1) = lua["dofs1"].get_or(dofs(1));
        dofs(2) = lua["dofs2"].get_or(dofs(2));
        dofs(3) = lua["dofs3"].get_or(dofs(3));
        dofs(4) = lua["dofs4"].get_or(dofs(4));
        dofs(5) = lua["dofs5"].get_or(dofs(5));
        dofs(6) = lua["dofs6"].get_or(dofs(6));


        dynamic_vector<double> R = gradrec.oper * dofs;

        //std::cout << cl << std::endl;

        //std::cout << "***********************************" << std::endl;
        //std::cout << dofs.transpose() << std::endl;
        //std::cout << "*****" << std::endl;
        //std::cout << R.transpose() << std::endl;

        //std::cout << "R" << std::endl;
        //std::cout << R.transpose() << std::endl;
        //std::cout << "Postprocess" << std::endl;


        auto inner_mesh = mlb.inner_mesh();

        dynamic_vector<double> z = dynamic_vector<double>::Zero(R.size()-1);
        z(0) = lua["z0"].get_or(z(0));
        z(1) = lua["z1"].get_or(z(1));
        z(2) = lua["z2"].get_or(z(2));
        z(3) = lua["z3"].get_or(z(3));
        z(4) = lua["z4"].get_or(z(4));
        z(5) = lua["z5"].get_or(z(5));
        z(6) = lua["z6"].get_or(z(6));

        R = R.head(R.size()-1);

        double avg = 0;
        for (size_t i = 0; i < R.size(); i++)
            avg += R(i);

        for (auto& icl : inner_mesh)
        {
            auto pts = make_test_points(inner_mesh, icl, 10);
            for (auto& pt : pts)
            {
                auto v = mlb.eval_functions(icl, pt);
                //auto g = mlb.eval_gradients(icl, pt);
                //assert(v.size() == g.rows());
                //ofs << pt.x() << " " << pt.y() << " ";

                //for (size_t zz = 0; zz < v.size(); zz++)
                //    ofs << v(zz) << " " << g(zz,0) << " " << g(zz,1) << " ";

                //ofs << std::endl;

                //std::cout << (gradrec.stiff_mat * z).transpose() << std::endl;

                auto val  = R.dot(v);
                auto val2 = disk::evaluate_dofs(msh, cl, ccb, dofs, pt);
                ofs << pt.x() << " " << pt.y() << " " << val << " " << val2 << std::endl;
            }
        }

        //mlp.assemble(msh, cl, ccb, cfb);

        if ( 0 == --howmany )
            break;
    }

    ofs.close();


    return 0;
}
