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

template<typename Mesh>
void test_gradient_reconstruction(sol::state& lua, const Mesh& msh)
{
    typedef Mesh                            mesh_type;
    typedef typename Mesh::scalar_type      scalar_type;

    size_t k_outer =    lua["k_outer"];
    size_t k_inner =    lua["k_inner"];
    size_t rl =         lua["refinement_levels"];

    std::cout << "K outer: " << k_outer << std::endl;
    std::cout << "K inner: " << k_inner << std::endl;
    std::cout << "Refinement levels: " << rl << std::endl;

    std::cout << k_outer << " " << k_inner << " " << rl << std::endl;

    auto basis = make_scaled_monomial_scalar_basis(msh, k_outer-1, k_outer);
    auto ccb = basis.first;
    auto cfb = basis.second;

    disk::multiscale_local_problem<mesh_type> mlp(k_inner, rl);


    std::string output_filename = lua["testgr_filename"].get_or(std::string("testgr.dat"));

    size_t testgr_eval_num = lua["testgr_eval_num"].get_or(5);

    std::ofstream ofs(output_filename);

    size_t howmany = lua["howmany"].get_or(9999999);
    for (auto& cl : msh)
    {
        disk::multiscale_local_basis<mesh_type> mlb(msh, cl, ccb, cfb, k_inner, rl);
        
        disk::gradient_reconstruction_multiscale<mesh_type, decltype(ccb), decltype(cfb)> gradrec(msh, cl, mlb, ccb, cfb);

        auto proj = disk::make_projector(msh, cl, ccb, cfb);

        auto f = [&](const typename mesh_type::point_type& pt) -> typename mesh_type::point_type::value_type
        {
            return lua["process_point"](pt.x(), pt.y());
        };

        auto dofs = proj.project(f);

        //auto bar = barycenter(msh, cl);

        //bar = bar+bar/4;

        //std::cout << dofs.transpose() << std::endl;
        //std::cout << f(bar) << std::endl;
        //std::cout << disk::evaluate_dofs(msh, cl, ccb, dofs, bar) << std::endl;

        dynamic_vector<double> R = gradrec.oper * dofs;
        R = R.head(R.size()-1);

        //std::cout << cl << std::endl;

        //std::cout << "***********************************" << std::endl;
        //std::cout << dofs.transpose() << std::endl;
        //std::cout << "*****" << std::endl;
        //std::cout << R.transpose() << std::endl;

        //std::cout << "R" << std::endl;
        //std::cout << R.transpose() << std::endl;
        //std::cout << "Postprocess" << std::endl;


        auto inner_mesh = mlb.inner_mesh();

        for (auto& icl : inner_mesh)
        {
            
            auto pts = make_test_points(inner_mesh, icl, testgr_eval_num);
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

                R = dynamic_vector<scalar_type>::Zero(R.size());

                R(0) = 0;

                R(1) = 0;
                R(2) = 0;
                
                R(3) = 1;
                R(4) = 0;
                
                R(5) = 0;
                R(6) = 0;

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
}


template<typename T>
auto
do_static_condensation(const dynamic_matrix<T>& mat,
                       const dynamic_vector<T>& vec,
                       size_t split)
{
    auto cell_size = split;
    assert(mat.rows() == mat.cols());
    auto face_size = mat.rows() - split;

    dynamic_matrix<T> K_TT = mat.topLeftCorner(cell_size, cell_size);
    dynamic_matrix<T> K_TF = mat.topRightCorner(cell_size, face_size);
    dynamic_matrix<T> K_FT = mat.bottomLeftCorner(face_size, cell_size);
    dynamic_matrix<T> K_FF = mat.bottomRightCorner(face_size, face_size);

    assert(K_TT.cols() == cell_size);
    assert(K_TT.cols() + K_TF.cols() == mat.cols());
    assert(K_TT.rows() + K_FT.rows() == mat.rows());
    assert(K_TF.rows() + K_FF.rows() == mat.rows());
    assert(K_FT.cols() + K_FF.cols() == mat.cols());

    auto K_TT_ldlt = K_TT.llt();
    dynamic_matrix<T> AL = K_TT_ldlt.solve(K_TF);
    dynamic_vector<T> bL = K_TT_ldlt.solve(vec);

    dynamic_matrix<T> AC = K_FF - K_FT * AL;
    dynamic_vector<T> bC = /* no projection on faces, eqn. 26*/ - K_FT * bL;

    return std::make_pair(AC, bC);
}

template<typename T>
dynamic_vector<T>
recover_static_condensation(const dynamic_matrix<T>& mat,
                            const dynamic_vector<T>& vec,
                            const dynamic_vector<T>& solF,
                            size_t split)
{
    auto cell_size = split;
    assert(mat.rows() == mat.cols());
    auto face_size = mat.rows() - split;

    dynamic_vector<T> ret( mat.rows() );

    dynamic_matrix<T> K_TT = mat.topLeftCorner(cell_size, cell_size);
    dynamic_matrix<T> K_TF = mat.topRightCorner(cell_size, face_size);

    dynamic_vector<T> solT = K_TT.llt().solve(vec - K_TF*solF);

    ret.head(cell_size) = solT;
    ret.tail(face_size) = solF;

    return ret;
}


template<typename Mesh>
void test_full_problem(sol::state& lua, const Mesh& msh)
{
    typedef Mesh                        mesh_type;
    typedef typename Mesh::scalar_type  scalar_type;
    typedef Eigen::SparseMatrix<scalar_type> sparse_matrix_type;

    size_t k_outer = lua["k_outer"];
    size_t k_inner = lua["k_inner"];
    size_t rl = lua["refinement_levels"];


    auto basis = make_scaled_monomial_scalar_basis(msh, k_outer-1, k_outer);
    auto ccb = basis.first;
    auto cfb = basis.second;

    disk::multiscale_local_problem<mesh_type> mlp(k_inner, rl);

    auto num_cells          = msh.cells_size();
    auto num_faces          = msh.faces_size();
    auto num_cell_dofs      = ccb.size();
    auto num_face_dofs      = cfb.size();
    auto num_boundary_faces = msh.boundary_faces_size();
    auto num_internal_faces = msh.internal_faces_size();

    auto system_size = /*num_cells*num_cell_dofs +*/ num_internal_faces*num_face_dofs;

    sparse_matrix_type              A(system_size, system_size);
    dynamic_vector<scalar_type>     b(system_size), x(system_size);

    b = dynamic_vector<scalar_type>::Zero(system_size);

    typedef Eigen::Triplet<scalar_type>         triplet_type;
    std::vector<triplet_type>                   triplets;

    std::vector<size_t> face_compress_map( num_faces );
    std::vector<size_t> face_expand_map( num_internal_faces );
    size_t fn = 0, fi = 0;
    for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++, fi++)
    {
        std::cout << *itor << std::endl;
        if( msh.is_boundary(*itor) )
            continue;

        face_compress_map.at(fi) = fn;
        face_expand_map.at(fn) = fi;
        fn++;
    }

    std::cout << face_compress_map << std::endl;
    std::cout << face_expand_map << std::endl;

    //std::ofstream mat_ofs("sysmat.dat");
    //std::ofstream rhs_ofs("sysrhs.dat");

    auto f = [](const typename mesh_type::point_type& pt) ->
                    typename mesh_type::point_type::value_type
    {
        //return 1.;
        //return 2 * M_PI * M_PI * sin(pt.x() * M_PI) * sin(pt.y() * M_PI);
        return sin(pt.x())*sin(pt.y());
    };

    size_t elemnum = 0;
    for (auto& cl : msh)
    {
        std::cout << "Assembly: " << elemnum << "/" << msh.cells_size() << std::endl;

        disk::multiscale_local_basis<mesh_type> mlb(msh, cl, ccb, cfb, k_inner, rl);
        disk::gradient_reconstruction_multiscale<mesh_type, decltype(ccb), decltype(cfb)> gradrec(msh, cl, mlb, ccb, cfb);
        
        auto proj = disk::make_projector(msh, cl, ccb, cfb);

        dynamic_vector<scalar_type> dofs = make_rhs(msh, cl, ccb, cfb, f);
        dynamic_vector<scalar_type> cdofs = dofs.head(num_cell_dofs);
        
        auto fcs = faces(msh, cl);
        //std::sort(fcs.begin(), fcs.end());
        std::vector<size_t> l2g(/*num_cell_dofs +*/ fcs.size() * num_face_dofs);

        //for (size_t cell_i = 0; cell_i < num_cell_dofs; cell_i++)
        //    l2g[cell_i] = elemnum * num_cell_dofs + cell_i;

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];

            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;
            auto face_offset = /*num_cells * num_cell_dofs +*/ face_compress_map.at(face_id) * num_face_dofs;
            auto pos = face_i * num_face_dofs;

            for (size_t i = 0; i < num_face_dofs; i++)
            {
                if ( msh.is_boundary(fc) )
                    l2g.at(/*num_cell_dofs+*/pos+i) = 0xDEADBEEF;
                else
                    l2g.at(/*num_cell_dofs+*/pos+i) = face_offset+i;
            }
        }

        auto sc = do_static_condensation(gradrec.data, cdofs, num_cell_dofs);

        //size_t dsz = gradrec.data.rows();
        size_t dsz = sc.first.rows();
        for (size_t i = 0; i < dsz; i++)
        {
            if (l2g[i] == 0xDEADBEEF)
                continue;

            for (size_t j = 0; j < dsz; j++)
            {
                if (l2g[j] == 0xDEADBEEF)
                    continue;

                triplets.push_back( triplet_type(l2g[i], l2g[j], sc.first(i,j)) );
                //mat_ofs << l2g[i]+1 << " " << l2g[j]+1 << " " << gradrec.data(i,j) << std::endl;
            }

            b(l2g.at(i)) += sc.second(i);
        }
        //b.block(elemnum * num_cell_dofs, 0, num_cell_dofs, 1) = dofs;

        elemnum++;

    }

    //mat_ofs.close();
    //rhs_ofs.close();

    A.setFromTriplets(triplets.begin(), triplets.end());

#ifdef HAVE_INTEL_MKL
    Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
#else
    Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
#endif

    solver.analyzePattern(A);
    solver.factorize(A);
    x = solver.solve(b);

    std::cout << b.transpose() << std::endl;
    std::cout << x.transpose() << std::endl;

    std::ofstream sol_ofs("solution.dat");
    std::ofstream sol_ofs_ms("solution_ms.dat");

    //std::cout << " *** POSTPROCESSING ***" << std::endl;
    elemnum = 0;
    for (auto& cl : msh)
    {
        std::cout << cl << std::endl;
        auto fcs = faces(msh, cl);
        //std::sort(fcs.begin(), fcs.end());

        dynamic_vector<scalar_type> rhs = make_rhs(msh, cl, ccb, cfb, f);

        dynamic_vector<scalar_type> dofs =
            dynamic_vector<scalar_type>::Zero(/*num_cell_dofs +*/ fcs.size()*num_face_dofs);

        //dofs.head(num_cell_dofs) = x.block(elemnum * num_cell_dofs, 0, num_cell_dofs, 1);

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
            auto face_offset = /*num_cells * num_cell_dofs +*/ face_compress_map.at(face_id) * num_face_dofs;
            auto pos = face_i * num_face_dofs;

            dofs.block(/*num_cell_dofs + */pos, 0, num_face_dofs, 1) =
                x.block(face_offset, 0, num_face_dofs, 1);
        }

        //auto tps = make_test_points(msh, cl, 5);
        //for(auto& tp : tps)
        //{
        //    auto phi = ccb.eval_functions(msh, cl, tp);
        //    auto val = dofs.head(num_cell_dofs).dot(phi);
        //    
        //    sol_ofs << tp.x() << " " << tp.y() << " " << val << std::endl;
        //}
//#if 0
        disk::multiscale_local_basis<mesh_type> mlb(msh, cl, ccb, cfb, k_inner, rl);
        disk::gradient_reconstruction_multiscale<mesh_type, decltype(ccb), decltype(cfb)> gradrec(msh, cl, mlb, ccb, cfb);
        
        dynamic_vector<scalar_type> rhs_c = rhs.head(num_cell_dofs);
        dynamic_vector<scalar_type> local_dofs = recover_static_condensation(gradrec.data, rhs_c, dofs, num_cell_dofs);

        dynamic_vector<scalar_type> R = gradrec.oper * local_dofs;
        R = R.head(R.size()-1);

        std::cout << "Elem dofs" << std::endl;
        std::cout << R.transpose() << std::endl;

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

int
main(int argc, char **argv)
{
    if ( !std::regex_match(argv[1], std::regex(".*\\.mesh2d$") ) )
    {
        std::cout << "Mesh format must be Netgen 2D" << std::endl;
        return 1;
    }

    auto msh = disk::load_netgen_2d_mesh<double>(argv[1]);

    std::cout << "AVG H: " << average_diameter(msh) << std::endl;

    sol::state lua;
    lua.do_file("params.lua");

    std::string mode = lua["mode"].get_or( std::string("testgr") );
    std::cout << "Mode: " << mode << std::endl;

    if (mode == "testgr")
        test_gradient_reconstruction(lua, msh);
    else if (mode == "testfull")
        test_full_problem(lua, msh);

    return 0;
}


