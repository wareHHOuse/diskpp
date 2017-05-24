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

//#undef NDEBUG

#include "common/eigen.hpp"

#include <map>

#include "colormanip.h"

#include "../../config.h"

#include "loaders/loader.hpp"
#include "hho/hho_multiscale.hpp"


#include "contrib/sol2/sol.hpp"

#include "mesh/mesh_hierarchy.hpp"

template<typename T>
std::ostream&
operator<<(std::ostream& os, const std::vector<T>& vec)
{
    for (auto& elem : vec)
        os << elem << " ";

    return os;
}



template<typename Mesh>
void plot_multiscale_basis_functions(sol::state& lua, const Mesh& msh)
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


    size_t pms_eval_num = lua["pms_eval_num"].get_or(5);

    

    size_t elemnum = 0;
    for (auto& cl : msh)
    {
        std::stringstream ss;
        ss << "./mb/multiscale_basis_C" << elemnum << ".dat";
        std::ofstream ofs(ss.str());
        if (!ofs.is_open())
            std::cout << "Problem opening " << ss.str() << std::endl;

        disk::multiscale_local_basis<mesh_type> mlb(msh, cl, ccb, cfb, k_inner, rl);

        auto inner_mesh = mlb.inner_mesh();

        for (auto& icl : inner_mesh)
        {
            auto tps = make_test_points(inner_mesh, icl, pms_eval_num);
            for (auto& tp : tps)
            {
                auto val = mlb.eval_functions(icl, tp);
                ofs << tp.x() << " " << tp.y() << " ";
                for (size_t i = 0; i < val.size(); i++)
                    ofs << val(i) << " ";
                ofs << std::endl;
            }
        }

        ofs.close();
        elemnum++;
    }
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
bool
conjugated_gradient(const Eigen::SparseMatrix<T>& A,
                    const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
                    Eigen::Matrix<T, Eigen::Dynamic, 1>& x)
{
    size_t                      N = A.cols();
    size_t                      iter = 0;
    T                           nr, nr0;
    T                           alpha, beta, rho;
    
    Eigen::Matrix<T, Eigen::Dynamic, 1> d(N), r(N), r0(N), y(N);
    x = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(N);
 

    r0 = d = r = b - A*x;
    nr = nr0 = r.norm();
    
    std::ofstream ofs("cocg_nopre_convergence.txt");
    
    while ( nr/nr0 > 1e-8 && iter < 4000 && nr/nr0 < 10000 )
    {
        std::cout << "                                                 \r";
        std::cout << " -> Iteration " << iter << ", rr = ";
        std::cout << nr/nr0 << "\b\r";
        std::cout.flush();
        
        ofs << nr/nr0 << std::endl;
        y = A*d;
        rho = r.dot(r);
        alpha = rho/d.dot(y);
        x = x + alpha * d;
        r = r - alpha * y;
        beta = r.dot(r)/rho;
        d = r + beta * d;
        
        nr = r.norm();
        iter++;
    }
    
    ofs << nr/nr0 << std::endl;
    ofs.close();
    
    std::cout << " -> Iteration " << iter << ", rr = " << nr/nr0 << std::endl;
    
    return true;
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











#define WITH_STATIC_CONDENSATION
#define PLOT_SOLUTION




template<typename Mesh>
void 
test_full_problem_error(sol::state& lua, const Mesh& msh, size_t rl,
                        disk::mesh_hierarchy<typename Mesh::scalar_type>& hierarchy,
                        const std::pair<size_t, dynamic_vector<typename Mesh::scalar_type>>& monoscale_solution)
{
    typedef Mesh                        mesh_type;
    typedef typename Mesh::scalar_type  scalar_type;
    typedef Eigen::SparseMatrix<scalar_type> sparse_matrix_type;

    size_t k_outer = lua["k_outer"];
    size_t k_inner = lua["k_inner"];
    size_t testpoints = lua["test_points"].get_or(3);


    auto basis = make_scaled_monomial_scalar_basis(msh, k_outer-1, k_outer);
    auto ccb = basis.first;
    auto cfb = basis.second;

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
        if( msh.is_boundary(*itor) )
            continue;

        face_compress_map.at(fi) = fn;
        face_expand_map.at(fn) = fi;
        fn++;
    }

    auto f = [](const typename mesh_type::point_type& pt) ->
                    typename mesh_type::point_type::value_type
    {
        //return 2.0 * M_PI * M_PI * sin(pt.x() * M_PI) * sin(pt.y() * M_PI);
        return sin(pt.x())*sin(pt.y());
    };

    /*
    auto f = [](const typename mesh_type::point_type& p) ->
                    typename mesh_type::point_type::value_type {
        auto eps = 0.1;
        return +2.0 * M_PI * M_PI * ( 100.0*sin(M_PI*p.y()/eps)*sin(M_PI*p.y()/eps)*cos(M_PI*p.x()/eps)*cos(M_PI*p.x()/eps) + 1.0) * sin(M_PI*p.x()) * sin(M_PI*p.y())
               - 200 * M_PI * M_PI * sin(M_PI*p.x()) * sin(M_PI*p.y()/eps) * cos(M_PI*p.y()) * cos(M_PI*p.x()/eps) * cos(M_PI*p.x()/eps) * cos(M_PI*p.y()/eps) / eps
               + 200 * M_PI * M_PI * sin(M_PI*p.y()) * sin(M_PI*p.x()/eps) * cos(M_PI*p.x()) * sin(M_PI*p.y()/eps) * sin(M_PI*p.y()/eps) * cos(M_PI*p.x()/eps) / eps;
    
    };
    */
    /*
    auto f = [](const typename mesh_type::point_type& p) ->
                    typename mesh_type::point_type::value_type {
        return
        + 4.0 * M_PI * M_PI * (100*sin(M_PI*p.y())*sin(M_PI*p.y())*cos(M_PI*p.x())*cos(M_PI*p.x()) + 1)*sin(M_PI*p.x())*sin(M_PI*p.y())
        + 400 * M_PI * M_PI * sin(M_PI*p.x()) * pow(sin(M_PI*p.y()),3) * cos(M_PI*p.x()) * cos(M_PI*p.x()) 
        - 400 * M_PI * M_PI * sin(M_PI*p.x()) * sin(M_PI*p.y()) * cos(M_PI*p.x()) * cos(M_PI*p.x()) * cos(M_PI*p.y()) * cos(M_PI*p.y());
    };
    */

    /*
    auto f = [](const typename mesh_type::point_type& p) ->
                    typename mesh_type::point_type::value_type {

        return 
        - 2*(100*pow(sin(10.0*M_PI*p.y()),2)*pow(cos(10.0*M_PI*p.x()),2) + 1)*(-M_PI*M_PI*sin(M_PI*p.x())*sin(M_PI*p.y())
        + 10.0*M_PI*M_PI*sin(10.0*M_PI*p.x())*sin(10.0*M_PI*p.y())) 
        - 2000.0*M_PI*(M_PI*sin(M_PI*p.x())*cos(M_PI*p.y()) 
        - M_PI*sin(10.0*M_PI*p.x())*cos(10.0*M_PI*p.y()))*sin(10.0*M_PI*p.y())*pow(cos(10.0*M_PI*p.x()),2)*cos(10.0*M_PI*p.y()) 
        + 2000.0*M_PI*(M_PI*sin(M_PI*p.y())*cos(M_PI*p.x()) 
        - M_PI*sin(10.0*M_PI*p.y())*cos(10.0*M_PI*p.x()))*sin(10.0*M_PI*p.x())*pow(sin(10.0*M_PI*p.y()),2)*cos(10.0*M_PI*p.x());
    };
    */

    std::vector<dynamic_matrix<scalar_type>> gr_opers, gr_datas;
    gr_opers.reserve( msh.cells_size() );
    gr_datas.reserve( msh.cells_size() );

    std::vector<disk::multiscale_local_basis<mesh_type>> multiscale_bases;
    multiscale_bases.reserve( msh.cells_size() );

    size_t elemnum = 0;
    for (auto& cl : msh)
    {
        std::cout << "Assembly: " << elemnum << "/" << msh.cells_size() << ": ";
        std::cout.flush();
        disk::multiscale_local_basis<mesh_type> mlb(msh, cl, ccb, cfb, k_inner, rl);
        multiscale_bases.push_back(mlb);
        std::cout << "B ";
        std::cout.flush();

        disk::gradient_reconstruction_multiscale<mesh_type, decltype(ccb), decltype(cfb)> gradrec(msh, cl, mlb, ccb, cfb);
        gr_opers.push_back( gradrec.oper );
        gr_datas.push_back( gradrec.data );
        std::cout << "GR ";
        std::cout.flush();

        auto proj = disk::make_projector(msh, cl, ccb, cfb);

        dynamic_vector<scalar_type> dofs = make_rhs(msh, cl, ccb, cfb, f);
        dynamic_vector<scalar_type> cdofs = dofs.head(num_cell_dofs);
        
        auto fcs = faces(msh, cl);
        std::vector<size_t> l2g(/*num_cell_dofs +*/ fcs.size() * num_face_dofs);

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
            }

            b(l2g.at(i)) += sc.second(i);
        }

        std::cout << "ASM\r";
        std::cout.flush();

        elemnum++;
    }

    std::cout << std::endl;


    A.setFromTriplets(triplets.begin(), triplets.end());

#ifdef HAVE_INTEL_MKL
    Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
#else
    Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
#endif

    solver.analyzePattern(A);
    solver.factorize(A);
    x = solver.solve(b);

    //std::cout << " *** POSTPROCESSING ***" << std::endl;

    typedef disk::quadrature<mesh_type, typename mesh_type::cell>      cell_quadrature_type;

    auto monoscale_degree   = monoscale_solution.first;
    auto monoscale_all_dofs = monoscale_solution.second;
    //std::cout << "MM: " << monoscale_all_dofs.size() << std::endl;
    auto monoscale_mesh = *std::next(hierarchy.meshes_begin(), hierarchy.meshes_size()-1);
    typedef disk::scaled_monomial_scalar_basis<mesh_type, typename mesh_type::cell> 
               monoscale_cell_basis_type;
    
    monoscale_cell_basis_type     monoscale_cell_basis(monoscale_degree);
    cell_quadrature_type          cell_quadrature(2*monoscale_degree+6); //XXX
    size_t monoscale_cbs = monoscale_cell_basis.size();
    
    scalar_type error = 0.0, l2_err = 0.0;

    std::ofstream ofs("comparison.dat");
    std::ofstream probl("probl.dat");

    elemnum = 0;
    for (auto& cl : msh)
    {
        std::cout << "Postprocess: " << elemnum << "/" << msh.cells_size() << "\r";
        std::cout.flush();
        auto fcs = faces(msh, cl);

        dynamic_vector<scalar_type> rhs = make_rhs(msh, cl, ccb, cfb, f);

        dynamic_vector<scalar_type> dofs =
            dynamic_vector<scalar_type>::Zero(/*num_cell_dofs +*/ fcs.size()*num_face_dofs);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            if (msh.is_boundary(fc))
                continue;

            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;
            auto face_offset = /*num_cells * num_cell_dofs +*/ face_compress_map.at(face_id) * num_face_dofs;
            auto pos = face_i * num_face_dofs;

            dofs.block(/*num_cell_dofs + */pos, 0, num_face_dofs, 1) =
                x.block(face_offset, 0, num_face_dofs, 1);
        }

        //disk::multiscale_local_basis<mesh_type> mlb(msh, cl, ccb, cfb, k_inner, rl);
        auto mlb = multiscale_bases.at(elemnum);
        dynamic_vector<scalar_type> rhs_c = rhs.head(num_cell_dofs);
        dynamic_vector<scalar_type> local_dofs = recover_static_condensation(gr_datas.at(elemnum), rhs_c, dofs, num_cell_dofs);

        dynamic_vector<scalar_type> R = gr_opers.at(elemnum) * local_dofs;
        R = R.head(R.size()-1);

        auto ms_inner_mesh = mlb.inner_mesh();
        for (auto& icl : ms_inner_mesh)
        {   
            scalar_type local_err = 0.0;
            scalar_type local_l2_err = 0.0;
            auto qps = cell_quadrature.integrate(ms_inner_mesh, icl);
            
            for(auto& qp : qps)
            {
                dynamic_vector<scalar_type> grad_multi = dynamic_vector<scalar_type>::Zero(2);
                dynamic_vector<scalar_type> grad_mono = dynamic_vector<scalar_type>::Zero(2);

                auto tp = qp.point();
                auto local_tensor = disk::make_material_tensor(tp);
                local_tensor(0,0) = sqrt(local_tensor(0,0));
                local_tensor(1,1) = sqrt(local_tensor(1,1));
                auto msphi = mlb.eval_functions(icl, tp);
                auto msval = R.dot(msphi);

                dynamic_matrix<scalar_type> msdphi = mlb.eval_gradients(icl, tp);
                for (size_t i = 0; i < R.size(); i++)
                    grad_multi += (msdphi.block(i, 0, 1, 2) * R(i)).transpose();

                try {
                    auto mono_cell_id = hierarchy.locate_point(tp, hierarchy.meshes_size()-1);
                    auto mono_cell = *std::next(monoscale_mesh.cells_begin(), mono_cell_id);
                    auto mono_dofs_offset = monoscale_cbs * mono_cell_id;
                    dynamic_vector<scalar_type> mono_dofs =
                        monoscale_all_dofs.block(mono_dofs_offset, 0, monoscale_cbs, 1);

                    auto mono_phi = monoscale_cell_basis.eval_functions(monoscale_mesh, mono_cell, tp);
                    auto mono_val = mono_dofs.dot(mono_phi);

                    dynamic_matrix<scalar_type> mono_dphi = monoscale_cell_basis.eval_gradients(monoscale_mesh, mono_cell, tp);
                    for (size_t i = 0; i < mono_dofs.size(); i++)
                        grad_mono += (mono_dphi.block(i, 0, 1, 2) * mono_dofs(i)).transpose();
                    
                    //auto eps = 0.1;
                    //grad_mono(0) = M_PI * cos(M_PI * tp.x()) * sin(M_PI * tp.y()) + M_PI * sin(M_PI*tp.y()/eps)*cos(M_PI*tp.x()/eps);
                    //grad_mono(1) = M_PI * sin(M_PI * tp.x()) * cos(M_PI * tp.y()) + M_PI * sin(M_PI*tp.x()/eps)*cos(M_PI*tp.y()/eps);

                    ofs << tp.x() << " " << tp.y() << " " << msval << " " << mono_val << std::endl;
                    local_l2_err += qp.weight() * (msval-mono_val) * (msval-mono_val);
                }

                catch(...)
                {
                    std::cout << "problem with: " << tp << std::endl;
                    probl << tp.x() << " " << tp.y() << " " << 1 << std::endl;
                }

                dynamic_vector<scalar_type> diff = local_tensor*(grad_multi-grad_mono);

                local_err += qp.weight() * diff.dot(diff);
            }
            error += local_err;
            l2_err += local_l2_err;
        }

        elemnum++;

        //std::cout << "Gradients: " << std::endl;
        //std::cout << grad_multi << std::endl;
        //std::cout << grad_mono << std::endl;
    }

    std::cout << std::endl;
    ofs.close();
    probl.close();

    std::cout << "Error: " << sqrt(error) << std::endl;
    std::cout << "Error: " << sqrt(l2_err) << std::endl;
}






template<typename mesh_type>
std::pair<size_t, dynamic_vector<typename mesh_type::scalar_type>>
load_monoscale_solution(const mesh_type& msh, const std::string& sol_fn)
{
    typedef typename mesh_type::cell                cell_type;
    typedef typename mesh_type::scalar_type         scalar_type;

    std::ifstream ifs(sol_fn, std::ifstream::binary);
    
    size_t sol_num_elements, cell_basis_deg, face_basis_deg;
    
    ifs.read(reinterpret_cast<char *>(&sol_num_elements), sizeof(size_t));
    ifs.read(reinterpret_cast<char *>(&cell_basis_deg), sizeof(size_t));
    ifs.read(reinterpret_cast<char *>(&face_basis_deg), sizeof(size_t));
    
    if (sol_num_elements != msh.cells_size())
    {
        std::cout << "Solution has a different number of elements than the mesh (";
        std::cout << sol_num_elements << " vs. " << msh.cells_size() << ")" << std::endl;
        throw std::invalid_argument("");
    }
        
    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    
    cell_basis_type     cell_basis_k1(cell_basis_deg+1);
    
    size_t local_dofs_size = cell_basis_k1.size();
    size_t dofs_vec_size = msh.cells_size() * local_dofs_size;
    dynamic_vector<scalar_type> cell_dofs = dynamic_vector<scalar_type>::Zero(dofs_vec_size);
    
    for (size_t i = 0; i < dofs_vec_size; i++)
        ifs.read(reinterpret_cast<char *>(&cell_dofs(i)), sizeof(scalar_type));
    
    ifs.close();

    return std::make_pair(cell_basis_deg+1,cell_dofs);
}


void run_multiscale_tests(sol::state& lua)
{
    using T = double;

    typedef disk::simplicial_mesh<T,2>          mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;
    typedef Eigen::SparseMatrix<scalar_type>    sparse_matrix_type;

    size_t k_outer = lua["k_outer"].get_or(1);
    size_t k_inner = lua["k_inner"].get_or(1);
    size_t hierarchy_levels = lua["hierarchy_levels"].get_or(5);

    size_t testpoints = lua["test_points"].get_or(3);
    std::string initial_mesh_filename = lua["initial_mesh_filename"];

    mesh_type initial_mesh;
    disk::netgen_mesh_loader<T, 2> loader;
    loader.verbose(true);
    loader.read_mesh(initial_mesh_filename);
    loader.populate_mesh(initial_mesh);

    disk::mesh_hierarchy<scalar_type> hierarchy(initial_mesh, hierarchy_levels);

    auto mn = hierarchy.meshes_size() - 1;
    auto finer_mesh = *std::next(hierarchy.meshes_begin(), mn);

    if ( lua["dump_last"].get_or(0) )
        dump_netgen_format(finer_mesh, "finer_mesh.mesh2d");

    auto monoscale_solution = 
        load_monoscale_solution(finer_mesh, "../diffusion/solution.bin");

    std::cout <<"done" << std::endl;

    std::array<size_t, 6> levels;

    levels[0] = lua["levels0"].get_or(4);
    levels[1] = lua["levels1"].get_or(4);
    levels[2] = lua["levels2"].get_or(4);
    levels[3] = lua["levels3"].get_or(4);
    levels[4] = lua["levels4"].get_or(4);
    levels[5] = lua["levels5"].get_or(4);

    if ( lua["do0"].get_or(0) )
        test_full_problem_error(lua, *(hierarchy.meshes_begin()+0), levels[0], hierarchy, monoscale_solution);
    
    if ( lua["do1"].get_or(0) )
        test_full_problem_error(lua, *(hierarchy.meshes_begin()+1), levels[1], hierarchy, monoscale_solution);
    
    if ( lua["do2"].get_or(0) )
        test_full_problem_error(lua, *(hierarchy.meshes_begin()+2), levels[2], hierarchy, monoscale_solution);
    
    if ( lua["do3"].get_or(0) )
        test_full_problem_error(lua, *(hierarchy.meshes_begin()+3), levels[3], hierarchy, monoscale_solution);
    
    if ( lua["do4"].get_or(0) )
        test_full_problem_error(lua, *(hierarchy.meshes_begin()+4), levels[4], hierarchy, monoscale_solution);
    
    //if ( lua["do5"].get_or(0) )
    //    test_full_problem_error(lua, *(hierarchy.meshes_begin()+5), levels[5], hierarchy, monoscale_solution);

    //test_full_problem_error(lua, *(hierarchy.meshes_begin()+0), 6, hierarchy, monoscale_solution);
    //test_full_problem_error(lua, *(hierarchy.meshes_begin()+1), 5, hierarchy, monoscale_solution);
    //test_full_problem_error(lua, *(hierarchy.meshes_begin()+2), 4, hierarchy, monoscale_solution);
    //test_full_problem_error(lua, *(hierarchy.meshes_begin()+3), 4, hierarchy, monoscale_solution);
    //test_full_problem_error(lua, *(hierarchy.meshes_begin()+4), 4, hierarchy, monoscale_solution);
    //test_full_problem_error(lua, *(hierarchy.meshes_begin()+5), 3, hierarchy, monoscale_solution);
    //test_full_problem_error(lua, *(hierarchy.meshes_begin()+6), 2, hierarchy, monoscale_solution);


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
    size_t testpoints = lua["test_points"].get_or(3);


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

#ifdef WITH_STATIC_CONDENSATION
    auto system_size = /*num_cells*num_cell_dofs +*/ num_internal_faces*num_face_dofs;
#else
    auto system_size = num_cells*num_cell_dofs + num_internal_faces*num_face_dofs;
#endif

    sparse_matrix_type              A(system_size, system_size);
    dynamic_vector<scalar_type>     b(system_size), x(system_size);

    b = dynamic_vector<scalar_type>::Zero(system_size);

    typedef Eigen::Triplet<scalar_type>         triplet_type;
    std::vector<triplet_type>                   triplets;

#ifdef WITH_STATIC_CONDENSATION
    std::vector<size_t> face_compress_map( num_faces );
    std::vector<size_t> face_expand_map( num_internal_faces );
    size_t fn = 0, fi = 0;
    for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++, fi++)
    {
        //std::cout << *itor << std::endl;
        if( msh.is_boundary(*itor) )
            continue;

        face_compress_map.at(fi) = fn;
        face_expand_map.at(fn) = fi;
        fn++;
    }

    //std::cout << face_compress_map << std::endl;
    //std::cout << face_expand_map << std::endl;
#endif

    //std::ofstream mat_ofs("sysmat.dat");
    //std::ofstream rhs_ofs("sysrhs.dat");

    auto f = [](const typename mesh_type::point_type& pt) ->
                    typename mesh_type::point_type::value_type
    {
        //return 1.;
        //return 2 * M_PI * M_PI * sin(pt.x() * M_PI) * sin(pt.y() * M_PI);
        return sin(pt.x())*sin(pt.y());
    };

    std::vector<dynamic_matrix<scalar_type>> gr_opers, gr_datas;
    gr_opers.reserve( msh.cells_size() );
    gr_datas.reserve( msh.cells_size() );

    size_t elemnum = 0;
    for (auto& cl : msh)
    {
        std::cout << "Assembly: " << elemnum << "/" << msh.cells_size() << "\r";
        std::cout.flush();

        //std::cout << "MLB" << std::endl;
        disk::multiscale_local_basis<mesh_type> mlb(msh, cl, ccb, cfb, k_inner, rl);
        //std::cout << "GR" << std::endl;
        disk::gradient_reconstruction_multiscale<mesh_type, decltype(ccb), decltype(cfb)> gradrec(msh, cl, mlb, ccb, cfb);
        gr_opers.push_back( gradrec.oper );
        gr_datas.push_back( gradrec.data );
        //std::cout << "XX" << std::endl;
        auto proj = disk::make_projector(msh, cl, ccb, cfb);

        dynamic_vector<scalar_type> dofs = make_rhs(msh, cl, ccb, cfb, f);
        dynamic_vector<scalar_type> cdofs = dofs.head(num_cell_dofs);
        
        auto fcs = faces(msh, cl);
        //std::sort(fcs.begin(), fcs.end());
#ifdef WITH_STATIC_CONDENSATION
        std::vector<size_t> l2g(/*num_cell_dofs +*/ fcs.size() * num_face_dofs);
#else
        std::vector<size_t> l2g(num_cell_dofs + fcs.size() * num_face_dofs);

        for (size_t cell_i = 0; cell_i < num_cell_dofs; cell_i++)
            l2g[cell_i] = elemnum * num_cell_dofs + cell_i;
#endif
        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];

            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;
#ifdef WITH_STATIC_CONDENSATION
            auto face_offset = /*num_cells * num_cell_dofs +*/ face_compress_map.at(face_id) * num_face_dofs;
#else
            auto face_offset = num_cells * num_cell_dofs + face_compress_map.at(face_id) * num_face_dofs;
#endif
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

    std::cout << std::endl;

    //mat_ofs.close();
    //rhs_ofs.close();

    A.setFromTriplets(triplets.begin(), triplets.end());
/*
#ifdef HAVE_INTEL_MKL
    Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
#else
    Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
#endif

    solver.analyzePattern(A);
    solver.factorize(A);
    x = solver.solve(b);
*/

    conjugated_gradient(A, b, x);

//    std::cout << b.transpose() << std::endl;
//    std::cout << x.transpose() << std::endl;

    //dynamic_vector<scalar_type> x1 = x;
    //x1(0) = 7;
    //dynamic_vector<scalar_type> b1 = A*x1;

    //std::cout << "NORM: " << (b1-b).norm() << std::endl;

#ifdef PLOT_SOLUTION
    std::ofstream sol_ofs("solution.dat");
    std::ofstream sol_ofs_ms("solution_ms.dat");
#endif


    //disk::simplicial_mesh<T,2> ref_sol_msh;
    //disk::netgen_mesh_loader<T, 2> ref_loader;
    //ref_loader.read_mesh("ref_sol_mesh.mesh");
    //ref_loader.populate_mesh(ref_sol_msh);


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
            //std::cout << " - " << fc << " ";
            if (msh.is_boundary(fc))
            {
                //std::cout << "skipped" << std::endl;
                continue;
            }
            //std::cout << std::endl;

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
        //disk::gradient_reconstruction_multiscale<mesh_type, decltype(ccb), decltype(cfb)> gradrec(msh, cl, mlb, ccb, cfb);

        dynamic_vector<scalar_type> rhs_c = rhs.head(num_cell_dofs);
        dynamic_vector<scalar_type> local_dofs = recover_static_condensation(gr_datas.at(elemnum), rhs_c, dofs, num_cell_dofs);

        dynamic_vector<scalar_type> R = gr_opers.at(elemnum) * local_dofs;
        R = R.head(R.size()-1);

        //std::cout << "Elem dofs" << std::endl;
        //std::cout << R.transpose() << std::endl;

        auto ms_inner_mesh = mlb.inner_mesh();
        for (auto& icl : ms_inner_mesh)
        {
            auto itps = make_test_points(ms_inner_mesh, icl, testpoints);
            for(auto& tp : itps)
            {
                auto msphi = mlb.eval_functions(icl, tp);
                auto msval = R.dot(msphi);
#ifdef PLOT_SOLUTION
                sol_ofs_ms << tp.x() << " " << tp.y() << " " << msval << std::endl;
#endif
            }
        }
//#endif
        elemnum++;
    }
#ifdef PLOT_SOLUTION
    sol_ofs.close();
    sol_ofs_ms.close();
#endif

}

int
main(int argc, char **argv)
{
    sol::state lua;
    lua.do_file("params.lua");
    
    /*
    if ( !std::regex_match(argv[1], std::regex(".*\\.mesh2d$") ) )
    {
        std::cout << "Mesh format must be Netgen 2D" << std::endl;
        return 1;
    }

    auto msh = disk::load_netgen_2d_mesh<double>(argv[1]);

    std::cout << "AVG H: " << average_diameter(msh) << std::endl;

    

    std::string mode = lua["mode"].get_or( std::string("testgr") );
    std::cout << "Mode: " << mode << std::endl;

    if (mode == "testgr")
        test_gradient_reconstruction(lua, msh);
    else if (mode == "testfull")
        test_full_problem(lua, msh);
    else if (mode == "pms")
        plot_multiscale_basis_functions(lua, msh);
    else if (mode == "run_ms_test")
    */
    
    run_multiscale_tests(lua);

    return 0;
}


