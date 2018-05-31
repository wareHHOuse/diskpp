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
 * Karol Cascavita (C) 2018                     karol.cascavita@enpc.fr
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
#include <iomanip>
#include <regex>

#include <unistd.h>

#include "revolution/bases"
#include "revolution/quadratures"
#include "revolution/methods/hho"

#include "core/loaders/loader.hpp"

#include "solvers/solver.hpp"

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"
#include "output/postMesh.hpp"

#include "output/silo.hpp"

// Return polar angle of p with respect to origin o
template<typename T>
auto
to_angle(const point<T,2>& p, const point<T,2>& o)
{
  return atan2(p.y() - o.y(), p.x() - o.x());
}
template<typename Mesh, typename Points>
auto
sort_by_polar_angle(const Mesh & msh,
                    const typename Mesh::cell &  cl,
                    const Points& pts)
{
    typedef point<typename Mesh::scalar_type,2> point_type;
    //Warningn this may work only on convex polygons
    auto h = barycenter(msh, cl);
    auto sorted_pts = pts;

    std::sort( sorted_pts.begin(), sorted_pts.end(),
                            [&](const point<typename Mesh::scalar_type,2>& va, const point_type& vb )
        {
            auto theta_a = to_angle(va, h);
            auto theta_b = to_angle(vb, h);
            return (theta_a < theta_b);
        }
    );
    return sorted_pts;
}
/*
* Copyright 2000 softSurfer, 2012 Dan Sunday
* This code may be freely used and modified for any purpose
* providing that this copyright notice is included with it.
* SoftSurfer makes no warranty for this code, and cannot be held
* liable for any real or imagined damage resulting from its use.
* Users of this code must verify correctness for their application.
*===================================================================
* isLeft(): tests if a point is Left|On|Right of an infinite line.
*    Input:  three points P0, P1, and P2
*    Return: >0 for P2 left of the line through P0 and P1
*            =0 for P2  on the line
*            <0 for P2  right of the line
*    See: Algorithm 1 "Area of Triangles and Polygons"
*===================================================================
*/
template<typename T>
auto
isLeft( const point<T,2>& P0, const point<T,2>& P1, const point<T,2>& P2 )
{
    auto ret = ( (P1.x() - P0.x()) * (P2.y() - P0.y())
            - (P2.x() -  P0.x()) * (P1.y() - P0.y()) );
    return ret;
}
/*===================================================================
* wn_PnPoly(): winding number test for a point in a polygon
*      Input:   P = a point,
*               vts = vertex points of a polygon vts[n+1] with vts[n]=V[0]
*      Return:  the winding number (= 0 only when P is outside)
*===================================================================
*/
template<typename Mesh>
bool
wn_PnPoly(const Mesh& msh,
          const typename Mesh::cell& cl,
          const point<typename Mesh::scalar_type,2>& P)
{
    auto vts = points(msh, cl);

    typedef  typename Mesh::scalar_type scalar_type;

    std::vector< point<scalar_type,2>> svts(vts.size()+1);
    auto svts_temp = sort_by_polar_angle(msh, cl, vts);
    std::copy(std::begin(svts_temp), std::end(svts_temp), std::begin(svts));
    svts[vts.size()] = svts_temp[0];

    int winding_number = 0;

    for (int i = 0 ; i < svts.size() - 1; i++)
    {
        if (svts.at(i).y() <= P.y() )
        {
            auto upward_crossing = svts[i+1].y()  > P.y();
            if (upward_crossing)
            {
                auto P_on_left = isLeft( svts[i], svts[i+1], P) > scalar_type(0);
                if (P_on_left)
                    ++winding_number; //valid up intersect
            }
        }
        else
        {
            // start y > P.y (no test needed)
            auto downward_crossing = svts[i+1].y()  <= P.y();
            if (downward_crossing)
            {
                auto P_on_right = isLeft( svts[i], svts[i+1], P) < scalar_type(0);
                if (P_on_right)
                    --winding_number;  //valid down intersect
            }
        }
    }
    if(winding_number != 0)
        return true;
    return false;
}
template<typename Mesh>
void
plot_over_line(const Mesh    & msh,
                const std::pair<point<typename Mesh::scalar_type,2>,
                                point<typename Mesh::scalar_type,2>>   & e,
                const Matrix<typename Mesh::scalar_type, Dynamic, 1> & sol,
                const typename revolution::hho_degree_info& hdi,
                const std::string & filename)
{
    typedef Matrix<typename Mesh::scalar_type, Dynamic, 1>  vector_type;
    typedef point<typename Mesh::scalar_type,2>             point_type;

    std::ofstream pfs(filename);
    if(!pfs.is_open())
        std::cout << "Error opening file :"<< filename <<std::endl;

    auto N = 400;
    auto h = (e.second - e.first)/N;
    size_t i = 0 ;
    auto pts = std::vector<point_type>(N);
    auto dim = Mesh::dimension;
    for(auto& p: pts)
    {
        p = e.first + i * h;

        for(auto cl: msh)
        {
            if(wn_PnPoly( msh, cl, p))
            {
                auto cbs   = revolution::vector_basis_size(hdi.cell_degree(), dim, dim);
                auto cell_ofs = revolution::priv::offset(msh, cl);

                vector_type s = sol.block(cell_ofs * cbs, 0, cbs, 1);
                auto cb = revolution::make_vector_monomial_basis(msh,
                                                        cl, hdi.cell_degree());
                auto phi = cb.eval_functions(p);
                vector_type vel = phi.transpose() * s;
                pfs<< p.x() << " "<< p.y() << " "<< vel(0) << " "<< vel(1)<< std::endl;
            }
        }
        i++;
    }
    pfs.close();
    return;
}

template<typename Mesh>
void
compute_discontinuous_velocity(const Mesh& msh,
                        const dynamic_vector< typename Mesh::scalar_type>& sol,
                        const typename revolution::hho_degree_info& hdi,
                        const std::string& filename)
{
    typedef Mesh mesh_type;
    typedef typename Mesh::scalar_type scalar_type;
    const size_t cell_degree   = hdi.cell_degree();
    const size_t cbs = revolution::vector_basis_size(cell_degree,
                                    Mesh::dimension, Mesh::dimension);
    // compute mesh for post-processing
    disk::PostMesh<mesh_type> post_mesh = disk::PostMesh<mesh_type>(msh);
    gmsh::Gmesh gmsh    = disk::convertMesh(post_mesh);
    auto        storage = post_mesh.mesh().backend_storage();

    const static_vector<scalar_type, Mesh::dimension> vzero =
                            static_vector<scalar_type, Mesh::dimension>::Zero();

    const size_t nb_nodes(gmsh.getNumberofNodes());

    // first(number of data at this node), second(cumulated value)
    std::vector<std::pair<size_t, static_vector<scalar_type, Mesh::dimension>>>
                                    value(nb_nodes, std::make_pair(0, vzero));

    size_t cell_i = 0;

    for (auto& cl : msh)
    {
        auto cell_ofs = revolution::priv::offset(msh, cl);
        Matrix<scalar_type, Dynamic, 1> x = sol.block(cell_ofs * cbs, 0, cbs, 1);

        const auto cell_nodes = post_mesh.nodes_cell(cell_i);
        auto  cbas = revolution::make_vector_monomial_basis(msh, cl, cell_degree);

       // Loop on the nodes of the cell
        for (auto& point_id : cell_nodes)
        {
            const auto pt = storage->points[point_id];

            const auto phi = cbas.eval_functions(pt);
            assert(phi.rows() == cbs);
            const auto depl = revolution::eval(x, phi);

            // Add displacement at node
            value[point_id].first++;
            value[point_id].second += depl;
        }
        cell_i++;
    }

    std::vector<gmsh::Data>    data;    // create data
    std::vector<gmsh::SubData> subdata; // create subdata
    data.reserve(nb_nodes);             // data has a size of nb_node

    // Compute the average value and save it
    for (size_t i_node = 0; i_node < value.size(); i_node++)
    {
        const static_vector<scalar_type, Mesh::dimension> depl_avr =
                            value[i_node].second / double(value[i_node].first);

        const gmsh::Data tmp_data(i_node + 1, disk::convertToVectorGmsh(depl_avr));
        data.push_back(tmp_data);
    }

    // Create and init a nodedata view
    gmsh::NodeData nodedata(3, 0.0, "depl_node_cont", data, subdata);
    // Save the view
    nodedata.saveNodeData(filename, gmsh);

    return;
}

template<typename Mesh>
class augmented_lagrangian_stokes
{
    typedef Mesh mesh_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;
    typedef typename mesh_type::scalar_type scalar_type;

    typedef disk::mechanics::BoundaryConditions<mesh_type> boundary_type;

    using point_type = typename mesh_type::point_type;

    typedef dynamic_vector<scalar_type>       vector_type;
    typedef dynamic_matrix<scalar_type>       matrix_type;

    typedef Matrix<scalar_type, 2, 2>         tensor_type;
    typedef Matrix<scalar_type, 2, 1>         vector2d_type;
    typedef std::function<vector2d_type (const point_type &)> vector_funtion_type;
    typedef std::function<scalar_type   (const point_type &)> scalar_funtion_type;

    vector_funtion_type     velocity, rhs_fun;
    scalar_funtion_type     pressure;
    vector_type             multiplier, auxiliar;

    typename revolution::hho_degree_info di;
    bool                    use_sym_grad;
    scalar_type             factor;
    scalar_type             viscosity;
    scalar_type             alpha;
    size_t                  cbs, fbs, pbs, sbs;

public:
    vector_type             sol, sol_old;
    scalar_type             convergence;

    augmented_lagrangian_stokes(const Mesh& msh,
                            const typename revolution::hho_degree_info & hdi,
                            const bool& use_symmetric_gradient):
                            di(hdi), use_sym_grad(use_symmetric_gradient)
    {
        rhs_fun  = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
            return Matrix<scalar_type, 2, 1>::Zero();
        };

        velocity = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
            return Matrix<scalar_type, 2, 1>::Zero();
        };
        pressure = [](const point_type& p) -> scalar_type {
            return 0;
        };

        factor = 1.;//(use_sym_grad)? 2. : 1.;
        viscosity = 100;
        alpha = 0.1;
        convergence = 0.;
        auto dim =  Mesh::dimension;

        cbs = revolution::vector_basis_size(di.cell_degree(), dim, dim);
        fbs = revolution::vector_basis_size(di.face_degree(), dim - 1, dim);
        pbs = revolution::scalar_basis_size(di.face_degree(), dim);
        sbs = revolution::matrix_basis_size(di.face_degree(), dim, dim);
    };

    auto
    define_assembler(const mesh_type& msh)
    {
        boundary_type bnd(msh);

        auto wall = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
            return Matrix<scalar_type, 2, 1>::Zero();
        };
        auto movingWall = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
            return Matrix<scalar_type, 2, 1>{1,0};
        };

        bnd.addDirichletBC(0, 1, movingWall);
        bnd.addDirichletBC(0, 2, wall);
        bnd.addDirichletBC(0, 3, wall);
        bnd.addDirichletBC(0, 4, wall);

        auto assembler = revolution::make_stokes_assembler_alg(msh, di, bnd);

        return assembler;
    }

    template<typename Assembler>
    auto
    initialize(const mesh_type& msh, const Assembler& assembler)
    {
        auto systsz = assembler.global_system_size();
        sol = vector_type::Zero(systsz);
        sol_old = vector_type::Zero(systsz);

        multiplier = vector_type::Zero(msh.cells_size() * sbs);
        auxiliar   = vector_type::Zero(msh.cells_size() * sbs);

        return;
    }

    template<typename Assembler>
    auto
    compute_errors( const mesh_type& msh,
                    const Assembler& assembler)
    {
        auto dim =  Mesh::dimension;

        scalar_type error_vel(0);

        for (auto& cl : msh)
        {
        	auto bar = barycenter(msh, cl);

            //energy error
            Matrix<scalar_type, Dynamic, 1> svel =  assembler.take_velocity(msh, cl, sol);
            Matrix<scalar_type, Dynamic, 1> svel_old =  assembler.take_velocity(msh, cl, sol_old);

            Matrix<scalar_type, Dynamic, 1> diff_vel = svel - svel_old;
            auto gr = revolution::make_hho_stokes(msh, cl, di, use_sym_grad);
            Matrix<scalar_type, Dynamic, Dynamic> stab;
            stab = make_hho_fancy_stabilization_vector(msh, cl, gr.first, di);
            auto G = revolution::make_hlow_stokes(msh, cl, di, use_sym_grad);

            Matrix<scalar_type, Dynamic, Dynamic> B = factor * viscosity * (G.second + stab);

            error_vel += diff_vel.dot(B * diff_vel);
        }

        //ofs.close();
        return std::make_pair(std::sqrt(error_vel), 0.);
    }

    template<typename Assembler>
    void
    save_auxiliar(  const mesh_type& msh,
                        const cell_type& cl,
                        const Assembler& assembler)
    {
        auto gamma = compute_auxiliar(msh, cl, assembler);
        auto cell_ofs  = revolution::priv::offset(msh, cl);
        auxiliar.block(cell_ofs * sbs, 0, sbs, 1) = gamma;
        return;
    }

    template<typename Assembler>
    Matrix<scalar_type, Dynamic, 1>
    compute_auxiliar(   const mesh_type& msh,
                        const cell_type& cl,
                        const Assembler& assembler)
    {
        auto u_TF  = assembler.take_velocity(msh, cl, sol_old);
        auto value = 1./ (factor *(viscosity + alpha));
        auto G = revolution::make_hlow_stokes(msh, cl, di, use_sym_grad);

        auto cell_ofs    = revolution::priv::offset(msh, cl);
        Matrix<scalar_type, Dynamic, 1> Gu = G.first * u_TF;
        Matrix<scalar_type, Dynamic, 1> stress = multiplier.block(cell_ofs * sbs, 0, sbs, 1);
        Matrix<scalar_type, Dynamic, 1> gamma = value * (stress +  factor * alpha * Gu);

        return gamma;
    }

    template<typename Assembler>
    auto
    update_multiplier(const mesh_type& msh, const Assembler& assembler)
    {
        convergence = 0.;
        auto dim = Mesh::dimension;

        for(auto cl: msh)
        {
            auto u_TF = assembler.take_velocity(msh, cl, sol);
            auto G = revolution::make_hlow_stokes(msh, cl, di, use_sym_grad);
            auto cell_ofs = revolution::priv::offset(msh, cl);

            Matrix<scalar_type, Dynamic, 1> Gu = G.first * u_TF;
            Matrix<scalar_type, Dynamic, 1> gamma = compute_auxiliar( msh,  cl, assembler);
            Matrix<scalar_type, Dynamic, 1> gamma_old = auxiliar.block(cell_ofs *sbs, 0, sbs, 1);

            Matrix<scalar_type, Dynamic, 1> diff_stress = factor * alpha * (Gu - gamma);
            Matrix<scalar_type, Dynamic, 1> diff_gamma  = alpha * (gamma - gamma_old);

            auto sb =  revolution::make_matrix_monomial_basis(msh, cl, di.face_degree());
            Matrix<scalar_type, Dynamic, Dynamic> mass = make_mass_matrix(msh, cl, sb);

            convergence += diff_stress.dot(mass * diff_stress);// + diff_gamma.dot(mass * diff_gamma);

            multiplier.block(cell_ofs * sbs, 0, sbs, 1) +=   diff_stress;
        }

        return;
    }

    template<typename Assembler>
    Matrix<scalar_type, Dynamic, 1>
    make_rhs_alg(   const mesh_type& msh,
                    const cell_type& cl,
                    const Assembler& assembler)
    {
        auto G = revolution::make_hlow_stokes(msh, cl, di, use_sym_grad);
        auto cb = revolution::make_vector_monomial_basis(msh, cl, di.cell_degree());
        auto sb = revolution::make_matrix_monomial_basis(msh, cl, di.face_degree());
        auto cell_ofs =  revolution::priv::offset(msh, cl);
        auto num_faces = howmany_faces(msh, cl);

        auto stress = multiplier.block( sbs * cell_ofs,  0, sbs, 1);
        auto gamma  = compute_auxiliar( msh,  cl, assembler);
        vector_type str_agam = stress - factor * alpha * gamma;

        Matrix<scalar_type, Dynamic, Dynamic> mm = revolution::make_mass_matrix(msh, cl, sb);

        Matrix<scalar_type, Dynamic, 1> rhs =
                    Matrix<scalar_type, Dynamic, 1>::Zero(cbs + fbs * num_faces);

        //(f, v_T)
        rhs.block( 0, 0, cbs, 1) = make_rhs(msh, cl, cb, rhs_fun);

        //(stress - factor * alpha * gamma, Gv)
        rhs -=  G.first.transpose() * mm * str_agam;

        return rhs;
    }

    template<typename Assembler>
    auto
    make_global_rhs(const mesh_type& msh, Assembler& assembler)
    {
        assembler.initialize_rhs();

        for (auto cl : msh)
        {
            Matrix<scalar_type, Dynamic, 1> local_rhs = make_rhs_alg(msh, cl, assembler);

            assembler.assemble_rhs(msh, cl, local_rhs);

            save_auxiliar(msh, cl, assembler);
        }
        assembler.finalize_rhs();

        return;
    }

    template<typename Assembler>
    auto
    make_global_matrix(const mesh_type& msh, Assembler& assembler)
    {
        auto celdeg = di.cell_degree();

        assembler.initialize_lhs();

        for (auto cl : msh)
        {
            auto cb  = revolution::make_vector_monomial_basis(msh, cl, celdeg);
            auto G = revolution::make_hlow_stokes(msh, cl, di, use_sym_grad);
            auto gr  = revolution::make_hho_stokes(msh, cl, di, use_sym_grad);
            Matrix<scalar_type, Dynamic, Dynamic> stab;
            stab = make_hho_fancy_stabilization_vector(msh, cl, gr.first, di);
            auto dr  = make_hho_divergence_reconstruction_stokes_rhs(msh, cl, di);

            Matrix<scalar_type, Dynamic, Dynamic> A = factor * (alpha * G.second + viscosity * stab);

            assembler.assemble_lhs(msh, cl, A, -dr.second);
        }

        assembler.finalize_lhs();

        return;
    }

    template<typename Assembler>
    auto
    run_stokes_like(const mesh_type& msh, Assembler& assembler, const size_t iter)
    {
        sol_old = sol;

        make_global_rhs(msh,assembler);

        if(iter == 0)
            make_global_matrix(msh, assembler);

        //dump_sparse_matrix(assembler.LHS, "stokes.txt");
        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        sol = vector_type::Zero(systsz);
        disk::solvers::pardiso_params<scalar_type> pparams;
        mkl_pardiso_ldlt(pparams, assembler.LHS, assembler.RHS, sol);
        //mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

        return;
    }

};


template<typename Mesh, typename T, typename Assembler>
void
quiver( const Mesh& msh, const dynamic_vector<T>& sol, const Assembler& assembler,
        const typename revolution::hho_degree_info & di)
{
    std::ofstream ofs("quiver_vel.data");

    if (!ofs.is_open())
        std::cout << "Error opening errors "<<std::endl;

    auto cbs = revolution::vector_basis_size(di.cell_degree(),
                                            Mesh::dimension, Mesh::dimension);
    for (auto& cl: msh)
    {
        auto cell_ofs = revolution::priv::offset(msh, cl);
        Matrix<T, Dynamic, 1>  s = assembler.take_velocity(msh, cl, sol);

        auto cb  = revolution::make_vector_monomial_basis(msh, cl, di.cell_degree());
        auto qps =  revolution::integrate(msh, cl, 2 * di.face_degree());
        for (auto& qp : qps)
        {
            Matrix<T, Dynamic, 2>  phi = cb.eval_functions(qp.point());
            Matrix<T,  1, Dynamic>  st = (s.block(0,0,cbs,1)).transpose();
            Matrix<T, Dynamic, 1>  vt =  st * phi;

            ofs << qp.point().x() << "  "<< qp.point().y()<< "   ";
            ofs << vt(0) << "  "<< vt(1) << "    "<< std::endl;
        }
    }
    ofs.close();
}

template< typename T>
std::string
tostr(const T & var)
{
    std::ostringstream  ostr;
    ostr << var;
    return ostr.str();
}

template<typename Mesh>
auto
run_alg_stokes(const Mesh& msh, size_t degree,// std::ofstream & ofs,
                bool use_sym_grad)
{
    using T = typename Mesh::coordinate_type;
    T tolerance = 1.e-9, Ninf = 10.e+10;
    size_t max_iters = 50000;

    typename revolution::hho_degree_info hdi(degree, degree);
    augmented_lagrangian_stokes<Mesh> als(msh, hdi, use_sym_grad);

    auto assembler = als.define_assembler(msh);
    als.initialize(msh, assembler);

    T error_old = 0.;

    for(size_t i = 0; i < max_iters; i++)
    {
        als.run_stokes_like(msh, assembler, i);
        als.update_multiplier(msh, assembler);
        auto error = als.compute_errors(msh, assembler);
        auto convergence = std::sqrt(als.convergence);

        if(i % 25 == 0)
            std::cout << "  i : "<< i<<"  - " << convergence<<std::endl;
        assert(convergence < Ninf);
        if( convergence < tolerance)
            break;
        //if(error.first < 10.e-8 && i > 1)
        //    break;
    }

    auto final_error = als.compute_errors(msh, assembler);

    auto dim = Mesh::dimension;
    auto cbs   = revolution::vector_basis_size(hdi.cell_degree(), dim, dim);
    Matrix< T, Dynamic, 1> cell_sol(cbs * msh.cells_size());
    for(auto cl : msh)
    {
        auto cell_ofs = revolution::priv::offset(msh, cl);
        Matrix<T, Dynamic, 1> svel =  assembler.take_velocity(msh, cl, als.sol);
        cell_sol.block(cell_ofs * cbs, 0, cbs, 1) = svel.block(0, 0, cbs, 1);
    }

    typedef point<T,2>             point_type;

    compute_discontinuous_velocity( msh, cell_sol, hdi, "driven2d.msh");
    auto p_x = std::make_pair(point_type({0.0, 0.5}), point_type({1.0, 0.5}));
    auto p_y = std::make_pair(point_type({0.5, 0.0}), point_type({0.5, 1.0}));
    plot_over_line(msh, p_x, cell_sol, hdi, "plot_over_line_x_alg.data");
    plot_over_line(msh, p_y, cell_sol, hdi, "plot_over_line_y_alg.data");

    //quiver( msh, als.sol, assembler, hdi);
    return final_error;
}

#if 0
void convergence_test_typ1(void)
{
    using T = double;
    bool use_sym_grad = false;
    std::vector<std::string> meshfiles;

    //meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_1.typ1");
    //meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_2.typ1");
    //meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_3.typ1");
    //meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_4.typ1");
    //meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_5.typ1");
    //meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_6.typ1");

    /*
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_1.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_2.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_3.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_4.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_5.typ1");
    */
    /*
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_1.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_2.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_3.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_4.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_5.typ1");
    */
    std::cout << "                   velocity H1-error";
    std::cout << "    -     pressure L2-error "<< std::endl;

    for (size_t k = 2; k < 3; k++)
    {
        std::cout << "DEGREE " << k << std::endl;

        std::ofstream ofs("errors_k" + tostr(k) + ".data");
        if (!ofs.is_open())
            std::cout << "Error opening errors "<<std::endl;

        std::vector<T> mesh_hs;
        std::vector<std::pair<T,T>> errors;

        //for (size_t i = 0; i < 1; i++)
        for (size_t i = 0; i < meshfiles.size(); i++)
        {
            typedef disk::generic_mesh<T, 2>  mesh_type;

            std::cout << " Mesh : "<< i << std::endl;
            mesh_type msh;
            disk::fvca5_mesh_loader<T, 2> loader;
            if (!loader.read_mesh(meshfiles.at(i)))
            {
                std::cout << "Problem loading mesh." << std::endl;
                continue;
            }
            loader.populate_mesh(msh);

            auto error = run_alg_stokes(msh, k, ofs, use_sym_grad);

            mesh_hs.push_back( disk::mesh_h(msh) );
            errors.push_back(error);
            ofs << " " << std::endl;
        }
        ofs.close();

        for (size_t i = 0; i < mesh_hs.size(); i++)
        {
            if (i == 0)
            {
                std::cout << "    ";
                std::cout << std::scientific << std::setprecision(4) << mesh_hs.at(i) << "    ";
                std::cout << std::scientific << std::setprecision(4) << errors.at(i).first;
                std::cout << "     -- " << "          ";
                std::cout << std::scientific << std::setprecision(4) << errors.at(i).second;
                std::cout << "     -- " << std::endl;
            }
            else
            {
                auto rate = std::log( errors.at(i).first/errors.at(i-1).first ) /
                            std::log( mesh_hs.at(i)/mesh_hs.at(i-1) );
                std::cout << "    ";
                std::cout << std::scientific  << std::setprecision(4) << mesh_hs.at(i) << "    ";
                std::cout << std::scientific  << std::setprecision(4) << errors.at(i).first << "    ";
                std::cout << std::fixed<< std::setprecision(2) << rate << "          ";

                auto pres_rate = std::log( errors.at(i).second/errors.at(i-1).second ) /
                            std::log( mesh_hs.at(i)/mesh_hs.at(i-1) );
                std::cout << std::scientific  << std::setprecision(4) << errors.at(i).second << "    ";
                std::cout << std::fixed << std::setprecision(2) << pres_rate << std::endl;
            }
        }
    }
}

int main(void)
{
    convergence_test_typ1();
}
#endif

//#if 0

int main(int argc, char **argv)
{
    using RealType = double;
    bool use_sym_grad = true;

    char    *filename       = nullptr;
    int ch;
    size_t degree = 0;

    while ( (ch = getopt(argc, argv, "k:")) != -1 )
    {
        switch(ch)
        {
            case 'k':
                degree = atoi(optarg);
                break;

            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    if (argc != 1)
    {
        std::cout << "Please specify a 2D mesh" << std::endl;

        return 0;
    }

    filename = argv[0];

    if (std::regex_match(filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;

        typedef disk::generic_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::fvca5_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        run_alg_stokes(msh, degree, use_sym_grad); // ;ofs, use_sym_grad);
    }

    if (std::regex_match(filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;

        typedef disk::simplicial_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::netgen_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        run_alg_stokes(msh, degree, use_sym_grad); //, ofs, use_sym_grad);
    }

    /* Medit 2d*/
    if (std::regex_match(filename, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        typedef disk::generic_mesh<RealType, 2>  mesh_type;
        mesh_type msh = disk::load_medit_2d_mesh<RealType>(filename);
        run_alg_stokes(msh, degree, use_sym_grad); //, ofs, use_sym_grad);
    }
/*
    if (std::regex_match(filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: Cartesian 2D" << std::endl;

        typedef disk::cartesian_mesh<RealType, 2>   mesh_type;

        mesh_type msh;
        disk::cartesian_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        test_bases(msh);
    }
    */

    return 0;
}
//#endif
