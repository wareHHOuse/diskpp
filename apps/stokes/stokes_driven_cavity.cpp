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

#include "bases/bases.hpp"
#include "quadratures/quadratures.hpp"
#include "methods/hho"

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
    typedef point<typename Mesh::coordinate_type,2> point_type;
    //Warningn this may work only on convex polygons
    auto h = barycenter(msh, cl);
    auto sorted_pts = pts;

    std::sort( sorted_pts.begin(), sorted_pts.end(),
                            [&](const point<typename Mesh::coordinate_type,2>& va, const point_type& vb )
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
          const point<typename Mesh::coordinate_type,2>& P)
{
    auto vts = points(msh, cl);

    typedef  typename Mesh::coordinate_type scalar_type;

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
                const std::pair<point<typename Mesh::coordinate_type,2>,
                                point<typename Mesh::coordinate_type,2>>   & e,
                const Matrix<typename Mesh::coordinate_type, Dynamic, 1> & sol,
                const typename disk::hho_degree_info& hdi,
                const std::string & filename)
{
    typedef Matrix<typename Mesh::coordinate_type, Dynamic, 1>  vector_type;
    typedef point<typename Mesh::coordinate_type,2>             point_type;

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
                auto cbs   = disk::vector_basis_size(hdi.cell_degree(), dim, dim);
                auto cell_ofs = disk::priv::offset(msh, cl);

                vector_type s = sol.block(cell_ofs * cbs, 0, cbs, 1);
                auto cb = disk::make_vector_monomial_basis(msh,
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
                        const dynamic_vector< typename Mesh::coordinate_type>& sol,
                        const typename disk::hho_degree_info& hdi,
                        const std::string& filename)
{
    typedef Mesh mesh_type;
    typedef typename Mesh::coordinate_type scalar_type;
    const size_t cell_degree   = hdi.cell_degree();
    const size_t cbs = disk::vector_basis_size(cell_degree,
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
        auto cell_ofs = disk::priv::offset(msh, cl);
        Matrix<scalar_type, Dynamic, 1> x = sol.block(cell_ofs * cbs, 0, cbs, 1);

        const auto cell_nodes = post_mesh.nodes_cell(cell_i);
        auto  cbas = disk::make_vector_monomial_basis(msh, cl, cell_degree);

       // Loop on the nodes of the cell
        for (auto& point_id : cell_nodes)
        {
            const auto pt = storage->points[point_id];

            const auto phi = cbas.eval_functions(pt);
            assert(phi.rows() == cbs);
            const auto depl = disk::eval(x, phi);

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

template<typename Mesh, typename Assembler>
void
post_processing(const Mesh& msh, const Assembler& assembler,
                const typename disk::hho_degree_info& di,
                const dynamic_vector<typename Mesh::coordinate_type>& sol,
                bool use_sym_grad = true)

{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, 1>  vector_type;
    typedef typename Mesh::point_type  point_type;

    size_t const dim = Mesh::dimension;
    auto rbs = disk::vector_basis_size(di.reconstruction_degree(), dim, dim);
    auto cbs = disk::vector_basis_size(di.cell_degree(), dim, dim);

    std::ofstream ofs("data.data");
    if (!ofs.is_open())
        std::cout << "Error opening file"<<std::endl;

    for(auto cl : msh)
    {

        auto gr  = disk::make_hho_stokes(msh, cl, di, use_sym_grad);
        vector_type svel =  assembler.take_velocity(msh, cl, sol);

        //this is only for k = 0 (for k >0, it is ok for plotting  but not for div u )
        auto bar = barycenter(msh, cl);

        //Velocity
        vector_type cell_vel = svel.block(0,0, cbs, 1);
        auto cb  = disk::make_vector_monomial_basis(msh, cl, di.cell_degree());
        auto phi = cb.eval_functions(bar);
        Matrix<T, Mesh::dimension, 1> ueval = disk::eval(cell_vel, phi);

        //Stress
        Matrix<T, Dynamic,  DIM> c_dphi;
        if(use_sym_grad)
        {
            auto cb = disk::make_sym_matrix_monomial_basis(msh, cl, di.reconstruction_degree());
            Matrix<T, Dynamic,  DIM> c_dphi_tmp = cb.eval_gradients(bar);
            c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
        }
        else
        {
            auto cb = disk::make_sym_matrix_monomial_basis(msh, cl, di.reconstruction_degree());
            Matrix<T, Dynamic,  DIM> c_dphi_tmp = cb.eval_gradients(bar);
            c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
        }

        auto gr = disk::make_hho_stokes(msh, cl, hdi, use_sym_grad);
        Matrix<T, Dynamic, 1> Gu = g.first * svel;
        Matrix<T, dim, dim> grad_eval = disk::eval(Gu, c_dphi);
        T divu = grad_eval(0,0) + grad_eval(1,1);
        ofs << ueval(0)   << " " << ueval(1) << " " << divu <<std::endl;
    }
    ofs.close();

    compute_discontinuous_velocity( msh, sol, di, "depl2d.msh");
    auto p_x = std::make_pair(point_type({0.0, 0.5}), point_type({1.0, 0.5}));
    auto p_y = std::make_pair(point_type({0.5, 0.0}), point_type({0.5, 1.0}));
    plot_over_line(msh, p_x, sol, di, "plot_over_line_x.data");
    plot_over_line(msh, p_y, sol, di, "plot_over_line_y.data");

    return;
}

template<typename Mesh>
auto
run_stokes(const Mesh& msh, size_t degree, bool use_sym_grad = true)
{
    typedef Mesh mesh_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    typedef dynamic_matrix<scalar_type>     matrix_type;
    typedef disk::vector_boundary_conditions<mesh_type> boundary_type;

    auto rhs_fun = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
        return Matrix<scalar_type, 2, 1>::Zero();
    };

    auto wall = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
        return Matrix<scalar_type, 2, 1>::Zero();
    };
    auto movingWall = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
        return Matrix<scalar_type, 2, 1>{1,0};
    };

    typename disk::hho_degree_info hdi(degree, degree);
    boundary_type bnd(msh);

    bnd.addDirichletBC(0, 1, movingWall);
    bnd.addDirichletBC(0, 2, wall);
    bnd.addDirichletBC(0, 3, wall);
    bnd.addDirichletBC(0, 4, wall);

    auto assembler = disk::make_stokes_assembler(msh, hdi, bnd);

    auto viscosity = scalar_type(100);
    scalar_type factor = (use_sym_grad)? 2. * viscosity : viscosity;

    for (auto cl : msh)
    {
        auto gr = disk::make_hho_stokes(msh, cl, hdi, use_sym_grad);
        Matrix<scalar_type, Dynamic, Dynamic> stab;
        stab = make_vector_hho_stabilization(msh, cl, gr.first, hdi);
        auto dr = make_hho_divergence_reconstruction_rhs(msh, cl, hdi);
        auto cb = disk::make_vector_monomial_basis(msh, cl, hdi.cell_degree());
        auto rhs = make_rhs(msh, cl, cb, rhs_fun);
        assembler.assemble(msh, cl, factor * (gr.second + stab), -dr, rhs);
    }

    assembler.finalize();

    //dump_sparse_matrix(assembler.LHS, "stokes.txt");

    size_t systsz = assembler.LHS.rows();
    size_t nnz = assembler.LHS.nonZeros();

    dynamic_vector<scalar_type> sol = dynamic_vector<scalar_type>::Zero(systsz);

    disk::solvers::pardiso_params<scalar_type> pparams;
    mkl_pardiso_ldlt(pparams, assembler.LHS, assembler.RHS, sol);

    post_processing(msh, assembler, hdi, sol, use_sym_grad);

    return;

}

//#if 0
int main(int argc, char **argv)
{
    using RealType = double;
    bool use_sym_grad = false;

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

        run_stokes(msh, degree, use_sym_grad);
        std::cout << "fini" << std::endl;
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

        run_stokes(msh, degree, use_sym_grad);
    }

    /* Medit 2d*/
    if (std::regex_match(filename, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        typedef disk::generic_mesh<RealType, 2>  mesh_type;
        mesh_type msh = disk::load_medit_2d_mesh<RealType>(filename);
        run_stokes(msh, degree, use_sym_grad);
    }

    return 0;
}
//#endif
